#include "TimeIntegration.hpp"
#include "MaterialPoints.hpp"
#include "Nodes.hpp"
#include "eig3.h"
#include "Matrix.h"
#include <math.h>
#include "omp.h"
#include <algorithm>
#include <stdio.h>

CPUGPU void updateMatPoint(Parameters *mpmParameters,MaterialPoints *mpts, Nodes *nds,double dt, std::vector<MaterialModels>* materials,int mp){
// Update the density, strain, stress, position and velocity
	// Store some nodal values for the loops
		int Nx          = nds->Nx-1;
		int Ny          = nds->Ny-1;
		int Nz          = nds->Nz-1;
		double cellsize = nds->cellsize;
		double xmin     = nds->position[0*3 + 0];
		double ymin     = nds->position[0*3 + 1];
		double zmin     = nds->position[0*3 + 2];
			// Update Density
			mpts->density[mp]              += mpts->densityRate[mp]*dt;
			mpts->volume[mp]               = mpts->mass[mp]/mpts->density[mp];
				// Update Velocity
				double dampingCoef =  mpmParameters->velocityDampingCoef;
				if(mpts->hold[mp]==false)
				{
				mpts->velocity[mp*3+0]         = mpts->velocity[mp*3+0] + (mpts->acceleration[mp*3+0] - dampingCoef*mpts->velocity[mp*3+0] )*dt;
				mpts->velocity[mp*3+1]         = mpts->velocity[mp*3+1] + (mpts->acceleration[mp*3+1] - dampingCoef*mpts->velocity[mp*3+1] )*dt;
				if(mpmParameters->is3D)
				{
					mpts->velocity[mp*3+2]     =mpts->velocity[mp*3+2]+ (mpts->acceleration[mp*3+2] - dampingCoef*mpts->velocity[mp*3+2] )* dt;
				}

					// Update Position
					mpts->position[mp*3+0]         += mpts->Gvelocity[mp*3+0] * dt;
					mpts->position[mp*3+1]         += mpts->Gvelocity[mp*3+1] * dt;
					if(mpmParameters->is3D)
					{
						mpts->position[mp*3+2]     += mpts->Gvelocity[mp*3+2] * dt;
					}
				}

			// Update Displacement
			mpts->displacement[mp*3+0] = mpts->position[mp*3+0] - mpts->pos0[mp*3+0];
			mpts->displacement[mp*3+1] = mpts->position[mp*3+1] - mpts->pos0[mp*3+1];
			mpts->displacement[mp*3+2] = mpts->position[mp*3+2] - mpts->pos0[mp*3+2];

			// Update Temperature
			mpts->temperature[mp]          += mpts->temperatureRate[mp] * dt;
			
			// Update Strain
			mpts->strain[mp*6+0]           += mpts->strainRate[mp*6+0]*dt;//exx
			mpts->strain[mp*6+1]           += mpts->strainRate[mp*6+1]*dt;//eyy
			mpts->strain[mp*6+2]           += mpts->strainRate[mp*6+2]*dt;//ezz
			mpts->strain[mp*6+3]           += mpts->strainRate[mp*6+3]*dt;//exy
			mpts->strain[mp*6+4]           += mpts->strainRate[mp*6+4]*dt;//exz
			mpts->strain[mp*6+5]           += mpts->strainRate[mp*6+5]*dt;//eyz
			// Update stresses using the constitutive law
			int matID = mpts->material[mp];
			bool isRigid = (*materials)[matID].isRigid;
			bool isDamage = (*materials)[matID].isDamage;
			if(!isRigid){ // only update stresses for non-rigid material points
				constitutiveModel(mpmParameters,mpts,dt,mp,materials);
				if(isDamage)
				{
					applyDamage(mpmParameters,mpts,dt,mp,materials);
				}
			}
			
        // if(mp==950)
        // {
        //     printf("%f",mpts->displacement[950*3 + 0]);
        // }
				



		// Periodic Corrections - Particle location
		double Lx, Ly, Lz;
		Lx = nds->L[0]-2*nds->cellsize; Ly = nds->L[1]-2*nds->cellsize; Lz = nds->L[2]-2*nds->cellsize;
		if(mpts->position[mp*3+0] > Lx)
		{
			mpts->position[mp*3+0] =  1.0*nds->cellsize +fmod((mpts->position[mp*3+0] + Lx ),Lx);//1.5*nds->cellsize +
		}
		if(mpts->position[mp*3+0] < 1*nds->cellsize)
		{
			mpts->position[mp*3+0] += Lx-1.0*nds->cellsize ;
		}
		if(mpts->position[mp*3+1] > Ly)
		{
			mpts->position[mp*3+1] = 1.0*nds->cellsize + fmod((mpts->position[mp*3+1] + Ly ),Ly);
		}
		if(mpts->position[mp*3+1] < 1*nds->cellsize)
		{
			mpts->position[mp*3+1] += Ly-1.0*nds->cellsize ;
		}
        if(mpmParameters->is3D)
        {
            if(mpts->position[mp*3+2] > Lz)
            {
                mpts->position[mp*3+2] = 1.0*nds->cellsize + fmod((mpts->position[mp*3+2] + Lz ),Lz);
            }
            if(mpts->position[mp*3+2] < 1*nds->cellsize)
            {
                mpts->position[mp*3+2] += Lz-1.0*nds->cellsize ;
            }
        }		
// Need to set the cell IDs for particles here to enable OpenMP extrapolation to nodes
		int cellid ;
		int idx = ( (mpts->position[mp*3+0]-xmin-0.25*cellsize) / cellsize );
		int idy = ( (mpts->position[mp*3+1]-ymin-0.25*cellsize) / cellsize );
		int idz = ( (mpts->position[mp*3+2]-zmin-0.25*cellsize) / cellsize );
		cellid = idx + idy*(Nx) + idz*(Nx)*(Ny);
		mpts->cellID[mp]=cellid;

};


CPUGPU void constitutiveModel(Parameters *mpmParameters, MaterialPoints *mpts, double dt, int mp, std::vector<MaterialModels>* materials){
	int matID   = mpts->material[mp];
    int matType = (*materials)[matID].model;
	

    switch (matType)
	{
		case 1:
			elastic(mpmParameters,mpts,dt,mp,materials); // in ConstitutiveModels.cpp/hpp
			break;
		case 2:
			newtonianLiquid(mpmParameters,mpts,dt,mp,materials); // in ConstitutiveModels.cpp/hpp
			break;
		case 3:
			perfectPlasticity(mpmParameters,mpts,dt,mp,materials); // in ConstitutiveModels.cpp/hpp
			break;
		case 4:
			druckerPrager_Damage(mpmParameters,mpts,dt,mp,materials); // in ConstitutiveModels.cpp/hpp
			break;
		case 5:
			turbulentLiquid(mpmParameters,mpts,dt,mp,materials); // in ConstitutiveModels.cpp/hpp
			break;
	};
}

CPUGPU void resetMaterialPointRates(MaterialPoints* mpts, int mp){
			//Now all state variables have been updated, reset the rate terms
		mpts->densityRate[mp]       = 0.0;
		mpts->acceleration[mp*3+0]  = 0.0;
		mpts->acceleration[mp*3+1]  = 0.0;
		mpts->acceleration[mp*3+2]  = 0.0;
		mpts->Gvelocity[mp*3+0]     = 0.0;
		mpts->Gvelocity[mp*3+1]     = 0.0;
		mpts->Gvelocity[mp*3+2]     = 0.0;
		mpts->strainRate[mp*6+0]    = 0.0;
		mpts->strainRate[mp*6+1]    = 0.0;
		mpts->strainRate[mp*6+2]    = 0.0;
		mpts->strainRate[mp*6+3]    = 0.0;
		mpts->strainRate[mp*6+4]    = 0.0;
		mpts->strainRate[mp*6+5]    = 0.0;
		mpts->temperatureRate[mp]   = 0.0;
};



CPUGPU void applyDamage(Parameters *mpmParameters,MaterialPoints* mpts,double dt,int mp,std::vector<MaterialModels>* materials){
	int matID                   = mpts->material[mp];
	int damageType = (*materials)[matID].damageModel;
	switch(damageType)
	{
		case 1:
			gradyKippDamage(mpmParameters,mpts,dt,mp,materials); // in DamageModels.cpp/hpp
	}
};



CPUGPU void enforceMatPointBC(MaterialPoints* mpts, int phase,int mp){

		int bc_type = mpts->BCTYPE[mp];
		switch(bc_type)
		{
			case BCTypes::mp_vx:
				bc_mp_vx(mpts,mp,phase);
				break;
			case BCTypes::mp_vy:
				bc_mp_vy(mpts,mp,phase);
				break;
			case BCTypes::mp_vz:
				bc_mp_vz(mpts,mp,phase);
				break;
			case BCTypes::mp_vx_vy:
				bc_mp_vx_vy(mpts,mp,phase);
                break;
            case BCTypes::mp_vxsin:
                bc_mp_vxsin(mpts,mp,phase);
                break;
            case BCTypes::mp_vysin:
                bc_mp_vysin(mpts,mp,phase);
                break;
            case BCTypes::mp_syysin:
                bc_mp_syysin(mpts,mp,phase);
				break;
			case BCTypes::mp_fx:
                bc_mp_fx(mpts,mp,phase);
                break;
			case BCTypes::mp_fy:
                bc_mp_fy(mpts,mp,phase);
				break;
			case BCTypes::mp_fz:
                bc_mp_fz(mpts,mp,phase);
                break;
			case BCTypes::mp_fx_fy:
                bc_mp_fx_fy(mpts,mp,phase);
                break;
			case BCTypes::mp_stress_xx_yy:
				bc_mp_stress_xx_yy(mpts,mp,phase);
				break;
	}
};

CPUGPU void bc_mp_fx_fy(MaterialPoints* mpts, int mp, int phase){
		mpts->acceleration[3*mp+0] += mpts->BCVAL[4*mp+0]/mpts->mass[mp];
		mpts->acceleration[3*mp+1] += mpts->BCVAL[4*mp+1]/mpts->mass[mp];

};

CPUGPU void bc_mp_stress_xx_yy(MaterialPoints* mpts, int mp, int phase){
	mpts->stress[6*mp + 0] = mpts->BCVAL[4*mp+0];
	mpts->stress[6*mp + 1] = mpts->BCVAL[4*mp+1];
	mpts->stresstemp[6*mp + 0] = mpts->BCVAL[4*mp+0];
	mpts->stresstemp[6*mp + 1] = mpts->BCVAL[4*mp+1];

};

CPUGPU void bc_mp_vx(MaterialPoints* mpts, int mp, int phase){
	mpts->velocity[3*mp+0] = mpts->BCVAL[4*mp+0];
};

CPUGPU void bc_mp_vxsin(MaterialPoints* mpts, int mp, int phase){
    mpts->velocity[3*mp+0] = mpts->BCVAL[4*mp+0]*sin(2.0*3.14159*mpts->BCVAL[4*mp+1]*mpts->t);
	
}

CPUGPU void bc_mp_vy(MaterialPoints* mpts, int mp, int phase){
	mpts->velocity[3*mp+0] = 0.0;
	mpts->velocity[3*mp+1] = mpts->BCVAL[4*mp+0];
	mpts->velocity[3*mp+2] = 0.0;

};
CPUGPU void bc_mp_vz(MaterialPoints* mpts, int mp, int phase){
	mpts->velocity[3*mp+0] = 0.0;
	mpts->velocity[3*mp+1] = 0.0;
   	mpts->velocity[3*mp+2] = mpts->BCVAL[4*mp+0];


};
CPUGPU void bc_mp_vysin(MaterialPoints* mpts, int mp, int phase){
    mpts->velocity[3*mp+1] = mpts->BCVAL[4*mp+0]*sin(2.0*3.14159*mpts->BCVAL[4*mp+1]*mpts->t);
	mpts->Gvelocity[3*mp+1] =0.0;
	mpts->acceleration[3*mp+1]=0.0;

}

CPUGPU void bc_mp_syysin(MaterialPoints* mpts, int mp, int phase){
    mpts->stress[6*mp+1] = mpts->BCVAL[4*mp+0]*sin(2.0*3.14159*mpts->BCVAL[4*mp+1]*mpts->t);
	mpts->velocity[3*mp+0] = 0.0;
	mpts->velocity[3*mp+1] = 0.0;
	mpts->velocity[3*mp+2] = 0.0;
	mpts->Gvelocity[3*mp+0] = 0.0;
	mpts->Gvelocity[3*mp+1] = 0.0;
	mpts->Gvelocity[3*mp+2] = 0.0;

}

CPUGPU void bc_mp_vx_vy(MaterialPoints* mpts, int mp, int phase){
	mpts->velocity[3*mp+0] = mpts->BCVAL[4*mp+0];
	mpts->velocity[3*mp+1] = mpts->BCVAL[4*mp+1];
	mpts->velocity[3*mp+2] = 0.0;

};

CPUGPU void bc_mp_fx(MaterialPoints* mpts, int mp, int phase){
    mpts->acceleration[3*mp+0] += mpts->BCVAL[4*mp+0]/mpts->mass[mp];
	mpts->acceleration[3*mp+1] += 0.0;
};

CPUGPU void bc_mp_fy(MaterialPoints* mpts, int mp, int phase){
    mpts->acceleration[3*mp+0] = 0.0;
	mpts->acceleration[3*mp+1] += mpts->BCVAL[4*mp+0]/mpts->mass[mp];
};

CPUGPU void bc_mp_fz(MaterialPoints* mpts, int mp, int phase){
    mpts->acceleration[3*mp+0] += 0.0;
	mpts->acceleration[3*mp+1] += 0.0;
    mpts->acceleration[3*mp+2] += mpts->BCVAL[4*mp+0]/mpts->mass[mp];
};

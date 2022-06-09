#include "Extrapolation.hpp"
#include "MaterialPoints.hpp"
#include "MaterialModels.hpp"
#include "Nodes.hpp"
#include <stdio.h>
#include <cmath>
#include "omp.h"
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <math.h>



class sort_indices
{
   private:
	   managed_vector<int> *mparr;
   public:
     	sort_indices(managed_vector<int>* parr)
		{
			mparr = parr;
	    }

		bool operator()(int i,int j) const
		{
			return mparr->at(i)<mparr->at(j); 
		}
};

CPUGPU void extrapToNodes(Parameters *mpmParameters, MaterialPoints* mpts, Nodes* nds,std::vector<MaterialModels>* materials, int n)
{
/////////////////////////////////////////////////////////////
        //------Step 2. Hash nodes to get the starting cell (based on front lower left cdts)
        /////////////////////////////////////////////////////////////
        int Nx          = nds->Nx-1;
        int Ny          = nds->Ny-1;
        int Nz          = nds->Nz-1;
        double cellsize = nds->cellsize;
        double xmin     = nds->position[0*3 + 0];
        double ymin     = nds->position[0*3 + 1];
        double zmin     = nds->position[0*3 + 2];
        int nidx_x = ( (nds->position[n*3+0]-xmin) / cellsize );
        int nidx_y = ( (nds->position[n*3+1]-ymin) / cellsize );
        int nidx_z = ( (nds->position[n*3+2]-zmin) / cellsize );
        int nodeCellID = nidx_x + nidx_y*Nx + nidx_z*Nx*Ny;
        /////////////////////////////////////////////////////////////
        //-------Step 3. Loop over each cell and use the sorted array to get the particle IDs
        /////////////////////////////////////////////////////////////
        if(mpmParameters->is3D)
        {
            for (int k = -2; k < 2; k++)
            {
                for (int j = -2; j < 2; j++)
                {
                    for (int i = -2; i < 2; i++)
                    {
                        int celli   = nidx_x+i;
                        int cellj   = nidx_y+j;
                        int cellk   = nidx_z+k;
                        int cell_ID = celli+cellj*Nx+cellk*Nx*Ny; // Calculate the cellID the node need to look for matpoints

                        if((cell_ID>=Nx*Ny*Nz) || (cell_ID<0) ) // If the cell is not actually visible to this node (edge detection)
                        {
                            cell_ID = -1;
                        }

                        if(cell_ID>=0)//Skip eroneous cell IDs
                        {
                            // Get the first particle and the number of particles in this cell
                            int p_index  = nds->cells[cell_ID];
                            int npInCell = nds->cell_pcount[cell_ID];

                            if(npInCell>0)
                            {
                                ///////////////////////////////////////////////////////////// 
                                //------- Step 4. Using particle IDs, gather information on Mass, Momenta, Stress, etc.
                                /////////////////////////////////////////////////////////////
                                for(int mp = p_index; mp<npInCell+p_index; mp++)
                                {
                                    // Calculate the weighting function between node n and mat point mp, and write mass, momenta etc. to node. 
                                    extrapolate_Matpnt_to_Node(mpmParameters,mpts,nds,materials,mp,n);// in Extrapolation.c/hpp
                                }; // Use all mat points in current cell and contribute to node
                            };// If this cell has any particles inside it
                        };//checking for non-negative cell IDs
                    };//x4 cells
                };//x4 cells
            };//x4 cells
        }
        else{ // in 2D
            for (int j = -2; j < 2; j++)
            {
                for (int i = -2; i < 2; i++)
                {
                    int celli   = nidx_x+i;
                    int cellj   = nidx_y+j;
                    int cell_ID = celli+cellj*Nx; // Calculate the cellID the node need to look for matpoints

                    if((cell_ID>=Nx*Ny) || (cell_ID<0) ) // If the cell is not actually visible to this node (edge detection)
                    {
                        cell_ID = -1;
                    }

                    if(cell_ID>=0)//Skip eroneous cell IDs
                    {
                        // Get the first particle and the number of particles in this cell
                        int p_index  = nds->cells[cell_ID];
                        int npInCell = nds->cell_pcount[cell_ID];

                        if(npInCell>0)
                        {
                            ///////////////////////////////////////////////////////////// 
                            //------- Step 4. Using particle IDs, gather information on Mass, Momenta, Stress, etc.
                            /////////////////////////////////////////////////////////////
                            for(int mp = p_index; mp<npInCell+p_index; mp++)
                            {
                                // Calculate the weighting function between node n and mat point mp, and write mass, momenta etc. to node. 
                                extrapolate_Matpnt_to_Node(mpmParameters,mpts,nds,materials,mp,n);// in Extrapolation.c/hpp
                            }; // Use all mat points in current cell and contribute to node
                        };// If this cell has any particles inside it
                    };//checking for non-negative cell IDs
                };//x4 cells
            };//x4 cells
		}
	};

CPUGPU void extrapToMatPoints(Parameters *mpmParameters, MaterialPoints* mpts, Nodes* nds,std::vector<MaterialModels>* materials,int mp)
{
	int N2D[9];
	int N3D[27];
	int ndNumbers= 9;
	mpts->temperatureGrad[3*mp+0] = 0.0; // have to zero the temp grad here.
	mpts->temperatureGrad[3*mp+1] = 0.0; // have to zero the temp grad here.
	mpts->temperatureGrad[3*mp+2] = 0.0; // have to zero the temp grad here.
	mpts->normals[mp*3 +0]        = 0.0;
	mpts->normals[mp*3 +1]        = 0.0;
	mpts->normals[mp*3 +2]        = 0.0;
	if (mpmParameters->is3D)
	{
		getGIMPNodes(mpmParameters, mpts, nds, mp, N3D);
		ndNumbers = 27;
	}
	else
	{
		getGIMPNodes(mpmParameters, mpts, nds, mp, N2D);
		ndNumbers = 9;
	}
	
	for(int nd=0;nd<ndNumbers;nd++)
	{
		int n;
		if (mpmParameters->is3D)
		{
			n =N3D[nd];
		}
		else
		{
			 n = N2D[nd];
		}
	 // Calculate the weighting function between node n and mat point mp, and write strainrate, densityrate etc. to material point mp. 		
		extrapolate_Node_to_Matpnt(mpmParameters,mpts,nds,materials,mp,n);// in Extrapolation.c/hpp
	}//Node loop
}

//------------------------------------HELPER FUNCTIONS------------------------------------------------//
CPUGPU void extrapolate_Matpnt_to_Node(Parameters *mpmParameters, MaterialPoints* mpts, Nodes* nds,std::vector<MaterialModels>* materials, int mp, int n){
// Find the shape function and gradient values for the material point/node pair
	double shape_x     = shapeFunctionGIMP(mpts,nds,mp,n,0,0);
	double shape_y     = shapeFunctionGIMP(mpts,nds,mp,n,0,1);
	double shape_z     = shapeFunctionGIMP(mpts,nds,mp,n,0,2);
	double gradShape_x = shapeFunctionGIMP(mpts,nds,mp,n,1,0);
	double gradShape_y = shapeFunctionGIMP(mpts,nds,mp,n,1,1);
	double gradShape_z = shapeFunctionGIMP(mpts,nds,mp,n,1,2);
	double mp_mass     = mpts->mass[mp];	
	int    matID       = mpts->material[mp];
	double mp_density  = mpts->density[mp];
	double mp_temp     = mpts->temperature[mp];
	double mp_tempGradX = mpts->temperatureGrad[3*mp + 0];
	double mp_tempGradY = mpts->temperatureGrad[3*mp + 1];
	double mp_tempGradZ = mpts->temperatureGrad[3*mp + 2];
	int    mpObjID	   = mpts->objID[mp];
	bool   isRigid     = (*materials)[matID].isRigid;
	if(isRigid) mp_mass = 1.0;
	if (!mpmParameters->is3D)
	{
		shape_z = 1;
		gradShape_z = 0.0;
	}
//////////////////////////////// Contact/Multibody ////////////////////////////////////////////////////////
	// Add each body to a separate velocity field first
	if(nds->bodyID[2*n + 0] == -1)
	{
		 nds->bodyID[2*n + 0] = mpObjID;
	}
	if((nds->bodyID[2*n + 1] == -1)&&(nds->bodyID[2*n + 0] != mpObjID))
	{
		nds->bodyID[2*n + 1] = mpObjID;
	}


	if(isRigid && (nds->hasRigid[2*n+0]==false) && (nds->bodyID[2*n + 0] == mpObjID))
	{
		nds->hasRigid[2*n+0]=true;
	}
	else if(isRigid && (nds->hasRigid[2*n+1]==false) && (nds->bodyID[2*n + 1] == mpObjID))
	{
		nds->hasRigid[2*n+1]=true;
	}




	if(nds->bodyID[2*n + 0]  == mpObjID)
	{

		nds->mass_body1[n]           += shape_x*shape_y*shape_z*mp_mass;
		nds->momentum_body1[3*n + 0] += shape_x*shape_y*shape_z*mp_mass*mpts->velocity[mp*3 + 0]; //px
		nds->momentum_body1[3*n + 1] += shape_x*shape_y*shape_z*mp_mass*mpts->velocity[mp*3 + 1]; //py
		nds->momentum_body1[3*n + 2] += shape_x*shape_y*shape_z*mp_mass*mpts->velocity[mp*3 + 2]; //pz
			// Internal forces
		nds->force_body1[n*3 + 0]	 -= (mp_mass/mp_density)*((gradShape_x*shape_y*shape_z*mpts->stress[mp*6+0]) + (gradShape_y*shape_x*shape_z*mpts->stress[mp*6+3]) + (gradShape_z*shape_x*shape_y*mpts->stress[mp*6+5]) );//Fint_x
		nds->force_body1[n*3 + 1] 	 -= (mp_mass/mp_density)*((gradShape_x*shape_y*shape_z*mpts->stress[mp*6+3]) + (gradShape_y*shape_x*shape_z*mpts->stress[mp*6+1]) + (gradShape_z*shape_x*shape_y*mpts->stress[mp*6+4]) );// Fint_y
		nds->force_body1[n*3 + 2]    -= (mp_mass/mp_density)*((gradShape_x*shape_y*shape_z*mpts->stress[mp*6+5]) + (gradShape_y*shape_x*shape_z*mpts->stress[mp*6+4]) + (gradShape_z*shape_x*shape_y*mpts->stress[mp*6+2]) );// Fint_z
if(!isRigid){
        // External forces
		nds->force_body1[n*3 + 0] 	 += shape_x*shape_y*shape_z*mp_mass*mpmParameters->gravity[0];// Fext_x
		nds->force_body1[n*3 + 1] 	 += shape_x*shape_y*shape_z*mp_mass*mpmParameters->gravity[1];// Fext_y
		nds->force_body1[n*3 + 2]    += shape_x*shape_y*shape_z*mp_mass*mpmParameters->gravity[2];// Fext_z
}
					// Calculate normals
		nds->normals_body1[3*n + 0] += (gradShape_x*shape_y*shape_z*mpts->volume[mp]);//;*(shape_x*shape_y*shape_z*mpts->volume[mp])/(nds->cellsize*nds->cellsize);
		nds->normals_body1[3*n + 1] += (shape_x*gradShape_y*shape_z*mpts->volume[mp]);//*(shape_x*shape_y*shape_z*mpts->volume[mp])/(nds->cellsize*nds->cellsize);
		nds->normals_body1[3*n + 2] += (shape_x*shape_y*gradShape_z*mpts->volume[mp]);//*(shape_x*shape_y*shape_z*mpts->volume[mp])/(nds->cellsize*nds->cellsize);
	}
	else if(nds->bodyID[2*n + 1]  == mpObjID)
	{
		nds->mass_body2[n]           += shape_x*shape_y*shape_z*mp_mass;
		nds->momentum_body2[3*n + 0] += shape_x*shape_y*shape_z*mp_mass*mpts->velocity[mp*3 + 0]; //px
		nds->momentum_body2[3*n + 1] += shape_x*shape_y*shape_z*mp_mass*mpts->velocity[mp*3 + 1]; //py
		nds->momentum_body2[3*n + 2] += shape_x*shape_y*shape_z*mp_mass*mpts->velocity[mp*3 + 2]; //pz
		
		nds->force_body2[n*3 + 0]	    -= (mp_mass/mp_density)*((gradShape_x*shape_y*shape_z*mpts->stress[mp*6+0]) + (gradShape_y*shape_x*shape_z*mpts->stress[mp*6+3]) + (gradShape_z*shape_x*shape_y*mpts->stress[mp*6+5]) );//Fint_x
		nds->force_body2[n*3 + 1] 	    -= (mp_mass/mp_density)*((gradShape_x*shape_y*shape_z*mpts->stress[mp*6+3]) + (gradShape_y*shape_x*shape_z*mpts->stress[mp*6+1]) + (gradShape_z*shape_x*shape_y*mpts->stress[mp*6+4]) );// Fint_y
		nds->force_body2[n*3 + 2]       -= (mp_mass/mp_density)*((gradShape_x*shape_y*shape_z*mpts->stress[mp*6+5]) + (gradShape_y*shape_x*shape_z*mpts->stress[mp*6+4]) + (gradShape_z*shape_x*shape_y*mpts->stress[mp*6+2]) );// Fint_z

if(!isRigid){
		// External forces
		nds->force_body2[n*3 + 0] 	 += shape_x*shape_y*shape_z*mp_mass*mpmParameters->gravity[0];// Fext_x
		nds->force_body2[n*3 + 1] 	 += shape_x*shape_y*shape_z*mp_mass*mpmParameters->gravity[1];// Fext_y
        nds->force_body2[n*3 + 2]    += shape_x*shape_y*shape_z*mp_mass*mpmParameters->gravity[2];// Fext_z
}
				// Calculate normals
		nds->normals_body2[3*n + 0] += (gradShape_x*shape_y*shape_z*mpts->volume[mp]);//*(shape_x*shape_y*shape_z*mpts->volume[mp])/(nds->cellsize*nds->cellsize);
		nds->normals_body2[3*n + 1] += (shape_x*gradShape_y*shape_z*mpts->volume[mp]);//*(shape_x*shape_y*shape_z*mpts->volume[mp])/(nds->cellsize*nds->cellsize);
		nds->normals_body2[3*n + 2] += (shape_x*shape_y*gradShape_z*mpts->volume[mp]);//*(shape_x*shape_y*shape_z*mpts->volume[mp])/(nds->cellsize*nds->cellsize);
	}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
	

	// Calculate the center of mass fields
	// mass, momenta and force to the nodes
	
	nds->mass[n]           	    += shape_x*shape_y*shape_z*mp_mass;
	nds->momentum[n*3 + 0] 	    += shape_x*shape_y*shape_z*mp_mass*mpts->velocity[mp*3 + 0];// px
	nds->momentum[n*3 + 1] 	    += shape_x*shape_y*shape_z*mp_mass*mpts->velocity[mp*3 + 1];// py
	nds->momentum[n*3 + 2]      += shape_x*shape_y*shape_z*mp_mass*mpts->velocity[mp*3 + 2];// pz

	// Internal forces
	nds->force[n*3 + 0]	        -= (mp_mass/mp_density)*((gradShape_x*shape_y*shape_z*mpts->stress[mp*6+0]) + (gradShape_y*shape_x*shape_z*mpts->stress[mp*6+3]) + (gradShape_z*shape_x*shape_y*mpts->stress[mp*6+5]) );//Fint_x
	nds->force[n*3 + 1] 	    -= (mp_mass/mp_density)*((gradShape_x*shape_y*shape_z*mpts->stress[mp*6+3]) + (gradShape_y*shape_x*shape_z*mpts->stress[mp*6+1]) + (gradShape_z*shape_x*shape_y*mpts->stress[mp*6+4]) );// Fint_y
	nds->force[n*3 + 2]         -= (mp_mass/mp_density)*((gradShape_x*shape_y*shape_z*mpts->stress[mp*6+5]) + (gradShape_y*shape_x*shape_z*mpts->stress[mp*6+4]) + (gradShape_z*shape_x*shape_y*mpts->stress[mp*6+2]) );// Fint_z

	// External forces
	nds->force[n*3 + 0] 	    += shape_x*shape_y*shape_z*mp_mass*mpmParameters->gravity[0];// Fext_x
	nds->force[n*3 + 1] 	    += shape_x*shape_y*shape_z*mp_mass*mpmParameters->gravity[1];// Fext_y
	nds->force[n*3 + 2]         += shape_x*shape_y*shape_z*mp_mass*mpmParameters->gravity[2];// Fext_z

	// Thermal Properties
	nds->Tmass[n]               += shape_x*shape_y*shape_z*mp_mass*(*materials)[matID].specificHeatCapacity;
	nds->temperature[n]         += shape_x*shape_y*shape_z*mp_temp*mp_mass*(*materials)[matID].specificHeatCapacity;
	nds->temperatureRate[n]     -= (mp_mass/mp_density)*((*materials)[matID].thermalConductivity)*((gradShape_x*shape_y*shape_z*mp_tempGradX)+(shape_x*gradShape_y*shape_z*mp_tempGradY)+(shape_x*shape_y*gradShape_z*mp_tempGradZ)); // Internal thermal
	

};

CPUGPU void extrapolate_Node_to_Matpnt(Parameters *mpmParameters, MaterialPoints* mpts, Nodes* nds,std::vector<MaterialModels>* materials, int mp, int n){
// Check for 0 nodal masses and aCPUGPU void
	if(nds->mass[n]!=0.0)
	{
		// Find the shape function and gradient values for the material point/node pair
		double shape_x     = shapeFunctionGIMP(mpts,nds,mp,n,0,0);
		double shape_y     = shapeFunctionGIMP(mpts,nds,mp,n,0,1);
		double shape_z     = shapeFunctionGIMP(mpts,nds,mp,n,0,2);
		double gradShape_x = shapeFunctionGIMP(mpts,nds,mp,n,1,0);
		double gradShape_y = shapeFunctionGIMP(mpts,nds,mp,n,1,1);
		double gradShape_z = shapeFunctionGIMP(mpts,nds,mp,n,1,2);
		int    matID       = mpts->material[mp];
		double mp_mass     = mpts->mass[mp];	
		double mp_density  = mpts->density[mp];
		int    mpObjID	   = mpts->objID[mp];
		bool   isRigid     = (*materials)[matID].isRigid;
		if (!mpmParameters->is3D)
		{
			shape_z = 1;
			gradShape_z = 0.0;
		}
		// Calculate density rate, acceleration, and strain rate on particles from nodes

		if(nds->bodyID[2*n+0]==mpObjID)
		{ 
			mpts->densityRate[mp]           -= mpts->density[mp]*((gradShape_x*shape_y*shape_z*nds->velocity_body1[n*3+0])+(gradShape_y*shape_z*shape_x*nds->velocity_body1[n*3+1])+(gradShape_z*shape_x*shape_y*nds->velocity_body1[n*3+2]));

			mpts->acceleration[mp*3+0]      += shape_x*shape_y*shape_z*nds->acceleration_body1[n*3 +0];// ax
			mpts->acceleration[mp*3+1]      += shape_x*shape_y*shape_z*nds->acceleration_body1[n*3 +1];// ay
			mpts->acceleration[mp*3+2]      += shape_x*shape_y*shape_z*nds->acceleration_body1[n*3 +2];// az

			mpts->Gvelocity[mp*3+0]         += shape_x*shape_y*shape_z*nds->velocity_body1[n*3 +0];// vx
			mpts->Gvelocity[mp*3+1]         += shape_x*shape_y*shape_z*nds->velocity_body1[n*3 +1];// vy
			mpts->Gvelocity[mp*3+2]         += shape_x*shape_y*shape_z*nds->velocity_body1[n*3 +2];// vz

			mpts->strainRate[mp*6 +0]       += (gradShape_x*shape_y*shape_z*nds->velocity_body1[n*3+0]);//exx
			mpts->strainRate[mp*6 +1]       += (gradShape_y*shape_x*shape_z*nds->velocity_body1[n*3+1]);//eyy
			mpts->strainRate[mp*6 +2]       += (gradShape_z*shape_x*shape_y*nds->velocity_body1[n*3+2]);//ezz
			mpts->strainRate[mp*6 +3]       += 0.5*(gradShape_x*shape_y*shape_z*nds->velocity_body1[n*3+1] + gradShape_y*shape_x*shape_z*nds->velocity_body1[n*3+0]);//exy
			mpts->strainRate[mp*6 +4]       += 0.5*(gradShape_y*shape_x*shape_z*nds->velocity_body1[n*3+2] + gradShape_z*shape_x*shape_y*nds->velocity_body1[n*3+1]);//eyz
			mpts->strainRate[mp*6 +5]       += 0.5*(gradShape_x*shape_y*shape_z*nds->velocity_body1[n*3+2] + gradShape_z*shape_x*shape_y*nds->velocity_body1[n*3+0]);//ezx
							// Contact Mechanics
			mpts->normals[3*mp + 0] += shape_x*shape_y*shape_z*nds->normals_body1[3*n + 0];
			mpts->normals[3*mp + 1] += shape_x*shape_y*shape_z*nds->normals_body1[3*n + 1];
			mpts->normals[3*mp + 2] += shape_x*shape_y*shape_z*nds->normals_body1[3*n + 2];
		}
		else if(nds->bodyID[2*n+1]==mpObjID)
		{
			mpts->densityRate[mp]           -= mpts->density[mp]*((gradShape_x*shape_y*shape_z*nds->velocity_body2[n*3+0])+(gradShape_y*shape_z*shape_x*nds->velocity_body2[n*3+1])+(gradShape_z*shape_x*shape_y*nds->velocity_body2[n*3+2]));
			mpts->acceleration[mp*3+0]      += shape_x*shape_y*shape_z*nds->acceleration_body2[n*3 +0];// ax
			mpts->acceleration[mp*3+1]      += shape_x*shape_y*shape_z*nds->acceleration_body2[n*3 +1];// ay
			mpts->acceleration[mp*3+2]      += shape_x*shape_y*shape_z*nds->acceleration_body2[n*3 +2];// az

			mpts->Gvelocity[mp*3+0]         += shape_x*shape_y*shape_z*nds->velocity_body2[n*3 +0];// vx
			mpts->Gvelocity[mp*3+1]         += shape_x*shape_y*shape_z*nds->velocity_body2[n*3 +1];// vy
			mpts->Gvelocity[mp*3+2]         += shape_x*shape_y*shape_z*nds->velocity_body2[n*3 +2];// vz

			mpts->strainRate[mp*6 +0]       += (gradShape_x*shape_y*shape_z*nds->velocity_body2[n*3+0]);//exx
			mpts->strainRate[mp*6 +1]       += (gradShape_y*shape_x*shape_z*nds->velocity_body2[n*3+1]);//eyy
			mpts->strainRate[mp*6 +2]       += (gradShape_z*shape_x*shape_y*nds->velocity_body2[n*3+2]);//ezz
			mpts->strainRate[mp*6 +3]       += 0.5*(gradShape_x*shape_y*shape_z*nds->velocity_body2[n*3+1] + gradShape_y*shape_x*shape_z*nds->velocity_body2[n*3+0]);//exy
			mpts->strainRate[mp*6 +4]       += 0.5*(gradShape_y*shape_x*shape_z*nds->velocity_body2[n*3+2] + gradShape_z*shape_x*shape_y*nds->velocity_body2[n*3+1]);//eyz
			mpts->strainRate[mp*6 +5]       += 0.5*(gradShape_x*shape_y*shape_z*nds->velocity_body2[n*3+2] + gradShape_z*shape_x*shape_y*nds->velocity_body2[n*3+0]);//ezx
							// Contact Mechanics
			mpts->normals[3*mp + 0] += shape_x*shape_y*shape_z*nds->normals_body2[3*n + 0];
			mpts->normals[3*mp + 1] += shape_x*shape_y*shape_z*nds->normals_body2[3*n + 1];
			mpts->normals[3*mp + 2] += shape_x*shape_y*shape_z*nds->normals_body2[3*n + 2];
		}
		

		
		mpts->temperatureGrad[3*mp + 0] += (gradShape_x*shape_y*shape_z)*nds->temperature[n];
		mpts->temperatureGrad[3*mp + 1] += (shape_x*gradShape_y*shape_z)*nds->temperature[n];
		mpts->temperatureGrad[3*mp + 2] += (shape_x*shape_y*gradShape_z)*nds->temperature[n];
		mpts->temperatureRate[mp]       += shape_x*shape_y*shape_z*nds->temperatureRate[n]; 
			


	}
};


CPUGPU void sortParticlesByCellID_CPU(Parameters* mpmParameters,MaterialPoints* mpts, Nodes* nds){
	
	int NcellsX = nds->Nx-1;
	int NcellsY = nds->Ny-1;
	int NcellsZ = nds->Nz-1;
	if(!mpmParameters->is3D)
	{
		NcellsZ = 1;
	}
	#pragma omp parallel for
	for(int g = 0; g < (NcellsX* NcellsY * NcellsZ); g++)
	{
		nds->cell_pcount[g]=0;
	}

	#pragma omp parallel for
	for(int p = 0; p < mpts->numParticles; p++) 
	{	
		mpts->pidx[p]=p;
	}
	// must be done in serial to aCPUGPU void data race
	for(int p = 0; p < mpts->numParticles; p++) 
	{	
		int gid = mpts->cellID[p];

		if(gid >= 0) 
		{
			nds->cell_pcount[gid]++;
		}
	}
	
	// sort particle index by cell ID to have particles grouped by cell
	std::sort(mpts->pidx.begin(), mpts->pidx.end(), sort_indices(&(mpts->cellID)) );
    // Need to sort actual particles
    // Copy particles to their sorted positions
	
    #pragma omp parallel for
    for (int i = 0; i < mpts->numParticles; i++) 
	{
		int idx_prev = mpts->pidx[i];
		mpts->postemp[i*3+0]     = mpts->position[idx_prev*3+0];
		mpts->postemp[i*3+1]     = mpts->position[idx_prev*3+1];
		mpts->postemp[i*3+2]     = mpts->position[idx_prev*3+2];
		mpts->pos0temp[i*3+0]    = mpts->pos0[idx_prev*3+0];
		mpts->pos0temp[i*3+1]    = mpts->pos0[idx_prev*3+1];
		mpts->pos0temp[i*3+2]    = mpts->pos0[idx_prev*3+2];
		mpts->disptemp[i*3+0]    = mpts->displacement[idx_prev*3+0];
		mpts->disptemp[i*3+1]    = mpts->displacement[idx_prev*3+1];
		mpts->disptemp[i*3+2]    = mpts->displacement[idx_prev*3+2];
		mpts->veltemp[i*3+0]     = mpts->velocity[idx_prev*3+0];
		mpts->veltemp[i*3+1]     = mpts->velocity[idx_prev*3+1];
		mpts->veltemp[i*3+2]     = mpts->velocity[idx_prev*3+2];
		mpts->gveltemp[i*3+0]    = mpts->Gvelocity[idx_prev*3+0];
		mpts->gveltemp[i*3+1]    = mpts->Gvelocity[idx_prev*3+1];
		mpts->gveltemp[i*3+2]    = mpts->Gvelocity[idx_prev*3+2];
		mpts->stresstemp[i*6+0]  = mpts->stress[idx_prev*6+0];
		mpts->stresstemp[i*6+1]  = mpts->stress[idx_prev*6+1];
		mpts->stresstemp[i*6+2]  = mpts->stress[idx_prev*6+2];
		mpts->stresstemp[i*6+3]  = mpts->stress[idx_prev*6+3];
		mpts->stresstemp[i*6+4]  = mpts->stress[idx_prev*6+4];
		mpts->stresstemp[i*6+5]  = mpts->stress[idx_prev*6+5];
		mpts->straintemp[i*6+0]  = mpts->strain[idx_prev*6+0];
		mpts->straintemp[i*6+1]  = mpts->strain[idx_prev*6+1];
		mpts->straintemp[i*6+2]  = mpts->strain[idx_prev*6+2];
		mpts->straintemp[i*6+3]  = mpts->strain[idx_prev*6+3];
		mpts->straintemp[i*6+4]  = mpts->strain[idx_prev*6+4];
		mpts->straintemp[i*6+5]  = mpts->strain[idx_prev*6+5];
		mpts->acc_plastic_strain_temp[i] = mpts->acc_plastic_strain[idx_prev];
		mpts->plastictemp[i*6+0] = mpts->plastic_strain[idx_prev*6+0];
		mpts->plastictemp[i*6+1] = mpts->plastic_strain[idx_prev*6+1];
		mpts->plastictemp[i*6+2] = mpts->plastic_strain[idx_prev*6+2];
		mpts->plastictemp[i*6+3] = mpts->plastic_strain[idx_prev*6+3];
		mpts->plastictemp[i*6+4] = mpts->plastic_strain[idx_prev*6+4];
		mpts->plastictemp[i*6+5] = mpts->plastic_strain[idx_prev*6+5];
		mpts->damagetemp[i]      = mpts->damage[idx_prev];
		mpts->smaxtemp[i]        = mpts->smax[idx_prev];
		mpts->tempGradtemp[i]    = mpts->temperatureGrad[idx_prev];
		mpts->temptemp[i]		 = mpts->temperature[idx_prev];
		mpts->masstemp[i]        = mpts->mass[idx_prev];
		mpts->densitytemp[i]     = mpts->density[idx_prev];
		mpts->volumetemp[i]      = mpts->volume[idx_prev];
		mpts->normalstemp[i*3+0] = mpts->normals[idx_prev*3+0];
		mpts->normalstemp[i*3+1] = mpts->normals[idx_prev*3+1];
		mpts->normalstemp[i*3+2] = mpts->normals[idx_prev*3+2];
		mpts->cellIDtemp[i]      = mpts->cellID[idx_prev];
		mpts->materialtemp[i]    = mpts->material[idx_prev];
		mpts->objIDtemp[i]       = mpts->objID[idx_prev];
		mpts->holdtemp[i]        = mpts->hold[idx_prev];
		mpts->BCVALtemp[i*4+0]   = mpts->BCVAL[idx_prev*4+0];
		mpts->BCVALtemp[i*4+1]   = mpts->BCVAL[idx_prev*4+1];
		mpts->BCVALtemp[i*4+2]   = mpts->BCVAL[idx_prev*4+2];
		mpts->BCVALtemp[i*4+3]   = mpts->BCVAL[idx_prev*4+3];
		mpts->BCTYPEtemp[i]      = mpts->BCTYPE[idx_prev];
    }

	// Now that the temp arrays have the sorted variables, flip points so that the main arrays are updated
		// After the sorting, flip the old and temporary arrays

		mpts->position.swap(mpts->postemp);
		mpts->pos0.swap(mpts->pos0temp);
		mpts->displacement.swap(mpts->disptemp);
		mpts->velocity.swap(mpts->veltemp);
		mpts->Gvelocity.swap(mpts->gveltemp);
		mpts->stress.swap(mpts->stresstemp);
		mpts->strain.swap(mpts->straintemp);
		mpts->plastic_strain.swap(mpts->plastictemp);
		mpts->acc_plastic_strain.swap(mpts->acc_plastic_strain_temp);
		mpts->damage.swap(mpts->damagetemp);
		mpts->smax.swap(mpts->smaxtemp);
		mpts->mass.swap(mpts->masstemp);
		mpts->density.swap(mpts->densitytemp);
		mpts->objID.swap(mpts->objIDtemp);
		mpts->volume.swap(mpts->volumetemp);
		mpts->normals.swap(mpts->normalstemp);
		mpts->temperature.swap(mpts->temptemp);
		mpts->temperatureGrad.swap(mpts->tempGradtemp);
		mpts->temperatureRate.swap(mpts->tempRatetemp);
		mpts->BCVAL.swap(mpts->BCVALtemp);
		mpts->material.swap(mpts->materialtemp);
		mpts->cellID.swap(mpts->cellIDtemp);
		mpts->hold.swap(mpts->holdtemp);
		mpts->BCTYPE.swap(mpts->BCTYPEtemp);
		
	// Now need to find the starting particle index for each cell
	int currStart = 0;
	//int cells[NcellsX* NcellsY * NcellsZ];// Index of cells with starting pidx
	nds->cells[0] = 0;
	for (int i = 1; i < NcellsX* NcellsY * NcellsZ; i++) 
	{
		currStart     += nds->cell_pcount[i-1];
		nds->cells[i] = currStart;
	}
	
};

CPUGPU void sortParticlesByCellID_GPU(Parameters* mpmParameters,MaterialPoints* mpts, Nodes* nds){

};


CPUGPU double shapeFunctionClassic(MaterialPoints *mpts,Nodes* nds ,int mpIdx, int ndIdx, int N_or_G, int xyz){
// matpoints - Material Points array holding all the persistant qualities
// nodes - Node array holding all properties relvant for nodes
// mpIdx - Index for current material point
// ndIdx - Index for current node
// N_or_G - Switch to either calculashapeFunctionClassicte Shape of Gradient function value

//Step 1 - get x,y,z positions of material point and node
	double mat_pos = mpts->position[mpIdx*3 +xyz];
	double node_pos = nds->position[ndIdx*3 +xyz];
	double cellsize = nds->cellsize;
	double dist = node_pos-mat_pos;
	double G=0.;
	double S=0.;

	if(N_or_G==0) // Calculate Shape Function
	{
		if(fabs(dist)>cellsize)
			{
				S=0.0;
			}
		else if(dist<=cellsize)
			{
				S= 1.- (fabs(dist)/cellsize);
			}
		return S;
	}
	else if(N_or_G==1)//Calculate Gradient Function
	{
		if(fabs(dist)>cellsize)
			{
				G=0.0;
			}
		else if(dist<=0.0)
			{
				G= 1./cellsize;
			}
		else if(dist>0.0)
			{
				G= -1./cellsize;
			}
		return G;
	}
	return G;
};

CPUGPU double shapeFunctionGIMP(MaterialPoints *mpts,Nodes* nds ,int mpIdx, int ndIdx, int N_or_G, int xyz){
// matpoints - Material Points array holding all the persistant qualities
// nodes - Node array holding all properties relvant for nodes
// mpIdx - Index for current material point
// ndIdx - Index for current node
// N_or_G - Switch to either calculate Shape of Gradient function value

//Step 1 - get x,y,z positions of material point and node
double mat_pos  = mpts->position[mpIdx*3 +xyz];
double node_pos = nds->position[ndIdx*3 +xyz];
double cellsize = nds->cellsize;
double lp       = 0.5*mpts->psep; // 0.5*sqrt(mpts->volume[mpIdx]);//TO DO need to fix this
double dist     = node_pos-mat_pos;
double G        = 0.0;
double S        = 0.0;

if(N_or_G==0) // Calculate Shape Function
{
	if (dist <= -cellsize-lp)
	{
	    S=0.;
	 }
	else if((-cellsize-lp < dist)&&(dist<= -cellsize+lp))
	{
	    S = pow((cellsize + lp + dist),2.0)/(4.*cellsize*lp);
	}
	   
	else if((-cellsize+lp<dist)&&(dist<=-lp))
	{
	    S = 1. + (dist/cellsize);
	}
	   
	else if((-lp<dist)&&(dist<=lp))
	{
	    S= 1. - (pow(dist,2.0) + pow(lp,2.0))/(2.*cellsize*lp);
	}
	
	else if((lp<dist)&&(dist<=cellsize-lp))
	{
	    S= 1. - (dist/cellsize);
	}
	 
	else if((cellsize-lp<dist)&&(dist<=cellsize+lp))
	{
	    S= pow((cellsize + lp - dist),2.0) / (4.*cellsize*lp);
	}
	 
	else if(cellsize+lp>dist)
	{
	    S=0.;
	}
	return S;    

}
else if(N_or_G==1)//Calculate Gradient Function
{

	if (dist < -cellsize-lp)
	  {
	    G =0.0;
	  }
	else if((-cellsize-lp < dist)&&(dist< -cellsize+lp))
	  {
	    G= (cellsize + lp + dist)/(2.*cellsize*lp);
	  }
	else if((-cellsize+lp<=dist)&&(dist<-lp))
	   {
	    G= 1./cellsize;
	   }
	else if((-lp<=dist)&&(dist<lp))
	  {
	    G= -dist/(cellsize*lp);
	  }
	else if((lp<=dist)&&(dist<cellsize-lp))
	  {
	    G= -1./cellsize;
	  }
	else if((cellsize-lp<=dist)&&(dist<=cellsize+lp))
	  {
	    G= -(cellsize+lp-dist)/(2.*cellsize*lp);
	  }
	else if(cellsize+lp>dist)
	 { 
	   	G=0.;
	 }

	G = G*mat_pos/(fabs(mat_pos));

	return -G;
	};
return G;
};

CPUGPU void getGIMPNodes(Parameters* mpmParameters, MaterialPoints *mpts, Nodes *nds, int mpIdx,int *N){

	if(mpmParameters->is3D)
	{
	// Get all 27 Nodes surrounding the material points in 3d
	int Nx = nds->Nx;
	int Ny = nds->Ny;
	int Nz = nds->Nz;
	// For each material point, the 8 nodes that it interacts with
	int index_x = (mpts->position[mpIdx*3 + 0]-(mpts->psep/2) - nds->position[0])/nds->cellsize;
	int index_y = (mpts->position[mpIdx*3 + 1]-(mpts->psep/2) - nds->position[1])/nds->cellsize;
	int index_z = (mpts->position[mpIdx*3 + 2]-(mpts->psep/2) - nds->position[2])/nds->cellsize;
	int nCount = 0;

	for (int k = 0; k < 3; k++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int i = 0; i < 3; i++)
			{
				int ni    = index_x+i;
				int nj    = index_y+j;
				int nk    = index_z+k;
				N[nCount] = ni+nj*Nx+nk*Nx*Ny;

				if( (N[nCount]>=(Nx*Ny*Nz)) || (N[nCount]<0))
				{
					N[nCount] = 0;
				}
				nCount++;
			}
		}
	}
}
else
{
	// Get all 9 Nodes surrounding the material points in 2d
	int Nx = nds->Nx;
	int Ny = nds->Ny;
	// For each material point, the 8 nodes that it interacts with
	int index_x = (mpts->position[mpIdx*3 + 0]-(mpts->psep/2) - nds->position[0])/nds->cellsize;
	int index_y = (mpts->position[mpIdx*3 + 1]-(mpts->psep/2) - nds->position[1])/nds->cellsize;
	int nCount = 0;
		for (int j = 0; j < 3; j++)
		{
			for (int i = 0; i < 3; i++)
			{
				int ni    = index_x+i;
				int nj    = index_y+j;
				N[nCount] = ni+nj*Nx;

				if( (N[nCount]>=(Nx*Ny)) || (N[nCount]<0))
				{
					N[nCount] = 0;
				}
				nCount++;
			}
		}
	}

};

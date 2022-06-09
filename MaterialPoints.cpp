// Class file for MaterialPoints
#include"MaterialPoints.hpp"
#include"Nodes.hpp"
#include <cmath>
#include <stdio.h>
#include <math.h>
#include <string.h>


MaterialPoints::MaterialPoints(Nodes* nodes) // constructor
{
	this->nodes  = nodes;
	int npMAX    = 8*nodes->numNodes; // Max case for material points is to assume entire domain is filled with 8 particles per cell
    N            = new int[27];
    t            = 0;
	numParticles = 0;
};
MaterialPoints::~MaterialPoints(){
	//delete position;
};



void MaterialPoints::addPyckGeom(double* pos, int numParticles){
	this->numParticles = (int)numParticles;
	for(int i=0; i<this->numParticles;i++)
	{
		this->add(i,pos[i*3+0],pos[i*3+1],pos[i*3+2]);
	}
};

void MaterialPoints::addPyckIntField(IntField* intField){
	if(intField->name.compare("material")==0)
	{
		#pragma omp parallel for
		for(int i=0; i<this->numParticles;i++)
		{
			this->material[i] = intField->data[i];
			this->materialtemp[i] = intField->data[i];
		}
		
	}
	else if(intField->name.compare("objectID")==0)
	{
		#pragma omp parallel for
		for(int i=0; i<this->numParticles;i++)
		{
			this->objID[i]     = intField->data[i];
			this->objIDtemp[i] = intField->data[i];
		}
		
	}

};
void MaterialPoints::addPyckDoubleField(DoubleField* doubleField){

};

// take an index and a value to set particle materialID 
void MaterialPoints::setMaterialByCell(int idx, int materialID){
	#pragma omp parallel for
	for(int i=0; i<numParticles; i++)
	{
		int cellID_val = cellID[i];
		if (cellID_val == idx)
		{
			// Set the material
			this->material[i] = materialID;
			this->materialtemp[i] = materialID;
		}
	}


};

void MaterialPoints::setMaterialByCell(std::vector<int> materialArray){
	#pragma omp parallel for
	for(int i=0; i<numParticles; i++)
	{
		int cellID_val = cellID[i];
		int material_ID = materialArray[cellID_val];
		this->material[i] = material_ID;
		this->materialtemp[i] = material_ID;
	}
};


void MaterialPoints::setBC(int p, BCTypes::function type, std::vector<double> val )
{
	this->BCTYPE[p]    = type;
	for(int i=0; i < val.size(); i++)
		this->BCVAL[4*p+i] = val[i];
	};

void MaterialPoints::add(int i, double x, double y,double z){

	position.push_back(x);
	postemp.push_back(x);
    position.push_back(y);
	postemp.push_back(y);
    position.push_back(z);
	postemp.push_back(z);
	pos0.push_back(x);
	pos0temp.push_back(x);
	pos0.push_back(y);
	pos0temp.push_back(y);
	pos0.push_back(z);
	pos0temp.push_back(z);
	displacement.push_back(0.0);
	displacement.push_back(0.0);
	displacement.push_back(0.0);
	disptemp.push_back(0.0);
	disptemp.push_back(0.0);
	disptemp.push_back(0.0);
    velocity.push_back(0.0);
	veltemp.push_back(0.0);
    velocity.push_back(0.0);
	veltemp.push_back(0.0);
    velocity.push_back(0.0);
	veltemp.push_back(0.0);
	Gvelocity.push_back(0.0);
	gveltemp.push_back(0.0);
	Gvelocity.push_back(0.0);
	gveltemp.push_back(0.0);
	Gvelocity.push_back(0.0);
	gveltemp.push_back(0.0);
	stress.push_back(0.0);
	stresstemp.push_back(0.0);
    stress.push_back(0.0);
	stresstemp.push_back(0.0);
    stress.push_back(0.0);
	stresstemp.push_back(0.0);
    stress.push_back(0.0);
	stresstemp.push_back(0.0);
    stress.push_back(0.0);
	stresstemp.push_back(0.0);
    stress.push_back(0.0);
	stresstemp.push_back(0.0);
	strain.push_back(0.0);
	straintemp.push_back(0.0);
    strain.push_back(0.0);
	straintemp.push_back(0.0);
    strain.push_back(0.0);
	straintemp.push_back(0.0);
    strain.push_back(0.0);
	straintemp.push_back(0.0);
    strain.push_back(0.0);
	straintemp.push_back(0.0);
    strain.push_back(0.0);
	straintemp.push_back(0.0);
	acc_plastic_strain.push_back(0.0);
	acc_plastic_strain_temp.push_back(0.0);
	plastic_strain.push_back(0.0);
	plastictemp.push_back(0.0);
    plastic_strain.push_back(0.0);
	plastictemp.push_back(0.0);
    plastic_strain.push_back(0.0);
	plastictemp.push_back(0.0);
    plastic_strain.push_back(0.0);
	plastictemp.push_back(0.0);
    plastic_strain.push_back(0.0);
	plastictemp.push_back(0.0);
    plastic_strain.push_back(0.0);
	plastictemp.push_back(0.0);
	damage.push_back(0.0);
	damagetemp.push_back(0.0);
	smax.push_back(0.0);
	smaxtemp.push_back(0.0);
	material.push_back(0);
	materialtemp.push_back(0);
	objID.push_back(0);
	objIDtemp.push_back(0);
	acceleration.push_back(0.0);
	acceleration.push_back(0.0);
	acceleration.push_back(0.0);
	Gvelocity.push_back(0.0);
	Gvelocity.push_back(0.0);
	Gvelocity.push_back(0.0);
	strainRate.push_back(0.0);
	strainRate.push_back(0.0);
	strainRate.push_back(0.0);
	strainRate.push_back(0.0);
	strainRate.push_back(0.0);
	strainRate.push_back(0.0);
	densityRate.push_back(0.0);
	temperature.push_back(0.0);
	temptemp.push_back(0.0);
	temperatureRate.push_back(0.0);
	tempRatetemp.push_back(0.0);
	temperatureGrad.push_back(0.0);
	temperatureGrad.push_back(0.0);
	temperatureGrad.push_back(0.0);
	tempGradtemp.push_back(0.0);
	tempGradtemp.push_back(0.0);
	tempGradtemp.push_back(0.0);
	density.push_back(0.0);
	densitytemp.push_back(0.0);
	mass.push_back(0.0);
	masstemp.push_back(0.0);
	volume.push_back(0.0);
	volumetemp.push_back(0.0);
	normals.push_back(0.0);
	normalstemp.push_back(0.0);
	normals.push_back(0.0);
	normalstemp.push_back(0.0);
	normals.push_back(0.0);
	normalstemp.push_back(0.0);
	hold.push_back(false);
	holdtemp.push_back(false);
	BCTYPE.push_back(BCTypes::none);
	BCTYPEtemp.push_back(BCTypes::none);
	BCVAL.push_back(0.0);
	BCVALtemp.push_back(0.0);
	BCVAL.push_back(0.0);
	BCVALtemp.push_back(0.0);
	BCVAL.push_back(0.0);
	BCVALtemp.push_back(0.0);
	BCVAL.push_back(0.0);
	BCVALtemp.push_back(0.0);
	int Nx          = this->nodes->Nx;
	int Ny          = this->nodes->Ny;
	int Nz          = this->nodes->Nz;
	double cellsize = this->nodes->cellsize;
	double xmin     = this->nodes->position[0*3 + 0];
	double ymin     = this->nodes->position[0*3 + 1];
	double zmin     = this->nodes->position[0*3 + 2];
	int cellIDval   = 0;
	int idx         = ( (position[i*3+0]-xmin) / cellsize );
	int idy         = ( (position[i*3+1]-ymin) / cellsize );
	int idz         = ( (position[i*3+2]-zmin) / cellsize );
	cellIDval       = idx + idy*(Nx-1) + idz*(Nx-1)*(Ny-1);
	cellID.push_back(cellIDval);
	cellIDtemp.push_back(cellIDval);
	pidx.push_back(0);
};

// Getters
double  MaterialPoints::getPosX(int p){return position[p*3 + 0];};
double  MaterialPoints::getPosY(int p){return position[p*3 + 1];};
double  MaterialPoints::getPosZ(int p){return position[p*3 + 2];};
double  MaterialPoints::getDispX(int p){return displacement[p*3 + 0];};
double  MaterialPoints::getDispY(int p){return displacement[p*3 + 1];};
double  MaterialPoints::getDispZ(int p){return displacement[p*3 + 2];};
double  MaterialPoints::getVelX(int p){return velocity[p*3 + 0];};
double  MaterialPoints::getVelY(int p){return velocity[p*3 + 1];};
double  MaterialPoints::getVelZ(int p){return velocity[p*3 + 2];};
double  MaterialPoints::getStressXX(int p){return stress[p*3 + 0];};
double  MaterialPoints::getStressYY(int p){return stress[p*3 + 1];};
double  MaterialPoints::getStressZZ(int p){return stress[p*3 + 2];};
double  MaterialPoints::getStressXY(int p){return stress[p*3 + 3];};
double  MaterialPoints::getStressYZ(int p){return stress[p*3 + 4];};
double  MaterialPoints::getStressZX(int p){return stress[p*3 + 5];};
int     MaterialPoints::getObjID(int p){return objID[p];};
int		MaterialPoints::getCellID(int p){return cellID[p];};
int     MaterialPoints::getMatID(int p){return material[p];};
double* MaterialPoints::getPosArr(){return this->position.data();};    // Casting here since the managed allocator and std::vectors had different VTPWrite behaviour
double* MaterialPoints::getDispArr(){return this->displacement.data();};
double* MaterialPoints::getMassArr(){return  this->mass.data();};
int*    MaterialPoints::getMatArr(){return this->material.data();};
int*    MaterialPoints::getCellIDArr(){return this->cellID.data();};
double* MaterialPoints::getVelArr(){return  this->velocity.data();};
double* MaterialPoints::getAccArr(){return this->acceleration.data();};
double* MaterialPoints::getStressArr(){return  this->stress.data();};
double* MaterialPoints::getStrainArr(){return this->strain.data();};
double* MaterialPoints::getTemperatureArr(){return  this->temperature.data();};
double* MaterialPoints::getPlasticStrainArr(){return  this->plastic_strain.data();};
double* MaterialPoints::getAccPlasticStrainArr(){return this->acc_plastic_strain.data();};
double* MaterialPoints::getDamageArr(){return this->damage.data();};
double* MaterialPoints::getSMaxArr(){return this->smax.data();};
double* MaterialPoints::getNormalsArr(){return this->normals.data();};
int*    MaterialPoints::getObjectIDArr(){return this->objID.data();};
double* MaterialPoints::getVolumeArr(){return this->volume.data();};
double* MaterialPoints::getDensityArr(){return this->density.data();};

void MaterialPoints::setRotatingBodyVel(int objectID, double centerX, double centerY, double centerZ, double omega){
	#pragma omp parallel for
	for(int p = 0; p < numParticles; p++)
	{
	        if(objID[p]==objectID)
			{
	            double x        = position[p*3 +0];
	            double y        = position[p*3 +1];
	            double r        = sqrt(pow(x-centerX,2.0) + pow(y-centerY,2.0));
	            double costheta = (y-centerY)/r;
	            double sintheta = (x-centerX)/r;
	            double vx       = -omega*r*costheta;
	            double vy       = omega*r*sintheta;
	            velocity[p*3+0] = vx;
		        velocity[p*3+1] = vy;
	            velocity[p*3+2] = 0.0;
				stress[p*6+0]   = 0.0;
				stress[p*6+1]   = 0.0;
				stress[p*6+2]   = 0.0;
				stress[p*6+3]   = 0.0;
				stress[p*6+4]   = 0.0;
				stress[p*6+5]   = 0.0;
		}
	}
};


void MaterialPoints::setVel(int p, double vx, double vy, double vz){
	this->velocity[p*3+0] = vx;
    this->velocity[p*3+1] = vy;
    this->velocity[p*3+2] = vz;
}

void MaterialPoints::setDamage(int p,double D){
	this->damage[p] = D;
}

void MaterialPoints::setStress(int p, double sxx,double syy,double szz,double sxy,double syz,double sxz){
	stress[p*6+0] = sxx;
    stress[p*6+1] = syy;
    stress[p*6+2] = szz;
	stress[p*6+3] = sxy;
	stress[p*6+4] = syz;
	stress[p*6+5] = sxz;
}
void MaterialPoints::setTemp(int p, double T){
	temperature[p] = T;
}

void MaterialPoints::setObjectID(int p, int objid)
{
	this->objID[p] = objid;
	this->objIDtemp[p] = objid;
}

void MaterialPoints::setMaterialID(int p, int materialID)
{
	this->material[p]=materialID;
	this->materialtemp[p]=materialID;
}

void MaterialPoints::holdParticle(int p)
{
	this->hold[p]=true;
	this->holdtemp[p]=true;
}
void MaterialPoints::Serialise(int step, char * filename){

FILE *fp_mpts;
    char filename_mpts[50]; 
    sprintf(filename_mpts,"output/%s_mpts_%d.csv",filename,step);
    fp_mpts=fopen(filename_mpts,"w+");
    double vm,P,T,x,y,z,vx,vy,vz,rho;
    fprintf(fp_mpts, "x,y,z,vm,P,T,vx,vy,vz,rho\n");
    for(int i=0; i<this->numParticles;i++)
    {
        vm  = pow(pow(this->stress[6*i+0],2) + pow(this->stress[6*i+1],2) + pow(this->stress[6*i+2],2) + pow(2.0*this->stress[6*i+3],2) +pow(2.0*this->stress[6*i+4],2) + pow(2.0*this->stress[6*i+5],2),0.5);
        x   = this->position[3*i+0];
        y   = this->position[3*i+1];
        z   = this->position[3*i+2];
        vx  = this->velocity[3*i+0];
        vy  = this->velocity[3*i+1];
        vz  = this->velocity[3*i+2];
	    rho = this->density[i];
		P   = 0.33*(this->stress[6*i+0] + this->stress[6*i+1] + this->stress[6*i+2]);
		T   = this->temperature[i];
        fprintf(fp_mpts, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",x,y,z,vm,P,T,vx,vy,vz,rho);
    }
    fclose(fp_mpts);
}

// Class file for Nodes
#include"Nodes.hpp"
#include "Parameters.hpp"
#include <cmath>
#include <stdio.h>
#include <math.h>
#include <string.h>
Nodes::Nodes(std::vector<double> L, double cellsize) // constructor
{
    this->L          = L;
    this->cellsize   = cellsize;
    this->numNodes   = round((1+(L[0]/cellsize))*(1+(L[1]/cellsize))*(1+(L[2]/cellsize)));
    this->numCells   = (L[0]/cellsize)*(L[1]/cellsize)*(L[2]/cellsize);
    this->t          = 0.0;
	CPUGPUMALLOC(position,			 double, 3 * numNodes);
	CPUGPUMALLOC(velocity,			 double, 3 * numNodes);
	CPUGPUMALLOC(velocity_body1,	 double, 3 * numNodes);
	CPUGPUMALLOC(velocity_body2,	 double, 3 * numNodes);
	CPUGPUMALLOC(acceleration,		 double, 3 * numNodes);
	CPUGPUMALLOC(acceleration_body1, double, 3 * numNodes);
	CPUGPUMALLOC(acceleration_body2, double, 3 * numNodes);
	CPUGPUMALLOC(force,				 double, 3 * numNodes);
	CPUGPUMALLOC(force_body1,		 double, 3 * numNodes);
	CPUGPUMALLOC(force_body2,		 double, 3 * numNodes);
	CPUGPUMALLOC(momentum,		     double, 3 * numNodes);
	CPUGPUMALLOC(momentum_body1,	 double, 3 * numNodes);
	CPUGPUMALLOC(momentum_body2,	 double, 3 * numNodes);
	CPUGPUMALLOC(mass,				 double, 1 * numNodes);
	CPUGPUMALLOC(mass_body1,		 double, 1 * numNodes);
	CPUGPUMALLOC(mass_body2,		 double, 1 * numNodes);
	CPUGPUMALLOC(normals_body1,		 double, 3 * numNodes);
	CPUGPUMALLOC(normals_body2,		 double, 3 * numNodes);
	CPUGPUMALLOC(temperature,		 double, 1 * numNodes);
	CPUGPUMALLOC(temperatureRate,	 double, 1 * numNodes);
	CPUGPUMALLOC(Tmass,				 double, 1 * numNodes);
	CPUGPUMALLOC(cell_pcount,		 int,	 1 * numNodes);
	CPUGPUMALLOC(cells,				 int,	 1 * numNodes);
	CPUGPUMALLOC(bodyID,			 int,	 2 * numNodes);
	CPUGPUMALLOC(hasRigid,			 bool,	 2 * numNodes);
	CPUGPUMALLOC(cellList,			 int,	 1 * numNodes);
	CPUGPUMALLOC(BCTYPE,			 int,	 1 * numNodes);
	CPUGPUMALLOC(BCVAL,				 double, 4 * numNodes);
	// position         = new double[3*numNodes];
	// velocity         = new double[3*numNodes];
	// velocity_body1   = new double[3*numNodes];
	// velocity_body2   = new double[3*numNodes];
	// acceleration     = new double[3*numNodes];
	// acceleration_body1     = new double[3*numNodes];
	// acceleration_body2     = new double[3*numNodes];
	// force            = new double[3*numNodes];
	// force_body1            = new double[3*numNodes];
	// force_body2            = new double[3*numNodes];
	// momentum         = new double[3*numNodes];
	// momentum_body1   = new double[3*numNodes];
	// momentum_body2   = new double[3*numNodes];
	// mass             = new double[numNodes];
	// mass_body1       = new double[numNodes];
	// mass_body2       = new double[numNodes];
	// normals_body1    = new double[3*numNodes]; // assumining only two bodies will be interacting per node (2 * vectors)
	// normals_body2	 = new double[3*numNodes]; // assumining only two bodies will be interacting per node (2 * vectors)
	// temperature      = new double[numNodes];
	// temperatureRate  = new double[numNodes];
	// Tmass		     = new double[numNodes];
	// cell_pcount      = new int[numNodes];// Number of particles per each cell
	// cells            = new int[numNodes];
	// bodyID            = new int[2*numNodes];
	// hasRigid		 = new bool[2*numNodes];
	// cellList         = new int[numNodes]; 
	// BCTYPE	         = new int[numNodes];
	// BCVAL	         = new double[4*numNodes];
	this->Nz         = (1+(L[2]/cellsize));
	this->Ny         = (1+(L[1]/cellsize));	
    this->Nx         = (1+(L[0]/cellsize));
	printf("Number of cells: \nIn X: %d\nIn Y: %d\nIn Z: %d\n",this->Nx, this->Ny,this->Nz);
    // Generate grid
    int ncount=0; // Dummy index for node arrays
	for (int k = 0; k < this->Nz; k++)
	{
		for (int j = 0; j < this->Ny; j++)
		{
			for (int i = 0; i < this->Nx; i++)
			{
				double x                   = i*cellsize;
				double y                   = j*cellsize;
				double z                   = k*cellsize;
				position[ncount*3 + 0]     = x;
				position[ncount*3 + 1]     = y;
				position[ncount*3 + 2]     = z;
				velocity[ncount*3 + 0]     = 0.0;
				velocity[ncount*3 + 1]     = 0.0;
				velocity[ncount*3 + 2]     = 0.0;
				velocity_body1[ncount*3 + 0]     = 0.0;
				velocity_body1[ncount*3 + 1]     = 0.0;
				velocity_body1[ncount*3 + 2]     = 0.0;
				velocity_body2[ncount*3 + 0]     = 0.0;
				velocity_body2[ncount*3 + 1]     = 0.0;
				velocity_body2[ncount*3 + 2]     = 0.0;
				momentum_body1[ncount*3 + 0]     = 0.0;
				momentum_body1[ncount*3 + 1]     = 0.0;
				momentum_body1[ncount*3 + 2]     = 0.0;
				momentum_body2[ncount*3 + 0]     = 0.0;
				momentum_body2[ncount*3 + 1]     = 0.0;
				momentum_body2[ncount*3 + 2]     = 0.0;
				acceleration[ncount*3 + 0] = 0.0;
				acceleration[ncount*3 + 1] = 0.0;
				acceleration[ncount*3 + 2] = 0.0;
				acceleration_body1[ncount*3 + 0] = 0.0;
				acceleration_body1[ncount*3 + 1] = 0.0;
				acceleration_body1[ncount*3 + 2] = 0.0;
				acceleration_body2[ncount*3 + 0] = 0.0;
				acceleration_body2[ncount*3 + 1] = 0.0;
				acceleration_body2[ncount*3 + 2] = 0.0;
				force[ncount*3 + 0]        = 0.0;
				force[ncount*3 + 1]        = 0.0;
				force[ncount*3 + 2]        = 0.0;
				force_body1[ncount*3 + 0]        = 0.0;
				force_body1[ncount*3 + 1]        = 0.0;
				force_body1[ncount*3 + 2]        = 0.0;
				force_body2[ncount*3 + 0]        = 0.0;
				force_body2[ncount*3 + 1]        = 0.0;
				force_body2[ncount*3 + 2]        = 0.0;
				momentum[ncount*3 + 0]     = 0.0;
				momentum[ncount*3 + 1]     = 0.0;
				momentum[ncount*3 + 2]     = 0.0;
				mass[ncount]               = 0.0;
				mass_body1[ncount]         = 0.0;
				mass_body2[ncount]         = 0.0;
				temperature[ncount]        = 0.0;
				temperatureRate[ncount]    = 0.0;
				Tmass[ncount]              = 0.0;
				BCTYPE[ncount]             = 0;
				bodyID[2*ncount + 0]       = -1;
				bodyID[2*ncount + 1]       = -1;
				hasRigid[2*ncount + 0]     = false;
				hasRigid[2*ncount + 1]     = false;
				ncount++;
			}
		}
	}
   
};

Nodes::~Nodes(){

}
// Getters
double* Nodes::getPosArr(){return position;}
double* Nodes::getMassArr(){return mass;};
double* Nodes::getVelArr(){return velocity;};
double* Nodes::getForceArr(){return force;};
double  Nodes::getPosX(int n){return position[n*3 + 0];};
double  Nodes::getPosY(int n){return position[n*3 + 1];};
double  Nodes::getPosZ(int n){return position[n*3 + 2];};

void Nodes::setBC(int nIdx, BCTypes::function type, std::vector<double> val )
{
	BCTYPE[nIdx]    = type;
	for(int i=0; i < val.size(); i++)
	{
		BCVAL[4*nIdx+i] = val[i];
	};
	
}

void Nodes::Serialise(int step, char * filename){
FILE *fp_nds;
    char filename_nds[50]; 
    sprintf(filename_nds,"output/%s_nds_%d.csv",filename,step);
    fp_nds=fopen(filename_nds,"w+");
    fprintf(fp_nds, "Pos_x,Pos_y,Pos_z\n");
    for(int i=0; i<this->numNodes;i++)
    {
        fprintf(fp_nds, "%f,%f,%f\n",this->position[i*3 +0],this->position[i*3 +1],this->position[i*3 +2]);
    }
    fclose(fp_nds);

}

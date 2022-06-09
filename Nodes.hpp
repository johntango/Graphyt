// Header file for Nodes class
#ifndef NODES_H
#define NODES_H

#include <vector>
#include "BoundaryConditions.hpp"
class Nodes
{
	public:
	Nodes(std::vector<double> L, double cellsize);
    ~Nodes();
	void                solveNodalEquations(Nodes *nds, double dt, int step);
	void                setBC(int nIdx, BCTypes::function type, std::vector<double> val );
	void                Serialise(int step, char* filename);
    double              getPosX(int n);
	double              getPosY(int n);
	double              getPosZ(int n);
	double*             getPosArr();
    double*             getMassArr();
    double*             getVelArr();
	double*             getForceArr();
    std::vector<double> L;// Lengths of the grid in each dimension (from 0)
    double              cellsize;
    double              t;
    int                 numNodes;
    int                 numCells;
   	double*             position;
	double*             velocity; // Centre of mass velocity field
	double*             velocity_body1;
	double*             velocity_body2;
	double*             momentum_body1; // body 1 velocity field
	double*             momentum_body2; // body 2 velocity field
	int*                bodyID;
	bool * 				hasRigid;
	double*             acceleration;
	double*             acceleration_body1;
	double*             acceleration_body2;
	double*             force;
	double*             force_body1;
	double*             force_body2;
	double*             momentum;
	double*             mass;
	double*             mass_body1;
	double*             mass_body2;
	double*				normals_body1;
	double*				normals_body2;
	double*             temperature;
	double*             temperatureRate;
	double*             Tmass;
	int*                cell_pcount;// Number of particles per each cell
	int*                cells;
	int*                cellList; 
	int*                BCTYPE;
	double*             BCVAL;
	int                 Nx,Ny,Nz;

};


#endif
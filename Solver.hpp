#ifndef SOLVER_H
#define SOLVER_H
#include "MaterialPoints.hpp"
#include "Nodes.hpp"
#include "Parameters.hpp"
#include "TimeIntegration.hpp"
#include "Extrapolation.hpp"
#include "Grid.hpp"
#include "BoundaryConditions.hpp"
#include "MaterialModels.hpp"




class Solver
{
	public:
	Solver(Nodes* nds, MaterialPoints* mpts,Parameters* mpmParameters, std::vector<MaterialModels> materials);
    ~Solver();
    void                        initialize();
    void                        iterate_GPU(double t,int step, Parameters* mpmParameters);
    void                        iterate_CPU_steps(int nsteps, Parameters* mpmParameters);
    void                        iterate_CPU(double t,int step, Parameters* mpmParameters);
    void                        setDt(double dt);
    void                        enforceBCs(MaterialPoints* mpts, Nodes* nds, int phase);
    void						Synchronise();
    void                        sync_GPU();
    Nodes*                      nds;
    MaterialPoints*             mpts;
    Parameters*                 mpmParameters;
    std::vector<MaterialModels> materials;
    double                      dt;
    double                      tmax;
    double                      startTimeLoop;
    double                      avgMillionParticleUpdatesPer_second;
    double                      avgtpcore;
};



#endif
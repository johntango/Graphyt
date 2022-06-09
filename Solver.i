%module Solver

%{
#include "Solver.hpp"    
#include "MaterialModels.hpp"    
%}

%template(vectorm) std::vector<MaterialModels>;

class Solver{
    public:
	Solver(Nodes* nds, MaterialPoints* mpts,Parameters* mpmParameters, std::vector<MaterialModels> materials);
    ~Solver();
    void iterate_CPU(double t, int step,Parameters* mpmParameters);
    void iterate_CPU_steps(int nsteps, Parameters* mpmParameters);
    void iterate_GPU(double t, int step,Parameters* mpmParameters);
    void setDt(double dt);
    double dt,tmax,startTimeLoop;
};
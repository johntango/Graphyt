%module Nodes

%{
#include "Nodes.hpp"    
#include "BoundaryConditions.hpp"
%}

class Nodes{
    public:
	Nodes(std::vector<double> L, double cellsize);
    ~Nodes();
    void   setBC(int nIdx, BCTypes::function type, std::vector<double> val );
    double getPosX(int n);
	double getPosY(int n);
	double getPosZ(int n);
    double * getPosArr();
    double * getMassArr();
    double * getVelArr();
    double * getForceArr();
	void Serialise(int step, char * filename);
    double cellsize;
    int numNodes;
};
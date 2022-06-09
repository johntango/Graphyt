%module Parameters

%{
#include "Parameters.hpp"    
%}

class Parameters{
    public:
	Parameters(double tmax,bool is3D);
	~Parameters();
	void setGravity(double gx, double gy, double gz);
	void    setDampingCoef(double c);
	void    setFriction(double mu);
	void    setCFL(double cfl);
};

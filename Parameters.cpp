#include "Parameters.hpp"
#include <cstdlib>
#include <vector>
#include "omp.h"
// Header file to define all data structs used in the simulation:
// - Simulation Parameters

Parameters::Parameters(double tmax,bool is3D)
{
	this->tmax                = tmax;
	this->is3D          	  = is3D;
	this->CFL				  = 0.1;
	this->mu                  = 0.0;
	this->velocityDampingCoef = 0.0;
};

Parameters::~Parameters(){};

void Parameters::setGravity(double gx, double gy, double gz){
	gravity[0] = gx;
	gravity[1] = gy;
	gravity[2] = gz;
};

void Parameters::setDampingCoef(double c){
	this->velocityDampingCoef = c;
};

void Parameters::setCFL(double cfl){
	this->CFL = cfl;
};

void Parameters::setFriction(double mu){
	this->mu = mu;
}

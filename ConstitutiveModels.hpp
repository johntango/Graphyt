#ifndef CONSTITUTIVE_H
#define CONSTITUTIVE_H

#include "MaterialPoints.hpp"
#include "MaterialModels.hpp"
#include "Parameters.hpp"
#include "Nodes.hpp"
#include "eig3.h"
#include "Matrix.h"
#include <math.h>
#include "omp.h"
#include <algorithm>
#include <stdio.h>


CPUGPU void elastic(Parameters *mpmParameters,MaterialPoints* mpts,double dt,int mp,std::vector<MaterialModels>* materials);
CPUGPU void newtonianLiquid(Parameters *mpmParameters,MaterialPoints* mpts,double dt,int mp,std::vector<MaterialModels>* materials);
CPUGPU void perfectPlasticity(Parameters *mpmParameters,MaterialPoints* mpts,double dt,int mp,std::vector<MaterialModels>* materials);
CPUGPU void druckerPrager_Damage(Parameters *mpmParameters,MaterialPoints* mpts,double dt,int mp,std::vector<MaterialModels>* materials);
CPUGPU void turbulentLiquid(Parameters *mpmParameters,MaterialPoints* mpts,double dt,int mp,std::vector<MaterialModels>* materials);

#endif
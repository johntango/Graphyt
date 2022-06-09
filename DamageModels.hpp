#ifndef DAMAGE_H
#define DAMAGE_H

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


CPUGPU void gradyKippDamage(Parameters *mpmParameters,MaterialPoints* mpts,double dt,int mp,std::vector<MaterialModels>* materials);

#endif
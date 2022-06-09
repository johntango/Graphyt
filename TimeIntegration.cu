#include "TimeIntegration.hpp"
#include "cuda.h"
#include "MaterialPoints.hpp"
#include "Nodes.hpp"
#include "eig3.h"
#include "Matrix.h"
#include <math.h>
#include "omp.h"
#include <algorithm>
#include "managed_allocator.hpp"
#include <stdio.h>

#include "TimeIntegration.cpp"
__global__ void updateMatPointKernel(Parameters *mpmParameters,MaterialPoints *mpts, Nodes *nds,double dt) 
{

}

void updateMatPoint_GPU(Parameters *mpmParameters,MaterialPoints *mpts, Nodes *nds,double dt, std::vector<MaterialModels>* materials){

}


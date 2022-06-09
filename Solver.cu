#include "cuda.h"
#include "Solver.cpp"


void Solver::sync_GPU() {

		cudaDeviceSynchronize();
};
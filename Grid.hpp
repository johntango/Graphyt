#ifndef GRID_H
#define GRID_H
#include "Parameters.hpp"
#include "MaterialPoints.hpp"
#include "Nodes.hpp"
// Header file for Grid.cpp
CPUGPU void solveNodalEquations(Parameters* mpmParameters,Nodes *nds, double dt,int step,int n);
CPUGPU void calculateColoumbFrictionContactNormals(double* friction_normals,double* contactNormals,double mu_prime,double* relVel_body, Nodes* nds,Parameters* parameters, int n,int body);
CPUGPU void collectPeriodicContributions(Parameters* mpmParameters,Nodes *nds);
CPUGPU void resetNode(Nodes* nds, int n);
CPUGPU void enforceGridBC(Nodes* nds,int phase, int n);
// BC Type Functions
CPUGPU void bc_grid_vx(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_vy(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_vz(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_temp(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_tflux_x(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_tflux_y(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_tflux_z(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_fx(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_fy(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_fz(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_vx_vy(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_vx_vz(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_vy_vz(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_fx_fy(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_fy_fz(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_fx_fz(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_fx_fy_fz(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_vx_vy_vz(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_vx_temp(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_vy_temp(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_vz_temp(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_fx_temp(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_fy_temp(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_fz_temp(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_vx_vy_temp(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_vx_vz_temp(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_vy_vz_temp(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_fx_fy_temp(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_fy_fz_temp(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_fx_fz_temp(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_fx_fy_fz_temp(Nodes* nds, int ndIdx,int phase);
CPUGPU void bc_grid_vx_vy_vz_temp(Nodes* nds, int ndIdx,int phase);
#endif
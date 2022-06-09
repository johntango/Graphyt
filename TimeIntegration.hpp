#ifndef TIMEINTEGRATE_H
#define TIMEINTEGRATE_H
#include "Parameters.hpp"
#include "MaterialPoints.hpp"
#include "BoundaryConditions.hpp"
#include "Nodes.hpp"
#include "MaterialModels.hpp"

#include "ConstitutiveModels.hpp"
#include "DamageModels.hpp"

// Header file for TimeIntegration.cpp
CPUGPU void updateMatPoint(Parameters *mpmParameters,MaterialPoints *mpts, Nodes *nds,double dt,std::vector<MaterialModels> *materials,int mp);

 void updateMatPoint_GPU(Parameters *mpmParameters,MaterialPoints *mpts, Nodes *nds,double dt, std::vector<MaterialModels>* materials);

CPUGPU void constitutiveModel(Parameters *mpmParameters, MaterialPoints *mpts, double dt, int mp, std::vector<MaterialModels>* materials);

CPUGPU void applyDamage(Parameters *mpmParameters,MaterialPoints* mpts,double dt,int mp,std::vector<MaterialModels>* materials);
CPUGPU void resetMaterialPointRates(MaterialPoints* mpts, int mp);
CPUGPU void enforceMatPointBC(MaterialPoints* mpts, int phase, int mp);
CPUGPU void bc_mp_vx(MaterialPoints* mpts, int mp, int phase);
CPUGPU void bc_mp_vxsin(MaterialPoints* mpts, int mp, int phase);
CPUGPU void bc_mp_vy(MaterialPoints* mpts, int mp, int phase);
CPUGPU void bc_mp_vysin(MaterialPoints* mpts, int mp, int phase);
CPUGPU void bc_mp_syysin(MaterialPoints* mpts, int mp, int phase);
CPUGPU void bc_mp_fx(MaterialPoints* mpts, int mp, int phase);
CPUGPU void bc_mp_fy(MaterialPoints* mpts, int mp, int phase);
CPUGPU void bc_mp_fz(MaterialPoints* mpts, int mp, int phase);
CPUGPU void bc_mp_stress_xx_yy(MaterialPoints* mpts, int mp, int phase);
CPUGPU void bc_mp_fx_fy(MaterialPoints* mpts, int mp, int phase);
CPUGPU void bc_mp_vz(MaterialPoints* mpts, int mp, int phase);

CPUGPU void mp_vz(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fx(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fy(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fz(MaterialPoints* mpts, int mp, int phase);
CPUGPU void bc_mp_vx_vy(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_vx_vz(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_vy_vz(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fx_fy(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fy_fz(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fx_fz(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fx_fy_fz(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_vx_vy_vz(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_vx_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_vy_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_vz_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fx_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fy_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fz_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_vx_vy_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_vx_vz_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_vy_vz_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fx_fy_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fy_fz_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fx_fz_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fx_fy_fz_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_vx_vy_vz_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_vx_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_vy_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_vz_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fx_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fy_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fz_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_vx_vy_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_vx_vz_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_vy_vz_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fx_fy_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fy_fz_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fx_fz_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_fx_fy_fz_temp(MaterialPoints* mpts, int mp, int phase);
CPUGPU void mp_vx_vy_vz_temp(MaterialPoints* mpts, int mp, int phase);

#endif

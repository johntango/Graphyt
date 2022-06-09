#include "MaterialModels.hpp"
#include <stdio.h>

MaterialModels::MaterialModels() // constructor
{
    // Intialise a material with defaul quantities
    this->bulkMod              = 1e7;
    this->shearMod             = 5000.0;
    this->density              = 1000.0;
    this->isRigid              = false;
    this->viscosity            = 1.0;
    this->thermalConductivity  = 391;  // Thermal properties for Copper
    this->specificHeatCapacity = 385.2;// Thermal properties for Copper
    this->yieldStress          = 1e7;
    this->phi                  = 40.0;
    this->psi                  = 10.0;
    this->cohesion             = 100000.0;
    this->criticalStrain       = 0.01;
    this->isDamage             = false;
    this->damageModel          = 1;
    this->Cg                   = 586.52;
    this->damageM              = 9;
    this->damageK              = 5.39e36;
    this->turbStress           = 1e10;
    this->model                = 1;
};

MaterialModels::~MaterialModels(){};

// Functions
void MaterialModels::setBulkMod(double bulkMod){this->bulkMod=bulkMod;};
void MaterialModels::setShearMod(double shearMod){this->shearMod=shearMod;};
void MaterialModels::setDensity(double density){this->density=density;};
void MaterialModels::setIsRigid(bool isRigid){this->isRigid=isRigid;};
void MaterialModels::setViscosity(double viscosity){this->viscosity=viscosity;};
void MaterialModels::setThermalConductivity(double thermalConductivity){this->thermalConductivity=thermalConductivity;};
void MaterialModels::setSpecificHeatCapacity(double specificHeatCapacity){this->specificHeatCapacity=specificHeatCapacity;};
void MaterialModels::setYieldStress(double yieldStress){this->yieldStress=yieldStress;};
void MaterialModels::setCohesion(double cohesion){this->cohesion = cohesion;};
void MaterialModels::setCriticalStrain(double criticalStrain){this->criticalStrain = criticalStrain;};
void MaterialModels::setFrictionAngle(double phi){this->phi = phi;};
void MaterialModels::setDilatancyAngle(double psi){this->psi = psi;};
void MaterialModels::setIsDamage(bool isDamage){this->isDamage = isDamage;};
void MaterialModels::setDamageModel(int model){this->damageModel = model;};
void MaterialModels::setCrackSpeed(double Cg){this->Cg = Cg;};
void MaterialModels::setDamageM(double m){this->damageM = m;};
void MaterialModels::setDamageK(double k){this->damageK = k;};
void MaterialModels::setTurbulentStress(double turbStress){this->turbStress = turbStress;};
void MaterialModels::setModel(int model){this->model = model;};



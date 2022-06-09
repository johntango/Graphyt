#ifndef MATMODELS_H
#define MATMODELS_H


class MaterialModels
{
	public:
            MaterialModels();
            ~MaterialModels();
            void   setBulkMod(double bulkMod);
            void   setShearMod(double shearMod);
            void   setDensity(double density);
            void   setIsRigid(bool isRigid);
            void   setViscosity(double viscosity);
            void   setThermalConductivity(double thermalConductivity);
            void   setSpecificHeatCapacity(double specificHeatCapacity);
            void   setYieldStress(double yieldStress);
            void   setCohesion(double cohesion);
            void   setCriticalStrain(double criticalStrain);
            void   setFrictionAngle(double phi);
            void   setDilatancyAngle(double psi);
            void   setIsDamage(bool isDamage);
            void   setDamageModel(int model);
            void   setCrackSpeed(double Cg);
            void   setDamageM(double m);
            void   setDamageK(double k);
            void   setTurbulentStress(double turbStress);
            void   setModel(int model);
            double bulkMod;
            double shearMod;
            double density;
            bool   isRigid;
            double viscosity;
            double thermalConductivity;
            double specificHeatCapacity;
            double yieldStress;
            double cohesion;
            double criticalStrain;
            double phi;
            double psi;
            bool   isDamage;
            int    damageModel;
            double Cg;
            double damageM;
            double damageK;
            double turbStress;
            int    model;
};


#endif
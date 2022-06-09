%module MaterialModels

%{
    #include "MaterialModels.hpp"    
%}

class MaterialModels{
    public:
            MaterialModels();
            ~MaterialModels();
            void setBulkMod(double bulkMod);
            void setShearMod(double shearMod);
            void setDensity(double density);
            void setIsRigid(double isRigid);
            void setViscosity(double viscosity);
            void setThermalConductivity(double thermalConductivity);
            void setSpecificHeatCapacity(double specificHeatCapacity);
            void setYieldStress(double yieldStress);
            void setCohesion(double cohesion);
            void setCriticalStrain(double criticalStrain);
            void setFrictionAngle(double phi);
            void setDilatancyAngle(double psi);
            void   setIsDamage(bool isDamage);
            void   setDamageModel(int model);
            void   setCrackSpeed(double Cg);
            void   setDamageM(double m);
            void   setDamageK(double k);
            void setTurbulentStress(double turbStress);
            void setModel(int model);

};
%module MaterialPoints

%{
	#include "MaterialPoints.hpp"    
	#include "BoundaryConditions.hpp"
%}

class MaterialPoints{
    public:
			MaterialPoints(Nodes*);
			~MaterialPoints();
			void    add(int i,double x, double y,double z);
			void    addPyckGeom(double* pos, int numParticle);
			void    addPyckIntField(IntField* intField);
			void    addPyckDoubleField(DoubleField* doubleField);
			double  getPosX(int p);
			double  getPosY(int p);
			double  getPosZ(int p);
            double  getDispX(int p);
            double  getDispY(int p);
            double  getDispZ(int p);
            double  getVelX(int p);
            double  getVelY(int p);
            double  getVelZ(int p);
            double  getStressXX(int p);
            double  getStressYY(int p);
            double  getStressZZ(int p);
            double  getStressXY(int p);
            double  getStressYZ(int p);
            double  getStressZX(int p);
			int     getObjID(int p);
			int     getMatID(int p);
			int		getCellID(int p);
			double* getPosArr();
			double* getMassArr();
			int*    getMatArr();
			int*    getCellIDArr();
			double* getVelArr();
			double* getAccArr();
			double* getDispArr();
			double* getStressArr();
			double* getStrainArr();
			double* getTemperatureArr();
			double* getPlasticStrainArr();
			double* getAccPlasticStrainArr();
			double* getDamageArr();
			double* getSMaxArr();
			double* getNormalsArr();
			double*  getVolumeArr();
			double*  getDensityArr();
			int*    getObjectIDArr();
			void    setBC(int mp, BCTypes::function type, std::vector<double> val );
			void    setDamage(int p,double D);
			void    setVel(int p, double vx, double vy, double vz);
			void    setStress(int p, double sxx,double syy,double szz,double sxy,double syz,double sxz);
			void    setTemp(int p, double T);
			void    setRotatingBodyVel(int objectID, double centerX, double centerY, double centerZ, double omega);
			void    setMaterialByCell(std::vector<int> materialArray);
			void    setMaterialByCell(int idx, int materialID);
			void    setMaterialID(int p, int materialID);
			void 	setObjectID(int p, int objid);
			void    holdParticle(int p);
			void    Serialise(int step, char * filename);
			int     numParticles;
			double  psep;
};

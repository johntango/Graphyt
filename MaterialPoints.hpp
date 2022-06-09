#ifndef MATERIALPOINT_H
#define MATERIALPOINT_H
#include "Nodes.hpp"
#include "MaterialModels.hpp"
#include "BoundaryConditions.hpp"
#include <vector>
#include "intField.h"
#include "doubleField.h"


#ifdef __CUDACC__
	#include "managed_allocator.hpp"
	template<class T>
	using managed_vector = std::vector<T, managed_allocator<T>>;
#else
	template<class T>
	using managed_vector = std::vector<T>;
#endif

class MaterialPoints
{
	public:
	MaterialPoints(Nodes*);
	~MaterialPoints();
	void     add(int i, double x, double y,double z);
	void     addPyckGeom(double* pos, int numParticle);
	void     addPyckIntField(IntField* intField);
	void     addPyckDoubleField(DoubleField* doubleField);
	double   getPosX(int p);
	double   getPosY(int p);
	double   getPosZ(int p);
    double   getDispX(int p);
	double   getDispY(int p);
	double   getDispZ(int p);
    double   getVelX(int p);
	double   getVelY(int p);
	double   getVelZ(int p);
    double   getStressXX(int p);
	double   getStressYY(int p);
	double   getStressZZ(int p);
    double   getStressXY(int p);
	double   getStressYZ(int p);
	int      getMatID(int p);
	double   getStressZX(int p);
	int      getObjID(int p);
	int		 getCellID(int p);

	double*  getPosArr();
	double*  getDispArr();
	double*  getMassArr();
	int*     getMatArr();
	int*     getCellIDArr();
	double*  getVelArr();
	double*  getAccArr();
	double*  getStressArr();
	double*  getStrainArr();
	double*  getTemperatureArr();
	double*  getPlasticStrainArr();
	double*  getAccPlasticStrainArr();
	double*  getDamageArr();
	double*  getSMaxArr();
	double*  getNormalsArr();
	double*  getSurfaceParticles();
	double*  getVolumeArr();
	double*  getDensityArr();
	int*     getObjectIDArr();
	void     setRotatingBodyVel(int objectID, double centerX, double centerY, double centerZ, double omega);
	void     setVelBC(int p, double vx, double vy, double vz);
	void     setBC(int p, BCTypes::function type, std::vector<double> val );
	void     setDamage(int p,double D);
	void     setVel(int p, double vx, double vy, double vz);
	void     setTemp(int p, double T);
	void     setStress(int p, double sxx,double syy,double szz,double sxy,double syz,double sxz);
	void     setMaterialByCell(std::vector<int> materialArray);
	void     setMaterialByCell(int idx, int materialID);
	void     setMaterialID(int idx, int materialID);
	void 	 setObjectID(int p, int objid);
	void     holdParticle(int p);
	void     Serialise(int ste, char* filename);
	Nodes*   nodes;
	int      numParticles;
    double   psep;
    double   t;
	managed_vector<double>  position;
	managed_vector<double>  pos0;
	managed_vector<double>  pos0temp;
	managed_vector<double>  displacement;
	managed_vector<double>  disptemp;
	managed_vector<double>  velocity;
	managed_vector<double>  Gvelocity;
	managed_vector<double>  acceleration;
	managed_vector<double>  strainRate;
	managed_vector<double>  densityRate;
	managed_vector<double>  stress;
	managed_vector<double>  plastic_strain;
	managed_vector<double>  acc_plastic_strain;
	managed_vector<double>  acc_plastic_strain_temp;
	managed_vector<double>  damage;
	managed_vector<double>  damagetemp;
	managed_vector<double>  smax;
	managed_vector<double>  smaxtemp;
	managed_vector<double>  strain;
	managed_vector<double>  temperature;
	managed_vector<double>  temperatureRate;
	managed_vector<double>  temperatureGrad;
	managed_vector<double>  mass;
	managed_vector<double>  density;
	managed_vector<double>  volume;
	managed_vector<double>  normals;
	managed_vector<double>  normalstemp;
	managed_vector<int>     cellID;
	managed_vector<int>		material;
	managed_vector<int>     materialtemp;
	managed_vector<int>     objID;
	managed_vector<int>     objIDtemp;
	managed_vector<bool>    hold;	
	managed_vector<bool>    holdtemp;
	managed_vector<double>  postemp;
	managed_vector<double>  veltemp;
  	managed_vector<double>  gveltemp;
  	managed_vector<double>  stresstemp;
  	managed_vector<double>  straintemp;
	managed_vector<double>  temptemp;
	managed_vector<double>  tempRatetemp;
	managed_vector<double>  tempGradtemp;
	managed_vector<double>  plastictemp;
  	managed_vector<double>  masstemp;
  	managed_vector<double>  densitytemp;
  	managed_vector<double>  volumetemp;
  	managed_vector<int>     cellIDtemp;
  	managed_vector<int>     cell_pidx;// Cell ID for correspoinding particle
	managed_vector<int>     pidx;  // Index of particle ids to be sorted by cellID
	int*                    N;
	managed_vector<int>     BCTYPE;
	managed_vector<int>     BCTYPEtemp;
	managed_vector<double>  BCVAL;
	managed_vector<double>  BCVALtemp;

	
};


#endif

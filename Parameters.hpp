#ifndef PARAMETERS_H
#define PARAMETERS_H



#ifndef CPUGPU
		#define CPUGPU  // prefix to functions so that they will run on host and/or device
		#define CPUGPUMALLOC(pointer, type, count) pointer = new type[count]; // memory managment for arrays using compatible allocators
#endif
#ifdef __CUDACC__ 
		#define CPUGPU __host__ __device__
		#define CPUGPUMALLOC(pointer, type, count) cudaMallocManaged(pointer, count*sizeof(type));
#endif


class Parameters
{
	public:
	Parameters(double tmax, bool is3D);
	~Parameters();
	void    setGravity(double gx, double gy, double gz);
	void    setDampingCoef(double c);
	void    setCFL(double cfl);
	void    setFriction(double mu);
	double  dt;
	double  tmax;
	double  CFL;
	bool    is3D;
	double  mu;
	double  gravity[3];
	double  velocityDampingCoef;
};

#endif

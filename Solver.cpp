
#include "Solver.hpp"
#include <cmath>
#include <stdio.h>
#include <math.h> 

Solver::Solver(Nodes* nds, MaterialPoints* mpts, Parameters* mpmParameters, std::vector<MaterialModels> materials){
    this->mpts          = mpts;
    this->nds           = nds;
    this->materials     = materials;
    this->mpmParameters    = mpmParameters;
    this->startTimeLoop = omp_get_wtime(); // Get start time of total timeloop
    this->avgMillionParticleUpdatesPer_second = 0.0;
    this->avgtpcore = 0.0;
    double maxK         = 1;
    double minrho       = 1e15;
    
    for(int i = 0; i<materials.size();i++)
    {
        
        if(materials[i].bulkMod > maxK)
        {
            maxK   = materials[i].bulkMod;
        }
        
        if(materials[i].density<minrho)
        {
            minrho = materials[i].density;
        }
    }
    
    this->dt   = mpmParameters->CFL*(nds->cellsize / sqrt(maxK/minrho));
    this->tmax = mpmParameters->tmax;
    
    printf("Time step(ms): %f\n",1000.0*this->dt);
    initialize(); // Initialize the material point masses and volumes
}
Solver::~Solver(){};

void Solver::initialize(){
    // Set the mass and volumes of matpoints based on their materials
    if(this->mpmParameters->is3D)
    {
        #pragma omp for
        for(int p = 0; p<this->mpts->numParticles; p++)
        {
            int matID = mpts->material[p];
            mpts->density[p] = materials[matID].density;
            mpts->volume[p]  = mpts->psep*mpts->psep*mpts->psep;
            mpts->mass[p]    = mpts->density[p]*mpts->volume[p];
            if(materials[matID].isRigid)
            {
                mpts->mass[p]    = 1.0;
            }
        }
    }
    else
    {
        #pragma omp for
        for(int p = 0; p<this->mpts->numParticles; p++)
        { // For 2D calculations, volume and mass is per unit depth
            int matID = mpts->material[p];
            mpts->density[p] = materials[matID].density;
            mpts->volume[p] = mpts->psep*mpts->psep;
            mpts->mass[p]   = mpts->density[p]*mpts->volume[p];
        }
    }

};
void Solver::setDt(double dt)
{
    this->dt = dt;
};

void Solver::iterate_CPU_steps(int nsteps, Parameters* mpmParameters)
{
    double t = 0.0;
    double dt = mpmParameters->dt;
    for(int i = 0; i < nsteps; i++)
    {
        iterate_CPU(t, i, mpmParameters);
        t = t + dt;
    }
}

void Solver::iterate_CPU(double t, int step,Parameters* mpmParameters)
{
    /////////////////////////////////////////// Sort Particles  ///////////////////////////////////////////////////////
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
        double startStepTime = omp_get_wtime();
        if((step%10==0)||(step==1))
       {
        sortParticlesByCellID_CPU(mpmParameters,mpts,this->nds);
       }
    //----------------------------------------------------------------------------------------------------------------//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //#
    //#
    //#
    //#
    //#
    /////////////////////////////////////////////// Reset Nodes ///////////////////////////////////////////////////////
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
        double timecheckpoint_1 = omp_get_wtime();
        //#pragma omp parallel for schedule(dynamic, this->nds->numNodes/(omp_get_num_threads()*16))
		#pragma omp parallel for 
        for(int n = 0; n < this->nds->numNodes; n++)
        {
            resetNode(this->nds,n);// in Grid.c/hpp
        }
        double timecheckpoint_2 = omp_get_wtime();
    //----------------------------------------------------------------------------------------------------------------//
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //#
    //#
    //#
    //#
    //#
    /////////////////////////////////////////////// Reset Material Point Rates ///////////////////////////////////////////////////////
        #pragma omp parallel for schedule(dynamic,mpts->numParticles/(omp_get_num_threads()*16))
		for(int mp=0; mp<mpts->numParticles; mp++)
        {
          resetMaterialPointRates(mpts,mp);
        }
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //#
    //#
    //#
    //#
    //#
    ////////////////////////////////// Extrapolate Material Points To Nodes ////////////////////////////////////////////
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
        // For each node, find the cells it needs to use and populate a list of cell IDs to use for particle info
        // 4x4x4 cells in 3D.
        #pragma omp parallel for schedule(dynamic, this->nds->numNodes/(omp_get_num_threads()*16))
        for(int n = 0; n<nds->numNodes;n++)
        {
            extrapToNodes(mpmParameters,mpts,nds,&(this->materials),n);
        };// Loop over each node (in parallel)
        double timecheckpoint_3 = omp_get_wtime();    
    //----------------------------------------------------------------------------------------------------------------//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //#
    //#
    //#
    //#
    //#
    //////////////////////////// Grid Boundary Conditions - Post Extrapolation To Nodes ////////////////////////////////
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    #pragma omp parallel for schedule(dynamic, this->nds->numNodes/(omp_get_num_threads()*16))
		for(int n=0; n< this->nds->numNodes; n++)
        {
            enforceGridBC(nds,1,n); // in Grid.cpp/hpp
        }
        double timecheckpoint_4 = omp_get_wtime();
    //----------------------------------------------------------------------------------------------------------------//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //#
    //#
    //#
    //#
    //#
    ////////////////////////////////////////// Solve Nodal Equations //////////////////////////////////////////////////
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    // Collect Periodic boundary contributions
        collectPeriodicContributions(mpmParameters,nds);
        #pragma omp parallel for schedule(dynamic, this->nds->numNodes/(omp_get_num_threads()*16))
		for(int n=0; n< this->nds->numNodes; n++)
        {
            solveNodalEquations(mpmParameters,this->nds,this->dt,step,n); // in Grid.c/hpp
        }
        double timecheckpoint_5 = omp_get_wtime();
    //----------------------------------------------------------------------------------------------------------------//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //#
    //#
    //#
    //#
    //#
    ////////////////////////////// Grid Boundary Conditions - Post Nodal Equations ////////////////////////////////////
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    #pragma omp parallel for schedule(dynamic, this->nds->numNodes/(omp_get_num_threads()*16))
		for(int n=0; n< this->nds->numNodes; n++)
        {
            enforceGridBC(nds,2,n); // in Grid.cpp/hpp
        }
        double timecheckpoint_6 = omp_get_wtime();
    //----------------------------------------------------------------------------------------------------------------//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //#
    //#
    //#
    //#
    //#
    /////////////////////////////////// Extrapolate Nodal Values to Material Points ///////////////////////////////////
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    #pragma omp parallel for schedule(dynamic,mpts->numParticles/(omp_get_num_threads()*16))
		for(int mp=0; mp<mpts->numParticles; mp++)
        {
          extrapToMatPoints(mpmParameters,mpts,nds,&(this->materials),mp);
        }
        double timecheckpoint_7 = omp_get_wtime();
    //----------------------------------------------------------------------------------------------------------------//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //#
    //#
    //#
    //#
    //#    
    //////////////////////////// Particle Boundary Conditions - Post Nodal Equations ///////////////////////////////////
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//  
    #pragma omp parallel for schedule(dynamic, this->nds->numNodes/(omp_get_num_threads()*16))
		for(int n=0; n< this->nds->numNodes; n++)
        {
            enforceGridBC(nds,3,n); // in Grid.cpp/hpp
        }

        #pragma omp parallel for //schedule(dynamic, this->mpts->numParticles/(omp_get_num_threads()*16))        
        for(int mp=0; mp<mpts->numParticles; mp++)
        {
            enforceMatPointBC(mpts,4,mp); // in TimeIntegration.cpp/hpp
        }
        double timecheckpoint_8 = omp_get_wtime();
    //----------------------------------------------------------------------------------------------------------------//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //#
    //#
    //#
    //#
    //#
    ////////////////////////////////////// Update Material Points /////////////////////////////////////////////////////
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//        
	    #pragma omp parallel for schedule(dynamic, mpts->numParticles/(omp_get_num_threads()*16))    
		for(int mp=0; mp<mpts->numParticles; mp++)
        {
            updateMatPoint(mpmParameters,mpts,this->nds,this->dt,&(this->materials),mp);// in TimeIntegration.c/hpp
        }
        double timecheckpoint_9 = omp_get_wtime();
    //----------------------------------------------------------------------------------------------------------------//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //#
    //#
    //#
    //#
    //#   
    //////////////////////////// Particle Boundary Conditions - Post Update of Particles ///////////////////////////////
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// 
    #pragma omp parallel for schedule(dynamic,mpts->numParticles/(omp_get_num_threads()*16))    	
		for(int mp=0; mp<mpts->numParticles; mp++)
        {
            enforceMatPointBC(mpts,5,mp); // in TimeIntegration.cpp/hpp
        }
        double timecheckpoint_10 = omp_get_wtime();
    //----------------------------------------------------------------------------------------------------------------//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //#
    //#
    //#
    //#
    //#        
    ////////////////////////////// General Timing and other final measures /////////////////////////////////////////////
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
        this->nds->t=t;
        mpts->t=t;
        double endStepTime = omp_get_wtime();
        double totStepTime = endStepTime - startStepTime;
        int threads;
        #pragma omp parallel
        threads = omp_get_num_threads();
		
        this->avgMillionParticleUpdatesPer_second += (mpts->numParticles / (1000000.0*totStepTime))/1000.0;
        this->avgtpcore += threads*(1000*totStepTime/(mpts->numParticles))/1000.0;
        if(step%1000==0)
        {
            printf("Step: %d, SimTime (ms) = %.1f |(%.2f)|MPups:%.3f|t/p/core:%.10f ms \n",step,t*1000,totStepTime,this->avgMillionParticleUpdatesPer_second,this->avgtpcore);
            this->avgMillionParticleUpdatesPer_second = 0.0;
            this->avgtpcore = 0.0;
        }
        double time_resetNodes, time_extraptoNodes, time_enforceBCs1, time_solveNodes, time_enforceBCs2, time_extraptoMP, time_enforceBCs3, time_updateMP, time_enforceBCs4;
        time_resetNodes    = timecheckpoint_2 - timecheckpoint_1;
        time_extraptoNodes = timecheckpoint_3 - timecheckpoint_2;
        time_enforceBCs1   = timecheckpoint_4 - timecheckpoint_3;
        time_solveNodes    = timecheckpoint_5 - timecheckpoint_4;
        time_enforceBCs2   = timecheckpoint_6 - timecheckpoint_5;
        time_extraptoMP    = timecheckpoint_7 - timecheckpoint_6;
        time_enforceBCs3   = timecheckpoint_8 - timecheckpoint_7;
        time_updateMP      = timecheckpoint_9 - timecheckpoint_8;
        time_enforceBCs4   = timecheckpoint_10 - timecheckpoint_9;
            //----------------------------------------------------------------------------------------------------------------//
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
#ifdef __CUDACC__
    void Solver::iterate_GPU(double t, int step,Parameters* mpmParameters)
    {
        updateMatPoint_GPU(mpmParameters,mpts,this->nds,this->dt,&(this->materials));// in TimeIntegration.cu
		cudaDeviceSynchronize();
        // Still in development
    }
#else
void Solver::iterate_GPU(double t, int step,Parameters* mpmParameters){};
#endif


void Solver::Synchronise(){};



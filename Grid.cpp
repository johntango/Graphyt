#include "Grid.hpp"
#include "MaterialPoints.hpp"
#include <algorithm>
#include "Nodes.hpp"
#include "omp.h"
#include "stdio.h"
#include "math.h"

CPUGPU void solveNodalEquations(Parameters* mpmParameters,Nodes *nds, double dt, int step,int n){
// for each node and solve for it's updated velocity and acceleration, if more than one body is present at a node, resolve the contact mechanics corrections
		if(nds->mass[n]!=0.0)
		{
				// Find current grid velocity from momenta
				nds->velocity[n*3 +0] = nds->momentum[n*3 + 0]/nds->mass[n];
				nds->velocity[n*3 +1] = nds->momentum[n*3 + 1]/nds->mass[n];
				nds->velocity[n*3 +2] = nds->momentum[n*3 + 2]/nds->mass[n];
				// Find current acceleration from forces
			// MAKE CORRECTIONS FOR NODES INVOLVED IN CONTACT
				// Normalise normal vectos
				double mag1 =sqrt(pow(nds->normals_body1[3*n + 0],2.0) + pow(nds->normals_body1[3*n + 1],2.0) + pow(nds->normals_body1[3*n + 2],2.0));
				double mag2 =sqrt(pow(nds->normals_body2[3*n + 0],2.0) + pow(nds->normals_body2[3*n + 1],2.0) + pow(nds->normals_body2[3*n + 2],2.0)) ;
				if(mag1>1e-8)
				{
					nds->normals_body1[3*n +0] /= mag1;
					nds->normals_body1[3*n +1] /= mag1;
					nds->normals_body1[3*n +2] /= mag1;
					
				}
				else{
					nds->normals_body1[3*n +0] = 0.0;
					nds->normals_body1[3*n +1] = 0.0;
					nds->normals_body1[3*n +2] = 0.0;
				}
				if(mag2>1e-8)
				{
					nds->normals_body2[3*n +0] /= mag2;
					nds->normals_body2[3*n +1] /= mag2;
					nds->normals_body2[3*n +2] /= mag2;
				}
				else{
					nds->normals_body2[3*n +0] = 0.0;
					nds->normals_body2[3*n +1] = 0.0;
					nds->normals_body2[3*n +2] = 0.0;
				}
					nds->normals_body2[3*n +0] = -nds->normals_body1[3*n+0];
					nds->normals_body2[3*n +1] = -nds->normals_body1[3*n+1];
					nds->normals_body2[3*n +2] = -nds->normals_body1[3*n+2];

				// Contact Detection
				if(nds->bodyID[2*n+1]!=-1) // only perform when two bodies are detected at a node
				{
				
					
						double approach1 = 0.0;
						double approach2 = 0.0;
						double relVel_body1[3];
						double relVel_body2[3];
							// Find out if they are approaching and correct
						if(nds->mass_body1[n]!=0.0)
						{
							
							relVel_body1[0] = (nds->momentum_body1[3*n+0]/nds->mass_body1[n]) - (nds->velocity[3*n+0]) 	;
							relVel_body1[1] = (nds->momentum_body1[3*n+1]/nds->mass_body1[n]) - (nds->velocity[3*n+1]) 	;
							relVel_body1[2] = (nds->momentum_body1[3*n+2]/nds->mass_body1[n]) - (nds->velocity[3*n+2]) 	;
							approach1 = (relVel_body1[0]*nds->normals_body1[3*n+0] + relVel_body1[1]*nds->normals_body1[3*n+1] + relVel_body1[2]*nds->normals_body1[3*n+2]);
						}
						if(nds->mass_body2[n]!=0.0)
						{
							
							relVel_body2[0] = (nds->momentum_body2[3*n+0]/nds->mass_body2[n]) - (nds->velocity[3*n+0]) 	;
							relVel_body2[1] = (nds->momentum_body2[3*n+1]/nds->mass_body2[n]) - (nds->velocity[3*n+1]) 	;
							relVel_body2[2] = (nds->momentum_body2[3*n+2]/nds->mass_body2[n]) - (nds->velocity[3*n+2]) 	;
							approach2 = (relVel_body2[0]*nds->normals_body2[3*n+0] + relVel_body2[1]*nds->normals_body2[3*n+1] + relVel_body2[2]*nds->normals_body2[3*n+2]);
						}
						bool isApproach = ((approach1 >0.01) || (approach2 > 0.01));
						if(isApproach) // two boddies are moving toward eachother
						{	
							if(nds->hasRigid[n*2+0]) // body 1 is rigid
							{
								if(nds->mass_body1[n]!=0.0){
								nds->velocity_body1[3*n+0]=nds->momentum_body1[3*n+0]/nds->mass_body1[n];
								nds->velocity_body1[3*n+1]=nds->momentum_body1[3*n+1]/nds->mass_body1[n];
								nds->velocity_body1[3*n+2]=nds->momentum_body1[3*n+2]/nds->mass_body1[n];
								if(nds->mass_body2[n]!=0.0)
								{
								nds->velocity_body2[3*n+0]=nds->momentum_body2[3*n+0]/nds->mass_body2[n] - (nds->momentum_body2[3*n+0] - nds->mass_body2[n]*nds->velocity_body1[3*n+0])/nds->mass_body2[n];
								nds->velocity_body2[3*n+1]=nds->momentum_body2[3*n+1]/nds->mass_body2[n] - (nds->momentum_body2[3*n+1] - nds->mass_body2[n]*nds->velocity_body1[3*n+1])/nds->mass_body2[n];
								nds->velocity_body2[3*n+2]=nds->momentum_body2[3*n+2]/nds->mass_body2[n] - (nds->momentum_body2[3*n+2] - nds->mass_body2[n]*nds->velocity_body1[3*n+2])/nds->mass_body2[n];
								}
								}
							}
							else if (nds->hasRigid[n*2+1]) // body 2 is rigid
							{
								if(nds->mass_body2[n]!=0.0){
								nds->velocity_body2[3*n+0]=nds->momentum_body2[3*n+0]/nds->mass_body2[n];
								nds->velocity_body2[3*n+1]=nds->momentum_body2[3*n+1]/nds->mass_body2[n];
								nds->velocity_body2[3*n+2]=nds->momentum_body2[3*n+2]/nds->mass_body2[n];
								if(nds->mass_body1[n]!=0.0)
								{
								nds->velocity_body1[3*n+0]=nds->momentum_body1[3*n+0]/nds->mass_body1[n] - (nds->momentum_body1[3*n+0] - nds->mass_body1[n]*nds->velocity_body2[3*n+0])/nds->mass_body1[n];
								nds->velocity_body1[3*n+1]=nds->momentum_body1[3*n+1]/nds->mass_body1[n] - (nds->momentum_body1[3*n+1] - nds->mass_body1[n]*nds->velocity_body2[3*n+1])/nds->mass_body1[n];
								nds->velocity_body1[3*n+2]=nds->momentum_body1[3*n+2]/nds->mass_body1[n] - (nds->momentum_body1[3*n+2] - nds->mass_body1[n]*nds->velocity_body2[3*n+2])/nds->mass_body1[n];
								}
								}
							}
							else{ // two non-rigid bodies approaching

							// BODY 1 
							// Correct velocity to avoid interpenetration
							if(nds->mass_body1[n]!=0.0){
								
								double contact_normals[3]; // normals for velocity correction
								double friction_normals[3];// normals for force corrections
								double mu_prime= 0.0; // effective coefficient of friction
								calculateColoumbFrictionContactNormals(contact_normals,friction_normals,mu_prime,relVel_body1,nds,mpmParameters,n,1);
								nds->velocity_body1[3*n+0] = (nds->momentum_body1[3*n+0]/nds->mass_body1[n]) - approach1*contact_normals[0] ; 
								nds->velocity_body1[3*n+1] = (nds->momentum_body1[3*n+1]/nds->mass_body1[n]) - approach1*contact_normals[1] ;
								nds->velocity_body1[3*n+2] = (nds->momentum_body1[3*n+2]/nds->mass_body1[n]) - approach1*contact_normals[2] ;
								nds->acceleration_body1[3*n+0] = (nds->force_body1[3*n+0]/nds->mass_body1[n]) - (approach1*contact_normals[0])/dt;// + (approach1*mu_prime*friction_normals[0])/dt;
								nds->acceleration_body1[3*n+1] = (nds->force_body1[3*n+1]/nds->mass_body1[n]) - (approach1*contact_normals[1])/dt;// + (approach1*mu_prime*friction_normals[1])/dt;
								nds->acceleration_body1[3*n+2] = (nds->force_body1[3*n+2]/nds->mass_body1[n]) - (approach1*contact_normals[2])/dt;// + (approach1*mu_prime*friction_normals[2])/dt;
								// nds->acceleration_body1[3*n+0] = (approach1*contact_normals[0])/dt -(approach1*mu_prime*friction_normals[0])/dt;
								// nds->acceleration_body1[3*n+1] = (approach1*contact_normals[1])/dt -(approach1*mu_prime*friction_normals[1])/dt;
								// nds->acceleration_body1[3*n+2] = (approach1*contact_normals[2])/dt -(approach1*mu_prime*friction_normals[2])/dt;
							
							if(step==1)
							{
								nds->velocity_body1[3*n+0] += nds->acceleration_body1[3*n+0]*0.5*dt; 
								nds->velocity_body1[3*n+1] += nds->acceleration_body1[3*n+1]*0.5*dt; 
								nds->velocity_body1[3*n+2] += nds->acceleration_body1[3*n+2]*0.5*dt; 
							}
							else{
								nds->velocity_body1[3*n+0] += nds->acceleration_body1[3*n+0]*dt; 
								nds->velocity_body1[3*n+1] += nds->acceleration_body1[3*n+1]*dt; 
								nds->velocity_body1[3*n+2] += nds->acceleration_body1[3*n+2]*dt; 
							}
							}
							// BODY 2
							// Correct velocity to avoid interpenetration
							if(nds->mass_body2[n]!=0.0)
							{
								double contact_normals[3]; // normals for velocity correction
								double friction_normals[3];// normals for force corrections
								double mu_prime=0.0; // effective coefficient of friction
								calculateColoumbFrictionContactNormals(contact_normals,friction_normals,mu_prime,relVel_body2,nds,mpmParameters,n,2);
								nds->velocity_body2[3*n+0] = (nds->momentum_body2[3*n+0]/nds->mass_body2[n]) - approach2*contact_normals[0] ; 
								nds->velocity_body2[3*n+1] = (nds->momentum_body2[3*n+1]/nds->mass_body2[n]) - approach2*contact_normals[1] ;
								nds->velocity_body2[3*n+2] = (nds->momentum_body2[3*n+2]/nds->mass_body2[n]) - approach2*contact_normals[2] ;
								nds->acceleration_body2[3*n+0] = (nds->force_body2[3*n+0]/nds->mass_body2[n]) - (approach2*contact_normals[0])/dt;// + (approach2*mu_prime*friction_normals[0])/dt;
								nds->acceleration_body2[3*n+1] = (nds->force_body2[3*n+1]/nds->mass_body2[n]) - (approach2*contact_normals[1])/dt;// + (approach2*mu_prime*friction_normals[1])/dt;
								nds->acceleration_body2[3*n+2] = (nds->force_body2[3*n+2]/nds->mass_body2[n]) - (approach2*contact_normals[2])/dt;// + (approach2*mu_prime*friction_normals[2])/dt;
								// nds->acceleration_body2[3*n+0] = (approach2*contact_normals[0])/dt -(approach2*mu_prime*friction_normals[0])/dt;
								// nds->acceleration_body2[3*n+1] = (approach2*contact_normals[1])/dt -(approach2*mu_prime*friction_normals[1])/dt;
								// nds->acceleration_body2[3*n+2] = (approach2*contact_normals[2])/dt -(approach2*mu_prime*friction_normals[2])/dt;
								
								if(step==1)
								{
									nds->velocity_body2[3*n+0] += nds->acceleration_body2[3*n+0]*0.5*dt; 
									nds->velocity_body2[3*n+1] += nds->acceleration_body2[3*n+1]*0.5*dt; 
									nds->velocity_body2[3*n+2] += nds->acceleration_body2[3*n+2]*0.5*dt; 
								}
								else
								{
									nds->velocity_body2[3*n+0] += nds->acceleration_body2[3*n+0]*dt; 
									nds->velocity_body2[3*n+1] += nds->acceleration_body2[3*n+1]*dt; 
									nds->velocity_body2[3*n+2] += nds->acceleration_body2[3*n+2]*dt; 
								}
							}
						}
						}
						else
						{ // Bodies arent approaching, so they can be updated separately
							// BODY 1 
							if(nds->mass_body1[n]!=0.0)
							{
								nds->velocity_body1[3*n+0] = (nds->momentum_body1[3*n+0]/nds->mass_body1[n]) ; 
								nds->velocity_body1[3*n+1] = (nds->momentum_body1[3*n+1]/nds->mass_body1[n]);
								nds->velocity_body1[3*n+2] = (nds->momentum_body1[3*n+2]/nds->mass_body1[n]) ;
								nds->acceleration_body1[3*n+0] = (nds->force_body1[3*n+0]/nds->mass_body1[n]) ;
								nds->acceleration_body1[3*n+1] = (nds->force_body1[3*n+1]/nds->mass_body1[n]);
								nds->acceleration_body1[3*n+2] = (nds->force_body1[3*n+2]/nds->mass_body1[n]) ;
								if(step==1)
								{
									nds->velocity_body1[3*n+0] += nds->acceleration_body1[3*n+0]*0.5*dt; 
									nds->velocity_body1[3*n+1] += nds->acceleration_body1[3*n+1]*0.5*dt; 
									nds->velocity_body1[3*n+2] += nds->acceleration_body1[3*n+2]*0.5*dt; 
								}
								else
								{
									nds->velocity_body1[3*n+0] += nds->acceleration_body1[3*n+0]*dt; 
									nds->velocity_body1[3*n+1] += nds->acceleration_body1[3*n+1]*dt; 
									nds->velocity_body1[3*n+2] += nds->acceleration_body1[3*n+2]*dt; 
								}
							}
									// BODY 2
							if(nds->mass_body2[n]!=0.0)
							{
								nds->velocity_body2[3*n+0] = (nds->momentum_body2[3*n+0]/nds->mass_body2[n]) ; 
								nds->velocity_body2[3*n+1] = (nds->momentum_body2[3*n+1]/nds->mass_body2[n])  ;
								nds->velocity_body2[3*n+2] = (nds->momentum_body2[3*n+2]/nds->mass_body2[n]) ;
								nds->acceleration_body2[3*n+0] = (nds->force_body2[3*n+0]/nds->mass_body2[n]);
								nds->acceleration_body2[3*n+1] = (nds->force_body2[3*n+1]/nds->mass_body2[n]);
								nds->acceleration_body2[3*n+2] = (nds->force_body2[3*n+2]/nds->mass_body2[n]);
									
								if(step==1)
								{
									nds->velocity_body2[3*n+0] += nds->acceleration_body2[3*n+0]*0.5*dt; 
									nds->velocity_body2[3*n+1] += nds->acceleration_body2[3*n+1]*0.5*dt; 
									nds->velocity_body2[3*n+2] += nds->acceleration_body2[3*n+2]*0.5*dt; 
								}
								else
								{
									nds->velocity_body2[3*n+0] += nds->acceleration_body2[3*n+0]*dt; 
									nds->velocity_body2[3*n+1] += nds->acceleration_body2[3*n+1]*dt; 
									nds->velocity_body2[3*n+2] += nds->acceleration_body2[3*n+2]*dt; 
								}
							}
						}
				}
				else // only one body detected at this node
				{
					if(nds->mass_body1[n]!=0.0){
						nds->velocity_body1[3*n+0] = (nds->momentum_body1[3*n+0]/nds->mass_body1[n]) ; 
						nds->velocity_body1[3*n+1] = (nds->momentum_body1[3*n+1]/nds->mass_body1[n]);
						nds->velocity_body1[3*n+2] = (nds->momentum_body1[3*n+2]/nds->mass_body1[n]) ;
						nds->acceleration_body1[3*n+0] = (nds->force_body1[3*n+0]/nds->mass_body1[n]) ;
						nds->acceleration_body1[3*n+1] = (nds->force_body1[3*n+1]/nds->mass_body1[n]);
						nds->acceleration_body1[3*n+2] = (nds->force_body1[3*n+2]/nds->mass_body1[n]) ;

						if(step==1)
						{
							nds->velocity_body1[3*n+0] += nds->acceleration_body1[3*n+0]*0.5*dt; 
							nds->velocity_body1[3*n+1] += nds->acceleration_body1[3*n+1]*0.5*dt; 
							nds->velocity_body1[3*n+2] += nds->acceleration_body1[3*n+2]*0.5*dt; 
						}
						else{
							nds->velocity_body1[3*n+0] += nds->acceleration_body1[3*n+0]*dt; 
							nds->velocity_body1[3*n+1] += nds->acceleration_body1[3*n+1]*dt; 
							nds->velocity_body1[3*n+2] += nds->acceleration_body1[3*n+2]*dt; 
						}
					}
				
				}
		}
	};
	
CPUGPU void calculateColoumbFrictionContactNormals(double* contactNormals,double* friction_normals,double mu_prime,double* relVel_body, Nodes* nds,Parameters* parameters, int n,int body){
		contactNormals[0] = 0.0;
		contactNormals[1] = 0.0;
		contactNormals[2] = 0.0;
		double mu = parameters->mu;
		double omega_mag;
		double omega[3];
		double omegan[3];
		double mu0;
	if(body==1)
	{
		omega[0] = relVel_body[1]*nds->normals_body1[3*n+2] - relVel_body[2]*nds->normals_body1[3*n+1];
		omega[1] = relVel_body[0]*nds->normals_body1[3*n+2] - relVel_body[2]*nds->normals_body1[3*n+0];
		omega[2] = relVel_body[0]*nds->normals_body1[3*n+1] - relVel_body[1]*nds->normals_body1[3*n+0];
		omega_mag = sqrt(pow(omega[0],2.)+pow(omega[1],2.)+pow(omega[2],2.)) + 0.001;
		omega[0] /= omega_mag;
		omega[1] /= omega_mag;
		omega[2] /= omega_mag;
		mu0 = omega_mag/(relVel_body[0]*nds->normals_body1[3*n+0] + relVel_body[1]*nds->normals_body1[3*n+1] + relVel_body[2]*nds->normals_body1[3*n+2]) +0.001;
		omegan[0] = omega[2]*nds->normals_body1[3*n+1] - omega[1]*nds->normals_body1[3*n+2];
		omegan[1] = omega[2]*nds->normals_body1[3*n+0] - omega[0]*nds->normals_body1[3*n+2];
		omegan[2] = omega[1]*nds->normals_body1[3*n+0] - omega[0]*nds->normals_body1[3*n+1];
		friction_normals[0] = omegan[0];
		friction_normals[1] = omegan[1];
		friction_normals[2] = omegan[2];
		mu = std::min(mu,mu0);
		mu_prime = mu;
		omegan[0]*=mu;
		omegan[1]*=mu;
		omegan[2]*=mu;

		contactNormals[0] = nds->normals_body1[3*n+0]+omegan[0];
		contactNormals[1] = nds->normals_body1[3*n+1]+omegan[1];
		contactNormals[2] = nds->normals_body1[3*n+2]+omegan[2];
	}
	if(body==2)
	{
		omega[0] = relVel_body[1]*nds->normals_body2[3*n+2] - relVel_body[2]*nds->normals_body2[3*n+1];
		omega[1] = relVel_body[0]*nds->normals_body2[3*n+2] - relVel_body[2]*nds->normals_body2[3*n+0];
		omega[2] = relVel_body[0]*nds->normals_body2[3*n+1] - relVel_body[1]*nds->normals_body2[3*n+0];
		omega_mag = sqrt(pow(omega[0],2.)+pow(omega[1],2.)+pow(omega[2],2.)) + 0.001;
		omega[0] /= omega_mag;
		omega[1] /= omega_mag;
		omega[2] /= omega_mag;
		mu0 = omega_mag/(relVel_body[0]*nds->normals_body2[3*n+0] + relVel_body[1]*nds->normals_body2[3*n+1] + relVel_body[2]*nds->normals_body2[3*n+2]) + 0.001;
		omegan[0] = omega[2]*nds->normals_body2[3*n+1] - omega[1]*nds->normals_body2[3*n+2];
		omegan[1] = omega[2]*nds->normals_body2[3*n+0] - omega[0]*nds->normals_body2[3*n+2];
		omegan[2] = omega[1]*nds->normals_body2[3*n+0] - omega[0]*nds->normals_body2[3*n+1];
		mu = std::min(mu,mu0);
		friction_normals[0] = omegan[0];
		friction_normals[1] = omegan[1];
		friction_normals[2] = omegan[2];
		mu_prime = mu;
		omegan[0]*=mu;
		omegan[1]*=mu;
		omegan[2]*=mu;

		contactNormals[0] = nds->normals_body2[3*n+0]+omegan[0];
		contactNormals[1] = nds->normals_body2[3*n+1]+omegan[1];
		contactNormals[2] = nds->normals_body2[3*n+2]+omegan[2];
	}
};

CPUGPU void collectPeriodicContributions(Parameters* mpmParameters,Nodes *nds){

if(mpmParameters->is3D)
{
int Nx      = nds->Nx;
int Ny      = nds->Ny;
int Nz      = nds->Nz;
int numNode = 3;

// X- Direction 
for (int i = 0; i < numNode; i++)
{
	#pragma omp parallel for
	for (int j = 0; j < Ny; j++)
	{
		for (int k = 0; k < Nz; k++)
		{
			// X - direction boundaries
			int node     = i+j*Nx+k*Nx*Ny;
			int antinode = (Nx-numNode)+i+j*Nx+k*Nx*Ny;
			double sum   = 0.0;
			
			sum                                           = nds->momentum[3*node+0] + nds->momentum[3*antinode + 0 ];
			nds->momentum[3*node+0]                       = sum ;
			nds->momentum[3*antinode+ 0 ]                 = sum;

			sum                                           = nds->momentum[3*node+1] + nds->momentum[3*antinode + 1];
			nds->momentum[3*node+1]                       = sum ;
			nds->momentum[3*antinode+ 1]                  = sum;

			sum                                           = nds->momentum[3*node+2] + nds->momentum[3*antinode + 2];
			nds->momentum[3*node+2]                       = sum ;
			nds->momentum[3*antinode+ 2]                  = sum;
			sum                                           = nds->momentum_body1[3*node+0] + nds->momentum_body1[3*antinode + 0 ];
			nds->momentum_body1[3*node+0]                       = sum ;
			nds->momentum_body1[3*antinode+ 0 ]                 = sum;

			sum                                           = nds->momentum_body1[3*node+1] + nds->momentum_body1[3*antinode + 1];
			nds->momentum_body1[3*node+1]                       = sum ;
			nds->momentum_body1[3*antinode+ 1]                  = sum;

			sum                                           = nds->momentum_body1[3*node+2] + nds->momentum_body1[3*antinode + 2];
			nds->momentum_body1[3*node+2]                       = sum ;
			nds->momentum_body1[3*antinode+ 2]                  = sum;
			sum                                           = nds->momentum_body2[3*node+0] + nds->momentum_body2[3*antinode + 0 ];
			nds->momentum_body2[3*node+0]                       = sum ;
			nds->momentum_body2[3*antinode+ 0 ]                 = sum;

			sum                                           = nds->momentum_body2[3*node+1] + nds->momentum_body2[3*antinode + 1];
			nds->momentum_body2[3*node+1]                       = sum ;
			nds->momentum_body2[3*antinode+ 1]                  = sum;

			sum                                           = nds->momentum_body2[3*node+2] + nds->momentum_body2[3*antinode + 2];
			nds->momentum_body2[3*node+2]                       = sum ;
			nds->momentum_body2[3*antinode+ 2]                  = sum;

			sum                                           = nds->mass[node] + nds->mass[antinode];
			nds->mass[antinode]                           = sum;
			nds->mass[node]                               = sum; 
			sum                                           = nds->mass_body1[node] + nds->mass_body1[antinode];
			nds->mass_body1[antinode]                           = sum;
			nds->mass_body1[node]                               = sum; 
			sum                                           = nds->mass_body2[node] + nds->mass_body2[antinode];
			nds->mass_body2[antinode]                           = sum;
			nds->mass_body2[node]                               = sum; 

			sum                                           = nds->force[3*node+0] + nds->force[3*antinode + 0];
			nds->force[3*antinode + 0]                    = sum;
			nds->force[3*node+0]                          = sum; 

			sum                                           = nds->force[3*node+1] + nds->force[3*antinode + 1];
			nds->force[3*antinode + 1]                    = sum;
			nds->force[3*node+1]                          = sum; 

			sum                                           = nds->force[3*node+2] + nds->force[3*antinode + 2];
			nds->force[3*antinode + 2 ]                   = sum;
			nds->force[3*node+2]                          = sum; 
			sum                                           = nds->force_body1[3*node+0] + nds->force_body1[3*antinode + 0];
			nds->force_body1[3*antinode + 0]                    = sum;
			nds->force_body1[3*node+0]                          = sum; 

			sum                                           = nds->force_body1[3*node+1] + nds->force_body1[3*antinode + 1];
			nds->force_body1[3*antinode + 1]                    = sum;
			nds->force_body1[3*node+1]                          = sum; 

			sum                                           = nds->force_body1[3*node+2] + nds->force_body1[3*antinode + 2];
			nds->force_body1[3*antinode + 2 ]                   = sum;
			nds->force_body1[3*node+2]                          = sum; 
			sum                                           = nds->force_body2[3*node+0] + nds->force_body2[3*antinode + 0];
			nds->force_body2[3*antinode + 0]                    = sum;
			nds->force_body2[3*node+0]                          = sum; 

			sum                                           = nds->force_body2[3*node+1] + nds->force_body2[3*antinode + 1];
			nds->force_body2[3*antinode + 1]                    = sum;
			nds->force_body2[3*node+1]                          = sum; 

			sum                                           = nds->force_body2[3*node+2] + nds->force_body2[3*antinode + 2];
			nds->force_body2[3*antinode + 2 ]                   = sum;
			nds->force_body2[3*node+2]                          = sum; 
		}
	}
}

// Y - Direction 
for (int j = 0; j < numNode; j++)
{
	#pragma omp parallel for
	for (int k = 0; k < Nz; k++)
	{

		for (int i = 0; i < Nx; i++)
		{
			// Y - direction boundaries
			int node     = i+j*Nx+k*Nx*Ny;
			int antinode = i+(j+Ny-numNode)*Nx+k*Nx*Ny;
			double sum   = 0.0;
			
			
			sum                                           = nds->momentum[3*node+0] + nds->momentum[3*antinode + 0 ];
			nds->momentum[3*node+0]                       = sum ;
			nds->momentum[3*antinode+ 0 ]                 = sum;

			sum                                           = nds->momentum[3*node+1] + nds->momentum[3*antinode + 1];
			nds->momentum[3*node+1]                       = sum ;
			nds->momentum[3*antinode+ 1]                  = sum;

			sum                                           = nds->momentum[3*node+2] + nds->momentum[3*antinode + 2];
			nds->momentum[3*node+2]                       = sum ;
			nds->momentum[3*antinode+ 2]                  = sum;
			sum                                           = nds->momentum_body1[3*node+0] + nds->momentum_body1[3*antinode + 0 ];
			nds->momentum_body1[3*node+0]                       = sum ;
			nds->momentum_body1[3*antinode+ 0 ]                 = sum;

			sum                                           = nds->momentum_body1[3*node+1] + nds->momentum_body1[3*antinode + 1];
			nds->momentum_body1[3*node+1]                       = sum ;
			nds->momentum_body1[3*antinode+ 1]                  = sum;

			sum                                           = nds->momentum_body1[3*node+2] + nds->momentum_body1[3*antinode + 2];
			nds->momentum_body1[3*node+2]                       = sum ;
			nds->momentum_body1[3*antinode+ 2]                  = sum;
			sum                                           = nds->momentum_body2[3*node+0] + nds->momentum_body2[3*antinode + 0 ];
			nds->momentum_body2[3*node+0]                       = sum ;
			nds->momentum_body2[3*antinode+ 0 ]                 = sum;

			sum                                           = nds->momentum_body2[3*node+1] + nds->momentum_body2[3*antinode + 1];
			nds->momentum_body2[3*node+1]                       = sum ;
			nds->momentum_body2[3*antinode+ 1]                  = sum;

			sum                                           = nds->momentum_body2[3*node+2] + nds->momentum_body2[3*antinode + 2];
			nds->momentum_body2[3*node+2]                       = sum ;
			nds->momentum_body2[3*antinode+ 2]                  = sum;

			sum                                           = nds->mass[node] + nds->mass[antinode];
			nds->mass[antinode]                           = sum;
			nds->mass[node]                               = sum; 
			sum                                           = nds->mass_body1[node] + nds->mass_body1[antinode];
			nds->mass_body1[antinode]                           = sum;
			nds->mass_body1[node]                               = sum; 
			sum                                           = nds->mass_body2[node] + nds->mass_body2[antinode];
			nds->mass_body2[antinode]                           = sum;
			nds->mass_body2[node]                               = sum; 

			sum                                           = nds->force[3*node+0] + nds->force[3*antinode + 0];
			nds->force[3*antinode + 0]                    = sum;
			nds->force[3*node+0]                          = sum; 

			sum                                           = nds->force[3*node+1] + nds->force[3*antinode + 1];
			nds->force[3*antinode + 1]                    = sum;
			nds->force[3*node+1]                          = sum; 

			sum                                           = nds->force[3*node+2] + nds->force[3*antinode + 2];
			nds->force[3*antinode + 2 ]                   = sum;
			nds->force[3*node+2]                          = sum; 
			sum                                           = nds->force_body1[3*node+0] + nds->force_body1[3*antinode + 0];
			nds->force_body1[3*antinode + 0]                    = sum;
			nds->force_body1[3*node+0]                          = sum; 

			sum                                           = nds->force_body1[3*node+1] + nds->force_body1[3*antinode + 1];
			nds->force_body1[3*antinode + 1]                    = sum;
			nds->force_body1[3*node+1]                          = sum; 

			sum                                           = nds->force_body1[3*node+2] + nds->force_body1[3*antinode + 2];
			nds->force_body1[3*antinode + 2 ]                   = sum;
			nds->force_body1[3*node+2]                          = sum; 
			sum                                           = nds->force_body2[3*node+0] + nds->force_body2[3*antinode + 0];
			nds->force_body2[3*antinode + 0]                    = sum;
			nds->force_body2[3*node+0]                          = sum; 

			sum                                           = nds->force_body2[3*node+1] + nds->force_body2[3*antinode + 1];
			nds->force_body2[3*antinode + 1]                    = sum;
			nds->force_body2[3*node+1]                          = sum; 

			sum                                           = nds->force_body2[3*node+2] + nds->force_body2[3*antinode + 2];
			nds->force_body2[3*antinode + 2 ]                   = sum;
			nds->force_body2[3*node+2]                          = sum; 
		}		
	}
}

// Z - Direction 
for (int k = 0; k <numNode; k++)
{
	#pragma omp parallel for
	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			// Z - direction boundaries
			int node     = i+j*Nx+k*Nx*Ny;
			int antinode = i+j*Nx+(k+Nz-numNode)*Nx*Ny;
			double sum   = 0.0;
			
			
			sum                                           = nds->momentum[3*node+0] + nds->momentum[3*antinode + 0 ];
			nds->momentum[3*node+0]                       = sum ;
			nds->momentum[3*antinode+ 0 ]                 = sum;

			sum                                           = nds->momentum[3*node+1] + nds->momentum[3*antinode + 1];
			nds->momentum[3*node+1]                       = sum ;
			nds->momentum[3*antinode+ 1]                  = sum;

			sum                                           = nds->momentum[3*node+2] + nds->momentum[3*antinode + 2];
			nds->momentum[3*node+2]                       = sum ;
			nds->momentum[3*antinode+ 2]                  = sum;
			sum                                           = nds->momentum_body1[3*node+0] + nds->momentum_body1[3*antinode + 0 ];
			nds->momentum_body1[3*node+0]                       = sum ;
			nds->momentum_body1[3*antinode+ 0 ]                 = sum;

			sum                                           = nds->momentum_body1[3*node+1] + nds->momentum_body1[3*antinode + 1];
			nds->momentum_body1[3*node+1]                       = sum ;
			nds->momentum_body1[3*antinode+ 1]                  = sum;

			sum                                           = nds->momentum_body1[3*node+2] + nds->momentum_body1[3*antinode + 2];
			nds->momentum_body1[3*node+2]                       = sum ;
			nds->momentum_body1[3*antinode+ 2]                  = sum;
			sum                                           = nds->momentum_body2[3*node+0] + nds->momentum_body2[3*antinode + 0 ];
			nds->momentum_body2[3*node+0]                       = sum ;
			nds->momentum_body2[3*antinode+ 0 ]                 = sum;

			sum                                           = nds->momentum_body2[3*node+1] + nds->momentum_body2[3*antinode + 1];
			nds->momentum_body2[3*node+1]                       = sum ;
			nds->momentum_body2[3*antinode+ 1]                  = sum;

			sum                                           = nds->momentum_body2[3*node+2] + nds->momentum_body2[3*antinode + 2];
			nds->momentum_body2[3*node+2]                       = sum ;
			nds->momentum_body2[3*antinode+ 2]                  = sum;

			sum                                           = nds->mass[node] + nds->mass[antinode];
			nds->mass[antinode]                           = sum;
			nds->mass[node]                               = sum; 
			sum                                           = nds->mass_body1[node] + nds->mass_body1[antinode];
			nds->mass_body1[antinode]                           = sum;
			nds->mass_body1[node]                               = sum; 
			sum                                           = nds->mass_body2[node] + nds->mass_body2[antinode];
			nds->mass_body2[antinode]                           = sum;
			nds->mass_body2[node]                               = sum; 

			sum                                           = nds->force[3*node+0] + nds->force[3*antinode + 0];
			nds->force[3*antinode + 0]                    = sum;
			nds->force[3*node+0]                          = sum; 

			sum                                           = nds->force[3*node+1] + nds->force[3*antinode + 1];
			nds->force[3*antinode + 1]                    = sum;
			nds->force[3*node+1]                          = sum; 

			sum                                           = nds->force[3*node+2] + nds->force[3*antinode + 2];
			nds->force[3*antinode + 2 ]                   = sum;
			nds->force[3*node+2]                          = sum; 
			sum                                           = nds->force_body1[3*node+0] + nds->force_body1[3*antinode + 0];
			nds->force_body1[3*antinode + 0]                    = sum;
			nds->force_body1[3*node+0]                          = sum; 

			sum                                           = nds->force_body1[3*node+1] + nds->force_body1[3*antinode + 1];
			nds->force_body1[3*antinode + 1]                    = sum;
			nds->force_body1[3*node+1]                          = sum; 

			sum                                           = nds->force_body1[3*node+2] + nds->force_body1[3*antinode + 2];
			nds->force_body1[3*antinode + 2 ]                   = sum;
			nds->force_body1[3*node+2]                          = sum; 
			sum                                           = nds->force_body2[3*node+0] + nds->force_body2[3*antinode + 0];
			nds->force_body2[3*antinode + 0]                    = sum;
			nds->force_body2[3*node+0]                          = sum; 

			sum                                           = nds->force_body2[3*node+1] + nds->force_body2[3*antinode + 1];
			nds->force_body2[3*antinode + 1]                    = sum;
			nds->force_body2[3*node+1]                          = sum; 

			sum                                           = nds->force_body2[3*node+2] + nds->force_body2[3*antinode + 2];
			nds->force_body2[3*antinode + 2 ]                   = sum;
			nds->force_body2[3*node+2]                          = sum; 
		}				
	}
}
}
else{
	int Nx      = nds->Nx;
	int Ny      = nds->Ny;
	int numNode = 3;
	
	// X- Direction 
	for (int i = 0; i < numNode; i++)
	{
		#pragma omp parallel for
		for (int j = 0; j < Ny; j++)
		{
				// X - direction boundaries
				int node     = i+j*Nx;
				int antinode = (Nx-numNode)+i+j*Nx;
				double sum   = 0.0;
				
				sum                                           = nds->momentum[3*node+0] + nds->momentum[3*antinode + 0 ];
				nds->momentum[3*node+0]                       = sum ;
				nds->momentum[3*antinode+ 0 ]                 = sum;
	
				sum                                           = nds->momentum[3*node+1] + nds->momentum[3*antinode + 1];
				nds->momentum[3*node+1]                       = sum ;
				nds->momentum[3*antinode+ 1]                  = sum;
	
				sum                                           = nds->momentum[3*node+2] + nds->momentum[3*antinode + 2];
				nds->momentum[3*node+2]                       = sum ;
				nds->momentum[3*antinode+ 2]                  = sum;
				sum                                           = nds->momentum_body1[3*node+0] + nds->momentum_body1[3*antinode + 0 ];
				nds->momentum_body1[3*node+0]                       = sum ;
				nds->momentum_body1[3*antinode+ 0 ]                 = sum;
	
				sum                                           = nds->momentum_body1[3*node+1] + nds->momentum_body1[3*antinode + 1];
				nds->momentum_body1[3*node+1]                       = sum ;
				nds->momentum_body1[3*antinode+ 1]                  = sum;
	
				sum                                           = nds->momentum_body1[3*node+2] + nds->momentum_body1[3*antinode + 2];
				nds->momentum_body1[3*node+2]                       = sum ;
				nds->momentum_body1[3*antinode+ 2]                  = sum;
				sum                                           = nds->momentum_body2[3*node+0] + nds->momentum_body2[3*antinode + 0 ];
				nds->momentum_body2[3*node+0]                       = sum ;
				nds->momentum_body2[3*antinode+ 0 ]                 = sum;
	
				sum                                           = nds->momentum_body2[3*node+1] + nds->momentum_body2[3*antinode + 1];
				nds->momentum_body2[3*node+1]                       = sum ;
				nds->momentum_body2[3*antinode+ 1]                  = sum;
	
				sum                                           = nds->momentum_body2[3*node+2] + nds->momentum_body2[3*antinode + 2];
				nds->momentum_body2[3*node+2]                       = sum ;
				nds->momentum_body2[3*antinode+ 2]                  = sum;
	
				sum                                           = nds->mass[node] + nds->mass[antinode];
				nds->mass[antinode]                           = sum;
				nds->mass[node]                               = sum; 
				sum                                           = nds->mass_body1[node] + nds->mass_body1[antinode];
				nds->mass_body1[antinode]                           = sum;
				nds->mass_body1[node]                               = sum; 
				sum                                           = nds->mass_body2[node] + nds->mass_body2[antinode];
				nds->mass_body2[antinode]                           = sum;
				nds->mass_body2[node]                               = sum; 
	
				sum                                           = nds->force[3*node+0] + nds->force[3*antinode + 0];
				nds->force[3*antinode + 0]                    = sum;
				nds->force[3*node+0]                          = sum; 
	
				sum                                           = nds->force[3*node+1] + nds->force[3*antinode + 1];
				nds->force[3*antinode + 1]                    = sum;
				nds->force[3*node+1]                          = sum; 
	
				sum                                           = nds->force[3*node+2] + nds->force[3*antinode + 2];
				nds->force[3*antinode + 2 ]                   = sum;
				nds->force[3*node+2]                          = sum; 
				sum                                           = nds->force_body1[3*node+0] + nds->force_body1[3*antinode + 0];
				nds->force_body1[3*antinode + 0]                    = sum;
				nds->force_body1[3*node+0]                          = sum; 
	
				sum                                           = nds->force_body1[3*node+1] + nds->force_body1[3*antinode + 1];
				nds->force_body1[3*antinode + 1]                    = sum;
				nds->force_body1[3*node+1]                          = sum; 
	
				sum                                           = nds->force_body1[3*node+2] + nds->force_body1[3*antinode + 2];
				nds->force_body1[3*antinode + 2 ]                   = sum;
				nds->force_body1[3*node+2]                          = sum; 
				sum                                           = nds->force_body2[3*node+0] + nds->force_body2[3*antinode + 0];
				nds->force_body2[3*antinode + 0]                    = sum;
				nds->force_body2[3*node+0]                          = sum; 
	
				sum                                           = nds->force_body2[3*node+1] + nds->force_body2[3*antinode + 1];
				nds->force_body2[3*antinode + 1]                    = sum;
				nds->force_body2[3*node+1]                          = sum; 
	
				sum                                           = nds->force_body2[3*node+2] + nds->force_body2[3*antinode + 2];
				nds->force_body2[3*antinode + 2 ]                   = sum;
				nds->force_body2[3*node+2]                          = sum; 
			}
	}
	
	// Y - Direction 
	for (int j = 0; j < numNode; j++)
	{
		#pragma omp parallel for

			for (int i = 0; i < Nx; i++)
			{
				// Y - direction boundaries
				int node     = i+j*Nx;
				int antinode = i+(j+Ny-numNode)*Nx;
				double sum   = 0.0;
				
				
				sum                                           = nds->momentum[3*node+0] + nds->momentum[3*antinode + 0 ];
				nds->momentum[3*node+0]                       = sum ;
				nds->momentum[3*antinode+ 0 ]                 = sum;
	
				sum                                           = nds->momentum[3*node+1] + nds->momentum[3*antinode + 1];
				nds->momentum[3*node+1]                       = sum ;
				nds->momentum[3*antinode+ 1]                  = sum;
	
				sum                                           = nds->momentum[3*node+2] + nds->momentum[3*antinode + 2];
				nds->momentum[3*node+2]                       = sum ;
				nds->momentum[3*antinode+ 2]                  = sum;
				sum                                           = nds->momentum_body1[3*node+0] + nds->momentum_body1[3*antinode + 0 ];
				nds->momentum_body1[3*node+0]                       = sum ;
				nds->momentum_body1[3*antinode+ 0 ]                 = sum;
	
				sum                                           = nds->momentum_body1[3*node+1] + nds->momentum_body1[3*antinode + 1];
				nds->momentum_body1[3*node+1]                       = sum ;
				nds->momentum_body1[3*antinode+ 1]                  = sum;
	
				sum                                           = nds->momentum_body1[3*node+2] + nds->momentum_body1[3*antinode + 2];
				nds->momentum_body1[3*node+2]                       = sum ;
				nds->momentum_body1[3*antinode+ 2]                  = sum;

				sum                                           = nds->momentum_body2[3*node+0] + nds->momentum_body2[3*antinode + 0 ];
				nds->momentum_body2[3*node+0]                       = sum ;
				nds->momentum_body2[3*antinode+ 0 ]                 = sum;
	
				sum                                           = nds->momentum_body2[3*node+1] + nds->momentum_body2[3*antinode + 1];
				nds->momentum_body2[3*node+1]                       = sum ;
				nds->momentum_body2[3*antinode+ 1]                  = sum;
	
				sum                                           = nds->momentum_body2[3*node+2] + nds->momentum_body2[3*antinode + 2];
				nds->momentum_body2[3*node+2]                       = sum ;
				nds->momentum_body2[3*antinode+ 2]                  = sum;
	
				sum                                           = nds->mass[node] + nds->mass[antinode];
				nds->mass[antinode]                           = sum;
				nds->mass[node]                               = sum; 
				sum                                           = nds->mass_body1[node] + nds->mass_body1[antinode];
				nds->mass_body1[antinode]                           = sum;
				nds->mass_body1[node]                               = sum; 
				
				sum                                           = nds->mass_body2[node] + nds->mass_body2[antinode];
				nds->mass_body2[antinode]                           = sum;
				nds->mass_body2[node]                               = sum; 
	
				sum                                           = nds->force[3*node+0] + nds->force[3*antinode + 0];
				nds->force[3*antinode + 0]                    = sum;
				nds->force[3*node+0]                          = sum; 
	
				sum                                           = nds->force[3*node+1] + nds->force[3*antinode + 1];
				nds->force[3*antinode + 1]                    = sum;
				nds->force[3*node+1]                          = sum; 
	
				sum                                           = nds->force[3*node+2] + nds->force[3*antinode + 2];
				nds->force[3*antinode + 2 ]                   = sum;
				nds->force[3*node+2]                          = sum; 
				sum                                           = nds->force_body1[3*node+0] + nds->force_body1[3*antinode + 0];
				nds->force_body1[3*antinode + 0]                    = sum;
				nds->force_body1[3*node+0]                          = sum; 
	
				sum                                           = nds->force_body1[3*node+1] + nds->force_body1[3*antinode + 1];
				nds->force_body1[3*antinode + 1]                    = sum;
				nds->force_body1[3*node+1]                          = sum; 
	
				sum                                           = nds->force_body1[3*node+2] + nds->force_body1[3*antinode + 2];
				nds->force_body1[3*antinode + 2 ]                   = sum;
				nds->force_body1[3*node+2]                          = sum; 
				sum                                           = nds->force_body2[3*node+0] + nds->force_body2[3*antinode + 0];
				nds->force_body2[3*antinode + 0]                    = sum;
				nds->force_body2[3*node+0]                          = sum; 
	
				sum                                           = nds->force_body2[3*node+1] + nds->force_body2[3*antinode + 1];
				nds->force_body2[3*antinode + 1]                    = sum;
				nds->force_body2[3*node+1]                          = sum; 
	
				sum                                           = nds->force_body2[3*node+2] + nds->force_body2[3*antinode + 2];
				nds->force_body2[3*antinode + 2 ]                   = sum;
				nds->force_body2[3*node+2]                          = sum; 
			}		
	}
}
};

CPUGPU void resetNode(Nodes* nds,int n){
	//Reset node to previous state
		nds->velocity[n*3 +0]    = 0.0;
		nds->velocity[n*3 +1]    = 0.0;
		nds->velocity[n*3 +2]    = 0.0;

		nds->velocity_body1[n*3 + 0]  = 0.0;
		nds->velocity_body1[n*3 + 1]  = 0.0;
		nds->velocity_body1[n*3 + 2]  = 0.0;
		nds->velocity_body2[n*3 + 0]  = 0.0;
		nds->velocity_body2[n*3 + 1]  = 0.0;
		nds->velocity_body2[n*3 + 2]  = 0.0;

		nds->acceleration[n*3+0] = 0.0;
		nds->acceleration[n*3+1] = 0.0;
		nds->acceleration[n*3+2] = 0.0;

		nds->acceleration_body1[n*3+0] = 0.0;
		nds->acceleration_body1[n*3+1] = 0.0;
		nds->acceleration_body1[n*3+2] = 0.0;
		nds->acceleration_body2[n*3+0] = 0.0;
		nds->acceleration_body2[n*3+1] = 0.0;
		nds->acceleration_body2[n*3+2] = 0.0;

		nds->force[n*3+0]        = 0.0;
		nds->force[n*3+1]        = 0.0;
		nds->force[n*3+2]        = 0.0;
		nds->force_body1[n*3+0]        = 0.0;
		nds->force_body1[n*3+1]        = 0.0;
		nds->force_body1[n*3+2]        = 0.0;
		nds->force_body2[n*3+0]        = 0.0;
		nds->force_body2[n*3+1]        = 0.0;
		nds->force_body2[n*3+2]        = 0.0;

		nds->momentum[n*3+0]     = 0.0;
		nds->momentum[n*3+1]     = 0.0;
		nds->momentum[n*3+2]     = 0.0;
		nds->momentum_body1[n*3 + 0]     = 0.0;
		nds->momentum_body1[n*3 + 1]     = 0.0;
		nds->momentum_body1[n*3 + 2]     = 0.0;
		nds->momentum_body2[n*3 + 0]     = 0.0;
		nds->momentum_body2[n*3 + 1]     = 0.0;
		nds->momentum_body2[n*3 + 2]     = 0.0;

		nds->mass[n]             = 0.0;
		nds->mass_body1[n]         = 0.0;
		nds->mass_body2[n]         = 0.0;
		
		nds->temperature[n]      = 0.0;
		nds->temperatureRate[n]  = 0.0;
		nds->Tmass[n]            = 0.0;

		nds->bodyID[2*n + 0]          = -1;
		nds->bodyID[2*n + 1]          = -1;
		nds->hasRigid[2*n + 0]     = false;
		nds->hasRigid[2*n + 1]     = false;

		nds->normals_body1[n*3+0]      =0.0;
		nds->normals_body1[n*3+1]      =0.0;
		nds->normals_body1[n*3+2]      =0.0;
		nds->normals_body2[n*3+0]      =0.0;
		nds->normals_body2[n*3+1]      =0.0;
		nds->normals_body2[n*3+2]      =0.0;
};

CPUGPU void enforceGridBC(Nodes* nds,int phase, int n){
		int bc_type = nds->BCTYPE[n];
		switch(bc_type)
		{
			case BCTypes::grid_vx:
				bc_grid_vx(nds,n,phase);
				break;
			case BCTypes::grid_vy:
				bc_grid_vy(nds,n,phase);
				break;
			case BCTypes::grid_vz:
				bc_grid_vx_vz(nds,n,phase);
				break;
			case BCTypes::grid_temp:
				bc_grid_temp(nds,n,phase);
				break;
			case BCTypes::grid_tflux_x:
				bc_grid_tflux_x(nds,n,phase);
				break;
			case BCTypes::grid_tflux_y:
				bc_grid_tflux_y(nds,n,phase);
				break;
			case BCTypes::grid_tflux_z:
				bc_grid_tflux_z(nds,n,phase);
				break;
			case BCTypes::grid_fx:
				bc_grid_fx(nds,n,phase);
				break;
			case BCTypes::grid_fy:
				bc_grid_fy(nds,n,phase);
				break;
			case BCTypes::grid_fz:
				bc_grid_fz(nds,n,phase);
				break;
			case BCTypes::grid_vx_vy:
				bc_grid_vx_vy(nds,n,phase);
				break;
			case BCTypes::grid_vx_vz:
				bc_grid_vx_vz(nds,n,phase);
				break;
			case BCTypes::grid_vy_vz:
				bc_grid_vy_vz(nds,n,phase);
				break;
			case BCTypes::grid_fx_fy:
				bc_grid_fx_fy(nds,n,phase);
				break;
			case BCTypes::grid_fy_fz:
				bc_grid_fy_fz(nds,n,phase);
				break;
			case BCTypes::grid_fx_fz:
				bc_grid_fx_fz(nds,n,phase);
				break;
			case BCTypes::grid_fx_fy_fz:
				bc_grid_fx_fy_fz(nds,n,phase);
				break;
			case BCTypes::grid_vx_vy_vz:
				bc_grid_vx_vy_vz(nds,n,phase);
				break;
			case BCTypes::grid_vx_temp:
				bc_grid_vx_temp(nds,n,phase);
				break;
			case BCTypes::grid_vy_temp:
				bc_grid_vy_temp(nds,n,phase);
				break;
			case BCTypes::grid_vz_temp:
				bc_grid_vz_temp(nds,n,phase);
				break;
			case BCTypes::grid_fx_temp:
				bc_grid_fx_temp(nds,n,phase);
				break;
			case BCTypes::grid_fy_temp:
				bc_grid_fy_temp(nds,n,phase);
				break;
			case BCTypes::grid_fz_temp:
				bc_grid_fz_temp(nds,n,phase);
				break;
			case BCTypes::grid_vx_vy_temp:
				bc_grid_vx_vy_temp(nds,n,phase);
				break;
			case BCTypes::grid_vx_vz_temp:
				bc_grid_vx_vz_temp(nds,n,phase);
				break;
			case BCTypes::grid_vy_vz_temp:
				bc_grid_vy_vz_temp(nds,n,phase);
				break;
			case BCTypes::grid_fx_fy_temp:
				bc_grid_fx_fy_temp(nds,n,phase);
				break;
			case BCTypes::grid_fy_fz_temp:
				bc_grid_fy_fz_temp(nds,n,phase);
				break;
			case BCTypes::grid_fx_fz_temp:
				bc_grid_fx_fz_temp(nds,n,phase);
				break;
			case BCTypes::grid_fx_fy_fz_temp:
				bc_grid_fx_fy_fz_temp(nds,n,phase);
				break;
			case BCTypes::grid_vx_vy_vz_temp:
				bc_grid_vx_vy_vz_temp(nds,n,phase);
				break;
			case BCTypes::none:
				break;
		}
};

CPUGPU void bc_grid_vx(Nodes* nds, int n,int phase){
	nds->velocity_body1[3*n+0]     = nds->BCVAL[4*n+0];
	nds->acceleration_body1[3*n+0] = 0.0;
	nds->velocity_body2[3*n+0]     = nds->BCVAL[4*n+0];
	nds->acceleration_body2[3*n+0] = 0.0;
	nds->velocity[3*n+0]     = nds->BCVAL[4*n+0];
	nds->acceleration[3*n+0] = 0.0;
};
CPUGPU void bc_grid_vy(Nodes* nds, int n,int phase){
	nds->velocity_body1[3*n+1]     = nds->BCVAL[4*n+0];
	nds->acceleration_body1[3*n+1] = 0.0;
	nds->velocity_body2[3*n+1]     = nds->BCVAL[4*n+0];
	nds->acceleration_body2[3*n+1] = 0.0;
	nds->velocity[3*n+1]     = nds->BCVAL[4*n+0];
	nds->acceleration[3*n+1] = 0.0;
};
CPUGPU void bc_grid_vz(Nodes* nds, int n,int phase){
	nds->velocity_body1[3*n+2]     = nds->BCVAL[4*n+0];
	nds->acceleration_body1[3*n+2] = 0.0;
	nds->velocity_body2[3*n+2]     = nds->BCVAL[4*n+0];
	nds->acceleration_body2[3*n+2] = 0.0;
	nds->velocity[3*n+2]     = nds->BCVAL[4*n+0];
	nds->acceleration[3*n+2] = 0.0;
};
CPUGPU void bc_grid_temp(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_tflux_x(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_tflux_y(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_tflux_z(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_fx(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_fy(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_fz(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_vx_vy(Nodes* nds, int n,int phase){
	nds->velocity_body1[3*n+0]     = nds->BCVAL[4*n+0];
	nds->velocity_body1[3*n+1]     = nds->BCVAL[4*n+1];
	
	nds->acceleration_body1[3*n+0] = 0.0;
	nds->acceleration_body1[3*n+1] = 0.0;

	nds->velocity_body2[3*n+0]     = nds->BCVAL[4*n+0];
	nds->velocity_body2[3*n+1]     = nds->BCVAL[4*n+1];
	
	nds->acceleration_body2[3*n+0] = 0.0;
	nds->acceleration_body2[3*n+1] = 0.0;

	nds->velocity[3*n+0]     = nds->BCVAL[4*n+0];
	nds->velocity[3*n+1]     = nds->BCVAL[4*n+1];
	
	nds->acceleration[3*n+0] = 0.0;
	nds->acceleration[3*n+1] = 0.0;
	
};
CPUGPU void bc_grid_vx_vz(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_vy_vz(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_fx_fy(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_fy_fz(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_fx_fz(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_fx_fy_fz(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_vx_vy_vz(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_vx_temp(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_vy_temp(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_vz_temp(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_fx_temp(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_fy_temp(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_fz_temp(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_vx_vy_temp(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_vx_vz_temp(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_vy_vz_temp(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_fx_fy_temp(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_fy_fz_temp(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_fx_fz_temp(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_fx_fy_fz_temp(Nodes* nds, int n,int phase){

};
CPUGPU void bc_grid_vx_vy_vz_temp(Nodes* nds, int n,int phase){

};

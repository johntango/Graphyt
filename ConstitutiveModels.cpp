#include "ConstitutiveModels.hpp"

CPUGPU void elastic(Parameters *mpmParameters,MaterialPoints* mpts,double dt,int mp,std::vector<MaterialModels>* materials){
    // Decompose strain rates into volumetric and deviatoric
double strainRate_vol[3]={0.0,0.0,0.0};
double strainRate_dev[6]={0.0,0.0,0.0,0.0,0.0,0.0};
int matID = mpts->material[mp];
double G  =  (*materials)[matID].shearMod;
double K  =  (*materials)[matID].bulkMod;

// Decompose strain rates into volumetric and deviatoric 
strainRate_vol[0] = 0.3333 * (mpts->strainRate[mp*6+0] + mpts->strainRate[mp*6+1] +mpts->strainRate[mp*6+2] );
strainRate_vol[1] = 0.3333 * (mpts->strainRate[mp*6+0] + mpts->strainRate[mp*6+1] +mpts->strainRate[mp*6+2] );
strainRate_vol[2] = 0.3333 * (mpts->strainRate[mp*6+0] + mpts->strainRate[mp*6+1] +mpts->strainRate[mp*6+2] );

strainRate_dev[0] = mpts->strainRate[mp*6+0] - strainRate_vol[0];
strainRate_dev[1] = mpts->strainRate[mp*6+1] - strainRate_vol[1];
strainRate_dev[2] = mpts->strainRate[mp*6+2] - strainRate_vol[2];
strainRate_dev[3] = mpts->strainRate[mp*6+3];
strainRate_dev[4] = mpts->strainRate[mp*6+4];
strainRate_dev[5] = mpts->strainRate[mp*6+5];
// Calculate temporary stresses assuming elastic step
mpts->stress[mp*6+0] += (2.0*G*strainRate_dev[0] + 3.0*K*strainRate_vol[0]) * dt;
mpts->stress[mp*6+1] += (2.0*G*strainRate_dev[1] + 3.0*K*strainRate_vol[1]) * dt;
mpts->stress[mp*6+2] += (2.0*G*strainRate_dev[2] + 3.0*K*strainRate_vol[2]) * dt;
mpts->stress[mp*6+3] += (2.0*G*strainRate_dev[3]                        ) * dt;
mpts->stress[mp*6+4] += (2.0*G*strainRate_dev[4]                        ) * dt;
mpts->stress[mp*6+5] += (2.0*G*strainRate_dev[5]                        ) * dt;

};


CPUGPU void newtonianLiquid(Parameters *mpmParameters,MaterialPoints* mpts,double dt,int mp,std::vector<MaterialModels>* materials){

double strainRate_vol[3]={0.0,0.0,0.0};
double strainRate_dev[6]={0.0,0.0,0.0,0.0,0.0,0.0};
int matID = mpts->material[mp];
double K =  (*materials)[matID].bulkMod;
double mu =  (*materials)[matID].viscosity;

strainRate_vol[0] = 0.33333 * (mpts->strainRate[mp*6+0] + mpts->strainRate[mp*6+1] +mpts->strainRate[mp*6+2] );
strainRate_vol[1] = 0.33333 * (mpts->strainRate[mp*6+0] + mpts->strainRate[mp*6+1] +mpts->strainRate[mp*6+2] );
strainRate_vol[2] = 0.33333 * (mpts->strainRate[mp*6+0] + mpts->strainRate[mp*6+1] +mpts->strainRate[mp*6+2] );

double P = K*(mpts->strain[mp*6+0]+mpts->strain[mp*6+1]+mpts->strain[mp*6+2]);
//P = (2.15e9/7.15)*( pow((1000.0/mpts->density[mp]),7.15) - 1.0) + 101325.0;

strainRate_dev[0] =   mpts->strainRate[mp*6+0] - strainRate_vol[0];
strainRate_dev[1] =   mpts->strainRate[mp*6+1] - strainRate_vol[1];
strainRate_dev[2] =   mpts->strainRate[mp*6+2] - strainRate_vol[2];
strainRate_dev[3] =   mpts->strainRate[mp*6+3];
strainRate_dev[4] =   mpts->strainRate[mp*6+4];
strainRate_dev[5] =   mpts->strainRate[mp*6+5];

// total stress = pI + 2*mueij - (2/3)mu*tr(eij)I
mpts->stress[mp*6+0] = P + mu*(2.-(2./3.))*(strainRate_dev[0])*dt;
mpts->stress[mp*6+1] = P + mu*(2.-(2./3.))*(strainRate_dev[1])*dt;
mpts->stress[mp*6+3] = 2.0*mu*(strainRate_dev[3])*dt;
if(mpmParameters->is3D){
mpts->stress[mp*6+2] = P + mu*(2.-(2./3.))*(strainRate_dev[2])*dt;
   mpts->stress[mp*6+4] = 2.0*mu*(strainRate_dev[4])*dt;
   mpts->stress[mp*6+5] = 2.0*mu*(strainRate_dev[5])*dt;
}
};

CPUGPU void turbulentLiquid(Parameters *mpmParameters,MaterialPoints* mpts,double dt,int mp,std::vector<MaterialModels>* materials){

double strainRate_vol[3]={0.0,0.0,0.0};
double strainRate_dev[6]={0.0,0.0,0.0,0.0,0.0,0.0};
int matID = mpts->material[mp];
double K =  (*materials)[matID].bulkMod;
double mu_mat =  (*materials)[matID].viscosity;
double turbStress = (*materials)[matID].turbStress;
double mu = mu_mat;


strainRate_dev[0] =   mpts->strainRate[mp*6+0] - strainRate_vol[0];
strainRate_dev[1] =   mpts->strainRate[mp*6+1] - strainRate_vol[1];
strainRate_dev[2] =   mpts->strainRate[mp*6+2] - strainRate_vol[2];
strainRate_dev[3] =   mpts->strainRate[mp*6+3];
strainRate_dev[4] =   mpts->strainRate[mp*6+4];
strainRate_dev[5] =   mpts->strainRate[mp*6+5];

// model 1 - Turbulence exhibited as a decrease in instantaneous viscosity
// Check for exisiting turbulescence
if(mpts->damage[mp]==1) // material cannot connect with neighbors so viscosity is effectivly reduced?
{
    mu=mu_mat/2.0;
}
//Check for onset of turbulescence
// Based on the criterion-> KE = VisEnergy + turblenceStress
//KE = 1/2 mass(vx^2+vy^2+vz^2)
//VE = 1/2 mu(e_rate_dev_xy^2+e_rate_dev_zy^2+e_rate_dev_xz^2) 
double KE = 0.5 * mpts->mass[mp] * ((mpts->velocity[mp*3+0] * mpts->velocity[mp*3+0]) + (mpts->velocity[mp*3+1] * mpts->velocity[mp*3+1]) + (mpts->velocity[mp*3+2] * mpts->velocity[mp*3+2]) );
double VE = 0.5*mu_mat*(strainRate_dev[3]*strainRate_dev[3] + strainRate_dev[4]*strainRate_dev[4] +strainRate_dev[5]*strainRate_dev[5]);
double turbulent_criteria = KE-VE-turbStress;
if(turbulent_criteria>0){
    mpts->damage[mp]= 1.0;
};
// Dissipation will cause the turbulence to fade away, this will restore the regular viscosity
if(turbulent_criteria<=0)
{
    mpts->damage[mp]= 0.0;
    mu =  mu_mat; // set mu back to material property
}



strainRate_vol[0] = 0.33333 * (mpts->strainRate[mp*6+0] + mpts->strainRate[mp*6+1] +mpts->strainRate[mp*6+2] );
strainRate_vol[1] = 0.33333 * (mpts->strainRate[mp*6+0] + mpts->strainRate[mp*6+1] +mpts->strainRate[mp*6+2] );
strainRate_vol[2] = 0.33333 * (mpts->strainRate[mp*6+0] + mpts->strainRate[mp*6+1] +mpts->strainRate[mp*6+2] );

double P = K*(mpts->strain[mp*6+0]+mpts->strain[mp*6+1]+mpts->strain[mp*6+2]);
//P = (2.15e9/7.15)*( pow((1000.0/mpts->density[mp]),7.15) - 1.0) + 101325.0;


// total stress = pI + 2*mueij - (2/3)mu*tr(eij)I
mpts->stress[mp*6+0] = P + mu*(2.-(2./3.))*(strainRate_dev[0])*dt;
mpts->stress[mp*6+1] = P + mu*(2.-(2./3.))*(strainRate_dev[1])*dt;
mpts->stress[mp*6+3] = 2.0*mu*(strainRate_dev[3])*dt;
if(mpmParameters->is3D){
mpts->stress[mp*6+2] = P + mu*(2.-(2./3.))*(strainRate_dev[2])*dt;
   mpts->stress[mp*6+4] = 2.0*mu*(strainRate_dev[4])*dt;
   mpts->stress[mp*6+5] = 2.0*mu*(strainRate_dev[5])*dt;
}
}




CPUGPU void perfectPlasticity(Parameters *mpmParameters,MaterialPoints* mpts,double dt,int mp,std::vector<MaterialModels>* materials){

double stress[6]            = {0.0,0.0,0.0,0.0,0.0,0.0};
double strainRate_vol[3]    = {0.0,0.0,0.0};
double strainRate_dev[6]    = {0.0,0.0,0.0,0.0,0.0,0.0};
int    matID                = mpts->material[mp];
double G                    = (*materials)[matID].shearMod;
double K                    = (*materials)[matID].bulkMod;
double yieldStress          = (*materials)[matID].yieldStress;
bool isDamage                = (*materials)[matID].isDamage;
// Decompose strain rates into volumetric and deviatoric 
strainRate_vol[0] = 0.33 * (mpts->strainRate[mp*6+0] + mpts->strainRate[mp*6+1] +mpts->strainRate[mp*6+2] );
strainRate_vol[1] = 0.33 * (mpts->strainRate[mp*6+0] + mpts->strainRate[mp*6+1] +mpts->strainRate[mp*6+2] );
strainRate_vol[2] = 0.33 * (mpts->strainRate[mp*6+0] + mpts->strainRate[mp*6+1] +mpts->strainRate[mp*6+2] );
strainRate_dev[0] = mpts->strainRate[mp*6+0] - strainRate_vol[0];
strainRate_dev[1] = mpts->strainRate[mp*6+1] - strainRate_vol[1];
strainRate_dev[2] = mpts->strainRate[mp*6+2] - strainRate_vol[2];
strainRate_dev[3] = mpts->strainRate[mp*6+3];
strainRate_dev[4] = mpts->strainRate[mp*6+4];
strainRate_dev[5] = mpts->strainRate[mp*6+5];
// Calculate temporary stresses assuming elastic step
stress[0] = mpts->stress[mp*6+0] + (2*G*strainRate_dev[0] + 3*K*strainRate_vol[0]) * dt;
stress[1] = mpts->stress[mp*6+1] + (2*G*strainRate_dev[1] + 3*K*strainRate_vol[1]) * dt;
if(mpmParameters->is3D){
stress[2] = mpts->stress[mp*6+2] + (2*G*strainRate_dev[2] + 3*K*strainRate_vol[2]) * dt;
stress[4] = mpts->stress[mp*6+4] + (2*G*strainRate_dev[4]                        ) * dt;
stress[5] = mpts->stress[mp*6+5] + (2*G*strainRate_dev[5]                        ) * dt;
}
stress[3] = mpts->stress[mp*6+3] + (2*G*strainRate_dev[3]                        ) * dt;

// Calculte VonMises Stress
double equivStress = sqrt( 0.5*pow((stress[0]-stress[1]),2.0) + 0.5*pow((stress[1]-stress[2]),2.0) + 0.5*pow((stress[2]-stress[0]),2.0) + pow(stress[3],2.0) + pow(stress[4],2.0) + pow(stress[5],2.0) );
// Check Yield Surface

if(equivStress > yieldStress)
{
    // Material has yielded, need to decompose elastic and plastic strains
    double plasticStrainRate[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    double elasticStrainRate[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    strainRate_dev[0] = 0.0;
    strainRate_dev[1] = 0.0;
    strainRate_dev[2] = 0.0;
    strainRate_dev[3] = 0.0;
    strainRate_dev[4] = 0.0;
    strainRate_dev[5] = 0.0;
    

    // Find the increment of plastic strain 
    double dL = ((equivStress-yieldStress) / (3.0*G) );



    // Use the flow direction to calculate plastic strain increment

    plasticStrainRate[0] = dL*(0.5/equivStress)*(2.0*stress[0] - stress[1] - stress[2]); 
    plasticStrainRate[1] = dL*(0.5/equivStress)*(2.0*stress[1] - stress[2] - stress[0]); 
    plasticStrainRate[2] = dL*(0.5/equivStress)*(2.0*stress[2] - stress[0] - stress[1]); 
    plasticStrainRate[3] = dL*(0.5/equivStress)*(6.0*stress[3]);
    plasticStrainRate[4] = dL*(0.5/equivStress)*(6.0*stress[4]);
    plasticStrainRate[5] = dL*(0.5/equivStress)*(6.0*stress[5]);
    mpts->acc_plastic_strain[mp] += dt*sqrt( 0.5*pow((plasticStrainRate[0]-plasticStrainRate[1]),2.0) + 0.5*pow((plasticStrainRate[1]-plasticStrainRate[2]),2.0) + 0.5*pow((plasticStrainRate[2]-plasticStrainRate[0]),2.0) + pow(plasticStrainRate[3],2.0) + pow(plasticStrainRate[4],2.0) + pow(plasticStrainRate[5],2.0) );
    // Increase the mat point plastic strain
    mpts->plastic_strain[mp*6+0] += plasticStrainRate[0];
    mpts->plastic_strain[mp*6+1] += plasticStrainRate[1];
    mpts->plastic_strain[mp*6+2] += plasticStrainRate[2];
    mpts->plastic_strain[mp*6+3] += plasticStrainRate[3];
    mpts->plastic_strain[mp*6+4] += plasticStrainRate[4];
    mpts->plastic_strain[mp*6+5] += plasticStrainRate[5];
    // Reduce the amount of available strain for elastic update
    elasticStrainRate[0] = mpts->strainRate[mp*6+0] - (plasticStrainRate[0]/dt);
    elasticStrainRate[1] = mpts->strainRate[mp*6+1] - (plasticStrainRate[1]/dt);
    elasticStrainRate[2] = mpts->strainRate[mp*6+2] - (plasticStrainRate[2]/dt);
    elasticStrainRate[3] = mpts->strainRate[mp*6+3] - (plasticStrainRate[3]/dt);
    elasticStrainRate[4] = mpts->strainRate[mp*6+4] - (plasticStrainRate[4]/dt);
    elasticStrainRate[5] = mpts->strainRate[mp*6+5] - (plasticStrainRate[5]/dt);
    // Now calculate the elastic strains and update the stress
    strainRate_dev[0] = elasticStrainRate[0] - strainRate_vol[0];
    strainRate_dev[1] = elasticStrainRate[1] - strainRate_vol[1];
    strainRate_dev[2] = elasticStrainRate[2] - strainRate_vol[2];
    strainRate_dev[3] = elasticStrainRate[3];
    strainRate_dev[4] = elasticStrainRate[4];
    strainRate_dev[5] = elasticStrainRate[5];
    // Now update stresses using hooke's law
    mpts->stress[mp*6+0] += (2.0*G*strainRate_dev[0] + 3.0*K*strainRate_vol[0]) * dt;
    mpts->stress[mp*6+1] += (2.0*G*strainRate_dev[1] + 3.0*K*strainRate_vol[1]) * dt;
    mpts->stress[mp*6+3] += (2.0*G*strainRate_dev[3]                          ) * dt;
    if(mpmParameters->is3D){
    mpts->stress[mp*6+2] += (2.0*G*strainRate_dev[2] + 3.0*K*strainRate_vol[2]) * dt;
    
    mpts->stress[mp*6+4] += (2.0*G*strainRate_dev[4]                          ) * dt;
    mpts->stress[mp*6+5] += (2.0*G*strainRate_dev[5]                          ) * dt;
    }
}
else{
    // Elastic so update stresses elastically
    mpts->stress[mp*6+0] = stress[0];
    mpts->stress[mp*6+1] = stress[1];
    mpts->stress[mp*6+3] = stress[3];
    mpts->stress[mp*6+2] = stress[2];
    mpts->stress[mp*6+4] = stress[4];
    mpts->stress[mp*6+5] = stress[5];
}

};



CPUGPU void druckerPrager_Damage(Parameters *mpmParameters,MaterialPoints* mpts,double dt,int mp,std::vector<MaterialModels>* materials){
double stress[6]            = {0.0,0.0,0.0,0.0,0.0,0.0};
double strainRate_vol[3]    = {0.0,0.0,0.0};
double strainRate_dev[6]    = {0.0,0.0,0.0,0.0,0.0,0.0};
double plasticStrainRate[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
int matID                   = mpts->material[mp];
double G                    = (*materials)[matID].shearMod;
double K                    = (*materials)[matID].bulkMod;
double c                    = (*materials)[matID].cohesion;
double phi                  = (*materials)[matID].phi * 3.14159/180.0;
double psi                  = (*materials)[matID].psi * 3.14159/180.0;
double eps_cr               = (*materials)[matID].criticalStrain;
bool isDamage                = (*materials)[matID].isDamage;
double F;
double J2;
double eta;
double eta_bar;
double xi;
double I1; 
double dsxx,dsyy,dszz,dsxy,dsyz,dsxz;
// Decompose strain rates into volumetric and deviatoric 
strainRate_vol[0] = 0.3333 * (mpts->strainRate[mp*6+0] + mpts->strainRate[mp*6+1] +mpts->strainRate[mp*6+2] );
strainRate_vol[1] = 0.3333 * (mpts->strainRate[mp*6+0] + mpts->strainRate[mp*6+1] +mpts->strainRate[mp*6+2] );
strainRate_vol[2] = 0.3333 * (mpts->strainRate[mp*6+0] + mpts->strainRate[mp*6+1] +mpts->strainRate[mp*6+2] );

strainRate_dev[0] = mpts->strainRate[mp*6+0] - strainRate_vol[0];
strainRate_dev[1] = mpts->strainRate[mp*6+1] - strainRate_vol[1];
strainRate_dev[2] = mpts->strainRate[mp*6+2] - strainRate_vol[2];
strainRate_dev[3] = mpts->strainRate[mp*6+3];
strainRate_dev[4] = mpts->strainRate[mp*6+4];
strainRate_dev[5] = mpts->strainRate[mp*6+5];
// Calculate temporary stresses assuming elastic step
stress[0] = mpts->stress[mp*6+0] + (2.0*G*strainRate_dev[0] + 3.0*K*strainRate_vol[0]) * dt;
stress[1] = mpts->stress[mp*6+1] + (2.0*G*strainRate_dev[1] + 3.0*K*strainRate_vol[1]) * dt;
stress[2] = mpts->stress[mp*6+2] + (2.0*G*strainRate_dev[2] + 3.0*K*strainRate_vol[2]) * dt;
stress[3] = mpts->stress[mp*6+3] + (2.0*G*strainRate_dev[3]                        ) * dt;
stress[4] = mpts->stress[mp*6+4] + (2.0*G*strainRate_dev[4]                        ) * dt;
stress[5] = mpts->stress[mp*6+5] + (2.0*G*strainRate_dev[5]                        ) * dt;


I1 = (1.0/3.0)*(stress[0]+stress[1]+stress[2]);
dsxx = stress[0] - I1;
dsyy = stress[1] - I1;
dszz = stress[2] - I1;
dsxy = stress[3];
dsyz = stress[4];
dsxz = stress[5];
J2 = (0.5) * (dsxx*dsxx + dsyy*dsyy + dszz*dszz + 2.0*dsxy*dsxy + 2.0*dsyz*dsyz + 2.0*dsxz*dsxz);

eta = 6.0*sin(phi) / (sqrt(3.0)*(3.0-sin(phi)));
xi  = 6.0*cos(phi) / (sqrt(3.0)*(3.0-sin(phi)));
if(phi==psi) //Associative
{
    eta_bar = eta;
}
else{
eta_bar = 6.0*sin(psi) / (sqrt(3.0)*(3.0 - sin(psi)));
}


F = sqrt(J2) + eta*I1 - xi*c;
bool criterion;
if(isDamage)
{
    criterion = (F>0.0)&&(I1< 0.0);
}
else
{
    criterion = F>0.0;
} 


if(criterion)
{
    // Plastic update
    
    // start N-R loop
    int iteration = 0;
    double error = 1.0;
    double dL_old = 0.0;
    double dL_new = 10.0;
    double Fp = 0.0;
    double equivpStrain = mpts->acc_plastic_strain[mp];
    while(error >1e-8)
    {
        iteration++;
        // Assuming linear cohesion
        double cval = c - (equivpStrain+(xi*dL_old))*(c/eps_cr);
        double dc_val = -c/eps_cr;
            if(cval<c/10.0)
            {
                cval = c/10.0;
                dc_val = 0.0;
            }
        
        Fp = sqrt(J2) - G*dL_old + eta*(I1 - eta_bar*K*dL_old) - xi*cval;
        double dFp = -G - K*eta*eta_bar - xi*xi*dc_val;
        
        dL_new = dL_old - (Fp/dFp);
        
        error = fabs(Fp/cval);
        dL_old = dL_new;
        if(iteration >200){
        break;			
        }
        
    }

    // Now need to check for the apex
    double G_apex =sqrt(J2) - G*dL_old;
    if ((G_apex < 0.0) && (!isDamage))
    {
        iteration =0;
        error = 1.0;
        double deps= 0.0;
        while(error>1e-8)
        {
            iteration ++;
            // Assuming linear cohesion
            double cval = c - (equivpStrain+((xi/eta)*deps))*(c/eps_cr);
            double dc_val = -c/eps_cr;
            if(cval<c/10.0)
            {
                cval = c/10.0;
                dc_val = 0.0;
            }
            double Gp = cval*(xi/eta_bar) - I1 + K*deps;
            double dGp = (xi/eta)*dc_val*(xi/eta_bar) + K;
            deps = deps- (Gp/dGp);
            error =fabs(Gp/cval );

            if(iteration >200){
        
                break;			
            }

        }

        if(deps <0.0){
                deps=0.0;
        }
        else{
            mpts->stress[6*mp+0] = I1 -  K*deps;
            mpts->stress[6*mp+1] = I1 -  K*deps;
            mpts->stress[6*mp+3] = 0.0;
            mpts->stress[6*mp+2] = I1 - K*deps;
            mpts->stress[6*mp+4] = 0.0;
            mpts->stress[6*mp+5] = 0.0;

            mpts->acc_plastic_strain[mp] += deps*(xi/eta);
            mpts->plastic_strain[mp*6+0] += deps*( eta_bar/3.0);
            mpts->plastic_strain[mp*6+1] += deps*( eta_bar/3.0);
            mpts->plastic_strain[mp*6+2] += deps*(eta_bar/3.0);
            mpts->plastic_strain[mp*6+3] += 0.0;
            mpts->plastic_strain[mp*6+4] += 0.0;
            mpts->plastic_strain[mp*6+5] += 0.0;
        }
    }
    else{

            mpts->stress[6*mp+0] = (1.0-(G*dL_old/sqrt(J2)) )*dsxx + (I1 - K*eta_bar*dL_old);
            mpts->stress[6*mp+1] = (1.0-(G*dL_old/sqrt(J2)) )*dsyy + (I1 - K*eta_bar*dL_old);
            mpts->stress[6*mp+3] = (1.0-(G*dL_old/sqrt(J2)) )*dsxy;
            mpts->stress[6*mp+2] = (1.0-(G*dL_old/sqrt(J2)) )*dszz + (I1 - K*eta_bar*dL_old);
            mpts->stress[6*mp+4] = (1.0-(G*dL_old/sqrt(J2)) )*dsyz;
            mpts->stress[6*mp+5] = (1.0-(G*dL_old/sqrt(J2)) )*dsxz;

            mpts->acc_plastic_strain[mp] += dL_old*xi;
            mpts->plastic_strain[mp*6+0] += dL_old*((0.5/J2)*dsxx + eta_bar/3.0);
            mpts->plastic_strain[mp*6+1] += dL_old*((0.5/J2)*dsyy + eta_bar/3.0);
            mpts->plastic_strain[mp*6+2] += dL_old*((0.5/J2)*dszz + eta_bar/3.0);
            mpts->plastic_strain[mp*6+3] += dL_old*(0.5/J2)*dsxy;
            mpts->plastic_strain[mp*6+4] += dL_old*(0.5/J2)*dsyz;
            mpts->plastic_strain[mp*6+5] += dL_old*(0.5/J2)*dsxz;
            
        
    }
}
else
{
    // Elastic so update stresses elastically
    mpts->stress[mp*6+0] = stress[0];
    mpts->stress[mp*6+1] = stress[1];
    mpts->stress[mp*6+3] = stress[3];
    mpts->stress[mp*6+2] = stress[2];
    mpts->stress[mp*6+4] = stress[4];
    mpts->stress[mp*6+5] = stress[5];

}


};
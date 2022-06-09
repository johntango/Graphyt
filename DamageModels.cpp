#include "DamageModels.hpp"

CPUGPU void gradyKippDamage(Parameters *mpmParameters,MaterialPoints* mpts,double dt,int mp,std::vector<MaterialModels>* materials){
	int matID                   = mpts->material[mp];
	double Cg                   = (*materials)[matID].Cg;
	double m                    = (*materials)[matID].damageM;
	double k                    = (*materials)[matID].damageK;
	double G                    = (*materials)[matID].shearMod;
	double K                    = (*materials)[matID].bulkMod;
        
        if(!mpmParameters->is3D)
        {
                double str[4];
        
        // 2D implementation
                str[0] = mpts->stress[6*mp + 0];
                str[1] = mpts->stress[6*mp + 1];
                str[2] = mpts->stress[6*mp + 3];
                str[3] = mpts->stress[6*mp + 2];
                double theta, s1, s2;
                if (fabs(str[0] - str[1]) < 1e-5) theta = 3.14159 / 4.0;
                else theta = 0.5 * atan(2.0 * str[2] / (str[0] - str[1]));
                double c = cos(theta);
                double s = sin(theta);
                s1 = str[0] * c * c + str[1] * s * s + 2.0 * c * s * str[2];
                s2 = str[0] * s * s + str[1] * c * c - 2.0 * c * s * str[2];
                double strmax = std::max(s1, std::max(s2, 0.0));
        
                mpts->smax[mp] = strmax;
        
                double eff_strain, emin;
                double alpha = (8.0 * 3.14159 * Cg * Cg * Cg * k) / (((m + 1.0) * (m + 2.0) * (m + 3.0)));
        //	double alpha1 = (8.0*3.14159*k*(m+3.0)*(m+3.0)) / ((m+1.0)*(m+2.0));
                if (strmax > 1e-3)
                {
                        eff_strain = strmax/(K + (4.0*G/3.0) );
                        double vol = mpts->mass[mp]/mpts->density[mp];
                        emin = pow(vol * k, -1.0 / m);
                        if (eff_strain > emin)
                        {
                                double dD13dt = ((m + 3.0) / 3.0) * pow(alpha, 1.0 / 3.0) * pow(eff_strain, m / 3.0);
                                //double dD13dt1 = (1.0/3.0)*Cg*pow(alpha1,(1.0/3.0))*pow(eff_strain,(m/3.0));
                                mpts->damage[mp] = pow(pow(mpts->damage[mp], 1.0 / 3.0) + (dt * dD13dt), 3.0);
                                if (mpts->damage[mp] > 1.0) mpts->damage[mp] = 1.0;
                                if (s1 > 0.0) s1 = (1.0 - mpts->damage[mp]) * s1;
                                if (s2 > 0.0) s2 = (1.0 - mpts->damage[mp]) * s2;
                                str[0] = s1 * c * c + s2 * s * s;
                                str[1] = s1 * s * s + s2 * c * c;
                                str[2] = (s1 - s2) * c * s;
                                str[3] = 0.0;
                        }

                    }
                if(mpts->damage[mp]==1.0)
                {
                        if (mpts->damage[mp] > 1.0) mpts->damage[mp] = 1.0;
                        if (s1 > 0.0) s1 = (1.0 - mpts->damage[mp]) * s1;
                        if (s2 > 0.0) s2 = (1.0 - mpts->damage[mp]) * s2;
                        str[0] = s1 * c * c + s2 * s * s;
                        str[1] = s1 * s * s + s2 * c * c;
                        str[2] = (s1 - s2) * c * s;
                        str[3] = 0.0;      
                }
        

                mpts->stress[6*mp + 0] = str[0]; //xx
                mpts->stress[6*mp + 1] = str[1]; //yy
                mpts->stress[6*mp + 2] = str[3]; //zz
                mpts->stress[6*mp + 3] = str[2]; //xy
        }
        else{ // Full 3D implementation
                double s[3][3];
                s[0][0] = mpts->stress[6*mp+ 0];//xx
                s[0][1] = mpts->stress[6*mp+ 3];//xy
                s[0][2] = mpts->stress[6*mp+ 5];//xz
                s[1][0] = mpts->stress[6*mp+ 3];//xy
                s[1][1] = mpts->stress[6*mp+ 1];//yy
                s[1][2] = mpts->stress[6*mp+ 4];//yz
                s[2][0] = mpts->stress[6*mp+ 5];//xz
                s[2][1] = mpts->stress[6*mp+ 4];//yz
                s[2][2] = mpts->stress[6*mp+ 2];//zz
                double checkSum = fabs(s[0][0]+s[1][1] +s[2][2] +s[0][1] +s[1][2] +s[0][2]);
                if(checkSum >1e-3)
                {
                        double eigenvectors[3][3];
                        double eigenvalues[3];
                        eigen_decomposition(s, eigenvectors, eigenvalues);
                        double s1 = eigenvalues[0];
                        double s2 = eigenvalues[1];
                        double s3 = eigenvalues[2];
                        double strmax = std::max(s1, std::max(s2, std::max(s3, 0.0)));
                
                        mpts->smax[mp] = strmax;
                
                        double eff_strain, emin;
                        double alpha = (8.0 * 3.14159 * Cg * Cg * Cg * k) / (((m + 1.0) * (m + 2.0) * (m + 3.0)));
                
                        if(strmax > 1e-3)
                        {
                                
                                eff_strain = strmax / (K + 4.0 * G / 3.0);
                                double vol = mpts->mass[mp]/mpts->density[mp];
                                emin = pow(vol * k, -1.0 / m);
                                if (eff_strain > emin)
                                {
                                    double dD13dt = ((m + 3.0) / 3.0) * pow(alpha, 1.0 / 3.0) * pow(eff_strain, m / 3.0);
                                    //double dD13dt1 = (1.0/3.0)*Cg*pow(alpha1,(1.0/3.0))*pow(eff_strain,(m/3.0));
                                    mpts->damage[mp] = pow(pow(mpts->damage[mp], 1.0 / 3.0) + (dt * dD13dt), 3.0);
                                    if (mpts->damage[mp] > 1.0) mpts->damage[mp] = 1.0;
                                    if (s1 > 0.0) s1 = (1.0 - mpts->damage[mp]) * s1;
                                    if (s2 > 0.0) s2 = (1.0 - mpts->damage[mp]) * s2;
                                    if (s3 > 0.0) s3 = (1.0 - mpts->damage[mp]) * s3;
                            }
                        //     if(mpts->damage[mp]>0.5)
                        //     {
                        //         if (mpts->damage[mp] > 1.0) mpts->damage[mp] = 1.0;
                        //         if (s1 > 0.0) s1 = (1.0 - mpts->damage[mp]) * s1;
                        //         if (s2 > 0.0) s2 = (1.0 - mpts->damage[mp]) * s2;
                        //         if (s3 > 0.0) s3 = (1.0 - mpts->damage[mp]) * s3;
                                
    
                        //     }
                        }
                        double diag[3][3];
                        diag[0][0] = s1;  diag[0][1] = 0.0; diag[0][2] = 0.0;
                        diag[1][0] = 0.0; diag[1][1] = s2;  diag[1][2] = 0.0;
                        diag[2][0] = 0.0; diag[2][1] = 0.0; diag[2][2] = s3;
                        double norm1 = sqrt(eigenvectors[0][0] * eigenvectors[0][0] + eigenvectors[1][0] * eigenvectors[1][0] + eigenvectors[2][0] * eigenvectors[2][0]);
                        double norm2 = sqrt(eigenvectors[0][1] * eigenvectors[0][1] + eigenvectors[1][1] * eigenvectors[1][1] + eigenvectors[2][1] * eigenvectors[2][1]);
                        double norm3 = sqrt(eigenvectors[0][2] * eigenvectors[0][2] + eigenvectors[1][2] * eigenvectors[1][2] + eigenvectors[2][2] * eigenvectors[2][2]);
                        double evnorm[3][3]; // probably not necessary since eigenvectors are already normalized
                        evnorm[0][0] = eigenvectors[0][0] / norm1; evnorm[0][1] = eigenvectors[0][1] / norm2; evnorm[0][2] = eigenvectors[0][2] / norm3;
                        evnorm[1][0] = eigenvectors[1][0] / norm1; evnorm[1][1] = eigenvectors[1][1] / norm2; evnorm[1][2] = eigenvectors[1][2] / norm3;
                        evnorm[2][0] = eigenvectors[2][0] / norm1; evnorm[2][1] = eigenvectors[2][1] / norm2; evnorm[2][2] = eigenvectors[2][2] / norm3;
                        double evnormtransposed[3][3];
                        evnormtransposed[0][0] = evnorm[0][0]; evnormtransposed[0][1] = evnorm[1][0]; evnormtransposed[0][2] = evnorm[2][0];
                        evnormtransposed[1][0] = evnorm[0][1]; evnormtransposed[1][1] = evnorm[1][1]; evnormtransposed[1][2] = evnorm[2][1];
                        evnormtransposed[2][0] = evnorm[0][2]; evnormtransposed[2][1] = evnorm[1][2]; evnormtransposed[2][2] = evnorm[2][2];
                        double temp[3][3];
                        double finalstress[3][3];
                        matrixMatrixMultiply(diag, evnormtransposed, temp);
                        matrixMatrixMultiply(evnorm, temp, finalstress);
                        mpts->stress[6*mp + 0] = finalstress[0][0] ;
                        mpts->stress[6*mp + 1] = finalstress[1][1];
                        mpts->stress[6*mp + 2] = finalstress[2][2];
                        mpts->stress[6*mp + 3] = finalstress[0][1] ;
                        mpts->stress[6*mp + 4] = finalstress[1][2];
                        mpts->stress[6*mp + 5] = finalstress[0][2];
                }
        }

};
#!/usr/bin/python
import numpy as np
import random
 
# Written by Joe Morris, jmorris4@slb.com
# References:
# Morris, J. P., "A numerical investigation of the scaling of fracture stiffness," ARMA 12-610, in 46th U.S. Rock Mechanics/Geomechanics Symposium, June 24-27, Chicago, 2012
# Morris, J. P., Jocker, J., and Prioul, R., "Exploring alternative characterizations of fracture stiffness and their respective scaling behaviors," ARMA 13-197, in the 47th U.S. Rock Mechanics/Geomechanics Symposium, June 23-26, San Francisco, 2013
def drazerKoplik(nx,ny,zeta,seed):
    # zeta is the Hurst exponent
    # zeta=0.0 # Low numbers lead to grainy spatial distributions
    # zeta=10.0 # High numbers lead to smooth surfaces
    # zeta=-1.0 # Will recover the original, spatially uncorrelated Gaussian random distribution
    #zeta=0.8 # The literature indicates that this is typical of fractures
 
    # To get reproducible results we can set the seed
    random.seed(seed)
    
    zrand = np.zeros([nx,ny],float)
    for j in range(0, ny):
        for i in range (0, nx):
            zrand[i,j] = random.normalvariate(0.0,1.0);
    Zorig = np.fft.fft2( zrand )
    Z=np.zeros([nx,ny],float)
    ifreqs=np.fft.fftfreq(nx)
    jfreqs=np.fft.fftfreq(ny)
    for j in range(0, ny):
        for i in range (0, nx):
            ifreq=ifreqs[i]
            jfreq=jfreqs[j]
            zfreq=complex(ifreq,jfreq)
            if ( (i != 0) or (j != 0) ):
                Z[i,j] = zfreq**(-zeta-1.0)*Zorig[i,j]
    z = np.fft.ifft2( Z )
    ap = z.real + z.imag
    mmax=ap.max()
    mmin=ap.min()
    ap=(ap-mmin)/(mmax-mmin)
    return ap
 

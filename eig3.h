// Eigen-decomposition for symmetric 3x3 real matrices.
// Public domain, copied from the public domain Java library JAMA

#ifndef _eig_h
#include "Parameters.hpp"
// Symmetric matrix A => eigenvectors in columns of V, corresponding eigenvalues in d
CPUGPU void eigen_decomposition(double A[3][3], double V[3][3], double d[3]);

#endif

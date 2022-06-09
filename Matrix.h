// Groups all the 3D matrix-vector related operations

#ifndef _Matrix_h
#include "Parameters.hpp"
// Matrix Matrix multiplication
CPUGPU void matrixMatrixMultiply(double m[3][3], double n[3][3], double result[3][3]);

// Matrix Vector Multiplication
CPUGPU void matrixVectorMultiply(double m[3][3], double v[3], double result[3]);

// Vector Magnitude
CPUGPU double vectorMagnitude(double v[3]);

// Vector Print
CPUGPU void vectorPrint(double v[3]);

// Matrix Print
CPUGPU void matrixPrint(double m[3][3]);

// Matrix Orthonormalization Gram-Schmidt
CPUGPU void matrixOrthonormalize(double m[3][3], double result[3][3]);

#endif
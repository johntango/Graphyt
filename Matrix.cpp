/**
    Spark / SOLIDSSPH
    Matrix.cpp
    Purpose: Groups all 3D matrix-related functions

    @author MIT Geonumerics
    @version 1.0 01/01/01
*/

#include "Matrix.h"
#include <math.h>
#include <stdio.h>

/**
    Computes the magnitude of a 3D vector

    @param v An array of 3 elements
    @return the magnitude of the vector
*/
CPUGPU double vectorMagnitude(double v[3])
{
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}
/**
    Computes the multiplication of two 3x3 matrix

    @param m An array of 3x3 elements
    @param n An array of 3x3 elements
    @param result An array of 3x3 elements
    @return Void.
*/
CPUGPU void matrixMatrixMultiply(double m[3][3], double n[3][3], double result[3][3])
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			result[i][j] = m[i][0] * n[0][j] + m[i][1] * n[1][j] + m[i][2] * n[2][j];
}
/**
    Computes the multiplication of a 3x3 matrix and a 3D vector

    @param m An array of 3x3 elements
    @param v An array of 3 elements
    @param result An array of 3x3 elements
    @return Void.
*/
CPUGPU void matrixVectorMultiply(double m[3][3], double v[3], double result[3])
{
	for (int i = 0; i < 3; i++)
		result[i] = m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2];
}
/**
    Computes the orthonormlization of a 3x3 matrix

    @param m An array of 3x3 elements
    @param result An array of 3x3 elements
    @return Void.
*/
CPUGPU void matrixOrthonormalize(double m[3][3], double result[3][3])
{
	// Gram-Schmidt Orthonormalization

	double v1[3], v2[3], v3[3];
	v1[0] = m[0][0];  v1[1] = m[1][0];  v1[2] = m[2][0];
	v2[0] = m[0][1];  v2[1] = m[1][1];  v2[2] = m[2][1];
	v3[0] = m[0][2];  v3[1] = m[1][2];  v3[2] = m[2][2];

	double magnitude;
	double proj1, proj2;

	// u1
	magnitude = vectorMagnitude(v1);
	v1[0] = v1[0] / magnitude;
	v1[1] = v1[1] / magnitude;
	v1[2] = v1[2] / magnitude;

	// u2
	proj1 = ( v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] ) / ( v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2] );
	v2[0] = v2[0] - v1[0] * proj1;
	v2[1] = v2[1] - v1[1] * proj1;
	v2[2] = v2[2] - v1[2] * proj1;
	magnitude = vectorMagnitude(v2);
	v2[0] = v2[0] / magnitude;
	v2[1] = v2[1] / magnitude;
	v2[2] = v2[2] / magnitude;

	// u3
	proj1 = ( v1[0] * v3[0] + v1[1] * v3[1] + v1[2] * v3[2] ) / ( v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2] );
	proj2 = ( v2[0] * v3[0] + v2[1] * v3[1] + v2[2] * v3[2] ) / ( v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2] );
	v3[0] = v3[0] - v1[0] * proj1 - v2[0] * proj2;
	v3[1] = v3[1] - v1[1] * proj1 - v2[1] * proj2;
	v3[2] = v3[2] - v1[2] * proj1 - v2[2] * proj2;
	magnitude = vectorMagnitude(v3);
	v3[0] = v3[0] / magnitude;
	v3[1] = v3[1] / magnitude;
	v3[2] = v3[2] / magnitude;


	result[0][0] = v1[0];  result[1][0] = v1[1];  result[2][0] = v1[2];
	result[0][1] = v2[0];  result[1][1] = v2[1];  result[2][1] = v2[2];
	result[0][2] = v3[0];  result[1][2] = v3[1];  result[2][2] = v3[2];

}
/**
    Prints a 3D vector

    @param v An array of 3 elements
    @return Void.
*/
CPUGPU void vectorPrint(double v[3])
{
	printf("    %.3e\t %.3e\t %.3e\n", v[0], v[1], v[2]);
}
/**
    Prints a 3x3 matrix

    @param v An array of 3x3 elements
    @return Void.
*/
CPUGPU void matrixPrint(double m[3][3])
{
	for (int i = 0; i < 3; i++)
		vectorPrint(m[i]);
}
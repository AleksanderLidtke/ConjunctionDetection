/*
A set of functions to perform basic vector and matrix operations on 1D vectors or 2D matrices.

@version 1.1.0
@since 23 Oct 2014 09:16:00
@author Aleksander Lidtke
@email al11g09@soton.ac.uk, alek_l@onet.eu

CHANGELOG:
25 Jul 2014 - 1.0.1 - Fixed a small bug in inverseThreeByThree that would make one term slightly off the actual value.
30 Jul 2014 - 1.0.2 - Added the functionality to compute statistics of vectors.
03 Aug 2014 - 1.0.3 - Begun to throw errors when attempting to inverse matrices with 0 determinants.
08 Aug 2014 - 1.0.4 - Changed some of the FatalErrors to MathWarnings.
23 Oct 2014 - 1.1.0 - Added namespace VectorOperations
*/
#pragma once

#include <cmath>
#include <vector>
#include <float.h>
#include "FatalError.h"
#include "MathWarning.h"

namespace VectorOperations{
	void linspace(std::vector<double>* outputVecPtr, double start, double stop, int* noPoints);
	void centralLinspace(std::vector<double>* outputVecPtr, double centre, double spacing, int* noPoints);
	void vectorDifference(std::vector<double>* outputVecPtr, std::vector<double>* vect1Ptr, std::vector<double>* vec2Ptr);
	double vectorMagnitudeSquared(std::vector<double>* vecPtr);
	double vectorMagnitude(std::vector<double>* vecPtr);
	void unitVector(std::vector<double>* vecPtr);
	double dotProduct(std::vector<double>* vec1Ptr, std::vector<double>* vec2Ptr);
	void crossProduct(std::vector<double>* outputVecPtr, std::vector<double>* vec1Ptr, std::vector<double>* vec2Ptr);
	std::vector<double> vectorMultiplyByScalar(std::vector<double>* vecPtr, double coefficient);
	std::vector< std::vector<double> > matrixMultiplyByScalar(std::vector< std::vector<double> >* matPtr, double coefficient);
	std::vector <double> vectorMultiplyByMatrix(std::vector<double>* vecPtr, std::vector< std::vector<double> >* matPtr);
	std::vector< std::vector<double> > matrixMultiplyByMatrix(std::vector< std::vector<double> >* mat1Ptr, std::vector< std::vector<double> >* mat2Ptr);
	std::vector< std::vector<double> > addMatrices(std::vector< std::vector<double> >* mat1Ptr, std::vector< std::vector<double> >* mat2Ptr);
	double twoByTwoDeterminant( std::vector< std::vector<double> >* matPtr );
	double threeByThreeDeterminant( std::vector< std::vector<double> >* matPtr );
	std::vector< std::vector<double> > inverseTwoByTwo( std::vector< std::vector<double> >* matPtr );
	std::vector< std::vector<double> > inverseThreeByThree( std::vector< std::vector<double> >* matPtr );
	std::vector< std::vector<double> > transposeMatrix( std::vector< std::vector<double> >* originalPtr );
	double mean(std::vector<double>* inputVec);
	double standardDeviation(std::vector<double>* inputVec);
}//END namespace VectorOperations

namespace EquationsSolving{
	std::vector< std::vector<double> > luDecomposition(std::vector< std::vector<double> > A, std::vector<int> &index, double &d);
	void luSubst(std::vector < std::vector<double> >* A, std::vector<int>* index, std::vector<double>* b);
}//END namespace EquationsSolving

/*
A class that provides a similar functionality as SpaceObject but for externally-generated ephemerides.

@version 1.0.3
@since 02 Dec 2014 18:13:00
@author Aleksander Lidtke
@email al11g09@soton.ac.uk, alek_l@onet.eu

CHANGELOG:
01 Aug 2014 - 1.0.1 - changed the check whether all ephemeris lines have been read not to throw warnings when the last line is empty (will happen when exporting ephemerides from STK 10.0.0).
23 Oct 2014 - 1.0.2 - added namespaces VectorOperations and EquationsSolving.
02 Dec 2014 - 1.0.3 - Fixed a bug where the transverseUnit_inertial in inertialToRTC would be calculated by crossing XxZ unit vectors instead of ZxX.
*/
#pragma once

#include <cmath>
#include <math.h>
#include <vector>
#include <string>
#include <random>
#include <sstream>
#include <fstream>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "sgp4ext.h"
#include "sgp4unit.h"
#include "VectorOperations.h"

#define PI 3.14159265358979323846
#define SECOND_IN_JDAYS 1.157401129603386e-05; // 1 second expressed in JDAYS. JDAY/s
#define JDAY_IN_SECONDS 86400.46863810098; // 1 JDAY expressed in seconds. s/JDAY

class ExternalSpaceObject{
private:
	std::vector< std::vector<double> > InertialToRTC(std::vector<double>* positionPtr, std::vector<double>* velocityPtr);
	void LoadEphemerisTable( std::string fullFileName );

	std::vector<double> rhsX, rhsY, rhsZ; // Right-hand sides of the equations to be used to interpolate the ephemeris table.
	int NO_INTERPOLATION_POINTS; // Number of points to be used to interpolate the ephemeris table from, currently fixed to two.
	std::vector<double> interpolationArguments;
	std::vector<int> luDecompositionIndices; // Indices of the rows that have been swapped while doing the LU decomposition (needed for backsubstitution).
	double luDecompositionNoSwaps; // Number of row swaps done during the LU decomposition.
	std::vector< std::vector<double> > A; // A matrix of coefficients to solve the system of equations for interpolation in dimensionless time.
	std::vector< std::vector<double> > LU; // L-U decomposed matrix A.

public:
	/* Methods. */
	ExternalSpaceObject(void);
	ExternalSpaceObject(std::string SSC, std::string ephemerisFileName, int COV=1, int wgsEllipsoid=84, double objectRadius=5.0);
	~ExternalSpaceObject(void);

	void PropagateJDAY(std::vector<double>* posPtr, std::vector<double>* vloPtr, double JDAY, bool updateCurrentState=false);
	void Propagate(std::vector<double>* posPtr, std::vector<double>* veloPtr, int year, int month, int day, int hour, int minute, double second,  bool updateCurrentState=false);
	void CalculateRadius(std::vector<double>* posPtr);
	void SetHardBodyRadius(double bodyRadius);
	void SetCovarianceMatrixRTC(std::vector< std::vector<double> > CovarianceMatrixRTC);

	/* Propagation-sepcific attrributes. */
	double tumin, mu, Re, xke, j2, j3, j4, j3oj2; // Gravity model properties.
	gravconsttype gravityModel; // Gravity model being used - WGS84, WGS72 or old WGS72.
	std::vector<double> EphemerisEpochs; // Julian Day epochs of all the ephemeris points of this object.
	std::vector< std::vector<double> > EphemerisPositions; // 3 Cartesian position components in the TEME reference frame of all the ephemeris points of this object, same order as EphemerisEpochs.
	std::vector< std::vector<double> > EphemerisVelocities; // 3 Cartesian velocity components in the TEME reference frame of all the ephemeris points of this object, same order as EphemerisEpochs.

	/* General orbital elements and other attributes - may be updated as propagation progresses but initialised with the values at TLE epoch. */
	std::vector<double> currentPos, currentVelo; // Current position and velocity.
	double CURRENT_RADIUS, currentEpochJDAY; // Corresponding osculating radius and epoch.
	double hardBodyRadius; // Physical radius of the object (used to compute the collision probability) and epoch of the TLE that was used to create it.
	int earthEllipsoid; // Year of the WGS ellipsoid used to propgate the object, 72 or 84.
	std::string NORAD_ID; // SSC/NORAD ID number of the object, may be empty if the object isn't in the catalgoue.
	std::string EPHEMERIS_FILE_NAME; // Name of the ephemeris file that holds the ephemeris table for this object.
	std::vector<double> r0, v0; // Cartesian position and velocity in TEME at the object's epoch.
	double EpochJDAY; // Epoch of this object's first recorded ephemeris point, Julian Days.

	/* Covariance matrix related attributes. */
	std::vector< std::vector<double> > FullCovarianceMatrixRTC; // 6x6 position and velocity covariance matrix in km and km/sec.
	std::vector< std::vector<double> > PositionCovarianceMatrixRTC; // 3x3 position covariance matrix in km.

};

/*
*---------------------------------------------------------------------
* HELPER FUNCTIONS FOR EQUATIONS' SOLVING.
*---------------------------------------------------------------------
*/
std::vector< std::vector<double> > luDecomposition(std::vector< std::vector<double> > A, std::vector<int> &index, double &d);
void luSubst(std::vector < std::vector<double> >* A, std::vector<int>* index, std::vector<double>* b);

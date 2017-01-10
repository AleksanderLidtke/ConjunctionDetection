/*
A class that provides a convenient interface to the SGP4 propagator.
Based on code by David Vallado freely available on clestrak.com (last accessed 03/02/2014).

@version 1.4.0
@since 09 Dec 2016 12:20:00
@author Aleksander Lidtke
@email al11g09@soton.ac.uk, alek_l@onet.eu

CHANGELOG:
23 Jul 2014 - 1.1.0 - Added the covariance parameter that enables a fixed covariance matrix to be set.
24 Jul 2014 - 1.1.1 - Fixed a bug where the NORAD_ID attribute would be changed when reading three line elements hence causing the entire object creation framework to fail.
25 Jul 2014 - 1.1.2 - Added the capability to read three-line elements from a directory.
29 Jul 2014 - 1.2.0 - Begun to assume position covariance for object for which only one TLE is available and hence cannot do this using V.P. Osweiler's approach. Used standard deviations for various orbital regimes from Flohrer, T., Assessment and categorisation of TLE orbit errors for the US SSN catalogue, 1 Jan 2008 snapshot.
29 Jul 2014 - 1.2.1 - Started drafting the method that spawns TLEs when setting a fixed covariance for an object in order to enable the covariance to be propagated in a Monte Carlo fashion.
					 Still need code that will create a TLE from a state vector and velocity covariance matrix to create Monte Carlo TLEs.
30 Jul 2014 - 1.2.2 - Added filtering based on orbital energy to remove outlying TLEs when reading them from hard drive to estimate the covariance. This is based on approach from Patera, R.P., Space Event Detection Method, 2008.
31 Jul 2014 - 1.2.3 - Fixed an error that would make compilation to fail under g++ 4.8.1 regarding the isnan function used in ReadThreeLineElements method of SpaceObject.
04 Aug 2014 - 1.2.4 - Begun to throw an exception when don't have enough TLEs to estimate the covariance for a given object - it may be temporary for AMOS not to obsure the effects of growing covariance.
05 Aug 2014 - 1.3.0 - Added the capability to simulate the covariance matrices that would likely have been produced had the TLE been generated using the European Space Surveillance System.
11 Aug 2014 - 1.3.1 - Updated the covariances in COV 3 mode to the ones given by T. Flohrer.
26 Aug 2014 - 1.3.2 - Increased MIN_TLES to 5 from to following the discussion with Roberto Armellin to have more reliable statitics of the covariance matrix estimates.
23 Oct 2014 - 1.3.3 - added namespaces VectorOperations amd EquationsSolving.
02 Dec 2014 - 1.3.4 - Fixed a bug where the transverseUnit_inertial in inertialToRTC would be calculated by crossing XxZ unit vectors instead of ZxX.
09 May 2016 - 1.4.0 - Added the option to set the position variances through command line arguments to facilitate catalogue accuracy studies.
*/
#pragma once

#include <iostream>
#include <stdio.h>
#include <cstring>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <random>

#include "sgp4ext.h"
#include "sgp4unit.h"
#include "VectorOperations.h"

#define PI 3.14159265358979323846
#define MIN_TLES 5 // Minimum number of TLEs that should be available to estimate the covariance using V.P. Osweiler's method.

class SpaceObject{
private:
	void ReadThreeLineElements( std::string CovarianceThreeLineElementsDirectory="./Covariance 3LEs/", double rejectionThreshold = 3.0 );
	std::vector< std::vector<double> > InertialToRTC(std::vector<double>* positionPtr, std::vector<double>* velocityPtr);
public:
	/* Methods. */
	SpaceObject(void);
	SpaceObject(std::vector<std::string> threeLineElement, int COV=1, std::string covpath="./Covariance 3LEs", double bstarMultiplier=1.0, int wgsEllipsoid=84, double objectRadius=5.0, char opsMode='i');
	~SpaceObject(void);

	void PropagateJDAY(std::vector<double>* posPtr, std::vector<double>* vloPtr, double JDAY, bool updateCurrentState=false);
	void Propagate(std::vector<double>* posPtr, std::vector<double>* veloPtr, int year, int month, int day, int hour, int minute, double second,  bool updateCurrentState=false);
	void CalculatePerigeeRadius(void);
	void CalculateApogeeRadius(void);
	void SetHardBodyRadius(double bodyRadius);
	void SetPositionVariances(double radialVariance, double inTrackVariance, double crossTrackVariance);
	void SetCovarianceMatrixRTC(std::vector< std::vector<double> > CovarianceMatrixRTC, int noTLEs=10);
	void ComputeCovarianceOSW( double epochJDAY );

	/* Propagation-sepcific attrributes. */
	char OPSMODE; // Improved operation mode of SGP4, should result in a smoother behaviour.
	double tumin, mu, Re, xke, j2, j3, j4, j3oj2; // Gravity model properties.
	gravconsttype gravityModel; // Gravity model being used - WGS84, WGS72 or old WGS72.
	elsetrec sgp4Sat; // SGP4 propagation object.
	char classification; // Object classification (all here will be U - unclassified.
	char intldesg[11]; // International designator - launch year, launch number, launch piece.
	double BstarMultiplier; // Multiplication factor of the TLEs' B*, used to simulate variations in solar activity.

	/* SGP4-specific attributes. */
	int cardnumb, numb;
    long revnum, elnum;
	int nexp, ibexp;

	/* General orbital elements and other attributes - may be updated as propagation progresses but initialised with the values at TLE epoch. */
	std::vector<double> currentPos, currentVelo; // Current position and velocity.
	double CURRENT_PERIGEE_RADIUS, CURRENT_APOGEE_RADIUS, currentEpochJDAY; // Corresponding radii and epoch.
	double hardBodyRadius, TLEepochJDAY; // Physical radius of the object (used to compute the collision probability) and epoch of the TLE that was used to create it.
	int earthEllipsoid; // Year of the WGS ellipsoid used to propgate the object, 72 or 84.
	std::string NORAD_ID; // SSC/NORAD ID number of the object.
	double SemiMajorAxis, Eccentricity, Inclination, MeanAnomaly, LongAscendingNode, ArgumentOfPerigee, AltitudeOfPerigee; // Classical orbital elements and altitudes at the epoch of the TLE.
	double r0[3]; double v0[3]; // Cartesian position and velocity in TEME at the TLE's epoch.

	/* Covariance matrix related attributes. */
	std::vector< std::vector<double> > FullCovarianceMatrixRTC; // 6x6 position and velocity covariance matrix in km and km/sec.
	std::vector< std::vector<double> > PositionCovarianceMatrixRTC; // 3x3 position covariance matrix in km.
	std::vector<elsetrec> CovarianceTLEs; // TLEs that are representative of the covariance matrix and are used to propagate it in a Monte Carlo fashion.

};

std::vector< std::vector<double> > luDecomposition(std::vector< std::vector<double> > A, std::vector<int> &index, double &d);
void luSubst(std::vector < std::vector<double> >* A, std::vector<int>* index, std::vector<double>* b);


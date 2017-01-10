/*
A class that provides a framework for looking for conjunctions between objects in Earth orbit. The analysis is
done on a pair-wise basis, i.e. pairs of objects are created from the pool of objects of interest read from a
three-line element file and conjunction between each pair are found.
Various types of analyses are supported (one object against multiple ones, several against many or all-on-all).
Multiple convenience methods for reading input, performing the analyses and writing output are included.

@version 1.8.1
@since 09 May 2016 19:05:00
@author Aleksander Lidtke
@email al11g09@soton.ac.uk, alek_l@onet.eu

CHANGELOG:
23 Apr 2014 - 1.1.0 - merged conjunction detection with object pairs' creation.
14 May 2014 - 1.2.0 - added a fix to deallocate temporary position vectors.
22 May 2014 - 1.2.1 - removed computing the relative positions at TCA to increase speed and displaying all the found conjunctions to make the logs more legible.
11 Jun 2014 - 1.3.0 - added the skipping of certain analysis steps in case of XYZ prefilter failure as this has been shown not to affect the overall precision.
23 Jul 2014 - 1.4.0 - Added the covariance parameter that enables a fixed covariance matrix to be set for space objects.
31 Jul 2014 - 1.4.1 - Removed the fine search code as it will never be used - the code is accurate as it is and adding additional steps to increase the accuracy would make the execution too slow.
02 Aug 2014 - 1.4.2 - Fixed a bug where the progress information would not be displayed in one-on-all analysis mode.
03 Aug 2014 - 1.5.0 - Begun not to record conjunctions with erroneous collision probabilities in covariance mode 2 (estimation using V.P. Oswieler's method).
04 Aug 2014 - 1.5.1 - Added the capability to set a radius for the primary object when using external ephemerides in support of mission analysis at desgin stages.
	     	  1.5.2 - Added a catch block around calculateMaximumSphericalProbability for COV 1 as a FatalError would sometimes be thrown when inverting matrices therein.
05 Aug 2014 - 1.6.0 - Added the capability to simulate the covariance matrices that would likely have been produced had the TLE been generated using the European Space Surveillance System.
23 Aug 2014 - 1.6.1 - Added the capability to enforce direct numerical integration from the command line thus removing the need to maintain different versions of the code.
23 Oct 2014 - 1.6.2 - Added namespaces VectorOperations and EquationsSolving.
08 Dec 2014 - 1.6.3 - Started to check whether an external radius file has been provided to this simulation.
	     	  1.6.4 - Started to output the title of the simulation when displaying this simulation's progress information.
29 Jan 2016 - 1.7.0 - Added the option to skip collision probability calculations.
09 May 2016 - 1.7.1 - Corrected some spelling mistakes in the docs.
 			- 1.8.0 - Added the option to set the position variances through command line arguments to facilitate catalogue accuracy studies.
 			- 1.8.1 - Implemented a fix that ensures that analyses can be run back to back without finding the same conjunctions twice due to slighlty overlapping analysis intervals.
*/
#pragma once

#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <time.h>
#include <map>
#include "sgp4unit.h"
#include "sgp4ext.h"
#include "FatalError.h"
#include "MathWarning.h"
#include "SpaceObject.h"
#include "Conjunction.h"
#include "VectorOperations.h"
#include "ExternalSpaceObject.h"

#define SECOND_IN_JDAYS 1.157401129603386e-05; // 1 second expressed in JDAYS. JDAY/s
#define JDAY_IN_SECONDS 86400.46863810098; // 1 JDAY expressed in seconds. s/JDAY

class Simulation{
	private:
		/*
		 *---------------------------------------------------------------------
		 * HELPER METHODS FOR INTERPOLATION AND CONJUNCTION DETECTION.
		 *---------------------------------------------------------------------
		 */
		void findConjunctionsBetweenObjectPair(std::pair<SpaceObject, SpaceObject>* objectPairPtr, std::vector<Conjunction>* conjunctionsPtr, double* analysisStart, double* analysisStop);
		void findConjunctionsBetweenObjectPair(std::pair<ExternalSpaceObject, SpaceObject>* objectPairPtr, std::vector<Conjunction>* conjunctionsPtr);
		double findTCA(std::vector< std::vector<double> >* relativeInterpolatingPolynomialCoefficients, double TCAguess, int maxiter=1000, double accuracy=1E-8);
		void findInterpolatingCoefficients(std::vector< std::vector<double> >* posesInterpPtr, std::vector< std::vector<double> >* velosInterpPtr, std::vector<double>* nodeEpochsPtr, std::vector< std::vector<double> >* interpCoeffsPtr);
		bool XYZpreFilter(std::vector<double>* posVecPtr, double* thresholdValuePtr);
		double getInterpolatedRelativeDistanceSquared(double s);
		double getInterpolatedRelativeRangeRateSquared(double s);
		double getInterpolatedRelativeAccelerationSquared(double s);

		/*
		 *---------------------------------------------------------------------
		 * ATTRIBUTES.
		 *---------------------------------------------------------------------
		 */

		/* Gravity model. */
		int wgsEllipsoidUsed; gravconsttype gravityModel;
		double tumin, GRAVITY_CONSTANT, MEAN_EARTH_RADIUS, xke, j2, j3, j4, j3oj2;

		/* Analysis settings. */
		std::string TITLE;
		double CONJUNCTION_THRESHOLD_DISTANCE; // km, all cases when objects get below this separation distance will be found.
		double CONJUNCTION_THRESHOLD_DISTANCE_SQUARED; // km^2, square of the threshold distance.
		double COARSE_TIME_STEP_S; // s, time step in which a number of interpolation points is created to generate an ephemeris table from which to look for conjunctions.
		double COARSE_TIME_STEP; // Julian days, time step in which a number of interpolation points is created to generate an ephemeris table from which to look for conjunctions.
		int NO_INTERPOLATION_POINTS; // Number of pooints used to interpolate trajectories.
		double PERIGEE_APOGEE_PAD; // km, a buffer to be used in perigee/apogee pre-filtering.
		double PERIGEE_APOGEE_ORBITAL_REGIME_PAD; // km, a buffer that is used to filter out object that are in different orbital regimes and hence will never conjunct.
		double bStarScalingFactor; // Used to tweak B* to emulate varying solar conditions.
		double defaultR_RB, defaultR_PL, defaultR_DEB, defaultR_Other; // Default radii in metres for objects no in the database.
		double defaultVarianceRad, defaultVarianceIn, defaultVarianceCross; // Default position variances in m^2 that we can use to simulate varying catalogue accuracy.
		time_t TIME_STARTED, TIME_FINISHED;
		int CovarianceType; // Type of covariance to be used in collision probability calculations. 1 for spherical maximum spherical collision probability, 2 for using estimates of the covariance matrix based on historic TLEs for every object and computing true and maximum collision probabilities based on those.
		bool EnforceDirectNumericalIntegration; // Whether to use direct numerical integration when estimating the collision probability.

		/* Pre-filters' settings. */
		double SURFACE_GRAVITY_ACCELERATION; // km/s^2
		double ESCAPE_VELOCITY; // km/s
		double THRESHOLD_RADIUS;
		double THRESHOLD_RADIUS_SQUARED; // Add a safety factor of 2.0 in order not to miss conjunctions.
		std::vector<double> DISTANCE_THRESHOLD;
		double ACCELERATION_RADIUS;
		double ACCELERATION_RADIUS_SQUARED;

		/* Interpolation parameters. */
		int polynomialOrder; // Order of the interpolating polynomial.
		std::vector<double> interpolationArguments; // Dimensionless "time" [0,1] in the interpolation interval.
		std::vector< std::vector<double> > A; // Coefficients that are used to compute the interpolating poilynomial (interpolation of each Cartesian coordinate w.r.t. time).
											  // The order of interpolating polynomials' coefficients is a_0, a_1, ..., a_N for a polynomial of degree N.
		std::vector<double> rhsX; // Right hand sides of the sets of equations to be solved to find the interpolation coefficients. The interpolating polynomial coefficients in every coordianate themselves will also be stored here temporarily.
		std::vector<double> rhsY;
		std::vector<double> rhsZ;
		std::vector<int> luDecompositionIndices; // Indices of the rows that have been swapped while doing the LU decomposition (needed for backsubstitution).
		double luDecompositionNoSwaps; // Number of row swaps done during the LU decomposition.
		std::vector< std::vector<double> > LU; // A LU decomposed matrix A.
		double m6, m5, m4, m3, m2, m1, m0; // Coefficients of the relative distance, velcoty and acceleration functions (a 6th order polynomial for cubic interpolation).

	public:
		/*
		 *--------------------------------------------------------------------- 
		 * CONSTRUCTORS AND DECONSTRUCTORS.
		 *---------------------------------------------------------------------
		 */
		Simulation(void);
		Simulation(std::string title, int wgsEllipsoidType, double thresholdDistance=50., double coarseTimeStep=300., int noInterpolationPoints=2, double perigeeApogeeFilterPad=50.);
		~Simulation(void);

		/*
		 *---------------------------------------------------------------------
		 * PRIMARY METHODS.
		 *---------------------------------------------------------------------
		 */
		void setDirectNumericalIntegration( bool flag );
		void performAnalysis(int startYear, int startMonth, int startDay, int startHour, int startMinute, double startSecond, int stopYear, int stopMonth, int stopDay, int stopHour, int stopMinute, double stopSecond, int mode, std::vector<std::string> primariesSSCs);
		void performAnalysis(double ANALYSIS_INTERVAL_START, double ANALYSIS_INTERVAL_STOP, int mode, std::vector<std::string> primariesSSCs);
		void performAnalysis(std::string primariesSSCs, std::string primaryEphemerisFileName, double primaryRadius);
		void createObjectsFromTLEfile(const char* fileName, double defaultRadius_RB, double defaultRadius_PL, double defaultRadius_DEB, double defaultRadius_Other, int COV, std::string covpath, double defaultVarRad=20.0, double defaultVarIn=100.0, double defaultVarCross=20.0, double bStarMultiplicationFactor=1.0, const char * radiusFileName="stkRadiusFile.rad");
		void saveConjunctions( const char * outFileName, std::string version);

		/*
		 *---------------------------------------------------------------------
		 * ATTRIBUTES.
		 *---------------------------------------------------------------------
		 */

		/* Objects and conjunctions. */
		std::map<std::string, SpaceObject>* objectsPtr; // Objects included in the analysis.
		std::vector<Conjunction>* conjunctionsFoundPtr; // Conjunctions that have been found in this simulation.
};

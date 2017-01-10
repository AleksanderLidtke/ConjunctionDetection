/*
A class that provides a convenient means of recording, analysing and handling orbital conjunctions.

@version 1.6.2
@since 11 Mar 2015 13:53:00
@author Aleksander Lidtke
@email al11g09@soton.ac.uk, alek_l@onet.eu

CHANGELOG:
28 Jul 2014 - 1.1.0 - Finished verifying calculateCollisionProbability and calculateMaximumCollisionProbability. calculateMaximumCollisionProbability has to use
	numerical integration of the PDF with at least 1000 terms in order to be reasonably accurate, so commented out the series expansion's implementation.
03 Aug 2014 - 1.1.1 - Added an iterations cap when looking for maximum collision probability with fixed covariance as sometimes the desired tolerance would never be reached.
04 Aug 2014 - 1.1.2 - Added an iterations cap when looking for maximum spherical probability as well. In some rare cases the iterations would continue until double floating point precision would be exceeded.
	This issue would, interestingly, only reveal itself when compiling in Release configuration.
23 Aug 2014 - 1.2.0 - Added parameters to calculate collision probability methods that enforce numerical integration to be used.
26 Aug 2014 - 1.3.0 - Returned the functionality to caluclate maximum collision probability using the series expansion, if desired.
26 Aug 2014 - 1.3.1 - Started using covarianceScalingFactorSquared in CalculateMaximumCollisionProbability for speed and increased default numer of terms to be used in Chan's series expansion to 1000 from 5.
29 Aug 2014 - 1.3.2 - Stopped throwing a FatalError when the maximum collision probability exceeded 1.0 as this is attributed to extremely close miss distances 
	that are caused by erroneous state vectors of objects that have re-entered and will be filtered out in post-processin.
29 Aug 2014 - 1.3.3 - Changed the number of Chan's terms back to 5 to avoid haiving NaNs in probabilities - the terms will get smaller and smaller the more there are and eventually become NaNs.
23 Oct 2014 - 1.3.4 - Added VectorOperations namespace everywhere.
12 Nov 2013 - 1.3.5 - Started using Simpson's rule to do the numerical PDF integration to get Pc.
19 Nov 2014 - 1.3.6 - Updated calculateMaximumSphericalProbability to match the Pc|MAX algorithm developed for AMOS (changed the while loop criterion and started looking at the NaNs).
		    - Moved hZ* in integrateProbabilityPDF from the for loop in k to the end to reduce the number of floating point operations.
02 Dec 2014 - 1.3.7 - Fixed a bug where the transverseUnit_inertial in inertialToRTC would be calculated by crossing XxZ unit vectors instead of ZxX.
04 Jan 2015 - 1.3.8 - Changed the golden-ratio search in calculateMaximumCollisionProbability to look for the scaling factor rather than scaling factor squared in hope that this will make it converge better.
05 Jan 2015 - 1.3.8 - Corrected a bug in calculateMaximumCollisionProbability where the covarianceScalingFactor would be uninitialised thus breaking the entire golden ratio search and causing the Pc|MAX to be erroneous.
28 Jan 2015 - 1.4.0 - Started to use diagonalised covariance matrix in calculateMaximumCollisionProbability even when using direct numerical integration to be consistent with the series expansion algorithm. It's OK to use non-diagonalised matrix in calculateCollisionProbability though.
					- calculateCollisionProbability will now throw a FatalError when unsuccessfully inversing the covariance matrix instead of trying to diagonalise it first.
29 Jan 2015 - 1.5.0 - Corrected errors in integrateProbabilityPDFTrapeze and integrateProbabilityPDFSimpson where only quarter of a circle was used to integrate PDF even though it isn't constant over that circle. This would result in collision probability to be higher for objects whose in-plane position was in the same quarter of the B-plane as this quarter of the circle.
11 Mar 2015 - 1.6.0 - Added a dedicated function to calculate collision probability using the series expansion by K. Chan.
			- 1.6.1 - cleared some unused variables and updated the docs.
			- 1.6.2 - Added a golden ratio search when looking for PcMAX using series expansion too.
*/
#pragma once

#include <cmath>
#include <vector>
#include <string>
#include <float.h>
#include <sstream>
#include <iostream>
#include "sgp4ext.h"
#include "dsyevv3.h"
#include "VectorOperations.h"

class Conjunction{
	private:

	public:
		/* Methods. */
		Conjunction(void);
		Conjunction(double TCA_JDAY, double MD_km, double relativeV, double collisionRadius, std::string object1NORAD_ID, std::string object2NORAD_ID);
		~Conjunction(void);
		std::string getReadableTCA( void );
		void calculateMaximumSphericalCollisionProbability( bool EnforceNumericalIntegration=false, int noTerms=5 );
		void calculateCollisionProbability( std::vector<double>* position1, std::vector<double>* velocity1, std::vector< std::vector<double> >* RTCcovMat1, std::vector<double>* position2, std::vector<double>* velocity2, std::vector< std::vector<double> >* RTCcovMat2, bool EnforceNumericalIntegration=false, int noTerms=5 );
		void calculateMaximumCollisionProbability( std::vector<double>* position1, std::vector<double>* velocity1, std::vector< std::vector<double> >* RTCcovMat1, std::vector<double>* position2, std::vector<double>* velocity2, std::vector< std::vector<double> >* RTCcovMat2, bool EnforceNumericalIntegration=false, int noTerms=5 );
		
		/* Declare a friend to be able to write to stdout neatly. */
		friend std::ostream& operator<<(std::ostream& os, const Conjunction& conj);

		/* Attributes. */
		double missDistance, TCA, relativeVelocity, trueProbability, maximumProbability, combinedCollisionRadius;
		std::string primaryNORAD_ID, secondaryNORAD_ID, readableTCA;
};

/* Helper functions. */
void shift2(double &a, double &b, const double c);
void shift3(double &a, double &b, double &c, const double d);
double integrateProbabilityPDFSeries(double combinedCollisionRadius, std::vector<double>* relativeProjectedPosition, std::vector< std::vector<double> >* Pd, int noTerms=50);
double integrateProbabilityPDFTrapeze( double combinedCollisionRadius, std::vector<double>* relativeProjectedPosition, std::vector< std::vector<double> >* Pd, std::vector<std::vector<double> >* PdInverse, int noIntervals=5000);
double integrateProbabilityPDFSimpson( double combinedCollisionRadius, std::vector<double>* relativeProjectedPosition, std::vector< std::vector<double> >* Pd, std::vector< std::vector<double> >* PdInverse, int noIntervals=5000);
std::vector< std::vector<double> > InertialToRTC(std::vector<double>* positionPtr, std::vector<double>* velocityPtr);
std::vector<double> EigenValues2x2( std::vector<std::vector<double> >* matPtr);

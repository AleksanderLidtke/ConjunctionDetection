#include "Conjunction.h"

Conjunction::Conjunction(void){
	/* Default constructor, does nothing special. */
};

Conjunction::Conjunction(double TCA_JDAY, double MD_km, double relativeV, double collisionRadius, std::string object1NORAD_ID, std::string object2NORAD_ID){
	/* Create an orbital conjunction object that stores all the relevant information, such as:
	@param TCA_JDAY - Time of Closest Approach in Julian days.
	@param MD_km - Miss Distance in km.
	@param relativeVelocity - relative velocity at TCA in km per second.
	@param collisionRadius - the combined radius of both conjuncting objects in km.
	@param object1NORAD_ID, object2NORAD_ID - NORAD IDs of both objects involved.
	*/
	missDistance=MD_km; TCA=TCA_JDAY; relativeVelocity=relativeV;  combinedCollisionRadius=collisionRadius; // Save the basic information.
	primaryNORAD_ID=object1NORAD_ID; secondaryNORAD_ID=object2NORAD_ID;

	readableTCA = std::string(); // Convert the TCA from Julian days to year/month/day hour:minute:second format.
	readableTCA = getReadableTCA();
	maximumProbability = 0.0; trueProbability = -1.0; // Won't be able to compute true collision probability with no covariance, need to initialise it with some characteristic value.
};

std::string Conjunction::getReadableTCA(void){
	/* Convert the epoch of closest approach to a readable format year/month/day hour:minute:second
	and write it to the outStr it.
	@param outStrPtr - TCA epoch will be wrriten to this string.
	*/
	int year, mon, day, hr, min; // Get the TCA epoch in a readable format.
	double sec;
	invjday(TCA, year, mon, day, hr, min, sec);
	std::stringstream ss;
	ss<< day<<"/"<<mon<<"/"<<year<<" "<<hr<<":"<<min<<":"<<sec; // Convert the TCA to a string.
	return ss.str();
};

void Conjunction::calculateMaximumSphericalCollisionProbability( bool EnforceNumericalIntegration, int noTerms ){
	/* Compute the maximum probability of collision between the two objects at the time of closest approach. 
	In order to limit the amount of computation that will be done the miss distance class attribute will be used
	to do this as the actual relative position and velocity has negligible effect on the maximum collision probability 
	when spherical uncertainty is assumed. 
	The maximum probability as given in N. Berend (1999) will be computed assuming a spherical or ellipsoidal uncertainty 
    distribution of both objects depending on the spherical flag. This may be further numerically 
    optimised to ensure that the true maximum probaility is computed. A 1:1:1 (RTC) aspect ratio of the 
    covariance matrices and uncorrelated errors in every dimension will be assumed.
    If the ratio of the combined objects' radius to the nominal miss distance is less or equal to 80% the probability will be computed
    by a series expansion of the PDF integral proposed by K. Chen (2008). This has been found accurate enough and not to
    require further numerical optimisation (the analytical estimate of the highest collision probability is accurate within several percent). 
    Otherwise the relative position PDF will be numerically integrated and optimised to yield a true worst case collision probability (can reach 1.0).
	Numerical integration of the PDF can be enforced in all cases by setting EnforceNumericalIntegration=true.
	The collision probability estiamtion approach assumes linear relative trajectory, hence the miss distance that has been
	provided to the contructor MUST BE those of the closest approach to make this assumption viable.
	@param EnforceNumericalIntegration - whether to enforce the numerical integration to be used instead of Chan's series expansion.
	@param noTerms - number of terms that will be used to compute the collision probability if the series expansion is used.
	@throws MathWarning - when the probability cannot be calculated due to uninversible covariance matrices.
	@return - the found maximum collision probability can be accessed through this.maximumProbability attribute.
	@deprrecated - use one of the newer methods that use the actual aspect ratio of the covariance matrix. Kept for compatibility with old simulations.
	*/
	std::vector<double> relativeProjectedPosition = std::vector<double>(2, 0.0);    
	// Synthesise the relative geometry to limit the amount of calculation without fidelity loss - the maximum probability for spherical uncertainties depends only on the miss distance, not the relative positions and velocities.
	relativeProjectedPosition.at(1) = missDistance;

    /* FIND THE FIRS APROXIMATE COMBINED COVARIANCE MATRIX TO YIELD WORST-CASE COLLISION PROBABILITY */
    double covarianceScalingFactor =  VectorOperations::vectorMagnitude(&relativeProjectedPosition)/std::sqrt(2.0); // Standard deviation which yields the maximum collision probability for spherical uncertainty distribution.
    std::vector< std::vector<double> > projectedCovarianceMatrixInPlane = std::vector< std::vector<double> >(2, std::vector<double>(2, 0.0) );
	projectedCovarianceMatrixInPlane.at(0).at(0) = 1.0; projectedCovarianceMatrixInPlane.at(1).at(1) = 1.0; // When looking for the worst-case can assume zero cross-term.
	std::vector< std::vector<double> > Pd = VectorOperations::matrixMultiplyByScalar( &projectedCovarianceMatrixInPlane, covarianceScalingFactor*covarianceScalingFactor ); // The worst-case covariance matrix in the collision plane. 
			
    /* COMPUTE THE COLLISION PROBABILITY */
    if( ( (combinedCollisionRadius > 0.8*VectorOperations::vectorMagnitude(&relativeProjectedPosition)) && (covarianceScalingFactor > DBL_MIN ) ) || EnforceNumericalIntegration ){ // Check whether to use direct integration or series expansion. In rare cases covariance scaling facotr will be 0, won't be able to inverse the covariance matrix then...
		// Use direct integration. Need to find the true worst-case Pd as in this parameter regime the analytcial estimate doesn't work very well.
		
		/* Do a golden ratio search for the true worst-case covariance scaling factor. */
		double TOL = 0.01; // Tolerance to within which the maximum will be found.
		int MAX_ITER = 100; // Maximum number of iterations allowed.
		int currentIter = 0; // Number of iterations done already.
		double fTEMP = 0.0; // -1*Current best estimate of Pc.

		std::vector< std::vector<double> > PdInverse = VectorOperations::inverseTwoByTwo( &Pd );
		
		const double R=0.61803399, C=1.0-R; // The golden ratios.
		double f1,f2,x0,x1,x2,x3; // At any time keep track of four points x0,x1,x2 and x3.
								  // Do a golden search for ax < bx < cx and f(bx) < f(ax) && f(bx) < f(cx).

		double ax = 0.0; // Starting points that bracket the root (normally less than covarianceScalingFactor). Sometimes it's 0.
		double bx = 0.5*covarianceScalingFactor;
		double cx = 1.5*covarianceScalingFactor;
		
		x0=ax; x3=cx;
		if( fabs(bx-cx) > fabs(cx-ax) ){ // Make x0 to x1 the smaller segment.
			x1=bx; // And fill in the new point to be tried.
			x2=bx+C*(cx-bx);
		}else{
			x2=bx;
			x1=bx-C*(bx-ax);
		};

		Pd = VectorOperations::matrixMultiplyByScalar( &projectedCovarianceMatrixInPlane, x1*x1 ); // Initial endpoint function evaluations.
		PdInverse = VectorOperations::inverseTwoByTwo( &Pd );
		f1 = -1.0*integrateProbabilityPDFSimpson( combinedCollisionRadius, &relativeProjectedPosition, &Pd, &PdInverse, 5000);

		Pd = VectorOperations::matrixMultiplyByScalar( &projectedCovarianceMatrixInPlane, x2*x2 );
		PdInverse = VectorOperations::inverseTwoByTwo( &Pd );
		f2 = -1.0*integrateProbabilityPDFSimpson( combinedCollisionRadius, &relativeProjectedPosition, &Pd, &PdInverse, 5000);

		// See if the probabilities are numbers and keep track of that, for ill-conditioned covariance matrices they'll be NaNs.
		#ifdef _MSC_VER // MSVS has a different _isnan function from <float.h>, Linux compilers use std::isnan from <cmath>
			int IS_NAN_f1, IS_NAN_f2;
			IS_NAN_f1 = _isnan(f1);
			IS_NAN_f2 = _isnan(f2);
		#else
			bool IS_NAN_f1, IS_NAN_f2;
			IS_NAN_f1 = std::isnan(f1);
			IS_NAN_f2 = std::isnan(f2);
		#endif

		while( (fabs(x3-x0) > TOL*( fabs(x1)+fabs(x2) )) && !IS_NAN_f1 && !IS_NAN_f2 ){
			if(f2 < f1){ // One possible outcome of the golden ratio root braceting.
				shift3(x0,x1,x2,R*x2+C*x3);
				Pd = VectorOperations::matrixMultiplyByScalar( &projectedCovarianceMatrixInPlane, x2*x2 );
				PdInverse = VectorOperations::inverseTwoByTwo( &Pd );
				fTEMP = -1.0*integrateProbabilityPDFSimpson( combinedCollisionRadius, &relativeProjectedPosition, &Pd, &PdInverse, 5000);
				shift2(f1,f2,fTEMP);
				currentIter += 1;
			}else{ // The other outcome.
				shift3(x3,x2,x1,R*x1+C*x0);
				Pd = VectorOperations::matrixMultiplyByScalar( &projectedCovarianceMatrixInPlane, x1*x1 );
				PdInverse = VectorOperations::inverseTwoByTwo( &Pd );
				fTEMP = -1.0*integrateProbabilityPDFSimpson( combinedCollisionRadius, &relativeProjectedPosition, &Pd, &PdInverse, 5000);
				shift2(f2,f1,fTEMP);
				currentIter += 1;
			};

			if(currentIter==MAX_ITER){
				std::cerr<<"Terminated iterations when looking for true maximum probability, current accuracy is "<<fabs(x3-x0)<<" and the desired "<<TOL*( fabs(x1)+fabs(x2) )<<". Using current best Pc_MAX estimate."<<std::endl;
				break;
			}
		}; //END while searching for the worst-case covariance scaling factor

		if( (f1 < f2) && !IS_NAN_f1){ // Output the best of the two current values.
			covarianceScalingFactor = x1;
			maximumProbability = -1.0*f1; // Scale to the actual probability value (looking for a minimum here).
		}else if( !IS_NAN_f2 ){
			covarianceScalingFactor = x2;
			maximumProbability = -1.0*f2;
		};

	}else if(covarianceScalingFactor < DBL_MIN ){ // Set probability to 1.0, there is 0 miss distance.
		maximumProbability = 1.0;
	}else{ // Use series expansion. No need to numerically find the worst-case covariance matrix here - in this regime the analytical guess for Pd is accurate enough.
		/* 11 mar 2015 - do the golden ratio search anyway to be able to prove that this works in the thesis. Otherwise run into risk that some super close conjunctions have higher PcMAX than the medioka ones. */
        /* Rotate the covariance matrix into the principal axes. */
		double rotationAngle;
		if( Pd.at(0).at(0)-Pd.at(1).at(1) <= DBL_MIN ){ // Both standard deviations are approximately equal, choose a fixed angle as given by K. Chan.
			rotationAngle = PI/4.0;
		}else{ // For cases where a specific covariance matrix will be specified, will never be true for diagonal matrices as used for maximum probability.
			rotationAngle = std::atan( 2.0*Pd.at(0).at(1)/(Pd.at(0).at(0) - Pd.at(1).at(1)) ); // Pd[0,1] is the cross term, Pd[0,0] and Pd[1,1] are the standard deviations squared.
		};
        
        std::vector< std::vector<double> >inPlaneRotationMatrix = std::vector< std::vector<double> >(2, std::vector<double>(2, 0.0) ); // Rotates the in-plane coordinate system into principal axes.
		inPlaneRotationMatrix.at(0).at(0)=std::cos(rotationAngle);inPlaneRotationMatrix.at(0).at(1)=-1.0*std::sin(rotationAngle);
		inPlaneRotationMatrix.at(1).at(0)=std::sin(rotationAngle);inPlaneRotationMatrix.at(1).at(1)=std::cos(rotationAngle);
        
		VectorOperations::vectorMultiplyByMatrix( &relativeProjectedPosition, &inPlaneRotationMatrix); // Rotate to principal axes.
        // Would need to find eigenvalues of Pd here if the matrix had a cross-term. But this is never the case for spherical probability.
        
		/* Do a golden ratio search for the true worst-case covariance scaling factor. */
		double TOL = 0.01; // Tolerance to within which the maximum will be found.
		int MAX_ITER = 100; // Maximum number of iterations allowed.
		int currentIter = 0; // Number of iterations done already.
		double fTEMP = 0.0; // -1*Current best estimate of Pc.

		std::vector< std::vector<double> > PdInverse = VectorOperations::inverseTwoByTwo( &Pd );
		
		const double R=0.61803399, C=1.0-R; // The golden ratios.
		double f1,f2,x0,x1,x2,x3; // At any time keep track of four points x0,x1,x2 and x3.
								  // Do a golden search for ax < bx < cx and f(bx) < f(ax) && f(bx) < f(cx).

		double ax = 0.0; // Starting points that bracket the root (normally less than covarianceScalingFactor). Sometimes it's 0.
		double bx = 0.5*covarianceScalingFactor;
		double cx = 1.5*covarianceScalingFactor;
		
		x0=ax; x3=cx;
		if( fabs(bx-cx) > fabs(cx-ax) ){ // Make x0 to x1 the smaller segment.
			x1=bx; // And fill in the new point to be tried.
			x2=bx+C*(cx-bx);
		}else{
			x2=bx;
			x1=bx-C*(bx-ax);
		};

		Pd = VectorOperations::matrixMultiplyByScalar( &projectedCovarianceMatrixInPlane, x1*x1 ); // Initial endpoint function evaluations.
		f1 = -1.0*integrateProbabilityPDFSeries(combinedCollisionRadius, &relativeProjectedPosition, &Pd);

		Pd = VectorOperations::matrixMultiplyByScalar( &projectedCovarianceMatrixInPlane, x2*x2 );
		f2 = -1.0*integrateProbabilityPDFSeries(combinedCollisionRadius, &relativeProjectedPosition, &Pd);

		// See if the probabilities are numbers and keep track of that, for ill-conditioned covariance matrices they'll be NaNs.
		#ifdef _MSC_VER // MSVS has a different _isnan function from <float.h>, Linux compilers use std::isnan from <cmath>
			int IS_NAN_f1, IS_NAN_f2;
			IS_NAN_f1 = _isnan(f1);
			IS_NAN_f2 = _isnan(f2);
		#else
			bool IS_NAN_f1, IS_NAN_f2;
			IS_NAN_f1 = std::isnan(f1);
			IS_NAN_f2 = std::isnan(f2);
		#endif

		while( (fabs(x3-x0) > TOL*( fabs(x1)+fabs(x2) )) && !IS_NAN_f1 && !IS_NAN_f2 ){
			if(f2 < f1){ // One possible outcome of the golden ratio root braceting.
				shift3(x0,x1,x2,R*x2+C*x3);
				Pd = VectorOperations::matrixMultiplyByScalar( &projectedCovarianceMatrixInPlane, x2*x2 );
				fTEMP = -1.0*integrateProbabilityPDFSeries(combinedCollisionRadius, &relativeProjectedPosition, &Pd);
				shift2(f1,f2,fTEMP);
				currentIter += 1;
			}else{ // The other outcome.
				shift3(x3,x2,x1,R*x1+C*x0);
				Pd = VectorOperations::matrixMultiplyByScalar( &projectedCovarianceMatrixInPlane, x1*x1 );
				fTEMP = -1.0*integrateProbabilityPDFSeries(combinedCollisionRadius, &relativeProjectedPosition, &Pd);
				shift2(f2,f1,fTEMP);
				currentIter += 1;
			};

			if(currentIter==MAX_ITER){
				std::cerr<<"Terminated iterations when looking for true maximum probability, current accuracy is "<<fabs(x3-x0)<<" and the desired "<<TOL*( fabs(x1)+fabs(x2) )<<". Using current best Pc_MAX estimate."<<std::endl;
				break;
			}
		}; //END while searching for the worst-case covariance scaling factor

		if( (f1 < f2) && !IS_NAN_f1){ // Output the best of the two current values.
			covarianceScalingFactor = x1;
			maximumProbability = -1.0*f1; // Scale to the actual probability value (looking for a minimum here).
		}else if( !IS_NAN_f2 ){
			covarianceScalingFactor = x2;
			maximumProbability = -1.0*f2;
		};
	};//END else (compute series-based probability).
};//END calculateMaximumCollisionProbability

void Conjunction::calculateCollisionProbability( std::vector<double>* position1, std::vector<double>* velocity1, std::vector< std::vector<double> >* RTCcovMat1, std::vector<double>* position2, std::vector<double>* velocity2, std::vector< std::vector<double> >* RTCcovMat2, bool EnforceNumericalIntegration, int noTerms ){
	/* Compute the true probability of collision between the two objects given their positions, velocities and position covariance matrices
    (must be in the same frame of reference and units) at the time of closest approach.
	If the ratio of the combined objects' radius to the nominal miss distance is less or equal to 80% the probability will be computed
    by a series expansion of the PDF integral proposed by K. Chen (2008). Otherwise, or if enforced through EnforceNumericalIntegration parameter,
	the relative position PDF will be numerically integrated.
    The collision probability estiamtion approach assumes linear relative trajectory, hence the provided positions and velocities
    MUST BE those of the closest approach to make this assumption viable.
    @param position1, position2 - positions of objects 1 and 2 in the form of 1x3 vectors in an inertial reference frame at the TCA.
    @param velocity1, velocity2 - velocities of those objects, also in the form of 1x3 vectors in an inertial reference frame at the TCA.
	@param RTCcovMat1, RTCcovMat2 - 3x3 position covariance matrices of the both objects at the TCA in the respective radial - in-track - cross-track frames of reference.
	@patam EnforceNumericalIntegration - whether to enforce the numerical integration to be used instead of Chan's series expansion.
	@param noTerms - number of terms that will be used to compute the collision probability if the series expansion is used.
	@throws FatalError - when the probability cannot be calculated due to uninversible covariance matrices or when one of the eigenvalues of the combined covariance matrix is zero.
	@return - the found maximum collision probability can be accessed through this.maximumProbability attribute.
	*/

	// Use a reference frame centred on object1 throughout the process (SCRF), the inputs are in the World Reference Frame (WRF, e.g. TEME or ECI) or RTC of respective objects for covariance matrices.
    
    /* RELATIVE POSITION AND VELOCITY VECTORS */
    std::vector<double> relPositionWRF = std::vector<double>(3,0.0); std::vector<double> xAxisVectorSCRF = std::vector<double>(3,0.0);
	VectorOperations::vectorDifference(&relPositionWRF, position1, position2); VectorOperations::vectorDifference(&xAxisVectorSCRF, velocity1, velocity2); // Object 2 in 1-fixed frame. Account for translation from WRF to SCRF.
    VectorOperations::unitVector(&xAxisVectorSCRF); // This is the X-axis of the relative velocity-based reference frame.
    
    /* CONSTRUCTION THE RELATIVE REFERENCE FRAME. */
    // Find the remaining components of S/C-centred relative velocity-aligned reference frame.
    std::vector<double> yAxisVectorSCRF = VectorOperations::vectorMultiplyByScalar( &xAxisVectorSCRF, VectorOperations::dotProduct(&xAxisVectorSCRF, &relPositionWRF) ); // Normal to X-axis and pointing towards the intersection of the collision plane with the relative velocity line.
	VectorOperations::vectorDifference( &yAxisVectorSCRF, &relPositionWRF, &yAxisVectorSCRF ); // relPositionWRF - the extra component along the X-axis = the component of relative position normal to the X-axis.

    std::vector<double> zAxisVectorSCRF = std::vector<double>(3,0.0);
	VectorOperations::crossProduct(&zAxisVectorSCRF, &yAxisVectorSCRF, &xAxisVectorSCRF);  // Normal to Y- and X- axes.
     
    // Non-dimensionalise the remaining components (in-plane).
	VectorOperations::unitVector( &yAxisVectorSCRF );
	VectorOperations::unitVector( &zAxisVectorSCRF );
    
    // Unit vectors of the World Reference Frame.
	std::vector<double> xAxisUnitVectorWRF = std::vector<double>(3,0.0); xAxisUnitVectorWRF.at(0) = 1.0;
    std::vector<double> yAxisUnitVectorWRF = std::vector<double>(3,0.0); yAxisUnitVectorWRF.at(1) = 1.0;
    std::vector<double> zAxisUnitVectorWRF = std::vector<double>(3,0.0); zAxisUnitVectorWRF.at(2) = 1.0;
    
    // Construct the rotation matrix from WRF to SCRF.
    std::vector< std::vector<double> > rotationMatrix = std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0) );
	rotationMatrix.at(0).at(0) = VectorOperations::dotProduct(&xAxisVectorSCRF, &xAxisUnitVectorWRF); rotationMatrix.at(0).at(1) = VectorOperations::dotProduct(&xAxisVectorSCRF, &yAxisUnitVectorWRF); rotationMatrix.at(0).at(2) = VectorOperations::dotProduct(&xAxisVectorSCRF, &zAxisUnitVectorWRF);
	rotationMatrix.at(1).at(0) = VectorOperations::dotProduct(&yAxisVectorSCRF, &xAxisUnitVectorWRF); rotationMatrix.at(1).at(1) = VectorOperations::dotProduct(&yAxisVectorSCRF, &yAxisUnitVectorWRF); rotationMatrix.at(1).at(2) = VectorOperations::dotProduct(&yAxisVectorSCRF, &zAxisUnitVectorWRF);
	rotationMatrix.at(2).at(0) = VectorOperations::dotProduct(&zAxisVectorSCRF, &xAxisUnitVectorWRF); rotationMatrix.at(2).at(1) = VectorOperations::dotProduct(&zAxisVectorSCRF, &yAxisUnitVectorWRF); rotationMatrix.at(2).at(2) = VectorOperations::dotProduct(&zAxisVectorSCRF, &zAxisUnitVectorWRF);
    
	std::vector< std::vector<double> > rotationMatrixInverse = VectorOperations::inverseThreeByThree(&rotationMatrix); // Will need these to rotate both covariance matrices to the B-plane.
	std::vector< std::vector<double> > rotationMatrixTranspose = VectorOperations::transposeMatrix(&rotationMatrix);
	//std::vector< std::vector<double> > rotationMatrixTransposeInverse = VectorOperations::inverseThreeByThree(&rotationMatrixTranspose);

    /* COMPUTE RELATIVE POSITIONS AND VELOCITIES */
    // Relative position expressed in the object1-fixed reference frame (SCRF).
    std::vector<double> temp = VectorOperations::vectorMultiplyByMatrix(&relPositionWRF, &rotationMatrix);
	std::vector<double> relativeProjectedPosition = std::vector<double>(2, 0.0);
	relativeProjectedPosition.at(0) = temp.at(1); relativeProjectedPosition.at(1) = temp.at(2); // Discard the velocity (X) component, only care about the components within the collision plane.
    
    /* FIND THE COMBINED COVARIANCE MATRIX IN THE B-PLANE. */
	std::vector< std::vector<double> > TEME2RTC_primary = InertialToRTC(position1, velocity1); // Rotation matrix between the TEME and RTC reference frames at the TCA for the primary.
	std::vector< std::vector<double> > TEME2RTC_secondary = InertialToRTC(position2, velocity2); // Rotation matrix between the TEME and RTC reference frames at the TCA for the secondary.

	// RTC to TEME rotation for the primary - formulate necessary matrices and then multiply out for rotation.
	std::vector< std::vector<double> > TEME2RTCinverse_primary = VectorOperations::inverseThreeByThree(&TEME2RTC_primary); // Need this to rotate from VNC to TEME.
	std::vector< std::vector<double> > TEME2RTCtranspose_primary = VectorOperations::transposeMatrix(&TEME2RTC_primary); // Also need the transpose when rotating covaraince matrices.
	std::vector< std::vector<double> > TEME2RTCtransposeInverse_primary = VectorOperations::inverseThreeByThree(&TEME2RTCtranspose_primary);
	
	std::vector< std::vector<double> > tempMat_primary = VectorOperations::matrixMultiplyByMatrix(RTCcovMat1, &TEME2RTCtransposeInverse_primary); // The operation is RotMat^-1 CovMat RotMat.T^-1.
	std::vector< std::vector<double> > posCovarianceTEME_primary = VectorOperations::matrixMultiplyByMatrix(&TEME2RTCinverse_primary, &tempMat_primary); // Finished rotations to TEME.

	// RTC to TEME rotation for the secondary - formulate necessary matrices and then multiply out for rotation.
	std::vector< std::vector<double> > TEME2RTCinverse_secondary = VectorOperations::inverseThreeByThree(&TEME2RTC_secondary);
	std::vector< std::vector<double> > TEME2RTCtranspose_secondary = VectorOperations::transposeMatrix(&TEME2RTC_secondary);
	std::vector< std::vector<double> > TEME2RTCtransposeInverse_secondary = VectorOperations::inverseThreeByThree(&TEME2RTCtranspose_secondary);
	
	std::vector< std::vector<double> > tempMat_secondary = VectorOperations::matrixMultiplyByMatrix(RTCcovMat2, &TEME2RTCtransposeInverse_secondary);
	std::vector< std::vector<double> > posCovarianceTEME_secondary = VectorOperations::matrixMultiplyByMatrix(&TEME2RTCinverse_secondary, &tempMat_secondary); // Finished rotations to TEME.

	// Add the matrices together in TEME - assume uncorrelated errors for both objects
	std::vector< std::vector<double> > combinedCovarianceMatrixTEME = VectorOperations::addMatrices( &posCovarianceTEME_primary, &posCovarianceTEME_secondary ); // Combined 3x3 position covariance defined in TEME at TCA.

	// Et viola, we have the combined covariance matrix rotated into the B-plane.
	std::vector< std::vector<double> > tempMat_combined = VectorOperations::matrixMultiplyByMatrix(&combinedCovarianceMatrixTEME, &rotationMatrixTranspose);
	std::vector< std::vector<double> > combinedCovarianceMatrixBplane = VectorOperations::matrixMultiplyByMatrix(&rotationMatrix, &tempMat_combined);

	std::vector< std::vector<double> > Pd = std::vector< std::vector<double> >(2, std::vector<double>(2, 0.0) ); // Covariance matrix defined in the B-plane, at the TCA and ignoring the velocity (X) components.
	Pd.at(0).at(0) = combinedCovarianceMatrixBplane.at(1).at(1); Pd.at(0).at(1) = combinedCovarianceMatrixBplane.at(1).at(2);
	Pd.at(1).at(0) = combinedCovarianceMatrixBplane.at(2).at(1); Pd.at(1).at(1) = combinedCovarianceMatrixBplane.at(2).at(2);

    /* COMPUTE THE COLLISION PROBABILITY */
    if( ( combinedCollisionRadius > 0.8*VectorOperations::vectorMagnitude(&relativeProjectedPosition) ) || EnforceNumericalIntegration ){ // Check whether to use direct integration or series expansion.
		// Use direct integration.
		std::vector< std::vector<double> > PdInverse;
		try{
			PdInverse = VectorOperations::inverseTwoByTwo( &Pd ); // At this point Pd0 is the same as Pd.
			std::vector<double> relPosMultipliedByPdInverse = VectorOperations::vectorMultiplyByMatrix(&relativeProjectedPosition, &PdInverse);
		}catch(MathWarning){
			std::string msg1 = "Impossible to inverse the covariance matrix in calculateCollisionProbability.";
			throw FatalError(msg1,__FILE__,__LINE__);
		}

		trueProbability = integrateProbabilityPDFSimpson( combinedCollisionRadius, &relativeProjectedPosition, &Pd, &PdInverse, 5000);

	}else{ // Use series expansion. No need to numerically find the worst-case covariance matrix here - in this regime the analytical guess for Pd is accurate enough.
        /* Rotate the covariance matrix and the in-plane relative position into the principal axes. */
		double rotationAngle; // FInd the angle by which to rotate.
		if( Pd.at(0).at(0)-Pd.at(1).at(1) <= DBL_MIN ){ // Both standard deviations are approximately equal, choose a fixed angle as given by K. Chan.
			rotationAngle = PI/4.0;
		}else{ // For cases where a finite covariance matrix will be specified, will never be true for diagonal matrices as used for maximum probability.
			rotationAngle = std::atan( 2.0*Pd.at(0).at(1)/(Pd.at(0).at(0) - Pd.at(1).at(1)) ); // Pd[0,1] is the cross term, Pd[0,0] and Pd[1,1] are the standard deviations squared.
		};
        
        std::vector< std::vector<double> >inPlaneRotationMatrix = std::vector< std::vector<double> >(2, std::vector<double>(2, 0.0) ); // Rotates the in-plane coordinate system into principal axes.
		inPlaneRotationMatrix.at(0).at(0)=std::cos(rotationAngle);inPlaneRotationMatrix.at(0).at(1)=-1.0*std::sin(rotationAngle);
		inPlaneRotationMatrix.at(1).at(0)=std::sin(rotationAngle);inPlaneRotationMatrix.at(1).at(1)=std::cos(rotationAngle);
        
		relativeProjectedPosition = VectorOperations::vectorMultiplyByMatrix( &relativeProjectedPosition, &inPlaneRotationMatrix); // Rotate the position into the principal axes.

        // Need to find eigenvalues here since the matrix has cross-terms. This will diagonalise it i.e. rotate into the principal axes.
		std::vector< std::vector<double> > PdDiagonalised = std::vector< std::vector<double> >(2, std::vector<double>(2, 0.0) );
		std::vector<double> PdEigenValues = EigenValues2x2(&Pd);
		PdDiagonalised.at(0).at(0) = PdEigenValues.at(0); PdDiagonalised.at(1).at(1) = PdEigenValues.at(1);
        
		trueProbability = integrateProbabilityPDFSeries(combinedCollisionRadius, &relativeProjectedPosition, &PdDiagonalised);
	};//END else (compute series-based probability).
};//END calculateCollisionProbability

void Conjunction::calculateMaximumCollisionProbability( std::vector<double>* position1, std::vector<double>* velocity1, std::vector< std::vector<double> >* RTCcovMat1, std::vector<double>* position2, std::vector<double>* velocity2, std::vector< std::vector<double> >* RTCcovMat2, bool EnforceNumericalIntegration, int noTerms ){
	/* Compute the maximum probability of collision between the two objects given their positions, velocities and position covariance matrices
    (must be in the same frame of reference and units) at the time of closest approach by scaling the covariance matrices to yield the worst-case collision probability
	according to the method given by N. Berend (1999). The relative position PDF will be either numerically integrated or K. Chan's series expansion will be used if desired
	or the combined objects' size is more than 80% of the miss distance. Direct numerical integration of the PDF can be enforced with the EnforceNumericalIntegration parameter.
	The collision probability estiamtion approach assumes linear relative trajectory, hence the provided positions and velocities
    MUST BE those of the closest approach to make this assumption viable.
    @param position1, position2 - positions of objects 1 and 2 in the form of 1x3 vectors in an inertial reference frame at the TCA.
    @param velocity1, velocity2 - velocities of those objects, also in the form of 1x3 vectors in an inertial reference frame at the TCA.
	@param RTCcovMat1, RTCcovMat2 - 3x3 position covariance matrices of the both objects at the TCA in the respective radial - in-track - cross-track frames of reference.
	@param EnforceNumericalIntegration - whether to enforce the numerical integration to be used instead of Chan's series expansion.
	@param noTerms - number of terms to use in Chan's series expansion it is's being used.
	@throws FatalError - when the probability cannot be calculated due to uninversible covariance matrices or when one of the eigenvalues of the combined covariance matrix is zero.
	@return - the found maximum collision probability can be accessed through this.maximumProbability attribute.
	*/

	// Use a reference frame centred on object1 throughout the process (SCRF), the inputs are in the World Reference Frame (WRF, e.g. TEME or ECI) or RTC of respective objects for covariance matrices.
    
    /* RELATIVE POSITION AND VELOCITY VECTORS */
    std::vector<double> relPositionWRF = std::vector<double>(3,0.0); std::vector<double> xAxisVectorSCRF = std::vector<double>(3,0.0);
	VectorOperations::vectorDifference(&relPositionWRF, position1, position2); VectorOperations::vectorDifference(&xAxisVectorSCRF, velocity1, velocity2); // Object 2 in 1-fixed frame. Account for translation from WRF to SCRF.
    VectorOperations::unitVector(&xAxisVectorSCRF); // This is the X-axis of the relative velocity-based reference frame.
    
    /* CONSTRUCTION THE RELATIVE REFERENCE FRAME. */
    // Find the remaining components of S/C-centred relative velocity-aligned reference frame.
    std::vector<double> yAxisVectorSCRF = VectorOperations::vectorMultiplyByScalar( &xAxisVectorSCRF, VectorOperations::dotProduct(&xAxisVectorSCRF, &relPositionWRF) ); // Normal to X-axis and pointing towards the intersection of the collision plane with the relative velocity line.
	VectorOperations::vectorDifference( &yAxisVectorSCRF, &relPositionWRF, &yAxisVectorSCRF ); // relPositionWRF - the extra component along the X-axis = the component of relative position normal to the X-axis.

    std::vector<double> zAxisVectorSCRF = std::vector<double>(3,0.0);
	VectorOperations::crossProduct(&zAxisVectorSCRF, &yAxisVectorSCRF, &xAxisVectorSCRF);  // Normal to Y- and X- axes.
     
    // Non-dimensionalise the remaining components (in-plane).
	VectorOperations::unitVector( &yAxisVectorSCRF );
	VectorOperations::unitVector( &zAxisVectorSCRF );
    
    // Unit vectors of the World Reference Frame.
	std::vector<double> xAxisUnitVectorWRF = std::vector<double>(3,0.0); xAxisUnitVectorWRF.at(0) = 1.0;
	std::vector<double> yAxisUnitVectorWRF = std::vector<double>(3,0.0); yAxisUnitVectorWRF.at(1) = 1.0;
	std::vector<double> zAxisUnitVectorWRF = std::vector<double>(3,0.0); zAxisUnitVectorWRF.at(2) = 1.0;
	    
    // Construct the rotation matrix from WRF to SCRF.
	std::vector< std::vector<double> > rotationMatrix = std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0) );
	rotationMatrix.at(0).at(0) = VectorOperations::dotProduct(&xAxisVectorSCRF, &xAxisUnitVectorWRF); rotationMatrix.at(0).at(1) = VectorOperations::dotProduct(&xAxisVectorSCRF, &yAxisUnitVectorWRF); rotationMatrix.at(0).at(2) = VectorOperations::dotProduct(&xAxisVectorSCRF, &zAxisUnitVectorWRF);
	rotationMatrix.at(1).at(0) = VectorOperations::dotProduct(&yAxisVectorSCRF, &xAxisUnitVectorWRF); rotationMatrix.at(1).at(1) = VectorOperations::dotProduct(&yAxisVectorSCRF, &yAxisUnitVectorWRF); rotationMatrix.at(1).at(2) = VectorOperations::dotProduct(&yAxisVectorSCRF, &zAxisUnitVectorWRF);
	rotationMatrix.at(2).at(0) = VectorOperations::dotProduct(&zAxisVectorSCRF, &xAxisUnitVectorWRF); rotationMatrix.at(2).at(1) = VectorOperations::dotProduct(&zAxisVectorSCRF, &yAxisUnitVectorWRF); rotationMatrix.at(2).at(2) = VectorOperations::dotProduct(&zAxisVectorSCRF, &zAxisUnitVectorWRF);
    
	std::vector< std::vector<double> > rotationMatrixInverse = VectorOperations::inverseThreeByThree(&rotationMatrix); // Will need these to rotate both covariance matrices to the B-plane.
	std::vector< std::vector<double> > rotationMatrixTranspose = VectorOperations::transposeMatrix(&rotationMatrix);

    /* COMPUTE RELATIVE POSITIONS AND VELOCITIES */
    // Relative position expressed in the object1-fixed reference frame (SCRF).
	std::vector<double> temp = VectorOperations::vectorMultiplyByMatrix(&relPositionWRF, &rotationMatrix);
	std::vector<double> relativeProjectedPosition = std::vector<double>(2, 0.0);
	relativeProjectedPosition.at(0) = temp.at(1); relativeProjectedPosition.at(1) = temp.at(2); // Discard the velocity (X) component, only care about the components within the collision plane.
    
    /* FIND THE COMBINED COVARIANCE MATRIX IN THE B-PLANE. */
	std::vector< std::vector<double> > TEME2RTC_primary = InertialToRTC(position1, velocity1); // Rotation matrix between the TEME and RTC reference frames at the TCA for the primary.
	std::vector< std::vector<double> > TEME2RTC_secondary = InertialToRTC(position2, velocity2); // Rotation matrix between the TEME and RTC reference frames at the TCA for the secondary.

	// RTC to TEME rotation for the primary - formulate necessary matrices and then multiply out for rotation.
	std::vector< std::vector<double> > TEME2RTCinverse_primary = VectorOperations::inverseThreeByThree(&TEME2RTC_primary); // Need this to rotate from VNC to TEME.
	std::vector< std::vector<double> > TEME2RTCtranspose_primary = VectorOperations::transposeMatrix(&TEME2RTC_primary); // Also need the transpose when rotating covaraince matrices.
	std::vector< std::vector<double> > TEME2RTCtransposeInverse_primary = VectorOperations::inverseThreeByThree(&TEME2RTCtranspose_primary);
	
	std::vector< std::vector<double> > tempMat_primary = VectorOperations::matrixMultiplyByMatrix(RTCcovMat1, &TEME2RTCtransposeInverse_primary); // The operation is RotMat^-1 CovMat RotMat.T^-1.
	std::vector< std::vector<double> > posCovarianceTEME_primary = VectorOperations::matrixMultiplyByMatrix(&TEME2RTCinverse_primary, &tempMat_primary); // Finished rotations to TEME.

	// RTC to TEME rotation for the secondary - formulate necessary matrices and then multiply out for rotation.
	std::vector< std::vector<double> > TEME2RTCinverse_secondary = VectorOperations::inverseThreeByThree(&TEME2RTC_secondary);
	std::vector< std::vector<double> > TEME2RTCtranspose_secondary = VectorOperations::transposeMatrix(&TEME2RTC_secondary);
	std::vector< std::vector<double> > TEME2RTCtransposeInverse_secondary = VectorOperations::inverseThreeByThree(&TEME2RTCtranspose_secondary);
	
	std::vector< std::vector<double> > tempMat_secondary = VectorOperations::matrixMultiplyByMatrix(RTCcovMat2, &TEME2RTCtransposeInverse_secondary);
	std::vector< std::vector<double> > posCovarianceTEME_secondary = VectorOperations::matrixMultiplyByMatrix(&TEME2RTCinverse_secondary, &tempMat_secondary); // Finished rotations to TEME.

	// Add the matrices together in TEME - assume uncorrrelated errors of both objects
	std::vector< std::vector<double> > combinedCovarianceMatrixTEME = VectorOperations::addMatrices( &posCovarianceTEME_primary, &posCovarianceTEME_secondary ); // Combined 3x3 position covariance defined in TEME at TCA.

	// Et viola, we have the combined covariance matrix rotated into the B-plane.
	std::vector< std::vector<double> > tempMat_combined = VectorOperations::matrixMultiplyByMatrix(&combinedCovarianceMatrixTEME, &rotationMatrixTranspose);
	std::vector< std::vector<double> > combinedCovarianceMatrixBplane = VectorOperations::matrixMultiplyByMatrix(&rotationMatrix, &tempMat_combined);

	std::vector< std::vector<double> > Pd0 = std::vector< std::vector<double> >(2, std::vector<double>(2, 0.0) ); // Reference covariance matrix defined in the B-plane, at the TCA and ignoring the velocity (X) components. This will never changed and will be used to find Pd by scaling with some scaling factor to yield worst-case collision probability.
	Pd0.at(0).at(0) = combinedCovarianceMatrixBplane.at(1).at(1); Pd0.at(0).at(1) = combinedCovarianceMatrixBplane.at(1).at(2); // Ignore the velocity i.e. X-component (along the relative veloctiy) of combinedCovarianceMatrixBplane.
	Pd0.at(1).at(0) = combinedCovarianceMatrixBplane.at(2).at(1); Pd0.at(1).at(1) = combinedCovarianceMatrixBplane.at(2).at(2);

	// Rotate the covariance matrix into the principal axes (diagonalise).
	std::vector<double> eigenValuesPd0 = EigenValues2x2( &Pd0 );
	Pd0.at(0).at(0) = eigenValuesPd0.at(0);Pd0.at(0).at(1) = 0.0;
	Pd0.at(1).at(1) = eigenValuesPd0.at(1); Pd0.at(1).at(0) = 0.0;

	// Rotate the relative position vector to the principal axes too.
	double rotationAngle;
	if( Pd0.at(0).at(0)-Pd0.at(1).at(1) <= DBL_MIN ){ // Both standard deviations are approximately equal, choose a fixed angle as given by K. Chan.
		rotationAngle = PI/4.0;
	}else{ // For cases where a particular covariance matrix will be specified, will never be true for diagonal matrices as used for maximum probability.
		rotationAngle = std::atan( 2.0*Pd0.at(0).at(1)/(Pd0.at(0).at(0) - Pd0.at(1).at(1)) ); // Pd0[0,1] is the cross term, Pd0[0,0] and Pd0[1,1] are the standard deviations squared.
	};
	std::vector< std::vector<double> >inPlaneRotationMatrix = std::vector< std::vector<double> >(2, std::vector<double>(2, 0.0) ); // Rotates the in-plane coordinate system into principal axes.
	inPlaneRotationMatrix.at(0).at(0)=std::cos(rotationAngle);inPlaneRotationMatrix.at(0).at(1)=-1.0*std::sin(rotationAngle);
	inPlaneRotationMatrix.at(1).at(0)=std::sin(rotationAngle);inPlaneRotationMatrix.at(1).at(1)=std::cos(rotationAngle);
	VectorOperations::vectorMultiplyByMatrix( &relativeProjectedPosition, &inPlaneRotationMatrix);

	// Scale Pd0 to worst-case to form Pd.
	double covarianceScalingFactorSquared;
	double covarianceScalingFactor; // Do the golden-ratio search on k not k^2. It converges better like this.
	std::vector< std::vector<double> > PdInverse(2, std::vector<double>(2, 0.0) );
	try{
		PdInverse = VectorOperations::inverseTwoByTwo( &Pd0 ); // At this point Pd0 is the same as Pd.
		std::vector<double> relPosMultipliedByPdInverse = VectorOperations::vectorMultiplyByMatrix(&relativeProjectedPosition, &PdInverse);
		covarianceScalingFactorSquared = VectorOperations::dotProduct(&relPosMultipliedByPdInverse, &relativeProjectedPosition)/2.0; // Scaling factor that yields the maximum collision probability according to N. Berend 1999.
		covarianceScalingFactor = std::sqrt( covarianceScalingFactorSquared );
	}catch(MathWarning){
		std::cerr<<"Impossible to inverse the covariance matrix in calculateMaximumCollisionProbability during the golden ratio search and cannot get an initial estimate of PcMAX."<<std::endl;
		std::string msg1 = "Impossible to inverse the covariance matrix in calculateMaximumCollisionProbability during the golden ratio search and cannot get an initial estimate of PcMAX.";
		throw FatalError(msg1,__FILE__,__LINE__);
	}
	// Sometimes matrix Pd0 will be badly conditioned when the covariance has been estimated using very few TLEs. This will cause the scaling factor to be undefined. Deal with it somehow.
	#ifdef _MSC_VER // MSVS has a different _isnan function from <float.h>, Linux compilers use std::isnan from <cmath>
		if( _isnan(covarianceScalingFactorSquared) ){
			covarianceScalingFactorSquared = 0.025; // Need to start the search from somewhere.
			covarianceScalingFactor = std::sqrt( covarianceScalingFactorSquared );
		}
	#else
		if( std::isnan(covarianceScalingFactorSquared) ){
			covarianceScalingFactorSquared = 0.025;
			covarianceScalingFactor = std::sqrt( covarianceScalingFactorSquared );
		}
	#endif

	std::vector< std::vector<double> > Pd = VectorOperations::matrixMultiplyByScalar( &Pd0, covarianceScalingFactorSquared ); // First guess of rescaling Pd to yield worst-case collision probability. Works when collision radius < ~80% miss distance.

	/* COMPUTE THE COLLISION PROBABILITY */
	if( ( combinedCollisionRadius > 0.8*VectorOperations::vectorMagnitude(&relativeProjectedPosition)) || EnforceNumericalIntegration ){ // Check whether to use direct integration or series expansion.
		/* Use direct integration and do a golden ratio search for the true worst-case covariance scaling factor. */
		double TOL = 0.01; // Tolerance to within which the maximum will be found.
		int MAX_ITER = 100; // Maximum number of golden ration search iterations after which we'll stop - without this some searches will never terminate and cause the algorithm to hang.
		int currentIter = 0; // Number of iterations done so far.
		double fTEMP = 0.0; // -1*Current best estimate of Pc.

		const double R=0.61803399, C=1.0-R; // The golden ratios.
		double f1,f2,x0,x1,x2,x3; // At any time keep track of four points x0,x1,x2 and x3.
									// Do a golden search for ax < bx < cx and f(bx) < f(ax) && f(bx) < f(cx).

		double ax = 0.0; // Starting points that bracket the root (normally less than covarianceScalingFactor). Sometimes it's 0.
		double bx = 0.5*covarianceScalingFactor;
		double cx = 1.5*covarianceScalingFactor;
	
		x0=ax; x3=cx;
		if( fabs(bx-cx) > fabs(cx-ax) ){ // Make x0 to x1 the smaller segment.
			x1=bx; // And fill in the new point to be tried.
			x2=bx+C*(cx-bx);
		}else{
			x2=bx;
			x1=bx-C*(bx-ax);
		};

		Pd = VectorOperations::matrixMultiplyByScalar( &Pd0, x1*x1 ); // Initial endpoint function evaluations.
		try{
			PdInverse = VectorOperations::inverseTwoByTwo( &Pd ); // Try to get inverse of Pd to estimate the Pc for x1 scaling factor, this may not work.
		}catch( MathWarning ){
			std::cerr<<"\tImpossible to inverse the covariance matrix in calculateMaximumCollisionProbability during the golden ratio search and cannot get an initial estimate of PcMAX."<<std::endl;
			std::string msg1 = "Impossible to inverse the covariance matrix in calculateMaximumCollisionProbability during the golden ratio search and cannot get an initial estimate of PcMAX.";
			throw FatalError(msg1,__FILE__,__LINE__);
		}
		f1 = -1.0*integrateProbabilityPDFSimpson( combinedCollisionRadius, &relativeProjectedPosition, &Pd, &PdInverse, 5000);

		Pd = VectorOperations::matrixMultiplyByScalar( &Pd0, x2*x2 );
		try{
			PdInverse = VectorOperations::inverseTwoByTwo( &Pd ); // Try to get inverse of Pd to estimate the Pc for x2 scaling factor, this may not work.
		}catch( MathWarning ){
			std::cerr<<"\tImpossible to inverse the covariance matrix in calculateMaximumCollisionProbability during the golden ratio search, cannot get an estimate of PcMAX."<<std::endl;
			std::string msg1 = "Impossible to inverse the covariance matrix in calculateMaximumCollisionProbability during the golden ratio searchcannot get an estimate of PcMAX.";
			throw FatalError(msg1,__FILE__,__LINE__);
		}
		f2 = -1.0*integrateProbabilityPDFSimpson( combinedCollisionRadius, &relativeProjectedPosition, &Pd, &PdInverse, 5000);

		// See if the probabilities are numbers and keep track of that, for ill-conditioned covariance matrices they'll be NaNs.
		#ifdef _MSC_VER // MSVS has a different _isnan function from <float.h>, Linux compilers use std::isnan from <cmath>
			int IS_NAN_f1, IS_NAN_f2;
			IS_NAN_f1 = _isnan(f1);
			IS_NAN_f2 = _isnan(f2);
		#else
			bool IS_NAN_f1, IS_NAN_f2;
			IS_NAN_f1 = std::isnan(f1);
			IS_NAN_f2 = std::isnan(f2);
		#endif

		while( ( fabs(x3-x0) > TOL*( fabs(x1)+fabs(x2) ) ) && !IS_NAN_f1 && !IS_NAN_f2 ){
			if(f2 < f1){ // One possible outcome of braceting the root.
				shift3(x0,x1,x2,R*x2+C*x3);
				Pd = VectorOperations::matrixMultiplyByScalar( &Pd0, x2*x2 );
				try{
					PdInverse = VectorOperations::inverseTwoByTwo( &Pd );
				}catch( MathWarning ){
					std::cerr<<"\tImpossible to inverse the covariance matrix in calculateMaximumCollisionProbability during the golden ratio search, using current best estimate of PcMAX and breaking."<<std::endl;
					break;
				};
				fTEMP = -1.0*integrateProbabilityPDFSimpson( combinedCollisionRadius, &relativeProjectedPosition, &Pd, &PdInverse, 5000);
				shift2(f1,f2,fTEMP);
				currentIter+=1;
			}else{ // The other outcome.
				shift3(x3,x2,x1,R*x1+C*x0);
				Pd = VectorOperations::matrixMultiplyByScalar( &Pd0, x1*x1 );
				try{
					PdInverse = VectorOperations::inverseTwoByTwo( &Pd );
				}catch( MathWarning ){
					std::cerr<<"\tImpossible to inverse the covariance matrix in calculateMaximumCollisionProbability during the golden ratio search, using current best estimate of PcMAX and breaking."<<std::endl;
					break;
				};
				fTEMP = -1.0*integrateProbabilityPDFSimpson( combinedCollisionRadius, &relativeProjectedPosition, &Pd, &PdInverse, 5000);
				shift2(f2,f1,fTEMP);
				currentIter+=1;
			};
			if(currentIter==MAX_ITER){
				std::cerr<<"\tTerminated iterations when looking for true maximum probability, current accuracy is "<<fabs(x3-x0)<<" and the desired "<<TOL*( fabs(x1)+fabs(x2) )<<std::endl;
				break;
			}

		}; //END while searching for the covarianceScalingFactor that yields true maximum of the collision probability.

		if( (f1 < f2) && !IS_NAN_f1 ){ // Output the best of the two current values. Sometimes they will be NaNs for ill-conditioned covariance matrices...
			covarianceScalingFactorSquared = x1*x1;
			covarianceScalingFactor = std::sqrt( covarianceScalingFactorSquared );
			maximumProbability = -1.0*f1; // Scale to the actual probability value (looking for a minimum here).
		}else if( !IS_NAN_f2 ){
			covarianceScalingFactorSquared = x2*x2;
			covarianceScalingFactor = std::sqrt( covarianceScalingFactorSquared );
			maximumProbability = -1.0*f2;
		}
	}else{ // Use series expansion. No need to numerically find the worst-case covariance matrix here - in this regime the analytical guess for Pd is accurate enough.
		/* 11 mar 2015 - do the golden ratio search anyway to be able to prove that this works in the thesis. Otherwise run into risk that some super close conjunctions have higher PcMAX than the medioka ones. */
        /* Rotate the covariance matrix into the principal axes. */
		double rotationAngle;
		if( Pd.at(0).at(0)-Pd.at(1).at(1) <= DBL_MIN ){ // Both standard deviations are approximately equal, choose a fixed angle as given by K. Chan.
			rotationAngle = PI/4.0;
		}else{ // For cases where a specific covariance matrix will be specified, will never be true for diagonal matrices as used for maximum probability.
			rotationAngle = std::atan( 2.0*Pd.at(0).at(1)/(Pd.at(0).at(0) - Pd.at(1).at(1)) ); // Pd[0,1] is the cross term, Pd[0,0] and Pd[1,1] are the standard deviations squared.
		};
        
        std::vector< std::vector<double> >inPlaneRotationMatrix = std::vector< std::vector<double> >(2, std::vector<double>(2, 0.0) ); // Rotates the in-plane coordinate system into principal axes.
		inPlaneRotationMatrix.at(0).at(0)=std::cos(rotationAngle);inPlaneRotationMatrix.at(0).at(1)=-1.0*std::sin(rotationAngle);
		inPlaneRotationMatrix.at(1).at(0)=std::sin(rotationAngle);inPlaneRotationMatrix.at(1).at(1)=std::cos(rotationAngle);
        
		VectorOperations::vectorMultiplyByMatrix( &relativeProjectedPosition, &inPlaneRotationMatrix); // Rotate to principal axes.
        // Would need to find eigenvalues of Pd here if the matrix had a cross-term. But this is never the case for spherical probability.
        
		/* Do a golden ratio search for the true worst-case covariance scaling factor. */
		double TOL = 0.01; // Tolerance to within which the maximum will be found.
		int MAX_ITER = 100; // Maximum number of iterations allowed.
		int currentIter = 0; // Number of iterations done already.
		double fTEMP = 0.0; // -1*Current best estimate of Pc.

		std::vector< std::vector<double> > PdInverse = VectorOperations::inverseTwoByTwo( &Pd );
		
		const double R=0.61803399, C=1.0-R; // The golden ratios.
		double f1,f2,x0,x1,x2,x3; // At any time keep track of four points x0,x1,x2 and x3.
								  // Do a golden search for ax < bx < cx and f(bx) < f(ax) && f(bx) < f(cx).

		double ax = 0.0; // Starting points that bracket the root (normally less than covarianceScalingFactor). Sometimes it's 0.
		double bx = 0.5*covarianceScalingFactor;
		double cx = 1.5*covarianceScalingFactor;
		
		x0=ax; x3=cx;
		if( fabs(bx-cx) > fabs(cx-ax) ){ // Make x0 to x1 the smaller segment.
			x1=bx; // And fill in the new point to be tried.
			x2=bx+C*(cx-bx);
		}else{
			x2=bx;
			x1=bx-C*(bx-ax);
		};

		Pd = VectorOperations::matrixMultiplyByScalar( &Pd0, x1*x1 ); // Initial endpoint function evaluations.
		f1 = -1.0*integrateProbabilityPDFSeries(combinedCollisionRadius, &relativeProjectedPosition, &Pd);

		Pd = VectorOperations::matrixMultiplyByScalar( &Pd0, x2*x2 );
		f2 = -1.0*integrateProbabilityPDFSeries(combinedCollisionRadius, &relativeProjectedPosition, &Pd);

		// See if the probabilities are numbers and keep track of that, for ill-conditioned covariance matrices they'll be NaNs.
		#ifdef _MSC_VER // MSVS has a different _isnan function from <float.h>, Linux compilers use std::isnan from <cmath>
			int IS_NAN_f1, IS_NAN_f2;
			IS_NAN_f1 = _isnan(f1);
			IS_NAN_f2 = _isnan(f2);
		#else
			bool IS_NAN_f1, IS_NAN_f2;
			IS_NAN_f1 = std::isnan(f1);
			IS_NAN_f2 = std::isnan(f2);
		#endif

		while( (fabs(x3-x0) > TOL*( fabs(x1)+fabs(x2) )) && !IS_NAN_f1 && !IS_NAN_f2 ){
			if(f2 < f1){ // One possible outcome of the golden ratio root braceting.
				shift3(x0,x1,x2,R*x2+C*x3);
				Pd = VectorOperations::matrixMultiplyByScalar( &Pd0, x2*x2 );
				fTEMP = -1.0*integrateProbabilityPDFSeries(combinedCollisionRadius, &relativeProjectedPosition, &Pd);
				shift2(f1,f2,fTEMP);
				currentIter += 1;
			}else{ // The other outcome.
				shift3(x3,x2,x1,R*x1+C*x0);
				Pd = VectorOperations::matrixMultiplyByScalar( &Pd0, x1*x1 );
				fTEMP = -1.0*integrateProbabilityPDFSeries(combinedCollisionRadius, &relativeProjectedPosition, &Pd);
				shift2(f2,f1,fTEMP);
				currentIter += 1;
			};

			if(currentIter==MAX_ITER){
				std::cerr<<"Terminated iterations when looking for true maximum probability, current accuracy is "<<fabs(x3-x0)<<" and the desired "<<TOL*( fabs(x1)+fabs(x2) )<<". Using current best Pc_MAX estimate."<<std::endl;
				break;
			}
		}; //END while searching for the worst-case covariance scaling factor

		if( (f1 < f2) && !IS_NAN_f1){ // Output the best of the two current values.
			covarianceScalingFactor = x1;
			maximumProbability = -1.0*f1; // Scale to the actual probability value (looking for a minimum here).
		}else if( !IS_NAN_f2 ){
			covarianceScalingFactor = x2;
			maximumProbability = -1.0*f2;
		};
	};//END else (compute series-based probability).

/* 29 Aug 2014 - conjunctions with PcMAX > 1.0 will now be filtered out in post-processing as they are caused by erroneous state vectors.
	if( maximumProbability>1.0){ // Throw an error if this happens.
		std::stringstream ss;
		ss<<"Pd:"<<Pd.at(0).at(0)<<","<<Pd.at(0).at(1)<<";"<<Pd.at(1).at(0)<<","<<Pd.at(1).at(1)<<" Collision radius:"<<combinedCollisionRadius<<" R:"<<relativeProjectedPosition.at(0)<<","<<relativeProjectedPosition.at(1)<<" Maximum Pc greater than 1.0.";
		throw FatalError(ss.str(),__FILE__,__LINE__);
	}*/

};//END calculateMaximumCollisionProbability

Conjunction::~Conjunction(void){
	 /* Deconstructor, does nothing special.*/
};

/* HELPER FUNCTIONS. */
std::ostream& operator<<(std::ostream& os, const Conjunction& conj){
	/* Used to pretty-print the Conjunction class. */
    os << conj.primaryNORAD_ID << ',' << conj.secondaryNORAD_ID << ',' << conj.readableTCA << ',' << conj.missDistance << ',' << conj.maximumProbability<< ',' << conj.trueProbability << ',' << conj.relativeVelocity <<','<< conj.combinedCollisionRadius;
    return os;
}

void shift2(double &a, double &b, const double c){
	/* Replace value of a with value of b, and the value of b with value of c.
	@reference - Numerical Recipes in C++, Second Edition.
	*/
	a=b;
	b=c;
};

void shift3(double &a, double &b, double &c, const double d){
	/* Replace the value of a with the value of b, the value of b with the value of c,
	and the value of c with the value of d.
	@reference - Numerical Recipes in C++, Second Edition.
	*/
	a=b;
	b=c;
	c=d;
};

double integrateProbabilityPDFTrapeze( double combinedCollisionRadius, std::vector<double>* relativeProjectedPosition, std::vector< std::vector<double> >* Pd, std::vector< std::vector<double> >* PdInverse, int noIntervals){
	/* Given the size of two objects, their position in the collision plane as well as the in-plane covariance matrix, which defines their
	relative position probability density function, compute the probability that the object will collide. This is done by numerically 
	integrating the PDF over a circle with the combinedCollisionRaius using a Trapezium method.

	@param combinedCollisionRadius - radii of the two objects, in the same units as other arguments, assumed km.
	@param relativeProjectedPosition - relative position of the objects in the collision plane.
	@param Pd - 2x2 position covariance matrix (full 3x3 matrix projected onto the collision plane and ignoring the velocity component).
	@param PdInverse - inverse of the covaraince matrix.
	@param noIntervals - number of intervals to be used to numerically integrate the PDF using the trapezium method. Around 50 000 is needed to get
					~10% accurate results for very low probabilities, around 1000 is enough for accurate maximum probabilities.
	@return - the collision probability of the two objects n this situation.
	*/
	/* Integration variables. */
	double hZ = 2.0*combinedCollisionRadius/((double)noIntervals + 2.0); // Interpolation steps. Need the extra 2 because we don't want to be looking at currentZ=combinedCollisionRadius as then currentY is undefined. In order to deal with this ignore two thin Z-strips close to the Z limits.
	double hY; // This depends on the other coordinate's value - recompute it at every step in the outer variable.
	double probability = 0.0;

	double yIntegral; // Temporary integral along Y (recomputed at each time step).
	double currentZ = -1.0*combinedCollisionRadius+hZ; // Start at the bottom of the circle and integrate up in Z and left to right in Y.
	double currentY, currentYlimit, expTerm; // Temporary variables.

	std::vector<double> discrepancyVector = std::vector<double>(2,0.0); // Instantenous [z y] (used in integration) minus the mean values.
	std::vector<double> discrepacyVectorMultipliedByPdInverse = std::vector<double>(2,0.0); // A temporary entity.

	/* 
	 * Use the trapezium rule to integrate the PDF over the circle with radius of collision radius.
	 * Only integrate over a quarter of the circle and then multiply by 4.
	 */
	for(int k=0; k<noIntervals; k++){ // Take steps along the Z-direction. Integrate along Y at each of them.
		hY = std::sqrt( std::pow(combinedCollisionRadius,2.0) - std::pow(currentZ,2.0) )/(double)noIntervals; // And the length of the step in Y at the current Z.
		currentYlimit = std::sqrt( std::pow(combinedCollisionRadius,2.0) - std::pow(currentZ,2.0) ); // Maximum value of Y at the current Z.
		currentY = -1.0*std::sqrt( std::pow(combinedCollisionRadius,2.0) - std::pow(currentZ,2.0) ); // Start at negative Y limit and go towards the positive limit.
		yIntegral = 0.0; // Reset this at every Y step.

		/* Integrate in Y along the current Z. */
		discrepancyVector.at(0) = currentZ - relativeProjectedPosition->at(0);
		discrepancyVector.at(1) = currentY - relativeProjectedPosition->at(1);

		discrepacyVectorMultipliedByPdInverse = VectorOperations::vectorMultiplyByMatrix(&discrepancyVector, PdInverse);
		expTerm = std::exp( -0.5 * VectorOperations::dotProduct( &discrepacyVectorMultipliedByPdInverse, &discrepancyVector ) ); // This is a rather lengthy expression so divide it into parts.
		yIntegral += 0.5*expTerm/( 2.0*PI*std::sqrt( VectorOperations::twoByTwoDeterminant(Pd) ) ); // First and last points where the function was evaluated, treat them slightly differently than the rest i.e. add 0.5.

		currentY += hY;
		while( currentY < currentYlimit ){ // Only go through the interior points of the trapezium rule.
			discrepancyVector.at(0) = currentZ - relativeProjectedPosition->at(0);
			discrepancyVector.at(1) = currentY - relativeProjectedPosition->at(1);

		    	discrepacyVectorMultipliedByPdInverse = VectorOperations::vectorMultiplyByMatrix(&discrepancyVector, PdInverse);
			expTerm = std::exp( -0.5 * VectorOperations::dotProduct( &discrepacyVectorMultipliedByPdInverse, &discrepancyVector ) ); // This is a rather lengthy expression so divide it into parts.
			yIntegral += expTerm/( 2.0*PI*std::sqrt( VectorOperations::twoByTwoDeterminant(Pd) ) );
			currentY += hY;
		};
		discrepancyVector.at(0) = currentZ - relativeProjectedPosition->at(0);
		discrepancyVector.at(1) = currentY - relativeProjectedPosition->at(1);

		discrepacyVectorMultipliedByPdInverse = VectorOperations::vectorMultiplyByMatrix(&discrepancyVector, PdInverse);
		expTerm = std::exp( -0.5 * VectorOperations::dotProduct( &discrepacyVectorMultipliedByPdInverse, &discrepancyVector ) ); // This is a rather lengthy expression so divide it into parts.
		yIntegral += 0.5*expTerm/( 2.0*PI*std::sqrt( VectorOperations::twoByTwoDeterminant(Pd) ) ); // Now Y is at the last point.
 
		yIntegral *= hY; // Limit the number of flops performed - only multiply once (trapezium rule).

		/* Add the contribution to the total integral for the current Z. */
		if(k==0 || k==noIntervals-1){
			probability += yIntegral/2.0; // Special treatment for the endpoints along the X-direction.
		}else{
			probability += yIntegral; // Interior points.
		};

		/* Take a step in Z, recompute/reset all the variables as necessary. */
		currentZ += hZ; // Take a step in Z.
	};

	probability *= hZ; // Finish the trapezium integration in Z.

	return probability;
};

double integrateProbabilityPDFSeries(double combinedCollisionRadius, std::vector<double>* relativeProjectedPosition, std::vector< std::vector<double> >* Pd, int noTerms){
	/* Given the size of two objects, their position in the collision plane as well as the diagonalised in-plane covariance matrix, which defines their
	relative position probability density function, compute the probability that the object will collide. This is done by expanding the PDF of their
	relative position into an infinite series as given by K. Chan 2009 and truncating after a given number of terms.

	@param combinedCollisionRadius - radii of the two objects, in the same units as other arguments, assumed km.
	@param relativeProjectedPosition - relative position of the objects in the collision plane.
	@param Pd - diagonalised 2x2 position covariance matrix (full 3x3 matrix projected onto the collision plane and ignoring the velocity component).
	@param noTerms - number of terms to be used to expand the PDF.
	@return - the collision probability of the two objects n this situation.
	*/
	double probability=0.0;
	
	// Pd is already diagonalised and scaled to worst-case, and the relativeProjectedPosition rotated to the principal axes.
	double u = combinedCollisionRadius*combinedCollisionRadius/( std::sqrt( Pd->at(0).at(0) ) * std::sqrt( Pd->at(1).at(1) ) ); // Dimensionless variable.
	double v = VectorOperations::dotProduct(relativeProjectedPosition, relativeProjectedPosition) * ( 1.0 + ( Pd->at(1).at(1)/Pd->at(0).at(0) - 1.0 )*( std::pow(relativeProjectedPosition->at(0),2.0) / (std::pow(relativeProjectedPosition->at(0),2.0) + std::pow(relativeProjectedPosition->at(1),2.0)) ) )/Pd->at(1).at(1); // Another dimensionless variable.
    
	/* Compute the collision probability with series expansion. */
	if(noTerms==1){
	    probability = std::exp(-v/2.0)*( 1.0-std::exp(-u/2.0) );
	}else if(noTerms==2){
	    probability = std::exp(-v/2.0)*( 1.0-std::exp(-u/2.0) + v/2.0*(1.0-(1.0+u/2.0)*std::exp(-u/2.0)) );
	}else{
		std::vector<double> mTerms = std::vector<double>();
		std::vector<double> kTerms = std::vector<double>();
		for(int m=0; m<noTerms; m++){
			kTerms.clear(); // Remove all the entries from the vector to calculate new kTerms.
			for(int k=0; k<m+1; k++){
				double factorialK=1.0;
				for(double b = 1.0; b <= k; b++) {
					factorialK *= b;
				};
//TODO: add a check here to look for cases when any of the kTerms turns to NaN as then the entire probability will. If it happens it's tiny - change it to DBL_MIN.
			    	kTerms.push_back( std::pow(u,k) / (std::pow(2.0,k) * factorialK) );
			};

			double factorialM=1.0;
			for(double b = 1.0; b <= m; b++) {
				factorialM *= b;
			};

			double sumKterms=0.0;
			for(std::vector<int>::size_type  i=0; i<kTerms.size(); i++){
				sumKterms += kTerms.at(i);
			};

			double expTermTemp = 1.0-std::exp(-1.0*u/2.0)*sumKterms;
			mTerms.push_back( std::pow(v, m)/( std::pow(2.0,m) * factorialM )*( expTermTemp ) );
		};

	double sumMterms=0.0;
	for(std::vector<int>::size_type  i=0; i<mTerms.size(); i++){
		sumMterms += mTerms.at(i);
	};

	probability = std::exp(-1.0*v/2.0)*sumMterms;
	};//END else (compute the probability for any number of terms noTerms).
	
	return probability;
}
		
double integrateProbabilityPDFSimpson( double combinedCollisionRadius, std::vector<double>* relativeProjectedPosition, std::vector< std::vector<double> >* Pd, std::vector< std::vector<double> >* PdInverse, int noIntervals){
	/* Given the size of two objects, their position in the collision plane as well as the in-plane covariance matrix, which defines their
	relative position probability density function, compute the probability that the object will collide. This is done by numerically 
	integrating the PDF over a circle with the combinedCollisionRaius using the Simpson's rule.

	@param combinedCollisionRadius - radii of the two objects, in the same units as other arguments, assumed km.
	@param relativeProjectedPosition - relative position of the objects in the collision plane.
	@param Pd - 2x2 position covariance matrix (full 3x3 matrix projected onto the collision plane and ignoring the velocity component).
	@param PdInverse - inverse of the covaraince matrix.
	@param noIntervals - number of intervals to be used to numerically integrate the PDF using the Simpson's rule (number of domain subdivisions).
	@return - the collision probability of the two objects n this situation.
	*/
	/* Integration variables. */
	double hZ = 2.0*combinedCollisionRadius/((double)noIntervals + 2.0); // Interpolation steps. Need the extra 2 because we don't want to be looking at currentZ=combinedCollisionRadius as then currentY is undefined. In order to deal with this ignore two thin Z-strips close to the Z limits.
	double hY; // This depends on the other coordinate's value - recompute it at every step in the outer variable.
	double probability=0.0;

	double yIntegral = 0.0; // Temporary integral along Y (recomputed at each time step).
	double currentZ = -1.0*combinedCollisionRadius+hZ; // Start at the bottom of the circle and integrate up in Z and left to right in Y.
	double currentY = -1.0*std::sqrt( std::pow(combinedCollisionRadius,2.0) - std::pow(currentZ,2.0) ); // Start at negative Y limit and go towards the positive limit.
	double expTerm; // Temporary variable.

	std::vector<double> discrepancyVector = std::vector<double>(2,0.0); // Instantenous [z y] (used in integration) minus the mean values.
	std::vector<double> discrepacyVectorMultipliedByPdInverse = std::vector<double>(2,0.0); // A temporary entity.

	/* 
	 * Use the Simpson's rule to integrate the PDF over the circle with radius of collision radius.
	 * Only integrate over a quarter of the circle and then multiply by 4.
	 */

	/* Get the contributions from the first step along Z with weight 1. */
	hY = 2.0*std::sqrt( std::pow(combinedCollisionRadius,2.0) - std::pow(currentZ,2.0) )/(double)noIntervals;
	discrepancyVector.at(0) = currentZ - relativeProjectedPosition->at(0);
	discrepancyVector.at(1) = currentY - relativeProjectedPosition->at(1);
	discrepacyVectorMultipliedByPdInverse = VectorOperations::vectorMultiplyByMatrix(&discrepancyVector, PdInverse);
	expTerm = std::exp( -0.5 * VectorOperations::dotProduct( &discrepacyVectorMultipliedByPdInverse, &discrepancyVector ) ); // This is a rather lengthy expression so divide it into parts.
	yIntegral += expTerm/( 2.0*PI*std::sqrt( VectorOperations::twoByTwoDeterminant(Pd) ) ); // First point has weight 1.

	for(int j=1; j<noIntervals-1; j++){ // Internal points are handled differently.
		currentY += hY;
		discrepancyVector.at(0) = currentZ - relativeProjectedPosition->at(0);
		discrepancyVector.at(1) = currentY - relativeProjectedPosition->at(1);
		discrepacyVectorMultipliedByPdInverse = VectorOperations::vectorMultiplyByMatrix(&discrepancyVector, PdInverse);
		expTerm = std::exp( -0.5 * VectorOperations::dotProduct( &discrepacyVectorMultipliedByPdInverse, &discrepancyVector ) );

		if(j % 2){ // j is odd.
			yIntegral += 4*expTerm/( 2.0*PI*std::sqrt( VectorOperations::twoByTwoDeterminant(Pd) ) );; // Add odd points with weitht 4...
		}else{
			yIntegral += 2*expTerm/( 2.0*PI*std::sqrt( VectorOperations::twoByTwoDeterminant(Pd) ) );; // ...even with weight 2.
		}
	}
	currentY += hY;
	discrepancyVector.at(0) = currentZ - relativeProjectedPosition->at(0);
	discrepancyVector.at(1) = currentY - relativeProjectedPosition->at(1);
	discrepacyVectorMultipliedByPdInverse = VectorOperations::vectorMultiplyByMatrix(&discrepancyVector, PdInverse);
	expTerm = std::exp( -0.5 * VectorOperations::dotProduct( &discrepacyVectorMultipliedByPdInverse, &discrepancyVector ) );
	yIntegral += expTerm; // Last point has weight 1.
	yIntegral = yIntegral * ( hY / 3.0 ); // Multiply by length of each step/3 to get the integral along this theta strip.

	probability += yIntegral; // 1.0 weight here as it's the first Z step.
	currentZ += hZ; // Go to the next Z.

	/* Get contributions for all the internal Z steps. */
	for(int i=1; i<noIntervals-1; i++){
		hY = 2.0*std::sqrt( std::pow(combinedCollisionRadius,2.0) - std::pow(currentZ,2.0) )/(double)noIntervals; // Need to re-compute this as it depends on the circumference at current Z.
		currentY = -1.0*std::sqrt( std::pow(combinedCollisionRadius,2.0) - std::pow(currentZ,2.0) ); yIntegral = 0.0; // Need to reset these every time.
		discrepancyVector.at(0) = currentZ - relativeProjectedPosition->at(0);
		discrepancyVector.at(1) = currentY - relativeProjectedPosition->at(1);
		discrepacyVectorMultipliedByPdInverse = VectorOperations::vectorMultiplyByMatrix(&discrepancyVector, PdInverse);
		expTerm = std::exp( -0.5 * VectorOperations::dotProduct( &discrepacyVectorMultipliedByPdInverse, &discrepancyVector ) ); // This is a rather lengthy expression so divide it into parts.
		yIntegral += expTerm/( 2.0*PI*std::sqrt( VectorOperations::twoByTwoDeterminant(Pd) ) ); // First point has weight 1.

		for(int j=1; j<noIntervals-1; j++){ // Internal points are handled differently.
			currentY += hY;
			discrepancyVector.at(0) = currentZ - relativeProjectedPosition->at(0);
			discrepancyVector.at(1) = currentY - relativeProjectedPosition->at(1);
			discrepacyVectorMultipliedByPdInverse = VectorOperations::vectorMultiplyByMatrix(&discrepancyVector, PdInverse);
			expTerm = std::exp( -0.5 * VectorOperations::dotProduct( &discrepacyVectorMultipliedByPdInverse, &discrepancyVector ) );

			if(j % 2){ // j is odd.
				yIntegral += 4*expTerm/( 2.0*PI*std::sqrt( VectorOperations::twoByTwoDeterminant(Pd) ) );; // Add odd points with weitht 4...
			}else{
				yIntegral += 2*expTerm/( 2.0*PI*std::sqrt( VectorOperations::twoByTwoDeterminant(Pd) ) );; // ...even with weight 2.
			}
		}
		currentY +=hY;
		discrepancyVector.at(0) = currentZ - relativeProjectedPosition->at(0);
		discrepancyVector.at(1) = currentY - relativeProjectedPosition->at(1);
		discrepacyVectorMultipliedByPdInverse = VectorOperations::vectorMultiplyByMatrix(&discrepancyVector, PdInverse);
		expTerm = std::exp( -0.5 * VectorOperations::dotProduct( &discrepacyVectorMultipliedByPdInverse, &discrepancyVector ) );
		yIntegral += expTerm;
		yIntegral = yIntegral * ( hY / 3.0 );

		// Find the contribution of this integral along theta to the combined 2D integral along R.
		if(i % 2){ // j is odd.
			probability += 4.*yIntegral; // Add odd points with weight 4...
		}else{
			probability += 2.*yIntegral; // ...even with weight 2.
		}
		currentZ += hZ; // Porceed to evaluate next Z location.
	}

	/* Get the contributions from the last step along Z with weight 1. */
	currentZ += hZ; // Go to the last Z.
	hY = 2.0*std::sqrt( std::pow(combinedCollisionRadius,2.0) - std::pow(currentZ,2.0) )/(double)noIntervals;
	currentY = -1.0*std::sqrt( std::pow(combinedCollisionRadius,2.0) - std::pow(currentZ,2.0) ); yIntegral = 0.0;
	discrepancyVector.at(0) = currentZ - relativeProjectedPosition->at(0);
	discrepancyVector.at(1) = currentY - relativeProjectedPosition->at(1);
	discrepacyVectorMultipliedByPdInverse = VectorOperations::vectorMultiplyByMatrix(&discrepancyVector, PdInverse);
	expTerm = std::exp( -0.5 * VectorOperations::dotProduct( &discrepacyVectorMultipliedByPdInverse, &discrepancyVector ) ); // This is a rather lengthy expression so divide it into parts.
	yIntegral += expTerm/( 2.0*PI*std::sqrt( VectorOperations::twoByTwoDeterminant(Pd) ) ); // First point has weight 1.

	for(int j=1; j<noIntervals-1; j++){ // Internal points are handled differently.
		currentY += hY;
		discrepancyVector.at(0) = currentZ - relativeProjectedPosition->at(0);
		discrepancyVector.at(1) = currentY - relativeProjectedPosition->at(1);
		discrepacyVectorMultipliedByPdInverse = VectorOperations::vectorMultiplyByMatrix(&discrepancyVector, PdInverse);
		expTerm = std::exp( -0.5 * VectorOperations::dotProduct( &discrepacyVectorMultipliedByPdInverse, &discrepancyVector ) );

		if(j % 2){ // j is odd.
			yIntegral += 4*expTerm/( 2.0*PI*std::sqrt( VectorOperations::twoByTwoDeterminant(Pd) ) );; // Add odd points with weitht 4...
		}else{
			yIntegral += 2*expTerm/( 2.0*PI*std::sqrt( VectorOperations::twoByTwoDeterminant(Pd) ) );; // ...even with weight 2.
		}
	}
	currentY +=hY;
	discrepancyVector.at(0) = currentZ - relativeProjectedPosition->at(0);
	discrepancyVector.at(1) = currentY - relativeProjectedPosition->at(1);
	discrepacyVectorMultipliedByPdInverse = VectorOperations::vectorMultiplyByMatrix(&discrepancyVector, PdInverse);
	expTerm = std::exp( -0.5 * VectorOperations::dotProduct( &discrepacyVectorMultipliedByPdInverse, &discrepancyVector ) );
	yIntegral += expTerm;
	yIntegral = yIntegral * ( hY / 3.0 );

	probability += yIntegral; // 1.0 weight here too.
	
	/* Correct for the step length (h/3) along Z. */
	probability = probability * ( hZ / 3.0 );

	return probability;
};

std::vector< std::vector<double> > InertialToRTC(std::vector<double>* positionPtr, std::vector<double>* velocityPtr){
    /* Compute a rotation matrix from an inertial frame of reference (TEME, ECI or similar) to 
    Radial - Transverse - In-track (also called Radial - Along-track - Cross-ctrack) frame.
    @param position, velcity - dimensional position and velocity in the inertial frame. Must have the same units, assumed km and km/sec.
    @return - 3x3 rotation matrix from inertial to RTC frames - DOES NOT INCLUDE TRANSLATION BETWEEN ORIGINS.
	*/
	/* Base vectors of the RTC frame expressed in the inertial frame. */
    std::vector<double> radialUnit_inertial = VectorOperations::vectorMultiplyByScalar( positionPtr, 1.0/VectorOperations::vectorMagnitude(positionPtr) ); // Unit radial vector expressed in inertial frame.
	
    std::vector<double> crossUnit_inertial = std::vector<double>(3, 0.0); // Cross-track unit vector expressed in the inertial reference frame.
	VectorOperations::crossProduct( &crossUnit_inertial, positionPtr, velocityPtr);
	VectorOperations::unitVector( &crossUnit_inertial ); // Make unit.
    
	std::vector<double> transverseUnit_inertial = std::vector<double>(3, 0.0); // Transverse unit vector expressed in the inertial reference frame.
	VectorOperations::crossProduct( &transverseUnit_inertial, &crossUnit_inertial, &radialUnit_inertial );
	VectorOperations::unitVector( &transverseUnit_inertial ); // Make unit.

	/* Base vectors of the inertial frame expressed in the inertial frame. */
    std::vector<double> xAxisUnitVectorInertial = std::vector<double>(3, 0.0); xAxisUnitVectorInertial.at(0)=1.0;
    std::vector<double> yAxisUnitVectorInertial = std::vector<double>(3, 0.0); yAxisUnitVectorInertial.at(1)=1.0;
    std::vector<double> zAxisUnitVectorInertial = std::vector<double>(3, 0.0); zAxisUnitVectorInertial.at(2)=1.0;

	/* Formulate the rotation matrix from the inertial to RTC frame. */
    std::vector< std::vector<double> > rotationMatrix = std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0) ); // Full 3x3 rotation matrix from inertial to RTC coordinate frames.
	rotationMatrix.at(0).at(0) = VectorOperations::dotProduct( &radialUnit_inertial, &xAxisUnitVectorInertial );
	rotationMatrix.at(0).at(1) = VectorOperations::dotProduct( &radialUnit_inertial, &yAxisUnitVectorInertial );
	rotationMatrix.at(0).at(2) = VectorOperations::dotProduct( &radialUnit_inertial, &zAxisUnitVectorInertial );

	rotationMatrix.at(1).at(0) = VectorOperations::dotProduct( &transverseUnit_inertial, &xAxisUnitVectorInertial );
	rotationMatrix.at(1).at(1) = VectorOperations::dotProduct( &transverseUnit_inertial, &yAxisUnitVectorInertial );
	rotationMatrix.at(1).at(2) = VectorOperations::dotProduct( &transverseUnit_inertial, &zAxisUnitVectorInertial );

	rotationMatrix.at(2).at(0) = VectorOperations::dotProduct( &crossUnit_inertial, &xAxisUnitVectorInertial );
	rotationMatrix.at(2).at(1) = VectorOperations::dotProduct( &crossUnit_inertial, &yAxisUnitVectorInertial );
	rotationMatrix.at(2).at(2) = VectorOperations::dotProduct( &crossUnit_inertial, &zAxisUnitVectorInertial );

	return rotationMatrix;
};

std::vector<double> EigenValues2x2( std::vector<std::vector<double> >* matPtr){
	/* Find eigenvalues of a 2x2 matrix for which an easy analytical solution is available. They will be returned in descending order.
	@param matPtr - pointer to a 2x2 matrix for which eigenvalues will be found.
	@return - vector of length 2 containing the two eigenvalues.
	*/
	std::vector<double> result = std::vector<double>(2, 0.0);

	double eig1 = 0.5*(-1.0*std::sqrt( matPtr->at(0).at(0)*matPtr->at(0).at(0) - 2.0*matPtr->at(0).at(0)*matPtr->at(1).at(1) + 4.0*matPtr->at(0).at(1)*matPtr->at(1).at(0) + matPtr->at(1).at(1)*matPtr->at(1).at(1) ) + matPtr->at(0).at(0) + matPtr->at(1).at(1));
	double eig2 = 0.5*(     std::sqrt( matPtr->at(0).at(0)*matPtr->at(0).at(0) - 2.0*matPtr->at(0).at(0)*matPtr->at(1).at(1) + 4.0*matPtr->at(0).at(1)*matPtr->at(1).at(0) + matPtr->at(1).at(1)*matPtr->at(1).at(1) ) + matPtr->at(0).at(0) + matPtr->at(1).at(1));
	
	if(eig1>eig2){
	result.at(0) = eig1;
	result.at(1) = eig2;
	} else {
	result.at(0) = eig2;
	result.at(1) = eig1;	
	}
	return result;
};

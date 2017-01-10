#include "ExternalSpaceObject.h"

ExternalSpaceObject::ExternalSpaceObject(void){
	/* Default constructor, does nothing special. */
};

ExternalSpaceObject::ExternalSpaceObject(std::string SSC, std::string ephemerisFileName, int COV, int wgsEllipsoid, double objectRadius){
	/* Create an object that represents something orbiting the Earth - an interface to handling external ephemeris files and performing
	conjunction assessment for them.
	@param SSC - NORAD ID catalogue number of the object, may be an empty string. It's used not to find conjunctions between this object and itself if it's recorded in the catalogue.
	@param ephemerisFileName - Name of a STK v 10 ephemeris file that contains the ephemeris of this object that will be read. Should contain Julian day epoch, Cartesian positions and velocities in meteres and metres per second in TEME reference frame.
	@param COV - type of covariance to be used for the objects. If 1 an uncertainty sphere 1:1:1 will be used and only maximum collision probability will be computed.
	@param wgsEllispoid - specifies which WGS Earth ellipsoid to use when computing altitudes etc., possible options are 72, 84 and 721 (old WGS72).
	@param objectRadius - radius of the object in metres, used to compute the probability of collisions.
	*/

	/* Define the desired EGS ellipsoid. */
	earthEllipsoid = wgsEllipsoid;

	if(wgsEllipsoid==84){
		gravityModel = wgs84;
	}else if(wgsEllipsoid==72){
		gravityModel = wgs84;
	}else if(wgsEllipsoid==721){
		gravityModel = wgs72;
	};
	getgravconst( gravityModel, tumin, mu, Re, xke, j2, j3, j4, j3oj2 ); // Get values for these variables.

	/* Load the ephemeris table. */
	NORAD_ID = SSC;
	EPHEMERIS_FILE_NAME = ephemerisFileName;
	LoadEphemerisTable( ephemerisFileName );

	/* Get other values of interest, mianly orbital elements. */
	EpochJDAY = EphemerisEpochs.at(0); // Epoch of the first ephemeris point.
	r0 = EphemerisPositions.at(0);
	v0 = EphemerisVelocities.at(0);

	CURRENT_RADIUS = VectorOperations::vectorMagnitude( &r0 ); // Don't compute osculating orbital elements and determine the periapsis and apoapsis radii hence. Use the current altitude band instead.

	hardBodyRadius = objectRadius/1000.; // Default value in km.

	//UNDONE: figure out how to handle real covariance for externally-generated ephemeris tables.
	if(COV==1){
		FullCovarianceMatrixRTC = std::vector< std::vector<double> >(6, std::vector<double>(6, 0.0)); // Initialise the covariance matrices.
		PositionCovarianceMatrixRTC = std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0));

		// Initialise the diagonal entries of the covariance matrices to simulate the desired aspect ratio.
		FullCovarianceMatrixRTC.at(0).at(0)=1.0; FullCovarianceMatrixRTC.at(1).at(1)=1.0; FullCovarianceMatrixRTC.at(2).at(2)=1.0; FullCovarianceMatrixRTC.at(3).at(3)=1.0; FullCovarianceMatrixRTC.at(4).at(4)=1.0; FullCovarianceMatrixRTC.at(5).at(5)=1.0;
		FullCovarianceMatrixRTC.at(0).at(0)=1.0; FullCovarianceMatrixRTC.at(1).at(1)=1.0; FullCovarianceMatrixRTC.at(2).at(2)=1.0;
	}

	/* INITIALISE INTERPOLATION PARAMETERS. */
	NO_INTERPOLATION_POINTS = 2;
	interpolationArguments = std::vector<double>( NO_INTERPOLATION_POINTS, 0.0 );
	interpolationArguments.at(0) = 0.0; // Currently interpolate from two points only.
	interpolationArguments.at(1) = 1.0;

	rhsX = std::vector<double>(2*NO_INTERPOLATION_POINTS,0.0); // Right hand sides of the sets of equations to be solved to find the interpolation coefficients..
	rhsY = std::vector<double>(2*NO_INTERPOLATION_POINTS,0.0);
	rhsZ = std::vector<double>(2*NO_INTERPOLATION_POINTS,0.0);

	/* Build the matrix of coefficients to solve the system of equations (also include rows corresponding to derivatives in the lower half of the matrix):
	[1, argumnet, argument^2, ..., argument^(N-1)] * [unknown polynomial coefficients] = [known data point values] */    
	A = std::vector< std::vector<double> >( 2*NO_INTERPOLATION_POINTS, std::vector<double> ( 2*NO_INTERPOLATION_POINTS, 0.0) );
	for(std::vector<int>::size_type j = 0; j < (std::vector<int>::size_type) NO_INTERPOLATION_POINTS; j++){ // First top half of the rows corresponds to the function values.
		for(std::vector<int>::size_type i=0; i < A.size(); i++){
			A.at(j).at(i) = pow( interpolationArguments.at(j), i );
		};
	};
	for(std::vector<int>::size_type j = 0; j < (std::vector<int>::size_type) NO_INTERPOLATION_POINTS; j++){ // The rest (bottom half) correpsonds to derivative values.
		for(std::vector<int>::size_type i=1; i < A.size(); i++){
			A.at(NO_INTERPOLATION_POINTS+j).at(i) = i*pow(interpolationArguments.at(j), i-1);
		};
	};

	/* Decompose the matrix in preparation for finding of the interpolation coefficients. */
	luDecompositionIndices = std::vector<int>();
	LU = EquationsSolving::luDecomposition(A, luDecompositionIndices, luDecompositionNoSwaps); // A LU decomposed matrix A.

};

void ExternalSpaceObject::PropagateJDAY(std::vector<double>* posPtr, std::vector<double>* veloPtr, double JDAY, bool updateCurrentState){
	/* Find the state of the satellite at the specified Julain day and output its position and velocity at that epoch to currentPosPtr
	and currentVeloPtr. Do this by cubic interpolation between two of the known ephemeris points (use position and velocity data)
	in a dimensionless time interval [0,1].
	@param posPtr, veloPtr - pointers to std::vector<double> of length 3 that contain current positions and velocities of the object.
	@param JDAY - Julian day at which the object's state is to be computed.
	@param updateCurrentState - whether to update currentPos, currentVelo and currentEpochJDAY with the new values.
	@throws FatalError - when the desired Julain Epoch is outside the period spanned by the ephemeris table.
	*/
	
	std::vector<double>::iterator UpperInterpolationBound = upper_bound(EphemerisEpochs.begin(), EphemerisEpochs.end(), JDAY); // Iterator to the first epoch that is greater than JDAY - this marks the upper bound of the interpolation domain.
	std::vector<double>::iterator LowerInterpolationBound = std::prev( UpperInterpolationBound ); // This must be the iterator to the epoch that marks the lower bound of the interpolation domain.
	
	std::vector<int>::size_type UPPER_INDEX = UpperInterpolationBound - EphemerisEpochs.begin(); // Need to get indices as will have to access position and velocity vectors at the same locations as the iterators are pointing to.
	std::vector<int>::size_type LOWER_INDEX = LowerInterpolationBound - EphemerisEpochs.begin();

	if( (JDAY<EphemerisEpochs.at(0)) || (JDAY>EphemerisEpochs.back()) ){
		std::string msg = "Requested state vector for the object from ";
		msg.append( EPHEMERIS_FILE_NAME );
		msg.append( " ephemeris file that is outside the ephemeris timespan." );
		throw FatalError(msg,__FILE__,__LINE__);
	}else{
		double dimensionalisationConstant = (*UpperInterpolationBound - *LowerInterpolationBound)*JDAY_IN_SECONDS; // Multiply velcioty in km/sec by this to get to km/JDAY and dived for opposite covnersion i.e. km/JDAY to km/sec.
		double dimensionlessEpoch = (JDAY - *LowerInterpolationBound)/(*UpperInterpolationBound - *LowerInterpolationBound); // Dimensionless location in the current interpolation interval.

		/* Get the right-hand side values for interpolation (new positions and velocities in the current interval). */
		// First entry on the RHS is the position at the beginning of the interpolation interval.
		rhsX.at(0) = EphemerisPositions.at(LOWER_INDEX).at(0);
		rhsY.at(0) = EphemerisPositions.at(LOWER_INDEX).at(1);
		rhsZ.at(0) = EphemerisPositions.at(LOWER_INDEX).at(2);
		// Second is the position at the end of the interval.
		rhsX.at(1) = EphemerisPositions.at(UPPER_INDEX).at(0);
		rhsY.at(1) = EphemerisPositions.at(UPPER_INDEX).at(1);
		rhsZ.at(1) = EphemerisPositions.at(UPPER_INDEX).at(2);
		// Third is the velocity at the beginning of the interval.
		rhsX.at(2) = EphemerisVelocities.at(LOWER_INDEX).at(0) * dimensionalisationConstant; // Change velocity to km/JDAY and non-dimensionalise.
		rhsY.at(2) = EphemerisVelocities.at(LOWER_INDEX).at(1) * dimensionalisationConstant;
		rhsZ.at(2) = EphemerisVelocities.at(LOWER_INDEX).at(2) * dimensionalisationConstant;
		// And finally the velocity at the end of the interval.
		rhsX.at(3) = EphemerisVelocities.at(UPPER_INDEX).at(0) * dimensionalisationConstant;
		rhsY.at(3) = EphemerisVelocities.at(UPPER_INDEX).at(1) * dimensionalisationConstant;
		rhsZ.at(3) = EphemerisVelocities.at(UPPER_INDEX).at(2) * dimensionalisationConstant;

		/* Solve the system of equations to find the interpolating polynomials' coefficients. */
		EquationsSolving::luSubst(&LU, &luDecompositionIndices, &rhsX); // The interpolating coefficients will be stored in the rhs vectors.
		EquationsSolving::luSubst(&LU, &luDecompositionIndices, &rhsY);
		EquationsSolving::luSubst(&LU, &luDecompositionIndices, &rhsZ);

		/* Get the interpolated state vectors. */
		posPtr->at(0) = rhsX.at(0) + rhsX.at(1)*dimensionlessEpoch + rhsX.at(2)*dimensionlessEpoch*dimensionlessEpoch + rhsX.at(3)*dimensionlessEpoch*dimensionlessEpoch*dimensionlessEpoch;
		posPtr->at(1) = rhsY.at(0) + rhsY.at(1)*dimensionlessEpoch + rhsY.at(2)*dimensionlessEpoch*dimensionlessEpoch + rhsY.at(3)*dimensionlessEpoch*dimensionlessEpoch*dimensionlessEpoch;
		posPtr->at(2) = rhsZ.at(0) + rhsZ.at(1)*dimensionlessEpoch + rhsZ.at(2)*dimensionlessEpoch*dimensionlessEpoch + rhsZ.at(3)*dimensionlessEpoch*dimensionlessEpoch*dimensionlessEpoch;

		veloPtr->at(0) = ( rhsX.at(1) + rhsX.at(2)*dimensionlessEpoch*2.0 + rhsX.at(3)*dimensionlessEpoch*dimensionlessEpoch*3.0 )/dimensionalisationConstant; // Can get interpolated velocities in every direction directly by differentiating the 3rd order interpolating polynomial.
		veloPtr->at(1) = ( rhsY.at(1) + rhsY.at(2)*dimensionlessEpoch*2.0 + rhsY.at(3)*dimensionlessEpoch*dimensionlessEpoch*3.0 )/dimensionalisationConstant; // Go back to km/sec.
		veloPtr->at(2) = ( rhsZ.at(1) + rhsZ.at(2)*dimensionlessEpoch*2.0 + rhsZ.at(3)*dimensionlessEpoch*dimensionlessEpoch*3.0 )/dimensionalisationConstant;


		if(updateCurrentState){ // Update the state variables if desired.
			currentEpochJDAY = JDAY;
			for(std::vector<int>::size_type i = 0; i < currentPos.size(); i++){ // Both vectors have the same length here as well - they are Cartesian components.
				currentPos.at(i) = posPtr->at(i);
				currentVelo.at(i) = veloPtr->at(i);
			};
		};
	}; // The desired epoch is within the timespan of the ephemeris table.
};

void ExternalSpaceObject::Propagate(std::vector<double>* posPtr, std::vector<double>* veloPtr, int year, int month, int day, int hour, int minute, double second, bool updateCurrentState){
	/* Find the state of the satellite at the specified Universal Time epoch and output its position and velocity at that epoch to currentPosPtr
	and currentVeloPtr. Do this by cubic interpolation between two of the known ephemeris points (use position and velocity data).
	Specifically will convert the desired UT epoch to Julian Days and call @see PropagateJDAY.
	@param posPtr, veloPtr - pointers to std::vector<double> of length 3 that contain current positions and velocities of the object.
	@param JDAY - Julian day at which the object's state is to be computed.
	@param updateCurrentState - whether to update currentPos, currentVelo and currentEpochJDAY with the new values.
	@throws FatalError - when the desired Julain Epoch is outside the period spanned by the ephemeris table.
	*/
	double JDAY = currentEpochJDAY; // Initialise with the last new location.
	jday(year, month, day, hour, minute, second, JDAY); // Convert desired epoch to Julian days.
	
	PropagateJDAY( posPtr, veloPtr, JDAY, updateCurrentState); // Call the method that accepts Julian Day arguments.
};

void ExternalSpaceObject::CalculateRadius(std::vector<double>* posPtr){
	/* Update the CURRENT_RADIUS based on the supplied Cartesian position. */
	CURRENT_RADIUS = VectorOperations::vectorMagnitude( posPtr ); // In km.
};

void ExternalSpaceObject::SetHardBodyRadius(double bodyRadius){
	/* Set the hardBodyRadius attribute to the specified value in m. */
	hardBodyRadius = bodyRadius/1000.; // Convert to km.
};

void ExternalSpaceObject::SetCovarianceMatrixRTC(std::vector< std::vector<double> > CovarianceMatrixRTC){
	/* Set the FullCovarianceMatrixRTC attribute. It contains a 6x6 position and velocity covariance matrix in km and km/sec and is defined in
	radial - in-track - cross-track reference frame at the epoch of this object.
	@param CovarianceMatrixRTC - 6x6 position and velocity covariance matrix in km and km/sec defined in the radial - in-track - cross-track reference frame at the same epoch as this ExternalSpaceObject, i.e. EpochJDAY.
	*/
	FullCovarianceMatrixRTC = CovarianceMatrixRTC; // Set the full 6x6 covariance matrix.

	// Also set the 3x3 position covariance that will be used to compute the collision probability.
	PositionCovarianceMatrixRTC = std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0) ); // Initialise.
	PositionCovarianceMatrixRTC.at(0).at(0) = CovarianceMatrixRTC.at(0).at(0); PositionCovarianceMatrixRTC.at(0).at(1) = CovarianceMatrixRTC.at(0).at(1); PositionCovarianceMatrixRTC.at(0).at(2) = CovarianceMatrixRTC.at(0).at(2);
	PositionCovarianceMatrixRTC.at(1).at(0) = CovarianceMatrixRTC.at(1).at(0); PositionCovarianceMatrixRTC.at(1).at(1) = CovarianceMatrixRTC.at(1).at(1); PositionCovarianceMatrixRTC.at(1).at(2) = CovarianceMatrixRTC.at(1).at(2);
	PositionCovarianceMatrixRTC.at(2).at(0) = CovarianceMatrixRTC.at(2).at(0); PositionCovarianceMatrixRTC.at(2).at(1) = CovarianceMatrixRTC.at(2).at(1); PositionCovarianceMatrixRTC.at(2).at(2) = CovarianceMatrixRTC.at(2).at(2);

	// Get 3x3 velocity covariance to seed TLEs from.
	std::vector< std::vector<double> > VelocityCovarianceMatrixRTC = std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0) );
	VelocityCovarianceMatrixRTC.at(0).at(0) = CovarianceMatrixRTC.at(3).at(3); VelocityCovarianceMatrixRTC.at(0).at(1) = CovarianceMatrixRTC.at(3).at(4); VelocityCovarianceMatrixRTC.at(0).at(2) = CovarianceMatrixRTC.at(3).at(5);
	VelocityCovarianceMatrixRTC.at(1).at(0) = CovarianceMatrixRTC.at(4).at(3); VelocityCovarianceMatrixRTC.at(1).at(1) = CovarianceMatrixRTC.at(4).at(4); VelocityCovarianceMatrixRTC.at(1).at(2) = CovarianceMatrixRTC.at(4).at(5);
	VelocityCovarianceMatrixRTC.at(2).at(0) = CovarianceMatrixRTC.at(5).at(3); VelocityCovarianceMatrixRTC.at(2).at(1) = CovarianceMatrixRTC.at(5).at(4); VelocityCovarianceMatrixRTC.at(2).at(2) = CovarianceMatrixRTC.at(5).at(5);
	//UNDONE: should find a way to propagate the covariance somehow? Maybe use Monte Carlo as in SpaceObject?
};

std::vector< std::vector<double> > ExternalSpaceObject::InertialToRTC(std::vector<double>* positionPtr, std::vector<double>* velocityPtr){
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
	VectorOperations::crossProduct( &transverseUnit_inertial, &crossUnit_inertial, &radialUnit_inertial);
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

void ExternalSpaceObject::LoadEphemerisTable( std::string fullFileName ){
	/* Read a STK v 10 ephemeris file that contains the ephemeris of this object. The epochs in Julian Days of all the ephemeris points will
	be stored in EphemerisEpochs attribute, the Cartesian positions in km in EphemerisPositions, and the Cartesian velocities
	in km/sec in EphemerisVelocities. The state vectors are defined in the True Equator Mean Equinox reference frame.
	@param - name of a file that contains Julian day epoch, Cartesian positions and velocities in meteres and metres per second in TEME reference frame.
	*/
	EphemerisEpochs = std::vector<double>(); // Initialise the vectors that will hold the ephemeris points.
	EphemerisPositions = std::vector< std::vector<double> >();
	EphemerisVelocities = std::vector< std::vector<double> >();

	/* Parse the ephemeris file. */
	std::ifstream TLEfileStream(fullFileName, std::ifstream::in);
	
	std::string EphemerisLine; // Currently read ephemeris point.

	int noEphemPoints; // Counter of the total number of ephemeris points read.
	double JDAYepoch;
	int EphemerisPointsRead = 0; // Counter.
	double temp; // Temporary number read from the file.

	while( std::getline(TLEfileStream, EphemerisLine) ){
		std::istringstream iss(EphemerisLine);

		/* Get information about the ephemeris file being parsed from the header and parse the ephemeris points. */
		if( EphemerisLine.find( "NumberOfEphemerisPoints" ) == 0 ){
			std::string noEphemPointsString = EphemerisLine.substr( EphemerisLine.find(" ")+1 ); // Don't want the white space before the number of points.
			std::basic_istringstream<char> noEphemPointsSS( noEphemPointsString ); // Convert from std::string to an int.
			noEphemPointsSS >> noEphemPoints;
			EphemerisEpochs.resize(noEphemPoints); // Allocate memory for the ephemeris points.
			EphemerisPositions.resize(noEphemPoints);
			EphemerisVelocities.resize(noEphemPoints);
			for(std::vector<int>::size_type i=0; i<EphemerisEpochs.size(); i++){
				EphemerisPositions.at(i).resize(3); // Resize every entry to be able to hold epoch, 3 Cartesian poisiton and 3 Cartesian velocity compontents (Julian Days, metres, and metres per second, respectively).
				EphemerisVelocities.at(i).resize(3);
			}
		}else if( EphemerisLine.find( "# Time of first point:" ) == 0 ){
			std::string JDAYepochString = EphemerisLine.substr( EphemerisLine.find("UTCG =")+6, 14 ); // Only get this part of the string that contains the epoch of the first point.
			std::basic_istringstream<char> JDAYepochSS( JDAYepochString );
			JDAYepochSS >> JDAYepoch; // Find the epoch of the first point - will get offset times from this epoch in seconds for every ephemeris point.
		}else if( EphemerisLine.size() != 0 ){ // Don't try to call .at(0) for empty strings.
			if( isdigit(EphemerisLine.at(0)) ){ // Only actual ephemeris points will start with digits.
				std::basic_istringstream<char> ephemPointSS( EphemerisLine ); // Parse this ephemeris point - always 7 numbers in scientific notation but don't preassume anyting.
				for(int counter=0; counter<7; counter++){ // Read the epoch, position, and velocity for this point.
					if(counter==0){
						ephemPointSS >> temp; // Change the epoch to Julian Days from the time offset in seconds from the first point.
						EphemerisEpochs.at(EphemerisPointsRead) = JDAYepoch + temp*SECOND_IN_JDAYS;
					}else if( (counter>0) && (counter<4) ){ // Now parsing position components.
						ephemPointSS >> temp;
						EphemerisPositions.at(EphemerisPointsRead).at(counter-1) = temp/1000.0; // Change to km.
					}else if( (counter>3) && (counter<7) ){ // Now parsing velocity compontents.
						ephemPointSS >> temp;
						EphemerisVelocities.at(EphemerisPointsRead).at(counter-4) = temp/1000.0;
					}
				}
				EphemerisPointsRead += 1; // Make sure this ephemeris point won't be overwritten.
			}
		}
	};

	int LastLineInfo = EphemerisLine.find("END Ephemeris");
	if( (LastLineInfo==0) || (EphemerisLine.empty()) ){ // Reached the last line of the file.
		std::cout<<"Parsed entre ephemeris file "<<fullFileName<<" and read "<<EphemerisEpochs.size()<<" ephemeris points."<<std::endl;
	}else if( (LastLineInfo!=0) && (!EphemerisLine.empty()) ){ // Haven't reached the last line of the file. Last line may be empty.
		std::cerr<<"Haven't reached the end of the ephemeris file while reading the ephemeris from "<<fullFileName<<"."<<std::endl;
	}else if(noEphemPoints != (int) EphemerisEpochs.size() ){
		std::cerr<<"Haven't read as many ephemeris points as quoted in the file header for "<<fullFileName<<"."<<std::endl;
	};
}

ExternalSpaceObject::~ExternalSpaceObject(void){
 /* Deconstructor, does nothing special.*/
};

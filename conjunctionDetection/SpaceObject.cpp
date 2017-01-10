#include "SpaceObject.h"

SpaceObject::SpaceObject(void){
	/* Default constructor, does nothing special. */
};

SpaceObject::SpaceObject(std::vector<std::string> threeLineElement, int COV, std::string covpath, double bstarMultiplier, int wgsEllipsoid, double objectRadius, char opsMode){
	/* Create an object that represents something orbiting the Earth - a convenience interface for SGP4 propagator by D. Vallado.
	@param threeLineElement - a vector of strings that cointains the object name (line 0) and the first and second lines of the TLE file for a given object. Can have trailing end of line characters and must have the leading 1 and 2.
	@param COV - type of covariance to be used for the objects.
		1 - an uncertainty sphere 1:1:1 will be used and only maximum collision probability will be computed.
		2 - a number of three-line elements will be read for every object and covariance will be estimated for those using the method of V.P. Osweiler. It will be assumed that
		3LEs are in the files called SSC.txt where SSC if the catalogue number of every object and that the 3LEs are sorted oldest to last (the most recent one is last).
		3 - fixed covariance matrices based on the predicted performance of the European Space Surveillance System (40x200x100 m^2 in each objects RTC) will be used.
	@param covpath - path from which the TLEs will be read to estimate when COV is set to 2, ignored otherwise.
	@param bstarMultiplier - a factor by which the B* coefficient of a satellite will be multiplied by, used to simulate atmospheric density changes.
	@param wgsEllispoid - specifies which WGS Earth ellipsoid to use, possible options are 72, 84 and 721 (old WGS72).
	@param objectRadius - radius of the object in metres, used to compute the probability of collisions.
	@param opsMode - operation mode of the SGP4, either i (improved) or a (afspc).
	*/
	char OPSMODE = opsMode; // Improved operation mode of SGP4, should result in a smoother behaviour.
	earthEllipsoid = wgsEllipsoid;
	BstarMultiplier = bstarMultiplier;

	/* Initialise SGP4-specific attributes. */
	revnum=0; elnum=0;
	
	// These are defined in SGP4EXT, will need to uncomment them when not using it.
	//const double deg2rad = PI/180.0;         // Conversion factor. 0.0174532925199433
	//const double xpdotp = 1440.0/(2.0*PI);  // Minutes in a dayper a full revolution = 229.1831180523293. Unit: min/rad
	sgp4Sat = elsetrec(); // SGP4 satellite struct that is used for propagation using the SGP4.
	
	/* Define the desired EGS ellipsoid. */
	if(wgsEllipsoid==84){
		gravityModel = wgs84;
	}else if(wgsEllipsoid==72){
		gravityModel = wgs84;
	}else if(wgsEllipsoid==721){
		gravityModel = wgs72;
	};
	getgravconst( gravityModel, tumin, mu, Re, xke, j2, j3, j4, j3oj2 ); // Get values for these variables.

	/* Convert the TLE into the expected format. */
	char TLE1[130];	char TLE2[130]; // TLE lines - char arrays of the appropriate length that are expected by sgp4init.

	#ifdef _MSC_VER // Depending on the compiler being used utilise the appropriate function to copy the TLE lines.
		strcpy_s(TLE1, threeLineElement.at(1).c_str());
		strcpy_s(TLE2, threeLineElement.at(2).c_str());
	#else
		std::strcpy(TLE1, threeLineElement.at(1).c_str());
		std::strcpy(TLE2, threeLineElement.at(2).c_str());
	#endif
	
	/* Set the implied decimal points since doing a formated read
     * fixes for bad input data values (missing, ...). */
    for (int j = 10; j <= 15; j++){
        if (TLE1[j] == ' ')
			TLE1[j] = '_';
	};
	if (TLE1[44] != ' ')
		TLE1[43] = TLE1[44];
	TLE1[44] = '.';
	if (TLE1[7] == ' ')
		TLE1[7] = 'U';
    if (TLE1[9] == ' ')
		TLE1[9] = '.';
	for (int j = 45; j <= 49; j++){
		if (TLE1[j] == ' ')
			TLE1[j] = '0';
	};
    if (TLE1[51] == ' ')
		TLE1[51] = '0';
    if (TLE1[53] != ' ')
		TLE1[52] = TLE1[53];
	TLE1[53] = '.';
    TLE2[25] = '.';
    for (int j = 26; j <= 32; j++){
        if (TLE2[j] == ' ')
			TLE2[j] = '0';
	};
    if (TLE1[62] == ' ')
        TLE1[62] = '0';
    if (TLE1[68] == ' ')
        TLE1[68] = '0';

	/* Parse the TLE. */
	#ifdef _MSC_VER // Depending on the compiler being used utilise the appropriate function to copy the TLE lines.
		sscanf_s(TLE1,"%2d %5ld %1c %10s %2d %12lf %11lf %7lf %2d %7lf %2d %2d %6ld ",&cardnumb,&sgp4Sat.satnum,&classification, sizeof(char),intldesg, 11*sizeof(char),&sgp4Sat.epochyr,&sgp4Sat.epochdays,&sgp4Sat.ndot,&sgp4Sat.nddot,&nexp,&sgp4Sat.bstar,&ibexp,&numb,&elnum);
		sscanf_s(TLE2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld %lf %lf %ld \n",&cardnumb,&sgp4Sat.satnum,&sgp4Sat.inclo,&sgp4Sat.nodeo,&sgp4Sat.ecco,&sgp4Sat.argpo,&sgp4Sat.mo,&sgp4Sat.no,&revnum);
	#else
		sscanf(TLE1,"%2d %5ld %1c %10s %2d %12lf %11lf %7lf %2d %7lf %2d %2d %6ld ",&cardnumb,&sgp4Sat.satnum,&classification,intldesg,&sgp4Sat.epochyr,&sgp4Sat.epochdays,&sgp4Sat.ndot,&sgp4Sat.nddot,&nexp,&sgp4Sat.bstar,&ibexp,&numb,&elnum);
		sscanf(TLE2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld %lf %lf %ld \n",&cardnumb,&sgp4Sat.satnum,&sgp4Sat.inclo,&sgp4Sat.nodeo,&sgp4Sat.ecco,&sgp4Sat.argpo,&sgp4Sat.mo,&sgp4Sat.no,&revnum);
	#endif
	
	std::stringstream strstream; // Record the NORAD ID as a string as well.
	strstream << sgp4Sat.satnum; strstream >> NORAD_ID;
		
	/* Iniitalise the SGP4 propagator variables and the propagator itself for this object. */
	// Find the TLE epoch in Julian days.
    int year = 2000 + sgp4Sat.epochyr; // N.B. this won't work for historic TLEs from before 2000.

	int epMon, epDay, epHr, epMin; double epSec; // TLE epoch components.
	days2mdhms(year, sgp4Sat.epochdays, epMon, epDay, epHr, epMin, epSec);
	jday(year, epMon, epDay, epHr, epMin, epSec, sgp4Sat.jdsatepoch); 

	// Find no, ndot, nddot.
    sgp4Sat.no   = sgp4Sat.no / xpdotp; // When using SGP4EXT this will already be in correct units, anyway multiply by: rad/min
    sgp4Sat.nddot= sgp4Sat.nddot * pow(10.0, nexp);
    sgp4Sat.bstar= sgp4Sat.bstar * pow(10.0, ibexp) * bstarMultiplier; // Multiply by the factor that allows variations in solar activity to be synthesised.
    
	// Convert to sgp4 units.
    sgp4Sat.a    = pow( sgp4Sat.no*tumin , (-2.0/3.0) );
    sgp4Sat.ndot = sgp4Sat.ndot  / (xpdotp*1440.0);  //* ? * minperday
    sgp4Sat.nddot= sgp4Sat.nddot / (xpdotp*1440.0*1440);
    
	// Find standard orbital elements.
    sgp4Sat.inclo = sgp4Sat.inclo  * deg2rad;
    sgp4Sat.nodeo = sgp4Sat.nodeo  * deg2rad;
    sgp4Sat.argpo = sgp4Sat.argpo  * deg2rad;
    sgp4Sat.mo    = sgp4Sat.mo     * deg2rad;

	// Iniitlaise the SGP4 for this TLE.
	sgp4init( gravityModel, OPSMODE, sgp4Sat.satnum, sgp4Sat.jdsatepoch-2433281.5, sgp4Sat.bstar, sgp4Sat.ecco, sgp4Sat.argpo, sgp4Sat.inclo, sgp4Sat.mo, sgp4Sat.no, sgp4Sat.nodeo, sgp4Sat);

	/* Call the propagator to get the initial state. */
    sgp4(gravityModel, sgp4Sat,  0.0, r0,  v0); // Third argument is time since epoch in minutes - 0.0 gives the initial position and velocity.
	currentPos.resize( sizeof(r0)/sizeof(r0[0]) ); currentVelo.resize( sizeof(v0)/sizeof(v0[0]) ); // Resize the vectors to be of appropriate length.
	currentPos = std::vector<double>(r0, r0+3); currentVelo = std::vector<double>(v0, v0+3); // Record the values of the initial position and velocity in the vectors.
	currentEpochJDAY = sgp4Sat.jdsatepoch;
	TLEepochJDAY = sgp4Sat.jdsatepoch; // In Julian Days.

	/* Get other values of interest, mianly orbital elements. */
	double SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper;
	rv2coe(r0, v0, mu, SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper );
	SemiMajorAxis = SMA; Eccentricity = ECC; Inclination = INCL/deg2rad; MeanAnomaly = MA/deg2rad; LongAscendingNode = LAN/deg2rad; ArgumentOfPerigee = ARGP/deg2rad; AltitudeOfPerigee = SMA*(1.0 - ECC) - Re; // Store the elements at the TLE epoch. Need them for assuming covariance (it varies depending on orbital regime).
	CURRENT_PERIGEE_RADIUS = SMA*(1.0 - ECC); // In km.
	CURRENT_APOGEE_RADIUS = SMA*(1.0 + ECC); // In km.

	hardBodyRadius = objectRadius/1000.; // Default value in km.

	switch(COV) {
	case 1: // Spherical covariance matrices, only care about the aspect ratio.
		FullCovarianceMatrixRTC = std::vector< std::vector<double> >(6, std::vector<double>(6, 0.0)); // Initialise the covariance matrices.
		PositionCovarianceMatrixRTC = std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0));

		// Initialise the diagonal entries of the covariance matrices to simulate the desired aspect ratio.
		FullCovarianceMatrixRTC.at(0).at(0)=1.0; FullCovarianceMatrixRTC.at(1).at(1)=1.0; FullCovarianceMatrixRTC.at(2).at(2)=1.0; FullCovarianceMatrixRTC.at(3).at(3)=1.0; FullCovarianceMatrixRTC.at(4).at(4)=1.0; FullCovarianceMatrixRTC.at(5).at(5)=1.0;
		PositionCovarianceMatrixRTC.at(0).at(0)=1.0; PositionCovarianceMatrixRTC.at(1).at(1)=1.0; PositionCovarianceMatrixRTC.at(2).at(2)=1.0;
		break;
	case 2: // Estimate covaraince for this SpaceObject using V.P. Osweiler's approach.
		ReadThreeLineElements( covpath ); // Read three line elements from a file - need them to evaluate covariance at any epoch.
		FullCovarianceMatrixRTC = std::vector< std::vector<double> >(6, std::vector<double>(6, 0.0)); // Initialise the covariance matrices.
		PositionCovarianceMatrixRTC = std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0));
		break;
	case 3: // ESSS performance, time-invariant covariance.
		FullCovarianceMatrixRTC = std::vector< std::vector<double> >(6, std::vector<double>(6, 0.0)); // Initialise the covariance matrices.
		PositionCovarianceMatrixRTC = std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0));
		// 40x200x100 m^2 in RTC.																															// Choose arbitrary values for velcoity covaraince, it isn't used.
		FullCovarianceMatrixRTC.at(0).at(0)=0.04*0.04;     FullCovarianceMatrixRTC.at(1).at(1)=0.2*0.2;     FullCovarianceMatrixRTC.at(2).at(2)=0.1*0.1;    FullCovarianceMatrixRTC.at(3).at(3)=1.0; FullCovarianceMatrixRTC.at(4).at(4)=1.0; FullCovarianceMatrixRTC.at(5).at(5)=1.0;
		PositionCovarianceMatrixRTC.at(0).at(0)=0.04*0.4; PositionCovarianceMatrixRTC.at(1).at(1)=0.2*0.2; PositionCovarianceMatrixRTC.at(2).at(2)=0.1*0.1;
		break;
	case 4: // For debugging only - set all the standard deviations to 1.0 to be able to compare with STK CAT (it also rotates the ellipsoids so they need to be spheres).
		FullCovarianceMatrixRTC = std::vector< std::vector<double> >(6, std::vector<double>(6, 0.0)); // Initialise the covariance matrices.
		PositionCovarianceMatrixRTC = std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0));
		FullCovarianceMatrixRTC.at(0).at(0)=1.0;     FullCovarianceMatrixRTC.at(1).at(1)=1.0;     FullCovarianceMatrixRTC.at(2).at(2)=1.0;    FullCovarianceMatrixRTC.at(3).at(3)=1.0; FullCovarianceMatrixRTC.at(4).at(4)=1.0; FullCovarianceMatrixRTC.at(5).at(5)=1.0;
		PositionCovarianceMatrixRTC.at(0).at(0)=1.0; PositionCovarianceMatrixRTC.at(1).at(1)=1.0; PositionCovarianceMatrixRTC.at(2).at(2)=1.0;
		break;
	case 5: // Don't use covaraince matrices - won't compute Pc at all.
		break;
	case 6: // The user has specified variances to use externally.
		FullCovarianceMatrixRTC = std::vector< std::vector<double> >(6, std::vector<double>(6, 0.0)); // Initialise the covariance matrices.
		PositionCovarianceMatrixRTC = std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0));
		// Will set the variances via a set method later.
		break;
	}//END switch covariance type
};

void SpaceObject::PropagateJDAY(std::vector<double>* posPtr, std::vector<double>* veloPtr, double JDAY, bool updateCurrentState){
	/* Propagate the satellite to the specified Julain day and output its position and velocity at that epoch to currentPosPtr
	and currentVeloPtr.
	@param posPtr, veloPtr - pointers to std::vector<double> of length 3 that contain current positions and velocities of the object.
	@param JDAY - Julian day at which the object's state is to be computed.
	@param updateCurrentState - whether to update currentPos, currentVelo and currentEpochJDAY with the new values.
	*/
	double r[3]; double v[3]; // Create arrays that SGP4 expects
	double minutesSinceEpoch = (JDAY-TLEepochJDAY)*1440.0;
    sgp4(gravityModel, sgp4Sat,  minutesSinceEpoch, r,  v);
	
	for(std::vector<int>::size_type i = 0; i < currentPos.size(); i++){ // Record the position and velocity in the vectors.
			posPtr->at(i) = r[i];
			veloPtr->at(i) = v[i];
		};

	if(updateCurrentState){ // Update the state variables if desired.
		currentEpochJDAY = JDAY;
		for(std::vector<int>::size_type i = 0; i < currentPos.size(); i++){ // Both vectors have the same length here as well - they are Cartesian components.
			currentPos.at(i) = r[i];
			currentVelo.at(i) = v[i];
		};
	};
};

void SpaceObject::Propagate(std::vector<double>* posPtr, std::vector<double>* veloPtr, int year, int month, int day, int hour, int minute, double second, bool updateCurrentState){
	/* Propagate the satellite to the specified Universal Time epoch and output its position and velocity at that epoch to currentPosPtr
	and currentVeloPtr.
	@param posPtr, veloPtr - pointers to std::vector<double> of length 3 that contain current positions and velocities of the object.
	@param year, month, day, hour, minute, second - Universal Time at which the object's state is to be computed.
	@param updateCurrentState - whether to update currentPos, currentVelo and currentEpochJDAY with the new values.
	*/
	double JDAY = currentEpochJDAY; // Initialise with the last new location.
	jday(year, month, day, hour, minute, second, JDAY); // Convert desired epoch to Julian days.
	double minutesSinceEpoch = (JDAY-TLEepochJDAY)*1440.0; // And minutes since TLE epoch.
	
	double r[3]; double v[3]; // Propagate using arrays that SGP4 expects
	sgp4(gravityModel, sgp4Sat,  minutesSinceEpoch, r,  v);

	for(std::vector<int>::size_type i = 0; i < currentPos.size(); i++){ // Record the position and velocity in the vectors.
			posPtr->at(i) = r[i];
			veloPtr->at(i) = v[i];
		};

	if(updateCurrentState){ // Update the state variables if desired.
		currentEpochJDAY = JDAY;
		for(std::vector<int>::size_type i = 0; i < currentPos.size(); i++){ // Both vectors have the same length - they are Cartesian components.
			currentPos[i] = r[i];
			currentVelo[i] = v[i];
		};
	};
};

void SpaceObject::CalculatePerigeeRadius(void){
	/* Update the currentPerigeeRadius based on currentPos and currentVelo values. */
	double SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper;
	rv2coe(&currentPos[0], &currentVelo[0], mu, SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper );
	CURRENT_PERIGEE_RADIUS = SMA*(1.0 - ECC); // In km.
};

void SpaceObject::CalculateApogeeRadius(void){
	/* Update the currentApogeeRadius based on currentPos and currentVelo values. */
	double SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper;
	rv2coe(&currentPos[0], &currentVelo[0], mu, SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper );
	CURRENT_APOGEE_RADIUS = SMA*(1.0 + ECC); // In km.
};

void SpaceObject::SetHardBodyRadius(double bodyRadius){
	/* Set the hardBodyRadius attribute to the specified value in m. */
	hardBodyRadius = bodyRadius/1000.; // Convert to km.
};

void SpaceObject::SetPositionVariances(double radialVariance, double inTrackVariance, double crossTrackVariance){
	/* Set position variances in the radial, in-track, cross-track covariance 
	matrices that are used to calculate the collision probabilities. Won't change
	the velocity covariance or the off-diagonal of the position covariance matrix.
	@param radialVariance - radial variance in m^2
	@param inTrackVariance - in-track variance in m^2
	@param radialVariance - cross-track variance in m^2
	 */
	FullCovarianceMatrixRTC.at(0).at(0)=radialVariance/1000./1000.; // Convert to km^2 from m^2.
	FullCovarianceMatrixRTC.at(1).at(1)=inTrackVariance/1000./1000.;
	FullCovarianceMatrixRTC.at(2).at(2)=crossTrackVariance/1000./1000.;
	// Don't touch the velocity covariance.
	//FullCovarianceMatrixRTC.at(3).at(3)=?;
	//FullCovarianceMatrixRTC.at(4).at(4)=?;
	//FullCovarianceMatrixRTC.at(5).at(5)=?;

	PositionCovarianceMatrixRTC.at(0).at(0)=radialVariance/1000./1000.;
	PositionCovarianceMatrixRTC.at(1).at(1)=inTrackVariance/1000./1000.;
	PositionCovarianceMatrixRTC.at(2).at(2)=crossTrackVariance/1000./1000.;
};

//TODO: work on creating TLEs in SetCovarianceMatrix to be able to propagate the set covariance to any epoch.
void SpaceObject::SetCovarianceMatrixRTC(std::vector< std::vector<double> > CovarianceMatrixRTC, int noTLEs){
	/* Set the FullCovarianceMatrixRTC attribute. It contains a 6x6 position and velocity covariance matrix in km and km/sec and is defined in
	radial - in-track - cross-track reference frame at the epoch of this object's TLE. Also initialise a set of TLEs stored in CovarianceTLEs
	that will be used to propagate the covariance to any epoch in a Monte Carlo fashion. Ignore the cross-terms when spawning these TLEs, only use
	standard deviations (i.e. assume uncorrelated errors at the epoch). Also ignore the velocity covariance when seeding the randomised TLEs.
	Assume that the epoch of the covariance is equal to the epoch of this SpaceObject.
	@param CovarianceMatrixRTC - 6x6 position and velocity covariance matrix in km and km/sec defined in the radial - in-track - cross-track reference frame at the same epoch as this SpaceObject, i.e. TLEepochJDAY.
	@param noTLEs - number of the TLEs that will be generated from CovarainceMatrixRTC's standard deviations and used to propagate the covariance in a Monte Carlo fashion.
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

	/* Rotate the position covariance matrix to TEME to make it easier to spawn random positions and hence also TLEs - use these to propagate the covariance in a Monte Carlo fashion. */
	std::vector<double> tempPos(3, 0.0); // Get the state vector of epoch of this object.
	std::vector<double> tempVelo(3, 0.0);
	PropagateJDAY(&tempPos, &tempVelo, TLEepochJDAY); // Here assume that the covariance epoch is to be set at the epoch of this space object - this defines the RTC-TEME rotation.

	std::vector< std::vector<double> > TEME2RTC = InertialToRTC(&tempPos, &tempVelo); // Rotation matrix from prime TLE's VNC to TEME axes.

	// RTC to TEME rotation - formulate necessary matrices and then multiply out for rotation.
	std::vector< std::vector<double> > TEME2RTCinverse = VectorOperations::inverseThreeByThree(&TEME2RTC); // Need this to rotate from VNC to TEME.
	std::vector< std::vector<double> > TEME2RTCtranspose = VectorOperations::transposeMatrix(&TEME2RTC); // Also need the transpose when rotating covaraince matrices.
	std::vector< std::vector<double> > TEME2RTCtransposeInverse = VectorOperations::inverseThreeByThree(&TEME2RTCtranspose);
	
	std::vector< std::vector<double> > tempMat = VectorOperations::matrixMultiplyByMatrix(&PositionCovarianceMatrixRTC, &TEME2RTCtransposeInverse); // The operation is RotMat^-1 CovMat RotMat.T^-1.
	std::vector< std::vector<double> > TEMEposCovaraince = VectorOperations::matrixMultiplyByMatrix(&TEME2RTCinverse, &tempMat); // Finished rotations to TEME.
	
	tempMat = VectorOperations::matrixMultiplyByMatrix(&VelocityCovarianceMatrixRTC, &TEME2RTCtransposeInverse); // Rotate the velcity covariance as well.
	std::vector< std::vector<double> > TEMEveloCovaraince = VectorOperations::matrixMultiplyByMatrix(&TEME2RTCinverse, &tempMat);

	/* Spawn random positions around the intiial position of this SpaceObject, hence generate TLEs that will be used to propagate the covariance matrix. */
	std::default_random_engine RandomGenerator; // Will generate random numbers used to seed the TLEs' positions.
	
	//TODO: get standard deviations from the diagonalised covariance matrix here if setting non-diagonal covariance matrices is to be supported.
	//double A[3][3]={ // Will be changed by the algorithm.
 //   {-1,2,2},
	//{2,2,-1},
	//{2,-1,2}
	//};
	//double Q[3][3]={ // Outputs for the eigenvalue and eigenvector finding algorithm. Eigenvectors are in columns.
 //   {0,0,0},
 //   {0,0,0},
 //   {0,0,0}
	//};
	//double w[3]={0.,0.,0.};

	//int retValue = dsyevv3(A, Q, w);

	std::normal_distribution<double> XDistribution(0.0, std::sqrt(TEMEposCovaraince.at(0).at(0)) ); // Seeds random position offsets in the X-direction according to the standard deviation of the set covaraince matrix.
	std::normal_distribution<double> YDistribution(0.0, std::sqrt(TEMEposCovaraince.at(1).at(1)) );
	std::normal_distribution<double> ZDistribution(0.0, std::sqrt(TEMEposCovaraince.at(2).at(2)) );

	std::normal_distribution<double> VxDistribution; // Seeds random velocity offsets in the X-direction according to the standard deviation of the set covaraince matrix.
	std::normal_distribution<double> VyDistribution;
	std::normal_distribution<double> VzDistribution;

	if( TEMEveloCovaraince.at(0).at(0) != 0.0 ){ // Trying to generate random numbers from a normal distribution with 0.0 standard deviation will cause problems.
		VxDistribution = std::normal_distribution<double>(0.0, std::sqrt(TEMEveloCovaraince.at(0).at(0)) );
	}
	if( TEMEveloCovaraince.at(1).at(1) != 0.0 ){
		VyDistribution = std::normal_distribution<double>(0.0, std::sqrt(TEMEveloCovaraince.at(1).at(1)) );
	}
	if( TEMEveloCovaraince.at(2).at(2) != 0.0 ){
		VzDistribution = std::normal_distribution<double>(0.0, std::sqrt(TEMEveloCovaraince.at(2).at(2)) );
	}

	double XOffset, YOffset, ZOffset, VxOffset=0.0, VyOffset=0.0, VzOffset=0.0; // Random distances and velocities by which to offset the current randomly seeded covariance TLE from the original one in the RTC frame.

	/* Find TLEs from the random positions that come from the covariance matrix to be set. */
	CovarianceTLEs = std::vector<elsetrec>(); // Re-set the TLEs - in case only one TLE was read from a file it's the TLE corresponding to this object and it should be the last entry of CovarianceTLEs. If the covariance is set in a different mode of operation don't have any TLEs in CovarainceTLEs anyway.
	for(int i=0; i<noTLEs; i++){
		//TODO: these randomised state vectors will be along the eigenvectors, not TEME's X, Y, and Z.
		XOffset = XDistribution( RandomGenerator ); // Get random perturbations to the position of the prime TLE derived based on the set covaraince matrix for this particular covariance TLE.
		YOffset = YDistribution( RandomGenerator );
		ZOffset = ZDistribution( RandomGenerator );
		if( TEMEveloCovaraince.at(0).at(0) != 0.0 ){ // Get random velocity perturbations if the covariance matrix for them is available.
			VxOffset = VxDistribution( RandomGenerator );
		}
		if( TEMEveloCovaraince.at(1).at(1) != 0.0 ){ 
			VyOffset = VyDistribution( RandomGenerator );
		}
		if( TEMEveloCovaraince.at(2).at(2) != 0.0 ){ 
			VzOffset = VzDistribution( RandomGenerator );
		}

		double randomisedTEMEposition[3] = {r0[0]+XOffset, r0[1]+YOffset, r0[2]+ZOffset}; // Slightly perturbed position of the prime TLE.
		double randomisedTEMEvelcotiy[3] = {v0[0]+VxOffset, v0[1]+VyOffset, v0[2]+VzOffset}; // Slightly perturbed position of the prime TLE.
		// elsetrec randomisedTLE = ReproduceTLE( randomisedTEMEposition, v0); // Use the velocity of this SpaceObject - approximation.
		//CovarianceTLEs.push_back(randomisedTLE);
		//TODO: generate a TLE  based on the randomised position vector here.
	}
	elsetrec temSGP4SAT = sgp4Sat; // Make a copy of this TLE's elsetrec, it's expected at the end of CovarainceTLEs.
	CovarianceTLEs.push_back( temSGP4SAT );
};

void SpaceObject::ReadThreeLineElements( std::string CovarianceThreeLineElementsDirectory, double rejectionThreshold ){
	/* Read three-line elements from a file called NORAD_ID.txt and save them to CovarianceTLEs. Use these to estimate the covariance using the method by
	V.P. Osweiler in method ComputeCovarianceOSW. Assume that the 3LEs are in the order oldest to newest i.e. that the most-recent one is last. If only
	one TLE is available for the object set covariance based on object's orbital regime of epoch according to data from 
	Flohrer, T. Assessment and categorisation of TLE orbit errors for the US SSN catalogue, 1 Jan 2008 snapshot.
	TLEs that have specific orbital energies very different from the linear least-squares fit will be rejected as suggested in 
	Patera, R.P., Space Event Detection Method, 2008.
	@param CovarianceThreeLineElementsDirectory - directory (relative or absolute) where the three-line elements' file from which the covariance will be
		estimated is located.
	@param rejectionThreshold - number of standard deviations that will be used to reject TLEs whose specific orbital energies are more than
		rejectionThreshold standard deviations away from the arithmetic mean.
	@throws FatalError - when too few TLEs are availbale to estimate the covariance for a given object.
	*/
	std::string tempString = CovarianceThreeLineElementsDirectory; // Read the file from this directory.
	tempString.append("/"); // Add the trailing backslash - Windows struggles to read it from input command line arguments.
	tempString.append(NORAD_ID); // Assume that the 3LEs from which to estimate the covariance are saved in files called NORAD_ID.txt. Make a copy of NORAD_ID attribute not to change it.
	const char* fileName = tempString.append(".txt").c_str();
	std::ifstream TLEfileStream(fileName, std::ifstream::in); // File streams that are necessary.

	CovarianceTLEs = std::vector<elsetrec>(); // TLEs that are representative of the covariance matrix and are used to propagate it in a Monte Carlo fashion.
	std::vector<elsetrec> TEMPCovarainceTLEs = std::vector<elsetrec>(); // Same as above but not filtered for outliers.
	std::vector<double> SpecificOrbitalEnergies = std::vector<double>(); // Specific orbital energies of non-filetered TLEs that are used to compute the interpolating polynomial.
	std::vector<double> TLEepochs = std::vector<double>(); // Epochs of the TLEs that are used to compute the interpolating polynomial.

	std::string TLEline; // Currently read TLE line.
	std::vector<std::string> currentTLE; currentTLE.resize(3); // Assembled TLE (second and third lines) and the object name (first line).

	int counterTLEs = 0; // Current line of the TLE, 0, 1 or 2.
	/* Read all the three-line elements from the file, those will be used to estimate and propagate the covariance matrix. */
	while( std::getline(TLEfileStream, TLEline) ){
		std::istringstream iss(TLEline);
		currentTLE.at(counterTLEs) = TLEline; // Add a new line.
		counterTLEs+=1;
		if(counterTLEs == 3){ // Read a whole TLE.
			counterTLEs = 0; // Re-set the counter.

			/* Convert the TLE into the expected format. */
			char TLE1[130];	char TLE2[130]; // TLE lines - char arrays of the appropriate length that are expected by sgp4init.
			elsetrec temporarySgp4Sat; // elsetrec for the curretnly-read TLE.

			#ifdef _MSC_VER // Depending on the compiler being used utilise the appropriate function to copy the TLE lines.
				strcpy_s(TLE1, currentTLE.at(1).c_str());
				strcpy_s(TLE2, currentTLE.at(2).c_str());
			#else
				std::strcpy(TLE1, currentTLE.at(1).c_str());
				std::strcpy(TLE2, currentTLE.at(2).c_str());
			#endif
	
			/* Set the implied decimal points since doing a formated read
			 * fixes for bad input data values (missing, ...). */
			for (int j = 10; j <= 15; j++){
				if (TLE1[j] == ' ')
					TLE1[j] = '_';
			};
			if (TLE1[44] != ' ')
				TLE1[43] = TLE1[44];
			TLE1[44] = '.';
			if (TLE1[7] == ' ')
				TLE1[7] = 'U';
			if (TLE1[9] == ' ')
				TLE1[9] = '.';
			for (int j = 45; j <= 49; j++){
				if (TLE1[j] == ' ')
					TLE1[j] = '0';
			};
			if (TLE1[51] == ' ')
				TLE1[51] = '0';
			if (TLE1[53] != ' ')
				TLE1[52] = TLE1[53];
			TLE1[53] = '.';
			TLE2[25] = '.';
			for (int j = 26; j <= 32; j++){
				if (TLE2[j] == ' ')
					TLE2[j] = '0';
			};
			if (TLE1[62] == ' ')
				TLE1[62] = '0';
			if (TLE1[68] == ' ')
				TLE1[68] = '0';

			/* Parse the TLE. */
			#ifdef _MSC_VER // Depending on the compiler being used utilise the appropriate function to copy the TLE lines.
				sscanf_s(TLE1,"%2d %5ld %1c %10s %2d %12lf %11lf %7lf %2d %7lf %2d %2d %6ld ",&cardnumb,&temporarySgp4Sat.satnum,&classification, sizeof(char),intldesg, 11*sizeof(char),&temporarySgp4Sat.epochyr,&temporarySgp4Sat.epochdays,&temporarySgp4Sat.ndot,&temporarySgp4Sat.nddot,&nexp,&temporarySgp4Sat.bstar,&ibexp,&numb,&elnum);
				sscanf_s(TLE2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld %lf %lf %ld \n",&cardnumb,&temporarySgp4Sat.satnum,&temporarySgp4Sat.inclo,&temporarySgp4Sat.nodeo,&temporarySgp4Sat.ecco,&temporarySgp4Sat.argpo,&temporarySgp4Sat.mo,&temporarySgp4Sat.no,&revnum);
			#else
				sscanf(TLE1,"%2d %5ld %1c %10s %2d %12lf %11lf %7lf %2d %7lf %2d %2d %6ld ",&cardnumb,&temporarySgp4Sat.satnum,&classification,intldesg,&temporarySgp4Sat.epochyr,&temporarySgp4Sat.epochdays,&temporarySgp4Sat.ndot,&temporarySgp4Sat.nddot,&nexp,&temporarySgp4Sat.bstar,&ibexp,&numb,&elnum);
				sscanf(TLE2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld %lf %lf %ld \n",&cardnumb,&temporarySgp4Sat.satnum,&temporarySgp4Sat.inclo,&temporarySgp4Sat.nodeo,&temporarySgp4Sat.ecco,&temporarySgp4Sat.argpo,&temporarySgp4Sat.mo,&temporarySgp4Sat.no,&revnum);
			#endif

			std::stringstream strstream; // Record the NORAD ID as a string as well.
			strstream << temporarySgp4Sat.satnum; strstream >> NORAD_ID;
		
			/* Iniitalise the SGP4 propagator variables and the propagator itself for this object. */
			// Find the TLE epoch in Julian days.
			int TLEyear = 2000 + temporarySgp4Sat.epochyr; // N.B. this won't work for historic TLEs from before 2000.

			int epMon, epDay, epHr, epMin; double epSec; // TLE epoch components.
			days2mdhms(TLEyear, temporarySgp4Sat.epochdays, epMon, epDay, epHr, epMin, epSec);
			jday(TLEyear, epMon, epDay, epHr, epMin, epSec, temporarySgp4Sat.jdsatepoch); 

			// Find no, ndot, nddot.
			temporarySgp4Sat.no   = temporarySgp4Sat.no / xpdotp; //* rad/min
			temporarySgp4Sat.nddot= temporarySgp4Sat.nddot * pow(10.0, nexp);
			temporarySgp4Sat.bstar= temporarySgp4Sat.bstar * pow(10.0, ibexp) * BstarMultiplier; // Multiply by the factor that allows variations in solar activity to be synthesised.
			
			// Convert to sgp4 units.
			temporarySgp4Sat.a    = pow( temporarySgp4Sat.no*tumin , (-2.0/3.0) );
			temporarySgp4Sat.ndot = temporarySgp4Sat.ndot  / (xpdotp*1440.0);  //* ? * minperday
			temporarySgp4Sat.nddot= temporarySgp4Sat.nddot / (xpdotp*1440.0*1440);
			
			// Find standard orbital elements.
			temporarySgp4Sat.inclo = temporarySgp4Sat.inclo  * deg2rad;
			temporarySgp4Sat.nodeo = temporarySgp4Sat.nodeo  * deg2rad;
			temporarySgp4Sat.argpo = temporarySgp4Sat.argpo  * deg2rad;
			temporarySgp4Sat.mo    = temporarySgp4Sat.mo     * deg2rad;

			// Iniitlaise the SGP4 for this TLE.
			sgp4init( gravityModel, OPSMODE, temporarySgp4Sat.satnum, temporarySgp4Sat.jdsatepoch-2433281.5, temporarySgp4Sat.bstar, temporarySgp4Sat.ecco, temporarySgp4Sat.argpo, temporarySgp4Sat.inclo, temporarySgp4Sat.mo, temporarySgp4Sat.no, temporarySgp4Sat.nodeo, temporarySgp4Sat);
			
			// Find the state vector of this TLE at the epoch of this SpaceObject.
			double r[3]; double v[3]; // Create arrays that SGP4 expects.
			double minutesSinceEpoch = (temporarySgp4Sat.jdsatepoch-TLEepochJDAY)*1440.0; // Minutes since epoch of this TLE to propagate it to the epoch of this SpaceObject.
			sgp4(gravityModel, temporarySgp4Sat, minutesSinceEpoch, r,  v);

			// Find classical orbital elements of this TLE at the epoch of this SpaceObject.
			double SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper;
			rv2coe(r, v, mu, SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper );

			SpecificOrbitalEnergies.push_back( 1.0*mu/(2.0*SMA) ); // Record the specific energy of this TLE, km^2 sec^-2.
			TLEepochs.push_back( temporarySgp4Sat.jdsatepoch ); 
			TEMPCovarainceTLEs.push_back( temporarySgp4Sat ); // Save this TLE for subsequent filtering.
		} // If read a full three-line element.
	} // While still reading the 3LE file.

	// Build up the matrix of coefficients that gives the over-defined system of equations.
	std::vector< std::vector<double> > A = std::vector< std::vector<double> >(SpecificOrbitalEnergies.size(), std::vector<double>(SpecificOrbitalEnergies.size(), 0.0)); // SpecificOrbitalEnergies.size() gives number of the points i.e. the size of the matrix of coefficients.
	for(std::vector<int>::size_type j = 0; j < TLEepochs.size(); j++){
		for(std::vector<int>::size_type i=0; i < TLEepochs.size(); i++){
			A.at(j).at(i) = pow( TLEepochs.at(j), i ); // The coefficients will be stored in the order a0, a1, a2, a3,...
		};
	};

	// Apply the linear least-squares condition to the matrix of coefficients and right-hand side vector.
	std::vector< std::vector<double> > Atranspose = VectorOperations::transposeMatrix( &A );
	std::vector< std::vector<double> > LHS = VectorOperations::matrixMultiplyByMatrix( &Atranspose, &A ); // LHS coefficients of the Ax = b system become (A'A)x = A'b
	std::vector<double> RHS = VectorOperations::vectorMultiplyByMatrix( &SpecificOrbitalEnergies, &Atranspose ); // Also change the RHS.

	double luDecompositionNoSwaps;
	std::vector<int> luDecompositionIndices = std::vector<int>();
	
	std::vector< std::vector<double> > LU = EquationsSolving::luDecomposition(LHS, luDecompositionIndices, luDecompositionNoSwaps);
	EquationsSolving::luSubst(&LU, &luDecompositionIndices, &RHS); // Results i.e. solution of (A'A)x = A'b will be stored in the RHS vector - these are the coefficients of the least-squares line.

	// Compute the residuals, i.e. the differences in orbital energies between the regression curve and actual values.
	std::vector<double> SpecificOrbitalEnergiesInterpolated(SpecificOrbitalEnergies.size(), 0.0);
	std::vector<double> SpecificOrbitalEnergyResiduals = SpecificOrbitalEnergies; // Start with the actual values, substract the interpolated value to get a residual.

	for(std::vector<int>::size_type i=0; i<SpecificOrbitalEnergies.size(); i++){ // Compute the residual in orbital energy for every TLE and store these to get statistics later.
		for(std::vector<int>::size_type j=0; j<RHS.size(); j++){ // Gather a contribution from all the terms of the regression polynomial.
			SpecificOrbitalEnergiesInterpolated.at(i) = SpecificOrbitalEnergiesInterpolated.at(i) - std::pow(TLEepochs.at(i), (double) j )*RHS.at(j);
			SpecificOrbitalEnergyResiduals.at(i) = SpecificOrbitalEnergyResiduals.at(i) - std::pow(TLEepochs.at(i), (double) j )*RHS.at(j);
		}
	}

	// Get the statistics of the residuals on which the rejection threshold is based.
	double meanOrbitalEnergyResiduals = VectorOperations::mean( &SpecificOrbitalEnergyResiduals );
	double STDEVofOrbitalEnergyResiduals = VectorOperations::standardDeviation( &SpecificOrbitalEnergyResiduals );

	// Filter out the outliers here - reject TLEs that have specific orbital energies that fall outside a number of standard deviations from the average.
	double lowerThreshold = meanOrbitalEnergyResiduals - rejectionThreshold*STDEVofOrbitalEnergyResiduals;
	double upperThreshold = meanOrbitalEnergyResiduals + rejectionThreshold*STDEVofOrbitalEnergyResiduals;
	for(std::vector<int>::size_type i=0; i<SpecificOrbitalEnergyResiduals.size(); i++){
		// Sometimes when there are many TLEs the LU decomposition will hit the double precision limit and hence the interpolated orbital energies and so also the residuals will be NaNs. Don't reject then.
		#ifdef _MSC_VER // MSVS has a different _isnan function from <float.h>, Linux compilers use std::isnan from <cmath>
			if( (SpecificOrbitalEnergyResiduals.at(i) <= upperThreshold) && (SpecificOrbitalEnergyResiduals.at(i) >= lowerThreshold) ){ // This TLE's orbital energy falls inside the bounds.
				CovarianceTLEs.push_back( TEMPCovarainceTLEs.at(i) );
			}else if( _isnan(SpecificOrbitalEnergyResiduals.at(i)) ){ // Have many TLEs in this set - cannot assess whether any of them falls outside the bounds as they are NaNs. It's safer to keep all of them.
				CovarianceTLEs.push_back( TEMPCovarainceTLEs.at(i) );
			}else{
				std::cerr<<"A TLE to estimate the covariance for "<<NORAD_ID<<" falls outside "<<rejectionThreshold<<" standard deviations in specific orbital energy, it will be rejected from the sample"<<std::endl;
			}
		#else
			if( (SpecificOrbitalEnergyResiduals.at(i) <= upperThreshold) && (SpecificOrbitalEnergyResiduals.at(i) >= lowerThreshold) ){ // This TLE's orbital energy falls inside the bounds.
				CovarianceTLEs.push_back( TEMPCovarainceTLEs.at(i) );
			}else if( std::isnan(SpecificOrbitalEnergyResiduals.at(i)) ){ // Have many TLEs in this set - cannot assess whether any of them falls outside the bounds as they are NaNs. It's safer to keep all of them.
				CovarianceTLEs.push_back( TEMPCovarainceTLEs.at(i) );
			}else{
				std::cerr<<"A TLE to estimate the covariance for "<<NORAD_ID<<" falls outside "<<rejectionThreshold<<" standard deviations in specific orbital energy, it will be rejected from the sample"<<std::endl;
			}
		#endif
	}

	/* 
	 * Handle special cases when too few TLEs are available for a given object.
	 * Need to assume some covariance in this case. Base this on old ESOC estimates
	 * of TLE covariance in various orbital regimes. Set covariance based on 
	 * Flohrer, T. Assessment and categorisation of TLE orbit errors for the US SSN catalogue, 1 Jan 2008 snapshot.
	 */

	if((int)CovarianceTLEs.size() < MIN_TLES){ // Seal with cases that have too few TLEs so cannot estimate the covariance for them.
		std::string msg = "Too few TLEs for current object to estimate the covariance from."; // Throw an exception for now, don't want to obscure the effects of growing covarinace by having some objects that covariance of which doesn't grow over time.
		throw FatalError(msg,__FILE__,__LINE__);

// Don't assume fixed covariance in the AMOS 2014 version of the code - want to see the growth of covariance over time and its impact.
/*		std::cerr<<"Don't have any TLEs to estimate the covariance from for "<<NORAD_ID<<". Assuming standard deviations based on orbital regime given by Flohrer and 0 covariance matrix's cross terms. Covariance will be assumed time-invariant."<<std::endl;
		
		std::vector< std::vector<double> > tempMatrix = std::vector< std::vector<double> >(6, std::vector<double>(6, 0.0) ); // Iniitlaise a covariance matrix to be set.
		//UNDONE: get velocity covariance to realistically seed representative TLEs to propagate the covariance.
		if( Eccentricity <= 0.1 ){ // Assign covaraince values depending on the orbital regime (eccentricity, inclination, and perigee altitude). Flohrer gives standard deviations, need to square them to get covariance matrix tersm.
			if( Inclination<30.0 ){
				if( AltitudeOfPerigee<800.0 ){
					tempMatrix.at(0).at(0) = 0.067*0.067;
					tempMatrix.at(1).at(1) = 0.118*0.118;
					tempMatrix.at(2).at(2) = 0.075*0.075;
				}else if( (AltitudeOfPerigee>=800.0) && (AltitudeOfPerigee<=25000.0) ){
					tempMatrix.at(0).at(0) = 0.191*0.191;
					tempMatrix.at(1).at(1) = 0.256*0.256;
					tempMatrix.at(2).at(2) = 0.203*0.203;
				}else{ // AltitudeOfPerigee > 25000 km.
					tempMatrix.at(0).at(0) = 0.357*0.357;
					tempMatrix.at(1).at(1) = 0.432*0.432;
					tempMatrix.at(2).at(2) = 0.083*0.083;
				}
			}else if( (Inclination>=30.0) && (Inclination<=60.0) ){
				if( AltitudeOfPerigee<800.0 ){
					tempMatrix.at(0).at(0) = 0.107*0.107;
					tempMatrix.at(1).at(1) = 0.308*0.308;
					tempMatrix.at(2).at(2) = 0.169*0.169;
				}else if( (AltitudeOfPerigee>=800.0) && (AltitudeOfPerigee<=25000.0) ){
					tempMatrix.at(0).at(0) = 0.071*0.071;
					tempMatrix.at(1).at(1) = 0.228*0.228;
					tempMatrix.at(2).at(2) = 0.095*0.095;
				}else{ // AltitudeOfPerigee > 25000 km.
					std::cerr<<"Encountered and object with no TLE to compute the covariance from and in an orbital regime: e<=0.1, 30<=i<=60, Hp>25000. Assuming 0 standard deviations."<<std::endl;
					tempMatrix.at(0).at(0) = 0.0; // Flohrer gives no information for this regime as objects are unlikely. Assume 0 covariance.
					tempMatrix.at(1).at(1) = 0.0;
					tempMatrix.at(2).at(2) = 0.0;
				}
			}else{ // Inclination > 60 degrees.
				if( AltitudeOfPerigee<800.0 ){
					tempMatrix.at(0).at(0) = 0.115*0.115;
					tempMatrix.at(1).at(1) = 0.517*0.517;
					tempMatrix.at(2).at(2) = 0.137*0.137;
				}else if( (AltitudeOfPerigee>=800.0) && (AltitudeOfPerigee<=25000.0) ){
					tempMatrix.at(0).at(0) = 0.091*0.091;
					tempMatrix.at(1).at(1) = 0.428*0.428;
					tempMatrix.at(2).at(2) = 0.114*0.114;
				}else{ // AltitudeOfPerigee > 25000 km.
					std::cerr<<"Encountered and object with no TLE to compute the covariance from and in an orbital regime: e<=0.1, 60<i, Hp>25000. Assuming 0 standard deviations."<<std::endl;
					tempMatrix.at(0).at(0) = 0.0; // Flohrer gives no information for this regime as objects are unlikely. Assume 0 covariance.
					tempMatrix.at(1).at(1) = 0.0;
					tempMatrix.at(2).at(2) = 0.0;
				}
			}
		}else{ // Eccentricity>0.1
			if( Inclination<30.0 ){
				if( AltitudeOfPerigee<800.0 ){
					tempMatrix.at(0).at(0) = 0.2252*0.2252;
					tempMatrix.at(1).at(1) = 0.427*0.427;
					tempMatrix.at(2).at(2) = 0.1421*0.1421;
				}else if( (AltitudeOfPerigee>=800.0) && (AltitudeOfPerigee<=25000.0) ){
					tempMatrix.at(0).at(0) = 0.1748*0.1748;
					tempMatrix.at(1).at(1) = 0.3119*0.3119;
					tempMatrix.at(2).at(2) = 0.971*0.971;
				}else{ // AltitudeOfPerigee > 25000 km.
					tempMatrix.at(0).at(0) = 0.402*0.402;
					tempMatrix.at(1).at(1) = 0.418*0.418;
					tempMatrix.at(2).at(2) = 0.083*0.083;
				}
			}else if( (Inclination>=30.0) && (Inclination<=60.0) ){
				if( AltitudeOfPerigee<800.0 ){
					tempMatrix.at(0).at(0) = 0.629*0.629;
					tempMatrix.at(1).at(1) = 0.909*0.909;
					tempMatrix.at(2).at(2) = 0.2057*0.2057;
				}else if( (AltitudeOfPerigee>=800.0) && (AltitudeOfPerigee<=25000.0) ){
					tempMatrix.at(0).at(0) = 0.1832*0.1832;
					tempMatrix.at(1).at(1) = 0.1878*0.1878;
					tempMatrix.at(2).at(2) = 0.1454*0.1454;
				}else{ // AltitudeOfPerigee > 25000 km.
					tempMatrix.at(0).at(0) = 0.4712*0.4712;
					tempMatrix.at(1).at(1) = 0.6223*0.6223;
					tempMatrix.at(2).at(2) = 0.1208*0.1208;
				}
			}else{ // Inclination > 60 degrees.
				if( AltitudeOfPerigee<800.0 ){
					tempMatrix.at(0).at(0) = 0.494*0.494;
					tempMatrix.at(1).at(1) = 0.814*0.814;
					tempMatrix.at(2).at(2) = 0.1337*0.1337;
				}else if( (AltitudeOfPerigee>=800.0) && (AltitudeOfPerigee<=25000.0) ){
					tempMatrix.at(0).at(0) = 0.529*0.529;
					tempMatrix.at(1).at(1) = 0.817*0.817;
					tempMatrix.at(2).at(2) = 0.1570*0.1570;
				}else{ // AltitudeOfPerigee > 25000 km.
					std::cerr<<"Encountered and object with no TLE to compute the covariance from and in an orbital regime: 0.1<e, 60<i, Hp>25000. Assuming 0 standard deviations."<<std::endl;
					tempMatrix.at(0).at(0) = 0.0; // Flohrer gives no information for this regime as objects are unlikely. Assume 0 covariance.
					tempMatrix.at(1).at(1) = 0.0;
					tempMatrix.at(2).at(2) = 0.0;
				}
			}
		}
		SetCovarianceMatrixRTC( tempMatrix ); // Set the attributes and spawn CovarianceTLEs for compatibility with other methods.
		*/
	}//END if(CovarianceTLEs.size()<MIN_TLES)
};

std::vector< std::vector<double> > SpaceObject::InertialToRTC(std::vector<double>* positionPtr, std::vector<double>* velocityPtr){
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

void SpaceObject::ComputeCovarianceOSW( double epochJDAY ){
	/* Compute full 6x6 covariance matrix (position and velocity) of Two-Line Elements at the desired epoch and save it in FullCovarianceMatrixRTC.
	Also save the 3x3 position covariance in PositionCovarianceMatrixRTC to use it later in true collision probability computations.
	This is done using the approach described by V. Osweiler 2006 i.e. by propagating all the TLEs to the desired epoch,
	computing residuals w.r.t. to the state vector corresponding to this SpaceObject at that epoch and estimating the covariance based on those residuals.
	Use TLEs saved in CovarianceTLEs to estimate the covariance (generated by @see ReadThreeLineELements).
	@param epochJDAY - Julian Day at which to compute the covariance.
	@return FullCovarianceMatrixRTC - 6x6 position (km) and velocity (km/sec) covariance matrix with variances on the diagonal and covariance coefficients off-diagonal.
		Positions are in rows 1 to 3, velocities 4 to 6. It is given in the velocity - normal - co-normal (cross-track) reference frame computed for this SpaceObject..
	@return PositionCovarianceMatrixRTC - 3x3 position (km) covariance matrix with variances on the diagonal and covariance coefficients off-diagonal.
		Positions are in rows 1 to 3. It is given in the velocity - normal - co-normal (cross-track) reference frame computed for this SpaceObject.
	*/
	//UNDONE: when generating TLEs from a covariance matrix is dealt with then won't need a safety check here as covarianceTLEs will never be empty and will be able to propagate the covariance to any epoch.
	if(CovarianceTLEs.size()>=MIN_TLES){ // Can propagate the covariance.
		/* Most-recent state and derived quantities. */
		std::vector<double> mostRecentPos = std::vector<double>(3, 0.0);
		std::vector<double> mostRecentVelo = std::vector<double>(3, 0.0);
		PropagateJDAY(&mostRecentPos, &mostRecentVelo, epochJDAY, false); // Get the position of the primary TLE (this SpaceObject) at the epoch of interest - it defines the local reference frame etc.

		std::vector< std::vector<double> > inertialToRTCrotation = InertialToRTC(&mostRecentPos, &mostRecentVelo); // Rotation matrix from the inertial frame to the RTC frame established by the state of this SpaceObject at the epoch of interest.

		/* Process all the TLEs that define the covariance matrix. */
		std::vector< std::vector<double> >* residuals = new std::vector< std::vector<double> >();
		residuals->resize(6); // Each entry of residuals is a vector with cartesian position and velocity components.S

		std::vector<double> tempPosResidual = std::vector<double>(3, 0.0); // 3 Cartesian position and 3 velocity components.
		std::vector<double> tempVeloResidual = std::vector<double>(3, 0.0);

		double meanX = 0.0, meanY = 0.0, meanZ = 0.0, meanVx=0.0, meanVy=0.0, meanVz=0.0; // Mean values of the residuals.
		int noResiduals = 0; // Number of the residuals read.

		/* Compute the residuals. */
		double r[3]; double v[3]; // Create arrays that SGP4 expects to find states of all the covariance-defining TLEs at the epoch of interest.
		for(std::vector<int>::size_type i=0; i<CovarianceTLEs.size()-1; i++){
			double minutesSinceEpoch = (epochJDAY-TLEepochJDAY)*1440.0; // Minutes since the prime TLE's epoch.
			sgp4(gravityModel, CovarianceTLEs.at(i), minutesSinceEpoch, r,  v); 
			for(std::vector<int>::size_type j=0; j<tempPosResidual.size(); j++){
				tempPosResidual.at(j) = r[j]-mostRecentPos.at(j);
				tempVeloResidual.at(j) = v[j]-mostRecentVelo.at(j);
			}

			tempPosResidual = VectorOperations::vectorMultiplyByMatrix(&tempPosResidual, &inertialToRTCrotation); // Rotate the new residual to VNC established by the prime TLE.
			tempVeloResidual = VectorOperations::vectorMultiplyByMatrix(&tempVeloResidual, &inertialToRTCrotation); 

			meanX += tempPosResidual.at(0); // Record all the residuals found to compute the average at the end.
			meanY += tempPosResidual.at(1);
			meanZ += tempPosResidual.at(2);
			meanVx += tempVeloResidual.at(0);
			meanVy += tempVeloResidual.at(1);
			meanVz += tempVeloResidual.at(2);

			for(std::vector<int>::size_type j=0; j<tempPosResidual.size(); j++){
				residuals->at(j).push_back( tempPosResidual.at(j) ); // Record this residual.
				residuals->at(j+3).push_back( tempVeloResidual.at(j) );
			}
			noResiduals += 1;
		}

		meanX = meanX/noResiduals; meanY = meanY/noResiduals; meanZ = meanZ/noResiduals;
		meanVx = meanVx/noResiduals; meanVy = meanVy/noResiduals; meanVz = meanVz/noResiduals;

		std::vector< std::vector<double> > residualsTranspose = VectorOperations::transposeMatrix( residuals ); // Transpose the matrix.

		std::vector< std::vector<double> > posVeloCovariance = VectorOperations::matrixMultiplyByMatrix(residuals, &residualsTranspose);
		FullCovarianceMatrixRTC = VectorOperations::matrixMultiplyByScalar(&posVeloCovariance, 1.0/noResiduals); // Divide through by the number of residuals to get the completed covariance matrix.

		// Extract only the position components of the covariance matrix.
		PositionCovarianceMatrixRTC.at(0).at(0) = FullCovarianceMatrixRTC.at(0).at(0); PositionCovarianceMatrixRTC.at(0).at(1) = FullCovarianceMatrixRTC.at(0).at(1); PositionCovarianceMatrixRTC.at(0).at(2) = FullCovarianceMatrixRTC.at(0).at(2);
		PositionCovarianceMatrixRTC.at(1).at(0) = FullCovarianceMatrixRTC.at(1).at(0); PositionCovarianceMatrixRTC.at(1).at(1) = FullCovarianceMatrixRTC.at(1).at(1); PositionCovarianceMatrixRTC.at(1).at(2) = FullCovarianceMatrixRTC.at(1).at(2);
		PositionCovarianceMatrixRTC.at(2).at(0) = FullCovarianceMatrixRTC.at(2).at(0); PositionCovarianceMatrixRTC.at(2).at(1) = FullCovarianceMatrixRTC.at(2).at(1); PositionCovarianceMatrixRTC.at(2).at(2) = FullCovarianceMatrixRTC.at(2).at(2);
	}//END have some covarianceTLEs to propagate the covariance.
};

SpaceObject::~SpaceObject(void){
 /* Deconstructor, does nothing special.*/
};

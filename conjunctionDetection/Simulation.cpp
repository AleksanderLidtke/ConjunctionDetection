#include "Simulation.h"

Simulation::Simulation(void){
	/* Default constrctor, does nothing. */
};

Simulation::~Simulation(void){
	/* Deconstructor, does nothing. */
};

Simulation::Simulation(std::string title, int wgsEllipsoidType, double thresholdDistance, double coarseTimeStep, int noInterpolationPoints, double perigeeApogeeFilterPad){
	/* GRAVITY MODEL PROPERTIES. */
	wgsEllipsoidUsed = 84; // Earth ellipsoid used.
	gravityModel = gravconsttype();
	getgravconst( gravityModel, tumin, GRAVITY_CONSTANT, MEAN_EARTH_RADIUS, xke, j2, j3, j4, j3oj2 );

	/* ANALYSIS SETTINGS. */
	TITLE = title;
	CONJUNCTION_THRESHOLD_DISTANCE = thresholdDistance; // km, all cases when objects get below this separation distance will be found.
	CONJUNCTION_THRESHOLD_DISTANCE_SQUARED = thresholdDistance*thresholdDistance; // km^2, square of the threshold distance.
    COARSE_TIME_STEP_S = coarseTimeStep; // s, time step in which a number of interpolation points is created to generate an ephemeris table from which to look for conjunctions.
    COARSE_TIME_STEP = COARSE_TIME_STEP_S * SECOND_IN_JDAYS; // Julian days, time step in which a number of interpolation points is created to generate an ephemeris table from which to look for conjunctions.
	NO_INTERPOLATION_POINTS = noInterpolationPoints; // Number of pooints used to interpolate trajectories.
    PERIGEE_APOGEE_PAD = perigeeApogeeFilterPad; // km, a buffer to be used in perigee/apogee pre-filtering.
	PERIGEE_APOGEE_ORBITAL_REGIME_PAD = 1000.0; // km, a buffer that is used to filter out object that are in different orbital regimes and hence will never conjunct.
	EnforceDirectNumericalIntegration = false; // Don't use direct integration by default.

	/* PRE-FILTERS' SETTINGS. */
	SURFACE_GRAVITY_ACCELERATION = GRAVITY_CONSTANT/(MEAN_EARTH_RADIUS*MEAN_EARTH_RADIUS); // km/s^2
	ESCAPE_VELOCITY = sqrt( 2.0*GRAVITY_CONSTANT/MEAN_EARTH_RADIUS ); // km/s
	THRESHOLD_RADIUS = CONJUNCTION_THRESHOLD_DISTANCE + ESCAPE_VELOCITY*COARSE_TIME_STEP_S;
	
    THRESHOLD_RADIUS_SQUARED = 2.0*THRESHOLD_RADIUS*THRESHOLD_RADIUS; // Add a safety factor of 2.0 in order not to miss conjunctions.
	double tempDistanceThreshold[3] = {THRESHOLD_RADIUS, THRESHOLD_RADIUS, THRESHOLD_RADIUS};
	DISTANCE_THRESHOLD = std::vector<double>( tempDistanceThreshold, tempDistanceThreshold+sizeof(tempDistanceThreshold)/sizeof(tempDistanceThreshold[0]) );
    ACCELERATION_RADIUS = CONJUNCTION_THRESHOLD_DISTANCE + SURFACE_GRAVITY_ACCELERATION*COARSE_TIME_STEP_S*COARSE_TIME_STEP_S;
    ACCELERATION_RADIUS_SQUARED = ACCELERATION_RADIUS*ACCELERATION_RADIUS;

	/* OBJECTS, OBJECT PAIRS AND CONJUNCTIONS. */
	objectsPtr = new std::map<std::string, SpaceObject>;
	conjunctionsFoundPtr = new std::vector<Conjunction>(); // Stores information about all the conjunctions.
	conjunctionsFoundPtr->reserve(1000000); // Reserve space for conjunctions to be stored.

	/* INTERPOLATION PARAMETERS. */
	polynomialOrder = 2*NO_INTERPOLATION_POINTS - 1;
	interpolationArguments = std::vector<double>(NO_INTERPOLATION_POINTS, 0.0);
    VectorOperations::linspace(&interpolationArguments, 0.0, 1.0, &NO_INTERPOLATION_POINTS); // Interpolate in interval [0,1]

	/* Build the matrix of coefficients to solve the system of equations (also include rows corresponding to derivatives in the lower half of the matrix):
    [1, argumnet, argument^2, ..., argument^(N-1)] * [unknown polynomial coefficients] = [known data point values] */
    A = std::vector< std::vector<double> >(polynomialOrder+1, std::vector<double>(polynomialOrder+1,0.0) );
    
	for(std::vector<int>::size_type j = 0; j < (std::vector<int>::size_type) int((polynomialOrder+1)/2.0); j++){ // First top half of the rows corresponds to the function values.
        for(std::vector<int>::size_type i=0; i < (std::vector<int>::size_type) polynomialOrder+1; i++){
            A.at(j).at(i) = pow( interpolationArguments.at(j), i );
		};
	};
    for(std::vector<int>::size_type j = 0; j < (std::vector<int>::size_type) int((polynomialOrder+1)/2.0); j++){ // The rest (bottom half) correpsonds to derivative values.
        for(std::vector<int>::size_type i=1; i < (std::vector<int>::size_type) polynomialOrder+1; i++){
            A.at(int((polynomialOrder+1)/2)+j).at(i) = i*pow(interpolationArguments.at(j), i-1);
		};
	};
	rhsX = std::vector<double>(2*NO_INTERPOLATION_POINTS,0.0); // Right hand sides of the sets of equations to be solved to find the interpolation coefficients..
	rhsY = std::vector<double>(2*NO_INTERPOLATION_POINTS,0.0);
	rhsZ = std::vector<double>(2*NO_INTERPOLATION_POINTS,0.0);

	luDecompositionIndices = std::vector<int>();
	LU = EquationsSolving::luDecomposition(A, luDecompositionIndices, luDecompositionNoSwaps);
	
	/* FINAL CONSTRUCTOR COMMANDS. */
	char buff[25]; time_t now = time (0);
    #ifdef _MSC_VER 
		struct tm timeinfo;
		localtime_s(&timeinfo, &now);
		strftime(buff, 25, "%Y-%m-%d %H:%M:%S.000", &timeinfo);
	#else
		strftime(buff, 25, "%Y-%m-%d %H:%M:%S.000", localtime (&now));
	#endif
	std::cout<<buff<<": Created "<<TITLE<<" simulation."<<std::endl;
};

void Simulation::setDirectNumericalIntegration( bool flag ){
	/* Set the EnforceDirectNumericalIntegration attribute that decides whether direct numerical integration of the relative position's probability
	density function (PDF) should be used at all times regardless of conjunction geomerty, object sizes etc.
	*/
	EnforceDirectNumericalIntegration = flag; // That's all that needs to happen.
};

void Simulation::performAnalysis(int startYear, int startMonth, int startDay, int startHour, int startMinute, double startSecond, int stopYear, int stopMonth, int stopDay, int stopHour, int stopMinute, double stopSecond, int mode, std::vector<std::string> primariesSSCs){
	/* Perform conjunction detection and analysis between the objects that have been created in @createObjectsFromTLEfile based on the pairs that will 
	be generated either in an one-on-all or all-on-all mode. Can also choose which objects to treat as primaries for the one-on-all mode.
	All the conjunctions' data will be stored in conjunctionsFoundPtr. A convenience wrapper for the same function using Julian days.
	@param start - UTCG epoch of the beginning of the analysis.
	@param stop - UTCG epoch of the end of the analysis.
	@param mode - 1 for one-on-all or multiple-on-all (several primaries), 2 for all-on-all.
	@param primariesSSCc - NORAD IDs of objects to be treated as primaries for which to find conjunctions against the catalogue that has been read. Ignored if mode==2 but at least 
		1 must be provided when mode==1.
	@throws FatalError - if mode is neither 1 nor 2. Also when primariesSSCs are empty when mode==1.
	*/
	double ANALYSIS_INTERVAL_START; jday(startYear,startMonth,startDay,startHour,startMinute,startSecond, ANALYSIS_INTERVAL_START);
    double ANALYSIS_INTERVAL_STOP; jday(stopYear,stopMonth,stopDay,stopHour,stopMinute,stopSecond, ANALYSIS_INTERVAL_STOP);
	
	// Call a method that accepts Julian day arguments.
	performAnalysis(ANALYSIS_INTERVAL_START, ANALYSIS_INTERVAL_STOP, mode, primariesSSCs);

};

void Simulation::performAnalysis(double ANALYSIS_INTERVAL_START, double ANALYSIS_INTERVAL_STOP, int mode, std::vector<std::string>primariesSSCs){
	/* Perform conjunction detection and analysis between the objects that have been created in @createObjectsFromTLEfile based on the pairs that will 
	be generated either in an one-on-all or all-on-all mode. Can also choose which objects to treat as primaries for the one-on-all mode.
	All the conjunctions' data will be stored in conjunctionsFoundPtr.
	@param ANALYSIS_INTERVAL_START - UTCG epoch of the beginning of the analysis in Julian days.
	@param ANALYSIS_INTERVAL_STOP - UTCG epoch of the end of the analysis in Julian days.
	@param mode - 1 for one-on-all or multiple-on-all (several primaries), 2 for all-on-all.
	@param primariesSSCc - NORAD IDs of objects to be treated as primaries for which to find conjunctions against the catalogue that has been read. Ignored if mode==2 but at least 
		1 must be provided when mode==1.
	@throws FatalError - if mode is neither 1 nor 2. Also when primariesSSCs are empty when mode==1.
	*/
	
	/* DISPLAY INITIAL INFORMATION. */
	TIME_STARTED = time(0);
	char buff[25];
	#ifdef _MSC_VER 
		struct tm timeinfo;
		localtime_s(&timeinfo, &TIME_STARTED);
		strftime(buff, 25, "%Y-%m-%d %H:%M:%S.000", &timeinfo);
	#else
		strftime(buff, 25, "%Y-%m-%d %H:%M:%S.000", localtime (&TIME_STARTED));
	#endif
	
	// Get information about the analysis period for display purposes.
	int yearStart, monStart, dayStart, hourStart, minuteStart;
	double secStart;
	invjday(ANALYSIS_INTERVAL_START, yearStart, monStart, dayStart, hourStart, minuteStart, secStart);
	int yearStop, monStop, dayStop, hourStop, minuteStop;
	double secStop;
	invjday(ANALYSIS_INTERVAL_STOP, yearStop, monStop, dayStop, hourStop, minuteStop, secStop);

	std::cout<<buff<<" Starting the analysis in mode "<<mode<<std::endl
		<<"Analysis interval start: "<<dayStart<<"/"<<monStart<<"/"<<yearStart<<" "<<hourStart<<":"<<minuteStart<<":"<<secStart<<std::endl
		<<"Analysis interval stop: "<<dayStop<<"/"<<monStop<<"/"<<yearStop<<" "<<hourStop<<":"<<minuteStop<<":"<<secStop<<std::endl;

	/* LOOK FOR CONJUNCTIONS */
	int counter = 0; // Count number of conjunctions that have been analysed.
	double totalPairs; // Number of object pairs to be analysed.
	if(mode==1){
		if(primariesSSCs.empty()){
			std::string msg1 = "Unspecified SSCs of the primaries in the one-on-all analysis mode that has been selected, which is meaningless.";
			throw FatalError(msg1,__FILE__,__LINE__);
		}else{
			totalPairs =(int)primariesSSCs.size()*( (int)objectsPtr->size() -1); // Number of object pairs to be analysed = number of primaries to analyse X number of pairs each.
			for(size_t i=0; i<primariesSSCs.size(); i++){
				std::cout<<"Primary SSC: "<<primariesSSCs.at(i)<<std::endl;
				for(std::map<std::string, SpaceObject>::iterator it = objectsPtr->begin(); it != objectsPtr->end(); ++it){
					if(primariesSSCs.at(i) != it->first){ // Don't include the primary in the secondaries.
						if(i>=1){
							if(primariesSSCs.at(i-1)!= it->first){ // Don't count the same pair twice (two primaries).
								std::pair<SpaceObject, SpaceObject> tempPair = std::pair<SpaceObject, SpaceObject>( objectsPtr->at(primariesSSCs.at(i)), it->second );
								findConjunctionsBetweenObjectPair( &tempPair, conjunctionsFoundPtr, &ANALYSIS_INTERVAL_START, &ANALYSIS_INTERVAL_STOP);
								counter += 1;
							};
						}else{
							std::pair<SpaceObject, SpaceObject> tempPair = std::pair<SpaceObject, SpaceObject>( objectsPtr->at(primariesSSCs.at(i)), it->second );
							findConjunctionsBetweenObjectPair( &tempPair, conjunctionsFoundPtr, &ANALYSIS_INTERVAL_START, &ANALYSIS_INTERVAL_STOP);
							counter += 1;
						};
					};

					if( fmod(counter,100.0) == 0 ){
						char tempBuff[25]; time_t now = time(0);
						#ifdef _MSC_VER 
							struct tm timeinfo;
							localtime_s(&timeinfo, &now);
							strftime(tempBuff, 25, "%Y-%m-%d %H:%M:%S.000", &timeinfo);
						#else
							strftime(tempBuff, 25, "%Y-%m-%d %H:%M:%S.000", localtime (&now));
						#endif
						std::cout<<tempBuff<<" "<<TITLE<<": Done "<<counter/totalPairs*100<<"% of object pairs out of total of "<<totalPairs<<" pairs to be analysed. Found "<<(int)conjunctionsFoundPtr->size()<<" conjunctions so far."<<std::endl;
					};
				};
			};
		};
	}else if(mode==2){
		totalPairs = 0.0;
		for(double n=1; n<(double)objectsPtr->size(); n++){
			totalPairs += n;
		};

		double counter = 0.0; // Number of objects analysed.
		for(std::map<std::string, SpaceObject>::iterator outerIt = objectsPtr->begin(); outerIt != objectsPtr->end(); ++outerIt){
			std::map<std::string, SpaceObject>::iterator innerIt = outerIt;
			for(std::advance(innerIt, 1); innerIt != objectsPtr->end(); ++innerIt){
				if(outerIt->first != innerIt->first){ // Don't include duplicates.
					std::pair<SpaceObject, SpaceObject> tempPair = std::pair<SpaceObject, SpaceObject>( outerIt->second, innerIt->second );
					findConjunctionsBetweenObjectPair( &tempPair, conjunctionsFoundPtr, &ANALYSIS_INTERVAL_START, &ANALYSIS_INTERVAL_STOP);
					counter += 1;
				};
			};
			if( fmod(counter,100.0) == 0 ){
				char tempBuff[25]; time_t now = time(0);
				#ifdef _MSC_VER 
					struct tm timeinfo;
					localtime_s(&timeinfo, &now);
					strftime(tempBuff, 25, "%Y-%m-%d %H:%M:%S.000", &timeinfo);
				#else
					strftime(tempBuff, 25, "%Y-%m-%d %H:%M:%S.000", localtime (&now));
				#endif
				std::cout<<tempBuff<<" "<<TITLE<<": Done "<<counter/totalPairs*100<<"% of object pairs out of total of "<<totalPairs<<" pairs to be analysed. Found "<<(int)conjunctionsFoundPtr->size()<<" conjunctions so far."<<std::endl;
			};
		};
	}else{
		std::string msg2 = "Supplied an unrecognised operations mode when creating pairs of objects for analysis.";
		throw FatalError(msg2,__FILE__,__LINE__);
	}
	TIME_FINISHED = time(0);
	char endBuff[25];
    #ifdef _MSC_VER 
		struct tm endTimeInfo;
		localtime_s(&endTimeInfo, &TIME_FINISHED);
		strftime(endBuff, 25, "%Y-%m-%d %H:%M:%S.000", &endTimeInfo);
	#else
		strftime(endBuff, 25, "%Y-%m-%d %H:%M:%S.000", localtime (&TIME_FINISHED));
	#endif
	std::cout<<endBuff<<": Analysed "<<totalPairs<<" pairs of objects in "<<difftime(TIME_FINISHED,TIME_STARTED)<<" seconds."<<std::endl;
};

void Simulation::performAnalysis(std::string primarySSC, std::string primaryEphemerisFileName, double primaryRadius){
	/* Perform conjunction detection and analysis between the objects that have been created in @createObjectsFromTLEfile based on the pairs that will 
	be generated either in an one-on-all or all-on-all mode. Can also choose which objects to treat as primaries for the one-on-all mode.
	All the conjunctions' data will be stored in conjunctionsFoundPtr.
	@param primarySSC - NORAD ID of the primary object for which to find conjunctions against the catalogue that has been read. Used to avoid finding conjunction between the same object.
	@param primaryEphemerisFileName - name of the file from which to read the ephemeris table for the primary.
	@param primaryRadius - radius of the primary object to be used in collision probabiity calculations.
	@throws FatalError - if mode is neither 1 nor 2. Also when primariesSSCs are empty when mode==1.
	*/

	/* DISPLAY INITIAL INFORMATION. */
	TIME_STARTED = time(0);
	char buff[25];
	#ifdef _MSC_VER 
		struct tm timeinfo;
		localtime_s(&timeinfo, &TIME_STARTED);
		strftime(buff, 25, "%Y-%m-%d %H:%M:%S.000", &timeinfo);
	#else
		strftime(buff, 25, "%Y-%m-%d %H:%M:%S.000", localtime (&TIME_STARTED));
	#endif
	
	/* Read the ephemeris of the primary and get information about the analysis period for display purposes. */
	ExternalSpaceObject primary = ExternalSpaceObject( primarySSC, primaryEphemerisFileName ); // Load the ephemeris for the primary.
	primary.SetHardBodyRadius( primaryRadius );

	int yearStart, monStart, dayStart, hourStart, minuteStart;
	double secStart;
	invjday(primary.EphemerisEpochs.at(0), yearStart, monStart, dayStart, hourStart, minuteStart, secStart);
	int yearStop, monStop, dayStop, hourStop, minuteStop;
	double secStop;
	invjday(primary.EphemerisEpochs.back(), yearStop, monStop, dayStop, hourStop, minuteStop, secStop);
	
	COARSE_TIME_STEP = primary.EphemerisEpochs.at(1) - primary.EphemerisEpochs.at(0); // Override this - the analysis intervals are coupled to ephemeris epochs for simplicity.
	COARSE_TIME_STEP_S = COARSE_TIME_STEP * JDAY_IN_SECONDS;

	std::cout<<buff<<" Starting the analysis for screening object no. "<<primarySSC<<" for conjunctions. Ephemeris is supplied from "<<primaryEphemerisFileName<<std::endl
		<<"Analysis interval start: "<<dayStart<<"/"<<monStart<<"/"<<yearStart<<" "<<hourStart<<":"<<minuteStart<<":"<<secStart<<std::endl
		<<"Analysis interval stop: "<<dayStop<<"/"<<monStop<<"/"<<yearStop<<" "<<hourStop<<":"<<minuteStop<<":"<<secStop<<std::endl;

	/* LOOK FOR CONJUNCTIONS */
	int counter = 0; // Count number of conjunctions that have been analysed.
	double totalPairs = ( (int)objectsPtr->size() -1); // Number of object pairs to be analysed = number of primaries to analyse X number of pairs each.
	for(std::map<std::string, SpaceObject>::iterator it = objectsPtr->begin(); it != objectsPtr->end(); ++it){
		if(primarySSC != it->first){ // Don't include the primary in the secondaries. primarySSC may be an empty string when it is not believed to be in the objectsPtr, but then this check will be passed as well.
			std::pair<ExternalSpaceObject, SpaceObject> tempPair = std::pair<ExternalSpaceObject, SpaceObject>( primary, it->second );
			findConjunctionsBetweenObjectPair( &tempPair, conjunctionsFoundPtr);
			counter += 1;
		};
	
		if( fmod(counter,100.0) == 0 ){
			char tempBuff[25]; time_t now = time(0);
			#ifdef _MSC_VER 
				struct tm timeinfo;
				localtime_s(&timeinfo, &now);
				strftime(tempBuff, 25, "%Y-%m-%d %H:%M:%S.000", &timeinfo);
			#else
				strftime(tempBuff, 25, "%Y-%m-%d %H:%M:%S.000", localtime (&now));
			#endif
			std::cout<<tempBuff<<" "<<TITLE<<": Done "<<counter/totalPairs*100<<"% of object pairs out of total of "<<totalPairs<<" pairs to be analysed. Found "<<(int)conjunctionsFoundPtr->size()<<" conjunctions so far."<<std::endl;
		};
	};

	TIME_FINISHED = time(0);
	char endBuff[25];
    #ifdef _MSC_VER 
		struct tm endTimeInfo;
		localtime_s(&endTimeInfo, &TIME_FINISHED);
		strftime(endBuff, 25, "%Y-%m-%d %H:%M:%S.000", &endTimeInfo);
	#else
		strftime(endBuff, 25, "%Y-%m-%d %H:%M:%S.000", localtime (&TIME_FINISHED));
	#endif
	std::cout<<endBuff<<": Analysed "<<counter<<" pairs of objects in "<<difftime(TIME_FINISHED,TIME_STARTED)<<" seconds."<<std::endl;
}

void Simulation::findConjunctionsBetweenObjectPair(std::pair<SpaceObject, SpaceObject>* objectPairPtr, std::vector<Conjunction>* conjunctionsPtr, double* analysisStart, double* analysisStop){
	/* Look for all the conjunctions that are predicted to happen between the given pair of objects in the time interval specified.
	Assume that objects that reach altitude of 100.0 km have re-entered the atmosphere.
	Use perigee/apogee sieve (computed at the end of the interpolation interval with a pad of PERIGEE_APOGEE_PAD km) with the caveat 
	that if the difference in altitudes	is 1000.0 km or more the pair will be ignored for speed.
	Also use the "smart sieve" as described in ROdriguez et. all 2002, namely:
		- XYZ
		- r^2
		- minimum distance
		- fine r^2
	First try to discard the pair at every coarse time step based on the prefilters. If all of them are passed find the conjunction
	based purely on the range measure i.e. find the epoch when the relative distance is minimum and below CONJUNCTION_THRESHOLD_DISTANCE.
	Once a conjunction is found its detailed information is refined by propagating the objects in a shorter interval around it.
	@param objectPairPtr - pointer to a pair of objects <object1, object2> between which the conjunctions are to be found.
	@param conujnctionsPtr - vector where all the found conjunctions will be recorded.
	@param analysisStart, analysisStop - Julian days marking the beginning and the end of the analysis intetrval, respectively.
	@param noInterpolationPoints - number of points to be used for interpolation of trajectories within a singe coarse analysis time step.
	@param COARSE_TIME_STEP_S - the duration of the analysis' time step (can, in principle, be different than the interpolation interval if more than 2 interpolation nodes are used).
	@param COARSE_TIME_STEP - same as COARSE_TIME_STEP_S but in Julian days.
	@param PERI_APOGEE_FILETR_PAD - buffer to be used when using the apogee/perigee prefilter, in km.
	@param THRESHOLD_RADIUS - the threshold radius as defined by Rodriguez 2002, in km.
	@param THRESHOLD_RADIUS_SQUARED - square of the threshold radius as defined by Rodriguez 2002, in km.
	@param ACCELERATION_RADIUS - the acceleration radius as defined by Rodriquez 2002, in km.
	@param ACCELERATION_RADIUS_SQUARED - square of the acceleration radius as defined by Rodriquez 2002, in km.
	*/
	double CURRENT_JDAY = *analysisStart;

	/* Positions, velocities and epochs of objects in Cartesian space used for interpolation. */
	std::vector< std::vector<double> > posesInterp1 = std::vector< std::vector<double> >(NO_INTERPOLATION_POINTS, std::vector<double>(3,0.0) ); // Each row has 3 Cartesian components.
	std::vector< std::vector<double> > velosInterp1 = std::vector< std::vector<double> >(NO_INTERPOLATION_POINTS, std::vector<double>(3,0.0) );
	std::vector< std::vector<double> > posesInterp2 = std::vector< std::vector<double> >(NO_INTERPOLATION_POINTS, std::vector<double>(3,0.0) );
	std::vector< std::vector<double> > velosInterp2 = std::vector< std::vector<double> >(NO_INTERPOLATION_POINTS, std::vector<double>(3,0.0) );
	std::vector<double> currentJDAYs = std::vector<double>(NO_INTERPOLATION_POINTS, 0.0); // Current interpolation points' epochs.

	/* Interpolating polynomials' coefficients for both objects - row indices correspond to dimensions, column indices to the particular coefficients. */
	std::vector< std::vector<double> >* interpCoeffsPointer1 = new std::vector< std::vector<double> >(3, std::vector<double>(2*NO_INTERPOLATION_POINTS, 0.0) );
	std::vector< std::vector<double> >* interpCoeffsPointer2 = new std::vector< std::vector<double> >(3, std::vector<double>(2*NO_INTERPOLATION_POINTS, 0.0) );
	std::vector< std::vector<double> >* interpCoeffsPointerRelativePointer = new std::vector< std::vector<double> >(3, std::vector<double>(2*NO_INTERPOLATION_POINTS, 0.0) );

	/* Positions and velocities of objects in Cartesian space used for pre-filtering. */
	std::vector<double> posesPreFiltering1 = std::vector<double>(3,0.0); std::vector<double> velosPreFiltering1 = std::vector<double>(3,0.0);
	std::vector<double> posesPreFiltering2 = std::vector<double>(3,0.0); std::vector<double> velosPreFiltering2 = std::vector<double>(3,0.0);

	/* Other pre-filtering associated variables. */
	std::vector<double> relPositionPreFiltering = std::vector<double>(3,0.0);
	std::vector<double> relVelocityPreFiltering = std::vector<double>(3,0.0);
	double relDistanceSquared, rMinSquared, distanceTravelledNormalToRmin, thresholdRadiusFine; // Look for the definitions of these in the pre-filters' implementation.

	/* Positions, velocities as well as relative distances and velocities at the Time of Closest Approach (TCA). */
	std::vector<double> posTCA1= std::vector<double>(3,0.0);
	std::vector<double> veloTCA1= std::vector<double>(3,0.0);
	std::vector<double> posTCA2= std::vector<double>(3,0.0);
	std::vector<double> veloTCA2 = std::vector<double>(3,0.0);
	std::vector<double> veloTCArelative= std::vector<double>(3,0.0);
	double coarseRelV, coarseMD; // From coarse interpolation.
	//double fineRelV, fineMD; // From fine interpolation.

	/* Temporary variables used for conjunction detection. */
	double coarseTCA, coarseTCA_JDAY; // From coarse interpolation.
	//double fineTCA, fineTCA_JDAY; // From fine interpolation. Don't use it in the current version.

	/*
	 *--------------------------------------------------------------------------------------------------------------------------------------------------
	 * PROPAGATE 1 COARSE STEP ITERATIVELY UNTIL THE END OF THE INTERVAL OF INTEREST.
	 *--------------------------------------------------------------------------------------------------------------------------------------------------
	 */
	while(CURRENT_JDAY < *analysisStop ){ // If CURRENT_JDAY is just a little less than analysisStop we'll find conjunctions outside of the analysis interval - need a check for this.
		VectorOperations::linspace(&currentJDAYs, CURRENT_JDAY, CURRENT_JDAY+COARSE_TIME_STEP, &NO_INTERPOLATION_POINTS); // Epochs that are interpolated between in the current interval.

		/*
		 *--------------------------------------------------------------------------------------------------------------------------------------------------
		 * RUNTIME PRE-FILTERING WITH SMART-SIEVE.
		 *--------------------------------------------------------------------------------------------------------------------------------------------------
		 */
		objectPairPtr->first.PropagateJDAY( &posesPreFiltering1, &velosPreFiltering1, currentJDAYs.at(0), true ); // Use the states and radii at the interval's beginning for pre-filtering.
		objectPairPtr->first.CalculateApogeeRadius(); objectPairPtr->first.CalculatePerigeeRadius();
		objectPairPtr->second.PropagateJDAY( &posesPreFiltering2, &velosPreFiltering2, currentJDAYs.at(0), true );
		objectPairPtr->second.CalculateApogeeRadius(); objectPairPtr->second.CalculatePerigeeRadius();
		
		// Check to see if neither object has re-entered. Use the state at interval's end.
		if( objectPairPtr->first.CURRENT_PERIGEE_RADIUS >= objectPairPtr->first.Re+100.0 || objectPairPtr->second.CURRENT_PERIGEE_RADIUS >= objectPairPtr->second.Re+100.0){
			// Add PERI_APOGEE_FILETR_PAD so as to keep more objects for further analysis.
			if( objectPairPtr->second.CURRENT_APOGEE_RADIUS+CONJUNCTION_THRESHOLD_DISTANCE+PERIGEE_APOGEE_PAD >= objectPairPtr->first.CURRENT_PERIGEE_RADIUS || objectPairPtr->first.CURRENT_APOGEE_RADIUS >= objectPairPtr->second.CURRENT_PERIGEE_RADIUS-PERIGEE_APOGEE_PAD-CONJUNCTION_THRESHOLD_DISTANCE){			
				// Perigee-apogee filter has been passed.
				VectorOperations::vectorDifference(&relPositionPreFiltering, &posesPreFiltering1, &posesPreFiltering2); // Position of object2 w.r.t. object1.
				if( !XYZpreFilter(&relPositionPreFiltering, &THRESHOLD_RADIUS) ){
					// XYZ filter has failed - skip some coarse time steps.
					CURRENT_JDAY = CURRENT_JDAY + COARSE_TIME_STEP* (int) 0.1*( VectorOperations::vectorMagnitude(&relPositionPreFiltering) - THRESHOLD_RADIUS )/2.0;
				}else{
					// XYZ filter has been passed.
					relDistanceSquared = VectorOperations::vectorMagnitudeSquared(&relPositionPreFiltering);
					if(relDistanceSquared <= THRESHOLD_RADIUS_SQUARED){
						// r^2 filter has been passed.
						VectorOperations::vectorDifference(&relVelocityPreFiltering, &velosPreFiltering1, &velosPreFiltering2); // True relative velocity.
						VectorOperations::unitVector(&relVelocityPreFiltering); // Make this vector unit (divide by length).
						distanceTravelledNormalToRmin = VectorOperations::dotProduct( &relPositionPreFiltering, &relVelocityPreFiltering); // Distance travelled by the secondary in the direction parallel to relative velocity (does not affect the miss distance).
						rMinSquared = relDistanceSquared - distanceTravelledNormalToRmin*distanceTravelledNormalToRmin; // Square of the actual miss distance.
						if(rMinSquared < ACCELERATION_RADIUS_SQUARED){
							// Minimum distance sieve has been passed.
							thresholdRadiusFine = ACCELERATION_RADIUS + 0.5*distanceTravelledNormalToRmin*COARSE_TIME_STEP_S; // Fine threshold radius (based on actual relative velocity) as defined by Rodriguez 2002, in km.
							if(relDistanceSquared <= thresholdRadiusFine*thresholdRadiusFine){
								/*
								 *--------------------------------------------------------------------------------------------------------------------------------------------------
								 * DO CONJUNCTION DETECTION IN THIS COARSE TIME STEP IF ALL FILTERS HAVE BEEN PASSED.
								 *--------------------------------------------------------------------------------------------------------------------------------------------------
								 */
								// Get the interpolation points of noth objects in the coarse time interval.
								// N.B. Cannot re-use positions from last propagation step as they will not be adjacent due to the pre-filters being implemented.
								for(size_t i=0; i<currentJDAYs.size(); i++){
									objectPairPtr->first.PropagateJDAY( &posesInterp1.at(i), &velosInterp1.at(i), currentJDAYs.at(i), false );
									objectPairPtr->second.PropagateJDAY( &posesInterp2.at(i), &velosInterp2.at(i), currentJDAYs.at(i), false );
								};
							
								// Find the interpolating polynomials' coefficients.
								findInterpolatingCoefficients(&posesInterp1, &velosInterp1, &currentJDAYs, interpCoeffsPointer1);
								findInterpolatingCoefficients(&posesInterp2, &velosInterp2, &currentJDAYs, interpCoeffsPointer2);
								for(size_t i=0; i<interpCoeffsPointerRelativePointer->size(); i++){
									VectorOperations::vectorDifference( &interpCoeffsPointerRelativePointer->at(i), &interpCoeffsPointer1->at(i), &interpCoeffsPointer2->at(i)); // Get the coefficients that give position of object 2 w.r.t. object 1.
								};

								// Find the TCA between this pairs of objects based on the interpolating polynomials' coefficients (TCA defined as the point where the relative velocity squared is 0.0).
								coarseTCA = findTCA(interpCoeffsPointerRelativePointer, 0.5);
							
								if( getInterpolatedRelativeDistanceSquared(coarseTCA) <= CONJUNCTION_THRESHOLD_DISTANCE_SQUARED && getInterpolatedRelativeDistanceSquared(coarseTCA)>=0 && 0.0<coarseTCA && coarseTCA<=1.0 ){ // If the distance is less than the threshold then need to find the closest approach precisely.
									/* 
									 * The conjunction passes the requirements in order to be recorded. Compute the information about it at the COARSE TCA. 
									 * Interpolate the positions and velocities rather than propagating again for speed.
									 */
									coarseTCA_JDAY = CURRENT_JDAY + coarseTCA*COARSE_TIME_STEP; // Julian Day of the TCA.
									
									// When running from a script there will be another simulation starting at the end epoch of this one. Need to avoid situations where conjunctions are found in that simulation's time domain. Cases like this may arise when the analysis intervals overlap ever so slightly.
									if(coarseTCA_JDAY < *analysisStop)
									{
										//Object 1.
										// The velocities aren't in km/s but in km/InterpolatinIntervalDuration => need to convert back to km/s.
										veloTCA1.at(0) = ( interpCoeffsPointer1->at(0).at(1) + interpCoeffsPointer1->at(0).at(2)*coarseTCA*2.0 + interpCoeffsPointer1->at(0).at(3)*coarseTCA*coarseTCA*3.0 ) / (currentJDAYs.back() - currentJDAYs.at(0))*SECOND_IN_JDAYS;
										veloTCA1.at(1) = ( interpCoeffsPointer1->at(1).at(1) + interpCoeffsPointer1->at(1).at(2)*coarseTCA*2.0 + interpCoeffsPointer1->at(1).at(3)*coarseTCA*coarseTCA*3.0 ) / (currentJDAYs.back() - currentJDAYs.at(0))*SECOND_IN_JDAYS;
										veloTCA1.at(2) = ( interpCoeffsPointer1->at(2).at(1) + interpCoeffsPointer1->at(2).at(2)*coarseTCA*2.0 + interpCoeffsPointer1->at(2).at(3)*coarseTCA*coarseTCA*3.0 ) / (currentJDAYs.back() - currentJDAYs.at(0))*SECOND_IN_JDAYS;

										//Object 2.
										veloTCA2.at(0) = ( interpCoeffsPointer2->at(0).at(1) + interpCoeffsPointer2->at(0).at(2)*coarseTCA*2.0 + interpCoeffsPointer2->at(0).at(3)*coarseTCA*coarseTCA*3.0 ) / (currentJDAYs.back() - currentJDAYs.at(0))*SECOND_IN_JDAYS;
										veloTCA2.at(1) = ( interpCoeffsPointer2->at(1).at(1) + interpCoeffsPointer2->at(1).at(2)*coarseTCA*2.0 + interpCoeffsPointer2->at(1).at(3)*coarseTCA*coarseTCA*3.0 ) / (currentJDAYs.back() - currentJDAYs.at(0))*SECOND_IN_JDAYS;
										veloTCA2.at(2) = ( interpCoeffsPointer2->at(2).at(1) + interpCoeffsPointer2->at(2).at(2)*coarseTCA*2.0 + interpCoeffsPointer2->at(2).at(3)*coarseTCA*coarseTCA*3.0 ) / (currentJDAYs.back() - currentJDAYs.at(0))*SECOND_IN_JDAYS;
	
										VectorOperations::vectorDifference( &veloTCArelative, &veloTCA1, &veloTCA2 );
										coarseRelV = VectorOperations::vectorMagnitude( &veloTCArelative ); // Need to know inertial velocities to compute the relative velocity as getInterpolatedRelativeRangeRateSquared will return 0 <=> TCA.
										coarseMD = std::sqrt( getInterpolatedRelativeDistanceSquared(coarseTCA) );

										// Create a Conjunction object that will store all the data about this event.
										Conjunction tempConj = Conjunction(coarseTCA_JDAY, coarseMD, coarseRelV, objectPairPtr->first.hardBodyRadius+objectPairPtr->second.hardBodyRadius, objectPairPtr->first.NORAD_ID, objectPairPtr->second.NORAD_ID);
										bool AppendConjunction = false; // Things may go badly when computing true collision probability due to zero eigenvalues, don't record such conjunctions.
										/* Choose how to compute the collision probability based on user-supplied argument. */
										switch( CovarianceType ){ 
											case 1: // Have no covariance information - assume a spherical uncertainty ellipsoid shape and compute the maximum probability for this.
												try{
													tempConj.calculateMaximumSphericalCollisionProbability(EnforceDirectNumericalIntegration);
													AppendConjunction = true;
												}catch(const MathWarning& MW){ // This will only happen when covariance scaling factor is very, very small and hence the scaled covariance matrix cannot be inversed.
													std::cerr<<"MathWarning when computing maximum spherical collision probability for "<<objectPairPtr->first.NORAD_ID<<" and "<<objectPairPtr->second.NORAD_ID<<", ignoring the conjunction."<<std::endl
														<<MW.what()<<std::endl;
												}catch(const FatalError& FE){
													std::cerr<<"FatalError when computing maximum spherical collision probability for "<<objectPairPtr->first.NORAD_ID<<" and "<<objectPairPtr->second.NORAD_ID<<", ignoring the conjunction."<<std::endl
													<<FE.what()<<std::endl;
												}
												break;
											case 2: // Estimate covariance from previous TLEs and compute true and maximum probabilities for that.
												// Propagate the covariance matrices of both objects to the epoch of the TCA.
												objectPairPtr->first.ComputeCovarianceOSW( coarseTCA_JDAY );
												objectPairPtr->second.ComputeCovarianceOSW( coarseTCA_JDAY );

												// For this CovarianceType also need to find the positions of the objects at the TCA.
												posTCA1.at(0) = interpCoeffsPointer1->at(0).at(0) + interpCoeffsPointer1->at(0).at(1)*coarseTCA + interpCoeffsPointer1->at(0).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(0).at(3)*coarseTCA*coarseTCA*coarseTCA;
	 											posTCA1.at(1) = interpCoeffsPointer1->at(1).at(0) + interpCoeffsPointer1->at(1).at(1)*coarseTCA + interpCoeffsPointer1->at(1).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(1).at(3)*coarseTCA*coarseTCA*coarseTCA;
	 											posTCA1.at(2) = interpCoeffsPointer1->at(2).at(0) + interpCoeffsPointer1->at(2).at(1)*coarseTCA + interpCoeffsPointer1->at(2).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(2).at(3)*coarseTCA*coarseTCA*coarseTCA;

												posTCA2.at(0) = interpCoeffsPointer2->at(0).at(0) + interpCoeffsPointer2->at(0).at(1)*coarseTCA + interpCoeffsPointer2->at(0).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(0).at(3)*coarseTCA*coarseTCA*coarseTCA;
	 											posTCA2.at(1) = interpCoeffsPointer2->at(1).at(0) + interpCoeffsPointer2->at(1).at(1)*coarseTCA + interpCoeffsPointer2->at(1).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(1).at(3)*coarseTCA*coarseTCA*coarseTCA;
	 											posTCA2.at(2) = interpCoeffsPointer2->at(2).at(0) + interpCoeffsPointer2->at(2).at(1)*coarseTCA + interpCoeffsPointer2->at(2).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(2).at(3)*coarseTCA*coarseTCA*coarseTCA;
										
												// Compute the true collision probability using the covariance data at TCA.
												try{
													tempConj.calculateCollisionProbability(&posTCA1, &veloTCA1, &objectPairPtr->first.PositionCovarianceMatrixRTC, &posTCA2, &veloTCA2, &objectPairPtr->second.PositionCovarianceMatrixRTC, EnforceDirectNumericalIntegration);
													AppendConjunction = true;
												}catch(const MathWarning& MW){
													std::cerr<<"MathWarning when computing true collision probability for "<<objectPairPtr->first.NORAD_ID<<" and "<<objectPairPtr->second.NORAD_ID<<", ignoring the conjunction."<<std::endl
														<<MW.what()<<std::endl;
												}catch(const FatalError& FE){
													std::cerr<<"FatalError when computing true collision probability for "<<objectPairPtr->first.NORAD_ID<<" and "<<objectPairPtr->second.NORAD_ID<<", ignoring the conjunction."<<std::endl
													<<FE.what()<<std::endl;
												}
												// Compute the maximum collision probability for this aspect ratio of the uncertainty ellipsoid.
												try{
													tempConj.calculateMaximumCollisionProbability(&posTCA1, &veloTCA1, &objectPairPtr->first.PositionCovarianceMatrixRTC, &posTCA2, &veloTCA2, &objectPairPtr->second.PositionCovarianceMatrixRTC, EnforceDirectNumericalIntegration);
													AppendConjunction = true;
												}catch(const MathWarning& MW){
													std::cerr<<"MathWarning when computing maximum collision probability for "<<objectPairPtr->first.NORAD_ID<<" and "<<objectPairPtr->second.NORAD_ID<<", ignoring the conjunction."<<std::endl
														<<MW.what()<<std::endl;
												}catch(const FatalError& FE){
													std::cerr<<"FatalError when computing maximum collision probability for "<<objectPairPtr->first.NORAD_ID<<" and "<<objectPairPtr->second.NORAD_ID<<", ignoring the conjunction."<<std::endl
													<<FE.what()<<std::endl;
												}
												break;
											case 3: case 4: // Use fixed European Space Surveiilance System performance or debugging 1 km STDEVs, i.e. don't propagate the covariance and go straight into probability calculations.
												// For this CovarianceType also need to find the positions of the objects at the TCA.
												posTCA1.at(0) = interpCoeffsPointer1->at(0).at(0) + interpCoeffsPointer1->at(0).at(1)*coarseTCA + interpCoeffsPointer1->at(0).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(0).at(3)*coarseTCA*coarseTCA*coarseTCA;
	 											posTCA1.at(1) = interpCoeffsPointer1->at(1).at(0) + interpCoeffsPointer1->at(1).at(1)*coarseTCA + interpCoeffsPointer1->at(1).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(1).at(3)*coarseTCA*coarseTCA*coarseTCA;
	 											posTCA1.at(2) = interpCoeffsPointer1->at(2).at(0) + interpCoeffsPointer1->at(2).at(1)*coarseTCA + interpCoeffsPointer1->at(2).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(2).at(3)*coarseTCA*coarseTCA*coarseTCA;

												posTCA2.at(0) = interpCoeffsPointer2->at(0).at(0) + interpCoeffsPointer2->at(0).at(1)*coarseTCA + interpCoeffsPointer2->at(0).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(0).at(3)*coarseTCA*coarseTCA*coarseTCA;
	 											posTCA2.at(1) = interpCoeffsPointer2->at(1).at(0) + interpCoeffsPointer2->at(1).at(1)*coarseTCA + interpCoeffsPointer2->at(1).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(1).at(3)*coarseTCA*coarseTCA*coarseTCA;
	 											posTCA2.at(2) = interpCoeffsPointer2->at(2).at(0) + interpCoeffsPointer2->at(2).at(1)*coarseTCA + interpCoeffsPointer2->at(2).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(2).at(3)*coarseTCA*coarseTCA*coarseTCA;
	
												// Compute the true collision probability using the covariance data at TCA.
												try{
													tempConj.calculateCollisionProbability(&posTCA1, &veloTCA1, &objectPairPtr->first.PositionCovarianceMatrixRTC, &posTCA2, &veloTCA2, &objectPairPtr->second.PositionCovarianceMatrixRTC, EnforceDirectNumericalIntegration);
													AppendConjunction = true;
												}catch(const MathWarning& MW){
													std::cerr<<"MathWarning when computing true collision probability for "<<objectPairPtr->first.NORAD_ID<<" and "<<objectPairPtr->second.NORAD_ID<<", ignoring the conjunction."<<std::endl
													<<MW.what()<<std::endl;
												}catch(const FatalError& FE){
													std::cerr<<"FatalError when computing true collision probability for "<<objectPairPtr->first.NORAD_ID<<" and "<<objectPairPtr->second.NORAD_ID<<", ignoring the conjunction."<<std::endl
													<<FE.what()<<std::endl;
												}
												// Compute the maximum collision probability for this aspect ratio of the uncertainty ellipsoid.
												try{
													tempConj.calculateMaximumCollisionProbability(&posTCA1, &veloTCA1, &objectPairPtr->first.PositionCovarianceMatrixRTC, &posTCA2, &veloTCA2, &objectPairPtr->second.PositionCovarianceMatrixRTC, EnforceDirectNumericalIntegration);
													AppendConjunction = true;
												}catch(const MathWarning& MW){
													std::cerr<<"MathWarning when computing maximum collision probability for "<<objectPairPtr->first.NORAD_ID<<" and "<<objectPairPtr->second.NORAD_ID<<", ignoring the conjunction."<<std::endl
														<<MW.what()<<std::endl;
												}catch(const FatalError& FE){
													std::cerr<<"FatalError when computing maximum collision probability for "<<objectPairPtr->first.NORAD_ID<<" and "<<objectPairPtr->second.NORAD_ID<<", ignoring the conjunction."<<std::endl
													<<FE.what()<<std::endl;
												}
												break;
											case 5: // Don't calculate PC at all, we only want to find the miss distances.
												tempConj.maximumProbability = -1.0; // Make it clear that we've not computed the PC.
												tempConj.trueProbability = -1.0;
												AppendConjunction = true; // Nothing can go wrong here.
												break;
											case 6: // Use time-invariant, externally set covariance matrices. Go straight into Pc calculations and skip PcMAX.
												// For this CovarianceType also need to find the positions of the objects at the TCA.
												posTCA1.at(0) = interpCoeffsPointer1->at(0).at(0) + interpCoeffsPointer1->at(0).at(1)*coarseTCA + interpCoeffsPointer1->at(0).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(0).at(3)*coarseTCA*coarseTCA*coarseTCA;
	 											posTCA1.at(1) = interpCoeffsPointer1->at(1).at(0) + interpCoeffsPointer1->at(1).at(1)*coarseTCA + interpCoeffsPointer1->at(1).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(1).at(3)*coarseTCA*coarseTCA*coarseTCA;
	 											posTCA1.at(2) = interpCoeffsPointer1->at(2).at(0) + interpCoeffsPointer1->at(2).at(1)*coarseTCA + interpCoeffsPointer1->at(2).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(2).at(3)*coarseTCA*coarseTCA*coarseTCA;

												posTCA2.at(0) = interpCoeffsPointer2->at(0).at(0) + interpCoeffsPointer2->at(0).at(1)*coarseTCA + interpCoeffsPointer2->at(0).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(0).at(3)*coarseTCA*coarseTCA*coarseTCA;
	 											posTCA2.at(1) = interpCoeffsPointer2->at(1).at(0) + interpCoeffsPointer2->at(1).at(1)*coarseTCA + interpCoeffsPointer2->at(1).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(1).at(3)*coarseTCA*coarseTCA*coarseTCA;
	 											posTCA2.at(2) = interpCoeffsPointer2->at(2).at(0) + interpCoeffsPointer2->at(2).at(1)*coarseTCA + interpCoeffsPointer2->at(2).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(2).at(3)*coarseTCA*coarseTCA*coarseTCA;
	
												// Compute the true collision probability using the covariance data at TCA.
												try{
													tempConj.calculateCollisionProbability(&posTCA1, &veloTCA1, &objectPairPtr->first.PositionCovarianceMatrixRTC, &posTCA2, &veloTCA2, &objectPairPtr->second.PositionCovarianceMatrixRTC, EnforceDirectNumericalIntegration);
													AppendConjunction = true;
												}catch( FatalError ){
													std::cerr<<"Impossible to compute true collision probability for "<<objectPairPtr->first.NORAD_ID<<" and "<<objectPairPtr->second.NORAD_ID<<", ignoring the conjunction."<<std::endl;
												}
											
												tempConj.maximumProbability = -1.0; // Make it clear that we've not computed PcMAX.
												break;
											
										}//END switch covariance type
										if(AppendConjunction){
											conjunctionsPtr->push_back( tempConj ); // Store in the results vector.
											//std::cout<<tempConj<<std::endl; // Print to a log file.
										}// No need to write again that this conjunction has been found but will be ignored.
									};//END if(coarseTCA_JDAY < *analysisStop)
								};//END if( getInterpolatedRelativeDistanceSquared(coarseTCA) <= CONJUNCTION_THRESHOLD_DISTANCE_SQUARED && getInterpolatedRelativeDistanceSquared(coarseTCA)>=0 && 0.0<coarseTCA && coarseTCA<=1.0 )
							};//END fine r^2 sieve
						};//END minimum distance sieve
					};//END r^2 sieve
				};//END XYZ sieve
			}else if( objectPairPtr->second.CURRENT_APOGEE_RADIUS+PERIGEE_APOGEE_PAD+PERIGEE_APOGEE_ORBITAL_REGIME_PAD >= objectPairPtr->first.CURRENT_PERIGEE_RADIUS && objectPairPtr->first.CURRENT_APOGEE_RADIUS >= objectPairPtr->second.CURRENT_PERIGEE_RADIUS-PERIGEE_APOGEE_PAD-PERIGEE_APOGEE_ORBITAL_REGIME_PAD ){
				break; // Perigee-apogee filter hasn't been passed by A LOT => objects are in differnet orbital regimes, can discard the pair straight away.
			};//END perigee/apogee sieve
		};//END check if neither object has re-entered.
		CURRENT_JDAY = CURRENT_JDAY + COARSE_TIME_STEP; // Go to the next coarse analysis interval.
	};//END while(CURRENT_JDAY <= *analysisStop )
delete interpCoeffsPointer1; // Deallocate the vectors.
delete interpCoeffsPointer2;
delete interpCoeffsPointerRelativePointer;
};//END findConjunctionsBetweenObjectPair

void Simulation::findConjunctionsBetweenObjectPair(std::pair<ExternalSpaceObject, SpaceObject>* objectPairPtr, std::vector<Conjunction>* conjunctionsPtr){
	/* Look for all the conjunctions that are predicted to happen between the given pair of objects in the time interval specified.
	Assume that objects that reach altitude of 100.0 km have re-entered the atmosphere.
	Use perigee/apogee sieve (computed at the end of the interpolation interval with a pad of PERIGEE_APOGEE_PAD km) with the caveat 
	that if the difference in altitudes	is 1000.0 km or more the pair will be ignored for speed.
	Also use the "smart sieve" as described in ROdriguez et. all 2002, namely:
		- XYZ
		- r^2
		- minimum distance
		- fine r^2
	First try to discard the pair at every coarse time step based on the prefilters. If all of them are passed find the conjunction
	based purely on the range measure i.e. find the epoch when the relative distance is minimum and below CONJUNCTION_THRESHOLD_DISTANCE.
	Once a conjunction is found its detailed information is refined by propagating the objects in a shorter interval around it.
	@param objectPairPtr - pointer to a pair of objects <object1, object2> between which the conjunctions are to be found.
	@param conujnctionsPtr - vector where all the found conjunctions will be recorded.
	@param analysisStart, analysisStop - Julian days marking the beginning and the end of the analysis intetrval, respectively.
	@param noInterpolationPoints - number of points to be used for interpolation of trajectories within a singe coarse analysis time step.
	@param COARSE_TIME_STEP_S - the duration of the analysis' time step (can, in principle, be different than the interpolation interval if more than 2 interpolation nodes are used).
	@param COARSE_TIME_STEP - same as COARSE_TIME_STEP_S but in Julian days.
	@param PERI_APOGEE_FILETR_PAD - buffer to be used when using the apogee/perigee prefilter, in km.
	@param THRESHOLD_RADIUS - the threshold radius as defined by Rodriguez 2002, in km.
	@param THRESHOLD_RADIUS_SQUARED - square of the threshold radius as defined by Rodriguez 2002, in km.
	@param ACCELERATION_RADIUS - the acceleration radius as defined by Rodriquez 2002, in km.
	@param ACCELERATION_RADIUS_SQUARED - square of the acceleration radius as defined by Rodriquez 2002, in km.
	*/

	/* Positions, velocities and epochs of objects in Cartesian space used for interpolation. */
	std::vector< std::vector<double> > posesInterp1 = std::vector< std::vector<double> >(NO_INTERPOLATION_POINTS, std::vector<double>(3,0.0) ); // Each row has 3 Cartesian components.
	std::vector< std::vector<double> > velosInterp1 = std::vector< std::vector<double> >(NO_INTERPOLATION_POINTS, std::vector<double>(3,0.0) );
	std::vector< std::vector<double> > posesInterp2 = std::vector< std::vector<double> >(NO_INTERPOLATION_POINTS, std::vector<double>(3,0.0) );
	std::vector< std::vector<double> > velosInterp2 = std::vector< std::vector<double> >(NO_INTERPOLATION_POINTS, std::vector<double>(3,0.0) );
	std::vector<double> currentJDAYs(NO_INTERPOLATION_POINTS, 0.0);

	/* Interpolating polynomials' coefficients for both objects - row indices correspond to dimensions, column indices to the particular coefficients. */
	std::vector< std::vector<double> >* interpCoeffsPointer1 = new std::vector< std::vector<double> >(3, std::vector<double>(2*NO_INTERPOLATION_POINTS, 0.0) );
	std::vector< std::vector<double> >* interpCoeffsPointer2 = new std::vector< std::vector<double> >(3, std::vector<double>(2*NO_INTERPOLATION_POINTS, 0.0) );
	std::vector< std::vector<double> >* interpCoeffsPointerRelativePointer = new std::vector< std::vector<double> >(3, std::vector<double>(2*NO_INTERPOLATION_POINTS, 0.0) );

	/* Positions and velocities of the second object in Cartesian space used for pre-filtering. Get positions of the primary (first object) from the ephemeris table directly. */
	std::vector<double> posesPreFiltering2 = std::vector<double>(3,0.0); std::vector<double> velosPreFiltering2 = std::vector<double>(3,0.0);

	/* Other pre-filtering associated variables. */
	std::vector<double> relPositionPreFiltering = std::vector<double>(3,0.0);
	std::vector<double> relVelocityPreFiltering = std::vector<double>(3,0.0);
	double relDistanceSquared, rMinSquared, distanceTravelledNormalToRmin, thresholdRadiusFine; // Look for the definitions of these in the pre-filters' implementation.

	/* Positions, velocities as well as relative distances and velocities at the Time of Closest Approach (TCA). */
	std::vector<double> posTCA1= std::vector<double>(3,0.0);
	std::vector<double> veloTCA1= std::vector<double>(3,0.0);
	std::vector<double> posTCA2= std::vector<double>(3,0.0);
	std::vector<double> veloTCA2 = std::vector<double>(3,0.0);
	std::vector<double> veloTCArelative= std::vector<double>(3,0.0);
	double coarseRelV, coarseMD; // From coarse interpolation.

	/* Temporary variables used for conjunction detection. */
	double coarseTCA, coarseTCA_JDAY; // From coarse interpolation.

	if( NO_INTERPOLATION_POINTS!=2 ){ // That's how it is.
		std::string msg = "Requested interpolation of the ephemeris table based on more than two points. This is not supported in the current version. Please use two interpolation points when supplying an external ephemeris table.";
		throw FatalError(msg,__FILE__,__LINE__);
	}

	/*
	 *--------------------------------------------------------------------------------------------------------------------------------------------------
	 * PROPAGATE 1 COARSE STEP ITERATIVELY UNTIL THE END OF THE INTERVAL OF INTEREST.
	 *--------------------------------------------------------------------------------------------------------------------------------------------------
	 */
	for( std::vector<int>::size_type analysisEpochID=0; analysisEpochID<objectPairPtr->first.EphemerisEpochs.size()-1; analysisEpochID++){ // Go through all the epochs from the ephemeris file and find the conjunctions between the objects in intervals given by them. Current epoch interval lies between i and i+1.
		/*
		 *--------------------------------------------------------------------------------------------------------------------------------------------------
		 * RUNTIME PRE-FILTERING WITH SMART-SIEVE.
		 *--------------------------------------------------------------------------------------------------------------------------------------------------
		 */
		objectPairPtr->first.CalculateRadius( &objectPairPtr->first.EphemerisPositions.at(analysisEpochID) ); // Do the prefiltering based on the states at the beginning of this interval given by epochs at i and i+1.
		objectPairPtr->second.PropagateJDAY( &posesPreFiltering2, &velosPreFiltering2, objectPairPtr->first.EphemerisEpochs.at(analysisEpochID), true );
		objectPairPtr->second.CalculateApogeeRadius(); objectPairPtr->second.CalculatePerigeeRadius();
		
		currentJDAYs.at(0) = objectPairPtr->first.EphemerisEpochs.at(analysisEpochID); // These are the bounds of the current analysis interval.
		currentJDAYs.at(1) = objectPairPtr->first.EphemerisEpochs.at(analysisEpochID+1);

		// Check to see if neither object has re-entered. Use the state at interval's end.
		if( objectPairPtr->first.CURRENT_RADIUS >= objectPairPtr->first.Re+100.0 || objectPairPtr->second.CURRENT_PERIGEE_RADIUS >= objectPairPtr->second.Re+100.0){
			// Add PERI_APOGEE_FILETR_PAD so as to keep more objects for further analysis.
			if( objectPairPtr->second.CURRENT_APOGEE_RADIUS+CONJUNCTION_THRESHOLD_DISTANCE+PERIGEE_APOGEE_PAD >= objectPairPtr->first.CURRENT_RADIUS || objectPairPtr->first.CURRENT_RADIUS >= objectPairPtr->second.CURRENT_PERIGEE_RADIUS-PERIGEE_APOGEE_PAD-CONJUNCTION_THRESHOLD_DISTANCE){
				// Perigee-apogee filter has been passed.
				VectorOperations::vectorDifference(&relPositionPreFiltering, &objectPairPtr->first.EphemerisPositions.at(analysisEpochID), &posesPreFiltering2); // Position of object2 w.r.t. object1.
				if( !XYZpreFilter(&relPositionPreFiltering, &THRESHOLD_RADIUS) ){
					// XYZ filter has failed - skip some coarse time steps.
					analysisEpochID = analysisEpochID + (std::vector<int>::size_type) std::floor( 0.1*( VectorOperations::vectorMagnitude(&relPositionPreFiltering) - THRESHOLD_RADIUS )/2.0 ); // Simply ignore the next few ephemeris intervals.
				}else{
					// XYZ filter has been passed.
					relDistanceSquared = VectorOperations::vectorMagnitudeSquared(&relPositionPreFiltering);
					if(relDistanceSquared <= THRESHOLD_RADIUS_SQUARED){
						// r^2 filter has been passed.
						VectorOperations::vectorDifference(&relVelocityPreFiltering, &objectPairPtr->first.EphemerisVelocities.at(analysisEpochID), &velosPreFiltering2); // True relative velocity.
						VectorOperations::unitVector(&relVelocityPreFiltering); // Make this vector unit (divide by length).
						distanceTravelledNormalToRmin = VectorOperations::dotProduct( &relPositionPreFiltering, &relVelocityPreFiltering); // Distance travelled by the secondary in the direction parallel to relative velocity (does not affect the miss distance).
						rMinSquared = relDistanceSquared - distanceTravelledNormalToRmin*distanceTravelledNormalToRmin; // Square of the actual miss distance.
						if(rMinSquared < ACCELERATION_RADIUS_SQUARED){
							// Minimum distance sieve has been passed.
							thresholdRadiusFine = ACCELERATION_RADIUS + 0.5*distanceTravelledNormalToRmin*COARSE_TIME_STEP_S; // Fine threshold radius (based on actual relative velocity) as defined by Rodriguez 2002, in km.
							if(relDistanceSquared <= thresholdRadiusFine*thresholdRadiusFine){
								// Fine r^2 sieve has been passed.
							
								/*
								 *--------------------------------------------------------------------------------------------------------------------------------------------------
								 * DO CONJUNCTION DETECTION IN THIS COARSE TIME STEP IF ALL FILTERS HAVE BEEN PASSED.
								 *--------------------------------------------------------------------------------------------------------------------------------------------------
								 */
								// Get the interpolation points of noth objects in the coarse time interval.
								// N.B. Cannot re-use positions from last propagation step as they will not be adjacent due to the pre-filters being implemented.
								for(size_t i=0; i<currentJDAYs.size(); i++){
									posesInterp1.at(i) = objectPairPtr->first.EphemerisPositions.at(analysisEpochID+i); // This could be easily adapted for support of more than 2 interpolation points. Currently will will always give the position at analysisEpochID and the next one i.e. the ones that mark  the ends of this analysis interval.
									velosInterp1.at(i) = objectPairPtr->first.EphemerisVelocities.at(analysisEpochID+i);
									objectPairPtr->second.PropagateJDAY( &posesInterp2.at(i), &velosInterp2.at(i), currentJDAYs.at(i), false );
								};
							
								// Find the interpolating polynomials' coefficients.
								findInterpolatingCoefficients(&posesInterp1, &velosInterp1, &currentJDAYs, interpCoeffsPointer1);
								findInterpolatingCoefficients(&posesInterp2, &velosInterp2, &currentJDAYs, interpCoeffsPointer2);
								for(size_t i=0; i<interpCoeffsPointerRelativePointer->size(); i++){
									VectorOperations::vectorDifference( &interpCoeffsPointerRelativePointer->at(i), &interpCoeffsPointer1->at(i), &interpCoeffsPointer2->at(i)); // Get the coefficients that give position of object 2 w.r.t. object 1.
								};

								// Find the TCA between this pairs of objects based on the interpolating polynomials' coefficients (TCA defined as the point where the relative velocity squared is 0.0).
								coarseTCA = findTCA(interpCoeffsPointerRelativePointer, 0.5);
							
								if( getInterpolatedRelativeDistanceSquared(coarseTCA) <= CONJUNCTION_THRESHOLD_DISTANCE_SQUARED && getInterpolatedRelativeDistanceSquared(coarseTCA)>=0 && 0.0<coarseTCA && coarseTCA<=1.0 ){ // If the distance is less than the threshold then need to find the closest approach precisely.
									/* 
									 * The conjunction passes the requirements in order to be recorded. Compute the information about it at the COARSE TCA. 
									 * Interpolate the positions and velocities rather than propagating again for speed.
									 */
									coarseTCA_JDAY = objectPairPtr->first.EphemerisEpochs.at(analysisEpochID) + coarseTCA*COARSE_TIME_STEP; // Julian Day of the TCA.
								    
									//Object 1.
									// The velocities aren't in km/s but in km/InterpolatinIntervalDuration => need to convert back to km/s.
									veloTCA1.at(0) = ( interpCoeffsPointer1->at(0).at(1) + interpCoeffsPointer1->at(0).at(2)*coarseTCA*2.0 + interpCoeffsPointer1->at(0).at(3)*coarseTCA*coarseTCA*3.0 ) / (currentJDAYs.back() - currentJDAYs.at(0))*SECOND_IN_JDAYS;
									veloTCA1.at(1) = ( interpCoeffsPointer1->at(1).at(1) + interpCoeffsPointer1->at(1).at(2)*coarseTCA*2.0 + interpCoeffsPointer1->at(1).at(3)*coarseTCA*coarseTCA*3.0 ) / (currentJDAYs.back() - currentJDAYs.at(0))*SECOND_IN_JDAYS;
									veloTCA1.at(2) = ( interpCoeffsPointer1->at(2).at(1) + interpCoeffsPointer1->at(2).at(2)*coarseTCA*2.0 + interpCoeffsPointer1->at(2).at(3)*coarseTCA*coarseTCA*3.0 ) / (currentJDAYs.back() - currentJDAYs.at(0))*SECOND_IN_JDAYS;

									//Object 2.
									veloTCA2.at(0) = ( interpCoeffsPointer2->at(0).at(1) + interpCoeffsPointer2->at(0).at(2)*coarseTCA*2.0 + interpCoeffsPointer2->at(0).at(3)*coarseTCA*coarseTCA*3.0 ) / (currentJDAYs.back() - currentJDAYs.at(0))*SECOND_IN_JDAYS;
									veloTCA2.at(1) = ( interpCoeffsPointer2->at(1).at(1) + interpCoeffsPointer2->at(1).at(2)*coarseTCA*2.0 + interpCoeffsPointer2->at(1).at(3)*coarseTCA*coarseTCA*3.0 ) / (currentJDAYs.back() - currentJDAYs.at(0))*SECOND_IN_JDAYS;
									veloTCA2.at(2) = ( interpCoeffsPointer2->at(2).at(1) + interpCoeffsPointer2->at(2).at(2)*coarseTCA*2.0 + interpCoeffsPointer2->at(2).at(3)*coarseTCA*coarseTCA*3.0 ) / (currentJDAYs.back() - currentJDAYs.at(0))*SECOND_IN_JDAYS;
	
									VectorOperations::vectorDifference( &veloTCArelative, &veloTCA1, &veloTCA2 );
									coarseRelV = VectorOperations::vectorMagnitude( &veloTCArelative ); // Need to know inertial velocities to compute the relative velocity as getInterpolatedRelativeRangeRateSquared will return 0 <=> TCA.
									coarseMD = std::sqrt( getInterpolatedRelativeDistanceSquared(coarseTCA) );

									// Create a Conjunction object that will store all the data about this event.
									Conjunction tempConj = Conjunction(coarseTCA_JDAY, coarseMD, coarseRelV, objectPairPtr->first.hardBodyRadius+objectPairPtr->second.hardBodyRadius, objectPairPtr->first.NORAD_ID, objectPairPtr->second.NORAD_ID);
									
									/* Choose how to compute the collision probability based on user-supplied argument. */
									bool AppendConjunction = false;  // Things may go badly when computing true collision probability due to zero eigenvalues, don't record such conjunctions.
									switch( CovarianceType ){ 
										case 1: // Have no covariance information - assume a spherical uncertainty ellipsoid shape and compute the maximum probability for this.
											tempConj.calculateMaximumSphericalCollisionProbability(EnforceDirectNumericalIntegration);
											AppendConjunction = true;
											break;
										case 2: // Estimate covariance from previous TLEs and compute true and maximum probabilities for that.
											// Propagate the covariance matrices of both objects to the epoch of the TCA.
											//UNDONE: figure out what to do when covariance of the ExternalSpaceObject is to be propagated to some epoch.
											//objectPairPtr->first.ComputeCovarianceOSW( coarseTCA_JDAY );
											objectPairPtr->second.ComputeCovarianceOSW( coarseTCA_JDAY );

											// For this CovarianceType also need to find the positions of the objects at the TCA.
											posTCA1.at(0) = interpCoeffsPointer1->at(0).at(0) + interpCoeffsPointer1->at(0).at(1)*coarseTCA + interpCoeffsPointer1->at(0).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(0).at(3)*coarseTCA*coarseTCA*coarseTCA;
 											posTCA1.at(1) = interpCoeffsPointer1->at(1).at(0) + interpCoeffsPointer1->at(1).at(1)*coarseTCA + interpCoeffsPointer1->at(1).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(1).at(3)*coarseTCA*coarseTCA*coarseTCA;
 											posTCA1.at(2) = interpCoeffsPointer1->at(2).at(0) + interpCoeffsPointer1->at(2).at(1)*coarseTCA + interpCoeffsPointer1->at(2).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(2).at(3)*coarseTCA*coarseTCA*coarseTCA;

											posTCA2.at(0) = interpCoeffsPointer2->at(0).at(0) + interpCoeffsPointer2->at(0).at(1)*coarseTCA + interpCoeffsPointer2->at(0).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(0).at(3)*coarseTCA*coarseTCA*coarseTCA;
 											posTCA2.at(1) = interpCoeffsPointer2->at(1).at(0) + interpCoeffsPointer2->at(1).at(1)*coarseTCA + interpCoeffsPointer2->at(1).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(1).at(3)*coarseTCA*coarseTCA*coarseTCA;
 											posTCA2.at(2) = interpCoeffsPointer2->at(2).at(0) + interpCoeffsPointer2->at(2).at(1)*coarseTCA + interpCoeffsPointer2->at(2).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(2).at(3)*coarseTCA*coarseTCA*coarseTCA;
										
											// Compute the true collision probability using the covariance data at TCA.
											try{
												tempConj.calculateCollisionProbability(&posTCA1, &veloTCA1, &objectPairPtr->first.PositionCovarianceMatrixRTC, &posTCA2, &veloTCA2, &objectPairPtr->second.PositionCovarianceMatrixRTC, EnforceDirectNumericalIntegration);
												AppendConjunction = true;
											}catch( FatalError ){
												std::cerr<<"Impossible to compute true collision probability for "<<objectPairPtr->first.NORAD_ID<<" and "<<objectPairPtr->second.NORAD_ID<<" at "<<coarseTCA_JDAY<<std::endl;
											}
											// Compute the maximum collision probability for this aspect ratio of the uncertainty ellipsoid.
											try{
												tempConj.calculateMaximumCollisionProbability(&posTCA1, &veloTCA1, &objectPairPtr->first.PositionCovarianceMatrixRTC, &posTCA2, &veloTCA2, &objectPairPtr->second.PositionCovarianceMatrixRTC, EnforceDirectNumericalIntegration);
												AppendConjunction = true;
											}catch( FatalError ){
												std::cerr<<"Impossible to compute maximum collision probability for "<<objectPairPtr->first.NORAD_ID<<" and "<<objectPairPtr->second.NORAD_ID<<" at "<<coarseTCA_JDAY<<std::endl;
											}

											break;
										case 3: case 4: // Use fixed European Space Surveiilance System performance or debugging 1 km STDEVs, i.e. don't propagate the covariance and go straight into probability calculations.
											// For this CovarianceType also need to find the positions of the objects at the TCA.
											posTCA1.at(0) = interpCoeffsPointer1->at(0).at(0) + interpCoeffsPointer1->at(0).at(1)*coarseTCA + interpCoeffsPointer1->at(0).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(0).at(3)*coarseTCA*coarseTCA*coarseTCA;
 											posTCA1.at(1) = interpCoeffsPointer1->at(1).at(0) + interpCoeffsPointer1->at(1).at(1)*coarseTCA + interpCoeffsPointer1->at(1).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(1).at(3)*coarseTCA*coarseTCA*coarseTCA;
 											posTCA1.at(2) = interpCoeffsPointer1->at(2).at(0) + interpCoeffsPointer1->at(2).at(1)*coarseTCA + interpCoeffsPointer1->at(2).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(2).at(3)*coarseTCA*coarseTCA*coarseTCA;

											posTCA2.at(0) = interpCoeffsPointer2->at(0).at(0) + interpCoeffsPointer2->at(0).at(1)*coarseTCA + interpCoeffsPointer2->at(0).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(0).at(3)*coarseTCA*coarseTCA*coarseTCA;
 											posTCA2.at(1) = interpCoeffsPointer2->at(1).at(0) + interpCoeffsPointer2->at(1).at(1)*coarseTCA + interpCoeffsPointer2->at(1).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(1).at(3)*coarseTCA*coarseTCA*coarseTCA;
 											posTCA2.at(2) = interpCoeffsPointer2->at(2).at(0) + interpCoeffsPointer2->at(2).at(1)*coarseTCA + interpCoeffsPointer2->at(2).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(2).at(3)*coarseTCA*coarseTCA*coarseTCA;
	
											// Compute the true collision probability using the covariance data at TCA.
											try{
												tempConj.calculateCollisionProbability(&posTCA1, &veloTCA1, &objectPairPtr->first.PositionCovarianceMatrixRTC, &posTCA2, &veloTCA2, &objectPairPtr->second.PositionCovarianceMatrixRTC, EnforceDirectNumericalIntegration);
												AppendConjunction = true;
											}catch( FatalError ){
												std::cerr<<"Impossible to compute true collision probability for "<<objectPairPtr->first.NORAD_ID<<" and "<<objectPairPtr->second.NORAD_ID<<", ignoring the conjunction."<<std::endl;
											}
											// Compute the maximum collision probability for this aspect ratio of the uncertainty ellipsoid.
											try{
												tempConj.calculateMaximumCollisionProbability(&posTCA1, &veloTCA1, &objectPairPtr->first.PositionCovarianceMatrixRTC, &posTCA2, &veloTCA2, &objectPairPtr->second.PositionCovarianceMatrixRTC, EnforceDirectNumericalIntegration);
												AppendConjunction = true;
											}catch( FatalError ){
												std::cerr<<"Impossible to compute maximum collision probability for "<<objectPairPtr->first.NORAD_ID<<" and "<<objectPairPtr->second.NORAD_ID<<", ignoring the conjunction."<<std::endl;
											}
											break;
										case 5: // Don't calculate PC at all, we only want to find the miss distances.
											tempConj.maximumProbability = -1.0; // Make it clear that we've not computed the PC.
											tempConj.trueProbability = -1.0;
											AppendConjunction = true; // Nothing can go wrong here.
											break;
										case 6: // Use time-invariant, externally set covariance matrices. Go straight into Pc calculations and skip PcMAX.
											// For this CovarianceType also need to find the positions of the objects at the TCA.
											posTCA1.at(0) = interpCoeffsPointer1->at(0).at(0) + interpCoeffsPointer1->at(0).at(1)*coarseTCA + interpCoeffsPointer1->at(0).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(0).at(3)*coarseTCA*coarseTCA*coarseTCA;
 											posTCA1.at(1) = interpCoeffsPointer1->at(1).at(0) + interpCoeffsPointer1->at(1).at(1)*coarseTCA + interpCoeffsPointer1->at(1).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(1).at(3)*coarseTCA*coarseTCA*coarseTCA;
 											posTCA1.at(2) = interpCoeffsPointer1->at(2).at(0) + interpCoeffsPointer1->at(2).at(1)*coarseTCA + interpCoeffsPointer1->at(2).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer1->at(2).at(3)*coarseTCA*coarseTCA*coarseTCA;

											posTCA2.at(0) = interpCoeffsPointer2->at(0).at(0) + interpCoeffsPointer2->at(0).at(1)*coarseTCA + interpCoeffsPointer2->at(0).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(0).at(3)*coarseTCA*coarseTCA*coarseTCA;
 											posTCA2.at(1) = interpCoeffsPointer2->at(1).at(0) + interpCoeffsPointer2->at(1).at(1)*coarseTCA + interpCoeffsPointer2->at(1).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(1).at(3)*coarseTCA*coarseTCA*coarseTCA;
 											posTCA2.at(2) = interpCoeffsPointer2->at(2).at(0) + interpCoeffsPointer2->at(2).at(1)*coarseTCA + interpCoeffsPointer2->at(2).at(2)*coarseTCA*coarseTCA + interpCoeffsPointer2->at(2).at(3)*coarseTCA*coarseTCA*coarseTCA;
	
											// Compute the true collision probability using the covariance data at TCA.
											try{
												tempConj.calculateCollisionProbability(&posTCA1, &veloTCA1, &objectPairPtr->first.PositionCovarianceMatrixRTC, &posTCA2, &veloTCA2, &objectPairPtr->second.PositionCovarianceMatrixRTC, EnforceDirectNumericalIntegration);
												AppendConjunction = true;
											}catch( FatalError ){
												std::cerr<<"Impossible to compute true collision probability for "<<objectPairPtr->first.NORAD_ID<<" and "<<objectPairPtr->second.NORAD_ID<<", ignoring the conjunction."<<std::endl;
											}
											
											tempConj.maximumProbability = -1.0; // Make it clear that we've not computed PcMAX.
											break;
											
									}//END switch covariance type

									if(AppendConjunction){
										conjunctionsPtr->push_back( tempConj ); // Store in the results vector.
										//std::cout<<tempConj<<std::endl; // Print to a log file.
									}else{
										std::cerr<<"A conjunction won't be recorded as at least one collision probability wasn't successfully computed."<<std::endl;
									}
								
								};//END if( getInterpolatedRelativeDistanceSquared(coarseTCA) <= CONJUNCTION_THRESHOLD_DISTANCE_SQUARED && getInterpolatedRelativeDistanceSquared(coarseTCA)>=0 && 0.0<coarseTCA && coarseTCA<=1.0 )
							};//END fine r^2 sieve
						};//END minimum distance sieve
					};//END r^2 sieve
				};//END XYZ sieve
			}else if( objectPairPtr->second.CURRENT_APOGEE_RADIUS+PERIGEE_APOGEE_PAD+PERIGEE_APOGEE_ORBITAL_REGIME_PAD >= objectPairPtr->first.CURRENT_RADIUS && objectPairPtr->first.CURRENT_RADIUS >= objectPairPtr->second.CURRENT_PERIGEE_RADIUS-PERIGEE_APOGEE_PAD-PERIGEE_APOGEE_ORBITAL_REGIME_PAD ){
				break; // Perigee-apogee filter hasn't been passed by A LOT => objects are in differnet orbital regimes, can discard the pair straight away.
			};//END perigee/apogee sieve
		};//END check if neither object has re-entered.
	};//END for all the epochs from the ephemeris table.
delete interpCoeffsPointer1; // Deallocate the vectors.
delete interpCoeffsPointer2;
delete interpCoeffsPointerRelativePointer;
};//END findConjunctionsBetweenObjectPair

double Simulation::findTCA(std::vector< std::vector<double> >* relativeInterpolatingPolynomialCoefficients, double TCAguess, int maxiter, double accuracy){
	/* Find the time of closest approach (TCA) between two objects given the polynomial coefficients that interpolate their positions in a dimensionless
	time interval [0,1]. This is done by finding the time t when the relative velocity squared is 0, i.e. at a local minimum of the relative distance.
	A newton-Rhapson method is employed in this purpose that will either pinpoint the root to within +-accuracy or until the maximum number of iterations, maxiter,
    have elapsed. A third-order interpolation is assumed.
	It is assumed that the interpolating coefficients will span an interval short enough for this to be a local minimum to also be global.
	@param relativeInterpolatingPolynomialCoefficients - RxN matrix, where R is a number of dimensions and N number of coefficients of the interpolating polynomial
		(its order + 1, i.e. 4). Each row should correspond to each dimension.
	@param TCAguess - argument value from which the Newton-Rhapson will be started, if no approximate TCA is known 0.5 is recommended.
	@param maxiter - maximum number of iterations that will be performed by the Newton-Rhapson root finder.
	@param accuracy - relative distance squared bound to within which the TCA will be found.
	@return - time of closest approach in a non-dimansional form. N.B. that it may be greater than the interpolation interval, i.e. may be less than 0.0 or greater than 1.0.
	*/

	/* Initialise or recompute the relative distance, velocity and acceleration functions squared coefficients.
	 * Relative distance magnitude is a sixth order polynomial if both trajectories are interpolated as cubics. Note that the relative coefficiencts have X, Y and Z comppoonents.
	 */
                //a*a
    m6 = relativeInterpolatingPolynomialCoefficients->at(0).at(3)*relativeInterpolatingPolynomialCoefficients->at(0).at(3) + relativeInterpolatingPolynomialCoefficients->at(1).at(3)*relativeInterpolatingPolynomialCoefficients->at(1).at(3) + relativeInterpolatingPolynomialCoefficients->at(2).at(3)*relativeInterpolatingPolynomialCoefficients->at(2).at(3);
                       //a                                                        b
    m5 = 2.0*( relativeInterpolatingPolynomialCoefficients->at(0).at(3)*relativeInterpolatingPolynomialCoefficients->at(0).at(2) + relativeInterpolatingPolynomialCoefficients->at(1).at(3)*relativeInterpolatingPolynomialCoefficients->at(1).at(2) + relativeInterpolatingPolynomialCoefficients->at(2).at(3)*relativeInterpolatingPolynomialCoefficients->at(2).at(2) );
                //2*a                            c                                                                                              b*b
    m4 = 2.0*relativeInterpolatingPolynomialCoefficients->at(0).at(3)*relativeInterpolatingPolynomialCoefficients->at(0).at(1) + relativeInterpolatingPolynomialCoefficients->at(0).at(2)*relativeInterpolatingPolynomialCoefficients->at(0).at(2) + 2.0*relativeInterpolatingPolynomialCoefficients->at(1).at(3)*relativeInterpolatingPolynomialCoefficients->at(1).at(1) + relativeInterpolatingPolynomialCoefficients->at(1).at(2)*relativeInterpolatingPolynomialCoefficients->at(1).at(2) + 2.0*relativeInterpolatingPolynomialCoefficients->at(2).at(3)*relativeInterpolatingPolynomialCoefficients->at(2).at(1) + relativeInterpolatingPolynomialCoefficients->at(2).at(2)*relativeInterpolatingPolynomialCoefficients->at(2).at(2);
               //2*b*c                                                                                                                        2*a*d
    m3 = 2.0*( relativeInterpolatingPolynomialCoefficients->at(0).at(2)*relativeInterpolatingPolynomialCoefficients->at(0).at(1) + relativeInterpolatingPolynomialCoefficients->at(0).at(3)*relativeInterpolatingPolynomialCoefficients->at(0).at(0) + relativeInterpolatingPolynomialCoefficients->at(1).at(2)*relativeInterpolatingPolynomialCoefficients->at(1).at(1) + relativeInterpolatingPolynomialCoefficients->at(1).at(3)*relativeInterpolatingPolynomialCoefficients->at(1).at(0) + relativeInterpolatingPolynomialCoefficients->at(2).at(2)*relativeInterpolatingPolynomialCoefficients->at(2).at(1) + relativeInterpolatingPolynomialCoefficients->at(2).at(3)*relativeInterpolatingPolynomialCoefficients->at(2).at(0) );
                 //2*b*d                                                                                                                   c*c
    m2 = 2.0*relativeInterpolatingPolynomialCoefficients->at(0).at(2)*relativeInterpolatingPolynomialCoefficients->at(0).at(0) + relativeInterpolatingPolynomialCoefficients->at(0).at(1)*relativeInterpolatingPolynomialCoefficients->at(0).at(1) +2.0*relativeInterpolatingPolynomialCoefficients->at(1).at(2)*relativeInterpolatingPolynomialCoefficients->at(1).at(0) + relativeInterpolatingPolynomialCoefficients->at(1).at(1)*relativeInterpolatingPolynomialCoefficients->at(1).at(1) +2.0*relativeInterpolatingPolynomialCoefficients->at(2).at(2)*relativeInterpolatingPolynomialCoefficients->at(2).at(0) + relativeInterpolatingPolynomialCoefficients->at(2).at(1)*relativeInterpolatingPolynomialCoefficients->at(2).at(1);
                 //2*c*d
    m1 = 2.0*( relativeInterpolatingPolynomialCoefficients->at(0).at(1)*relativeInterpolatingPolynomialCoefficients->at(0).at(0) + relativeInterpolatingPolynomialCoefficients->at(1).at(1)*relativeInterpolatingPolynomialCoefficients->at(1).at(0) + relativeInterpolatingPolynomialCoefficients->at(2).at(1)*relativeInterpolatingPolynomialCoefficients->at(2).at(0) );
                //d*d
    m0 = relativeInterpolatingPolynomialCoefficients->at(0).at(0)*relativeInterpolatingPolynomialCoefficients->at(0).at(0) + relativeInterpolatingPolynomialCoefficients->at(1).at(0)*relativeInterpolatingPolynomialCoefficients->at(1).at(0) + relativeInterpolatingPolynomialCoefficients->at(2).at(0)*relativeInterpolatingPolynomialCoefficients->at(2).at(0);

	/* Find the TCA using a Newton-Rhapson search. */
	double currentRoot = TCAguess; // Current best estimate of the root.
	double dT; // Current step in the timensionless time.
	for(int i=0;i<=maxiter;i++){
		dT = getInterpolatedRelativeRangeRateSquared(currentRoot)/getInterpolatedRelativeAccelerationSquared(currentRoot); // Current step in time.
		currentRoot = currentRoot - dT; // Update the root.
		if(fabs(dT)<=accuracy)
			return currentRoot; // Reached the desired accuracy, no need to iterate further.
	};
	return currentRoot; // Return the best estimate of the root.
};

double Simulation::getInterpolatedRelativeDistanceSquared(double s){
	/* Given a dimensionless epoch return the relative distance between objects based on the relative polynomials interpolating
	polynomials coefficients m6, m5, m4, m3, m2, m1 and m0 that are this class' attributes.
	@return square of the relative distance.
	*/
	return m6*std::pow(s,6.0) + m5*std::pow(s,5.0) + m4*std::pow(s,4.0) + m3*std::pow(s,3.0) + m2*std::pow(s,2.0) + m1*s + m0;
};

double Simulation::getInterpolatedRelativeRangeRateSquared(double s){
	/* Given a dimensionless epoch return the relative velocity between objects based on the relative polynomials interpolating
	polynomials coefficients m6, m5, m4, m3, m2, m1 and m0 that are this class' attributes.
	@return square of the relative velocity.
	*/
	return 6.0*m6*std::pow(s,5.0) + 5.0*m5*std::pow(s,4.0) + 4.0*m4*std::pow(s,3.0) + 3.0*m3*std::pow(s,2.0) + 2.0*m2*s + m1;
};

double Simulation::getInterpolatedRelativeAccelerationSquared(double s){
	/* Given a dimensionless epoch return the relative distance between objects based on the relative polynomials interpolating
	polynomials coefficients m6, m5, m4, m3, m2, m1 and m0 that are this class' attributes.
	@return square of the relative distance.
	*/
    return 30.0*m6*std::pow(s,4.0) + 20.0*m5*std::pow(s,3.0) + 12.0*m4*std::pow(s,2.0) + 6.0*m3*s + 2.0*m2;
};

void Simulation::findInterpolatingCoefficients(std::vector< std::vector<double> >* posesInterpPtr, std::vector< std::vector<double> >* velosInterpPtr, std::vector<double>* nodeEpochsPtr, std::vector< std::vector<double> >* interpCoeffsPtr){
	/* Use the function value and derivative information at N points (size of posersInterpPtr and velosInterpPtr) to fit a polynomial of order 2N-1 through those points.
	Can, in principle, work for any number of dimensions R but was only tested for 3, i.e. Cartesian position interpolation.
	@param posesInterpPtr - pointer to a vector of the size NxR (Nx3 for Cartesian position) containing function values (Cartesian positions at N points).
	@param velosInterpPtr - pointer to a vector of the size NxR (Nx3 for Cartesian velocities) containing function derivatives (Cartesian velocities at N points).
							The derivatives will be assumed to be in the same order as values.
	@param nodeEpochs - Julian day epochs corresponding to the function values and derivatives of size N.
	@param interpCoeffsPtr - pointer to an Rx(2N-1) vector that will be filled with interpolation polynomials' coefficient in every out of R dimensions.
	*/
	
	/* Get the new right-hand side values for interpolation (new positions and velocities in the current interval). */
	for(size_t i=0;i<(size_t)NO_INTERPOLATION_POINTS;i++){
		rhsX.at(i) = posesInterpPtr->at(i).at(0);
		rhsY.at(i) = posesInterpPtr->at(i).at(1);
		rhsZ.at(i) = posesInterpPtr->at(i).at(2);
		rhsX.at(i+NO_INTERPOLATION_POINTS) = velosInterpPtr->at(i).at(0) * (nodeEpochsPtr->back() - nodeEpochsPtr->at(0))*JDAY_IN_SECONDS; // Change velocity to km/JDAY and non-dimensionalise.
		rhsY.at(i+NO_INTERPOLATION_POINTS) = velosInterpPtr->at(i).at(1) * (nodeEpochsPtr->back() - nodeEpochsPtr->at(0))*JDAY_IN_SECONDS;
		rhsZ.at(i+NO_INTERPOLATION_POINTS) = velosInterpPtr->at(i).at(2) * (nodeEpochsPtr->back() - nodeEpochsPtr->at(0))*JDAY_IN_SECONDS;
	};
	
	/* Solve the system of equations to find the interpolating polynomials' coefficients.
	 * Use the pre-computed LU decomposition of the dimensionless time arguments for speed 
	 * as it will never change since it's always [0,1]. */
	EquationsSolving::luSubst(&LU, &luDecompositionIndices, &rhsX);
	EquationsSolving::luSubst(&LU, &luDecompositionIndices, &rhsY);
	EquationsSolving::luSubst(&LU, &luDecompositionIndices, &rhsZ);

	/* Combine the interpolating polynomials' coefficients into one 2-D vector for ease of handling. */
	for(size_t i=0;i<rhsX.size();i++){ // Go through all the coefficients.
		interpCoeffsPtr->at(0).at(i) = rhsX.at(i); // Each row of interpCoeffsPtr corresponds to a dimension.
		interpCoeffsPtr->at(1).at(i) = rhsY.at(i);
		interpCoeffsPtr->at(2).at(i) = rhsZ.at(i);
	}
	
};

bool Simulation::XYZpreFilter(std::vector<double>* posVecPtr, double* thresholdValuePtr){
	/* Returns true if the XYZ filter as given in Rodriguez et. all 2002 has been passed.
	@param posVecPtr - current relative position of two objects.
	@param thresholdValuePtr - threshold radius.
	@return - if any of the coordinates is below the threshold return true, if all of them are above the treshold return false.
	*/
	for(std::vector<int>::size_type  i=0; i<posVecPtr->size(); i++){
		if( posVecPtr->at(i)<*thresholdValuePtr){
			return true;
		};
	};
	return false; // All the relative position coordinates are above the threshold value.
};

void Simulation::createObjectsFromTLEfile(const char* fileName, double defaultRadius_RB, double defaultRadius_PL, double defaultRadius_DEB, double defaultRadius_Other, int COV, std::string covpath, double defaultVarRad, double defaultVarIn, double defaultVarCross, double bStarMultiplicationFactor, const char * radiusFileName){
	/* Read a file that contains multiple three line elements and create SpaceObjects from them. Return them to a std::map
	where the keys will be objects' NORAD IDs (SSCs).
	@param fileName - name of the file from which to read the TLE data.
	@param defaultR_RB - default object to be used for object that have no entry in the radiusFileName in metres and the 1st line of the three-line element contains R/B.
	@param defaultR_PL - default object to be used for object that have no entry in the radiusFileName in metres and the 1st line of the three-line element contains P/L.
	@param defaultR_DEB - default object to be used for object that have no entry in the radiusFileName in metres and the 1st line of the three-line element contains DEB.
	@param defaultR_Other - default object to be used for object that have no entry in the radiusFileName in metres and the 1st line of the three-line element does not contain R/B, P/L or DEB.
	@param COV - type of covariance to be used for the objects.
		1 - an uncertainty sphere 1:1:1 will be used and only maximum collision probability will be computed.
		2 - a number of three-line elements will be read for every object and covariance will be estimated for those using the method of V.P. Osweiler. It will be assumed that
		3LEs are in the files called SSC.txt where SSC if the catalogue number of every object and that the 3LEs are sorted oldest to last (the most recent one is last).
		3 - fixed covariance matrices based on the predicted performance of the European Space Surveillance System (20x100x20 m^2 in each objects RTC) will be used.
		4 - time-invariant 1x1x1 km^2 covariances in every object's RTC reference frame will be used to enable comparison to STK.
		5 - skip collision probability calculation and write -1 to the PC and PCMAX fields of all conjunctions.
		6 - use time-invariant covariance matrices defined by defaultvarRad, defaultVarIn, defaultVarCross.
	@param covpath - path from which the TLEs will be read to estimate when COV is set to 2, ignored otherwise.
	@param defaultVarRad, defaultVarIn, defaultVarCross - radial, in-track, cross-track variances to be used by all the objects. Only used when COV=6.
	@param bStarMultiplicationFactor - a factor by which the B* coefficient of a satellite will be multiplied by, used to simulate atmospheric density changes.
	@param radiusFileName - name of the file from which to read the hard body radii of objects.
	*/
	/* Save the data with which the population was created in this run. */
	bStarScalingFactor = bStarMultiplicationFactor;
	defaultR_RB=defaultRadius_RB; defaultR_PL=defaultRadius_PL; defaultR_DEB=defaultRadius_DEB; defaultR_Other=defaultRadius_Other;
	CovarianceType =  COV;
	defaultVarianceRad=defaultVarRad; defaultVarianceIn=defaultVarIn; defaultVarianceCross=defaultVarCross;

	std::ifstream TLEfileStream(fileName); // File streams that are necessary.
	std::ifstream radiusFileStream(radiusFileName);
	if(!radiusFileStream.is_open()){
		std::cerr<<"No external radius file provided to "<<TITLE<<"."<<std::endl;
	};
	
	/* Read the hard body radii to assign to objects when they are being created. */
	std::map<std::string, double>* hardBodyRadiiPtr = new std::map<std::string, double>();
	std::string radiusLine; // Currently read line with object's NORAD ID and hard body radius in m.

	int counterRadii = 0; // Current line of the TLE, 1 or 2.
	while( std::getline(radiusFileStream, radiusLine) ){
		std::istringstream radStringStream(radiusLine);
		if(counterRadii>5){ // Skip the header.
			int NORAD_ID; double radius; std::string NORAD_ID_STR;
			if (!(radStringStream >> NORAD_ID >> radius)) { break; } // Read the NORAD ID and the corresponding radius.
			
			std::stringstream strstream; // Convert the NORAD ID to a string.
			strstream << NORAD_ID; strstream>>NORAD_ID_STR;
			
			hardBodyRadiiPtr->insert ( std::pair<std::string,double>( NORAD_ID_STR, radius ) );
		}else{
			counterRadii+=1;
		};
	}

	/* Read the TLEs and create the SpaceObjects. */
	std::string TLEline; // Currently read TLE line.
	std::vector<std::string> currentTLE; currentTLE.resize(3); // Assembled TLE (second and third lines) and the object name (first line).

	int counterTLEs = 0; // Current line of the TLE, 0, 1 or 2.
	while( std::getline(TLEfileStream, TLEline) ){
		std::istringstream iss(TLEline);
		currentTLE.at(counterTLEs) = TLEline; // Add a new line.
		counterTLEs+=1;
		if(counterTLEs == 3){ // Read a whole TLE.
			counterTLEs = 0;
			try{
				SpaceObject tempSO = SpaceObject( currentTLE, COV, covpath, bStarMultiplicationFactor );
				objectsPtr->insert ( std::pair<std::string,SpaceObject>( tempSO.NORAD_ID, tempSO) );
				try{
					objectsPtr->at(tempSO.NORAD_ID).SetHardBodyRadius( hardBodyRadiiPtr->at(tempSO.NORAD_ID) ); // Assign the radius value.
				}catch( std::out_of_range ){ // Radius for this object does not exist in the database, use a default value for this object type.
					if(currentTLE.at(0).find("R/B") != std::string::npos){
						objectsPtr->at(tempSO.NORAD_ID).SetHardBodyRadius( defaultR_RB );
					}else if(currentTLE.at(0).find("P/L") != std::string::npos){
						objectsPtr->at(tempSO.NORAD_ID).SetHardBodyRadius( defaultR_PL );
					}else if(currentTLE.at(0).find("DEB") != std::string::npos){
						objectsPtr->at(tempSO.NORAD_ID).SetHardBodyRadius( defaultR_DEB );
					}else{
						objectsPtr->at(tempSO.NORAD_ID).SetHardBodyRadius( defaultR_Other );
					};
				};
				
				if(COV==6) // Are we supposed to set the position variances for all the objects?
				{ // Will set all the covariance matrices nicely and convert to km^2.
					objectsPtr->at(tempSO.NORAD_ID).SetPositionVariances(defaultVarRad, defaultVarIn, defaultVarCross);
				};
			}catch(FatalError){}; // Too few TLEs for tempObject exist to estimate the covariance, don't worry about that, it's planned.
		};
	}

	char buff[25]; time_t now = time (0);
    #ifdef _MSC_VER 
		struct tm timeinfo;
		localtime_s(&timeinfo, &now);
		strftime(buff, 25, "%Y-%m-%d %H:%M:%S.000", &timeinfo);
	#else
		strftime(buff, 25, "%Y-%m-%d %H:%M:%S.000", localtime (&now));
	#endif
	std::cout<<buff<<": Generated "<<objectsPtr->size()<<" objects."<<std::endl;
};

void Simulation::saveConjunctions( const char * outFileName, std::string version){
	/* Write all the conjunctions' data stored in conjunctionsFoundPtr to the specified file. Will overwrite the file if one with the same name exists in the output directory.
	Will also write a header explaining the file's collumn information.
	@param outFileName - name of the file which will be created and written to.
	@param version - version of the program that'll be written in the header file.
	*/
	std::ofstream outputFile; // Open the output file stream.
	outputFile.open( outFileName );	

	if (outputFile!=NULL){
		/* Write the header. */
		outputFile<<"Simulation: "<<TITLE<<std::endl;		
		outputFile<<"Version: "<<version<<std::endl;
		outputFile<<"Default radius for R/B (m): "<<defaultR_RB<<", Default radius for P/L (m): "<<defaultR_PL<<", Default radius for DEB (m): "<<defaultR_DEB<<", Default radius for others (m): "<<defaultR_Other<<std::endl;
		outputFile<<"Default radial variance (m^2): "<<defaultVarianceRad<<", Default in-track variance (m^2): "<<defaultVarianceIn<<", Default cross-track variance (m^2): "<<defaultVarianceCross<<std::endl;
		outputFile<<"B* scaling factor: "<<bStarScalingFactor<<std::endl;
		outputFile<<"Coarse time step (s): "<<COARSE_TIME_STEP_S<<std::endl;
		outputFile<<"Number of points used for interpolation: "<<NO_INTERPOLATION_POINTS<<std::endl;
		outputFile<<"Perigee/Apogee pre-filter pad (km): "<<PERIGEE_APOGEE_PAD<<std::endl;
		outputFile<<"Conjunction threshold distance (km): "<<CONJUNCTION_THRESHOLD_DISTANCE<<std::endl;
		outputFile<<"Type of covariance used: "<<CovarianceType<<std::endl;
		outputFile<<"Direct integration used for all Pc estimations: "<<EnforceDirectNumericalIntegration<<std::endl;
		outputFile<<"Execution time (s): "<<difftime(TIME_FINISHED,TIME_STARTED)<<std::endl;
		outputFile<<"Conjnctions found: "<<conjunctionsFoundPtr->size()<<std::endl;
		
		/* Write data about all the founc donjunctions. */
		outputFile<<"Primary SSC, Secondary SSC, TCA, Miss distance (km), Maximum collision probability, True collision probability, Relative velocity (km/s), Collision radius (km)"<<std::endl;
		for(std::vector<int>::size_type  i = 0; i<conjunctionsFoundPtr->size(); i++){
			outputFile << conjunctionsFoundPtr->at(i)<<std::endl;
		}
		outputFile.close(); // Close the file when done.
		std::cout<<"Saved the data for "<<conjunctionsFoundPtr->size()<<" conjunctions."<<std::endl;
	} else {
		std::cerr<<"Saving the conjunctions' data failed - file stream could not be opened."<<std::endl;
	}
};

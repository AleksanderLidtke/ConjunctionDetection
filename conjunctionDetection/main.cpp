#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "sgp4unit.h"
#include "sgp4ext.h"
#include "Simulation.h"
#include "SpaceObject.h"
#include "Conjunction.h"

#define SECOND_IN_JDAYS 1.157401129603386e-05; // 1 second expressed in Julian Days. JDAY/s
#define JDAY_IN_SECONDS 86400.46863810098; // 1 Julian Day expressed in seconds. s/JDAY

/* HELPER FUNCTIONS THAT HANDLE COMMAND LINE ARGUMENTS. */
char* getCmdOption(char ** begin, char ** end, const std::string & option){
	/* Get the value of the command line option that has been entered after the specified flag.
	@param begin, end - pointers to the beginning and the end of the char array of the command line arguments, i.e. argv, argv+argc, respectively.
	@param option - the flag that is to be found.
	@return - value corresponding to the given flag or 0 if nothing was supplied after the flag.
	*/
    char ** itr = std::find(begin, end, option);
    if(itr != end && ++itr != end){
        return *itr;
    };
    return 0;
};

bool cmdOptionExists(char** begin, char** end, const std::string& option){
	/* Check if a given command line option has been supplied.
	@param begin, end - pointers to the beginning and the end of the char array of the command line arguments, i.e. argv, argv+argc, respectively.
	@param option - the flag that is to be found.
	@return - true or false depending on whether the option has or has not been supplied to the program, respectively.
	*/
    return std::find(begin, end, option) != end;
};

/* MAIN PROGRAM THAT HANDLES INITIALISATION OF EVERYTHING AS WELL AS MANAGES CONJUNCTION DETECTION AND RECORDING. */
int main(int argc, char *argv[]){

	std::string VERSION = "1.8.0"; // Version of the entire program.
	/* Check whether to display help. Terminate the program afterwards. */
	if( cmdOptionExists(argv, argv+argc, "-h") || cmdOptionExists(argv, argv+argc, "--help") ){
		std::cout<<"Help for Conjunction Detection and Assessment v. "<<VERSION<<std::endl;
		std::cout<<"Argument flag / expected value type following the flag / default value / meaning"<<std::endl
			<<"-h || --help / NONE / NONE / Display this help information and terminate the program."<<std::endl
			<<"-v || --version / NONE / NONE / Display the current version of the program."<<std::endl
			<<"-u || --usage / NONE / NONE / Display information about how to use the program."<<std::endl
			<<"-t / const char* / cppRun / Change the title of the run that will be written to the output file."<<std::endl
			<<"-m / int / 1 / Mode of analysis: 1 for One-on-All (requires at least one primary to be supplied via the SSC argument and optionally also the PrimaryEphemeris file) or 2 for All-on-All."<<std::endl
			<<"-COV / int / 1 / Type of covariance to be used."<<std::endl
			<<"\t1 - a 1:1:1 uncertainty sphere will be assumed and only maximu collision probability will be computed."<<std::endl
			<<"\t2 - estimated covariance matrices will be estimated for every object using method of V.P. Osweiler 2006 (it is expected that the 3LEs will be written to files called SSC.txt where SSC is every object's SSC and the files will contain all the 3LEs with the most-recent one last)."<<std::endl
			<<"\t3 - time-invariant 20x100x20 m^2 covariances in every object's RTC reference frame will be used to simulate the desired performance of the European Space Surveillance System."<<std::endl
			<<"\t4 - time-invariant 1x1x1 km^2 covariances in every object's RTC reference frame will be used to enable comparison to STK."<<std::endl
			<<"\t5 - skip collision probability calculation and write -1 to the PC and PCMAX fields of all conjunctions."<<std::endl
			<<"\t6 - accept radial, in-track, cross-track variances in m^2 and skip maximum collision probability calculation for speed (write -1 to the PCMAX field of all conjunctions)."<<std::endl
			<<"-covpath / const char* / ./Covariance 3LEs / Directory from which to read the files with TLEs that will be used to estimate the covariance if COV is set to 2."<<std::endl
			<<"-EnforceDI / NONE / NONE / Whether to enforce direct numerical integration to be performed in all collision probability estimates. It will not be used in all cases, by default, unless object size is comparable to the miss distance (80% of it). K. Chan's series expansion will be used otherwise."<<std::endl
			<<"-TH / double / 20.0 / Threshold distance in km below which conjunctions will be recorded."<<std::endl
			<<"-coarsedT / double / 600.0 / Time step used to interpolate trajectories and find conjunctions."<<std::endl
			<<"-TLE / const char * / spaceTrack3LEs_23_10_2013.txt / File from which the three line elements will be read and used in the analysis."<<std::endl
			<<"-rRB / double / 1.7691 / Default radius for the objects marked as R/B (rocket bodies) in the three-line element file. Default value comes from MASTER 2009 population average for this type of object. In metres."<<std::endl
			<<"-rPL / double / 1.0350 / Default radius for the objects marked as P/L (payloads) in the three-line element file. Default value comes from MASTER 2009 population average for this type of object. In metres."<<std::endl
			<<"-rDEB / double / 0.1558 / Default radius for the objects marked as DEB (debris objects) in the three-line element file. Default value comes from MASTER 2009 population average for this type of object. In metres."<<std::endl
			<<"-rOther / double / 0.3470 / Default radius for the objects not marked as R/B, P/L or DEB in the three-line element file. Default value comes from MASTER 2009 population average for this type of object. In metres."<<std::endl
			<<"-bstar / double / 1.0 / Factor by which to multiply all objects B* factors to simulate variations in solar activity."<<std::endl
			<<"-varRad / double / 40.0 / Position variance in m^2 set for all the objects in the radial direction."<<std::endl
			<<"-varIn / double / 200.0 / Position variance in m^2 set for all the objects in the in-track direction."<<std::endl
			<<"-varCross / double / 100.0 / Position variance in m^2 set for all the objects in the cross-track direction."<<std::endl
			<<"-SSC / char* / NONE / SSC ID of the primary for which the conjunctions will be found. Ignored in All-on-All mode. When an ephemeris parameter argument is given it should be specified in order to ignore the primary from the group of secondary objects."<<std::endl
			<<"-PrimaryEphemeris / char* / NONE / STK v 10 ephemeris file that contains the ephemeris of the primary satellite that will be screened against conjunctions against the objects from TLE file. Should contain Julian day epoch, Cartesian positions and velocities in meteres and metres per second in TEME reference frame. The Julian Days override the coarseDt, jdayStart and jdayStop arguments."<<std::endl
			<<"-PrimaryRadius / double / 5.0 / Radius of the primary object in meteres that will be used for collision probability calculations, only used when using an external ephemeris file."<<std::endl
			<<"-o / const char* / cppRunOutFileName / Name of the output file that will contain information about the found conjunctions."<<std::endl
			<<"-wgs / integer 72, 721, or 84 / 84 / Change the WGS Earth ellipsoid type used."<<std::endl
			<<"-jdayStart / double / 23-10-2013 3:58:20.413924 / Julian day from which to start the analysis."<<std::endl
			<<"-jdayStop / double / 24-10-2013 3:58:20.413924 / Julian day at which to stop the analysis."<<std::endl;
		return 0; // Don't execute the program after displaying the help.
	};

	/* Display the usage information only. */
	if( cmdOptionExists(argv, argv+argc, "-u") || cmdOptionExists(argv, argv+argc, "--usage") ){
		std::cout<<"Conjunction Detection and Assessment v. "<<VERSION<<std::endl
			<<"The basic use cases intended and supported by this software are:"<<std::endl
			<<"1) Perform one-on-all conjunction screening (mode=1) whereby one primary object (either a TLE or an ephemeris table) will be screened for conjunctions against a supplied list of TLEs for secondaries. Various types of covariances may be assumed in the process:"<<std::endl
			<<"\t- only maximum probabilities may be computed when assuming spherical shape of the covariance ellipsoids (COV=1)"<<std::endl
			<<"\t- covariance of each individual TLE may be estimated based on previous TLEs for every object (COV=2). Then true probability using the covariance propagated to the epoch of the closest approach and maximum probability using the true aspect ratio of the covariance ellipsoid will be computed. Note that this covariance estimation cannot be performed for the ephemeris tables."<<std::endl
			<<"\t- covariance matrices that mimic the behaviour of the European Space Surveillance System (COV=3) can be set for all objects (40x200x100 m^2 in individual RTC reference frames)"<<std::endl
			<<"\t- arbitrary time-invariant covariance matrices may be used (COV=6) for all objects in individual RTC reference  frames using the varRad, varIn, varCross arguments."<<std::endl
			<<"\t- only the miss distances for all the conjunctions can be found (COV=5)"
			<<"2) Perform all-on-all conjunction analyses (mode=2) where conjunctions between all the supplied TLEs will be found. Use of ephemeris tables for all-on-all assessments is not supported. The supported covariance modes are:"<<std::endl
			<<"\t- only maximum probabilities may be computed when assuming spherical shape of the covariance ellipsoids (COV=1)"<<std::endl
			<<"\t- covariance of each individual TLE may be estimated based on previous TLEs for every object (COV=2). Then true probability using the covariance propagated to the epoch of the closest approach and maximum probability using the true aspect ratio of the covariance ellipsoid will be computed."<<std::endl
			<<"\t- covariance matrices that mimic the behaviour of the European Space Surveillance System (COV=3) can be set for all objects (40x200x100 m^2 in individual RTC reference frames)"<<std::endl
			<<"\t- arbitrary time-invariant covariance matrices may be used (COV=6) for all objects in individual RTC reference  frames using the varRad, varIn, varCross arguments."<<std::endl
			<<"\t- only the miss distances for all the conjunctions can be found (COV=5)"<<std::endl;
		return 0;
	}

	/* Display version only. */
	if( cmdOptionExists(argv, argv+argc, "-v") || cmdOptionExists(argv, argv+argc, "--version") ){
		std::cout<<"Conjunction Detection and Assessment v. "<<VERSION<<std::endl;
		return 0;
	}

	/* Go through the command line arguments and set the supplied values or the default ones. */
	std::string TITLE; // Title.
	if( cmdOptionExists(argv, argv+argc, "-t") ){
		std::basic_istringstream<char> titleSS( getCmdOption(argv, argv + argc, "-t") );
		titleSS >> TITLE;
	}else{
		TITLE = "cppRun";
	};

	int wgs; // WGS ellipsoid.
	if( cmdOptionExists(argv, argv+argc, "-wgs") ){
		std::basic_istringstream<char> wgsSS( getCmdOption(argv, argv + argc, "-wgs") );
		wgsSS >> wgs;
	}else{
		wgs=84;
	};

	int analysisMode; // One on all or all on all.
	if( cmdOptionExists(argv, argv+argc, "-m") ){
		std::basic_istringstream<char> modeSS( getCmdOption(argv, argv + argc, "-m") );
		modeSS >> analysisMode;
	}else{
		analysisMode=1;
	};

	int COV; // Type of covariance to be used.
	if( cmdOptionExists(argv, argv+argc, "-COV") ){
		std::basic_istringstream<char> covSS( getCmdOption(argv, argv + argc, "-COV") );
		covSS >> COV; // 1 for spherical ellipsoid, 2 for estimated.
	}else{
		COV = 1; // Use a sphere with maximum probability by default.
	};

	std::string COVPATH; // Where to take the covariance TLEs from.
	if( cmdOptionExists(argv, argv+argc, "-covpath") ){
		std::basic_istringstream<char> covpathSS( getCmdOption(argv, argv + argc, "-covpath") );
		covpathSS >> COVPATH;
		std::cout<<"Using "<<COVPATH<<" folder to read 3LEs from in order to estimate the covariance."<<std::endl;
	}else{
		COVPATH = "./Covariance 3LEs";
	};

	double threshold; // Threshold distance in km.
	if( cmdOptionExists(argv, argv+argc, "-TH") ){
		std::basic_istringstream<char> threshSS( getCmdOption(argv, argv + argc, "-TH") );
		threshSS >> threshold;
	}else{
		threshold=20.0;
	};

	double coarseDT; // Coarse time step.
	if( cmdOptionExists(argv, argv+argc, "-coarsedT") ){
		std::basic_istringstream<char> coarsedTSS( getCmdOption(argv, argv + argc, "-coarsedT") );
		coarsedTSS >> coarseDT;
	}else{
		coarseDT=600.0;
	};

	const char* TLEfile; // Objects' file.
	if( cmdOptionExists(argv, argv+argc, "-TLE") ){
		TLEfile = getCmdOption(argv, argv + argc, "-TLE");
		std::cout<<"Using object TLEs from "<<TLEfile<<std::endl;
	}else{
		TLEfile = "spaceTrack3LEs_23_10_2013.txt";
	};

	double rRB, rPL, rDEB, rOther; // Default object radii in metres.
	if( cmdOptionExists(argv, argv+argc, "-rRB") ){
		std::basic_istringstream<char> radRB_SS( getCmdOption(argv, argv + argc, "-rRB") );
		radRB_SS >> rRB;
	}else{
		rRB = 1.7691;
	};
	if( cmdOptionExists(argv, argv+argc, "-rPL") ){
		std::basic_istringstream<char> radPL_SS( getCmdOption(argv, argv + argc, "-rPL") );
		radPL_SS >> rPL;
	}else{
		rPL = 1.0350;
	};
	if( cmdOptionExists(argv, argv+argc, "-rDEB") ){
		std::basic_istringstream<char> radDEB_SS( getCmdOption(argv, argv + argc, "-rDEB") );
		radDEB_SS >> rDEB;
	}else{
		rDEB = 0.1558;
	};
	if( cmdOptionExists(argv, argv+argc, "-rOther") ){
		std::basic_istringstream<char> radOther_SS( getCmdOption(argv, argv + argc, "-rOther") );
		radOther_SS >> rOther;
	}else{
		rOther = 0.3470;
	};

	double bStarMultiplicationFactor; // B*
	if( cmdOptionExists(argv, argv+argc, "-bstar") ){
		std::basic_istringstream<char> bStarSS( getCmdOption(argv, argv + argc, "-bstar") );
		bStarSS >> bStarMultiplicationFactor;
	}else{
		bStarMultiplicationFactor = 1.0;
	};
	
	double varRad, varIn, varCross; // RIC variances to use for all the objects in m^2.
	if( cmdOptionExists(argv, argv+argc, "-varRad") ){
		std::basic_istringstream<char> varRad_SS( getCmdOption(argv, argv + argc, "-varRad") );
		varRad_SS >> varRad;
	}else{
		varRad = 20.0;
	};
	if( cmdOptionExists(argv, argv+argc, "-varIn") ){
		std::basic_istringstream<char> varIn_SS( getCmdOption(argv, argv + argc, "-varIn") );
		varIn_SS >> varIn;
	}else{
		varIn = 100.0;
	};
	if( cmdOptionExists(argv, argv+argc, "-varCross") ){
		std::basic_istringstream<char> varCross_SS( getCmdOption(argv, argv + argc, "-varCross") );
		varCross_SS >> varCross;
	}else{
		varCross = 20.0;
	};
	
	const char* primarySSC; // Primary objects.
	if( cmdOptionExists(argv, argv+argc, "-SSC") ){
		primarySSC = getCmdOption(argv, argv + argc, "-SSC");
	}else{
		primarySSC = "";
	};

	std::string primaryEphemerisFile; // Ephemeris table for the primary object.
	if( cmdOptionExists(argv, argv+argc, "-PrimaryEphemeris") ){
		primaryEphemerisFile = getCmdOption(argv, argv + argc, "-PrimaryEphemeris");
	}else{
		primaryEphemerisFile = std::string();
	};

	
	double primaryRadius; // Radius for the primary object to be used when using an external ephemeris file.
	if( cmdOptionExists(argv, argv+argc, "-PrimaryRadius") ){
		std::basic_istringstream<char> primaryRadiusSS(  getCmdOption(argv, argv + argc, "-PrimaryRadius") );
		primaryRadiusSS >> primaryRadius;
	}else{
		primaryRadius = 5.0;
	};

	const char* outFileName; // Output file name.
	if( cmdOptionExists(argv, argv+argc, "-o") ){
		outFileName = getCmdOption(argv, argv + argc, "-o");
	}else{
		outFileName = "cppRunOutFileName";
	};

	double ANALYSIS_INTERVAL_START; // Analysis start
	if( cmdOptionExists(argv, argv+argc, "-jdayStart") ){
		std::basic_istringstream<char> startSS( getCmdOption(argv, argv + argc, "-jdayStart") );
		startSS >> ANALYSIS_INTERVAL_START;
	}else{
		jday(2013,10,23,3,58,20.413924, ANALYSIS_INTERVAL_START);
	};

	double ANALYSIS_INTERVAL_STOP; // Analysis stop
	if( cmdOptionExists(argv, argv+argc, "-jdayStop") ){
		std::basic_istringstream<char> stopSS( getCmdOption(argv, argv + argc, "-jdayStop") );
		stopSS >> ANALYSIS_INTERVAL_STOP;
	}else{
		jday(2013,10,24,3,58,20.413924, ANALYSIS_INTERVAL_STOP);
	};

	/* Create the simulation and start it. */
	Simulation SimulationToRun = Simulation(TITLE, wgs, threshold, coarseDT);

	if( cmdOptionExists(argv, argv+argc, "-EnforceDI") ){ // User has requested direct numerical integration to be used at all times.
		SimulationToRun.setDirectNumericalIntegration( true );
	};

	SimulationToRun.createObjectsFromTLEfile(TLEfile, rRB, rPL, rDEB, rOther, COV, COVPATH, varRad, varIn, varCross, bStarMultiplicationFactor);

	//char* envisat_NORADID="27386"; char* zenit2_NORADID="27006"; char* strand1_NORADID="39090"; char* delta1RB_NORADID="862";
	//char* case1p="3669"; char* case2p = "7005";
	std::vector<std::string> primaries = std::vector<std::string>(); // List of primaries to be assessed against the catalogue.
	primaries.push_back(primarySSC); // This could be a list of primaries but this wouldn't work with many ephemeris tables and there are no input argument options to pass several primaries' SSCs.

	if( primaryEphemerisFile.empty() ){ // Find conjunctions between SGP4-generated ephemerides.
		SimulationToRun.performAnalysis(ANALYSIS_INTERVAL_START, ANALYSIS_INTERVAL_STOP, analysisMode, primaries); // Find the conjunctions.
		//SimulationToRun.performAnalysis(2013,10,23,3,58,20.413924,2014,10,23,3,58,20.413924, 2, primaries); // An alternative overloaded method.
	}else if( analysisMode==1 ){ // Here use the ephemeris table to get the state vectors of the primary.
		SimulationToRun.performAnalysis(primarySSC, primaryEphemerisFile, primaryRadius);
	}else if( (primaryEphemerisFile.empty()==false) && (analysisMode==2) ){ // User wants to run an all-on-all simulation but has specified an ephemeris table for the primary. This is meaningless.
		std::string msg = "Specified an all-on-all analysis mode AND an ephemeris table for the primary, which is ambiguous. Refine the input arguments list.";
		throw FatalError(msg,__FILE__,__LINE__);
	}
	
	SimulationToRun.saveConjunctions(outFileName, VERSION); // Save the results.

	std::cout<<"Found "<<SimulationToRun.conjunctionsFoundPtr->size()<<" conjunctions."<<std::endl;
	return 0;
}

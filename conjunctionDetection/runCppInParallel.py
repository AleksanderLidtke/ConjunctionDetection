#!/local/software/python/2.7.5/bin/python
# -*- coding: utf-8 -*-
"""
Run conjunction detection an assessment code in C++ by using Python to parallelise
it into several processes.

Created on Sun May  8 16:07:42 2016

@author: alek
@version: 1.0.0
@since: Sun May  8 16:07:42 2016

CHANGELOG:
Sun May  8 16:07:42 2016 - 1.0.0 - alek - Issued the first draft version.
"""

import sys, multiprocessing, subprocess, numpy, copy, getopt
from mpi4py import MPI

def runCpp(startStopEpochsAndID):
    """ Run the conjunction detection and assessment code with the arguments
    specified in the global namespace between the chosen start and end epochs.
    Write STDOUT and STDERR to log.logFileName with the run ID.
    
    Arguments
    ----------
    startStopEpochsAndID - numpy.ndarray of shape (1,3) containing floats with Julian Day
        epochs that mark beginnning and end epochs of the simulation to be run and a unique ID
        to record in output and title, in the respective entries.
    
    Returns
    ----------
    int output code of the C++ code or -1 if it hasn't finished.
    
    Global variables used
    ----------
    * arguments - list of arguments for the Cpp code
    * executable - path to the executable to use
    * logFileName - name of the file which to append with the run ID and 
      prepend with log. and where to write the STDERR and STDOUT logs.
    """
    args = copy.deepcopy(arguments) # Different processes will modify this. Make sure we don't mess it up.
    args[args.index("-t")+1]=args[args.index("-t")+1]+str(int(startStopEpochsAndID[2])) # Add simulation ID here. Cast it to an int (it's a float for interface reasons).
    args[args.index("-o")+1]=args[args.index("-o")+1]+str(int(startStopEpochsAndID[2]))
    args[args.index("-jdayStart")+1]=str(startStopEpochsAndID[0])
    args[args.index("-jdayStop")+1]=str(startStopEpochsAndID[1])
    
    ret=-1 # Output status.
    with open("log.{}_{}".format(logFileName,int(startStopEpochsAndID[2])),"w") as logF:
        proc = subprocess.Popen(args,executable=executable, stdout=logF, stderr=logF) # Dump stderr and stdout into the same file - there'll be many of them anyway.
        ret=proc.wait() # Wait for C++ to finish, particularly close the output files.
    
    return ret

def printHelp():
	""" Print the help message for the programme. """
	print """Use:
		        * --tle to specify the three line element file,
		        * -i to specify which ephemeris accuracy from the sampling plan to us,
		        * -t to set the simulation title,
		        * -o to set the output name,
		        * -l to set the log file name,
		        * -s how many separate simulations to split the desired time interval into, ideally no. nodes X no. cores per node,
		        * -p how many processes per node to use for the simulations,
		        * --jdayStart for the start Julian Day epoch,
		        * --jdayStop for the end Julian Day epoch,
				* --plan= to specify the sampling plan from which to take ephemeris uncertainties."""

if __name__=="__main__":
    " Start the MPI connection. "
    comm=MPI.COMM_WORLD # Connect to the MPI communicator.
    size=comm.Get_size() # Number of processes in the commmunicator (no. MPI parallel processes in this job).
    rank=comm.Get_rank() # The rank of this process in the communicator.
    name=MPI.Get_processor_name() # Name of the processor running this process, e.g. purple0009.
    print "Hello from {} with rank {}".format(name,rank) # Good to know that the nodes have started.

    " Default simulation settings. "
    logFileName="testLog" # Will write individual STDOUT and STDERR to log.logFileName_# where # is the integer ID of the run.
    outFileName="testOutFile" # Each parallel run will dump data to this file, appended with the simulation ID.
    simTitle="testTitle" # Title for all the simulations.
    noSims=4 # Number of separate simulations we'll run.
    noProcesses=4 # Number of processes we'll use.
    simulationIDs=range(noSims) # Integer IDs of every simulation, will be written to log files and output files.
    samplingPlanIDX=0 # Which point from the sampling plan to run the simulation for.
    samplingPlanName="JCASamplingPlan_mTH1-4585_PCTH_sR_sI_sC_50pts.npy" # Which sampling plan to use.

    tleSnapshot="spaceTrack3LEs_07_11_2013.txt" # Three line element snapshot to use.
    
    jdStart=2456603.5 # Epoch bounds between which conjunctions will be found.
    jdStop=2456603.5185185187
    # Dafaults are 07 Nov 2013 0:0:0.0 to 07 Nov 2013 00:26:40.000
    # With two simulations that gives a mid-point epoch of 07 Nov 2013 00:13:20.000
    
    executable="/home/al11g09/ConjunctionDetectionCpp/conjunctionDetection_PA_XYZ_R2_MD_fR2_CoarseSearch" # Cpp conjunction detecction executable to use.
    
    " Get commend line arguments from the bash script, if any. "
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:e:o:l:s:p:",["help","tle=","jdayStart=","jdayStop=","plan="]) # Skip the name of the script, i.e. the first argv
    except getopt.GetoptError:
        print "GetOpt error."
        printHelp()
        sys.exit(-1) # We're done, no idea what the user wants.
    
    for opt, arg in opts: # Go through the options and set the settings accordingly.
        if opt in ("-h", "--help"):
			printHelp()
			sys.exit(0) # Done
        elif opt in ("--tle"):
            tleSnapshot=arg
        elif opt in ("-i"):
            samplingPlanIDX=int(arg)
        elif opt in ("-e"):
            simTitle=arg
        elif opt in ("-o"):
            outFileName=arg
        elif opt in ("-l"):
            logFileName=arg
        elif opt in ("-s"):
            noSims=int(arg)
            simulationIDs=range(noSims) 
        elif opt in ("-p"):
            noProcesses=int(arg)
        elif opt in ("--jdayStart "):
            jdStart=float(arg)
        elif opt in ("--jdayStop"):
            jdStop=float(arg)
        elif opt in ("--plan"):
            samplingPlanName=arg

    " Read the chosen sampling plan to get the position uncertainty to simulate. "
    samplingPlan = numpy.load(samplingPlanName)
    radialVariances = numpy.power(samplingPlan[:,2],2) # RIC variances in m^2, as expected by C++.
    inTrackVariances = numpy.power(samplingPlan[:,3],2)
    crossTrackVariances = numpy.power(samplingPlan[:,4],2)

    " Pack the arguments of this script into arguments for C++. "
    # Start and end Julian Day epochs of separate simulations. C++ will handle
    # making sure that the conjunctions from overlapping analysis intervals
    # are not found more than once.
    epochs = numpy.linspace(jdStart,jdStop,noSims+1)
    
    # List of 3-tuples, each tuple contains the start end epoch of the given simulation, respectively.
    simulationEpochBounds=zip(epochs[:-1],epochs[1:],simulationIDs) # The only varying settings for all the simulations - start and end epochs, and corresponding IDs.
    simulationEpochBounds=numpy.array(simulationEpochBounds) # We'll be using arrays of indices on this, so it also has to be an array.
    
    # Arguments to supply to the Cpp conjunction detection and assessment code.
    arguments = ["-EnforceDI", # With small uncertainties we have to use precise integration to get the right PCs.
            "-m","2", # All on all.
            "-TH","1", # With accurate ephemerides we want only the really close conjunctions, rest is trash.
            "-COV","6", # A mode that accepts the position variances as arguments.
            "-TLE",tleSnapshot, # Use this TLE file - successfully used it many times before.
            "-rRB", "5.0",#"1.7691", # Default radii of the objects in metres.
            "-rPL", "5.0",#"1.0350",
            "-rDEB", "5.0",#"0.1558",
            "-rOther", "5.0",#"0.3470",
            # Use coarsedT less than 400 s - now the smallest differences in MD will play a role for the Pc.
            "-coarsedT", "100", # The default value in the program isn't the one I decided to use in the thesis.
            "-t", simTitle, # From here on the arguments will be changed by the runCpp function.
            "-o", outFileName,
            "-jdayStart", None,
            "-jdayStop", None,
            "-varRad", str(radialVariances[samplingPlanIDX]), # Variances in m^2.
            "-varIn", str(inTrackVariances[samplingPlanIDX]),
            "-varCross", str(crossTrackVariances[samplingPlanIDX])
                 ]
    
    " Check using MPI which simulationEpochBounds to run on this node. "
    simsPerNode=noSims/size # We want this many simulations per node. If it isn't equal to the no. nodex X no. cores so be it - the last node will do extra work.

    if rank==0: # We're on the lead node. Let it decide what the other nodes will do, i.e. which simulations to run.
        for i in range(1,size): # Send the indices of simsulations to run to other nodes.
            if i!=size-1: # Not the last node.
                simsOnNode=numpy.arange(i*simsPerNode,(i+1)*simsPerNode,1)
            else: # This one has to do extra work if there are too many simulations requested.
                simsOnNode=numpy.arange(i*simsPerNode,noSims,1)
            comm.Send([simsOnNode, MPI.INT], dest=i, tag=77) # Tell every non-lead node which simulations to run.
        simsOnNode=numpy.arange(0,simsPerNode,1) # The lead node will run these simulations itself.
    elif rank==size-1: # Last note - it might have to run more simulations.
        simsOnNode=numpy.empty(simsPerNode+(noSims-size*simsPerNode), dtype=numpy.int64) # This node might have to run some extra simulations - it'll get more indices than the other nodes.
        comm.Recv([simsOnNode, MPI.INT], source=0, tag=77)
    else: # Neither the lead or last note - it's told what simulations to run and we know how many.
        simsOnNode=numpy.empty(simsPerNode, dtype=numpy.int64)
        comm.Recv([simsOnNode, MPI.INT], source=0, tag=77)

    simulationEpochBounds=simulationEpochBounds[simsOnNode] # After all, this node will run these simulations.
    print "===============\nNode {} size {}, rank {}, args:\n{}\nsimsOnNode: {}\n===============".format(name,size,rank,sys.argv,simsOnNode)

    " Run the simulations in parallel. "
    processPool = multiprocessing.Pool(noProcesses)
    rets = processPool.map(runCpp,simulationEpochBounds) # Run all the simulation in parallel.
    processPool.close() # Free up the memory.
    processPool.join() # Make sure these are done before we try to access their outputs.
    #TODO this node is done. Could get some of the remaining simulations from the nodes that are still running to cut down the total wall time.
    
    " Collect the data for this series of runs into a single output file. "
    wroteHeader=False # Have we added the header tot he main file yet?
    with open(outFileName,"w") as mainOutF: # Collect all the data here.
        for intermediateFName in [outFileName+str(i) for i in simulationIDs]: # Output files of all the batch runs.
            with open(intermediateFName,"r") as intFile:
                intLines=intFile.readlines()
                # First row with the actual conjunction data.
                try:
                    dataStartIDX=[i  for i in range(len(intLines)) if intLines[i].startswith("Primary SSC")][0]+1
                    if wroteHeader: # Don't worry about the header, it's tehre already; go for the data straight away.
                        mainOutF.writelines(intLines[dataStartIDX:]) # Write the output lines from this file to the main one.
                    else: # Write the header and the data.
                        mainOutF.writelines(intLines)
                        wroteHeader=True
                except IndexError: # This shouldn't happen, but you never know.
                    print "No conjunctions in file {}.".format(intermediateFName)

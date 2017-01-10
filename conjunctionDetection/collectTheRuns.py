#!/local/software/python/2.7.5/bin/python
# -*- coding: utf-8 -*-
"""
Collect the output of C++ conjunction detection and assessment code from several runs into a single
output file.

Created on Mon May 23 08:57:42 2016

@author: alek
@version: 1.0.0
@since: Mon May 23 08:57:42 2016

CHANGELOG:
Mon May 23 08:57:42 2016 - 1.0.0 - alek - Issued the first draft version.
"""
import sys, os, numpy, copy, getopt

if __name__=="__main__":
	" Default values for the arguments. "
	dirName=os.getcwd()
	baseName="testName{}"
	noRuns=30
	outFileName="testOutFName"
    
	" Get commend line arguments from the bash script, if any. "
	try:
		opts, args = getopt.getopt(sys.argv[1:],"hi:d:n:i:o:",["help"]) # Skip the name of the script, i.e. the first argv
	except getopt.GetoptError:
		print "GetOpt error."
		sys.exit(-1) # We're done, no idea what the user wants.
    
	for opt, arg in opts: # Go through the options and set the settings accordingly.
		if opt in ("-h", "--help"):
			print """Use:
            * -d to specify the directory where the runs are located,
            * -n to specify the base file name of every run,
            * -i to specify how many runs are present in the directory,
            * -o to set the name of the output file that will contain all the information from individual runs."""
			sys.exit(0) # We're done, we've printed the help message.
		elif opt in ("-d"):
			dirName=arg
		elif opt in ("-n"):
			baseName=arg
		elif opt in ("-i"):
			noRuns=int(arg)
		elif opt in ("-o"):
			outFileName=arg

	runFNames=[os.path.join(dirName,baseName.format(i)) for i in range(noRuns)] # We're to collect data from these files.

	" Collect the data for this series of runs into a single output file. "
	wroteHeader=False # Have we added the header tot he main file yet?
	with open(outFileName,"w") as mainOutF: # Collect all the data here.
		for runFName in runFNames:
			with open(runFName,"r") as intFile:
				intLines=intFile.readlines()
                # First row with the actual conjunction data.
				try:
					dataStartIDX=[i  for i in range(len(intLines)) if intLines[i].startswith("Primary SSC")][0]+1
					if wroteHeader: # Don't worry about the header, it's there already; go for the data straight away.
						mainOutF.writelines(intLines[dataStartIDX:]) # Write the output lines from this file to the main one.
					else: # Write the header and the data.
						mainOutF.writelines(intLines)
						wroteHeader=True
				except IndexError: # This shouldn't happen, but you never know.
					print "No conjunctions in file {}.".format(intermediateFName)

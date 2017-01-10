# -*- coding: utf-8 -*-
"""
Parse and analyse orbital conjunctions' data.
Created on Wed Apr 23 15:12:47 2014
TO be used to analyse one on all data for Envisat and compare those to STK CAT results.

@author: Aleksander Lidtke
@version 1.0.0
@since 22/05/2014 15:42:00

CHANGELOG:

"""
import datetime, pylab, scipy.integrate, matplotlib.pyplot, matplotlib, numpy, csv, math, utilities, sys

def computeMaximumCollisionProbabilityRelative(relativePosition, collisionRadius):
    """ Compute the maximum probability of collision between the two objects given their relative position within the 
    collision plane (must be in the same units) at the time of closest approach. The maximum probability as given in
    N. Berend (1999) will be computed assuming a spherical uncertainty distribution of both objects and further numerically 
    optimised to ensure that the true maximum probaility is computed. This has been found to
    accurately estimate the collision risk even though it doesn't account for velocity uncertainty. Assuming a spherical
    uncertainty distribution ensures that no encounter geometries yield a higher risk than others.\n
    @param relativePosition - relative position of objects 1 and 2 within the collision plane in the form of 1x2 numpy arrays.\n
    @param collisionRadius - combined radius that contains extremities of the objects in the same units as position.\n
    @return (probability, error) - 2-tuple containig the probability of collision and estimated numerical integration error."""
    
    """ FIND THE COMBINED COVARIANCE MATRIX TO YIELD WORST-CASE COLLISION PROBABILITY """
    covarianceScalingFactor = numpy.linalg.norm(relativePosition)/math.sqrt(2) # Standard deviation which yields the maximum collision probability for spherical uncertainty distribution.
    projectedCovarianceMatrixInPlane = numpy.array([[1.,0.],[0.,1.]])
    Pd = covarianceScalingFactor*covarianceScalingFactor*projectedCovarianceMatrixInPlane # The worst-case covariance matrix in the collision plane.     
    PdInverse = numpy.linalg.inv(Pd) # And its inverse.
                
    """ COMPUTE THE COLLISION PROBABILITY """
    def integrand(y, x):
        """ Returns the probability of a successful attempt as defined by the probability density function
        in N. Berend (1999). Based upon covariance matrix Pd. """
        discrepancyVector = numpy.array([[x], [y]]) - relativePosition.reshape(-1,1) # Instantenous [x y] minus the mean values.
        expTerm = math.exp( -0.5*float( discrepancyVector.T.dot(PdInverse).dot(discrepancyVector) ) ) # This will always be a length-1 array.
        return expTerm/( 2*math.pi*math.sqrt( numpy.linalg.det( Pd ) ) ) # A rather lengthy expression so divide it into parts.
        
    # Numerically integrate the PDF
    probability, error = scipy.integrate.dblquad(integrand, -collisionRadius, collisionRadius, lambda x: -math.sqrt(collisionRadius**2-x**2), lambda x: math.sqrt(collisionRadius**2-x**2), epsabs=1E-40)
    
    def probabilityFunction(scalingFactor):
        Pd = scalingFactor*scalingFactor*projectedCovarianceMatrixInPlane # Covariance matrix to yield the worst-case collision probability...
        PdInverse = numpy.linalg.inv(Pd) # ...and its inverse.
        def integrand(y, x):
            """ Returns the probability of a successful attempt as defined by the probability density function
            in N. Berend (1999). Based upon covariance matrix Pd. """
            discrepancyVector = numpy.array([[x], [y]]) - relativePosition.reshape(-1,1) # Instantenous [x y] minus the mean values.
            expTerm = math.exp( -0.5*float( discrepancyVector.T.dot(PdInverse).dot(discrepancyVector) ) ) # This will always be a length-1 array.
            return expTerm/( 2*math.pi*math.sqrt( numpy.linalg.det( Pd ) ) ) # A rather lengthy expression so divide it into parts.
        
        probability, error = scipy.integrate.dblquad(integrand, -collisionRadius, collisionRadius, lambda x: -math.sqrt(collisionRadius**2-x**2), lambda x: math.sqrt(collisionRadius**2-x**2), epsabs=1E-40)
        return -probability # Actually we're looking for a maximum.
    
    actualCovarianceScalingFactor = scipy.optimize.fmin(probabilityFunction, covarianceScalingFactor) # Actual factor that gives the true local maximum of the collison probability.
    probability = -probabilityFunction(actualCovarianceScalingFactor) # - because the optimiser has to minimise a function so it returns - probaility.
        
    if float(error)/float(probability) >= 0.1: # 10% relative error, i.e plenty.
        utilities.customPrint('Relative error when computing the collision probability exceeded 10% with actual value of: '+str(float(error)/float(probability)), utilities.RED)

    return probability, error
    
class Conjunction(object):
    __doc__="""A class that is used to stroe data about a particular orbital conjunction."""

    def __init__(self, conjunctionRecord):
        """ A class that stores information about a particular orbital conjunction.\n
        Parameters
        ----------
        conjunctionRecord : string
            string (with the trailing end line character) that contains comma separated:
            SSC1, SSC2, TCA, Miss Distance (km), Pc, relative V at TCA (km/s), collision radius (km)
        """
        temp = conjunctionRecord.strip("\n").split(",") # List of strings with all the data about the particular conjunction.
        
        if len(temp[0])==4: # Brute-force method to cope with the varying lengths of SSCs.
            self.PRIMARY_SSC = '0'+temp[0]
        elif len(temp[0])==3:
            self.PRIMARY_SSC = '00'+temp[0]
        elif len(temp[0])==2:
            self.PRIMARY_SSC = '000'+temp[0]
        elif len(temp[0])==1:
            self.PRIMARY_SSC = '0000'+temp[0]
        else:
            self.PRIMARY_SSC = temp[0]
            
        if len(temp[1])==4:
            self.SECONDARY_SSC = '0'+temp[1]
        elif len(temp[1])==3:
            self.SECONDARY_SSC = '00'+temp[1]
        elif len(temp[1])==2:
            self.SECONDARY_SSC = '000'+temp[1]
        elif len(temp[1])==1:
            self.SECONDARY_SSC = '0000'+temp[1]
        else:
            self.SECONDARY_SSC = temp[1]
            
        try:
            # N.B. need temp[2].split(".") as in some format of the output TCA microseconds will be included.
            self.TCA = datetime.datetime.strptime( temp[2].split(".")[0], "%d/%m/%Y %H:%M:%S" ) + datetime.timedelta(microseconds=int(temp[2].split(".")[1])) # Hack to handle decimal fractions of a second with a leading 0.
        except IndexError: # Sometimes there will be no microseconds information, that isn't a problem.
            self.TCA = datetime.datetime.strptime( temp[2].split(".")[0], "%d/%m/%Y %H:%M:%S" )
        except ValueError: # Sometimes the microseconds will be a very small float in exponential notation, approximate to 0.
            if temp[2].split(".")[0][-2:]=='60': # datetime wants seconds to be [0,60), workaround.
                self.TCA = datetime.datetime.strptime( temp[2].split(".")[0][:-2], "%d/%m/%Y %H:%M:" ) # Simply ignore seconds here.
                self.TCA = self.TCA+datetime.timedelta(minutes=1) # Manuallly increase this.
            else:
                self.TCA = datetime.datetime.strptime( temp[2].split(".")[0], "%d/%m/%Y %H:%M:%S" )
                
        self.MISS_DISTANCE = float( temp[3] )
        if 'nan' in temp[4] or '-1.#IND' in temp[4] or '1.#INF' in temp[4]:
            raise Warning("NaN in max probability for {} and {}.".format(self.PRIMARY_SSC, self.SECONDARY_SSC), self.PRIMARY_SSC, self.SECONDARY_SSC)
        elif 'nan' in temp[5] or '-1.#IND' in temp[5] or '1.#INF' in temp[5]: # Set a very low true collision probability. Maximum goes unchanged as it's OK.
            self.MAX_COLLISION_PROBABILITY = float( temp[4] )
            self.NO_COLLISION_PROBABILITY_MAX = 1.0 - self.MAX_COLLISION_PROBABILITY
            self.TRUE_COLLISION_PROBABILITY = sys.float_info.min
            self.NO_COLLISION_PROBABILITY_TRUE = 1.0 - self.TRUE_COLLISION_PROBABILITY
        else:
            if float( temp[4] )>1.0:
                raise Warning("Maximum probability for conjunction between {} and {} is greater than 1.0".format(self.PRIMARY_SSC, self.SECONDARY_SSC), self.PRIMARY_SSC, self.SECONDARY_SSC)
            elif float( temp[5] )>1.0:
                raise Warning("True probability for conjunction between {} and {} is greater than 1.0".format(self.PRIMARY_SSC, self.SECONDARY_SSC), self.PRIMARY_SSC, self.SECONDARY_SSC)
            else:
                self.MAX_COLLISION_PROBABILITY = float( temp[4] )
                self.NO_COLLISION_PROBABILITY_MAX = 1.0 - self.MAX_COLLISION_PROBABILITY
                self.TRUE_COLLISION_PROBABILITY = float( temp[5] )
                self.NO_COLLISION_PROBABILITY_TRUE = 1.0 - self.TRUE_COLLISION_PROBABILITY

        self.RELATIVE_V = float( temp[6] )
        self.COLLISION_RADIUS = float( temp[7] )
        
    def __str__(self):
        return "{},{},{},{},{},{}\n".format(self.PRIMARY_SSC, self.SECONDARY_SSC, self.TCA, self.MISS_DISTANCE, self.MAX_COLLISION_PROBABILITY, self.TRUE_COLLISION_PROBABILITY)

class SimpleConjunction( object ):
    def __init__(self, primarySSC, secondarySSC, TCA, MD, MaxPc, TruePc):
        """ A conjunction object created from already parsed information.
        Parameters
        ----------
            primarySSC : string
                NORAD ID of the primary object.
            secondarySSC : string
                NORAD ID of the secondary object.
            TCA : datetime.Datetime
                Time of Closest Approach of the conjunction.
            MD : float
                Miss distance in km
            MaxPc : float
                Maximum collision probability.
            TruePc : float
                True collision probability.
        """
        
        if primarySSC.endswith("    "): # Deal with STK CAT output format. Go from the most to fewest whitespaces to avoid problems (something with two whitespaces at the end also endswith(" ") ).
            self.PRIMARY_SSC  = "0000"+primarySSC.strip("    ")  
        elif primarySSC.endswith("  "):
            self.PRIMARY_SSC  = "000"+primarySSC.strip("   ")
        elif primarySSC.endswith("  "):
            self.PRIMARY_SSC  = "00"+primarySSC.strip("  ")
        elif primarySSC.endswith(" "):
            self.PRIMARY_SSC = "0"+primarySSC.strip(" ")
        else:
            self.PRIMARY_SSC = primarySSC
            

        if secondarySSC.endswith("    "):
            self.SECONDARY_SSC = "0000"+secondarySSC.strip("    ")
        elif secondarySSC.endswith("   "):
            self.SECONDARY_SSC = "000"+secondarySSC.strip("   ")
        elif secondarySSC.endswith("  "):
            self.SECONDARY_SSC = "00"+secondarySSC.strip("  ")
        elif secondarySSC.endswith(" "):
            self.SECONDARY_SSC = "0"+secondarySSC.strip(" ")
        else:
            self.SECONDARY_SSC = secondarySSC
        
        self.TCA = TCA
        self.MISS_DISTANCE = MD
        self.MAX_COLLISION_PROBABILITY = MaxPc
        self.TRUE_COLLISION_PROBABILITY = TruePc
        self.NO_COLLISION_PROBABILITY_MAX = 1.0 - self.MAX_COLLISION_PROBABILITY
        self.NO_COLLISION_PROBABILITY_TRUE = 1.0 - self.TRUE_COLLISION_PROBABILITY
    
    def __str__(self):
        return "{},{},{},{},{},{}\n".format(self.PRIMARY_SSC, self.SECONDARY_SSC, self.TCA, self.MISS_DISTANCE, self.MAX_COLLISION_PROBABILITY, self.TRUE_COLLISION_PROBABILITY)
        
"""
============================================================================================================================
    IMPORT DATA FROM C++ FILE.
    Each file should contain conjunctions information and a header. The conjunction format is expected as a CSV file
    with the following collumns:
    SSC1, SSC2, TCA, Miss Distance (km), Max Pc, True Pc, relative V at TCA (km/s), collision radius (km)
============================================================================================================================
"""
decayEpochs = {} # Get the reentry epoch of all the objects and store accoring to the SSC.
decayedInSimulation = 0 # Decayed in the duration of the simulation.
decayedOutsideSimulation  = 0
with open("decayEpochs_23Oct2013","r") as deFile:
    deLines = deFile.readlines()
    for line in deLines:
        try:
            splitLine = line.split()
            decayEpochs[ splitLine[0] ] = datetime.datetime.strptime(splitLine[1]+" "+splitLine[2], '%Y-%m-%d %H:%M:%S')
            if datetime.datetime.strptime(splitLine[1]+" "+splitLine[2], '%Y-%m-%d %H:%M:%S') <= datetime.datetime(2013, 11, 23, 3, 58, 21):
                decayedInSimulation += 1
            else:
                decayedOutsideSimulation += 1
        except IndexError: # This object hasn't decayed yet.
            decayEpochs[ splitLine[0] ] = datetime.datetime(2200,1,1,0,0,0) # Very distant future just to have something.

fileName = 'envisat_1year_COV3_1kmStandardDeviations_DI'
lines=[]; # Read lines of the datafiles, one conjunction per line.
totalConjunctionsLoaded=0; # Nmber of conjunctions loaded from the files.

utilities.customPrint("Reading data from {}.".format(fileName), utilities.GREEN)
conjunctionsRead=0 # Initialise this
" Read the data. "
with open(fileName,"r") as inFile:
    tempLines = inFile.readlines(); # Lines from the given file only.

" Split the output header and conjunctions' information of generic length. "
header = [];
for i in range(0,len(tempLines)): # this was from 0 to len(...)+1
    try:
        int( tempLines[i][3] ) # This won't work for the header, only when applyint int to a conjunction (SSC is recorded first).
        header = tempLines[:i] # Record the header info here...
        tempLines = tempLines[i:] #...and the conjunctions info here.
        break # Found end of the header, stop now.
    except ValueError:
        if not tempLines[i].startswith("Conjnctions found:"):
            pass # Will get this when trying to convert a string to int i.e. when still parsing the header.
        else: # See if this line contains info about the number of conjunctions to be found in this file.
            conjunctionsRead = int( tempLines[i][19:] )
            totalConjunctionsLoaded = totalConjunctionsLoaded + conjunctionsRead
            
        
" Only save the conjunctions, discard the rest. "
utilities.customPrint("Read {} conjunctions from {}.".format(conjunctionsRead,fileName), utilities.BLUE)
lines.extend( tempLines );
del header, tempLines;

"""
============================================================================================================================
    IMPORT DATA FROM STK CAT FILE.
    Each file should contain conjunctions information and a header. The conjunction format is expected as a CSV file
    with the following collumns:
    SSC1, SSC2, TCA, Miss Distance (km), Max Pc, True Pc, relative V at TCA (km/s), collision radius (km)
============================================================================================================================
"""

STKconjunctions = [] # Different ellipsoid size (standard deviation squared c.f. standard deviation).
STK_SSCs=[]; STK_TCAs=[]; STKmissDistances=[]; STKtrueProbabilities=[]; STKmaxProbabilities=[];
with open('sgp4ESSS verificationReport COV1.csv', 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    rows = []
    for tempRow in reader:
        rows.append(tempRow)
    for row in rows[30:-1]: # Don't go through the header and the last line.
        tempConj = SimpleConjunction("27386", row[0].replace('"',''), datetime.datetime.strptime(row[1], "%d %b %Y %H:%M:%S.%f"), float(row[2]), float(row[3]), float(row[5]))
        if len(tempConj.SECONDARY_SSC)==4: # Correct for the leading 0s that STK does not output.
            tempConj.SECONDARY_SSC = '0'+tempConj.SECONDARY_SSC
        elif len(tempConj.SECONDARY_SSC)==3:
            tempConj.SECONDARY_SSC = '00'+tempConj.SECONDARY_SSC
        elif len(tempConj.SECONDARY_SSC)==2:
            tempConj.SECONDARY_SSC = '000'+tempConj.SECONDARY_SSC
        elif len(tempConj.SECONDARY_SSC)==1:
            tempConj.SECONDARY_SSC = '0000'+tempConj.SECONDARY_SSC
            
        if not decayEpochs[ tempConj.PRIMARY_SSC ] <= tempConj.TCA  and not decayEpochs[ tempConj.SECONDARY_SSC ] <= tempConj.TCA:
            if True: #tempConj.MISS_DISTANCE >= 5.0:
                STKconjunctions.append( tempConj )
                STK_SSCs.append( tempConj.SECONDARY_SSC )
                STK_TCAs.append( tempConj.TCA )
                STKmissDistances.append( tempConj.MISS_DISTANCE )
                STKmaxProbabilities.append( tempConj.MAX_COLLISION_PROBABILITY )
                STKtrueProbabilities.append( tempConj.TRUE_COLLISION_PROBABILITY )
        
STKconjunctions.sort(key=lambda x: x.TCA)
ACCUMULATED_PCs_MAX_STK = []; ACCUMULATED_PCs_TRUE_STK = []; accumulatedPcMAX_STK2 = 1.0; accumulatedPcTRUE_STK2 = 1.0;
for conjSTK in STKconjunctions: # Save the conjunctions' data once they've been sorted by TCA ascending.
    accumulatedPcMAX_STK2 = accumulatedPcMAX_STK2*conjSTK.NO_COLLISION_PROBABILITY_MAX
    accumulatedPcTRUE_STK2 = accumulatedPcTRUE_STK2*conjSTK.NO_COLLISION_PROBABILITY_TRUE
    
    ACCUMULATED_PCs_MAX_STK.append( 1.0-accumulatedPcMAX_STK2 )
    ACCUMULATED_PCs_TRUE_STK.append( 1.0-accumulatedPcTRUE_STK2 )
       
"""
============================================================================================================================
    PARSE THE CONJUNCTIONS INFORMATION TO GET GLOBAL STATISTICS ABOUT ALL OF THEM.
============================================================================================================================
"""
conjunctions = []; missDistances=[]; maxProbabilities = []; trueProbabilities = []; TCAs = [];
NO_BAD_PROBABILITY_CONJUNCTIONS = 0
NO_NANS_IN_CONJUNCTIONS = 0
NO_maxPC_LOWER_THAN_truePC = 0; badPCsEpochs=[]; badMaxPCs=[]; badTruePCs=[]; badPCsPrimaries=[]; badPCsSecondaries=[] # Epochs and probabilities for this case.
NO_TOO_LOW_maxPC = 0; tooLowMaxPCsPrimaries=[]; tooLowMaxPCsSecondaries=[]; # SSCs of the objets that have conjunctions with Pc MAX less than 1.0 even though it should be 1.0.

for conjLine in lines: # Create conjunctions from the records for ease of processing
    try:
        tempConj = Conjunction(conjLine)
        if not decayEpochs[ tempConj.PRIMARY_SSC ] <= tempConj.TCA  and not decayEpochs[ tempConj.SECONDARY_SSC ] <= tempConj.TCA: # Filter out the TLEs that have decayed already so this conjunction can't have taken place.
            conjunctions.append( tempConj )
            # Do consistency checks here
            if tempConj.COLLISION_RADIUS >= tempConj.MISS_DISTANCE and tempConj.MAX_COLLISION_PROBABILITY<1.0:
                utilities.customPrint("Miss distance lower than collision radius and Pc max less than 1.0", utilities.RED)
                utilities.customPrint( str(tempConj), utilities.GRAY)
                NO_TOO_LOW_maxPC += 1
                tooLowMaxPCsPrimaries.append(tempConj.PRIMARY_SSC)
                tooLowMaxPCsSecondaries.append(tempConj.SECONDARY_SSC)
                
            if tempConj.MAX_COLLISION_PROBABILITY < tempConj.TRUE_COLLISION_PROBABILITY:
                utilities.customPrint("Maximum collision probability lower than true collision probability", utilities.RED)
                utilities.customPrint( str(tempConj), utilities.GRAY)
                NO_maxPC_LOWER_THAN_truePC += 1
                badPCsEpochs.append(tempConj.TCA) # Store info about these erroneous cases for analysis.
                badMaxPCs.append(tempConj.MAX_COLLISION_PROBABILITY)
                badTruePCs.append(tempConj.TRUE_COLLISION_PROBABILITY)
                badPCsPrimaries.append(tempConj.PRIMARY_SSC)
                badPCsSecondaries.append(tempConj.SECONDARY_SSC)
            
    except Warning as wrng: # There seems to be something wrong with this conjunction. But maybe it involved an ignored TLE.
        if "NaN" in wrng.args[0] and not decayEpochs[ wrng.args[1] ] <= tempConj.TCA  and not decayEpochs[ wrng.args[2] ] <= tempConj.TCA: # wrng will have primary and secondary objects' SSCs - don't care about bad values for objects from the ignored list.
            utilities.customPrint( str(wrng.args), utilities.RED) # Print the error here as don't want this to be displayed if it came from one of the ignored TLEs.
            utilities.customPrint("\t"+conjLine, utilities.GRAY)
            NO_NANS_IN_CONJUNCTIONS += 1
        elif not decayEpochs[ wrng.args[1] ] <= tempConj.TCA  and not decayEpochs[ wrng.args[2] ] <= tempConj.TCA:
            utilities.customPrint(wrng.args[0], utilities.RED) # Print the error here as don't want this to be displayed if it came from one of the ignored TLEs.
            utilities.customPrint("\t"+conjLine, utilities.GRAY)
            NO_BAD_PROBABILITY_CONJUNCTIONS += 1

conjunctions.sort(key=lambda x: x.TCA) # Sort by TCA ascending.
utilities.customPrint("Found {} conjunctions out of {} that have collision probability values of either collision probability > 1.0".format(NO_BAD_PROBABILITY_CONJUNCTIONS, len(conjunctions)), utilities.RED)
utilities.customPrint("Found {} conjunctions out of {} that have NaN values of maximum collision probability".format(NO_NANS_IN_CONJUNCTIONS, len(conjunctions)), utilities.RED)
utilities.customPrint("Found {} conjunctions out of {} that have maximum Pc <1.0 when collision radius exceeds the miss distance".format(NO_TOO_LOW_maxPC, len(conjunctions)), utilities.RED)
utilities.customPrint("Found {} conjunctions out of {} that have maximum Pc lower than true Pc".format(NO_maxPC_LOWER_THAN_truePC, len(conjunctions)), utilities.RED)

ACCUMULATED_PCs_MAX = []; ACCUMULATED_PCs_TRUE = []; accumulatedPcMAX = 1.0; accumulatedPcTRUE = 1.0;
for conj in conjunctions: # Save the conjunctions' data once they've been sorted by TCA ascending.
    missDistances.append( conj.MISS_DISTANCE )
    maxProbabilities.append( conj.MAX_COLLISION_PROBABILITY )
    trueProbabilities.append( conj.TRUE_COLLISION_PROBABILITY )
    TCAs.append( conj.TCA )

    accumulatedPcMAX = accumulatedPcMAX*conj.NO_COLLISION_PROBABILITY_MAX
    accumulatedPcTRUE = accumulatedPcTRUE*conj.NO_COLLISION_PROBABILITY_TRUE
    
    ACCUMULATED_PCs_MAX.append( 1.0-accumulatedPcMAX )
    ACCUMULATED_PCs_TRUE.append( 1.0-accumulatedPcTRUE )

#" Make a histogram of misss distances. "
#binsMD = [0.5, 1.0, 5.0, 10.0, 15.0, 20.0, 25.0]
#pylab.figure()
#nMD, binsMD, patchesMD = pylab.hist( [missDistances, STKmissDistances], binsMD, rwidth=0.8, histtype='bar',color=['r', 'b'],label=['C++', 'STK CAT'])
#pylab.xticks(binsMD, rotation='vertical')
#pylab.xlabel('Miss distance bin (km)')
#pylab.ylabel('Conjunction count')
#pylab.legend()
#pylab.show()

#" Make a histogram of maximum probabilities. "
#binsPC = numpy.logspace(-5, 2, num=10)
#PChistFigure=pylab.figure()
#PChistAxes=PChistFigure.gca()
#nPC, binsPC, patchesPC = pylab.hist( [maxProbabilities, STKmaxProbabilities], binsPC, rwidth=0.8, histtype='bar')
#PChistAxes.set_xscale("log")
#PChistAxes.set_xlabel('Maximum collision probability')
#PChistAxes.set_ylabel('Conjunction count')
#pylab.show()

#" Make a histogram of true probabilities. "
#binsPCtrue = numpy.logspace(-30, 1, num=10)
#PChistFigureTRUE=pylab.figure()
#PChistAxesTRUE=PChistFigureTRUE.gca()
#PChistAxesTRUE.hist( [trueProbabilities, STKtrueProbabilities], binsPCtrue, rwidth=0.8, histtype='bar')
#PChistAxesTRUE.set_xscale("log")
#PChistAxesTRUE.set_xlabel('True collision probability')
#PChistAxesTRUE.set_ylabel('Conjunction count')
#pylab.show()


"""
============================================================================================================================
    PLOT EVOLUTION OF ACCUMULATED PROBABILITIES.
============================================================================================================================
"""
    
# FORMAT THE PLOT
ticksFontSize = 18
labelsFontSize = 30
titleFontSize = 34

matplotlib.rc('xtick', labelsize=ticksFontSize) 
matplotlib.rc('ytick', labelsize=ticksFontSize)  

fig, ax = matplotlib.pyplot.subplots(1,figsize=(12,8))
fig.suptitle(r'$Envisat,\ one\ year,\ conjunctions\ closer\ than\ 20\ km$', fontsize=titleFontSize) #\ external\ ephemeris

matplotlib.pyplot.grid(linewidth=2)
ax.tick_params(axis='both',reset=False,which='both',length=5,width=1.5)

ax.set_xlabel(r'$Epoch\ (UTCG)$',fontsize=labelsFontSize)
ax.set_ylabel(r'$Accumulated\ collision\ probability$',fontsize=labelsFontSize)

matplotlib.pyplot.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)

# DATA PLOTTING
ax.plot(TCAs, ACCUMULATED_PCs_MAX, c="0.5", ls='-', lw=5., label=r'$Maximum\ probability$')
ax.plot(TCAs, ACCUMULATED_PCs_TRUE, c="0.5", ls='--', lw=5., label=r'$True\ probability$')

ax.plot(STK_TCAs, ACCUMULATED_PCs_MAX_STK, c='k', ls='-', lw=5., label=r'$Maximum\ probability,\ STK\ CAT$')
ax.plot(STK_TCAs, ACCUMULATED_PCs_TRUE_STK, c='k', ls='--', lw=5., label=r'$True\ probability,\ STK\ CAT$')
    
# Shink current axis' height by 15% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.85])

matplotlib.pyplot.setp( ax.xaxis.get_majorticklabels(), rotation=30 )

# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.17), fancybox=True, shadow=True, ncol=4, prop={'size':18})
    
fig.show()

"""
============================================================================================================================
    PLOT HISTORY OF INDIVIDUAL MAXIMUM PROBABILITIES.
============================================================================================================================
"""

matplotlib.rc('xtick', labelsize=ticksFontSize) 
matplotlib.rc('ytick', labelsize=ticksFontSize)  

fig2, ax2 = matplotlib.pyplot.subplots(1,figsize=(12,8))
ax2.set_yscale("log")
fig2.suptitle(r'$Envisat,\ one\ year$', fontsize=titleFontSize)

matplotlib.pyplot.grid(linewidth=2)
ax2.tick_params(axis='both',reset=False,which='both',length=5,width=1.5)

ax2.set_xlabel(r'$Epoch\ (UTCG)$',fontsize=labelsFontSize)
ax2.set_ylabel(r'$Maximum\ collision\ probability$',fontsize=labelsFontSize)

matplotlib.pyplot.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)

# DATA PLOTTING
ax2.scatter(TCAs, maxProbabilities, c='r', marker='+', label=r'$Maximum\ probability,\ C++$')

ax2.scatter(STK_TCAs, STKmaxProbabilities, c='b', marker='x', label=r'$Maximum\ probability,\ STK\ CAT$')
    
# Shink current axis' height by 15% on the bottom
box = ax2.get_position()
ax2.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.85])

matplotlib.pyplot.setp( ax2.xaxis.get_majorticklabels(), rotation=30 )

# Put a legend below current axis
ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.17), fancybox=True, shadow=True, ncol=4, prop={'size':18})
    
fig2.show()

"""
============================================================================================================================
    PLOT HISTORY OF INDIVIDUAL TRUE PROBABILITIES.
============================================================================================================================
"""

matplotlib.rc('xtick', labelsize=ticksFontSize) 
matplotlib.rc('ytick', labelsize=ticksFontSize)  

fig3, ax3 = matplotlib.pyplot.subplots(1,figsize=(12,8))
ax3.set_yscale("log")
fig3.suptitle(r'$Envisat,\ one\ year,\ external\ ephemeris$', fontsize=titleFontSize)

matplotlib.pyplot.grid(linewidth=2)
ax3.tick_params(axis='both',reset=False,which='both',length=5,width=1.5)

ax3.set_xlabel(r'$Epoch\ (UTCG)$',fontsize=labelsFontSize)
ax3.set_ylabel(r'$True collision\ probability$',fontsize=labelsFontSize)

matplotlib.pyplot.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)

# DATA PLOTTING
ax3.scatter(TCAs, trueProbabilities, c='r', marker='+', label=r'$True probability,\ C++$')

ax3.scatter(STK_TCAs, STKtrueProbabilities, c='b', marker='x', label=r'$True\ probability,\ STK\ CAT$')
    
# Shink current axis' height by 15% on the bottom
box = ax3.get_position()
ax3.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.85])

matplotlib.pyplot.setp( ax3.xaxis.get_majorticklabels(), rotation=30 )

# Put a legend below current axis
ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.17), fancybox=True, shadow=True, ncol=4, prop={'size':18})
    
fig3.show()

"""
============================================================================================================================
    PLOT HISTORY OF INDIVIDUAL MISS DISTANCES.
============================================================================================================================
"""

matplotlib.rc('xtick', labelsize=ticksFontSize) 
matplotlib.rc('ytick', labelsize=ticksFontSize)  

fig4, ax4 = matplotlib.pyplot.subplots(1,figsize=(12,8))
fig4.suptitle(r'$Envisat,\ one\ year$', fontsize=titleFontSize)

matplotlib.pyplot.grid(linewidth=2)
ax4.tick_params(axis='both',reset=False,which='both',length=5,width=1.5)

ax4.set_xlabel(r'$Epoch\ (UTCG)$',fontsize=labelsFontSize)
ax4.set_ylabel(r'$Miss\ distance\ (km)$',fontsize=labelsFontSize)

matplotlib.pyplot.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)

# DATA PLOTTING
ax4.scatter(TCAs, missDistances, c='r', marker='+', label=r'$C++$')

ax4.scatter(STK_TCAs, STKmissDistances, c='b', marker='x', label=r'$STK\ CAT$')
    
# Shink current axis' height by 15% on the bottom
box = ax4.get_position()
ax4.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.85])

matplotlib.pyplot.setp( ax4.xaxis.get_majorticklabels(), rotation=30 )

# Put a legend below current axis
ax4.legend(loc='upper center', bbox_to_anchor=(0.5, -0.17), fancybox=True, shadow=True, ncol=4, prop={'size':18})
    
fig4.show()

"""
============================================================================================================================
    PLOT EVOLUTION OF ACCUMULATED TRUE PROBABILITIES ONLY.
    True probabilities are computed using different algorithms, maximum one is
    computed using S. Alfano so will be the same no matter what size of the
    covaraince ellipsoid.
============================================================================================================================
"""
matplotlib.rc('xtick', labelsize=ticksFontSize) 
matplotlib.rc('ytick', labelsize=ticksFontSize)  

fig5, ax5 = matplotlib.pyplot.subplots(1,figsize=(12,8))
fig5.suptitle(r'$Envisat,\ one\ year$', fontsize=titleFontSize)
#ax5.set_yscale("log")

matplotlib.pyplot.grid(linewidth=2)
ax5.tick_params(axis='both',reset=False,which='both',length=5,width=1.5)

ax5.set_xlabel(r'$Epoch\ (UTCG)$',fontsize=labelsFontSize)
ax5.set_ylabel(r'$Accumulated\ true\ collision\ probability$',fontsize=labelsFontSize)

matplotlib.pyplot.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)

# DATA PLOTTING
ax5.plot(TCAs, ACCUMULATED_PCs_TRUE, c='r', ls='--', lw=3., label=r'$C++$')
#ax5.plot(STK_TCAs, ACCUMULATED_PCs_TRUE_STK, c='g', ls='--', lw=3., label=r'$STK\ CAT ellipsoid\ size\ is\ \sigma^2$')
ax5.plot(STK_TCAs, ACCUMULATED_PCs_TRUE_STK, c='b', ls='--', lw=3., label=r'$STK\ CAT$')
    
# Shink current axis' height by 15% on the bottom
box = ax5.get_position()
ax5.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.85])

matplotlib.pyplot.setp( ax.xaxis.get_majorticklabels(), rotation=30 )

# Put a legend below current axis
ax5.legend(loc='upper center', bbox_to_anchor=(0.5, -0.17), fancybox=True, shadow=True, ncol=4, prop={'size':18})
    
fig5.show()
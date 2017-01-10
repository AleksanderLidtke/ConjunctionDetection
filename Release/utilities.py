# -*- coding: utf-8 -*-
"""
Created on Fri Nov 08 09:46:06 2013

Various plotting and output utilities useful for creating 3d Cartesian plots, priting to stdout in colour etc.
Also contains a generic class that can be used to handle SGP4 object creation and propagation.

@author: Aleksander Lidtke
@version 1.6.4
@since 12/08/2014 10:10:00
@email al11g09@soton.ac.uk alek_l@onet.eu

CHANGELOG:
30/01/2014 - 1.0.0 - removed the previous position storing from the SpaceObject - this will be retained externally.
03/03/2014 - 1.0.1 - added the calculatePerigeeRadiusWGS84 and calculateApogeeRadiusWGS84 methods to the SpaceObject to save the radii data during propagation (necessary for pre-filters' implementation).
05/03/2014 - 1.1.0 - added CURRENT_POSES and CURRENT_VELOS to the SpaceObject class attributes.
06/03/2014 - 1.2.0 - added PowerSeriesPolynomialInterpCoefficients3d function that is faster to execute than the original interpolating function PowerSeriesPolynomialInterp3d.
07/03/2014 - 1.3.0 - added a fix to the Hermite interpolation - need to nondimensionalise the gradients when provided in the form of velocity (that is a gradient w.r.t. time, not the lengthwise coordinate).
13/03/2014 - 1.4.0 - added the perigee/apogee filtering routines.
           - 1.4.0 - changed the SpaceObject.calculateRadii methods to use the CURRENT_POS and CURRENT_V for calculation and to return the radii.
14/03/2014 - 1.5.6 - SpaceObject now uses 2 points for interpolation of its trajectory as default - compatible with PowerSeriesPolynomialInterpGradientsCoefficients3d that uses function value and derivative.
           - 1.5.6 - Changed the divisions being done in conversion between JDAYs and seconds to multiplication for speed.
25/03/2014 - 1.6.0 - Added the conjunction class to stroe information about conjunctions.
28/03/2014 - 1.6.1 - Added the capability to change the default radius of objects.
31/03/2014 - 1.6.2 - Added a multiplier that can be applied to B* coefficients of objects in order to simulate variations in solar activity.
05/04/2014 - 1.6.3 - Added the capability to use different default radii for R/B, P/L, DEB and other types of objects.
12/08/2014 - 1.6.4 - Removed sgp4. from some imports to work on all machines without SGP4 added to PythonPath.
"""
import sys, numpy, math, datetime

from earth_gravity import wgs84;
from propagation import getgravconst
from sgp4io import twoline2rv
from ext import jday, invjday
from ext import rv2coe

GRAVITY_MODEL = 'wgs84';
MEAN_EARTH_RADIUS = getgravconst(GRAVITY_MODEL)[2] # km
GRAVITY_CONSTANT = getgravconst(GRAVITY_MODEL)[1] # km^3/s^2
SURFACE_GRAVITY_ACCELERATION = GRAVITY_CONSTANT/MEAN_EARTH_RADIUS**2 # km/s^2
ESCAPE_VELOCITY = math.sqrt( 2*GRAVITY_CONSTANT/MEAN_EARTH_RADIUS ) # km/s

threeLE_FILE_NAME='spaceTrack3LEs_23_10_2013.txt'

TLE1="""0 VANGUARD 1\n1 00005U 58002B   13293.59932806  .00000150  00000-0  16823-3 0  1542\n2 00005 034.2538 157.6066 1849678 212.0446 135.3229 10.84362834941283""" # Reference TLE for testing.
TLE_Envisat="""0 Envisat\n1 27386U 02009A   13296.16551405  .00000100  00000-0  47533-4 0  8339\n2 27386 098.4160 001.9229 0001126 079.5538 280.5787 14.37637337609529"""

GRAY, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8) # Colours for pretty printing.
SECOND_IN_JDAYS = 1.157401129603386e-05 # 1 second expressed in JDAYS. JDAY/s
JDAY_IN_SECONDS = 1./SECOND_IN_JDAYS # 1 JDAY expressed in seconds. s/JDAY

def customPrint(text, colourCode=GRAY):
    """ Prints some text to stdout in the given colour. Adds a new line at the end of text.\n
    Parameters
    ----------
        @param text - the string to be printed.
        @param colourCode - code of the colour to be used (e.g. 0,1,2,3,4) """
    output = "\x1b[1;%dm" % (30+colourCode) + text + "\x1b[0m\n"
    sys.stdout.write(output)
    
def centralLinspace(centre, spacing, noPoints):
    """ Given a number at the centre of an interval create a set of evenly spaced
    numbers centred on the given cetre.\n
    Parameters
    ----------
        @param centre - float giving the centre of the interval.\n
        @param spacing - the separation between the numbers in the interval.\n
        @params noPoints - number of points in the desired interval. Must be even.\n
    Retruns
    ----------
        @return - numpy.array 1XnoPoints containing evenly spaced points, separated by spacing and centred on centre. Contains the endpoints.\n
    Raises
    ----------
        @raises ValueError if noPoints is not an even number.
    """
    if noPoints%2 != 0:
        raise ValueError("Requested {} points to be seeded around a central location but the number of points must be even!".format(noPoints))
    
    points = numpy.linspace( centre-1.5*spacing, centre+1.5*spacing, noPoints)
    
    return points

def HermiteCubicInterp3d(X, Y, Z, noInterpolationPoints=500, dXdT=None, dYdT=None, dZdT=None, nodeEpochs=None):
    """ Interpolate a set of points in space using piecewise cubic Hermite polynomials which coserve the
    first derivative of the position (X, Y, Z) along the curve traced by the provided points. If not provided,
    his gradient will be consistent from interval to interval. If unsupplied it will be computed using the finite
    diffrerence method: central in interior, backwards or forwards at the endpoints).\n
    Parameters
    ----------
        @param X, Y, Z - coordinates of the know points that will be interpolated.\n
        @param times - values of the independtent variable that governs X, Y and Z corresponding to the supplied interpolation points.\n
        @param noInterpolationPoints - number of points to be returned that are evenly placed along the curve formed by the supplied points.\n
        @param dXdT, dYdT, dZdT - gradients of the coordinates (X, Y, Z) computed at each of the
            provided points w.r.t. time.\n
        @param nodeEpochs - epochs corresponding to the nodes at which the function values and gradients are supplied. MUST be provided if any 
            of the gradients is provided as it's necessary for nondimensionalisation.\n
    Returns
    ----------
        @return - noInterpolationPoints X 3 numpy.array containing X, Y and Z coordinates (in columns) of the interpolated points."""
    S = numpy.linspace(0., 1., len(X)) # Lengthwise coordinate along the curve.

    " Compute the gradients if not supplied. "
    if dXdT is None:
        dXdT = numpy.zeros(len(X))
        for i in range(len(X)):
            """ dXdS """
            if i == 0: # First point <=> FWD difference.
                dXdT[i] = (X[i+1]-X[i])/(S[i+1]-S[i])  
            elif i == len(X)-1: # Last point <=> BWD difference.
                dXdT[i] = (X[i]-X[i-1])/(S[i]-S[i-1])
            else: # Central difference in interior points.
                dXdT[i] = (X[i+1]-X[i-1])/(S[i+1]-S[i-1])
    else: # Nondimensionalise the gradient.
        dXdT = (nodeEpochs[-1]-nodeEpochs[0])*dXdT*JDAY_IN_SECONDS
            
    if dYdT is None:
        dYdT = numpy.zeros(len(Y))
        for i in range(len(Y)):
            if i == 0: # First point <=> FWD difference.
                dYdT[i] = (Y[i+1]-Y[i])/(S[i+1]-S[i])  
            elif i == len(Y)-1: # Last point <=> BWD difference.
                dYdT[i] = (Y[i]-Y[i-1])/(S[i]-S[i-1])
            else: # Central difference in interior points.
                dYdT[i] = (Y[i+1]-Y[i-1])/(S[i+1]-S[i-1])
    else: # Nondimensionalise the gradient.
        dYdT = (nodeEpochs[-1]-nodeEpochs[0])*dYdT*JDAY_IN_SECONDS
            
    if dZdT is None:
        dZdT = numpy.zeros(len(Z))
        for i in range(len(Z)):
            if i == 0: # First point <=> FWD difference.
                dZdT[i] = (Z[i+1]-Z[i])/(S[i+1]-S[i])  
            elif i == len(Z)-1: # Last point <=> BWD difference.
                dZdT[i] = (Z[i]-Z[i-1])/(S[i]-S[i-1])
            else: # Central difference in interior points.
                dZdT[i] = (Z[i+1]-Z[i-1])/(S[i+1]-S[i-1])
    else: # Nondimensionalise the gradient.
        dZdT = (nodeEpochs[-1]-nodeEpochs[0])*dZdT*JDAY_IN_SECONDS
    
    " Interpolate "
    sDesired = numpy.linspace(0., 1., noInterpolationPoints) # Points at which we want to interpolate the function.
    
    desiredPoints = numpy.zeros( (noInterpolationPoints,3) ) # X, Y and Z interpolated along S.
    
    for j in range(noInterpolationPoints):
        Sk = S[0]; Sk_1 = S[1]; # Known arguments binding the currently sought argument value.
        indLower=0; indUpper=1; # Indices of the binding known points.
        for k in range(len(X)): # Find the binding interval for the current point to be interpolated.
            if S[k] < sDesired[j] <= S[k+1]:
                Sk=S[k]
                indLower=k
                Sk_1=S[k+1]
                indUpper=k+1
                
        t = (sDesired[j]-Sk)/(Sk_1-Sk) # Dimensionless parameter - location of the point to be interpolated in its interval.
    
        # Hermite basis functions in the expanded form.
        h00 =  2*t*t*t  -  3*t*t  +  1
        h10 =    t*t*t  -  2*t*t  +  t
        h01 = -2*t*t*t  +  3*t*t
        h11 =    t*t*t  -    t*t
        
        # Compute the interpolated values of the function at each desired location from the basis functions and other parameters.
        desiredPoints[j,0] = h00*X[indLower] + h10*(Sk_1-Sk)*dXdT[indLower] + h01*X[indUpper] + h11*(Sk_1-Sk)*dXdT[indUpper]
        desiredPoints[j,1] = h00*Y[indLower] + h10*(Sk_1-Sk)*dYdT[indLower] + h01*Y[indUpper] + h11*(Sk_1-Sk)*dYdT[indUpper]
        desiredPoints[j,2] = h00*Z[indLower] + h10*(Sk_1-Sk)*dZdT[indLower] + h01*Z[indUpper] + h11*(Sk_1-Sk)*dZdT[indUpper]

    return desiredPoints


def LagrangeInterp3d(X, Y, Z, noInterpolationPoints=500):
    """ Interpolate a set of points in space using Lagrange polynomials using all the points (X, Y, Z) that 
    are provided.\n
    Parameters
    ----------
        @param X, Y, Z - coordinates of the know points that will be interpolated.\n
        @param noInterpolationPoints - number of points to be returned that are evenly placed along the curve formed by the supplied points.\n
    Returns
    ----------
        @return - noInterpolationPoints X 3 numpy.array containing X, Y and Z coordinates (in columns) of the interpolated points."""    
    S = numpy.linspace(0., 1., len(X)) # Lengthwise coordinate along the curve.
    
    " Interpolate "
    sDesired = numpy.linspace(0., 1., noInterpolationPoints) # Points at which we want to interpolate the function.
    
    desiredPoints=numpy.zeros( (noInterpolationPoints,3) ) # X, Y and Z interpolated along S.
    
    def lagrangeCoeff(termIndex, desiredArgument):
        """ Calculate the Lagrange polynomial coefficient for the given point.\n
        Parameters
        ----------
            @param temIndex - index of the term for which the coefficient is to be ealuated.\n
            @param desiredX - value of the argument (x) for which the value is being sought.\n
        Returns
        ----------
            @return - float being the value of the Lagrange interpolating polynomial's coefficient for the given term (y)
            when looking for the value at a particular x.
        See Also
        ----------
            @see http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html for more details."""
        v = 1.0 # Each coefficient is a product of many terms. Initialise with value of 1.
        for k in xrange(len(S)): # Each coefficient has a term that includes argument values at other points.
            if k != termIndex:
                v *= (desiredArgument-S[k]) / (S[termIndex]-S[k])
        return v
    
    for i,xi in enumerate(sDesired): # Evaluate value of the function at each desired point.
        for j in xrange(len(S)): # Value at each interpolated point is partially affected by the values at all the known points.
            desiredPoints[i,0] += X[j]*lagrangeCoeff(j,xi)
            desiredPoints[i,1] += Y[j]*lagrangeCoeff(j,xi)
            desiredPoints[i,2] += Z[j]*lagrangeCoeff(j,xi)

    return desiredPoints
    
def PowerSeriesPolynomialInterp3d(X, Y, Z, interpolationArguments, returnCoefficients=False):
    """ Interpolate a set of points in space using a power series polynomial expansion (finding coefficients for the N-1 order
    polynomial to interpolate the N data points). Use all the points (X, Y, Z) that are provided.
    Can return the polynomial coefficients if desired.\n
    Parameters
    ----------
        @param X, Y, Z - coordinates of the know points that will be interpolated.\n
        @param interpolationArguments - points evenly placed along the curve formed by the supplied points at which the interpoalting function will be evaluated.
                                        Must be in the same order as X, Y and Z.\n
        @param returnCoefficients - if True polynomial coefficients will be returned.\n
    Returns
    ----------
        @return - by default noInterpolationPoints X 3 numpy.array containing X, Y and Z coordinates (in columns) of the interpolated points
                    can also be a 2-tuple with the polynomial coefficients in the form of a numpy array (X, Y and Z interpolating polynomial's 
                    coefficients in rows 0, 1 and 2, respecitvely) as the second entry.
                    The coefficients are given starting at the zeroth order term up to the one with the highest order."""
    desiredPoints=numpy.zeros( (len(interpolationArguments),3) )

    """ Build the matrix of coefficients to solve the system of equations:
    [1, argumnet, argument^2, ..., argument^(N-1)] * [unknown polynomial coefficients] = [known data point values]"""
    numberOfPointsProvided = len(X)
    A = numpy.zeros( (numberOfPointsProvided, numberOfPointsProvided) )
    for j in range(numberOfPointsProvided):
        A[:,j] = numpy.linspace(interpolationArguments[0], interpolationArguments[-1], numberOfPointsProvided)**j # Points that have been provided must be evenly spaced through the interval.
    
    XpolyCoeffs=numpy.linalg.solve(A, X) # Coefficients of the interpolating polynomials in X, Y and Z
    YpolyCoeffs=numpy.linalg.solve(A, Y)
    ZpolyCoeffs=numpy.linalg.solve(A, Z)
    
    for i in range(len(interpolationArguments)):
        for k in range(len(XpolyCoeffs)): # Don't a priori know how many coefficients there will be.
            desiredPoints[i,0] = desiredPoints[i,0] + XpolyCoeffs[k]*(interpolationArguments[i]**k)
            desiredPoints[i,1] = desiredPoints[i,1] + YpolyCoeffs[k]*(interpolationArguments[i]**k)
            desiredPoints[i,2] = desiredPoints[i,2] + ZpolyCoeffs[k]*(interpolationArguments[i]**k)
    if not returnCoefficients:
        return desiredPoints
    elif returnCoefficients:
        return desiredPoints, numpy.vstack( (XpolyCoeffs, YpolyCoeffs, ZpolyCoeffs) )
        
def PowerSeriesPolynomialInterpCoefficients3d(X, Y, Z):
    """ Interpolate a set of points in space using a power series polynomial expansion (finding coefficients for the N-1 order
    polynomial to interpolate the N data points). Use all the points (X, Y, Z) that are provided.
    Return the interpolating polynomials' coefficients.\n
    Parameters
    ----------
        @param X, Y, Z - coordinates of the know points that will be interpolated.\n
    Returns
    ----------
        @return - polynomial coefficients in the form of a numpy array (X, Y and Z interpolating polynomial's 
                    coefficients in rows 0, 1 and 2, respecitvely). The coefficients are given starting at 
                    the zeroth order term up to the one with the highest order."""

    """ Build the matrix of coefficients to solve the system of equations:
    [1, argumnet, argument^2, ..., argument^(N-1)] * [unknown polynomial coefficients] = [known data point values]"""
    numberOfPointsProvided = len(X)
    A = numpy.zeros( (numberOfPointsProvided, numberOfPointsProvided) )
    for j in range(numberOfPointsProvided):
        A[:,j] = numpy.linspace(0.0, 1.0, numberOfPointsProvided)**j # Points that have been provided must be evenly spaced through the interval.
    
    XpolyCoeffs=numpy.linalg.solve(A, X) # Coefficients of the interpolating polynomials in X, Y and Z
    YpolyCoeffs=numpy.linalg.solve(A, Y)
    ZpolyCoeffs=numpy.linalg.solve(A, Z)
    
    return numpy.vstack( (XpolyCoeffs, YpolyCoeffs, ZpolyCoeffs) )
    
def PowerSeriesPolynomialInterpGradientsCoefficients3d(X, Y, Z, Xdot, Ydot, Zdot, nodeEpochs):
    """ Interpolate a set of points in space using a power series polynomial expansion (finding coefficients for the 3rd order
    polynomial to interpolate the 2 data points when given the function values and derivatives at those points).
    Returns the polynomial coefficients.\n
    Args:
        @param X, Y, Z - coordinates of the know points that will be interpolated.\n
        @param Xdot, Ydot, Zdot - corresponding rates of change, i.e. velocities.\n
        
        @param nodeEpochs - JDAY epochs corresponding to X, Y and Z as well as Xdot, Ydot and Zdot.\n
    Returns:
        @return - by default noInterpolationPoints X 3 numpy.array containing X, Y and Z coordinates (in columns) of the interpolated points
                    can also be a 2-tuple with the polynomial coefficients in the form of a numpy array (X, Y and Z interpolating polynomial's 
                    coefficients in rows 0, 1 and 2, respecitvely) as the second entry."""

    """ Build the matrix of coefficients to solve the system of equations:
    [1, argumnet, argument^2, ..., argument^(N-1)] * [unknown polynomial coefficients] = [known data point values]"""
    polynomialOrder = len(X)+len(Xdot)-1
    t=numpy.linspace(0., 1., len(X)); # Interpolate in interval [0,1]
    A = numpy.zeros( (polynomialOrder+1, polynomialOrder+1) );
        
    for j in range( int((polynomialOrder+1)/2) ): # First half of the rows corresponds to the function values.
        for i in range( int(polynomialOrder+1) ):
            A[j,i] = t[j]**i
        
    for j in range( int((polynomialOrder+1)/2) ): # The rest correpsonds to derivative values.
        for i in range( 1, polynomialOrder+1 ):
            A[int((polynomialOrder+1)/2)+j,i] = i*t[j]**(i-1)
    
    Xdot = (nodeEpochs[-1]-nodeEpochs[0])*Xdot*JDAY_IN_SECONDS # Change velocity to km/JDAY and non-dimensionalise.
    Ydot = (nodeEpochs[-1]-nodeEpochs[0])*Ydot*JDAY_IN_SECONDS
    Zdot = (nodeEpochs[-1]-nodeEpochs[0])*Zdot*JDAY_IN_SECONDS

    rhsX = numpy.hstack( [X, Xdot] )
    rhsY = numpy.hstack( [Y, Ydot] )
    rhsZ = numpy.hstack( [Z, Zdot] )
    
    XpolyCoeffs=numpy.linalg.solve(A, rhsX) # Coefficients of the interpolating polynomials in X, Y and Z
    YpolyCoeffs=numpy.linalg.solve(A, rhsY)
    ZpolyCoeffs=numpy.linalg.solve(A, rhsZ)

    return numpy.vstack( (XpolyCoeffs, YpolyCoeffs, ZpolyCoeffs) )

def PiecewiseCubicInterp3d(X, Y, Z, Xdots, Ydots, Zdots, noInterpolationPoints=500):
    """ Fit a piecewise cubic polymonial spline to nodes given the corresponding positions and gradients
    to interpolate the data in three dimensions.\n
    Parameters
    ----------
        @param X, Y, Z - coordinates of the nodes in the form of 1-D numpy.arrays.\n
        @param Xdots, Ydots, Zdots - gradients of the function at the nodes (must be in the same order as X, Y and Z).\n
        @param noInterpolationPoints - number of points to be returned that are evenly placed along the curve formed by the supplied points.\n
    Returns
    ----------
        @return - noInterpolationPoints X 3 numpy.array containing X, Y and Z coordinates (in columns) of the interpolated points.\n
    Raises
    ----------
        @throws ValueError when the desired number of points cannot be evenly split between the number of intervals supplied."""
    if noInterpolationPoints%(len(X)-1):
        raise ValueError("Requested {} points to be evenly spaced amongst {} intervals. This is not possible.".format(noInterpolationPoints,(len(X)-1)))
    
    sDesired = numpy.linspace(-1., 1., noInterpolationPoints/(len(X)-1)) # Number of points to be computed in each interval bound by two end points.
    
    desiredPoints = numpy.zeros( (noInterpolationPoints,3) ) # X, Y and Z interpolated along S.
    pointCounter = 0 # Index of the point that will be added to the results next.
    
    for intervalID in range(len(X)-1): # Go through each interval.
        " Find the interpolating polynomial in each interval. "
        aX =  0.25*X[intervalID] -0.25*X[intervalID+1] +0.25* Xdots[intervalID] +0.25* Xdots[intervalID+1]
        bX =                                           -0.25* Xdots[intervalID] +0.25* Xdots[intervalID+1]
        cX = -0.75*X[intervalID] +0.75*X[intervalID+1] -0.25* Xdots[intervalID] -0.25* Xdots[intervalID+1]
        dX =  0.5* X[intervalID] +0.5* X[intervalID+1] +0.25* Xdots[intervalID] -0.25* Xdots[intervalID+1]
        
        aY =  0.25*Y[intervalID] -0.25*Y[intervalID+1] +0.25* Ydots[intervalID] +0.25* Ydots[intervalID+1]
        bY =                                           -0.25* Ydots[intervalID] +0.25* Ydots[intervalID+1]
        cY = -0.75*Y[intervalID] +0.75*Y[intervalID+1] -0.25* Ydots[intervalID] -0.25* Ydots[intervalID+1]
        dY =  0.5* Y[intervalID] +0.5* Y[intervalID+1] +0.25* Ydots[intervalID] -0.25* Ydots[intervalID+1]

        aZ =  0.25*Z[intervalID] -0.25*Z[intervalID+1] +0.25* Zdots[intervalID] +0.25* Zdots[intervalID+1]
        bZ =                                           -0.25* Zdots[intervalID] +0.25* Zdots[intervalID+1]
        cZ = -0.75*Z[intervalID] +0.75*Z[intervalID+1] -0.25* Zdots[intervalID] -0.25* Zdots[intervalID+1]
        dZ =  0.5* Z[intervalID] +0.5* Z[intervalID+1] +0.25* Zdots[intervalID] -0.25* Zdots[intervalID+1]
        
        for j in range(len(sDesired)): # Find the values at all the points in this interval.
            desiredPoints[pointCounter,0] = aX*sDesired[j]*sDesired[j]*sDesired[j] + bX*sDesired[j]*sDesired[j] + cX*sDesired[j] + dX # Interpolate the values of all the coordinates for this point.
            desiredPoints[pointCounter,1] = aY*sDesired[j]*sDesired[j]*sDesired[j] + bY*sDesired[j]*sDesired[j] + cY*sDesired[j] + dY
            desiredPoints[pointCounter,2] = aZ*sDesired[j]*sDesired[j]*sDesired[j] + bZ*sDesired[j]*sDesired[j] + cZ*sDesired[j] + dZ
            pointCounter=pointCounter+1 # Mark the point as added.

    return desiredPoints

def dateTimeToJDAY(datetimeObj):
    """ Convert a datetime.datetime object into a Julian day.\n
    Prameters
    ----------
        @param datetime.datetime that has year, month, day, hour, minute and microsecond attributes.\n
    Returns
    ----------
        @return Julian day as a float corresponding to the input datetime."""
        
    return jday(datetimeObj.year, datetimeObj.month, datetimeObj.day, datetimeObj.hour, datetimeObj.minute, datetimeObj.second+datetimeObj.microsecond/1E6)

def JDAYtoDateTime(JDAY):
    """ Convert a Julian day into a datetime.datetime object.\n
    Prameters
    ----------
        @param Julian day as a float corresponding to the output datetime.\n
    Returns
    ----------
        @return datetime.datetime that has year, month, day, hour, minute and microsecond attributes."""
    year,month,day,hour,minute,second = invjday(JDAY) # This will return second as a decimal number. Want microseconds for datetime.
    microsecond, second = math.modf(second) # Separate the integer and decimal parts.
    return datetime.datetime(year,month,day,hour,minute,int(second), microsecond=int(microsecond*1E6))
    
def generateObjectsFrom3LE(filename=threeLE_FILE_NAME, defaultRadiusRB=1.7691, defaultRadiusPL=1.0350, defaultRadiusDEB=0.1558, defaultRadiusOther=0.3470, BstarMultiplier=1.0):
    """ Read three-line element data (3LE) in the NORAD format (standard, elements are mean and time given in UTC).
    Add a hard body radius to the objects if it can be loceted in the stkRadiusFile.rad and, if not, a default value.\n
    Parameters
    ----------
    @param filename - name of the filename (with extension) which contains 3LE data, line by line\n
    @param defaultRadiusRB - default radius to be set to all the satellites whos radius is not in the radii file and who have R/B in their name, m.\n
    @param defaultRadiusPL - default radius to be set to all the satellites whos radius is not in the radii file and who have P/L in their name, m.\n
    @param defaultRadiusDEB - default radius to be set to all the satellites whos radius is not in the radii file and who have DEB in their name, m.\n
    @param defaultRadiusOther - default radius to be set to all the satellites whos radius is not in the radii file and who don't have R/B, P/L or DEB in their name, m.\n
    @param BstarMultiplier - factor by which to multiply the B* coefficient. Can be used in order to simulate solar activity variations.\n
    Returns
    ----------
    @return - dictionary of SpaceObjects sorted according to their NORAD ID. """
    with open(filename, 'r') as tleFile:
        tleLines=tleFile.readlines()  
        customPrint('Imported 3LEs for '+str(int(len(tleLines)/3.0))+' objects.', GREEN)
    
    with open('stkRadiusFile.rad', 'r') as radiiFile: # Contains hard body radii for some objects.
        tempLines=radiiFile.readlines()       
    
    lines = [line for line in tempLines if not line.startswith('#') and not line.startswith('stk')] # Remove the surplus lines.
    satRadii={}
    for radiusData in lines:
        satRadii[radiusData.split()[0]]=float(radiusData.split()[1])
        
    satellites={}
    linesRead=0
    for tleLine in tleLines:
        if tleLine.startswith('0'):
            line0 = tleLine; linesRead+=1;
        elif tleLine.startswith('1'):
            line1=tleLine; linesRead+=1;
        elif tleLine.startswith('2'):
            line2=tleLine; linesRead+=1;
        if linesRead==3: # Read the whole 3LE.
            tempSat = SpaceObject(line0+line1+line2, bstarModifier=BstarMultiplier) # Create the satellite.
            try:
                tempSat.setRadius(satRadii[str(tempSat.NORAD_ID)]) # Radii are stored according to the NORAD catalogue number.
            except KeyError:
                if 'R/B' in line0:
                    tempSat.setRadius(defaultRadiusRB)
                elif 'P/L' in line0:
                    tempSat.setRadius(defaultRadiusPL)
                elif 'DEB' in line0:
                    tempSat.setRadius(defaultRadiusDEB)
                else:
                    tempSat.setRadius(defaultRadiusOther)
#                customPrint('No radius data exists for {}: {}, using a default value of {} m.'.format(tempSat.SGP4_SATELLITE.satnum, line0, tempSat.HARD_BODY_RADIUS*1000))
                
            satellites[str(tempSat.NORAD_ID)] = tempSat # Use the in-built functionality of sgp4 to extract the satnum (NORAD ID) that will be the key to this object.
            linesRead=0 # Reset the counter.
    return satellites
    
def filterObjectsByPerigeeRadius(satelliteDictionary, radius):
    """ Find the objects from the supplied satelliteDictionary which have perigee radius larger than the specified upper limit and
    remove these objects from the initial dictionary at the self.CURRENT_EPOCH_JDAY epoch (TLE epoch by default).\n
    Parameters
    ----------
    @param satelliteDisctionary- dictionary of SpaceObjects objects as returned by @see generateObjectsFrom3LE.
        These objects need to be propagataed to contain CURRENT_POS and CURRENT_V attributes - these will be 
        the position and velocity at the TLE epoch after object initialisation.\n
    @param radius - perigee radius below which objects will be kept."""
    for satID in satelliteDictionary.keys():
        if satelliteDictionary[satID].calculatePerigeeRadiusWGS84() >= radius:
            del satelliteDictionary[satID]
            
def filterObjectsByApogeeRadius(satelliteDictionary, radius):
    """ Find the objects from the supplied satelliteDictionary which have apogee radius smalle than the specified upper limit and
    remove these objects from the initial dictionary at the self.CURRENT_EPOCH_JDAY epoch (TLE epoch by default).\n
    Parameters
    ----------
    @param satelliteDisctionary- dictionary of SpaceObjects objects as returned by @see generateObjectsFrom3LE.
        These objects need to be propagataed to contain CURRENT_POS and CURRENT_V attributes - these will be 
        the position and velocity at the TLE epoch after object initialisation.\n
    @param radius - apogee radius above which objects will be kept."""
    for satID in satelliteDictionary.keys():
        if satelliteDictionary[satID].calculateApogeeRadiusWGS84() <= radius:
            del satelliteDictionary[satID]
 
class SpaceObject(object):
    def __init__(self, threeLE, bstarModifier=1.0, noInterpolationPoints=2):
        """ An object which can be propagated in space using the SGP4 propagator. The propagator object can be accessed via 
        self.SGP4_SATELLITE. Its name (1st 3LE line) is saved under self.SAT_NAME. The initial cartesian position and velocity 
        in TEME frame are set as the CURRENT_POS and CURRENT_V attributes. Epoch of the TLE, from which the object was created,
        is kept as self.TLE_EPOCH attribute in form of a tuple as returned by @see sgp4.ext.invjday. The same epoch is set as
        the current epoch self.CURRENT_EPOCH_JDAY (in julian days). TLE epoch is also available in datetime.datetime forat under TLE_EPOCH_DATETIME.
        NORAD_ID is the NORAD catalogue number of the given object including the leading 0, if it's present.\n
        \n
        The object also supports interpolation of its position and velocity - it has attributes self.CURRENT_POSES and
        self.CURRENT_VELOS that are numpy.Arrays of length noInterpolationPoints containing position and velocity data
        in cartesian space (size (noInterpolationPoints,3) ).\n
        Parameters
        ----------
            @param threeLE - NORAD 3 line element used to initialise the satellite object.\n
            @param bstarModifier - a factor by which the B* coefficient of a satellite will be multiplied by, used to simulate atmospheric density changes.\n"""
                
        for tleLine in threeLE.split('\n'):
            if not tleLine.startswith('1') and not tleLine.startswith('2'): # It's easy to forget to add 0 in some TLEs.
                if tleLine.startswith('0'):
                    line0 = tleLine[2:] # Don't want the 0 and spacebar after it in the SAT_NAME.
                else:
                    line0 = tleLine
                    
            elif tleLine.startswith('1'):
                line1=tleLine
                self.NORAD_ID=tleLine[2:7] # Include the leading 0, if present.
            elif tleLine.startswith('2'):
                line2=tleLine
                self.SGP4_SATELLITE = twoline2rv(line1, line2, wgs84, bstarModifier) # Append when read the second line of the TLE - will always come after 1st so next line will be a new object.
                self.SAT_NAME = line0.rstrip(); # Remove the trailing whitespaces.
                self.TLE_EPOCH = invjday(self.SGP4_SATELLITE.jdsatepoch)
                self.TLE_EPOCH_JDAY = self.SGP4_SATELLITE.jdsatepoch
                self.CURRENT_EPOCH_JDAY = self.SGP4_SATELLITE.jdsatepoch
                self.TLE_EPOCH_DATETIME = datetime.datetime( self.TLE_EPOCH[0],self.TLE_EPOCH[1],self.TLE_EPOCH[2],self.TLE_EPOCH[3],self.TLE_EPOCH[4],int(self.TLE_EPOCH[5]),microsecond=int((self.TLE_EPOCH[5]-int(self.TLE_EPOCH[5]))*1E6) )

                self.CURRENT_POS, self.CURRENT_V = self.SGP4_SATELLITE.propagate(self.TLE_EPOCH[0],self.TLE_EPOCH[1],self.TLE_EPOCH[2],self.TLE_EPOCH[3],self.TLE_EPOCH[4],self.TLE_EPOCH[5])
                
                self.HARD_BODY_RADIUS = 5.0/1000. # Default values.
                self.CURRENT_PERIGEE_RADIUS = 0.0
                self.CURRENT_APOGEE_RADIUS = 0.0
                
                self.CURRENT_POSES = numpy.zeros( (noInterpolationPoints,3) )
                self.CURRENT_VELOS = numpy.zeros( (noInterpolationPoints,3) )
    
    def __str__(self):
        """ Defines how the object will be printed. """
        return self.NORAD_ID+" "+self.SAT_NAME+".\n\tRadius (km): "+str(self.HARD_BODY_RADIUS)+".\n\tTLE epoch (yyyy, MM, dd, hh, m, s.0):\n\t"+str(self.TLE_EPOCH)
    
    def propagate(self, year, month, day, hour, minute, second, updatePositionAndVelocity=False):
        """ Propagate the object to the given UTC time. If desired save its current position and volocity in the
        True Equator Mean Equinox frame true of date under self.CURRENT_POS and self.CURRENT_V.\n
        Parameters
        ----------
            @param updatePositionAndVelocity - whether to save the positions and velocities in object's attributes.\n
        Returns
        ----------
            @return 2-tuple of 3-lists with: currentPosition, current velocity."""
        p,v=self.SGP4_SATELLITE.propagate(year, month, day, hour, minute, second)
        if updatePositionAndVelocity:
            self.CURRENT_POS=p; self.CURRENT_V = v; self.CURRENT_EPOCH_JDAY=jday(year, month, day, hour, minute, second)
        return p,v
        
    def propagateDatetime(self, datetimeObj, updatePositionAndVelocity=False):
        """ Propagate the object to given UTC time. If deisred its current position and volocity in the
        True Equator Mean Equinox frame true of date under self.CURRENT_POS and self.CURRENT_V.\n
        Parameters
        ----------
            @param datetimeObj - datetime.datetime to which object is to be propagated.\n
        Returns
        ----------
            @return 2-tuple of 3-lists with: currentPosition, current velocity."""
        p,v=self.SGP4_SATELLITE.propagate(datetimeObj.year,datetimeObj.month,datetimeObj.day,datetimeObj.hour,datetimeObj.minute,datetimeObj.second+datetimeObj.microsecond/1E6)
        if updatePositionAndVelocity:
            self.CURRENT_POS=p; self.CURRENT_V = v; self.CURRENT_EPOCH_JDAY=jday(datetimeObj.year, datetimeObj.month, datetimeObj.day, datetimeObj.hour, datetimeObj.minute, datetimeObj.second+datetimeObj.microsecond/1E6)
        return p,v
        
    def propagateJDAY(self, JDAY, updatePositionAndVelocity=False):
        """ Propagate the object to given UTC time given as a Julian Day. If desired and save object's
        current position and volocity in the True Equator Mean Equinox frame true of date under 
        self.CURRENT_POS and self.CURRENT_V.\n
        Parameters
        ----------
            @param datetimeObj - datetime.datetime to which object is to be propagated.\n
        Returns
        ----------
            @return 2-tuple of 3-lists with: currentPosition, current velocity."""
        p,v=self.SGP4_SATELLITE.propagateJDAY( JDAY )
        if updatePositionAndVelocity:
            year,month,day,hour,minute,second = invjday(JDAY)
            self.CURRENT_POS=p; self.CURRENT_V = v; self.CURRENT_EPOCH_JDAY=jday(year, month, day, hour, minute, second)
        return p,v
        
    def setRadius(self, radius=5.0):
        """ Set the radius of the object which will be used for assessing collison probability assessments.
        This radius can be accessed via self.HARD_BODY_RADIUS. A default value of 5.0 m is assumed to guarantee
        compatibility with STG Collision Assessment Toll (CAT) and SOCRATES data.\n
        Parameters
        ----------
            @param radius of the object in metres."""
        self.HARD_BODY_RADIUS = radius/1000.
        
    def findRadius(self, defaultRadius=5.0, radiusFile="stkRadiusFile.rad"):
        """ Set the radius of the object which will be used for assessing collison probability assessments.
        Use the object database file to find the radius of the object. It the radius does not exist in the
        database for the given object a default value of 5.0 m is assumed to guarantee
        compatibility with STG Collision Assessment Toll (CAT) and SOCRATES data.\
        This radius can be accessed via self.HARD_BODY_RADIUS (in km).\n
        Parameters
        ----------
            @param defaultRadius of the object in metres."""
        with open(radiusFile, 'r') as radiiFile: # Contains hard body radii for some objects.
            tempLines=radiiFile.readlines()       
    
        lines = [line for line in tempLines if not line.startswith('#') and not line.startswith('stk')] # Remove the surplus lines.
        satRadii={}
        for radiusData in lines:
            satRadii[radiusData.split()[0]]=radiusData.split()[1]
        
        try:
            self.HARD_BODY_RADIUS=float(satRadii[self.NORAD_ID])/1000. # Radii are stored according to the NORAD catalogue number.
        except KeyError:
            customPrint('No radius data exists for {}, using a default value of {} m.'.format(self.NORAD_ID,defaultRadius), RED)
            self.HARD_BODY_RADIUS = defaultRadius/1000.

    def calculatePerigeeRadiusWGS84(self):
        """ Calculate the perigee radius of the given object w.r.t. WGS84 geoid. Store it under self.CURRENT_PERIGEE_RADIUS in km.
        Use the self.CURRENT_POS and self.CURRRENT_V (position and velocity at self.CURRENT_EPOCH_JDAY) for calculation.\n
        Returns
        ----------
            @return the perigee radius in km."""
        p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = rv2coe(self.CURRENT_POS, self.CURRENT_V, GRAVITY_CONSTANT)
        del p, incl, omega, argp, nu, m, arglat, truelon, lonper;
        
        self.CURRENT_PERIGEE_RADIUS = a*(1.-ecc)
        return a*(1.-ecc)
                
    def calculateApogeeRadiusWGS84(self):
        """ Calculate the apogee radius of the given object w.r.t. WGS84 geoid.  Store it under self.CURRENT_APOGEE_RADIUS in km.\n
        Use the self.CURRENT_POS and self.CURRRENT_V (position and velocity at self.CURRENT_EPOCH_JDAY) for calculation.
        Returns
        ----------
            @return the apogee radius in km."""
        p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = rv2coe(self.CURRENT_POS, self.CURRENT_V, GRAVITY_CONSTANT)
        del p, incl, omega, argp, nu, m, arglat, truelon, lonper;
        
        self.CURRENT_APOGEE_RADIUS = a*(1.+ecc)
        return a*(1.+ecc)

class Conjunction(object):
    __doc__="""A class that is used to stroe data about a particular orbital conjunction."""

    def __init__(self, primaryNORADID, secondaryNORADID, epoch, missDistance, maxProbability, trueProbability=0.0):
        """ A class that stores information about a particular orbital conjunction.\n
        Parameters
        ----------
            @param epoch - time of closest approach in JDAYs.\n
            @param missDistance - closest approach distance in km.\n
            @param maxProbability - maximum collision probability between the objects during this conjunction.\n
            @param trueProbability - the true probability that the objects will collide.\n
        """
        self.PRIMARY_NORAD_ID = primaryNORADID
        self.SECONDARY_NORAD_ID = secondaryNORADID
        self.TCA = JDAYtoDateTime( epoch )
        self.MISS_DISTANCE = missDistance
        self.MAX_COLLISION_PROBABILITY = maxProbability
        self.TRUE_COLLISION_PROBABILITY = trueProbability
        
    def __str__(self):
        return "Conjunction between {} and {} at {} UTCG.\nMiss distance (km): {}\nMax collision probability: {}".format(self.PRIMARY_NORAD_ID, self.SECONDARY_NORAD_ID, self.TCA, self.MISS_DISTANCE, self.MAX_COLLISION_PROBABILITY)

import sys
import numpy as np
import sympy as sp
import numpy.linalg as lg   
import math as m
import SupportingFiles.engineBuilder as eb

##################################################################
######################## Flight Mechanics ########################
##################################################################
def cordTransform( vector, phi, theta, psi, toSystem ):
    phiRad = np.deg2rad(phi)
    thetaRad = np.deg2rad(theta)
    psiRad = np.deg2rad(psi)
    
    R1 = np.array([[1, 0, 0], 
                   [0, np.cos(phiRad), np.sin(phiRad)], 
                   [0, -np.sin(phiRad), np.cos(phiRad)]])
    
    R2 = np.array([[np.cos(thetaRad), 0, -np.sin(thetaRad)], 
                   [0, 1, 0], 
                   [np.sin(thetaRad), 0, np.cos(thetaRad)]])
    
    R3 = np.array([[np.cos(psiRad), np.sin(psiRad), 0], 
                   [-np.sin(psiRad), np.cos(psiRad), 0], 
                   [0, 0, 1]])
    
    transformedVector = 0

    if toSystem.lower() == "body": transformedVector = R1@R2@R3@vector
    else: transformedVector = R3.T@R2.T@R1.T@vector

    return transformedVector.reshape(-1, 1)


def inertialAcceleration( v, vDot, omega ):
    return np.array([[vDot[0] + omega[1]*v[2] - omega[2]*v[1]], 
                  [vDot[1] + omega[2]*v[0] - omega[0]*v[2]], 
                  [vDot[2] + omega[0]*v[1] - omega[1]*v[0]]])


##################################################################
#################### Air-Breathing Propulsion ####################
##################################################################

def engineCalculations(gammaA=1.4, gammaG=1.333, cpa=1.005, cpg=1.148):
    station = 0
    selection = 1
    tOut = 0
    pOut = 0

    pa = float(input("Enter ambient pressure: "))
    Ta = float(input("Enter ambient temp: "))
    mEta = float(input("Enter Mechanical Efficiency: "))

    engine = eb.Engine(mEta)

    pIn = pa
    tIn = Ta

    while selection != 0:
        station += 1

        print("\nSelect a station from the list")
        selection = int(input("(1) Compressor\n(2) Combuster\n(3) Turbine\n(0) Stop\n"))

        # Compressor
        if selection == 1:
            comp = eb.Compressor(pIn, tIn, station, gammaA, cpa, engine)
            pOut = comp.getPout()
            tOut = comp.getTout()

        # Combuster
        # TODO: assuing no heat exchanger
        elif selection == 2:
            comb = eb.Combustor(pIn, pOut, station, engine)
            pOut = comb.getPout()
            tOut = comb.getTout()

        #Turbine
        elif selection == 3:
            turb = eb.Turbine(pIn, tIn, station, gammaG, cpg, engine)
            pOut = turb.getPout()
            tOut = turb.getTout()

        pIn = pOut
        tIn = tOut

    print("\n")
    print(engine)

#TODO: DOES NOT WORK
def machFromAreaRatio(A, Astar, gamma=1.4):
    sol = A/Astar
    M = sp.symbols('M')

    equation = sp.Eq(1/M*((2/(gamma+1))*(1+(gamma-1)/2*M**2))**((gamma+1)/(2*(gamma-1))), sol)

    return sp.solve(equation)

def compIsoEfficiency(T01, T02, pRatio, gamma=1.4):
    exp = (gamma-1)/gamma
    eta = T01*(pRatio**(gamma) - 1)/(T02 - T01)
    return eta

def compPolyEfficiency(pRatio, tRatio, gamma=1.4):
    exp = (gamma-1)/gamma
    eta = (np.log(pRatio)**exp)/np.log(tRatio)
    return eta

def inletTstag(T, Ca, cp):
    return T+(Ca**2/(2*cp))

def inletPratio(T, Ca, cp, gamma, etaI):
    exp = gamma/(gamma-1)
    pRatio = (1+(etaI*Ca**2/(2*cp*T)))**exp
    
    return pRatio

def compTstag(T, eta, pRatio, gamma):
    exp = (gamma-1)/gamma
    delT = (T/eta)*((pRatio**exp)-1)
    tOut = delT+T

    return tOut, delT

def compPolyTstag(T, etaC, pRatio, gamma):
    exp = (gamma-1)/(etaC*gamma)
    tRatio = pRatio**exp
    tOut = T*tRatio
    delT = tOut - T

    return tOut, delT

def compWork(cp, delT, etaM):
    return cp*delT/etaM

def turbTstag(T, etaT, pRatio, gamma):
    exp = (gamma-1)/gamma
    delT = etaT*T*(1-((1/pRatio)**exp))
    tOut = T-delT

    return tOut, delT

def turbPolyTstag(T, etaT, pRatio, gamma):
    exp = etaT*(gamma-1)/gamma
    tRatio = pRatio**exp
    tOut = T/tRatio
    delT = T - tOut

    return tOut, delT

def turbPratio(T, delT, etaT, gamma):
    return (etaT*T/(etaT*T-delT))**(gamma/(gamma-1))

def turbWork(cp, delT, etaM):
    return cp*delT*etaM

##################################################################
########################## Gas Dynamics ##########################
##################################################################

def linterp( y1, y2, x1, x2, x ):
    y = y1 + (y2-y1)/(x2-x1)*(x-x1)
    return y

def massFluxM( pStag, R, TStag, M, gamma=1.4 ):
    flux = (pStag/np.sqrt(R*TStag))*(np.sqrt(gamma)*M/((1+(gamma-1)/2*M**2)**((gamma+1)/(2*(gamma-1)))))
    return flux

def maxMassFlow( AStar, pStag, R, TStag, gamma=1.4 ):
    flow = AStar*pStag*np.sqrt(gamma)*(1+((gamma-1)/2))**((gamma+1)/(2*(gamma-1)))/np.sqrt(R*TStag)
    return flow

def areaRatio( M, gamma=1.4 ):
    ratio = (1/M)*((1+((gamma-1)/2)*M**2)/((gamma+1)/2))**((gamma+1)/(2*(gamma-1)))
    return ratio

def stagT_T( M, gamma=1.4 ):
    ratio = 1 + (gamma-1)/2*M**2
    return ratio

def stagP_P( M, gamma=1.4 ):
    tRatio = stagT_T(M, gamma)
    pRatio = tRatio**(gamma/(gamma-1))
    return pRatio

def T2_T1( M1, M2, gamma=1.4):
    ratio = (1 + (gamma-1)/2*M1**2) / (1 + (gamma-1)/2*M2**2)
    return ratio

def stagPshock( M1, gamma=1.4 ):
    pRatio = ((M1**2*(gamma+1)/2)/(1+M1**2*(gamma-1)/2))**(gamma/(gamma-1))*(2*gamma*M1**2/(gamma+1)-((gamma-1)/(gamma+1)))**(1/(1-gamma))
    return pRatio

def ObliqueM2( M1, theta, beta, gamma=1.4 ):
    thetaRad = np.deg2rad(theta)
    betaRad = np.deg2rad(beta)
    M2 = np.sqrt((M1**2*np.sin(betaRad)**2 + 2/(gamma-1))/(np.sin(betaRad-thetaRad)**2*(2*gamma/(gamma-1)*M1**2*np.sin(betaRad)**2 - 1)))
    return M2

def ObliqueVnRelation( M1, beta, gamma=1.4 ):
    betaRad = np.deg2rad(beta)
    v1n_v2n = ((gamma+1)*M1**2*np.sin(betaRad)**2) / ((gamma-1)*M1**2*np.sin(betaRad)**2 + 2)
    return v1n_v2n

def ObliqueTratio( M1, beta, gamma=1.4 ):
    betaRad = np.deg2rad(beta)
    tr = (1 + (gamma-1)/2*M1**2*np.sin(betaRad)**2) * (2*gamma/(gamma-1)*M1**2*np.sin(betaRad)**2 - 1) / (M1**2*np.sin(betaRad)**2*(gamma+1)**2 / (2*(gamma-1))) 
    return tr

def ObliquePratio( M1, beta, gamma=1.4 ):
    betaRad = np.deg2rad(beta)
    pr = (2*gamma*M1**2*((np.sin(betaRad))**2)/(gamma+1)) - ((gamma-1)/(gamma+1))
    return pr

def ObliqueStagPratio( M1, beta, gamma=1.4 ):
    betaRad = np.deg2rad(beta)
    pr = ObliquePratio( M1, beta )
    prStag = pr**(1/(1-gamma)) * ((gamma+1)/2*M1**2*np.sin(betaRad)**2 / (1 + (gamma-1)/2*M1**2*np.sin(betaRad)**2))**(gamma/(gamma-1))
    return prStag

##################################################################
######################### Aero Computing #########################
##################################################################

# Error catch 2
# Makes sure correct number of value are entered
def ec2(nvalues, input_string, error_statement):
    import re

    while True:
        val = str(input(input_string))

        val = re.findall(r'[+]?\d*\.?\d+|[-+]?\d+', val)

        if len(val) == nvalues:
            print('Counted ' + str(len(val)) +' floats')
            break
        else:
            print(error_statement)
    val = [float(i) for i in val]

    return val

# Root Search
# Searches for a zero of a function f
# Uses interval (a, b) w/ increments dx
# Returns bounds (x1, x2) --> if no roots x1 = x2 = None
# Can be called again with returned bounds
def rootsearch(f,a,b,dx):
    from numpy import sign
    
    x1 = a; f1 = f(a)
    x2 = a + dx; f2 = f(x2)
    
    while sign(f1) == sign(f2):
        if x1 >= b: return None,None
        
        x1 = x2; f1 = f2
        x2 = x1 + dx; f2 = f(x2)
    
    else:
        return x1,x2
    
# Newton-Raphson (safe)
# Finds root of f = 0 using f and df (derivative)
# Derivitave can be found with sympy (check cherry tree)
# Assumes root to find is in (a, b)
# Midpoint of (a, b) is used as initial guess
# Brackets are changed with each iteration
# If an iteration is outside brackets 
# --> it is disreguarded and bisection is used
def newtonRaphson(f,df,a,b,tol=1.0e-9):
    import SupportingFiles.error as error
    from numpy import sign
    
    fa = f(a)
    
    if fa == 0.0: return a
    
    fb = f(b)
    
    if fb == 0.0: return b
    
    if sign(fa) == sign(fb): error.err("Root is not bracketed")
    
    x = 0.5*(a + b)
    
    for i in range(30):
        fx = f(x)
        
        if fx == 0.0: return x
        
        # Tighten the brackets on the root
        if sign(fa) != sign(fx): b = x
        
        else: a = x
        
        # Try a Newton-Raphson step
        dfx = df(x)
        
        # If division by zero, push x out of bounds
        try: dx = -fx/dfx
        except ZeroDivisionError: dx = b - a
        
        x = x + dx
        
        # If the result is outside the brackets, use bisection
        if (b - x)*(x - a) < 0.0:
            dx = 0.5*(b - a)
            x = a + dx
        
        # Check for convergence
        if abs(dx) < tol*max(abs(b),1.0): return x
    print("Too many iterations in Newton-Raphson")

# Theta Beta M
# Sets variables
def thetaBetaM():
    import numpy as np
    import math
    import sympy as sp
    c = 0
    r = 0

    while c == 0:
        # Gathers what needs to be found
        find = input("Do you need to find (deflection angle(T), shock-wave angle(B) or mock number(M)): ")

        # Comments are the same for each if statement
        if find == "T":
            # Gathers what the known variables are
            M, beta = ec2(2, "Input mock number then shock-wave angle: ", "Please input 2 numbers: ")
            
            # Creates a symbolic variable and changes degrees to radians
            beta = math.radians(beta)
            theta = sp.symbols('theta')

            # Creates equation and finds derivative
            f = ((2/sp.tan(beta))*((M**2 * (sp.sin(beta))**2 - 1) / (M**2 * (1.4 + sp.cos(2*beta)) + 2))) - sp.tan(theta)
            df = sp.diff(f, theta)

            # turns equations into usable functions
            newF = sp.lambdify(theta, f)
            newDf = sp.lambdify(theta, df)

            # Computes answer
            x1, x2 = rootsearch(newF, 0, 1, 0.1)
            r = newtonRaphson(newF, newDf, x1, x2)

            #Changes radians back to degrees
            r = math.degrees(r)
            beta = math.degrees(beta)

            # Prints solution 
            print("The deflection angle with a mach number of {:.1f} and shock-wave = {:.1f}° is {:.1f}°".format(M, beta, r))

            # Changes variable to exit while loop
            c = 1
        elif find == "B":
            M, theta = ec2(2, "Input mock number then deflection angle: ", "Please input 2 numbers: ")
            
            theta = math.radians(theta)
            beta = sp.symbols('beta')

            f = ((2/sp.tan(beta))*((M**2 * (sp.sin(beta))**2 - 1) / (M**2 * (1.4 + sp.cos(2*beta)) + 2))) - sp.tan(theta)
            df = sp.diff(f, beta)

            newF = sp.lambdify(beta, f)
            newDf = sp.lambdify(beta, df)

            x1, x2 = rootsearch(newF, 0.001, 2, 0.1)
            r = newtonRaphson(newF, newDf, x1, x2)

            r = math.degrees(r)
            theta = math.degrees(theta)

            print("The shock-wave angle with a mach number of {:.1f} and deflection = {:.1f}° is {:.1f}°".format(M, theta, r))

            c = 1
        elif find == "M":
            beta, theta = ec2(2, "Input shock-wave angle then deflection angle: ", "Please input 2 numbers: ")
            
            theta = math.radians(theta)
            beta = math.radians(beta)
            M = sp.symbols('M')

            f = ((2/sp.tan(beta))*((M**2 * (sp.sin(beta))**2 - 1) / (M**2 * (1.4 + sp.cos(2*beta)) + 2))) - sp.tan(theta)
            df = sp.diff(f, M)

            newF = sp.lambdify(M, f)
            newDf = sp.lambdify(M, df)

            x1, x2 = rootsearch(newF, 1, 11, 1)
            r = newtonRaphson(newF, newDf, x1, x2)
            
            theta = math.degrees(theta)
            beta = math.degrees(beta)

            print("The mach number with shock-wave = {:.1f}° and deflection = {:.1f}° is {:.1f}".format(beta, theta, r))

            c = 1
        else:
            print("Please enter either T, B, or M")

        return r
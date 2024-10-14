import AeroEqs as ae
import numpy as np

class Engine: 

    #variables
    components = []
    compressorWork = 0
    turbineWork = 0
    Ca = 0
    mTotal = 0
    mCold = 0
    mHot = 0

    def __init__(self, etaM, pa, Ta, startStation):
        self.etaM = etaM
        self.pa = pa
        self.Ta = Ta
        self.station = startStation

    #printing
    def __str__(self) -> str:

        for i in self.components: print(i)

        station = self.components[-1].getStation() + 1
        pressure = self.components[-1].getPout()
        temp = self.components[-1].getTout()

        print(f"({station}): P={pressure} T={temp}")
        return("")

    #methods
    def addComponent(self, component):
        self.station += 1
        self.components.append(component)

    #set
    def setCompressorWork(self, work):
        if self.compressorWork == 0: self.compressorWork = work
        else: self.compressorWork += work

    def setTurbineWork(self, work):
        if self.turbineWork == 0: self.turbineWork = work
        else: self.turbineWork += work

        if self.compressorWork != 0: self.totalWork = work - self.compressorWork

    def setCa(self, Ca): self.Ca = Ca

    def setMassFlow(self, m, BPR = 0):
        self.mTotal = m
        
        if BPR != 0:
            self.mCold = m*BPR/(BPR+1)
            self.mHot = m/(BPR+1)

    #get
    def getEtaM(self): return self.etaM
    def getWorkC(self): return self.compressorWork
    def getTa(self): return self.Ta
    def getPa(self): return self.pa
    def getStation(self): return self.station
    def getCa(self): return self.Ca


class Inlet:
    
    #Constructor
    def __init__(self, gamma, cp, engine: Engine):
        #set variables from params
        self.gamma = gamma
        self.cp = cp
        self.stationStart = engine.getStation()

        #get speed (engine) and efficiency (inlet)
        #assuming inlet means airplane
        engine.setCa(float(input("\nEnter aircraft speed (m/s): ")))
        self.eta = float(input("Enter Intake Isentropic Efficiency: "))

        #set p and t in
        self.pIn = engine.getPa()
        self.tIn = engine.getTa()

        #calculate T and P
        pRatio = self.inletPratio(engine.getCa())
        self.pOut = pRatio*engine.getPa()
        self.tOut = self.inletTstag(engine.getCa())

        print(f'\nP out = {self.pOut}')
        print(f'Temp Out = {self.tOut}')
        print("\n")

    #Printing
    def __str__(self) -> str:
        print(f"({self.stationStart}): P={self.pIn} T={self.tIn}")
        print("    ()                  ()")
        print("    ||                  ||")
        print("    ||                  ||")
        print("    ||                  ||")
        print("    ()                  ()")
        return("")

    # Get
    def getTout(self): return self.tOut
    def getPout(self): return self.pOut
    def getStation(self): return self.stationStart

    #math
    def inletTstag(self, Ca):
        #assuming Ca = m/s
        T = self.tIn
        cp = self.cp

        Tstag = T+(Ca**2/(2*cp*1000))

        return Tstag

    def inletPratio(self, Ca):
        T = self.tIn
        cp = self.cp
        etaI = self.eta
        gamma = self.gamma

        exp = gamma/(gamma-1)
        pRatio = (1+(etaI*Ca**2/(2*cp*T*1000)))**exp
        
        return pRatio
    

class Fan:

    #attributes
    compPolyEta = 0

    #Constructor
    def __init__(self, pIn, tIn, gamma, cp, engine: Engine):
        #setting attributes
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = engine.getStation()

        self.compPolyEta = float(input("Enter the Compressor Polytropic Efficiency: "))
        self.FPR = float(input("Enter Fan Pressure Ratio: "))
        self.BPR = float(input("Enter Engine BPR: "))

        self.pOut = self.fanPstag()
        self.tOut = self.fanTstag()

        

        #TODO: add fan nozzle
        

    def __str__(self) -> str:
        #TODO: add station, p, T
        print(f"(): P= T=")
        print(r"  //   ____     ____    \\")
        print(r" ||   /    \   /    \    ||")
        print(r" ||  (-------X-------)   ||")
        print(r" ||   \____/   \____/    ||")
        print("                           ")
        print(r"//\\                    //\\")
        return("")
    
    def fanPstag(self): 
        pOut = self.FPR*self.pIn
        return pOut

    def fanTstag(self):
        gamma = self.gamma

        exp = (gamma-1)/(self.compPolyEta*gamma)
        tRatio = self.FPR**exp
        tOut = tRatio*self.tIn

        return tOut
        
    #TODO: add choke test
        

class Compressor:

    #TODO: consider moving more into functions
    # Constructor
    def __init__(self, pIn, tIn, gamma, cp, engine: Engine):
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = engine.getStation()

        efficiencyChoice = 0

        self.pRatio = float(input("\nEnter Compressor pressure ratio: "))

        print("\nWhat Efficiency is available?")
        while efficiencyChoice == 0:

            #TODO: still has input taken even if wrong choice is given
            efficiencyChoice = int(input("(1) Isentropic\n(2) Polytropic\n"))
            eta = float(input("\nEnter Compressor Efficiency: "))
            
            #isentropic
            if efficiencyChoice == 1: 
                self.etaIsen = eta
                self.tOut, self.delT = self.compTstag()
            #polytropic
            elif efficiencyChoice == 2: 
                self.etaPoly = eta
                self.tOut, self.delT = self.compPolyTstag()
            else: 
                efficiencyChoice = 0
                print("Select from the list given")

        self.work = self.compWork(engine.getEtaM())

        self.pOut = self.pRatio*pIn

        #adding work to engine
        engine.setCompressorWork(self.work)

        print(f'\nP out = {self.pOut}')
        print(f'Temp Out = {self.tOut}')
        print(f'Change in Temp = {self.delT}')
        print(f'Compressor Work = {self.work}')
        print("\n")

    def __str__(self) -> str:
        print(f"({self.stationStart}): P={self.pIn} T={self.tIn}")
        print("     ____________________")
        print("    | ------------------ |")
        print("    | ------------------ |")
        print(r"     \ ---------------- /")
        print(r"      \ -------------- /")
        print(r"       \ ____________ /")
        print("        |            |")
        return("")
    
    # Get
    def getWork(self): return self.work
    def getTout(self): return self.tOut
    def getPout(self): return self.pOut
    def getStation(self): return self.stationStart

    #TODO: consider not returning but setting directly
    #Math
    #TODO: not used
    def compIsoEfficiency(self, T01, T02, pRatio, gamma=1.4):
        exp = (gamma-1)/gamma
        eta = T01*(pRatio**(gamma) - 1)/(T02 - T01)
        return eta

    #TODO: not used
    def compPolyEfficiency(self, pRatio, tRatio, gamma=1.4):
        exp = (gamma-1)/gamma
        eta = (np.log(pRatio)**exp)/np.log(tRatio)
        return eta

    def compTstag(self):
        #assuming only used for isentropic
        T = self.tIn
        eta = self.etaIsen
        pRatio = self.pRatio
        gamma = self.gamma
        
        exp = (gamma-1)/gamma
        delT = (T/eta)*((pRatio**exp)-1)
        tOut = delT+T

        return tOut, delT

    def compPolyTstag(self):
        #assuming only used for polytropic
        T = self.tIn
        etaC = self.etaPoly
        pRatio = self.pRatio
        gamma = self.gamma

        exp = (gamma-1)/(etaC*gamma)
        tRatio = pRatio**exp
        tOut = T*tRatio
        delT = tOut - T

        return tOut, delT

    def compWork(self, etaM):
        cp = self.cp
        delT = self.delT
        
        return cp*delT/etaM


class Combustor:

    # Constructor
    def __init__(self, pIn, tIn, engine: Engine):
        self.pIn = pIn
        self.tIn = tIn
        self.stationStart = engine.getStation()

        #assuming pDrop ONLY --> not pDrop/pIn
        self.pDrop = float(input("\nEnter pressure drop across combustor: "))
        self.pOut = pIn - self.pDrop
        
        #TODO: assuming turbine always next
        self.tOut = float(input("Enter the Turbine Inlet Temp: "))

        print(f'\nP out = {self.pOut}\n')

    def __str__(self) -> str:
        print(f"({self.stationStart}): P={self.pIn} T={self.tIn}")
        print("        |_____||_____|")
        print(r"       //\    ||    /\\")
        print(r"      //  \   ||   /  \\")
        print(r"      \\  /   ||   \  //")
        print(r"       \\/____||____\//")
        print("        |     ||     |")
        return("")

    # Get
    def getTout(self): return self.tOut
    def getPout(self): return self.pOut
    def getStation(self): return self.stationStart


class Turbine:

    # Constructor
    def __init__(self, pIn, tIn, gamma, cp, engine: Engine):
        #initialize params
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = engine.getStation()

        #for checking
        efficiencyType = "Polytropic"
        pRatioCheck = True

        #Check for pressure ratio
        print("\nIs the turbine pressure ratio known?")
        choice = int(input("(1) yes\n(2) no\n"))

        #TODO: turbine pressure can be unknown but still use polytropic --> HW 2b

        #get pressure ratio
        if choice == 1: 
            #assuming always polytropic
            self.pRatio = float(input("\nEnter Turbine pressure ratio: "))
        
        #get work to find pressure ratio
        else:
            #assuming always isentropic
            cWork = engine.getWorkC()

            self.delT = cWork/cp
            self.tOut = tIn-self.delT

            efficiencyType = "isentropic"

            pRatioCheck = False

        eta = float(input(f"\nEnter Turbine {efficiencyType} Efficiency: "))

        #polytropic
        if pRatioCheck:
            self.etaPoly = eta
            self.tOut, self.delT = self.turbPolyTstag()

        #isentropic
        else:
            self.etaIsen = eta
            self.pRatio = self.turbPratio()

        # get pout and work
        self.pOut = pIn/self.pRatio
        self.work = self.turbWork(engine.getEtaM())

        #adding work to engine
        engine.setTurbineWork(self.work)

        print(f'\nP out = {self.pOut}')
        print(f'delta T = {self.delT}')
        print(f'Temp Out = {self.tOut}')
        print(f'Turbine Work = {self.work}\n')

    #printing
    def __str__(self) -> str:
        print(f"({self.stationStart}): P={self.pIn} T={self.tIn}")
        print("        |     ||     |")
        print(R"       /      ||      \ ")
        print(R"      / -------------- \ ")
        print(R"     /        ||        \ ")
        print(R"    / ------------------ \ ")
        print("    |         ||          |")
        return("")
    
    # Get
    def getWork(self): return self.work
    def getTout(self): return self.tOut
    def getPout(self): return self.pOut
    def getStation(self): return self.stationStart

    #Math
    #TODO: not used
    def turbTstag(self, T, etaT, pRatio, gamma):
        exp = (gamma-1)/gamma
        delT = etaT*T*(1-((1/pRatio)**exp))
        tOut = T-delT

        return tOut, delT

    def turbPolyTstag(self):
        T = self.tIn
        etaT = self.etaPoly
        pRatio = self.pRatio
        gamma = self.gamma
        
        exp = etaT*(gamma-1)/gamma
        tRatio = pRatio**exp
        tOut = T/tRatio
        delT = T - tOut

        return tOut, delT

    def turbPratio(self):
        #assuming only for isentropic
        T = self.tIn
        delT = self.delT
        etaT = self.etaIsen
        gamma = self.gamma
        
        pRatio = (etaT*T/(etaT*T-delT))**(gamma/(gamma-1))

        return pRatio

    def turbWork(self, etaM):
        cp = self.cp
        delT = self.delT

        w = cp*delT*etaM

        return w
    
class Nozzle:

    def __init__(self, pIn, tIn, gamma, cp, engine: Engine):
        #initialize params
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = engine.getStation()

        self.eta = float(input("\nEnter Nozzle Isentropic Efficiency: "))

        self.pRatio = pIn/engine.getPa()
        self.critPR = self.critPRatio()

        #choked
        if self.pRatio > self.critPR:
            print("\nNozzle is choked")
            self.R = float(input("Enter R: "))

            self.pOut = pIn*self.nozzlePratio()
            self.tOut = tIn*self.nozzleTratio()
            self.rhoOut = self.nozzleRho()
            #TODO: consider being owned by engine
            self.Cj = self.nozzleV()

        print(f'\nP out = {self.pOut}')
        print(f'Temp Out = {self.tOut}')
        print(f'Rho Out = {self.rhoOut}')
        print(f'Cj = {self.Cj}')

    #printing
    def __str__(self) -> str:
        print(f"({self.stationStart}): P={self.pIn} T={self.tIn}")
        print("    |                     |")
        print(r"     \                   /")
        print(r"      \                 /")
        print(r"       \               /")
        print("        @             @")
        return("")
    
    #Math
    def critPRatio(self):
        etaJ = self.eta
        gamma = self.gamma

        den = (1-((1/etaJ)*(gamma-1)/(gamma+1)))**(gamma/(gamma-1))
        pRatio = 1/den
        
        return pRatio
    
    def nozzlePratio(self):
        #Pout/Pin
        return (1/self.critPR)

    def nozzleTratio(self):
        #Tout/Tin
        return (2/(self.gamma+1))
    
    def nozzleRho(self):
        return (self.pOut/(self.R*self.tOut))
    
    def nozzleV(self):
        #assuming R needs to be /1000
        return np.sqrt(self.gamma*self.R*self.tOut*1000)
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
    etaJ = 0
    R = 0
    BPR = 0

    def __init__(self, gammaA, gammaG, cpa, cpg):
        self.gammaA = gammaA
        self.gammaG = gammaG
        self.cpa = cpa
        self.cpg = cpg

        self.pa = float(input("Enter ambient pressure: "))
        self.Ta = float(input("Enter ambient temp: "))
        self.etaM = float(input("Enter Mechanical Efficiency: "))
        self.station = int(input("Enter the first station number: "))

        self.pIn = self.pa
        self.tIn = self.Ta

    #printing
    def __str__(self) -> str:

        for i in self.components: print(i)

        station = self.components[-1].getStation() + 1
        pressure = self.components[-1].getPout()
        temp = self.components[-1].getTout()

        print(f"({station}): P={pressure} T={temp}")
        return("")

    #methods
    def addComponent(self, selection):
        gammaA = self.gammaA
        gammaG = self.gammaG
        cpa = self.cpa
        cpg = self.cpg
        pIn = self.pIn
        tIn = self.tIn

        # Inlet
        if selection == 1:
            component = Inlet(gammaA, cpa, self)

        #fan
        elif selection == 2:
            component = Fan(pIn, tIn, gammaA, cpa, self)
            pass

        # Compressor
        elif selection == 3:
            component = Compressor(pIn, tIn, gammaA, cpa, self)

        # Combuster
        # assuing no heat exchanger
        elif selection == 4:
            component = Combustor(pIn, tIn, self)

        #Turbine
        elif selection == 5:
            tIn = float(input("\nEnter the Turbine Inlet Temp: "))
            component = Turbine(pIn, tIn, gammaG, cpg, self)

        elif selection == 6: 
            component = Nozzle(pIn, tIn, gammaG, cpg, self)

        self.components.append(component)
        self.pIn = component.pOut
        self.tIn = component.tOut

        self.station += 1

        #update fan nozzle exit station number
        if any(isinstance(component, Fan) for component in self.components): component.nozzle.stationStart = self.station+1

    #post assembly


    #set
    def setCompressorWork(self, work):
        if self.compressorWork == 0: self.compressorWork = work
        else: self.compressorWork += work

    def setTurbineWork(self, work):
        if self.turbineWork == 0: self.turbineWork = work
        else: self.turbineWork += work

        if self.compressorWork != 0: self.totalWork = work - self.compressorWork

    def setMassFlow(self, m, BPR = 0):
        self.mTotal = m
        
        if BPR != 0:
            self.mCold = m*BPR/(BPR+1)
            self.mHot = m/(BPR+1)

    #TODO: add thrust 


class Inlet:
    
    #Constructor
    def __init__(self, gamma, cp, engine: Engine):
        #set variables from params
        self.gamma = gamma
        self.cp = cp
        self.stationStart = engine.station

        #get speed (engine) and efficiency (inlet)
        #assuming inlet means airplane
        engine.Ca = float(input("\nEnter aircraft speed (m/s): "))
        eta = float(input("Enter Intake Isentropic Efficiency: "))

        #set p and t in
        self.pIn = engine.pa
        self.tIn = engine.Ta

        #calculate T and P
        pRatio = self.inletPratio(engine.Ca, eta)
        self.pOut = pRatio*self.pIn
        self.tOut = self.inletTstag(engine.Ca)

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

    def inletPratio(self, Ca, etaI):
        T = self.tIn
        cp = self.cp
        gamma = self.gamma

        exp = gamma/(gamma-1)
        pRatio = (1+(etaI*Ca**2/(2*cp*T*1000)))**exp
        
        return pRatio
    

class Fan:

    #Constructor
    def __init__(self, pIn, tIn, gamma, cp, engine: Engine):
        #setting attributes
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = engine.station

        self.etaPoly = float(input("Enter the Fan Polytropic Efficiency: "))
        self.FPR = float(input("Enter Fan Pressure Ratio: "))
        engine.BPR = float(input("Enter Engine BPR: "))

        self.pOut = self.fanPstag()
        self.tOut = self.fanTstag()

        self.nozzle = FanNozzle(self.pOut, self.tOut, gamma, cp, engine)

    #printing
    def __str__(self) -> str:
        print(f"({self.stationStart}): P={self.pIn} T={self.tIn}")
        print(r"  //   ____     ____    \\")
        print(r" ||   /    \   /    \    ||")
        print(r" ||  (-------X-------)   ||")
        print(r" ||   \____/   \____/    ||")
        print("                           ")
        print(r"//\\                    //\\")
        print(self.nozzle)
        return("")
    
    #math
    def fanPstag(self): 
        pOut = self.FPR*self.pIn
        return pOut

    def fanTstag(self):
        gamma = self.gamma

        exp = (gamma-1)/(self.etaPoly*gamma)
        tRatio = self.FPR**exp
        tOut = tRatio*self.tIn

        return tOut
    
    #get
    def getTout(self): return self.tOut
    def getPout(self): return self.pOut
    def getStation(self): return self.stationStart
        

class Compressor:

    # Constructor
    def __init__(self, pIn, tIn, gamma, cp, engine: Engine):
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = engine.station

        efficiencyChoice = 0

        self.pRatio = float(input("\nEnter Compressor pressure ratio: "))

        print("\nWhat Efficiency is available?")
        while efficiencyChoice == 0:

            efficiencyChoice = int(input("(1) Isentropic\n(2) Polytropic\n"))

            if efficiencyChoice != 1 and efficiencyChoice != 2:
                efficiencyChoice = 0
                print("Select from the list given")
            else:
                self.eta = float(input("\nEnter Compressor Efficiency: "))
                
                #isentropic
                if efficiencyChoice == 1: self.tOut, self.delT = self.compTstag()
                #polytropic
                else: self.tOut, self.delT = self.compPolyTstag()

        work = self.compWork(engine.etaM)

        self.pOut = self.pRatio*pIn

        #adding work to engine
        engine.setCompressorWork(work)

        print(f'\nP out = {self.pOut}')
        print(f'Temp Out = {self.tOut}')
        print(f'Change in Temp = {self.delT}')
        print(f'Compressor Work = {work}')
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
    def getTout(self): return self.tOut
    def getPout(self): return self.pOut
    def getStation(self): return self.stationStart

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
        eta = self.eta
        pRatio = self.pRatio
        gamma = self.gamma
        
        exp = (gamma-1)/gamma
        delT = (T/eta)*((pRatio**exp)-1)
        tOut = delT+T

        return tOut, delT

    def compPolyTstag(self):
        #assuming only used for polytropic
        T = self.tIn
        etaC = self.eta
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
        self.stationStart = engine.station

        #assuming pDrop ONLY --> not pDrop/pIn
        self.pDrop = float(input("\nEnter pressure drop across combustor: "))
        self.pOut = pIn - self.pDrop
        
        #placeholder for turbine inlet temp
        self.tOut = tIn

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
        #TODO: need to account for multiple turb/comp --> will alwasy take total comp work
        #initialize params
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = engine.station

        #for checking
        efficiencyType = "Isentropic"
        pRatioCheck = False

        #Check for pressure ratio
        print("\nIs the turbine pressure ratio known?")
        pChoice = int(input("(1) yes\n(2) no\n"))

        #TODO: turbine pressure can be unknown but still use polytropic --> HW 2b

        #get pressure ratio
        if pChoice == 1: 
            self.pRatio = float(input("\nEnter Turbine pressure ratio: "))
            pRatioCheck = True
            efficiencyChoice = 0

            print("\nWhat Efficiency is available?")
            while efficiencyChoice == 0:

                efficiencyChoice = int(input("(1) Isentropic\n(2) Polytropic\n"))

                if efficiencyChoice != 1 and efficiencyChoice != 2:
                    efficiencyChoice = 0
                    print("Select from the list given")
                else:
                    #polytropic
                    if efficiencyChoice == 2: 
                        efficiencyType = "Polytropic"
        
        #get work to find pressure ratio
        else:
            #assuming always isentropic
            cWork = engine.compressorWork

            self.delT = cWork/cp
            self.tOut = tIn-self.delT

        self.eta = float(input(f"\nEnter Turbine {efficiencyType} Efficiency: "))

        #polytropic
        if pRatioCheck:
            if efficiencyChoice == 2: self.tOut, self.delT = self.turbPolyTstag()
            else: self.tOut, self.delT = self.turbTstag()

        #isentropic
        else: self.pRatio = self.turbPratio()

        # get pout and work
        self.pOut = pIn/self.pRatio
        work = self.turbWork(engine.etaM)

        #adding work to engine
        engine.setTurbineWork(work)

        print(f'\nP out = {self.pOut}')
        print(f'delta T = {self.delT}')
        print(f'Temp Out = {self.tOut}')
        print(f'Turbine Work = {work}\n')

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
    def getTout(self): return self.tOut
    def getPout(self): return self.pOut
    def getStation(self): return self.stationStart

    #Math
    def turbTstag(self):
        exp = (self.gamma-1)/self.gamma
        delT = self.eta*self.tIn*(1-((1/self.pRatio)**exp))
        tOut = self.tIn-delT

        return tOut, delT

    def turbPolyTstag(self):
        T = self.tIn
        etaT = self.eta
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
        etaT = self.eta
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
        self.stationStart = engine.station

        if engine.etaJ == 0: engine.etaJ = float(input("\nEnter Nozzle Isentropic Efficiency: "))

        self.chokeTest(engine)

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
    
    # Get
    def getTout(self): return self.tOut
    def getPout(self): return self.pOut
    def getStation(self): return self.stationStart

    #Math
    def critPRatio(self, etaJ):
        gamma = self.gamma

        den = (1-((1/etaJ)*(gamma-1)/(gamma+1)))**(gamma/(gamma-1))
        pRatio = 1/den
        
        return pRatio
    
    def chokedPratio(self):
        #Pout/Pin
        return (1/self.critPR)

    def chokedTratio(self):
        #TStaticOut/Tin
        return (2/(self.gamma+1))
    
    def nozzleRho(self, R):
        return (self.pOut/(R*self.tOut))
    
    def nozzleV(self, R):
        #assuming R needs to be /1000
        return np.sqrt(self.gamma*R*self.tOut*1000)
    
    def chokeTest(self, engine: Engine):
        self.critPR = self.critPRatio(engine.etaJ)
        self.pRatio = self.pIn/engine.pa

        if engine.R == 0: engine.R = float(input("Enter R: "))

        #choked
        if self.pRatio > self.critPR:
            print("\nNozzle is Choked")

            self.pOut = self.pIn*self.chokedPratio()
            self.tOut = self.tIn*self.chokedTratio()
            self.rhoOut = self.nozzleRho(engine.R)
            self.Cj = self.nozzleV(engine.R)
        #not choked
        else: 
            print("\nNozzle is not Choked")

            gamma = self.gamma

            self.pOut = engine.pa
            self.tOut = self.tIn

            delTstagTstatic = engine.etaJ*self.tOut*(1-(engine.pa/self.pIn)**((gamma-1)/gamma))
            tStaticOut = self.tOut-delTstagTstatic
            tStagStaticRatio = self.tOut/tStaticOut

            pStagStaticRatio = tStagStaticRatio**(gamma/(gamma-1))
            
            Mach = np.sqrt((pStagStaticRatio**((gamma-1)/gamma)-1)**(2/(gamma-1)))
            a = np.sqrt(gamma*engine.R*tStaticOut)
            c = Mach*a

    
class FanNozzle(Nozzle):
    
    #constructor
    def __init__(self, pIn, tIn, gamma, cp, engine: Engine):
        #initialize params
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = 0

        if engine.etaJ == 0: engine.etaJ = float(input("\nEnter Nozzle Isentropic Efficiency: "))

        self.chokeTest(engine)

        print(f'\nP out = {self.pOut}')
        print(f'Temp Out = {self.tOut}')
        print(f'Rho Out = {self.rhoOut}')

    def __str__(self) -> str:
        print(f"({self.stationStart}): P={self.pOut} T={self.tOut}")
        return("")
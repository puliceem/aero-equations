import AeroEqs as ae
import numpy as np

class Engine: 

    #variables
    components = []
    compressorWork = 0
    turbineWork = 0
    Ca = 0
    mTotal = 0
    etaJ = 0
    R = 0
    BPR = 0
    power = 0
    engineC = []
    engineP = []
    mDot = []
    engineA = []
    engineRho = []
    thrust = 0
    f = 0
    sfc = 0

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
        for component in self.components:
            if isinstance(component, Fan): component.nozzle.stationStart = self.station+1

    #post assembly
    #TODO: consider adding mass flow and power to the try methods
    def checkMassFlow(self):
        #mass flow
        massFlowChoice = 0
        print("Is mass flow known? ")
        while massFlowChoice == 0:

            massFlowChoice = int(input("(1) Yes\n(2) No\n"))

            #reprompt
            if massFlowChoice != 1 and massFlowChoice != 2:
                massFlowChoice = 0
                print("Select from the list given")
            #Yes - mass flow known - set power too
            elif massFlowChoice == 1: 
                self.setMassFlow(float(input("Enter Total Mass Flow: ")))
                self.power = self.totalWork*self.mTotal
            #No - check if power known
            else:
                self.checkPower()

    def checkPower(self):
        #power
        powerChoice = 0
        print("Is power known? ")
        while powerChoice == 0:

            powerChoice = int(input("(1) Yes\n(2) No\n"))

            #reprompt
            if powerChoice != 1 and powerChoice != 2:
                powerChoice = 0
                print("Select from the list given")
            #Yes - power known - set massflow
            elif powerChoice == 1: 
                self.power = float(input("Enter Power: "))
                massFlow = self.power/self.totalWork
                self.setMassFlow(massFlow)
            #No - power and massflow not known
            else:
                pass

    #trying to calculate
    #TODO: use try except (maybe somethign else) to find prompts
    def Force(self):
        index = 0
        for item in self.engineC:
            F = ae.thrust(self.mDot[index], item, self.Ca, self.engineA[index], self.engineP[index], self.pa)
            self.setThrust(F)

            index += 1

    def Area(self):
        index = 0
        for item in self.engineC:
            self.engineA.append(ae.nozzleArea(self.mDot[index], self.engineRho[index], item))

            index += 1

    def fuelFlow(self):
        for component in self.components: 
            if isinstance(component, Turbine):
                turb = component
            elif isinstance(component, Combustor):
                comb = component

        self.f = ae.fuelFlow(turb.tIn, comb.tIn, comb.eta, self.cpa, self.cpg)

    def SFC(self): self.sfc = ae.SFC(self.f, self.totalWork)

    #set
    def setCompressorWork(self, work):
        if self.compressorWork == 0: self.compressorWork = work
        else: self.compressorWork += work

    def setTurbineWork(self, work):
        if self.turbineWork == 0: self.turbineWork = work
        else: self.turbineWork += work

        if self.compressorWork != 0: self.totalWork = self.turbineWork - self.compressorWork

    def setMassFlow(self, m):
        self.mTotal = m
        BPR = self.BPR

        if BPR != 0:
            mCold, mHot = ae.fanMassFlow(m, BPR)
            self.mDot.append(mCold)
            self.mDot.append(mHot)

    def setThrust(self, F):
        if self.thrust == 0: self.thrust = F
        else: self.thrust += F


class Inlet:
    
    #properties
    pRatio = 0
    pOut = 0
    tOut = 0

    #Constructor
    def __init__(self, gamma, cp, engine: Engine):
        #set variables from params
        self.gamma = gamma
        self.cp = cp
        self.stationStart = engine.station

        #get speed (engine) and efficiency (inlet)
        #assuming inlet means airplane
        engine.Ca = float(input("\nEnter aircraft speed (m/s): "))
        self.eta = float(input("Enter Intake Isentropic Efficiency: "))

        #set p and t in
        self.pIn = engine.pa
        self.tIn = engine.Ta

        #calculate T and P
        self.inletPout(engine.Ca)
        self.inletTstag

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
    def inletPout(self, Ca): self.pRatio, self.pOut = ae.inletPstagStatic(self.tIn, Ca, self.cp, self.gamma, self.eta, self.pIn)

    def inletTstag(self, Ca):
        #assuming Ca = m/s
        self.tOut = ae.inletTstag(self.tIn, Ca, self.cp)
    

class Fan:

    #properties
    pOut = 0
    tOut = 0
    tRatio = 0

    #Constructor
    def __init__(self, pIn, tIn, gamma, cp, engine: Engine):
        #setting attributes
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = engine.station

        self.eta = float(input("Enter the Fan Polytropic Efficiency: "))
        self.FPR = float(input("Enter Fan Pressure Ratio: "))
        engine.BPR = float(input("Enter Engine BPR: "))

        self.fanPstag()
        self.fanTstag()

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
    def fanPstag(self): self.pOut = ae.fanPstag(self.FPR, self.pIn)

    def fanTstag(self): self.tRatio, self.tOut = ae.fanTstag(self.tIn, self.gamma, self.eta, self.FPR)
    
    #get
    def getTout(self): return self.tOut
    def getPout(self): return self.pOut
    def getStation(self): return self.stationStart
        

class Compressor:

    #properties
    tOut = 0
    tRatio = 0
    work = 0

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
                if efficiencyChoice == 1: self.compTstag()
                #polytropic
                else: self.compPolyTstag()

        self.compWork(engine.etaM)

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
    def getTout(self): return self.tOut
    def getPout(self): return self.pOut
    def getStation(self): return self.stationStart

    #Math
    def compTstag(self): 
        #assuming only isentropic
        self.tOut, self.delT = ae.compTstag(self.tIn, self.eta, self.pRatio, self.gamma)

    def compPolyTstag(self):
        #assuming only polytropic
        self.tOut, self.delT = ae.compPolyTstag(self.tIn, self.eta, self.pRatio, self.gamma)

    def compWork(self, etaM): self.work = ae.compWork(self.delT, self.cp, etaM)


class Combustor:

    # Constructor
    def __init__(self, pIn, tIn, engine: Engine):
        self.pIn = pIn
        self.tIn = tIn
        self.stationStart = engine.station

        #assuming pDrop ONLY --> not pDrop/pIn
        self.pDrop = float(input("\nEnter pressure drop across combustor: "))
        self.pOut = pIn - self.pDrop

        #efficiency
        #TODO: determine if this is necessary here or only in engine.fuelFlow
        self.eta = float(input("\nEnter the Combustion Isentropic Efficiency: "))
        
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

    #properties
    tOut = 0
    delT = 0
    pRatio = 0
    work = 0

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

                #reprompt until correct choice given
                if efficiencyChoice != 1 and efficiencyChoice != 2:
                    efficiencyChoice = 0
                    print("Select from the list given")
                else:
                    #polytropic
                    if efficiencyChoice == 2: efficiencyType = "Polytropic"
        
        #get work to find pressure ratio
        else:
            #assuming always isentropic
            cWork = engine.compressorWork

            self.delT = cWork/cp
            self.tOut = tIn-self.delT

        self.eta = float(input(f"\nEnter Turbine {efficiencyType} Efficiency: "))

        #polytropic
        if pRatioCheck:
            if efficiencyChoice == 2: self.turbPolyTstag()
            else: self.turbTstag()

        #isentropic
        else: self.turbPratio()

        # get pout and work
        self.pOut = pIn/self.pRatio
        self.turbWork(engine.etaM)

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
    def getTout(self): return self.tOut
    def getPout(self): return self.pOut
    def getStation(self): return self.stationStart

    #Math
    def turbTstag(self):
        #assuming only isentropic
        self.tOut, self.delT = ae.turbTstag(self.tIn, self.eta, self.gamma, self.pRatio)

    def turbPolyTstag(self):
        #assuming only polytropic
        self.tOut, self.delT = ae.turbPolyTstag(self.tIn, self.eta, self.gamma, self.pRatio)

    def turbPratio(self):
        #assuming only for isentropic
        self.pRatio = ae.turbPratio(self.tIn, self.delT, self.eta, self.gamma)

    def turbWork(self, etaM): self.work = ae.turbWork(self.delT, self.cp, etaM)
    
class Nozzle:

    #properties
    pOut = 0
    tOut = 0
    rhoOut = 0
    pRatio = 0
    tRatio = 0
    critPR = 0
    Cj = 0

    # constructor
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
    def critPRatio(self, etaJ): self.critPR = ae.critPRatio(self.gamma, etaJ)
    
    def chokedPratio(self):
        #Pout/Pin
        self.pRatio, self.pOut = ae.chokedPratio(self.critPR, self.pIn)

    def chokedTratio(self):
        #TStaticOut/Tin
        self.tRatio, self.tOut = ae.chokedTratio(self.gamma, self.tIn)
    
    def nozzleRho(self, R): self.rhoOut = ae.nozzleRho(self.pOut, self.tOut, R)
    
    def nozzleV(self, R):
        #assuming R needs to be /1000
        self.Cj = ae.nozzleV(self.tOut, self.gamma, R)
    
    def chokeTest(self, engine: Engine):
        self.critPRatio(engine.etaJ)
        self.pRatio = self.pIn/engine.pa

        if engine.R == 0: engine.R = float(input("Enter R: "))

        #choked
        if self.pRatio > self.critPR:
            print("\nNozzle is Choked")

            self.chokedPratio()
            self.chokedTratio()
            self.nozzleRho(engine.R)
            self.nozzleV(engine.R)
        #not choked
        else: 
            #TODO: put equations in AeroEqs
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
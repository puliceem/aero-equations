import AeroEqs as ae

class Engine: 

    #variables
    components = []
    compressorWork = 0
    turbineWork = 0
    totalWork = 0
    Ca = 0
    mTotal = 0
    mRatio = 0
    etaJ = 0
    R = 0
    BPR = 0
    power = 0
    engineC = []
    engineP = []
    mDot = []
    engineA = []
    engineRho = []
    engineT = []
    thrustTotal = 0
    f = 0
    sfc = 0
    extraAnalysis = {}

    def __init__(self, gammaA=1.4, gammaG=1.333, cpa=1.005, cpg=1.148, pa=None, Ta=None, Ca=None, etaM=None, station=None):
        self.gammaA = gammaA
        self.gammaG = gammaG
        self.cpa = cpa
        self.cpg = cpg

        if pa == None: self.pa = float(input("Enter ambient pressure: "))
        else: self.pa = pa
        
        if Ta == None: self.Ta = float(input("Enter ambient temp: "))
        else: self.Ta = Ta
        
        if Ca == None: self.Ca = float(input("Enter aircraft speed (m/s): "))
        else: self.Ca = Ca
        
        if etaM == None: self.etaM = float(input("Enter Mechanical Efficiency: "))
        else: self.etaM = etaM
        
        if station == None: self.station = int(input("Enter the first station number: "))
        else: self.station = station

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
    def addComponent(self, selection, *args):
        gammaA = self.gammaA
        gammaG = self.gammaG
        cpa = self.cpa
        cpg = self.cpg
        pIn = self.pIn
        tIn = self.tIn

        if args: params = list(args)
        #TODO: may be a better way of doing this
        else: params = [None]*5

        # Inlet
        if selection == 1:
            component = Inlet(gammaA, cpa, self, eta=params[0])

        #Fan
        elif selection == 2:
            component = Fan(pIn, tIn, gammaA, cpa, self, eta=params[0], FPR=params[1], BPR=params[2], etaJ=params[3], R=params[4])

        # Compressor
        elif selection == 3:
            component = Compressor(pIn, tIn, gammaA, cpa, self, pRatio=params[0], eta=params[1], etaChoice=params[2])

        # Combuster
        # assuing no heat exchanger
        elif selection == 4:
            component = Combustor(pIn, tIn, self, pRatio=params[0], eta=params[1])

        #Turbine
        elif selection == 5:
            if not self.findComponent(Turbine) and not params: tIn = float(input("\nEnter the Turbine Inlet Temp: "))
            component = Turbine(pIn, tIn, gammaG, cpg, self, eta=params[0], pRatio=params[1], etaChoice=params[2])

        #Nozzle
        elif selection == 6: 
            component = Nozzle(pIn, tIn, gammaG, cpg, self, eta=params[0], R=params[1], typeChoice=params[2])

        self.components.append(component)
        self.pIn = component.pOut
        self.tIn = component.tOut

        self.station += 1

        #update fan nozzle exit station number
        fan = self.findComponent(Fan)
        if fan: fan.nozzle.stationStart = self.station+1

    #Methods
    def removeComponent(self):
        #remove latest entry from list
        self.components.pop()

        if self.components:
            #reset pIn tIn and station
            component = self.components[-1]
            self.pIn = component.pOut
            self.tIn = component.tOut
        else: 
            self.pIn = self.pa
            self.tIn = self.Ta

        self.station -= 1

    def returnVariable(self, index, *args):
        output = []
        
        if args: outputVariables = list(args)
        else: 
            count = 0
            variables = []
            outputVariables = []

            print("\nEnter (0) to quit or select from list")
            for variable in vars(self.components[index]):
                count+=1
                variables.append(variable)
                print(f"({count}) {variable}")

            while True:
                selection = int(input("Return: "))
                
                if selection == 0: break
                elif selection > count: print("\nPlease select a number form the list\n")
                else: outputVariables.append(variables[selection-1])
        
        for var in outputVariables: 
            output.append(getattr(self.components[index], var))
        
        return output

    def findComponent(self, type, first=True):
        #assuming each component has a max of 2
        foundComponent = False

        for component in self.components:
            if isinstance(component, type): 
                #assign index and iterate again if needed from 1 after
                foundComponent = component

                #break if first is to be found
                if first: break

        return foundComponent

    #post assembly
    def checkMassFlow(self):
        #mass flow
        if self.mTotal == 0:
            massFlowChoice = 0
            print("\nIs mass flow known? ")
            while massFlowChoice == 0:

                massFlowChoice = int(input("(1) Yes\n(2) No\n"))

                #reprompt
                if massFlowChoice != 1 and massFlowChoice != 2:
                    massFlowChoice = 0
                    print("Select from the list given")
                #Yes - mass flow known - set power too
                elif massFlowChoice == 1: 
                    self.setMassFlow(float(input("\nEnter Total Mass Flow: ")))

                    if self.totalWork: self.power = self.totalWork*self.mTotal

                    return True
                #No - check if power known
                else:
                    return self.checkPower()
        else: return True

    #TODO: does not check for mass flow if used first
    def checkPower(self):
        #power
        powerChoice = 0
        print("\nIs power known? ")
        while powerChoice == 0:

            powerChoice = int(input("(1) Yes\n(2) No\n"))

            #reprompt
            if powerChoice != 1 and powerChoice != 2:
                powerChoice = 0
                print("Select from the list given")
            #Yes - power known - set massflow
            elif powerChoice == 1: 
                self.power = float(input("\nEnter Power: "))
                massFlow = self.power/self.totalWork
                self.setMassFlow(massFlow)

                return True
            #No - power and massflow not known
            else:
                return False

    #trying to calculate
    #TODO: add any missing params
    def analysis(self):
        #check mass
         if self.checkMassFlow():

            #area
            try:
                self.Area()

                if len(self.engineA) == 2:
                    print(f"Fan Nozzle Area: {self.engineA[0]}")
                
                print(f"Core Nozzle Area: {self.engineA[1]}")

            except: print("Nozzle area cannot be found")

            #thrust
            try: 
                self.Force()

                if len(self.engineT) == 2:
                    print(f"Cold Stream Thrust: {self.engineT[0]}")
                    print(f"Hot Stream Thrust: {self.engineT[1]}")

                print(f"Total Engine Thrust: {self.thrustTotal}")

            except: print("Engine thrust cannot be found")
            
            #fuel flow
            try:
                self.fuelFlow()
                print(f"Fuel Flow: {self.f}")
            
            except: print("Fuel Flow cannot be found")

            #sfc
            try: 
                self.SFC()
                print(f"SFC: {self.sfc}")

            except: print("SFC cannot be found")


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

    #TODO: fuel flow and SFC may be calculated wrong --> wrong answers on hw
    def fuelFlow(self):
        for component in self.components: 
            if isinstance(component, Turbine):
                turb = component
            elif isinstance(component, Combustor):
                comb = component

        self.f = ae.fuelFlow(turb.tIn, comb.tIn, comb.eta, self.cpa, self.cpg)

    #TODO: figure out when which one needs to be used
    def SFC(self): self.sfc = ae.SFC(self.f, self.totalWork)
    def TSFC(self): self.tsfc = ae.TSFC(self.Ca, self.etaT, self.etaP)

    def overallEfficiency(self):
        mDotC = self.mDot[0]
        mDotH = self.mDot[1]
        CjC = self.engineC[0]
        CjH = self.engineC[1]

        self.etaT = ae.etaT(mDotC, mDotH, self.mDot, CjC, CjH, self.Ca, self.f)
        self.etaP = ae.etaP(mDotC, mDotH, self.mDot, CjC, CjH, self.Ca)
        self.etaO = self.etaT*self.etaP

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
        self.engineT.append(F)

        if self.thrustTotal == 0: self.thrustTotal = F
        else: self.thrustTotal += F


class Inlet:
    
    #properties
    pRatio = 0
    pOut = 0
    tOut = 0

    #Constructor
    def __init__(self, gamma, cp, engine: Engine, eta=None):
        #set variables from params
        self.gamma = gamma
        self.cp = cp
        self.stationStart = engine.station

        #get speed (engine) and efficiency (inlet)
        #TODO: check if inletPout needs changed
        if eta == None: self.eta = float(input("Enter Intake Isentropic Efficiency: "))
        else: self.eta = eta

        #set p and t in
        self.pIn = engine.pa
        self.tIn = engine.Ta

        #calculate T and P
        self.inletPout(engine.Ca)
        self.inletTstag(engine.Ca)

        print(f'\nP out = {self.pOut}')
        print(f'Temp Out = {self.tOut}')

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
    def __init__(self, pIn, tIn, gamma, cp, engine: Engine, eta=None, FPR=None, BPR=None, etaJ=None, R=None):
        #setting attributes
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = engine.station

        if eta == None: self.eta = float(input("\nEnter the Fan Polytropic Efficiency: "))
        else: self.eta = eta
        
        if FPR == None: self.FPR = float(input("Enter Fan Pressure Ratio: "))
        else: self.FPR = FPR
        
        if BPR == None: engine.BPR = float(input("Enter Engine BPR: "))
        else: self.BPR = BPR

        self.fanPstag()
        self.fanTstag()

        self.nozzle = FanNozzle(self.pOut, self.tOut, gamma, cp, engine, etaJ, R=R)

        print(f'\nP out = {self.pOut}')
        print(f'Temp Out = {self.tOut}')

    #printing
    def __str__(self) -> str:
        print(f"({self.stationStart}): P={self.pIn} T={self.tIn}")
        print(r"  //   ____             \\")
        print(r" ||   /    \             ||")
        print(r" ||  (-------X-------)   ||")
        print(r" ||            \____/    ||")
        print("                           ")
        print(r"//\\                    //\\")
        #TODO: looks horrible
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
    power = 0

    # Constructor
    def __init__(self, pIn, tIn, gamma, cp, engine: Engine, pRatio=None, eta=None, etaChoice=None):
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = engine.station

        efficiencyChoice = 0

        if pRatio == None: self.pRatio = float(input("\nEnter Compressor pressure ratio: "))
        else: self.pRatio = pRatio

        if eta == None: 

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
        
        else: 
            self.eta == eta
            
            if etaChoice == 1: self.compTstag()
            else: self.compPolyTstag()

        self.compWork(engine.etaM)

        self.pOut = self.pRatio*pIn

        #adding work to engine
        engine.setCompressorWork(self.work)

        print(f'\nP out = {self.pOut}')
        print(f'Temp Out = {self.tOut}')
        print(f'Compressor Work = {self.work}')

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
    def __init__(self, pIn, tIn, engine: Engine, pDrop=None, eta=None):
        self.pIn = pIn
        self.tIn = tIn
        self.stationStart = engine.station

        #assuming pDrop ONLY --> not pDrop/pIn
        if pDrop == None: self.pDrop = float(input("\nEnter pressure drop across combustor: "))
        else: self.pDrop = pDrop
        
        self.pOut = pIn - self.pDrop

        #efficiency
        if eta == None: self.eta = float(input("Enter the Combustion Isentropic Efficiency: "))
        else: self.eta = eta
        
        #placeholder for turbine inlet temp
        self.tOut = tIn

        print(f'\nP out = {self.pOut}')

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
    def __init__(self, pIn, tIn, gamma, cp, engine: Engine, eta=None, pRatio=None, etaChoice=None):
        #initialize params
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = engine.station

        if eta == None:

            #for checking
            efficiencyType = "Isentropic"
            efficiencyChoice = 0

            print("\nWhat Efficiency is available?")
            while efficiencyChoice == 0:

                efficiencyChoice = int(input("(1) Isentropic\n(2) Polytropic\n"))

                #reprompt until correct choice given
                if efficiencyChoice != 1 and efficiencyChoice != 2:
                    efficiencyChoice = 0
                    print("Select from the list given")
                #polytropic
                elif efficiencyChoice == 2: efficiencyType = "Polytropic"

            #get turbine efficiency
            self.eta = float(input(f"\nEnter Turbine {efficiencyType} Efficiency: "))

            #Check for pressure ratio
            print("\nIs the turbine pressure ratio known?")
            pChoice = int(input("(1) yes\n(2) no\n"))

            #TODO: should this default to powerBalance if possible?
            #get pressure ratio
            if pChoice == 1: 
                self.pRatio = float(input("\nEnter Turbine pressure ratio: "))
                #poly
                if efficiencyChoice == 2: self.turbPolyTstag()
                #isen
                else: self.turbTstag()
            
            #do power balance to get t
            else: self.powerBalance(engine)

        else: 
            self.eta = eta

            #no pRatio supplied
            if pRatio == None: self.powerBalance(engine)
            #pRatio supplied
            else: 
                self.pRatio = pRatio

                if etaChoice == 2: self.turbPolyTstag()
                else: self.turbTstag()
        
        # get pout and work
        self.pOut = pIn/self.pRatio
        self.turbWork(engine.etaM)

        #adding work to engine
        engine.setTurbineWork(self.work)

        print(f'\nP out = {self.pOut}')
        print(f'Temp Out = {self.tOut}')
        print(f'Turbine Work = {self.work}')

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
    
    #methods
    def powerBalance(self, engine: Engine):
        #initilizing
        component = False

        #assumign max 2 compressors 1 fan
        if engine.findComponent(Turbine):
            #other turbines exist
            
            if engine.findComponent(Fan):
                #fan exists
                fan: Fan = engine.findComponent(Fan)
                
                #check/set mdot
                if engine.checkMassFlow():
                    #calculate T values
                    self.fanPowerBalance(engine.mTotal, engine.mDot[1], fan.cp, self.cp, engine.etaM, fan.tOut, fan.tIn, self.tIn)

                #no mdot given --> mratio used instead
                else:
                    engine.mRatio = engine.BPR+1

                    self.delT = engine.mRatio*fan.cp*(fan.tOut-fan.tIn)/(engine.etaM*self.cp)
                    self.tOut = self.tIn-self.delT
                
                print(self.tOut)

            #find first compressor
            else: component = engine.findComponent(Compressor)
        
        #first turbine
        else:
            #find last compressor
            component = engine.findComponent(Compressor, first=False)

        #calculate T values (compressor)
        if component: self.compPowerBalance(component.cp, engine.etaM, component.tOut, component.tIn, self.tIn)

        self.powerBalancePratio()

    def compressorWork(self, engine: Engine):
        comp = engine.findComponent(Compressor, first=False)
        if comp: return comp.work
        else: 
            #TODO: if pass happens what happens next?
            pass

    # Get
    def getTout(self): return self.tOut
    def getPout(self): return self.pOut
    def getStation(self): return self.stationStart

    #Math
    def turbTstag(self):
        #assuming only isentropic
        self.tOut, self.delT = ae.turbTstag(self.tIn, self.eta, self.pRatio, self.gamma)

    def turbPolyTstag(self):
        #assuming only polytropic
        self.tOut, self.delT = ae.turbPolyTstag(self.tIn, self.eta, self.gamma, self.pRatio)

    def turbPratio(self):
        #assuming only for isentropic
        self.pRatio = ae.turbPratio(self.tIn, self.delT, self.eta, self.gamma)

    def turbWork(self, etaM): self.work = ae.turbWork(self.delT, self.cp, etaM)

    def fanPowerBalance(self, mTotal, mh, cpa, etaM, fanTout, fanTin, turbTin):
        self.tOut, self.delT = ae.fanPowerBalance(mTotal, mh, cpa, self.cp, etaM, fanTout, fanTin, turbTin=turbTin)

    def compPowerBalance(self, cpa, etaM, compTout, compTin, turbTin):
        self.tOut, self.delT = ae.compPowerBalance(cpa, self.cp, etaM, compTout, compTin, turbTin=turbTin)

    def powerBalancePratio(self):
        self.pRatio = ae.powerBalancePratio(self.tOut, self.tIn, self.eta, self.gamma)


class Nozzle:

    #properties
    pOut = 0
    tOut = 0
    tStaticOut = 0
    pStaticOut = 0
    rhoOut = 0
    pRatio = 0
    pOutRatio = 0
    tRatio = 0
    tOutRatio = 0
    critPR = 0

    # constructor
    def __init__(self, pIn, tIn, gamma, cp, engine: Engine, eta=None, R=None, typeChoice=None):
        #initialize params
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = engine.station

        if engine.etaJ == 0: 
            if eta == None: engine.etaJ = float(input("\nEnter Nozzle Isentropic Efficiency: "))
            else: engine.etaJ = eta

        if typeChoice == None: type = int(input("Is the Nozzle CD?\n(1) Yes\n(2)\n"))
        else: type = typeChoice

        if type == 1: pass
        #TODO: only accounts for converging
        else: self.chokeTest(engine, R=R)

        print(f'\nP(static) out = {self.pStaticOut}')
        print(f'Temp(static) Out = {self.tStaticOut}')
        #assuming last or only C calculated
        print(f'Cj = {engine.engineC[-1]}')

    #printing
    def __str__(self) -> str:
        print(f"({self.stationStart}): P={self.pIn} T={self.tIn}")
        print("    |                     |")
        print(r"     \                   /")
        print(r"      \                 /")
        print(r"       \               /")
        print("        @             @")
        return("")
    
    #methods
    def chokeTest(self, engine: Engine, R=None):
        self.critPRatio(engine.etaJ)
        self.pRatio = self.pIn/engine.pa

        if engine.R == 0: 
            if R == None: engine.R = float(input("Enter R: "))
            else: engine.R = R

        #choked
        if self.pRatio > self.critPR:
            print("\nNozzle is Choked")

            self.chokedPratio()
            self.chokedTratio()
            self.nozzleRho(engine.R)
            engine.engineC.append(self.nozzleV(engine.R))
        #not choked
        else: 
            print("\nNozzle is not Choked")

            self.pOut = engine.pa
            self.tOut = self.tIn

            self.notChokedToutRatio(engine.etaJ, engine.pa)
            self.notChokedPoutRatio()
            engine.engineC.append(self.notChokedCj(engine.R))

    # Get
    def getTout(self): return self.tOut
    def getPout(self): return self.pOut
    def getStation(self): return self.stationStart

    #Math
    def critPRatio(self, etaJ): self.critPR = ae.critPRatio(self.gamma, etaJ)
    
    def chokedPratio(self):
        #Pout/Pin
        self.pRatio, self.pStaticOut = ae.chokedPratio(self.critPR, self.pIn)

    def chokedTratio(self):
        #TStaticOut/Tin
        self.tRatio, self.tStaticOut = ae.chokedTratio(self.gamma, self.tIn)
    
    def nozzleRho(self, R): self.rhoOut = ae.nozzleRho(self.pStaticOut, self.tStaticOut, R)
    
    def nozzleV(self, R, Mach=1):
        #assuming R needs to be /1000
        return ae.nozzleV(self.tStaticOut, self.gamma, R, Mach)

    def notChokedToutRatio(self, etaJ, pa):
        #tStaticOut/Tin, T
        self.tOutRatio, self.tStaticOut = ae.notChokedToutRatio(etaJ, self.tOut, pa, self.pIn, self.gamma)

    def notChokedPoutRatio(self): self.pOutRatio = ae.noChokedPoutRatio(self.tOutRatio, self.gamma)

    def notChokedCj(self, R):
        Mach = ae.nozzleMach(self.pOutRatio, self.gamma)
        return self.nozzleV(R, Mach)
    
class FanNozzle(Nozzle):
    
    #constructor
    def __init__(self, pIn, tIn, gamma, cp, engine: Engine, eta=None, R=None):
        #initialize params
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = 0

        if engine.etaJ == 0: 
            if eta == None: engine.etaJ = float(input("\nEnter Nozzle Isentropic Efficiency: "))
            else: engine.etaJ = eta

        self.chokeTest(engine, R=R)

        #TODO: p(static) may not work for not choked
        print(f'\nP(static) out = {self.pStaticOut}')
        print(f'Temp(static) Out = {self.tStaticOut}')
        #assuming alwasy first to have C calculated
        print(f'Cj = {engine.engineC[0]}')

    def __str__(self) -> str:
        print(f"({self.stationStart}): P={self.pOut} T={self.tOut}")
        return("")
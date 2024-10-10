import AeroEqs as ae

class Engine: 

    #variables
    components = []
    compressorWork = 0
    turbineWork = 0
    Ca = 0

    def __init__(self, etaM):
        self.etaM = etaM

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

    #get
    def getEtaM(self): return self.etaM
    def getWorkC(self): return self.compressorWork


class Inlet:
    def __init__(self, pIn, tIn, stationStart, gamma, cp, engine: Engine):
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = stationStart

        self.Ca = float(input("\nEnter aircraft speed: "))
        self.eta = float(input("Enter Intake Isentropic Efficiency: "))

        
        

class Compressor:

    # Constructor
    def __init__(self, pIn, tIn, stationStart, gamma, cp, engine: Engine):
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = stationStart

        efficiencyChoice = 0

        self.pRatio = float(input("\nEnter Compressor pressure ratio: "))

        print("\nWhat Efficiency is available?")
        while efficiencyChoice == 0:

            #TODO: still has input taken even if wrong choice is given
            efficiencyChoice = int(input("(1) Isentropic\n(2) Polytropic\n"))
            eta = float(input("\nEnter Compressor Efficiency: "))
            
            if efficiencyChoice == 1: 
                self.etaIsen = eta
                self.tOut, self.delT = ae.compTstag(tIn, eta, self.pRatio, gamma)
            elif efficiencyChoice == 2: 
                self.etaPoly = eta
                self.tOut, self.delT = ae.compPolyTstag(tIn, eta, self.pRatio, gamma)
            else: 
                efficiencyChoice = 0
                print("Select from the list given")

        self.work = ae.compWork(cp, self.delT, engine.getEtaM())

        self.pOut = self.pRatio*pIn

        #adding to engine
        engine.addComponent(self)
        engine.setCompressorWork(self.work)

        print(f'\nP out = {self.pOut}')
        print(f'Temp Out = {self.tOut}')
        print(f'Change in Temp = {self.delT}')
        print(f'Compressor Work = {self.work}')
        print("\n")

    def __str__(self) -> str:
        print(f"({self.stationStart}): P={self.pIn} T={self.tIn}")
        print(" ____________________")
        print("| ------------------ |")
        print("| ------------------ |")
        print(r" \ ---------------- /")
        print(r"  \ -------------- /")
        print(r"   \ ____________ /")
        print("    |            |")
        return("")
    
    # Get
    def getWork(self): return self.work
    def getTout(self): return self.tOut
    def getPout(self): return self.pOut
    def getStation(self): return self.stationStart


class Combustor:

    # Constructor
    def __init__(self, pIn, tIn, stationStart, engine: Engine):
        self.pIn = pIn
        self.tIn = tIn
        self.stationStart = stationStart

        #assuming pDrop ONLY --> not pDrop/pIn
        self.pDrop = float(input("\nEnter pressure drop across combustor: "))
        self.pOut = pIn - self.pDrop
        
        #TODO: assuming turbine always next
        self.tOut = float(input("Enter the Turbine Inlet Temp: "))

        #adding to engine
        engine.addComponent(self)

        print(f'\nP out = {self.pOut}\n')

    def __str__(self) -> str:
        print(f"({self.stationStart}): P={self.pIn} T={self.tIn}")
        print("    |_____||_____|")
        print(r"   //\    ||    /\\")
        print(r"  //  \   ||   /  \\")
        print(r"  \\  /   ||   \  //")
        print(r"   \\/____||____\//")
        print("    |     ||     |")
        return("")

    # Get
    def getTout(self): return self.tOut
    def getPout(self): return self.pOut
    def getStation(self): return self.stationStart


class Turbine:

    # Constructor
    def __init__(self, pIn, tIn, stationStart, gamma, cp, engine: Engine):
        #initialize params
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = stationStart

        #for checking
        efficiencyType = "Polytropic"
        pRatioCheck = True

        #Check for pressure ratio
        print("\nIs the turbine pressure ratio known?")
        choice = int(input("(1) yes\n(2) no\n"))

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
            self.tOut, self.delT = ae.turbPolyTstag(tIn, eta, self.pRatio, gamma)

        #isentropic
        else:
            self.etaIsen = eta
            self.pRatio = ae.turbPratio(tIn, self.delT, eta, gamma)

        # get pout and work
        self.pOut = pIn/self.pRatio
        self.work = ae.turbWork(cp, self.delT, engine.getEtaM())

        #adding to engine
        engine.addComponent(self)
        engine.setTurbineWork(self.work)

        print(f'\nP out = {self.pOut}')
        print(f'delta T = {self.delT}')
        print(f'Temp Out = {self.tOut}')
        print(f'Turbine Work = {self.work}\n')

    def __str__(self) -> str:
        print(f"({self.stationStart}): P={self.pIn} T={self.tIn}")
        print("    |     ||     |")
        print(R"   /      ||      \ ")
        print(R"  / -------------- \ ")
        print(R" /        ||        \ ")
        print(R"/ ------------------ \ ")
        print("|         ||          |")
        return("")
    
    # Get
    def getWork(self): return self.work
    def getTout(self): return self.tOut
    def getPout(self): return self.pOut
    def getStation(self): return self.stationStart
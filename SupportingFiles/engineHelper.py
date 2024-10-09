import AeroEqs as ae

class Engine: 

    #variables
    components = []
    compressorWork = 0
    turbineWork = 0
    totalWork = 0

    #printing
    def __str__(self) -> str:
        for i in self.components: print(i)
        return("")

    #methods
    def addComponent(self, component):
        self.components.append(component)

    #set
    def setCompressorWork(self, work):
        self.compressorWork = work

    def setTurbineWork(self, work):
        self.turbineWork = work
        if self.compressorWork != 0: self.totalWork = work - self.compressorWork

def compressor(pIn, tIn, gamma, cp):
    etaM = 0
    work = 0

    pRatio = float(input("\nEnter Compressor pressure ratio: "))

    print("What Efficiency is available?")
    efficiencyChoice = int(input("(1) Isentropic\n(2) Polytropic"))
    eta = float(input("Enter Compressor Efficiency: "))

    #TODO: consider assuming 1
    print("Is the Mechanical Efficiency known?")
    choice = int(input("(1) yes\n(2) no\n"))
    if choice == 1: etaM = float(input("\nEnter Mechanical Efficiency: "))

    #assuming cold air
    tOut, delT = ae.compTstag(tIn, eta, pRatio, gamma)

    pOut = pRatio*pIn

    print(f'\nP out = {pOut}')
    print(f'Tout - Tin = {delT}')
    print(f'Temp Out = {tOut}')
    if choice == 1: 
        work = ae.compWork(cp, delT, etaM)
        print(f'Compressor Work = {work}')
    print("\n")

    return pOut, tOut, work

def combustor(pIn):

    #assuming pDrop ONLY --> not pDrop/pIn
    pDrop = float(input("\nEnter pressure drop across combustor: "))

    pOut = pIn - pDrop

    print(f'\nP out = {pOut}\n')

    return pOut

def turbine(pIn, tIn, gamma, cp, cWork):
    pRatio = 0
    delT = 0
    tOut = 0

    eta = float(input("\nEnter Turbine Efficiency: "))
    print("\nIs the turbine pressure ratio known?")
    choice = int(input("(1) yes\n(2) no\n"))

    if choice == 1: 
        pRatio = float(input("\nEnter Turbine pressure ratio: "))
        tOut, delT = ae.turbTstag(tIn, eta, pRatio, gamma)
    
    else: 
        if cWork == 0: cWork = float(input("\nEnter the compressor work: "))

        delT = cWork/cp
        tOut = tIn-delT

        pRatio = ae.turbPratio(tIn, delT, eta, gamma)

    pOut = pRatio*pIn
    work = ae.turbWork(cp, delT, eta)

    print(f'\nP out = {pOut}')
    print(f'delta T = {delT}')
    print(f'Temp Out = {tOut}')
    print(f'Turbine Work = {work}\n')

    return pOut, tOut, work
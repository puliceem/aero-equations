class Compressor:

    # Variables
    pOut = 0
    tOut = 0
    work = 0

    # Constructor
    def __init__(self, pIn, tIn, stationStart, gamma=1.4, cp=1.005):
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = stationStart

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
    
    # Set
    def setP_Out(self, pOut):
        self.pOut = pOut

    def setT_Out(self, tOut):
        self.tOut - tOut

    def setWork(self, work):
        self.work = work

class Combustor:

    # Variables
    pOut = 0

    # Constructor
    def __init__(self, pIn, tIn, stationStart):
        self.pIn = pIn
        self.tIn = tIn
        self.stationStart = stationStart

    def __str__(self) -> str:
        print(f"({self.stationStart}): P={self.pIn} T={self.tIn}")
        print("    |_____||_____|")
        print(r"   //\    ||    /\\")
        print(r"  //  \   ||   /  \\")
        print(r"  \\  /   ||   \  //")
        print(r"   \\/____||____\//")
        print("    |     ||     |")
        return("")
    
    # Set
    def setP_Out(self, pOut):
        self.pOut = pOut

class Turbine:

    # Variables
    pOut = 0
    tOut = 0
    work = 0

    # Constructor
    def __init__(self, pIn, tIn, stationStart, gamma=1.4, cp=1.005):
        self.pIn = pIn
        self.tIn = tIn
        self.gamma = gamma
        self.cp = cp
        self.stationStart = stationStart

    def __str__(self) -> str:
        print(f"({self.stationStart}): P={self.pIn} T={self.tIn}")
        print("    |     ||     |")
        print(R"   /      ||      \ ")
        print(R"  / -------------- \ ")
        print(R" /        ||        \ ")
        print(R"/ ------------------ \ ")
        print("|         ||          |")
        return("")
    
    # Set
    def setP_Out(self, pOut):
        self.pOut = pOut

    def setT_Out(self, tOut):
        self.tOut - tOut

    def setWork(self, work):
        self.work = work
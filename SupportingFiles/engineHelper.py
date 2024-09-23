def compressor(pIn, tIn, gamma, cp):
    etaM = 0
    work = 0

    pRatio = float(input("\nEnter Compressor pressure ratio: "))
    eta = float(input("Enter Compressor Efficiency: "))

    print("Is the Mechanical Efficiency known?")
    choice = int(input("(1) yes\n(2) no\n"))
    if choice == 1: etaM = float(input("\nEnter Mechanical Efficiency: "))

    #assuming cold air
    exp = (gamma-1)/(gamma)
    delT = (tIn/eta)*((pRatio**exp)-1)
    tOut = delT+tIn

    pOut = pRatio*pIn

    print(f'\nP out = {pOut}')
    print(f'Tout - Tin = {delT}')
    print(f'Temp Out = {tOut}')
    if choice == 1: 
        work = cp*delT/etaM
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

    eta = float(input("\nEnter Turbine Efficiency: "))
    print("\nIs the turbine pressure ratio known?")
    choice = int(input("(1) yes\n(2) no\n"))

    if choice == 1: 
        pRatio = float(input("\nEnter Turbine pressure ratio: "))
        exp = (gamma-1/gamma)
        delT = eta*tIn*(1-(1/pRatio)**exp)
    
    else: 
        if cWork == 0: cWork = float(input("\nEnter the compressor efficiency: "))

        delT = cWork/cp
        pRatio = (eta*tIn/(eta*tIn-delT))**(gamma/(gamma-1))

    tOut = delT*tIn

    pOut = pRatio*pIn
    work = cp*delT*eta

    print(f'\nP out = {pOut}')
    print(f'delta T = {delT}')
    print(f'Temp Out = {tOut}')
    print(f'Turbine Work = {work}\n')

    return pOut, tOut, work
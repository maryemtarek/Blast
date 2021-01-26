
DictPAM = {}
Proteins = ['G', 'A', 'V', 'L', 'I', 'P', 'S', 'T', 'D', 'E', 'N', 'Q', 'K', 'R', 'H', 'F', 'Y', 'W', 'M', 'C', 'B', 'Z', 'X']

def PAMMatrix():
    with open("PAM.txt", 'r') as file:
        lines = file.read().splitlines()
        PAM = []
        for line in lines:
            PAMList = [int(i) for i in line.split()] 
            PAM.append(PAMList)
        return PAM

def generateDic():
    PAM_MATRIX = PAMMatrix()

    for i in range(0, len(Proteins)):
        DictPAM[Proteins[i]] = {}
        Helper_Dict = {}
        for j in range(0, len(Proteins)):
            Helper_Dict[Proteins[j]] = PAM_MATRIX[i][j]

        DictPAM[Proteins[i]] = Helper_Dict
    
    return DictPAM

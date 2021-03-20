import LoadingPAM

wordLength = 3
Proteins = LoadingPAM.Proteins
DictPAM=LoadingPAM.DictPAM

def words(Query):
    wordsList = []
    #print (len(Query))

    for i in range(0, len(Query) - (wordLength - 1)):
        word = Query[i:wordLength + i]
        wordsList.append([word, i])
    return wordsList


def Neighboorhood(Words_List):
    Threshold = 20
    Seeds = []
    for i in range(0, len(Words_List)):
        CompareString = list(Words_List[i][0])
        s = list(Words_List[i][0])
        for k in range(0, len(s)):
            for j in range(0, len(Proteins)):
                CompareString[k] = Proteins[j]
                Score = DictPAM[s[0]][CompareString[0]] + DictPAM[s[1]][CompareString[1]] + DictPAM[s[2]][CompareString[2]]

                if Score >= Threshold:
                    Seeds.append(["".join(CompareString), Words_List[i][1], Score]) # Seeds 0--> word itself 1--> start ind in query  2--> Score

            CompareString = s

    return Seeds

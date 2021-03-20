
import Read
import LoadingPAM

seq=Read.Query
minNumOfRepeats = 3
minRepeatLength = 2

def repeats():
    #seq = Query
    maxRepeatLength = int(len(Read.Query) / minNumOfRepeats)
    for x in range(minRepeatLength, maxRepeatLength + 1):
        ProteinRepeatedSeq(x)
    return seq


def ProteinRepeatedSeq(repSeqLength):
    repeatedDict = {}
    repeat = 0
    start = -1
    ctr = 0  # counts the number of repeats
    i = 0
    while i < len(seq) - repSeqLength:
        j = i + repSeqLength
        if seq[i] == seq[j] or int(LoadingPAM.DictPAM[seq[i]][seq[i]]) == int(LoadingPAM.DictPAM[seq[j]][seq[j]]):
            if ctr == 0:
                ctr = 1
                start = i
            if j - start == repSeqLength * minNumOfRepeats - 1:
                ctr = minNumOfRepeats
                remove = seq[start:start + repSeqLength]

        elif ctr < minNumOfRepeats:
            start = -1
            ctr = 0

        else:
            repeatedDict[start] = j - 1
            repeat = 1
            start = -1
            ctr = 0

        i = i + 1

    if ctr >= minNumOfRepeats:
        repeatedDict[start] = len(seq) - 1
    if repeat == 1:
        RemoveRepeated(repeatedDict)


def RemoveRepeated(repeatedDict):
    global seq
    x = 0
    y = 0
    cnt = 0
    for key in repeatedDict:
        cnt = cnt + (y - x)
        seq = seq[:key - cnt] + seq[repeatedDict[key] + 1 - cnt:]
        x = key
        y = repeatedDict[key] + 1


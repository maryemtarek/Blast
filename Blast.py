import LoadingPAM
DictPAM = LoadingPAM.generateDic()
Proteins = LoadingPAM.Proteins

with open("Sequences.txt", 'r') as file:
     Query = file.read()

Genes = ["PYGPWGETPCGWFRRQGEHADGKRGEFQWEAAAAAPWG"]
def words(Query):
    wordLength = 3
    minNumOfRepeats = 3
    minRepeatLength = 2
    wordsList = []
    print (len(Query))

    for i in range(0, len(Query) - (wordLength - 1)):
        word = Query[i:wordLength + i]
        wordsList.append([word, i])
    return wordsList




def Neighboorhood(Words_List):
    Threshold = 21
    Seeds = []
    for i in range(0, len(Words_List)):
        CompareString = list(Words_List[i][0])
        s = list(Words_List[i][0])
        for k in range(0, len(s)):
            for j in range(0, len(Proteins)):          
                CompareString[k] = Proteins[j]            
                Score = DictPAM[s[0]][CompareString[0]] + DictPAM[s[1]][CompareString[1]] + DictPAM[s[2]][CompareString[2]]                      
            
                if Score >= Threshold:
                    Seeds.append(["".join(CompareString), Words_List[i][1], Score])

            CompareString = s

    return Seeds




 # Compare seeds with gene seq and add it to hit list  
def hit(seeds,gene):
    hitList=[]
    for i in range(0, len(seeds)):
        myWord = seeds[i][0]
        for j in range(0, len(gene)):
            if gene[j:j+3]==myWord: 
                hitList.append([i, j]) #i index in seeds list, j start hit in Gene
    return  hitList
 

#print(hit_list)

def Extend(hit_list,gene):
    hsp=[]
    # dec accuracy to inc time 
    threshold=4 
    for i in range(0, len(hit_list)):
        #hit_list[i][j]  j--> 0 index in seeds list, 1 start hit in Gene
        #      j=1
        #   3   (4   5   6)   7
        #  left (   hit   )      right
        leftIndGene =hit_list[i][1]-1 #first index to be extended from left
        rightIndGene = leftIndGene+4 #first index to be extended from right

        seedInd = hit_list[i][0] #Get index of seed in seeds list
        leftIndQuery = seeds[seedInd][1]-1
        rightIndQuery = leftIndQuery +4
        
        score = seeds[seedInd][2]
        mxScore = score
        while(mxScore-score<threshold): #terminate when score decrese 
    
            if leftIndGene>=0 and leftIndQuery>=0:
                score+=DictPAM[gene[leftIndGene]][Query[leftIndQuery]]
                leftIndQuery-=1
                leftIndGene-=1
            if rightIndGene<len(gene) and rightIndQuery<len(Query):
                score+=DictPAM[gene[rightIndGene]] [Query[rightIndQuery]]
                rightIndGene+=1
                rightIndQuery+=1           
            #cornern cases:
            if leftIndQuery<0 :
                # rightIndQuery make index in query boundries 
                # rightIndGene : query --> PQGE......
                #                gene  --> .......PQG
                if rightIndQuery>=len(Query) or rightIndGene>=len(gene):
                    break
            if leftIndGene<0 :
                # rightIndQuery make index in query boundries 
                #  rightIndQuery: query--> .......PQG
                #                 gene --> PQGE......
                if rightIndGene>=len(Query) or rightIndQuery>=len(Query):
                    break
            if mxScore<score:
                mxScore=score
        #         left                  right
        # finish: -1  0  1  ....... len  len+1
        hsp.append([rightIndGene-1,leftIndGene+1,rightIndQuery-1,leftIndQuery+1,score]) # +1,-1 for list boundries
    return hsp



def CheckOverlap(hsp,gene):
    newHspListFinal=[]
    ind=0
    size =len(hsp)
    while (ind<size): # iterate on items to remove common overlap
        # To add last element in list if not added
        if ind == len(hsp)-1 and newHspListFinal.count(hsp[ind])==0: # .count() check if element exists or not
            print (newHspListFinal.count(hsp[ind])==0)
            newHspListFinal.append(hsp[ind])
            break

        # F: First 
        # S: Sec
        FrightIndGene=hsp[ind][0]
        FlefttIndGene = hsp[ind][1]

        SrightIndGene = hsp[ind + 1][0]
        SleftIndGene = hsp[ind + 1][1]

        # print(FrightIndGene, SleftIndGene, ind)
        # if it overlap in gene
        if FrightIndGene >=SleftIndGene:

            FrightIndQuery = hsp[ind][2]
            FleftIndQuery = hsp[ind][3]
    
            SrightIndQuery = hsp[ind + 1][2]
            SleftIndQuery = hsp[ind + 1][3]
        
            #print(FrightIndGene, SleftIndGene, ind)
            # if it overlap in query
            if FrightIndQuery >= SleftIndQuery:
                # seq1--> PQGE
                # seq2--> PQGE
                if FleftIndQuery==SleftIndQuery and FrightIndQuery==SrightIndQuery:
                    newHspListFinal.append(hsp[ind])
                else:
                  
                    if FleftIndQuery == SleftIndQuery:
                        # seq1--> PQGEWQ
                        # seq2--> PQGE
                        if FrightIndQuery > SrightIndQuery:
                            newHspListFinal.append(hsp[ind])
                        # seq1--> PQGE
                        # seq2--> PQGEWQ
                        else:
                            newHspListFinal.append(hsp[ind+1])
                    # seq1--> PQGEWQ
                    # seq2-->    EWQPQF
                    else:        
                        newScore= hsp[ind][4]+hsp[ind+1][4]
                        geneCounter=SleftIndGene
                        # Remove duplicated scores from overlap
                        for j in range(SleftIndQuery, FrightIndQuery):
                            newScore -= DictPAM[gene[geneCounter]][Query[j]]
                            #print(newScore)
                            geneCounter+=1
                        
                        newHspListFinal.append([SrightIndGene, FlefttIndGene, SrightIndQuery, FleftIndQuery , newScore])                
                #   ind   ind+1   ind+2
                #    2      3       4
                # Move to seq2(ind+1) in list
                ind=ind+1 
            else: # No overlap in query
                newHspListFinal.append(hsp[ind])
        else: # No overlap in gene
            newHspListFinal.append(hsp[ind])
        ind+=1 # Move to nex element / (ind+2) in overlap
    return newHspListFinal



def PrintFinalList(Final_HSP_List,gene):
    for i in range (0,len(Final_HSP_List)):
        index1= Final_HSP_List[i][0]
        index2 =Final_HSP_List[i][1]
        print ("**************************************************************")
        print ("Gene :",gene[index2: index1+1])        
        print ("Query:", Query[Final_HSP_List[i][3]:Final_HSP_List[i][2]+1])    
        print("Score:",Final_HSP_List[i][4])
        print ("**************************************************************")
    

word_list = words(Query)
seeds = Neighboorhood(word_list)
for gene in Genes:
    hit_list = hit(seeds,gene)
    hsp_List = Extend(hit_list,gene)
    hsp_List.sort(key = lambda x: x[1]) # sort hsp according to leftIndQuery
    Final_HSP_List = CheckOverlap(hsp_List,gene)
    PrintFinalList(Final_HSP_List,gene)



#Refnces
# https://w...content-available-to-author-only...s.org/python-sort-list-according-second-element-sublist/

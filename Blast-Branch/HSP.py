import LoadingPAM
DictPAM=LoadingPAM.DictPAM

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

def Extend(hit_list,gene,Query,seeds):
    hsp=[]
    # dec accuracy to inc time
    threshold=4
    for i in range(0, len(hit_list)):
        #hit_list[i][j]  j--> 0 index in seeds list(word, start ind,score), 1 start hit in Gene
        #      j=1
        #   3   (4   5   6)   7
        #  left (   hit   )   right
        leftIndGene =hit_list[i][1]-1 #first index to be extended from left
        rightIndGene = leftIndGene+4 #first index to be extended from right

        seedInd = hit_list[i][0] #Get index of seed in seeds list
        leftIndQuery = seeds[seedInd][1]-1
        rightIndQuery = leftIndQuery +4

        score = seeds[seedInd][2]
        mxScore = score
        while(mxScore-score<threshold): #terminate when score decrease less than 4 of maxScore

            if leftIndGene>=0 and leftIndQuery>=0:
                score+=DictPAM[gene[leftIndGene]][Query[leftIndQuery]]
                leftIndQuery-=1
                leftIndGene-=1
            if rightIndGene<len(gene) and rightIndQuery<len(Query):
                score+=DictPAM[gene[rightIndGene]] [Query[rightIndQuery]]
                rightIndGene+=1
                rightIndQuery+=1
            #corner cases:
            if leftIndQuery<0 :
                # rightIndQuery make index in query boundaries
                # rightIndGene : query --> -1 PQGE -1 # Score all query
                #                gene  --> ...PQGE...
                # rightIndGene : query -->        PQGE......
                #                gene  --> .......PQG
                if rightIndQuery>=len(Query) or rightIndGene>=len(gene):
                    break
            if leftIndGene<0 :
                # rightIndQuery make index in query boundaries
                #  rightIndQuery: query--> .......PQG
                #                 gene -->        PQGE......
                if rightIndGene>=len(Query) or rightIndQuery>=len(Query):
                    break
            if mxScore<score:
                mxScore=score
        #         left                  right
        # finish: -1  0  1  ....... len  len+1
        hsp.append([rightIndGene-1,leftIndGene+1,rightIndQuery-1,leftIndQuery+1,score]) # +1,-1 for list boundaries
    return hsp



def CheckOverlap(hsp,gene,Query):
    newHspListFinal=[]
    ind=0
    size =len(hsp)
    while (ind<size): # iterate on items to remove common overlap
        # To add last element in list if not added
        if ind == len(hsp)-1 and newHspListFinal.count(hsp[ind])==0: # .count() check if element exists or not
            newHspListFinal.append(hsp[ind])
            break

        # F: First
        # S: Sec
        FrightIndGene=hsp[ind][0]
        FleftIndGene = hsp[ind][1]

        SrightIndGene = hsp[ind + 1][0]
        SleftIndGene = hsp[ind + 1][1]

        print(FrightIndGene, SleftIndGene, ind)

        # if it overlap in gene
        # seq1--> PQGEUSPQ
        # seq2-->     RLQGESPQ
        
        if FrightIndGene >=SleftIndGene :


            hsp1 = hsp[ind][2]
            hsp2 = hsp[ind+1][2]
            if hsp1 < hsp2:
                FleftQuery = hsp1
                FrightQuery = hsp[ind][3]
                SleftQuery = hsp2
                SrightQuery = hsp[ind+1][3]
            else: 
                SleftQuery = hsp1
                SrightQuery = hsp[ind][3]
                FleftQuery = hsp2
                FrightQuery = hsp[ind+1][3]



            if FrightQuery >= SleftQuery:
                leftQuery = min (FleftQuery,SleftQuery)
                rightQuery = max (FrightQuery,SrightQuery)
                leftGene = min (FleftIndGene,SleftIndGene)
                rightGene = max(FrightIndGene,SrightIndGene)
                geneCounter = rightGene
                newScore = hsp[ind][4]+hsp[ind+1][4]
                # Remove duplicated scores from overlap
                for j in range(leftQuery, rightQuery):
                    newScore -= DictPAM[gene[geneCounter]][Query[j]]
                    geneCounter +=1
                newHspListFinal.append ([rightGene,leftGene,rightQuery,leftQuery,newScore])
                ind +=1

            elif FleftQuery <=SleftQuery and FrightQuery >=SrightQuery:
                leftGene = min (FleftIndGene,SleftIndGene)
                rightGene = max(FrightIndGene,SrightIndGene)
                newHspListFinal.append()
                newHspListFinal.append ([rightGene,leftGene,FleftQuery,FrightQuery,hsp[ind][4]])
                ind +=1
                
            else : #No Overlap in Query
                newHspListFinal.append(hsp[ind])

        else: # No Overlap in gene
            newHspListFinal.append(hsp[ind])
        ind+=1 # Move to nex element / (ind+2) in overlap




    return newHspListFinal



def PrintFinalList(Final_HSP_List,gene,Query):
    for i in range (0,len(Final_HSP_List)):
        index1= Final_HSP_List[i][0]
        index2 =Final_HSP_List[i][1]
        print ("**************************************************************")
        print ("Gene :",gene[index2: index1+1])
        print ("Query:", Query[Final_HSP_List[i][3]:Final_HSP_List[i][2]+1])
        print("Score:",Final_HSP_List[i][4])
        print ("**************************************************************")

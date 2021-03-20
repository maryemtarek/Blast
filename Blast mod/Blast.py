
import LoadingPAM
import RemoveRepeats
import Read
import GeneratingSeeds
import HSP

DictPAM = LoadingPAM.generateDic()
Proteins = LoadingPAM.Proteins

print("Query before removing repeats :",Read.Query)
Query=RemoveRepeats.repeats()
print("Query after removing repeats :",Query)
print ("----------------------------------------------------------------")


word_list = GeneratingSeeds.words(Query)
seeds = GeneratingSeeds.Neighboorhood(word_list)
cnt =1
for gene in Read.Genes:
    print ("In Gene #",cnt," in Database: ")
    hit_list = HSP.hit(seeds,gene)
    hsp_List = HSP.Extend(hit_list,gene,Query,seeds)
    hsp_List.sort(key = lambda x: x[1]) # sort hsp according to leftIndGene
    Final_HSP_List = HSP.CheckOverlap(hsp_List,gene,Query)
    #print(Final_HSP_List)
    HSP.PrintFinalList(Final_HSP_List,gene,Query)
    cnt+=1



#Refnces
# https://w...content-available-to-author-only...s.org/python-sort-list-according-second-element-sublist/

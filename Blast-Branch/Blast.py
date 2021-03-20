
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

x="MECWYWQCKKNKGNPCGINFDWGENTNISIRWQISFIRWYLVTQFHRNYQFIEPQSTAPYSMRRTNYEIFGFRAKMHPHKCTCSVWHLVHRLVQETSNIAMYSIRSPVNERNMYWVQNHLHAYVIDNHRDWAAEHVKCKRNYSLNSNSIWIGGLEINIMECQAQMCCYMEPPMIVTHCTWQTQQDPLISNDAGPCRGHTTCQRQADRVQSPQLIKIPRMWGLCLCCYRLTDCHMGPDPHKGEHFAFYIEHASMKDAGMCPTKGLPVRRYMWDPYMALKSHQDTDDNDGTQDPEYWVPAFPDTTYQFPRDRRTTIHDQHSFMVWAMGCKWGLISDIIYEKYHCANRGIYVFKGAPHLCHPKAEDICGHYSFYRSQSRWHKSYTKVLDWNMPHYGDNFNMAWDPYGGAERQDFTEFKEWSISFYQDWKNSFDEDYRCKWTMERRLTYEPAGSSQQTCTRGAVTAPNFFPGHEDKHYWAMTPCRKPSETEKVLYLGCGERPCE"
for gene in range(410,420):
    print(gene," ",x[gene])



word_list = GeneratingSeeds.words(Query)
seeds = GeneratingSeeds.Neighboorhood(word_list)
print(seeds)
for gene in Read.Genes:

    hit_list = HSP.hit(seeds,gene)
    hsp_List = HSP.Extend(hit_list,gene,Query,seeds)

    print (gene[415],gene[416],gene[417])
    hsp_List.sort(key = lambda x: x[1]) # sort hsp according to leftIndGene
    print(hsp_List)
    Final_HSP_List = HSP.CheckOverlap(hsp_List,gene,Query)

    print(Final_HSP_List)
    HSP.PrintFinalList(Final_HSP_List,gene,Query)



#Refnces
# https://w...content-available-to-author-only...s.org/python-sort-list-according-second-element-sublist/

from itertools import product
import numpy as np
import pdb
import operator
from math import comb


## This script is to find out the best alignment for multiple DNA/RNA
## sequences using dynamic programming. This script allows the user 
## to define cost of mismatch(pxy), gap opening(pgap) and gain of matching(mRd). 

def neighbors(index):
    for relative_index in product((0, -1), repeat=len(index)):
        if not all(i == 0 for i in relative_index):
            yield tuple(i + i_rel for i, i_rel in zip(index, relative_index))

# Returnc index of neighbouring nucleotides in all N input sequences.
# For N input sequences, an index would be a tuple of N elements. 
# Get index of neighbouring bp of all sequence. 
def getNeighbour(tmp):
    # filter out indexes that are out of bound
    temp = list(neighbors(tmp))
    output = []

    l0, l1 = 1,0
    while(l0!=l1):
        for tp in temp:
            remov = False
            l0 = len(temp)
            # check if there is a base pair at tp in ith sequence 
            if any( tp[i] <0 for i in range(0,len(list(tmp)))):
                remov = True
            if remov:
                temp.remove(tp)
            else :
                output.insert(0,tuple(tp))
            l1 = len(temp)
    return(temp)

def dpMSA(sqs, pxy, pgap, mRd):
    # build score matrix twice as gene length to avoid 
    lengths = np.array([len(sequence) for sequence in sqs]);
    dim = []
    # define the matix row/column number twice of the longest input sequence
    for length in lengths:
        dim.append(max(lengths)*2)
    # scoreMatrix = np.zeros([max(lengths)*2,max(lengths) *2]);
    scoreMatrix = np.zeros(dim)

    
    # initialise the matrix by adding default gap cost
    for j in range(0,len(lengths)):
        for i in range(1,len(dim)+1):
            tupleIndexList = [0]*len(lengths)
            tupleIndexList[j] = i
    scoreMatrix[ tuple(tupleIndexList) ] = 0-i * pgap;
    scoreMatrix[ tuple([0,0,0]) ] = 0
    # scoreMatrix = np.full(scoreMatrix0.shape, np.NaN)
    
    # for the easy of access, add* in front of sequence to correct for index
    sqs_original = sqs
    sqs = []
    for sq in sqs_original:
        sqs.append("*"+sq)

   # iterate through index at all dimension and calculate the score
    lengths = tuple(lengths)
    ranges = [tuple(range(1, length+1)) for length in lengths[::-1]];
    for tupleIndex in product(*ranges):
        # if tupleIndex != lengths :
            tupleIndex = tuple(tupleIndex[::-1])
            # generate relative movement by 1 in all dimension 
            # neighborIndexes = getNeighbour(tupleIndex,sqs)
            neighborIndexes = list(neighbors(tupleIndex))
            costByNeighbor = []
            # print("Number of possible neighbor : " + str(len(neighborIndexes)))
            # looping over all indexes to find maximum score increase
            for neighborIndex in neighborIndexes:
                neighborIndex = tuple(neighborIndex)

                temp = 0
                relative_index =[]
                for i in range(0,len(tupleIndex)):
                    relative_index.append(list(tupleIndex)[i]-list(neighborIndex)[i])
                # get alignment by base pair according to index
                compare = 0
                gap_call = 0
                # print("#### calculating scores ###########")
                for i in range(0,len(tupleIndex)):
                    for j in range(0,len(tupleIndex)):
                        bp= sqs[i][list(neighborIndex)[i] -1]            
                        bp2=sqs[j][list(neighborIndex)[j] -1]

                        # since all comparison is done twice 
                        # cut the score by half before adding
                        if i == j:
                            # print("skip self-comparing")
                            print("",end="")
                        # no indel
                        elif relative_index[j] == 1 :
                            compare = compare +1                            
                            if sqs[i][neighborIndex[i]-1] == sqs[j][neighborIndex[j]-1] :
                                # print("    match |bp : " + bp + " bp2 : "+bp2)
                                temp = temp + mRd/2
                            else:
                                # print(" mismatch |bp : " + bp + " bp2 : "+bp2)
                                temp = temp - pxy/2
                        # indel
                        elif relative_index[j] == 0 :
                            # print(" gap  : " + bp + " bp2 : "+bp2)
                            gap_call = gap_call + 1
                            temp = temp - pgap/2
                        else:
                            print("relative index !=1, WRONG")
                
                # print("Compare : " + str(compare) + " Gap_call : " + str(gap_call))
                # for all element in neighborIndex calculated twice
                costByNeighbor.append(temp)
            # update scoreMatrix by all 
            chooseMax = []
            for i in range(0,len(neighborIndexes)):
                pS = scoreMatrix[ tuple(neighborIndexes[i])]
                if np.isnan(pS) : 
                    pS = 0
                chooseMax.append( pS + costByNeighbor[i])
            # print("chooseMax : ",end="")
            # print(chooseMax)
            i = chooseMax.index(max(chooseMax))
            # update scoreMatrix
            scoreMatrix[tupleIndex] = max(chooseMax)
    return (scoreMatrix)

# This function print the numeric score matrix that define the best alignment.
# The result would be printed to console as vectors of N integers. The number
# of vectors is the length of aligment. 
def printTrack(sqs, scoreMatrix,pxy,pgap,mRd):
    gene_index_mov = np.array([len(sequence)-1 for sequence in sqs])

    tmp = np.transpose(np.nonzero(score_matrix))
    tupleIndex = tuple(tmp[len(tmp)-1])
    zeroTuple = tuple([0]*len(sqs) )
    # print("tupleIndex  & zeroTuple : ", end="")
    # print(tupleIndex, end=" & ")
    # print(zeroTuple)

    tuple_array = []
    sqs_original = sqs
    sqs = []
    for sq in sqs_original:
        sqs.append("*"+sq)
    while( tuple(tupleIndex) != tuple(zeroTuple) ):
        tuple_array.append(tupleIndex)
        # print("tupleIndex\t",end="")
        # print(tupleIndex)
        neighborIndexes = list()
        neighborScores = list()
        neighborIndexes = getNeighbour(tupleIndex)
        for neighborIndex in neighborIndexes:
            neighborScores.append(scoreMatrix[neighborIndex])
        
        which_neighbour = neighborScores.index(max(neighborScores))
        # print("neighborIndex at : " + str(which_neighbour))
        
        tupleIndex = neighborIndexes[which_neighbour]
    # print(tuple_array)
    return(tuple_array)
    
### test run
sequences =  ["GCTGGAAGGCAT","GCAGAGCACG","GCTGGAAGCAT"]

# (sqs, pxy, pgap, mRd)
score_matrix = dpMSA(sequences,3,4,5)
aPath = printTrack(sequences,score_matrix,3,4,5)
# back track the path to print sequences
track = []
for item in aPath :
    track.insert(0,item)
track.insert(0,tuple([0]*len(sequences)))

sqs = []
for sq in sequences:
    sqs.append("*"+sq)

max_seq_len = (max(list(max(track))))

# print( len(track))  # bp
# print(len(track[0]))# seq #
out = np.empty([len(track[0]),len(track)], dtype=str)

# print(track)
track_bk = track
track_wk = []
print(" track length : " , len(track))

# for i in range(0,len(track[0])):
#     print("For the ", i,"th sequence")
i = 0
for j in range(0,len(track)-1):
    # print("At base pair: ", j)
    arr = np.asarray(track[j]).tolist()
    if (list(track)[j][i]  ==  list(track)[j+1][i]) :
        arr[i] = (-1)
    track_wk.append(arr)

for i in range(1,len(track[0])):
    print("For the ", i,"th sequence")
    for j in range(0,len(track)-1):
        # print("At base pair: ", j)
        arr = np.asarray(track_wk[j]).tolist()
        # if ( j== len(track) ):            
        #     if (list(track)[j][i]  ==  list(track)[j-1][i]) :
        #         arr[i] = (-1)
        # else :
        if (list(track)[j][i]  ==  list(track)[j+1][i]) :
            arr[i] = (-1)
        track_wk[j] = arr

for arr in track_wk:
    print(arr)


# i for index of sequence in input sequence list, j for base pair index
for i in range(0,len(track[0])):
    for j in range(1,len(track)-1):
        nucleotide = sqs[i][ list(track)[j][i] ]   
        out[i,j] = nucleotide     
        if (list(track)[j][i]  ==  list(track)[j+1][i]) :
            out[i,j] = "-"  

# print all sequences 
for o in out:
    print("".join(o))



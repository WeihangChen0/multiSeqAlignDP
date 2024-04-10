# multiSeqAlignDP
This is a python script that does multiple sequence alignment using dynamic programming. The dynamic programming method for multiple alignment is costly and hard to implement but it will return the best solution given the paramenter setting of cost and reward. Therefore, it can be used to test out performance of other alignment algorithm using limited number of short sequences. 

The problem can be divided into several steps:
1)	Build a zero matrix in N dimensions to store scores of alignment.
2)	Find the coordinates of a neighbor that is 1 base pair away in all dimensions.
3)	Calculate the cost by adding total match reward and mThe problem can be divided into several steps. 
4)	Put the score of the neighbor, i.e score of this position + corresponding cost in the neighbor’s coordinate.
5)	Find the max score among the neighbor and coordinate of that neighbor.
6)	Move to that neighbor and repeat step 2-5 until reaching the end of scoring matrix.
7)	Report the movement of coordinates of the path as profile.

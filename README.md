# Similarity-Search
<!--### Works:
1. Verify that the performance of HNSW (KNN Graph) is better
	* The distribution of the number of items in a PQ bucket (may be unbalanced)	* The distribution of the distance from an item to its bucket (for KNN graph, for ending node to start node) 	* The distribution of the distance between the items in the same bucket	* Measure how recall changes after a neighbor in top 50 has been found for PQ and KNN graph2. Verify our idea of improving the connectivity of KNN graphWalk a KNN graph and PQ IMI simultaneously
	* choose the one with the smallest distance to proceed. 
	* Build a PQ that contains approximately (n/k) buckets connect each PQ bucket to its k nearest neighbors in the dataset as if the buckets are artificial data introduced in the dataset
	* find those items in the dataset that are not connected to any PQ code, connect them to their nearest bucket.
	* Question: how many items need to be processed in the last stage? How does this improve time recall?
3. Show the challenge of using KNN graph for MIPSPlot the in-degree distribution of the items when using MIP and L2 distance for graph building  -->  
<!--### 10.15
Fisrt of all, I have read some papers related to KNN graph, such as HNSW and NNDES. I also read papers about PQ. I learned about the idea of the Norm-layered KNN graph for MIPS and recently I worked on some experienments to compare performance of L2 metric  and Innerproduct on HNSW and NNDES. Liu jie is programming for this algorithm and I will also try to help with him.

### 11.7
- [X] 把没用的code删掉
- [X] 跑yahoo和imagenet   
- [ ] 改stop condition
- [ ] 保存pq code to file
- [ ] test diff # of bridge vectors(1000)
	- netflix
		- non empty: 1632/16384
		- 500
		- 1000
		- 2000
	- yahoomusic
		- non empty: 4046/65536
		- 500 
		- 1000
		- 2000
	- imagenet
		- 423846/1048576
		- 500
		- 1000
		- 2000
		- 4000
		- 8000-->
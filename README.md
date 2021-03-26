# Personalized PageRank

## Tested Environment
- Ubuntu 20.04
- GCC 9.3.0
- SFMT 1.5.1 
- cmake 3.16

## Compile
Please make a folder with the name *build* first, if it does not exist. 
In Ubuntu, you can make it with the following command: 
```sh
$ mkdir build
```

Then you can compile the code as follows: 
```sh
$ cmake -B build
$ cd build && make -j && cd ..
```

In what follows, we assume that the program is complied and there is an execution file "*SpeedPPR*" in the folder "*build/*".

## Data

Our program accepts raw-datasets in the following format:

1. It begins with zero or more comment lines, each of which starts with a symbol '#'.   
2. Each line following the initial comment lines contains a pair of endpoint ids of an edge.  

Be careful with the format. We don't provide sanity check of the input graph format. 

## Clean the graphs

We don't use the raw-datasets directly. This function removes vertices without incoming nor outgoing edge, and relabels the ids of the rest vertices, starting from integer 0.

```sh
$ build/SpeedPPR -algo CLEAN_GRAPH \
    -graph <raw-dataset-path> \
    -is_undirected <yes or no> \
    -output_folder <output-folder-path>
```

For undirected graphs, our program duplicates each undirected edge into two directed edges. We assume that each edge appears only once in the original graph.

### *Example*

For example, if we have an *undirected* graph DBLP, stored in *"DataSetRaw/dblp/com-dblp.ungraph.txt"*, and we want to keep the output graph in the folder "*DataSetCleaned/dblp*", then we can run the following command.

```sh
$ build/SpeedPPR -algo CLEAN_GRAPH \
    -graph DataSetRaw/dblp/com-dblp.ungraph.txt \
    -is_undirected yes \
    -output_folder DataSetCleaned/dblp  
```
After the program returns, there should be 3 output files in "*DataSetCleaned/dblp*". 

1. **attribute.txt** : a meta file that records the number of vertices and edges in the output graph. 
   
2. **graph.txt** : the output graph in txt form. The first line contains a number, which is the number of vertices. Each of the following lines represents an edge by using the ids of its endpoints.  This is outputted only for reference.
   
3. **graph.bin** : the output graph in binary form. Each pair of unsigned integers represents the ids of an edge's endpoints. Note that the number of vertices is not kept as the first number in this file. It is the format that would be used by our program.

If we have a directed graph, e.g., "*DataSetRaw/web-Stanford/web-Stanford.txt*", we can run the following command. 

```sh
$ build/SpeedPPR -algo CLEAN_GRAPH \
    -graph DataSetRaw/web-Stanford/web-Stanford.txt \
    -is_undirected no \
    -output_folder DataSetCleaned/web-Stanford 
```

## Generate queries
Generate random source ids for single source PPR queries. Each line in the output file contains a node id sampled uniformly at random.

```sh
$ build/SpeedPPR -algo GEN_QUERY \
    -meta <graph-attribute-file> \
    -graph_binary <binary-graph-file> \
    -query <output-file-name> \
    -query_size <query count> 
```

### *Example*

```sh
# sample queries for dblp
$ build/SpeedPPR -algo GEN_QUERY \
    -meta DataSetCleaned/dblp/attribute.txt \
    -graph_binary DataSetCleaned/dblp/graph.bin \
    -query DataSetCleaned/dblp/dblp.query \
    -query_size 30 

# sample queries for web-Stanford
$ build/SpeedPPR -algo GEN_QUERY \
    -meta DataSetCleaned/web-Stanford/attribute.txt \
    -graph_binary DataSetCleaned/web-Stanford/graph.bin \
    -query DataSetCleaned/web-Stanford/web-Stanford.query \
    -query_size 30  
```

## Multi-thread Ground Truth  

We are ready to compute the ground truth PPR for the given query ids, in multi-thread. The solutions may be needed for later testing.  

```sh
$ build/SpeedPPR -algo GROUND_TRUTH \
    -meta <graph-attribute-file> \
    -graph_binary <binary-graph-file> \
    -query <query-file> \
    -query_size <query count> \
    -answer_folder <folder to save ground truth> \
    -alpha <alpha> \
```

### *Example*

In this example, the solutions are saved in the folder *answer/*. Create such a folder if necessary. Also, we might need to create sub-folders *dblp/* and *web-Stanford/* under the folder *answer/*. 

```sh
# crate the folder answer only if it does not exist
$ mkdir answer

# create the folders answer/dblp and answer/web-Stanford only if they do not exist
$ mkdir answer/dblp && mkdir answer/web-Stanford

# compute ground truth for dblp
# we omit the option -alpha <alpha>, by default alpha=0.2
$ build/SpeedPPR -algo GROUND_TRUTH \
    -meta DataSetCleaned/dblp/attribute.txt \
    -graph_binary DataSetCleaned/dblp/graph.bin \
    -query DataSetCleaned/dblp/dblp.query \
    -query_size 30 \
    -answer_folder answer/dblp 

# compute ground truth for web-Stanford
$ build/SpeedPPR -algo GROUND_TRUTH \
    -meta DataSetCleaned/web-Stanford/attribute.txt \
    -graph_binary DataSetCleaned/web-Stanford/graph.bin \
    -query DataSetCleaned/web-Stanford/web-Stanford.query \
    -query_size 30 \
    -answer_folder answer/web-Stanford   
```


## High-Precision PPR

Single thread high precision PPR computation. This function returns a PPR vector whose $\ell_1$ error (with respect to the true PPR vector) is within a specified threshold. If this threshold is not specified, by default it is set to *1 / {\# edges in the graph}*.

We provide three algorithms, namely

1. *PowForPush*: the *PowerPush* algorithm proposed in our paper.  
2. *PowItr*: the vanilla power iteration.  
3. *FwdPush*: forward push with First-In-First-Out queue.  

The solution output folder is specified by the option *-estimation_folder*, to distinguish the output option of *GROUND_TRUTH*.

```sh
$ build/SpeedPPR -algo <PowForPush / PowItr / FwdPush> \
    -meta <graph-attribute-file> \
    -graph_binary <binary-graph-file> \
    -query <query-file> \
    -query_size <query count> \
    -estimation_folder <result-output-folder> \
    -alpha <alpha> \
    -l1_error <l1-error>
```



### *Example*

```sh
# crate the folder estimation only if it does not exist
$ mkdir estimation

# crate the folder estimation/exact only if it does not exist
$ mkdir estimation/exact

# create the folders estimation/exact/dblp and estimation/exact/web-Stanford only if they do not exist
$ mkdir estimation/exact/dblp && mkdir estimation/exact/web-Stanford

# compute high precision PPR for dblp. PowForPush can be replaced by PowItr or FwdPush
$ build/SpeedPPR -algo PowForPush \
    -meta DataSetCleaned/dblp/attribute.txt \
    -graph_binary DataSetCleaned/dblp/graph.bin \
    -query DataSetCleaned/dblp/dblp.query \
    -estimation_folder estimation/exact/dblp \
    -query_size 30 \
    -alpha 0.2 \
    -l1_error 1e-8 

# compute high precision PPR for web-Stanford. PowForPush can be replaced by PowItr or FwdPush
$ build/SpeedPPR -algo PowForPush \
    -meta DataSetCleaned/web-Stanford/attribute.txt \
    -graph_binary DataSetCleaned/web-Stanford/graph.bin \
    -query DataSetCleaned/web-Stanford/web-Stanford.query \
    -estimation_folder estimation/exact/web-Stanford \
    -query_size 30 \
    -alpha 0.2 \
    -l1_error 1e-8 
```


## Index
Generate index for approximate PPR queries, used by the index version of SpeedPPR. 

```sh
$ build/SpeedPPR -algo BUILD_INDEX \
    -meta <graph-attribute-file> \
    -graph_binary <binary-graph-file> \
    -index_file <output-index-file-name> \
    -alpha <alpha> 
```

### *Example*

```sh
# build index for dblp
$ build/SpeedPPR -algo BUILD_INDEX \
    -meta DataSetCleaned/dblp/attribute.txt \
    -graph_binary DataSetCleaned/dblp/graph.bin \
    -index_file DataSetCleaned/dblp/index.bin \
    -alpha 0.2 

# build index for web-Stanford
$ build/SpeedPPR -algo BUILD_INDEX \
    -meta DataSetCleaned/web-Stanford/attribute.txt \
    -graph_binary DataSetCleaned/web-Stanford/graph.bin \
    -index_file DataSetCleaned/web-Stanford/index.bin \
    -alpha 0.2 
```

## Approximate PPR

This function answers approximate personalized PageRank queries. Given a source node *s*, for any node *t* in the graph, if the PPR *π(s, t) ≥ 1 / {\# vertices in the graph}*, it returns an estimate of *π(s, t)*, with relative error at most ε with probability *1 - 1 / {\# vertices in the graph}*.  

```sh
$ build/SpeedPPR -algo SpeedPPR \
    -meta <graph-attribute-file> \
    -graph_binary <binary-graph-file> \
    -index_file <index-file> \
    -query <query-file> \
    -query_size <query count> \
    -estimation_folder <result-output-folder> \
    -alpha <alpha> \
    -epsilon <relative-error> \
    -with_idx <no or yes> 
```

### *Example*  

```sh
# crate the folder estimation/index_free only if it does not exist
$ mkdir estimation/index_free

# create the folders estimation/index_free/dblp and estimation/index_free/web-Stanford only if they do not exist
$ mkdir estimation/index_free/dblp && mkdir estimation/index_free/web-Stanford

# without index dblp, the option "-with_idx no" can be dropped.
$ build/SpeedPPR -algo SpeedPPR \
    -meta DataSetCleaned/dblp/attribute.txt \
    -graph_binary DataSetCleaned/dblp/graph.bin \
    -query DataSetCleaned/dblp/dblp.query \
    -estimation_folder estimation/index_free/dblp \
    -query_size 30 \
    -alpha 0.2 \
    -epsilon 0.5

# without index web-Stanford, the option "-with_idx no" can be dropped.
$ build/SpeedPPR -algo SpeedPPR \
    -meta DataSetCleaned/web-Stanford/attribute.txt \
    -graph_binary DataSetCleaned/web-Stanford/graph.bin \
    -query DataSetCleaned/web-Stanford/web-Stanford.query \
    -estimation_folder estimation/index_free/web-Stanford \
    -query_size 30 \
    -alpha 0.2 \
    -epsilon 0.5 

# crate the folder estimation/with_index only if it does not exist
$ mkdir estimation/with_index

# create the folders estimation/with_index/dblp and estimation/with_index/web-Stanford only if they do not exist
$ mkdir estimation/with_index/dblp && mkdir estimation/with_index/web-Stanford

# with index dblp
$ build/SpeedPPR -algo SpeedPPR \
    -meta DataSetCleaned/dblp/attribute.txt \
    -graph_binary DataSetCleaned/dblp/graph.bin \
    -index_file DataSetCleaned/dblp/index.bin \
    -query DataSetCleaned/dblp/dblp.query \
    -estimation_folder estimation/with_index/dblp \
    -query_size 30 \
    -alpha 0.2 \
    -epsilon 0.5 \
    -with_idx yes
 
# without index web-Stanford
$ build/SpeedPPR -algo SpeedPPR \
    -meta DataSetCleaned/web-Stanford/attribute.txt \
    -graph_binary DataSetCleaned/web-Stanford/graph.bin \
    -index_file DataSetCleaned/web-Stanford/index.bin \
    -query DataSetCleaned/web-Stanford/web-Stanford.query \
    -estimation_folder estimation/with_index/web-Stanford \
    -query_size 30 \
    -alpha 0.2 \
    -epsilon 0.5 \
    -with_idx yes 
```

## Citation  

```c++
@article{WGWZ21,  
    title={Unifying the Global and Local Approaches: An Efficient Power Iteration with Forward Push},
    author={Wu, Hao and Gan, Junhao and Wei, Zhewei and Zhang, Rui},  
    journal={arXiv preprint arXiv:2101.03652}  
}
```




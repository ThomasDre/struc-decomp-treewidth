# Tree Decomposition Using Heuristics

## Requirements
- Boost Library
- R (for generating plots)

## Usage
treewidth (-v | -g <graph> [-e <min-deg,min-fill,max-card>] -t) -h
h ... Help
v ... run evaluation mode
g <graph> ... specify input graph
h <heuristic> ... specify either one of the available heuristics "min-deg", "min-fill", "max-card"
t ... generate a tree decomposition from the generated ordering

## How to use it (extended)

### Evaluation Mode
Using the -v flag, the program is run in the evaluation mode. 
It creates 100 random graphs using the Erdos-Renyi Model using the values n=[10,100,100] and p=[0.25,0.50,0.75]. 
All 3 heuristics are applied on the generated graphs and a tree decomposition is genereated for each. The aggregated treewidths for each pair of settings and each heuristic is saved
into a file which is stored in the results/eval folder.
The R file reads all results from this folder and generates the set of plots which are then saved in the plots folder.

Note: the R file needs to be executed manually

### Standard mode
Using a file that describes an input graph (file format https://github.com/ciaranm/glasgow-subgraph-solver?tab=readme-ov-file#file-formats) the following can be done:
- Generate a elimination ordering using on the available heuristics
- Also generate the tree decomposion for that graph and the specified elimination ordering heuristic

The results are stored under:
- results/eliminationOrderings: this will always be generated
- results/decompositions: the csv containing the tree decomposition will only be generated if the -t flag is set

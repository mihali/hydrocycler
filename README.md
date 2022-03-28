# hydrocycler

Usage: **hydrocycler.py file.xyz**

A program to find cycles of hydrogen bonding within a molecular cluster and to generate derivative molecular clusters by reversing the direction of the molecular- and H-bonding. These derivative molecular clusters serve as input for subsequent energy minimization. 

It views the covalent bonding and H-bonding network as a directed graph in the direction of H-bonding donation; that is 

 **[O-H> - - [O-H> - - [O-H> - - **

is viewed as a graph with three nodes in the left to right direction.  The result of a reversal of the above graph would look like

 **O] - - <H-O] - - <H-O] - - <H **

The input is a cartesian coordinate file and the outputs are cartesian coordinate files as well.                             
                               
## Installation

Procedure 
1. git clone https://github.com/mihali/hydrocycler.git 
2. cd hydrocycler
3. conda create --name hydrocycler python=3.7 
4. conda activate hydrocycler
5. conda install numpy scipy

## Test

There are files to try in the testfiles directory, e.g.

% python hydrocycler.py testfiles/si-4oh-27h2o.xyz

## References

1. Uses johnson.py from https://github.com/qpwo/python-simple-cycles
2. Algorithm employed is from Johnson (1975) https://doi.org/10.1137/0204007


  

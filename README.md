<img width="400" alt="Screen Shot 2022-03-30 at 9 59 51 AM" src="https://user-images.githubusercontent.com/10854556/161368022-17f12de8-4c98-483b-9340-e24c68c5e86b.png">


# hydrocycler

Usage: **python hydrocycler.py file.xyz**

An interactive script to find cycles of hydrogen bonding within an optimized molecular cluster and to generate derivative molecular clusters by reversing the direction of the molecular- and H-bonding. These derivative molecular clusters then serve as input for subsequent energy minimization. 

It views the covalent bonding and H-bonding network as a directed graph in the direction of H-bonding donation; that is 

 \[O-H> - - \[O-H> - - \[O-H> - - 

is viewed as a graph with three nodes in the rightward direction.  The result of a reversal of the graph above would look like

 O] - - <H-O] - - <H-O] - - <H 

The input is a cartesian coordinate file and the outputs are cartesian coordinate files as well.   
                               
**pymocyler.py** is a pymol helper script to facilitate visualization of H-bonding cycles.
                               
**hydrocycler_engine.py** is an automated script to generate .xyz files and a shell script to execute Gaussian .com files corresponding to these generated .xyz files. (To create .xyz files into .com files quickly, see **gt-xyz2com.py** in https://github.com/mihali/gt-x.)
 
                               
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

## Method

H-bonds can be determined from the cartesian coordinates of a cluster by measuring O-O distances between nearest neighbors and then taking H-O-O angles. A quick algorithm to perform a nearest neighbor search is by KDtree. After H-bonds have been identified in the cluster, the information can be defined as a "trio" consisting of the proton-donor oxygen, the acceptor oxygen, and the hydrogen. These can be imagined as forming a triangle. By pre-computing all the hydrogen positions where the donor and acceptor roles are swapped, reversals of cycles can easily be performed. 
                               
<img width="400" alt="Screen Shot 2022-03-31 at 10 59 07 AM" src="https://user-images.githubusercontent.com/10854556/161396405-5fb370f5-307f-4430-8e63-5c5a25d56633.png">
                               
## References

1. Uses johnson.py from https://github.com/qpwo/python-simple-cycles
2. Algorithm employed is from Johnson (1975) https://doi.org/10.1137/0204007


  

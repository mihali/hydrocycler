# hydrocycler

Usage: hydrocycler.py file.xyz

A program to find cycles of hydrogen bonding within a molecular cluster and to generate derivative molecular clusters by reversing the direction of the molecular- and H-bonding. It views the covalent bonding and H-bonding network as a directed graph in the direction of H-bonding donation; that is 

 [O-H> - - [O-H> - - [O-H> - -

is viewed as a graph with three nodes in the left to right direction. The input is a cartesian coordinate file and the output are cartesian coordinate files as well. The result of a reversal of the above graph would look like

 O] - - <H-O] - - <H-O] - - <H

## Installation

Installation procedure 
1. git clone https://github.com/mihali/hydrocycler.git 
2. cd hydrocycler
3. conda create --name hydrocycler python=3.7 
4. conda activate hydrocycler
5. conda install numpy scipy

## Test

There are files to try in the testfiles directory, e.g.

% python hydrocycler.py testfiles/si-4oh-27h2o.xyz

## Reference

1. Uses johnson.py from https://github.com/qpwo/python-simple-cycles



  

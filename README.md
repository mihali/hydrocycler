# hydrocycler

A program to find cycles of hydrogen bonding within a molecular cluster and to generate new molecular clusters by reversing the direction of the molecular bonding. It views the covalent bonding and H-bonding network as a directed graph in the direction of H-bonding, that is 
[O-H> -- [O-H> -- is viewed as a graph in the left to right direction.

The input is a cartesian coordinate file and the output are cartesian coordinate files as well. 

Usage: hydrocycler.py file.xyz

Needs johnson.py from https://github.com/qpwo/python-simple-cycles
Copy of johnson.py may be downloaded from the current repository. 

Proven installation procedure (Mac with conda)
1. Download hydrocycler.py
2. Make a directory (mkdir hydrocycler)
3. Move downloaded file to newly created hydrocycler directory. (mv \<srcpath\>/hydrocycler.py \<srcpath\>/johnson.py hydrocycler/ )
4. Create conda environment. (conda create --name hydrocycler python=3.7)
5. Activate environment. (conda activate hydrocycler)
6. Install numpy. (conda install numpy)
7. Install scipy. (conda install scipy)





  

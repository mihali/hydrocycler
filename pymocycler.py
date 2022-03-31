import sys, os
from pymol import cmd
import copy
from datetime import datetime
import numpy as np
from scipy.spatial import KDTree as kd
from collections import defaultdict

ts = datetime.now().strftime("%Y%m%d%H%m%s")
sys.setrecursionlimit(10000)

xyzfile = '/Users/mihali/Programming/hydrocycler/water_clusters_a_001_to_050/012.xyz'
#xyzfile = '/Users/mihali/Downloads/8h2oa.xyz'
#xyzfile = '/Users/mihali/Downloads/si-4oh-27h2o.xyz'


#==============

# cutoffs
oobondlim = 3.2
hobondlim = 1.2
hooangle = 0.523599    # 30 degrees
# hooangle = 0.785398    #  45 degrees


# moddir='/Users/mihali/opt/anaconda3/envs/python27/share/pymol/modules'
# sys.path.insert(0, moddir)
# os.environ['PYMOL_PATH'] = os.path.join(moddir, 'pymol/pymol_path')

# print ("XYX file to load:", end="")
# xyzfile = input()


def simple_cycles(G):
    # Yield every elementary cycle in python graph G exactly once
    # Expects a dictionary mapping from vertices to iterables of vertices
    def _unblock(thisnode, blocked, B):
        stack = set([thisnode])
        while stack:
            node = stack.pop()
            if node in blocked:
                blocked.remove(node)
                stack.update(B[node])
                B[node].clear()
    G = {v: set(nbrs) for (v,nbrs) in G.items()} # make a copy of the graph
    sccs = strongly_connected_components(G)
    while sccs:
        scc = sccs.pop()
        startnode = scc.pop()
        path=[startnode]
        blocked = set()
        closed = set()
        blocked.add(startnode)
        B = defaultdict(set)
        stack = [ (startnode,list(G[startnode])) ]
        while stack:
            thisnode, nbrs = stack[-1]
            if nbrs:
                nextnode = nbrs.pop()
                if nextnode == startnode:
                    yield path[:]
                    closed.update(path)
                elif nextnode not in blocked:
                    path.append(nextnode)
                    stack.append( (nextnode,list(G[nextnode]))  )
                    closed.discard(nextnode)
                    blocked.add(nextnode)
                    continue
            if not nbrs:
                if thisnode in closed:
                    _unblock(thisnode,blocked,B)
                else:
                    for nbr in G[thisnode]:
                        if thisnode not in B[nbr]:
                            B[nbr].add(thisnode)
                stack.pop()
                path.pop()
        remove_node(G, startnode)
        H = subgraph(G, set(scc))
        sccs.extend(strongly_connected_components(H))

def strongly_connected_components(graph):
    # Tarjan's algorithm for finding SCC's
    # Robert Tarjan. "Depth-first search and linear graph algorithms." SIAM journal on computing. 1972.
    # Code by Dries Verdegem, November 2012
    # Downloaded from http://www.logarithmic.net/pfh/blog/01208083168

    index_counter = [0]
    stack = []
    lowlink = {}
    index = {}
    result = []
    
    def _strong_connect(node):
        index[node] = index_counter[0]
        lowlink[node] = index_counter[0]
        index_counter[0] += 1
        stack.append(node)
    
        successors = graph[node]
        for successor in successors:
            if successor not in index:
                _strong_connect(successor)
                lowlink[node] = min(lowlink[node],lowlink[successor])
            elif successor in stack:
                lowlink[node] = min(lowlink[node],index[successor])

        if lowlink[node] == index[node]:
            connected_component = []

            while True:
                successor = stack.pop()
                connected_component.append(successor)
                if successor == node: break
            result.append(connected_component[:])
    
    for node in graph:
        if node not in index:
            _strong_connect(node)
    
    return result

def remove_node(G, target):
    # Completely remove a node from the graph
    # Expects values of G to be sets
    del G[target]
    for nbrs in G.values():
        nbrs.discard(target)

def subgraph(G, vertices):
    # Get the subgraph of G induced by set vertices
    # Expects values of G to be sets
   return {v: G[v] & vertices for v in vertices}

#=============
def createarrays (inp):

    hinput = []
    oinput = []
    xinput = []
    xyzdict = {}
    ocount = 0
    hcount = 0
    xcount = 0
    linenum=0
    
    for line in inp:
        toks = line.split()
        count = 0
        if len(toks) >= 3:
          if toks[0] == 'H':
            hinput.append ([float(toks[1]), float(toks[2]), float(toks[3])])
            count = hcount
            hcount = hcount + 1
          elif toks[0] == 'O':
            oinput.append([float(toks[1]), float(toks[2]), float(toks[3])])
            count = ocount 
            ocount = ocount + 1
          else:
            xinput.append([float(toks[1]), float(toks[2]), float(toks[3])])     
            count = xcount
            xcount = xcount + 1        
          xyzdict [(toks[0], count)] = [linenum, float(toks[1]), float(toks[2]), float(toks[3])]  
        linenum = linenum +1

    ocoords = np.array(oinput)
    hcoords = np.array(hinput)
    xcoords = np.array(xinput)

    return [ ocoords, hcoords, xcoords, xyzdict ]

#==============
def findcycles (ocoords,  hcoords):
#def findcycles ():

    okdt = kd(ocoords)                    # kdtree to find neighbors
    hkdt = kd(hcoords)

    count = 0
    ograph = {}
    trio = {}
    for oatom in ocoords:                 # go through each O coordinate
      oidx = okdt.query_ball_point(oatom, r=oobondlim )  # an O atom has neighbors
      for neighoidx in oidx:               # oxygen pair pair to test H-bond connection
        dist = np.linalg.norm(oatom - ocoords[neighoidx])
        if dist > 0.01:                     # exclude itself from generating the trio 
          hidx = hkdt.query_ball_point(oatom, r=hobondlim)
          for neighhidx in hidx:            
            hbondtest = isanhbond(hcoords[neighhidx], oatom, ocoords[neighoidx]) # isanhbond returns [flag, theta]
            if hbondtest[0] == 1:           # found an H-bond 
                trio [(count,neighoidx)] = [neighhidx, hbondtest[1]]
                if count in ograph:
                  ograph[count].append(neighoidx)
                else:
                  ograph[count] = [neighoidx]
      if count not in ograph:                   # terminal OH need keys
        ograph[count] = []
      count = count + 1
    return [ograph, trio]

#==============
def isanhbond(hatom, oatom1, oatom2):
    h=np.array (hatom)
    o1=np.array (oatom1)
    o2=np.array (oatom2)
      
    a = h - o1
    ia = a/np.linalg.norm(a)
    b = o2 - o1
    ib = b/np.linalg.norm(b)
    theta = np.arccos(np.dot (ia, ib))
    if theta < hooangle:
        return [1, theta]
    else:
        return [0, 0.0]     

#==============
def displaycycles (ograph, xyzdict, trio):

    cmd.do ( "set dash_gap, 0.5")
    cmd.do ( "set dash_radius, 0.05")
    for cycle in tuple(simple_cycles(ograph)):
      previous = -1
      for x in cycle:
          if (previous, x) in trio:
            numbero = xyzdict[('O', previous)][0] + 1 
            cmd.do ( "show spheres, resi %s"%numbero) 
            numbero2 = xyzdict[('O', x)][0] + 1
            numberh = xyzdict[('H',trio[(previous, x)][0])][0] + 1
            cmd.do ( "show spheres, resi %s"%numberh) 
            cmd.do ( "distance d, resi %s, resi %s"%(numberh,numbero2) )
            cmd.do ( "hide labels, d")
          previous = x      
      
      numbero = xyzdict[('O', previous)][0] + 1
      cmd.do ( "show spheres, resi %s"%numbero) 
      numbero2 = xyzdict[('O', cycle[0])][0] + 1
      numberh = xyzdict[('H',trio[(previous, cycle[0])][0])][0] + 1
      cmd.do ( "show spheres, resi %s"%numberh) 
      cmd.do ( "distance d, resi %s, resi %s"%(numberh,numbero2) )
      cmd.do ( "hide labels, d")

      return      

#==============
def fn ( inp ):
 
    [ocoords, hcoords, xcoords, xyzdict ] = createarrays (inp) 
    [ograph, trio] = findcycles (ocoords,  hcoords)
    
    numcycles = displaycycles(ograph, xyzdict, trio)

#=============
# BEGIN

cmd.load(xyzfile)
string = open ( xyzfile, 'r' ).readlines()[2:]
fn (string)

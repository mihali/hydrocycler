#!/usr/bin/env  python  
# This version takes one file input

import sys, os
import copy, time
from datetime import datetime
import numpy as np
from scipy.spatial import KDTree as kd
import johnson
from random import randint

ts = datetime.now().strftime("%Y%m%d%H%m%s")
sys.setrecursionlimit(10000)

#xyzfile = '/Users/mihali/Programming/hydrocycler/water_clusters_a_001_to_050/012.xyz'
#xyzfile = '/Users/mihali/Downloads/8h2oa.xyz'
xyzfile = '/Users/mihali/Downloads/si-4oh-27h2o.xyz'

#==============

# cutoffs
oobondlim = 3.2
hobondlim = 1.2
hooangle = 0.523599    # 30 degrees
# hooangle = 0.785398    #  45 degrees

#==============
# Log all output 

class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("hydrocycler-%s.out"%ts, "a")
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message) 
    def flush(self): 
        self.terminal.flush(); 
        self.log.flush()
sys.stdout = Logger()

#==============
def exportcartesian (xyzdict, filename):

  nts = datetime.now().strftime("%Y%m%d%H%m%s")
  filen = str(os.path.splitext(filename)[0]).split("/")[-1] 
  fd = open("%s-%s.xyz"%(filen,nts), "w")
  outputs = [sys.stdout, fd ]

  for output in outputs: 
   print(len(xyzdict.keys()),file=output)
   print("This was generated by Hydrocycler (Felipe, 2022)", file=output)
   for key in xyzdict:
     print(key[0], *(x for x in xyzdict[key]), file=output)
   print('\n\n', file=output)
  
  fd.flush()

#==============
def exportbatchjob (filename):
  global batchjob
  filen = str(os.path.splitext(filename)[0]).split("/")[-1] 
  print ("g16 %s-%s.com & \nwait"%(filen,ts), file=batchjob)

#==============
def isanhbond(h, o1, o2):
 
  a = h - o1
  ia = a/np.linalg.norm(a)
  b = o2 - o1
  ib = b/np.linalg.norm(b)
  theta = np.arccos(np.dot (ia, ib))
  if theta < hooangle:
    return [1, theta]
  else:
    return [0, 0.0]     

#=============
def createarrays (inp):

  hinput = []
  oinput = []
  xinput = []
  xyzdict = {}
  ocount = 0
  hcount = 0
  xcount = 0
  
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
      xyzdict [(toks[0], count)] = [float(toks[1]), float(toks[2]), float(toks[3])]  

  ocoords = np.array(oinput)
  hcoords = np.array(hinput)
  xcoords = np.array(xinput)

  return [ ocoords, hcoords, xcoords, xyzdict ]

#==============
def reverseit (o1, o2, h, theta):
  dist = np.absolute(o2-o1)
  unitv = (o2-o1)/dist
  displacement = dist - 2*np.absolute(h-o1)*np.cos(theta)
  return h + unitv * displacement

#==============
def printhistory(history):

    print ("Cycle Reversal History:")
    for x in history:
        print ([z+1 for z in x])  

#==============
def fn ( inp, file ):
 
 #--------------
  def modifycycle(choice): 
      [ograph, trio] = findcycles (ocoords,  hcoords)
      cycle = tuple(johnson.simple_cycles(ograph))[choice - 1]
      history.append(cycle)                                     # note it down
      print(cycle)
      previous = -1
      hdict = {}
      for node in cycle:
            if (previous, node) in trio:
              o1 = np.array(ocoords[previous])
              o2 = np.array(ocoords[node])
              h  = np.array(hcoords[trio[(previous, node)][0]])
              hdict[('H',trio[(previous, node)][0])] = h
              theta = trio[(previous, node)][1]
              newh = reverseit (o1, o2, h, theta) 
              hcoords[trio[(previous, node)][0]] = newh                  
              xyzdict[('H',trio[(previous, node)][0])] = newh          
            previous = node
      o1 = np.array(ocoords[previous])
      o2 = np.array(ocoords[cycle[0]])
      h  = np.array(hcoords[trio[(previous, cycle[0])][0]])
      hdict[('H',trio[(previous, cycle[0])][0])] = h
      theta = trio[(previous, cycle[0])][1]
      newh = reverseit (o1, o2, h, theta) 
      hcoords[trio[(previous, cycle[0])][0]] = newh
      xyzdict[('H',trio[(previous, cycle[0])][0])] = newh
      return hdict

 #--------------
  def restorecycle(hdict): 
      for (a,b) in hdict:
          if (a,b) in xyzdict:
              hcoords[b] = hdict[(a,b)]
              xyzdict[(a,b)] = hdict[(a,b)]          

  #--------------
  def sideways (ocoords, hcoords, xcoords):

      [ograph, trio] = findcycles (ocoords, hcoords)  
      numcycles = displaycycles(ograph)       # this is where we begin
      i = 0
      while i < numcycles:
        restoredict = modifycycle(i)
        exportcartesian (xyzdict, file)
        time.sleep(1)
        exportbatchjob (file)
        restorecycle (restoredict)
        [ograph, trio] = findcycles (ocoords,  hcoords)    
        i = i + 1

  #--------------
  def traverse (depth):
 
      if depth == 0:
        return
      [ograph, trio] = findcycles (ocoords,  hcoords)        
      numcycles = displaycycles(ograph)       # this is where we begin
      choicei = randint (1, numcycles)
      modifycycle(choicei)  
      exportcartesian (xyzdict, file)
      time.sleep(1)
      exportbatchjob (file)
      traverse (depth - 1)

  #--------------
  def displaycycles (ograph):

    print("\n")
    count = 1
    for cycle in tuple(johnson.simple_cycles(ograph)):
      print (str(count) + ": ", end='')
      print ([x+1 for x in cycle])        # incremented to coincide with residue number
      count = count + 1
    return count - 1

  #--------------
  def findcycles (ocoords,  hcoords):

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
  
  [ocoords, hcoords, xcoords, xyzdict ] = createarrays (inp) 
  [ograph, trio] = findcycles (ocoords,  hcoords)
  history = []

#  traverse(2)
  sideways(ocoords, hcoords, xcoords)

#=============
print ("\nWelcome to Hydrocycler (c) by Mihali Felipe (2022)\n")
print (" Calculating directed graphs with |O-H> nodes and H-bond edges.")
print (" H-bond parameters:")
print ("   H--O distance: %sA-%sA"%(hobondlim, oobondlim))
print ("   H-O--O angle:  <%s degrees "%int((hooangle*180/3.14159)))

batchjob = open("batchjob-%s.sh"%ts, "a") 
string = open ( xyzfile, 'r' ).readlines()[2:]

fn (string, xyzfile)

 #####################
 # Traversal plan
 # 1. We go through each level 1 config
 # 2. We go through all combinations (nuts)
 # 2. We modify random amounts in, random amounts across
 # each config made would go to bash script and have own com file
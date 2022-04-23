#!/usr/bin/env  python  
# This version takes one file input

import sys, os
import copy, time
from datetime import datetime
import numpy as np
from scipy.spatial import KDTree as kd
import johnson
from random import randint

ts = datetime.now().strftime("%y%m%d%H%M%S%f")
sys.setrecursionlimit(10000)

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
        self.log = open("hydrocycler_counter-%s.csv"%ts, "a")
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message) 
    def flush(self): 
        self.terminal.flush(); 
        self.log.flush()

sys.stdout = Logger()

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


#==============
def fn ( files ):  # filename

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

  for file in files:    
    string = open ( file, 'r' ).readlines()[2:] 
    [ocoords, hcoords, xcoords, xyzdictmaster ] = createarrays (string) 
    [ograph, trio] = findcycles (ocoords, hcoords)
    cycles = sorted(tuple(johnson.simple_cycles(ograph)))
    cycles_grp = tuple([ x for x in cycles ])
    cycles_max = max(map(len,cycles))
    cycles_sig = tuple([ x for sublist in cycles for x in sublist ])
    cycles_count = str(len(cycles_sig))
    print ("%s,\t%s,\t%s,\t%s"%(cycles_count, cycles_max, str(len(tuple([ x for x in cycles ]))),file))
#=============
print (" Calculating directed graphs with |O-H> nodes and H-bond edges.")
print (" H-bond parameters:")
print ("   H--O distance: %sA-%sA"%(hobondlim, oobondlim))
print ("   H-O--O angle:  <%s degrees \n"%int((hooangle*180/3.14159)))

argv = sys.argv[1:]
argc  = len(sys.argv)
command = sys.argv[0].split('/')[len(sys.argv[0].split('/')) - 1] 
if not argv:
   print("Usage: %s [file]" % command )
   print ("No stdin option for this command. Use file argument")
else:
  counter = 0
  for file in argv:
    fn ( [file] )
print ("Sig ln\tMax ln\tCyc ct\File\n\n")
print ("       Please see hydrocycler_counter-%s.csv"%ts)


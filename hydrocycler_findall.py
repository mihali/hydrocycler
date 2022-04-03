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
        self.log = open("hydrocycler_engine-%s.out"%ts, "a")
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message) 
    def flush(self): 
        self.terminal.flush(); 
        self.log.flush()

sys.stdout = Logger()

#==============
def exportcartesian (xyzdict, fd):

  outputs = [sys.stdout, fd ]

  for output in outputs: 
   print(len(xyzdict.keys()),file=output)
   print("This was generated by Hydrocycler (Felipe, 2022)", file=output)
   for key in xyzdict:
     print(key[0], *(x for x in xyzdict[key]), file=output)
   print('\n\n', file=output)

#==============
def exportbatchjob (filename):
  global batchjob
#  filen = str(os.path.splitext(filename)[0]).split("/")[-1] 
#  print ("g16 %s-%s.com & \nwait"%(filen,nts), file=batchjob)
  print ("g16 %s.com & \nwait"%(filename), file=batchjob)

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
def reverseit (o1, o2, h, theta):

  dist = np.linalg.norm(o2-o1)
  unitv = (o2-o1)/dist
  displacement = dist - 2*np.linalg.norm(h-o1)*np.cos(theta)
  return h + unitv * displacement

#==============
def fn ( files ):  # list of filenames

 #--------------
  def modifycycle(cycle, filename): 

      history.append(cycle)                                     # note it down
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
              hcoords[trio[(previous, node)][0]] = np.ndarray.tolist(newh)                  
              xyzdict[('H',trio[(previous, node)][0])] = np.ndarray.tolist(newh)          
            previous = node
      o1 = np.array(ocoords[previous])
      o2 = np.array(ocoords[cycle[0]])
      h  = np.array(hcoords[trio[(previous, cycle[0])][0]])
      hdict[('H',trio[(previous, cycle[0])][0])] = h
      theta = trio[(previous, cycle[0])][1]
      newh = reverseit (o1, o2, h, theta) 
      hcoords[trio[(previous, cycle[0])][0]] = np.ndarray.tolist(newh)
      xyzdict[('H',trio[(previous, cycle[0])][0])] = np.ndarray.tolist(newh)

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
  
  #-------------
  def recreatearrays (xyzdict):

    hinput = []
    oinput = []
    xinput = []
  
    for (toks, count) in xyzdict:
        if toks == 'H':
          hinput.append (xyzdict[toks, count])
        elif toks[0] == 'O':
          oinput.append (xyzdict[toks, count])
        else:
          xinput.append (xyzdict[toks, count]) 
  
    ocoords = np.array(oinput)
    hcoords = np.array(hinput)
    xcoords = np.array(xinput)

    return [ ocoords, hcoords, xcoords ]

  #--------------
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

# ------- fn starts here -----------#

  fileslocal = files.copy()
  fileschild =[]

  for file in fileslocal:                                     # open up current file
    string = open ( file, 'r' ).readlines()[2:] 
    [ocoords, hcoords, xcoords, xyzdictmaster ] = createarrays (string) 
    [ograph, trio] = findcycles (ocoords, hcoords)
    cycles = sorted(tuple(johnson.simple_cycles(ograph)))            # determine signature
    cycles_sig = tuple([ x for sublist in cycles for x in sublist ]) # if we've seen cycles before we don't
                                                                     # follow the cycle

    if (cycles_sig) not in cyclesdict:

      cyclesdict[cycles_sig]="DONE"                             # we've gone through this file

      for cycle in cycles:
        xyzdict = copy.deepcopy(xyzdictmaster)                  # each cyle corresponds to a new file
        modifycycle(cycle, xyzdict)  
        [ocoords2, hcoords2, xcoords2] = recreatearrays (xyzdict)
        [ograph2, trio2] = findcycles (ocoords2, hcoords2)
        cycles2 = sorted(tuple(johnson.simple_cycles(ograph2)))
        cycles_sig2 = tuple([ x for sublist in cycles2 for x in sublist ])
        if (cycles_sig2) not in cyclesdict:
          nts = datetime.now().strftime("%y%m%d%H%M%S%f")         # generate a name
          filen = str(os.path.splitext(filename)[0]).split("/")[-1] 
          fd = open("%s-%s.xyz"%(filen,nts), "w")
          fileschild.append("%s-%s.xyz"%(filen,nts))              # pass these to fn 
          exportcartesian(xyzdict, fd)
          exportbatchjob("%s-%s.xyz"%(filen,nts))
          fd.flush()
    else:
      pass
    
    if len(fileschild) > 0:
        print (fileschild)
        fn (fileschild) 
  


#=============
print (" Calculating directed graphs with |O-H> nodes and H-bond edges.")
print (" H-bond parameters:")
print ("   H--O distance: %sA-%sA"%(hobondlim, oobondlim))
print ("   H-O--O angle:  <%s degrees \n"%int((hooangle*180/3.14159)))

batchjob = open("batchjob-%s.sh"%ts, "a") 

files = []
history = []
cyclesdict = {}

argv = sys.argv[1:]
argc  = len(sys.argv)
command = sys.argv[0].split('/')[len(sys.argv[0].split('/')) - 1] 
if not argv:
   print("Usage: %s [file]" % command )
   print ("No stdin option for this command. Use file argument")
elif argc == 2:
   for filename in argv:
     files.append(filename) 
     fn ( files )
     
else:
   print("Usage: %s [file]" % command )


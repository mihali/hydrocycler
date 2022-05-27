#!/usr/bin/env  python  

import sys, os
import copy, time
from datetime import datetime
import numpy as np
from scipy.spatial import KDTree as kd
import johnson
import pickle

ts = datetime.now().strftime("%y%m%d%H%M%S%f")

#==============
# cutoffs
oobondlim = 3.2
hobondlim = 1.2
hooangle = 0.523599    # 30 degrees
# hooangle = 0.785398    #  45 degrees


#==============
def reverseit (o1, o2, h, theta):

    dist = np.linalg.norm(o2-o1)
    unitv = (o2-o1)/dist
    displacement = dist - 2*np.linalg.norm(h-o1)*np.cos(theta)
    return h + unitv * displacement

#==============
def exportcartesian (xyzdict, fd):

  outputs = [sys.stdout, fd ]
  for output in outputs: 
    print(len(xyzdict.keys()),file=output)
    print("This was generated by Hydrocycler (Felipe, 2022)", file=output)
    for key in xyzdict:
      print(key[0], *(x for x in xyzdict[key]), file=output)
    print('\n\n', file=output)
  fd.flush()

#==============
def displaycycles (ograph):

  print("\n")
  count = 1
  for cycle in tuple(johnson.simple_cycles(ograph)):
    print (str(count) + ": ", end='')
    print ([x+1 for x in cycle])        # incremented to coincide with residue number
    count = count + 1
  return count - 1


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
def gettrio (inp):

    [ ocoords, hcoords, xcoords, xyzdict1 ] = createarrays (inp)
    [ograph, trio1] = findcycles (ocoords,  hcoords)
    return trio1

#==============
def printtrio ():

    print ("Processing with trios:")
    for key in trio:
        print ("\t%s\t%s"%(str(key),str(trio[key])) , )

#===============
def createarrays (inp):

    hinput = []
    oinput = []
    xinput = []
    xyzdict1 = {}
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
            xyzdict1 [(toks[0], count)] = np.array([ float(toks[1]), float(toks[2]), float(toks[3]) ])  

    ocoords = np.array(oinput)
    hcoords = np.array(hinput)
    xcoords = np.array(xinput)

    return [ ocoords, hcoords, xcoords, xyzdict1 ]

#=============
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

#=============
def findcycles (ocoords,  hcoords):

    okdt = kd(ocoords)                    # kdtree to find neighbors
    hkdt = kd(hcoords)

    count = 0
    ograph = {}
    trio1 = {}

    for oatom in ocoords:                 # go through each O coordinate
        oidx = okdt.query_ball_point(oatom, r=oobondlim )  # an O atom has neighbors
        for neighoidx in oidx:               # oxygen pair pair to test H-bond connection
            dist = np.linalg.norm(oatom - ocoords[neighoidx])
            if dist > 0.01:                     # exclude itself from generating the trio1 
                hidx = hkdt.query_ball_point(oatom, r=hobondlim)
                for neighhidx in hidx:    
                    hbondtest = isanhbond(hcoords[neighhidx], oatom, ocoords[neighoidx]) 
                    if hbondtest[0] == 1:                                     # found an H-bond 
                        trio1 [(count,neighoidx)] = [neighhidx, hbondtest[1], hcoords[neighhidx]] 
                        o2 = ocoords[neighoidx]
                        o1 = ocoords[count]
                        h  = hcoords[neighhidx]
                        newh = reverseit (o1, o2, h, hbondtest[1])
                        trio1 [(neighoidx, count)] = [neighhidx, hbondtest[1], newh]
                        if count in ograph:
                            ograph[count].append(neighoidx)
                        else:
                            ograph[count] = [neighoidx]
        if count not in ograph:                   # terminal OH need keys
            ograph[count] = []
        count = count + 1
    return [ograph, trio1]



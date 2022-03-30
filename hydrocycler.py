#!/usr/bin/env  python  
# This version takes one file input

import sys, os
import copy
from datetime import datetime
import numpy as np
from scipy.spatial import KDTree as kd
import johnson

ts = datetime.now().strftime("%Y%m%d%H%m%s")

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
def displaycycles (ograph):
 print("\n")
 count = 1
 for cycle in tuple(johnson.simple_cycles(ograph)):
   print (str(count) + ": ", end='')
   print ([x+1 for x in cycle])        # incremented to coincide with residue number
   count = count + 1
 return count - 1

#==============
def findcycles (ocoords,  hcoords):
#def findcycles ():

 okdt = kd(ocoords)                    # kdtree to find neighbors
 hkdt = kd(hcoords)

 hootriplet = []
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

#=============
def createarrays (inp):

 hinput = []
 oinput = []
 xinput = []
 xyzdict = {}
 ocount = 0
 hcount = 0
 xcount = 0
 odict ={}
 
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

def print_menu(): 
    print ("\n============ H Y D R O C Y C L E R =================")      
    print (" 1. Show H-bond cycles")
    print (" 2. Show cycles and choose cycle for reversal")
    print (" 3. Save configuration and continue making changes")
    print (" 4. Save configuration and start from the beginning")
    print (" 5. Save configuration and exit")
    print (" 6. Show cycle reversal history for current configuration")
    print (" 7. Exit")
    print ("====================================================\n")      
    print ("Enter your choice [1-7]: ", end='') 


#==============
def fn ( inp, file ):
 
 [ocoords, hcoords, xcoords, xyzdict ] = createarrays (inp) 
 history = []
 ocoords_orig = copy.deepcopy(ocoords)
 hcoords_orig = copy.deepcopy(hcoords)
 xyzdict_orig = copy.deepcopy(xyzdict)

 [ograph, trio] = findcycles (ocoords,  hcoords)
 
 def modifycycles(numcycles): 
  print ("\nEnter a cycle to reverse. Enter nothing to get back to menu: ", end='')
  loop=True       
  while loop:
    choice = input()
    choicei = int(choice)
    if choice=="":
      loop=False
      return
    elif choice.isnumeric() and choicei > 0 and choicei <= numcycles:
      cycle = tuple(johnson.simple_cycles(ograph))[choicei - 1]
      history.append(cycle)                                     # note it down
      previous = -1
      for node in cycle:
        if (previous, node) in trio:
          o1 = ocoords[previous]
          o2 = ocoords[node]
          h  = hcoords[trio[(previous, node)][0]]
          theta = trio[(previous, node)][1]
          newh = reverseit (o1, o2, h, theta) 
          hcoords[trio[(previous, node)][0]] = newh                  
          xyzdict[('H',trio[(previous, node)][0])] = newh          
        previous = node
      o1 = ocoords[previous]
      o2 = ocoords[cycle[0]]
      h  = hcoords[trio[(previous, cycle[0])][0]]
      theta = trio[(previous, cycle[0])][1]
      newh = reverseit (o1, o2, h, theta) 

      hcoords[trio[(previous, cycle[0])][0]] = newh
      xyzdict[('H',trio[(previous, cycle[0])][0])] = newh

      loop=False
      return
    else:
      print ("Please enter a cycle to reverse. Enter nothing to get back to menu: ", end='')

 ######################
 # Ye Old-School menu #
 # ####################
 loop=True       
 while loop:          
    print_menu()    
    choice = input()

    if choice =='1':
        numcycles = displaycycles(ograph)
    elif choice=='2':   #"Show cycles (Johnson, 1975) and choose cycle for reversal")
        numcycles = displaycycles(ograph)
        modifycycles(numcycles)   # Change H positions
        [ograph, trio] = findcycles (ocoords,  hcoords)
        numcycles = displaycycles(ograph)
    elif choice=='3': #"Save configuration and continue making changes"
        exportcartesian (xyzdict, file)
    elif choice=='4': #"Save configuration and start from the beginning"
        exportcartesian (xyzdict, file)
        history = []
        ocoords = ocoords_orig
        hcoords = hcoords_orig
        xyzdict = xyzdict_orig
        [ograph, trio] = findcycles (ocoords,  hcoords)
    elif choice=='5': #"Save configuration and exit"
        exportcartesian (xyzdict, file)
        print("Session log file: hydrocycler-%s.out"%ts)
        print("Thank you for using Hydrocycler!\n")
        loop=False
    elif choice=='6':
        printhistory(history)
    elif choice=='7':
        print("Session log file: hydrocycler-%s.out"%ts)
        print("Thank you for using Hydrocycler!\n")
        loop=False 
    else:
        print("Enter any key between 1-7")

#=============
print ("\nWelcome to Hydrocycler (c) by Mihali Felipe (2022)\n")
print (" Calculating directed graphs with |O-H> nodes and H-bond edges.")
print (" H-bond parameters:")
print ("   H--O distance: %sA-%sA"%(hobondlim, oobondlim))
print ("   H-O--O angle:  <%s degrees "%int((hooangle*180/3.14159)))

files = sys.argv[1:]
argc  = len(sys.argv)
command = sys.argv[0].split('/')[len(sys.argv[0].split('/')) - 1] 

if not files:
   print("Usage: %s [file]" % command )
   print ("No stdin option for this command. Use file argument")
elif argc == 2:
   for file in files:
     string = open ( file, 'r' ).readlines()[2:]
     fn ( string, file )
else:
   print("Usage: %s [file]" % command )


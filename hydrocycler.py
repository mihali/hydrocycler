#!/usr/bin/env  python  
# This version takes one file input

import sys, os
import copy
from datetime import datetime
import numpy as np
from scipy.spatial import KDTree as kd
import johnson
from hydrocycler_utils import *
 
ts = datetime.now().strftime("%y%m%d%H%M%S%f")
history = []

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
def printhistory(history):

    print ("\nCycle Reversal History. Null ([]) means reset.:")
    for x in history:
        print ([z+1 for z in x])  

#==============
def print_menu(): 

    print ("\n============ H Y D R O C Y C L E R =================")      
    print (" 1. Show H-bond cycles")
    print (" 2. Show cycles and choose cycle for reversal")
    print (" 3. Save configuration and continue making changes")
    print (" 4. Save configuration and start from the beginning")
    print (" 5. Save configuration and exit")
    print (" 6. Show cycle reversal history for this session")
    print (" 7. Exit")
    print ("====================================================\n")      
    print ("Enter your choice [1-7]: ", end='') 

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
def modifycycles(xyzdict): 

  [ocoords, hcoords, xcoords] = recreatearrays (xyzdict)
  [ograph, trio2] = findcycles (ocoords, hcoords)  
  numcycles = displaycycles(ograph)
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
              newh = trio[(node, previous)][2]                  
              xyzdict[('H',trio[(previous, node)][0])] = newh          
          previous = node
      newh = trio[(cycle[0], previous)][2] 
      xyzdict[('H',trio[(previous, cycle[0])][0])] = newh
      loop=False
      return xyzdict
    else:
      print ("Please enter a cycle to reverse. Enter nothing to get back to menu: ", end='')


#==============
def exportxyz(xyzdict, filename): 

                    nts = datetime.now().strftime("%y%m%d%H%M%S%f")         
                    filen = str(os.path.splitext(filename)[0]).split("/")[-1] 
                    fd = open("%s-%s.xyz"%(filen,nts), "w")
                    exportcartesian(xyzdict, fd)

#==============
def fn ( inp, file ):
 
  global trio
  [ocoords, hcoords, xcoords, xyzdictmaster ] = createarrays (inp) 
  [ograph, trio] = findcycles (ocoords,  hcoords)
  xyzdict = copy.deepcopy(xyzdictmaster)                  

  # Ye Old-School menu #
  loop=True       
  while loop:          
      print_menu()    
      choice = input()

      if choice =='1':
          numcycles = displaycycles(ograph)
      elif choice=='2':   #"Show cycles (Johnson, 1975) and choose cycle for reversal")
          xyzdict = modifycycles(xyzdict)   
          [ocoords, hcoords, xcoords] = recreatearrays (xyzdict)
          [ograph, trio1] = findcycles (ocoords, hcoords)
          numcycles = displaycycles(ograph)
      elif choice=='3': #"Save configuration and continue making changes"
          exportxyz (xyzdict, file)
      elif choice=='4': #"Save configuration and start from the beginning"
          exportxyz (xyzdict, file)
          xyzdict = copy.deepcopy(xyzdictmaster)
          [ocoords, hcoords, xcoords] = recreatearrays (xyzdict)
          [ograph, trio1] = findcycles (ocoords, hcoords)
          history.append([])
      elif choice=='5': #"Save configuration and exit"
          exportxyz (xyzdict, file)
          print("Session log file: hydrocycler-%s.out"%ts)
          loop=False
      elif choice=='6':
          printhistory(history)
      elif choice=='7':
          print("Session log file: hydrocycler-%s.out"%ts)
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
        trio = gettrio (string)
        fn ( string, file )
        print ("\n\nThank you for using Hydrocycler!\n\n")
else:
   print("Usage: %s [file]" % command )


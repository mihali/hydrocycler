#!/usr/bin/env  python  

import sys, os
import copy, time
from datetime import datetime
import numpy as np
from scipy.spatial import KDTree as kd
import johnson
import pickle
from hydrocycler_utils import *

ts = datetime.now().strftime("%y%m%d%H%M%S%f")

#==============
# Log all output 

class Logger(object):

    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("hydrocycler_findall-%s.out"%ts, "a")
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message) 
    def flush(self): 
        self.terminal.flush(); 
        self.log.flush()

sys.stdout = Logger()

#==============
def exportbatchjob (filename):

  global batchjob
  print ("g16 %s.com & \nwait"%(filename), file=batchjob)

#==============
def modifycycle(cycle, xyzdict): 

    global trio
    previous = -1
    for node in cycle:
        if (previous, node) in trio:
            newh = trio[(node, previous)][2]                  
            xyzdict[('H',trio[(previous, node)][0])] = newh          
        previous = node
    newh = trio[(cycle[0], previous)][2] 
    xyzdict[('H',trio[(previous, cycle[0])][0])] = newh
    return xyzdict

#==============
def fn_recursive ( files ):  
  
    global cyclesdict

    fileslocal = files.copy()
    fileschild =[]
    
    for file in fileslocal:                                     # open up current file
        print ("Processing %s for cycles..."%file)
        string = open ( file, 'r' ).readlines()[2:] 
        [ocoords1, hcoords1, xcoords1, xyzdictmaster] = createarrays (string) 
        [ograph1, trio1] = findcycles (ocoords1, hcoords1)
        cycles = tuple(johnson.simple_cycles(ograph1))            # determine signature using
        cycles_srt = sorted(tuple(johnson.simple_cycles(ograph1)))            # determine signature using
        cycles_sig = tuple([ x for sublist in cycles_srt for x in sublist ]) # cycles for uniqueness 
        print ("This configuration has signature: %s"%str(cycles_sig))
        if (cycles_sig) not in cyclesdict:
            cyclesdict[cycles_sig]="DONE"                             # mark this file read

            for cycle in cycles:                                      # we go through each cycle
                xyzdict = copy.deepcopy(xyzdictmaster)                  # to create a file
                xyzdict = modifycycle(cycle, xyzdict)  
                [ocoords2, hcoords2, xcoords2] = recreatearrays (xyzdict)
                [ograph2, trio2] = findcycles (ocoords2, hcoords2)
                cycles2_srt = sorted(tuple(johnson.simple_cycles(ograph2)))
                cycles_sig2 = tuple([ x for sublist2 in cycles2_srt for x in sublist2 ])

                if (cycles_sig2) not in cyclesdict:                      # export this cartesian
                    nts = datetime.now().strftime("%y%m%d%H%M%S%f")         
                    filen = str(os.path.splitext(filename)[0]).split("/")[-1] 
                    fd = open("%s-%s.xyz"%(filen,nts), "w")
                    print ("%s-%s.xyz"%(filen,nts))
                    fileschild.append("%s-%s.xyz"%(filen,nts))           # pass these to fn 
                    exportcartesian(xyzdict, fd)
                    exportbatchjob("%s-%s"%(filen,nts))
                    fd.close()
                    if len(fileschild) > 0:
                        fn (fileschild) 

#=============

class XYZCycleDat:

    def __init__(self, pfilename ):
       global cyclesdict
       global cyclesseen
       self.pfilename = pfilename                                                # start from saved file
       string = open ( pfilename, 'r' ).readlines()[2:] 
       [ocoords1, hcoords1, xcoords1, xyzdictmaster] = createarrays (string) 
       [ograph1, trio1] = findcycles (ocoords1, hcoords1)
       cycles = tuple(johnson.simple_cycles(ograph1))            
       cycles_srt = sorted(tuple(johnson.simple_cycles(ograph1))) 
       cycles_sig = tuple([ x for sublist in cycles_srt for x in sublist ])      # get its signature from cycles
       if (cycles_sig) not in cyclesdict:
            cyclesdict[cycles_sig]="DONE"                 

            for cycle in cycles:                                                 # go through each cycle
                xyzdict = copy.deepcopy(xyzdictmaster)   
                xyzdict = modifycycle(cycle, xyzdict)                            # create child from parent
                [ocoords2, hcoords2, xcoords2] = recreatearrays (xyzdict)
                [ograph2, trio2] = findcycles (ocoords2, hcoords2)
                cycles2_srt = sorted(tuple(johnson.simple_cycles(ograph2)))
                cycles_sig2 = tuple([ x for sublist2 in cycles2_srt for x in sublist2 ]) # get child signature

                if (cycles_sig2) not in cyclesdict and (cycles_sig2) not in cyclesseen:  # child is unique ...
                    cyclesseen[cycles_sig2]="DONE"                 
                    print ("This configuration has signature: %s"%str(cycles_sig2))
                    nts = datetime.now().strftime("%y%m%d%H%M%S%f")                      # ... so save child
                    filen = str(os.path.splitext(filename)[0]).split("/")[-1] 
                    fd = open("%s-%s.xyz"%(filen,nts), "w")
                    print ("%s-%s.xyz"%(filen,nts))
                    files.append("%s-%s.xyz"%(filen,nts))                                # and put in stack
                    exportcartesian(xyzdict, fd)
                    exportbatchjob("%s-%s"%(filen,nts))
                    fd.close()

#=============

def fn_iter ( files ):  
      
    while files:                                     
        file=files.pop()
        XYZCycleDat(file)

#=============
print (" Calculating directed graphs with |O-H> nodes and H-bond edges.")
print (" H-bond parameters:")
print ("   H--O distance: %sA-%sA"%(hobondlim, oobondlim))
print ("   H-O--O angle:  <%s degrees \n"%int((hooangle*180/3.14159)))

batchjob = open("batchjob-%s.sh"%ts, "a") 

files = []
cyclesdict = {}
cyclesseen = {}

argv = sys.argv[1:]
argc  = len(sys.argv)
command = sys.argv[0].split('/')[len(sys.argv[0].split('/')) - 1] 
if not argv:
   print("Usage: %s [file]" % command )
   print ("No stdin option for this command. Use file argument")
elif argc == 2:
    for filename in argv:
        files.append(filename) 
        inp = open ( filename, 'r' ).readlines()[2:] 
        trio = gettrio (inp)
        print ("Processing with trios:")
        for key in trio:
            print ("\t%s\t%s"%(str(key),str(trio[key])) , )
        fn_iter ( files )                                            # recurses on list
        print ("\n\nThank you for using Hydrocycler!\n\n")
        filen = str(os.path.splitext(filename)[0]).split("/")[-1] 
        sigdict_fd = open ("%s-%s.sig"%(filen,ts), "wb")
        pickle.dump (cyclesdict, sigdict_fd, protocol=pickle.HIGHEST_PROTOCOL )

else:
   print("Usage: %s [file]" % command )


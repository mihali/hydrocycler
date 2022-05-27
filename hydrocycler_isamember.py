#=============
print ("\nWelcome to Hydrocycler (c) by Mihali Felipe (2022)\n")
print ("Checking against a previously generated signature file: ", end="")

import sys
from datetime import datetime
import pickle
import numpy as np
from scipy.spatial import KDTree as kd
import johnson
from hydrocycler_utils import *

ts = datetime.now().strftime("%y%m%d%H%M%S%f")

#==============
# Log all output 

class Logger(object):

    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("hydrocycler_isamember-%s.out"%ts, "a")
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message) 
    def flush(self): 
        self.terminal.flush(); 
        self.log.flush()

sys.stdout = Logger()


#==============
 
def fn (sig_dict, xyzfile):

    string = open (xyzfile, 'r' ).readlines()[2:] 
    [ocoords, hcoords, xcoords, xyzdictmaster ] = createarrays (string) 
    [ograph, trio] = findcycles (ocoords, hcoords)
    cycles = sorted(tuple(johnson.simple_cycles(ograph)))
    cycles_sig = tuple([ x for sublist in cycles for x in sublist ])
    if cycles_sig in sig_dict:
        print ("%s is in the H-bond family defined by signature file."%xyzfile)
    else:
        print ("%s is NOT in H-bond family defined by signature file."%xyzfile)

files = sys.argv[1:]
argc  = len(sys.argv)
command = sys.argv[0].split('/')[len(sys.argv[0].split('/')) - 1] 
if not files:
   print("Usage: %s signaturefile.sig cartesianfiles.xyz ..." % command )
   print ("No stdin option for this command. Use file argument")
elif argc > 2:
   print("%s\n"%files[0])
   with open(files[0], 'rb') as signature:
       sig_db = signature.read()
   sig_dict = pickle.loads(sig_db)  
   for file in files[1:]:
        fn ( sig_dict, file )
   print ("\n\nThank you for using Hydrocycler!\n\n")
else:
   print("Usage: %s signaturefile.sig cartesianfiles.xyz ..." % command )



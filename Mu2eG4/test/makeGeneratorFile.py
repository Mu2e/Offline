from ROOT import TFile, TTree, TDirectoryFile, TNtuple, TChain
import getopt
import sys
import array
import time
import math

######### Parse command line ###########

USAGE = """
Usage: python makeGeneratorFile.py [options]
Options:
   -i, --input  {name} : input ROOT file
   -o, --output {name} : output data file
   -d, --debug         : debug printout
   -p, --pdg {pdgCode} : PDG code of particles to be extracted
       --pmin {MeV/c}  : minimum momentum to be extracted (default:0)
   -v, --vd {VD_id}    : id of the virtual detector to be saved (default:8)
   -h, --help          : this message
"""
if len(sys.argv) < 2:
    print USAGE
    sys.exit(0)

try:
    opts,args=getopt.getopt(sys.argv[1:],'i:o:p:v:hd',
                            ['help','input=','output=','pdg=','vd=',
                             'pmin=','debug'])
except getopt.GetoptError:
    print USAGE
    sys.exit(2)

debug_print = False
input_file = None
output_file = None
pdg = None
vd = 8
pmin = 0.0

for o,a in opts:
    if o in ('-h','--help'):
        print USAGE
        sys.exit()
    elif o in ('-i','--input'):
        input_file=a
    elif o in ('-o','--output'):
        output_file=a
    elif o in ('-p','--pdg'):
        try:
            pdg = int(a)
        except:
            print "Wrong PDG code"
    elif o in ('-v','--vd'):
        try:
            vd = int(a)
        except:
            vd = None
            print "Wrong virtual detector id"
    elif o in ('--pmin'):
        try:
            pmin = float(a)
        except:
            pmin = None
            print "Wrong minimum momentum: ",a
    elif o in ('-d','--debug'):
        debug_print = True

all_defined = True

if input_file is None :
    print "Input file is not defined"
    all_defined = False

if output_file is None :
    print "Output file is not defined"
    all_defined = False

if pdg is None :
    print "Particle PDG code is not specified"
    all_defined = False

if vd is None :
    print "Virtual detector id is not specified"
    all_defined = False

if pmin is None :
    print "Minimum momentum is not specified"
    all_defined = False

if not all_defined :
    print USAGE
    sys.exit(2)

########## Main part ###############

# Open output file

fout = open(output_file,"w")

# Open ROOT file

chain = TChain("readvd/ntvd")
chain.Add(input_file)

nevents = chain.GetEntries()
print nevents," events found in ntuples"

tnow = time.time()
nwrite = 0

for i in xrange(nevents):

    if i%1000==0 and time.time()>(tnow+5) :
        print i," event read, ",nwrite," events written"
        tnow = time.time()
        
    # get the next tree in the chain and verify
    ientry = chain.LoadTree( i )
    if ientry < 0:
        break
    
    # copy next entry into memory and verify
    nb = chain.GetEntry( i )
    if nb <= 0:
        continue

    if( chain.pdg == pdg and chain.sid == vd ) :
        ptot = math.sqrt(chain.px*chain.px+chain.py*chain.py+chain.pz*chain.pz)
        if ptot<pmin : continue
        str = " %.2f %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d %.2f" % (
            chain.x, chain.y, chain.z,
            chain.px, chain.py, chain.pz,
            chain.time,
            int(chain.pdg), int(chain.evt), int(chain.trk), 0, 1.0)
        fout.write(str+"\n")
        nwrite+=1

fout.close()

    

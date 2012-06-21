from ROOT import TFile, TTree, TDirectoryFile, TNtuple, TChain
import getopt
import sys
import array
import time
import math

######### Parse command line ###########

USAGE = """
Generate file with particles at particular virtual detector, which can be used
as input of the event generator.

Usage: python makeGeneratorFile.py [options] [file list]
Options:
   -i, --input  {name} : input ROOT file
   -o, --output {name} : output data file
   -d, --debug         : debug printout
   -p, --pdg {pdgCode} : PDG code of particles to be extracted
       --pmin {MeV/c}  : minimum momentum to be extracted (default:0)
   -v, --vd {VD_id}    : id of the virtual detector to be saved (default:8)
   -t, --tau {lifetime}: use this lifetime in (ns) to calculate weight from
                       : proper time. This is used in the cases, when decays
                       : are switched off.
       --sin {sin(th)} : minimum sin(Theta) to be extracted (default:0)
   -h, --help          : this message

If -i option is used, only one input file is read. If -i option is not used,
all input files specified at the end are read.

Examples:

1. Save muons just in front of collimator 3

python makeGeneratorFile.py -i beamline_02.root -o muons_before_coll3.txt -p 13 -v 3

2. Save muons just after collimator 5 (DS entrance), using many input files

python makeGeneratorFile.py -o muons_before_coll3.txt -p 13 -v 8 /mu2e/data/users/logash/out/59919/59919_*/beamline_p.root

"""
if len(sys.argv) < 2:
    print USAGE
    sys.exit(0)

try:
    opts,args=getopt.getopt(sys.argv[1:],'i:o:p:v:hdt:',
                            ['help','input=','output=','pdg=','vd=',
                             'pmin=','debug','tau=', 'sin='])
except getopt.GetoptError:
    print USAGE
    sys.exit(2)

debug_print = False
input_file = None
output_file = None
pdg = None
vd = 8
pmin = 0.0
tau = None
sinth = 0.0

for o,a in opts:
    if o in ('-h','--help'):
        print USAGE
        sys.exit()
    elif o in ('-i','--input'):
        input_file=[ a ]
    elif o in ('-o','--output'):
        output_file=a
    elif o in ('-p','--pdg'):
        try:
            pdg = int(a)
        except:
            print "Wrong PDG code: ", a
    elif o in ('-t','--tau'):
        try:
            tau = float(a)
        except:
            print "Wrong lifetime: ", a
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
    elif o in ('--sin'):
        try:
            sinth = float(a)
        except:
            sinth = None
            print "Wrong minimum sin(th): ",a
    elif o in ('-d','--debug'):
        debug_print = True

all_defined = True

if input_file is None : input_file = args

if input_file is None :
    print "Input file is not defined"
    all_defined = False

if output_file is None :
    print "Output file is not defined"
    all_defined = False

if pdg is None :
    print "Particle PDG code is not specified"
    all_defined = False

if vd is None or vd<=0 :
    print "Virtual detector id is not specified or wrong number"
    all_defined = False

if pmin is None :
    print "Minimum momentum is not specified"
    all_defined = False

if sinth is None :
    print "Minimum sin(th) is not specified"
    all_defined = False

if not all_defined :
    print USAGE
    sys.exit(2)

########## Main part ###############

# Open output file

fout = open(output_file,"w")

# Open ROOT file

chain = TChain("readvd/ntpart")
for fname in input_file:
    chain.Add(fname)

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

    if chain.pdg != pdg : continue
    if chain.pvd[vd-1]<1.0 : continue
    
    ptot = chain.pvd[vd-1]
    if ptot<pmin : continue
    
    if sinth>0 :
        sth = math.sqrt(1-chain.pzvd[vd-1]**2/ptot**2)
        if sth<sinth : continue
        
    weight = chain.g4bl_weight
    if not (tau is None) and tau>0:
        weight = math.exp(-chain.gtvd[vd-1]/tau)
    str = " %.2f %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d %g" % (
        chain.xvd[vd-1], chain.yvd[vd-1], chain.zvd[vd-1],
        chain.pxvd[vd-1], chain.pyvd[vd-1], chain.pzvd[vd-1],
        chain.tvd[vd-1],
        int(chain.pdg), int(chain.evt), int(chain.trk), 0, weight)
    fout.write(str+"\n")
    nwrite+=1

fout.close()



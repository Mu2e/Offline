# Tracker Alignment with Millepede-II

Ryunosuke O'Neil
([@ryuwd](https://github.com/ryuwd))

A Tracker Alignment utility for performing Track-based alignment with the [Millepede-II software package](https://www.desy.de/~kleinwrt/MP2/doc/html/index.html). Only no-field cosmic tracks are currently supported.

### Alignment set-up
![Alignment flow(2)](https://user-images.githubusercontent.com/56410978/82936768-fa2e6500-9f86-11ea-81fe-b9f0bf20e842.png)
- 'Digis' refers to art files containing 'digi' data products. I am using DS-cosmic-nofield or DS-cosmic-nofield-alignselect
- Track Reco uses Richie's Time Fit ( [Mu2e-doc-33162](https://mu2e-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=33162) )
- 'MILLE' means to write input files for the 'PEDE' executable. This is implemented in `MilleDataWriter`.
- 'Align Tracker' refers to the stage carried out by `AlignedTrackerMaker` (TrackerConditions). It uses alignment constants provided to the Proditions interface to change tracker straw positions accordingly.
- 'Track Collection' is where the `AlignTrackCollector` module selects tracks, calculates the needed quantities for alignment, and writes the data to file using `MilleDataWriter`.
- 'PEDE' refers to the millepede executable that performs the alignment fit given the output file(s) produced by the `MilleDataWriter` class during the `AlignTrackCollector` job(s).


### MC Alignment test
1. Set up this base release.
```bash
source <Offline release>/setup.sh
```
2. Set up a TrackerAlignment helper environment.
Note that this will install some necessary python packages in your home directory with pip
```bash
source ${MU2E_BASE_RELEASE}/TrackerAlignment/scripts/setup.sh
```

2. Choose a working directory and setup your starting misalignment
```
cd /mu2e/data/users/$USER/
mkdir alignment-test
cd alignment-test

mu2ealign new MT_MDC2018_PlaneTranslationOnly 

# generates job.fcl and alignconstants_in.txt - the Tracker is misaligned according to the 'MT_MDC2018_PlaneTranslationOnly' configuration.
# (the corresp. DbService text file in ${MU2E_BASE_RELEASE}/TrackerAlignment/test/misalignments/ is copied)

```

It's important at this stage to choose which planes to fix in the alignment, in order to suppress weak modes. 
Fixing planes 5, 30 to zero in all DOFs, for example. This means changing the values in alignconstants_in.txt for Plane 5 + 30, and ensuring that these parameters are set in job.fcl:
```
physics.analyzers.AlignTrackCollector.ConstrainStrategy : "Fix"
physics.analyzers.AlignTrackCollector.FixPlane : [ 5, 30 ]

```


3. Run Track Collection + PEDE
You can run one process to collect the track data like this:
```bash
# running one job only
mu2ealign run
```

Or, you may run multiple processes:
```bash
# running multiple jobs e.g. 4 jobs processing 8 input files, 
# or 8 / 4 = 2 input art files per job
mu2ealign_genparallel 4 8
mu2ealign run
```

Once the jobs finish, you can run:
```bash
mu2ealign pede
```
Now you should make a new directory, and import the produced alignment constants:
```bash
mkdir iter1 && cd iter1
mu2ealign new ../alignconstants_out.txt
```
Then return to the beginning of this step.

Or, you can run multiple iterations automatically (easier, recommended)
```bash
# 
mu2ealign autorun 5 # for 5 alignment iterations
```

4. Examine output
```
# examine track diagnostics if needed
aligntrack_display TrackDiag.root < other trackdiag.root files to compare against >


# plot and show shifts and pulls
python ${MU2E_BASE_RELEASE}/TrackerAlignment/scripts/make_shiftplot.py alignconstants_out.txt iter*/alignconstants_out.txt
```
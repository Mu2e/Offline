# Tracker Alignment with Millepede-II

Write me!

### Expected alignment workflow
![Alignment flow(2)](https://user-images.githubusercontent.com/56410978/82936768-fa2e6500-9f86-11ea-81fe-b9f0bf20e842.png)
- 'Digis' refers to art files containing 'digi' data products. I am using DS-cosmic-nofield or DS-cosmic-nofield-alignselect
- Track Reco uses Richie's Time Fit ( [Mu2e-doc-33162](https://mu2e-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=33162) )
- `Mille` is a helper class provided by the Millepede alignment package to write input files for the `PEDE` executable, which performs the alignment fit.
- 'Align Tracker' refers to the stage carried out by `AlignedTrackerMaker` (TrackerConditions). It uses alignment constants provided to the Proditions interface to change tracker straw positions accordingly.
- 'Track Collection' is where the `AlignTrackCollector` module selects tracks, calculates the needed quantities for alignment, and writes the data to file using `Mille`.
- 'PEDE' refers to the millepede executable that performs the alignment fit given the output file(s) produced by the `Mille` class during the `AlignTrackCollector` job(s).


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
# (uses the corresp. DbService text file in ${MU2E_BASE_RELEASE}/TrackerAlignment/test/misalignments/)
mu2ealign new PlaneTranslationOnly 
```

3. Run Track Collection
```bash
mu2e -c job.fcl -S ${DS_COSMIC_NOFIELD_ALIGNSELECT} -n 25000

# examine track diagnostics if needed
aligntrack_display TrackDiag.root
```

4. Run PEDE step (alignment fit)
```bash
# adjust steering file if needed ...

mu2ealign pede


```

5. Start again from step 2, using the generated alignment constants to create a new job config
```bash
mkdir iter1 && cd iter1
mu2ealign new ../alignconstants_out.txt
```

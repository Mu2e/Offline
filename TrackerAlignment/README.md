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
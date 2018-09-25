#!/bin/bash
FCOUNT=0
#for mod in "RSD" "TTmakeSH" "TTmakePH" "TTflagBkgHits" "TTtimeClusterFinder"  "TThelixFinder" "TTKSFDeM" "TTKSFDeP" "CalHelixFinderDem" "CalHelixFinderDep" "CalSeedFitDem" "CalSeedFitDep"; do
for mod in "RSD" "makeSH" "makePH" "FlagBkgHits" "TimeClusterFinder"  "HelixFinder" "KSFDeM" "KSFDeP" "CaloCluserFast" "CalTimePeakFinder" "CalHelixFinderDem" "CalHelixFinderDep" "CalSeedFitDem" "CalSeedFitDep"; do
echo Processing module $mod
FCOUNT=`expr $FCOUNT + 1`
sqlite3  -separator "," $1 "SELECT Event, Time, ModuleLabel FROM TimeModule WHERE ModuleLabel=\"${mod}\";" > $FCOUNT$mod.csv
done


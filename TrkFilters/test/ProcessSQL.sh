#!/bin/bash
FCount = 0
for mod in "makeSH" "MakeStereoHits" "FlagBkgHits" "TimeClusterFinder"  "HelixFinder" "KSFDeM" ; do
echo Processing module $mod
FCOUNT=`expr $FCOUNT + 1`
sqlite3  -separator "," $1 "SELECT Event, Time, ModuleLabel FROM TimeModule WHERE ModuleLabel=\"${mod}\";" > $FCOUNT$mod.csv
done


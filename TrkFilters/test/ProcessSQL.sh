#!/bin/bash
for mod in "makeSH" "FSHPreStereo" "MakeStrawHitPositions" "MakeStereoHits" "FlagStrawHits" "FlagBkgHits" "TimeClusterFinder" "TCFilter" "PosHelixFinder" "PosHelixFilter" "NegHelixFinder" "NegHelixFilter" "KSFDeM" "DeMSeedFilter" "KSFUeMi" "UeMSeedFilter"; do
echo Processing module $mod
sqlite3  -separator "," $1 "SELECT Event, Time, ModuleLabel FROM TimeModule WHERE ModuleLabel=\"${mod}\";" > $mod.csv
done


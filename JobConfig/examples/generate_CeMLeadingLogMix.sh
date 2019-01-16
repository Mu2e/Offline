# Note: the jobs created by this script need at least 8GBytes of memory
generate_fcl --description CeMLeadingLogMix --dsconf MDC2018e --dsowner mu2e --include JobConfig/mixing/CeMLeadingLogMix.fcl \
--run-number 1002 --events-per-job 500 --njobs 1000 --max-engines 30 \
--auxinput=3:physics.filters.ootMixerCRV.fileNames:oot-CRV-recat.txt \
--auxinput=1:physics.filters.neutronMixerCRV.fileNames:neutron-CRV-cat.txt \
--auxinput=1:physics.filters.dioMixerCRV.fileNames:dio-CRV-cat.txt \
--auxinput=1:physics.filters.photonMixerCRV.fileNames:photon-CRV-cat.txt \
--auxinput=5:physics.filters.PSMixerCRV.fileNames:PS-CRV-recat.txt \
--auxinput=5:physics.filters.TSMixerCRV.fileNames:TS-CRV-recat.txt \
--auxinput=5:physics.filters.DSMixerCRV.fileNames:DS-CRV-recat.txt \
--auxinput=3:physics.filters.ootMixerTrkCal.fileNames:oot-TrkCal-recat.txt \
--auxinput=1:physics.filters.neutronMixerTrkCal.fileNames:neutron-TrkCal-cat.txt \
--auxinput=1:physics.filters.dioMixerTrkCal.fileNames:dio-TrkCal-cat.txt \
--auxinput=1:physics.filters.photonMixerTrkCal.fileNames:photon-TrkCal-cat.txt \
--auxinput=2:physics.filters.flashMixerTrkCal.fileNames:DS-flash-TrkCal-recat.txt \
--auxinput=1:physics.filters.protonMixerTrkCal.fileNames:proton-TrkCal.txt \
--auxinput=1:physics.filters.deuteronMixerTrkCal.fileNames:deuteron-TrkCal.txt 
rm -rf CeMLeadingLogMix
mv 000 CeMLeadingLogMix

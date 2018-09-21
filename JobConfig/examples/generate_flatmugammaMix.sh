generate_fcl --description flatmugammaMix --dsconf MDC2018a --dsowner mu2e --include JobConfig/mixing/flatmugammaMix.fcl \
--run-number 1002 --events-per-job 200000 --njobs 5000 --max-engines 30 \
--auxinput=15:physics.filters.ootMixerCRV.fileNames:oot-CRV-cat.txt \
--auxinput=1:physics.filters.neutronMixerCRV.fileNames:neutron-CRV-cat.txt \
--auxinput=1:physics.filters.dioMixerCRV.fileNames:dio-CRV-cat.txt \
--auxinput=1:physics.filters.photonMixerCRV.fileNames:photon-CRV-cat.txt \
--auxinput=25:physics.filters.PSMixerCRV.fileNames:PS-CRV-cut.txt \
--auxinput=25:physics.filters.TSMixerCRV.fileNames:TS-CRV-cat.txt \
--auxinput=25:physics.filters.DSMixerCRV.fileNames:DS-CRV-cut.txt \
--auxinput=15:physics.filters.ootMixerTrkCal.fileNames:oot-TrkCal-cat.txt \
--auxinput=1:physics.filters.neutronMixerTrkCal.fileNames:neutron-TrkCal-cat.txt \
--auxinput=1:physics.filters.dioMixerTrkCal.fileNames:dio-TrkCal-cat.txt \
--auxinput=1:physics.filters.photonMixerTrkCal.fileNames:photon-TrkCal-cat.txt \
--auxinput=10:physics.filters.flashMixerTrkCal.fileNames:DS-flash-TrkCal-cut.txt \
--auxinput=1:physics.filters.protonMixerTrkCal.fileNames:proton-TrkCal.txt \
--auxinput=1:physics.filters.deuteronMixerTrkCal.fileNames:deuteron-TrkCal.txt 
rm -rf flatmugammaMix
mv 000 flatmugammaMix

generate_fcl --description NoPrimaryFlat --dsconf MDC2018e --dsowner mu2e --embed JobConfig/mixing/NoPrimaryFlat.fcl \
--run-number 1002 --events-per-job 100 --njobs 10 --max-engines 30 \
--auxinput=45:physics.filters.ootMixerCRV.fileNames:oot-CRV-cat.txt \
--auxinput=3:physics.filters.neutronMixerCRV.fileNames:neutron-CRV-cat.txt \
--auxinput=3:physics.filters.dioMixerCRV.fileNames:dio-CRV-cat.txt \
--auxinput=3:physics.filters.photonMixerCRV.fileNames:photon-CRV-cat.txt \
--auxinput=25:physics.filters.PSMixerCRV.fileNames:PS-CRV-cut.txt \
--auxinput=25:physics.filters.TSMixerCRV.fileNames:TS-CRV-cat.txt \
--auxinput=25:physics.filters.DSMixerCRV.fileNames:DS-CRV-cut.txt \
--auxinput=45:physics.filters.ootMixerTrkCal.fileNames:oot-TrkCal-cat.txt \
--auxinput=3:physics.filters.neutronMixerTrkCal.fileNames:neutron-TrkCal-cat.txt \
--auxinput=3:physics.filters.dioMixerTrkCal.fileNames:dio-TrkCal-cat.txt \
--auxinput=3:physics.filters.photonMixerTrkCal.fileNames:photon-TrkCal-cat.txt \
--auxinput=30:physics.filters.flashMixerTrkCal.fileNames:DS-flash-TrkCal-cut.txt \
--auxinput=3:physics.filters.protonMixerTrkCal.fileNames:proton-TrkCal.txt \
--auxinput=3:physics.filters.deuteronMixerTrkCal.fileNames:deuteron-TrkCal.txt 
rm -rf NoPrimaryFlat
mv 000 NoPrimaryFlat

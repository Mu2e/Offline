generate_fcl --description DS-cosmicMix --dsconf MDC2018i --dsowner mu2e --include JobConfig/mixing/DS-cosmicMix.fcl \
--inputs=DS-cosmic.txt --merge-factor 5 --max-engines 30 \
--auxinput=3:physics.filters.ootMixerCRV.fileNames:oot-CRV-cat-cat.txt \
--auxinput=1:physics.filters.neutronMixerCRV.fileNames:neutron-CRV-cat.txt \
--auxinput=1:physics.filters.dioMixerCRV.fileNames:dio-CRV-cat.txt \
--auxinput=1:physics.filters.photonMixerCRV.fileNames:photon-CRV-cat.txt \
--auxinput=5:physics.filters.PSMixerCRV.fileNames:PS-CRV-cut-cat.txt \
--auxinput=5:physics.filters.TSMixerCRV.fileNames:TS-CRV-cat-cat.txt \
--auxinput=5:physics.filters.DSMixerCRV.fileNames:DS-CRV-cut-cat.txt \
--auxinput=3:physics.filters.ootMixerTrkCal.fileNames:oot-TrkCal-cat-cat.txt \
--auxinput=1:physics.filters.neutronMixerTrkCal.fileNames:neutron-TrkCal-cat.txt \
--auxinput=1:physics.filters.dioMixerTrkCal.fileNames:dio-TrkCal-cat.txt \
--auxinput=1:physics.filters.photonMixerTrkCal.fileNames:photon-TrkCal-cat.txt \
--auxinput=2:physics.filters.flashMixerTrkCal.fileNames:DS-flash-TrkCal-cut-cat.txt \
--auxinput=1:physics.filters.protonMixerTrkCal.fileNames:proton-TrkCal.txt \
--auxinput=1:physics.filters.deuteronMixerTrkCal.fileNames:deuteron-TrkCal.txt 
for dirname in 000 001 002 003 004 005 006 007 008 009; do
 if test -d $dirname; then
  echo "found dir $dirname"
  rm -rf DS-cosmicMix_$dirname
  mv $dirname DS-cosmicMix_$dirname
 fi
done

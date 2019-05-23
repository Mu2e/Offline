# Note: the jobs created by this script need at least 8GBytes of memory
generate_fcl --description ${name} --dsconf MDC2018h --dsowner mu2e ${includeOrEmbed} \
--run-number 1002 --events-per-job ${perjob} --njobs ${njobs} --max-engines 50 --first-subrun=1 \
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
rm -rf ${dirname}/${name}
rm -rf ${dirname}/${name}-001
rm -rf ${dirname}/${name}-002
rm -rf ${dirname}/${name}-003
rm -rf ${dirname}/${name}-004
rm -rf ${dirname}/${name}-005
rm -rf ${dirname}/${name}-006
rm -rf ${dirname}/${name}-007
mv 000 ${dirname}/${name}
mv 001 ${dirname}/${name}-001
mv 002 ${dirname}/${name}-002
mv 003 ${dirname}/${name}-003
mv 004 ${dirname}/${name}-004
mv 005 ${dirname}/${name}-005
mv 006 ${dirname}/${name}-006
mv 007 ${dirname}/${name}-007

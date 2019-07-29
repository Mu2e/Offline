# Note: the jobs created by this script need at least 8GBytes of memory
generate_fcl --description ${name} --dsconf MDC2018h --dsowner mu2e ${includeOrEmbed} \
--run-number 1002 --events-per-job ${perjob} --njobs ${njobs} --max-engines 50 --first-subrun=1 
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

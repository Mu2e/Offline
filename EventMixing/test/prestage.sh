#! /bin/bash
#
# Prestage mix-in files.
#
# Use the section numbers of the mix-in files in a deterministic way that
# depends only on the process number mod 7.
#

echo "Starting prestage: "  > prestage.log
echo "process: " $process  >> prestage.log
bgversion=$(($(($process%7))+1))
echo "background file version: " $bgversion >> prestage.log

type ifdh  >> prestage.log

echo "/mu2e/data/users/kutschke/Backgrounds/130922/Merged/dio/Tracker/dioBG_${bgversion}_data.root prestage/dioBG_data.root"             >> prestage.txt
echo "/mu2e/data/users/kutschke/Backgrounds/130922/Merged/proton/Tracker/protonBG_${bgversion}_data.root prestage/protonBG_data.root"    >> prestage.txt
echo "/mu2e/data/users/kutschke/Backgrounds/130922/Merged/neutron/Tracker/neutronBG_${bgversion}_data.root prestage/neutronBG_data.root" >> prestage.txt
echo "/mu2e/data/users/kutschke/Backgrounds/130922/Merged/photon/Tracker/photonBG_${bgversion}_data.root prestage/photonBG_data.root"    >> prestage.txt

mkdir prestage

tstart=$(date +%s)


# Get the lock and copy files
ifdh cp -f prestage.txt >> prestage.log
ret=$?

t2=$(date +%s)
echo "$(date) # Total stage-in time: $((t2-tstart)) seconds, status $ret"  >> prestage.log

ls -sh prestage >> prestage.log

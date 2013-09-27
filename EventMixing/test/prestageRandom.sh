#! /bin/bash
#
# Prestage mix-in files.
# Randomize the section numbers from each background species.
#

nSections=7

s1=$(($(($RANDOM%${nSections}))+1))
s2=$(($(($RANDOM%${nSections}))+1))
s3=$(($(($RANDOM%${nSections}))+1))
s4=$(($(($RANDOM%${nSections}))+1))

echo "Starting prestage: "                 > prestage.log
echo "process: " $process                 >> prestage.log
echo "Number of sections:   "  $nSections >> prestage.log
echo "Dio file section:     "  $s1        >> prestage.log
echo "Proton file section:  "  $s2        >> prestage.log
echo "Neutron file section: "  $s3        >> prestage.log
echo "Photon file section:  "  $s4        >> prestage.log

type ifdh  >> prestage.log


echo "/mu2e/data/users/kutschke/Backgrounds/130922/Merged/dio/Tracker/dioBG_${s1}_data.root prestage/dioBG_data.root"             >> prestage.txt
echo "/mu2e/data/users/kutschke/Backgrounds/130922/Merged/proton/Tracker/protonBG_${s2}_data.root prestage/protonBG_data.root"    >> prestage.txt
echo "/mu2e/data/users/kutschke/Backgrounds/130922/Merged/neutron/Tracker/neutronBG_${s3}_data.root prestage/neutronBG_data.root" >> prestage.txt
echo "/mu2e/data/users/kutschke/Backgrounds/130922/Merged/photon/Tracker/photonBG_${s4}_data.root prestage/photonBG_data.root"    >> prestage.txt

mkdir prestage

cat prestage.txt

tstart=$(date +%s)


# Get the lock and copy files
ifdh cp -f prestage.txt >> prestage.log
ret=$?

t2=$(date +%s)
echo "$(date) # Total stage-in time: $((t2-tstart)) seconds, status $ret"  >> prestage.log

ls -sh prestage >> prestage.log

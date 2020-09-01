#!/bin/bash

#Runs the TEve Event Display for Mu2e

echo "Running TEve Event Display for Mu2e..."

fcl=$1
art=$2
numofevts=$3

for var in "$@"
  do
	if [ "$var" == '-3DOnly' ]; then
	    sed -i 's/show2D.*/show2D : false/' TEveEventDisplay/fcl/prolog.fcl
	fi

	if [ "$var" == '-2Dand3D' ]; then
	    sed -i 's/show2D.*/show2D : true/' TEveEventDisplay/fcl/prolog.fcl
	fi

	if [ "$var" == '-DSOnly' ]; then
	    sed -i 's/showDSOnly :.*/showDSOnly : true/' TEveEventDisplay/fcl/prolog.fcl
	    sed -i 's/showCRV :.*/showCRV : false/' TEveEventDisplay/fcl/prolog.fcl
	fi
  
	if [ "$var" == '-CRVOnly' ]; then
	    sed -i 's/showDSOnly :.*/showDSOnly : false/' TEveEventDisplay/fcl/prolog.fcl
	    sed -i 's/showCRV :.*/showCRV : true/' TEveEventDisplay/fcl/prolog.fcl
	fi

	if [ "$var" == '-DSandCRV' ]; then
	    sed -i 's/showDSOnly :.*/showDSOnly : true/' TEveEventDisplay/fcl/prolog.fcl
	    sed -i 's/showCRV :.*/showCRV : true/' TEveEventDisplay/fcl/prolog.fcl
	fi

	if [ "$var" == 'hits' ]; then
	    sed -i 's/addHits.*/addHits : true/' TEveEventDisplay/fcl/prolog.fcl
	fi

	if [ "$var" == 'clusters' ]; then
	    sed -i 's/addClusters.*/addClusters : true/' TEveEventDisplay/fcl/prolog.fcl
	fi

	if [ "$var" == 'tracks' ]; then
	    sed -i 's/addTracks.*/addTracks : true/' TEveEventDisplay/fcl/prolog.fcl
	fi

	if [ "$var" == 'crvinfo' ]; then
	    sed -i 's/addCrvHits.*/addCrvHits : true/' TEveEventDisplay/fcl/prolog.fcl
	fi
done

mu2e -c ${fcl} ${art} --nevts ${numofevts}


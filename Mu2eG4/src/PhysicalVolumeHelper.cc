//
// A utility class to do indexolgy related to persistence of
// physical volume information.
//
// $Id: PhysicalVolumeHelper.cc,v 1.1 2010/03/23 20:58:29 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/03/23 20:58:29 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) This should be instantiated after the G4 initialization phase is complete.
//    ( After  _runManager->BeamOnBeginRun(); ).
//    At this time the information in PhysicalVolumeStore is static until 
//    the end of the run; in particular, the addresses of all physical
//    volumes are valid until the end of the run.
//
// 2) One purpose of the class is to construct the _persistentInfo object
//    that will be copied into the run-data.  This provides persistent
//    information about G4 volumes.
//
// 3) The second purpose is, given an G4VPhysicalVolume*, find the
//    index into _persistentInfo that corresponds to the given volume.
//
// 4) The implementation of 3) is a map that uses G4VPhysicalVolume*
//    as its key.  This is safe because of 1).  And this is faster
//    than the alternate solution which use a key that is a std::pair of
//    volume name and copy number.
//
// 5) The implementation will always correctly satisfy the request or 
//    will throw. The client code need not check for errors.

// Framework includes
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"

// G4 includes
#include "G4PhysicalVolumeStore.hh"
#include "G4Track.hh"

using namespace std;

namespace mu2e {
  
  PhysicalVolumeHelper::PhysicalVolumeHelper():
    _persistentInfo(),
    _volumeMap(){
  }

  // Return the index into _persistentInfo for the volume attached to this track.
  int PhysicalVolumeHelper::index( const G4Track* track ) const{
    G4VPhysicalVolume* physVol  = track->GetVolume();
    VolMapType_const_iterator physVolIter = _volumeMap.find(physVol);
    if ( physVolIter == _volumeMap.end() ){
      throw cms::Exception("RANGE")
        << "Cannot find physical volume in the map: "
        << physVol->GetName()
        << " Copy number: "
        << physVol->GetCopyNo()
        << " Size of map: "
        << _volumeMap.size()
        << "\n";
    }
    return physVolIter->second;
  }

  // Return the index into _persistentInfo for the volume attached to this track.
  int PhysicalVolumeHelper::index( G4VPhysicalVolume* physVol ) const{

    VolMapType_const_iterator physVolIter = _volumeMap.find(physVol);
    if ( physVolIter == _volumeMap.end() ){
      string message;
      if ( _volumeMap.empty() ){
        message = "The map is empty.  Is it not yet intialized?";
      } else{
        message = "The map is not empty: something has been corrupted.";
      }
      throw cms::Exception("RANGE")
        << "Cannot find physical volume in the map. "
        << message
        << "\n";
    }
    return physVolIter->second;

  }


  // Build _persistentInfo and _volumeMap from G4PhysicalVolumeStore.
  void PhysicalVolumeHelper::beginRun(){

    _volumeMap.clear();
    _persistentInfo.clear();

    // Loop over physical volume store.
    G4PhysicalVolumeStore* pstore = G4PhysicalVolumeStore::GetInstance();
    for ( std::vector<G4VPhysicalVolume*>::const_iterator i=pstore->begin(); i!=pstore->end(); ++i){

      // Make sure that this volume is not yet in the map.
      G4VPhysicalVolume* vpv = *i;
      VolMapType_iterator iter = _volumeMap.find(vpv);
      if ( iter != _volumeMap.end() ){
        throw cms::Exception("RANGE")
          << "Error building the persistent volume list.  Volume: "
          << vpv->GetName()
          << " Copy number: "
          << vpv->GetCopyNo()
          << " already exisist!\n";
      }

      // Add volume to the map and to the persistent info.
      _volumeMap.insert( std::make_pair(vpv,_volumeMap.size()) );
      _persistentInfo.push_back( PhysicalVolumeInfo( vpv->GetName(), vpv->GetCopyNo() ));
    }

  }

  // Clear information at the end of a run.
  void PhysicalVolumeHelper::endRun(){
    _volumeMap.clear();
    _persistentInfo.clear();
  }

}

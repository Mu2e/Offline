//
// A utility class to do indexolgy related to persistence of
// physical volume information.
//
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
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Offline/Mu2eG4/inc/PhysicalVolumeHelper.hh"

// G4 includes
#include "Geant4/G4PhysicalVolumeStore.hh"
#include "Geant4/G4Track.hh"

using namespace std;

namespace mu2e {

  PhysicalVolumeHelper::PhysicalVolumeHelper():
    _volumeMap(),
    m_helperInitialized(false)
    {}

  // Return the index into _persistentInfo for the volume attached to this track.
  int PhysicalVolumeHelper::index( const G4Track* track ) const{
    G4VPhysicalVolume* physVol  = track->GetVolume();
    VolMapType_const_iterator physVolIter = _volumeMap.find(physVol);
    if ( physVolIter == _volumeMap.end() ){
      throw cet::exception("RANGE")
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
      throw cet::exception("RANGE")
        << "Cannot find physical volume in the map. "
        << message
        << "\n";
    }
    return physVolIter->second;

  }


  // Build _persistentInfo and _volumeMap from G4PhysicalVolumeStore.
  void PhysicalVolumeHelper::beginRun(){

    _volumeMap.clear();
    _pSingleStage.clear();

    // Loop over physical volume store.
    G4PhysicalVolumeStore* pstore = G4PhysicalVolumeStore::GetInstance();
    unsigned current(0u);
    for (auto i=pstore->begin(); i!=pstore->end(); ++i, ++current){

      // Add volume to the map; it's an error if its already there.
      G4VPhysicalVolume* vpv = *i;
      pair<VolMapType_iterator,bool> ret = _volumeMap.insert( std::make_pair(vpv, current) );
      if ( !ret.second ){
        throw cet::exception("RANGE")
          << "Error building the persistent volume list.  Volume: "
          << vpv->GetName()
          << " Copy number: "
          << vpv->GetCopyNo()
          << " already exisist!\n";
      }

      //scorer volumes have no material and should be omitted
      if (!vpv->GetLogicalVolume()->GetMaterial()) {
        std::cout<<vpv->GetName()<<" "<<vpv->GetTranslation()<<" "
                 <<vpv->IsReplicated()<<" "<<vpv->IsParameterised()<<" "
                 <<vpv->GetLogicalVolume()->GetName()<<" "
                 <<std::endl;
      }
      if (!vpv->GetLogicalVolume()->GetMaterial()) continue;

      _pSingleStage[cet::map_vector_key(current)] =
        PhysicalVolumeInfo(vpv->GetName(), vpv->GetCopyNo(), vpv->GetLogicalVolume()->GetMaterial()->GetName() );
    }

    m_helperInitialized = true;

  }

  // Clear information at the end of a run.
  void PhysicalVolumeHelper::endRun(){
    _volumeMap.clear();
  }

}

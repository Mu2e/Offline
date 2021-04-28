// Shared code to write physical volume info by G4, G4MT
// into subrun or event (the mixing use case).
//
// Andrei Gaponenko, 2021

#ifndef Mu2eG4_writePhysicalVolumes_hh
#define Mu2eG4_writePhysicalVolumes_hh

#include <string>
#include <optional>

#include "canvas/Utilities/InputTag.h"

#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"

namespace mu2e {

  // Returns the current simStage index derived from the
  // PhysicalVolumeInfoMultiCollection size.
  template <class PRINCIPAL>
  unsigned writePhysicalVolumes(PRINCIPAL& store,
                                const std::optional<art::InputTag>& input,
                                const PhysicalVolumeInfoSingleStage& vi,
                                const std::string& outInstanceName);


}

#endif/*Mu2eG4_writePhysicalVolumes_hh*/

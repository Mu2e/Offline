//A. Edmonds based on RMCWeight_module.cc by Hasung Song (2018)

// C++ includes
#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

// cetlib includes
#include "cetlib_except/exception.h"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

// Mu2e includes
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/EventWeight.hh"

#include "Offline/TrkDiag/inc/WeightModule.hh"
#include "Offline/Mu2eUtilities/inc/BinnedSpectrumWeightPhys.hh"
namespace mu2e {
  typedef WeightModule<BinnedSpectrumWeightPhys> BinnedSpectrumWeight;
}
DEFINE_ART_MODULE(mu2e::BinnedSpectrumWeight)

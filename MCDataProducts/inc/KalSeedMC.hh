#ifndef MCDataProducts_KalSeedMC_hh
#define MCDataProducts_KalSeedMC_hh
///
//  MC truth analog to a KalSeed (Kalman fit track)
//
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/CaloClusterMC.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include <Rtypes.h>
#include <utility>
#include <vector>

namespace mu2e {
// define some structs and types
  struct MCDigiRel {
    int _sindex; // index into TrkStrawHitSeed of KalSeed, -1 if the associated digi wasn't used in the track fit
    art::Ptr<StrawDigiMC> _mcd; // associated DigiMC
    MCRelationship _prel; // relationship of this digiMC to the primary particle
  };
  typedef std::pair<SPPtr,Float_t> SPW; // Sim Particle weighted by the fraction of active diAgis
  typedef art::Ptr<CaloClusterMC> CCMCPtr;
  typedef std::pair<CCMCPtr,Float_t >CCW; // calo cluster weighted by primary energy fraction
  typedef art::Ptr<StepPointMC> SPP; // used only for virual detectors
  struct KalSeedMC {
    SPW const& primaryParticle() const { return _parts.size() > 0 ? _parts[0] : _nullptr; }
    SPW const& simParticle(size_t index) const { return _parts.at(index); }
// data products
    std::vector<SPW> _parts; // particles associated with digis on this track, sorted by fraction of active track hits p
    std::vector<MCDigiRel> _mcdigis; // digi MC information for digis associated with this track
    // StepPoints sampled at the virtual detectors
    std::vector<SPP> _vds;
    CCW _calo; // MC truth of associated calorimeter cluster, with weight to primary
   // null ptr to satisfy interface
    static SPW _nullptr;
  };
  typedef std::vector<mu2e::KalSeedMC> KalSeedMCCollection;
}
#endif

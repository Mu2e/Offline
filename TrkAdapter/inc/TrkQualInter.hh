//
// Interface for computing the track quality MVA output
// Dave Brown (LBNL) 2/12/17
//
#ifndef TrkAdapter_TrkQualInter_hh
#define TrkAdapter_TrkQualInter_hh
// framework
#include "fhiclcpp/ParameterSet.h"
// Utilities
#include "Mu2eUtilities/inc/MVATools.hh"
// TrkQual struct
#include "RecoDataProducts/inc/TrkQual.hh"

namespace mu2e {
  class TrkInter;

  class TrkQualInter {
    public:
// construct from a parameter set
    explicit TrkQualInter(fhicl::ParameterSet const&);
// compute the MVA value in the input/output struct.  This will
// OVERWRITE any existing MVA value
    void fillMVA(TrkQual& trkqual) const;
// fill the MVA input struct and compute the MVA value given an interface to a track
    void fillAll(TrkInter const& tinter, TrkQual& trkqual) const;
    private:
// configuration data
    int _debug;
// track quality computation
    std::unique_ptr<MVATools> _trkqualmva;
  };

}
#endif

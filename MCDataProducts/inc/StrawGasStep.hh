#ifndef MCDataProducts_StrawGasStep_hh
#define MCDataProducts_StrawGasStep_hh
//
// Class to summarize the passage of a single particle through a single straw's gas volume
// This consolidates the G4 steps and insulates the downstream straw response ionization simulation 
// from details of the G4 model
//
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "cetlib/map_vector.h"

#include "MCDataProducts/inc/ProcessCode.hh"
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/XYZVec.hh"
#include <Rtypes.h>

namespace mu2e {
  class StrawGasStep {
    public:
      typedef cet::map_vector_key key_type;

      StrawGasStep() : _eIon(0.0), _pathLen(0.0), _mom(0.0), _time(0.0) {}
      StrawGasStep( key_type   simParticleKey, StrawId    strawId, 
	Float_t  ionizingEdep, Float_t    pathLength, Float_t width, Float_t    momentum, Double_t   time, 
	XYZVec const& startPosition, XYZVec const& endPosition) : _simpart(simParticleKey),
	_strawId(strawId), _eIon(ionizingEdep), _pathLen(pathLength), _width(width), _mom(momentum), _time(time),
	_startpos(startPosition), _endpos(endPosition) {}

      key_type   simParticleKey() const { return _simpart; }
      StrawId    strawId()    const { return _strawId;}
      Float_t    ionizingEdep()    const { return _eIon; }
      Float_t    pathLength()   const { return _pathLen; }
      Float_t    radialWidth()   const { return _width; }
      Float_t    momentum()     const { return _mom; }
      Double_t   time()         const { return _time; }
      XYZVec const& startPosition() const { return _startpos; }
      XYZVec const& endPosition() const { return _endpos; }

    private:

      key_type      _simpart; // key to the particle generating this edep
      StrawId       _strawId; // straw
      Float_t       _eIon;  // ionization energy deposited in this volue
      Float_t       _pathLen;  // Length this particle traveled in this straw gas: this is NOT necessarily the end-start distance
      Float_t       _width; // transverse radial width of the charge cloud
      Float_t	    _mom; // scalar momentum of the particle in the middle of this gas volume
      Double_t      _time; // time particle enters this gas volume; must be double to allow for long-lived particles
      XYZVec	    _startpos, _endpos; //entrance and exit to the gas volume
  };

  // this class's collection is a map_vector
  typedef cet::map_vector<StrawGasStep> StrawGasStepCollection;
  typedef art::Assns<StrawGasStep,StepPointMC> StrawGasStepAssns;
}


#endif


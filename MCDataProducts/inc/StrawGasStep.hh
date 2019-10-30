#ifndef MCDataProducts_StrawGasStep_hh
#define MCDataProducts_StrawGasStep_hh
//
// Class to summarize the passage of a single particle through a single straw's gas volume
// This consolidates the G4 steps and insulates the downstream straw response ionization simulation 
// from details of the G4 model
//
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"

#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/ProcessCode.hh"
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/XYZVec.hh"
#include <Rtypes.h>

namespace mu2e {
  class StrawGasStep {
    public:
      struct StepType {
	constexpr static uint8_t _smsk = 0xF; // mask for shape field
	constexpr static uint8_t _ssft = 0; // shift for shape field
	constexpr static uint8_t _imsk = 0xF0; // mask for ionization field
	constexpr static uint8_t _isft = 4; // shift for ionization field
	enum shape {line=0,arc,curl,point }; // shape of the trajectory within the straw
	enum ionization { minion=0, highion, neutral }; // type of ionization
	uint16_t _stype;

	StepType() : _stype(0) {}
	StepType(StepType::shape shp, StepType::ionization ion) :
	  _stype( shp | (ion << _isft)) {}
	shape() const { return _stype & _smsk; }
	ionization() const { return (_stype & _imsk) >> _isft; }

      }; 

      StrawGasStep() : _eIon(0.0), _pathLen(0.), _mom(0.0), _time(0.0) {}
      StrawGasStep( StrawId    strawId, StepType stype,
	Float_t  Edep, Float_t    pathLength, Float_t width, Double_t   time, 
	XYZVec const& startPosition, XYZVec const& endPosition) :
	_strawId(strawId), _stype(stype), _eIon(Edep),
	_pathLen(pathLength), _width(width), _time(time),
	_startpos(startPosition), _endpos(endPosition) {}

      StrawId    strawId()    const { return _strawId;}
      StepType   stepType()    const { return _stype; }
      Float_t    ionizingEdep()    const { return _eIon; }
      Float_t    pathLength()   const { return _pathLen; }
      Float_t    radialWidth()   const { return _width; } 
      Double_t   time()         const { return _time; }
      XYZVec const& startPosition() const { return _startpos; }
      XYZVec const& endPosition() const { return _endpos; }
    private:
      StrawId       _strawId; // straw
      StepType	    _stype; // type of step: used downstream in response simulation
      Float_t       _eIon;  // ionization energy deposited in this straw by this particle
      Float_t       _pathLen;  // Length the primary particle traveled in this straw gas: this is NOT necessarily the end-start distance
      Float_t       _width; // transverse RMS of the charge cloud WRT the wire
      Double_t      _time; // time particle enters this gas volume; must be double to allow for long-lived particles
      XYZVec	    _startpos, _endpos; //entrance and exit to the gas volume
  };

  typedef std::vector<StrawGasStep> StrawGasStepCollection;
  typedef art::Assns<StrawGasStep,StepPointMC> StrawGasStepAssns;
}

#endif


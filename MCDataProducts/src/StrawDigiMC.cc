//
//  Summary of MC information used to create a StrawDigi
//
// Original author David Brown, LBNL
//
// Mu2e includes
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// Framework includes.
#include "cetlib_except/exception.h"
// C++ includes
#include <ostream>

using namespace std;

namespace mu2e {

  // Default constructor is required for persistable classes
  StrawDigiMC::StrawDigiMC()
    : _strawid(StrawId::_invalid)
  {}

  StrawDigiMC::StrawDigiMC(StrawId sid, double wetime[2],
      CLHEP::HepLorentzVector cpos[2], 
      art::Ptr<StepPointMC> stepMC[2], vector<art::Ptr<StepPointMC> > const& stepMCs) :
    _strawid(sid), _stepMCs(stepMCs)
  {
    for(size_t strawend=0;strawend<2;++strawend){
      _wetime[strawend] = wetime[strawend];
      _cpos[strawend] = cpos[strawend];
      _stepMC[strawend] = stepMC[strawend];
    }
  }

  StrawDigiMC::StrawDigiMC(StrawDigiMC const& other ) : 
    _strawid(other._strawid), _stepMCs(other._stepMCs) {
    for(size_t strawend=0;strawend<2;++strawend){
      _wetime[strawend] = other._wetime[strawend];
      _cpos[strawend] = other._cpos[strawend];
      _stepMC[strawend] = other._stepMC[strawend];
    }
  }

  double StrawDigiMC::driftDistance(StrawEnd strawend) const {
    double retval = -100.0;
    if(!_stepMC[strawend].isNull()){
      const Tracker& tracker = getTrackerOrThrow();
      // use the MC true sid, not the straws sid (digi could be from x-talk)
      Straw const& straw = tracker.getStraw(_stepMC[strawend]->strawIndex());
      retval = (_cpos[strawend] - straw.getMidPoint()).perp(straw.getDirection());
    }
    return retval;
  }

  double StrawDigiMC::distanceToMid(StrawEnd strawend) const {
    double retval = -100.0;
    if(!_stepMC[strawend].isNull()){
      const Tracker& tracker = getTrackerOrThrow();
      Straw const& straw = tracker.getStraw(_stepMC[strawend]->strawIndex());
      retval =  (_cpos[strawend] - straw.getMidPoint()).dot(straw.getDirection());
    }
    return retval;
  }

  bool StrawDigiMC::isCrossTalk(StrawEnd strawend) const {
    bool retval(false);
    if(!_stepMC[strawend].isNull()){
    // this is currently broken till we replace starwIndex with strawId FIXME!!
      //retval = _strawid == _stepMC[strawend]->strawIndex();
    }
    return retval;
  }

  double StrawDigiMC::energySum() const {
    double esum(0.0);
    for(auto imcs = _stepMCs.begin(); imcs!= _stepMCs.end(); ++ imcs){
      esum += (*imcs)->eDep();
    }
    return esum;
  }


  double StrawDigiMC::triggerEnergySum(StrawEnd strawend) const {
    double esum(0.0);
    if(!_stepMC[(size_t)strawend].isNull()){
      for(auto imcs = _stepMCs.begin(); imcs!= _stepMCs.end(); ++ imcs){
	// if the simParticle for this step is the same as the one which fired the discrim, add the energy
	if( (*imcs)->simParticle() == _stepMC[(size_t)strawend]->simParticle() ||
	    (*imcs)->simParticle()->parent() == _stepMC[(size_t)strawend]->simParticle() ||
	    (*imcs)->simParticle() == _stepMC[(size_t)strawend]->simParticle()->parent() )
	  esum += (*imcs)->eDep();
      }
    }
    return esum;
  }

  // Print the information found in this hit.
  void StrawDigiMC::print( ostream& ost, bool doEndl ) const {

    ost << "Straw Digi MC Truth for straw ends " << StrawEnd(TrkTypes::cal) << " : " << StrawEnd(TrkTypes::hv)
      << " cluster times : "      << _cpos[0].t() << " : " << _cpos[1].t()
      << " drift distance: "      << driftDistance(TrkTypes::cal) << " : " << driftDistance(TrkTypes::hv)
      << " distance to wire center: "     << distanceToMid(TrkTypes::cal) << " : " << distanceToMid(TrkTypes::hv)
      << " Energy: " << energySum();

    if ( doEndl ){
      ost << endl;
    }

  }

}

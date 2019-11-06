//
//  Summary of MC information used to create a StrawDigi
//
// Original author David Brown, LBNL
//
// Mu2e includes
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "DataProducts/inc/StrawEnd.hh"
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
      art::Ptr<StrawGasStep> sgs[2], art::Ptr<StepPointMC> stepmc[2]) :
    _strawid(sid)
  {
    for(size_t strawend=0;strawend<2;++strawend){
      _wetime[strawend] = wetime[strawend];
      _cpos[strawend] = cpos[strawend];
      _sgs[strawend] = sgs[strawend];
      _stepMC[strawend] = stepmc[strawend];
    }
  }

// legacy constructor: StrawGasSteps will be empty!
  StrawDigiMC::StrawDigiMC(StrawId sid, double wetime[2], 
	CLHEP::HepLorentzVector cpos[2], art::Ptr<StepPointMC> stepMC[2], std::vector<art::Ptr<StepPointMC> > const& stepmcs) :_strawid(sid)
  {
    for(size_t strawend=0;strawend<2;++strawend){
      _wetime[strawend] = wetime[strawend];
      _cpos[strawend] = cpos[strawend];
      _stepMC[strawend] = stepMC[strawend];
    }
  }

  StrawDigiMC::StrawDigiMC(const StrawDigiMC& rhs) : _strawid(rhs.strawId()) {
    for(int i_end=0;i_end<StrawEnd::nends;++i_end){
      StrawEnd::End end = static_cast<StrawEnd::End>(i_end);
      _wetime[end] = rhs.wireEndTime(end);
      _cpos[end] = rhs.clusterPosition(end);
      _sgs[end] = rhs.strawGasStep(end);
      _stepMC[end] = rhs.stepPointMC(end);
    }
  }

  StrawDigiMC::StrawDigiMC(const StrawDigiMC& rhs, art::Ptr<StepPointMC> stepMC[2] ) : _strawid(rhs.strawId()) {
    for(int i_end=0;i_end<StrawEnd::nends;++i_end){
      StrawEnd::End end = static_cast<StrawEnd::End>(i_end);
      _wetime[end] = rhs.wireEndTime(end);
      _cpos[end] = rhs.clusterPosition(end);
      _sgs[end] = rhs.strawGasStep(end);
      _stepMC[end] = stepMC[i_end];
    }
  }

  double StrawDigiMC::driftDistance(StrawEnd strawend) const {
    double retval = -100.0;
    if(!_stepMC[strawend].isNull()){
      const Tracker& tracker = *GeomHandle<Tracker>();
      // use the MC true sid, not the straws sid (digi could be from x-talk)
      Straw const& straw = tracker.getStraw(_stepMC[strawend]->strawId());
      retval = (_cpos[strawend] - straw.getMidPoint()).perp(straw.getDirection());
    }
    return retval;
  }

  double StrawDigiMC::distanceToMid(StrawEnd strawend) const {
    double retval = -100.0;
    if(!_stepMC[strawend].isNull()){
      const Tracker& tracker = *GeomHandle<Tracker>();
      Straw const& straw = tracker.getStraw(_stepMC[strawend]->strawId());
      retval =  (_cpos[strawend] - straw.getMidPoint()).dot(straw.getDirection());
    }
    return retval;
  }

  bool StrawDigiMC::isCrossTalk(StrawEnd strawend) const {
    bool retval(false);
    if(!_stepMC[strawend].isNull()){
      retval = _strawid == _stepMC[strawend]->strawId();
    }
    return retval;
  }

  double StrawDigiMC::energySum() const {
    return _sgs[0]->ionizingEdep();
  }


  double StrawDigiMC::triggerEnergySum(StrawEnd strawend) const {
    return _sgs[strawend]->ionizingEdep();
  }

  // Print the information found in this hit.
  void StrawDigiMC::print( ostream& ost, bool doEndl ) const {

    ost << "Straw Digi MC Truth for straw ends " << StrawEnd(StrawEnd::cal) << " : " << StrawEnd(StrawEnd::hv)
      << " cluster times : "      << _cpos[0].t() << " : " << _cpos[1].t()
      << " drift distance: "      << driftDistance(StrawEnd::cal) << " : " << driftDistance(StrawEnd::hv)
      << " distance to wire center: "     << distanceToMid(StrawEnd::cal) << " : " << distanceToMid(StrawEnd::hv)
      << " Energy: " << energySum();

    if ( doEndl ){
      ost << endl;
    }

  }

}

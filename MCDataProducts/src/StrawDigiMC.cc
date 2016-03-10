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
#include "cetlib/exception.h"
// C++ includes
#include <ostream>

using namespace std;

namespace mu2e {

  // Default constructor is required for persistable classes
  StrawDigiMC::StrawDigiMC()
    : _strawIndex(StrawIndex::NO_STRAW)
  {}

  StrawDigiMC::StrawDigiMC(StrawIndex index, double wetime[2], CLHEP::HepLorentzVector cpos[2], 
      art::Ptr<StepPointMC> stepMC[2], vector<art::Ptr<StepPointMC> > const& stepMCs) :
    _strawIndex(index), _stepMCs(stepMCs)
  {
    for(int itdc=0;itdc<StrawEnd::nends;++itdc){
      size_t jtdc = static_cast<size_t>(itdc);
      _wetime[jtdc] = wetime[jtdc];
      _cpos[jtdc] = cpos[jtdc];
      _stepMC[jtdc] = stepMC[jtdc];
    }
  }


  double StrawDigiMC::driftDistance(StrawDigi::TDCChannel itdc) const {
    size_t jtdc = static_cast<size_t>(itdc);
    double retval = -100.0;
    if(!_stepMC[jtdc].isNull()){
      const Tracker& tracker = getTrackerOrThrow();
      // use the MC true index, not the straws index (digi could be from x-talk)
      Straw const& straw = tracker.getStraw(_stepMC[jtdc]->strawIndex());
      retval = (_cpos[itdc] - straw.getMidPoint()).perp(straw.getDirection());
    }
    return retval;
  }

  double StrawDigiMC::distanceToMid(StrawDigi::TDCChannel itdc) const {
    size_t jtdc = static_cast<size_t>(itdc);
    double retval = -100.0;
    if(!_stepMC[jtdc].isNull()){
      const Tracker& tracker = getTrackerOrThrow();
      Straw const& straw = tracker.getStraw(_stepMC[jtdc]->strawIndex());
      retval =  (_cpos[itdc] - straw.getMidPoint()).dot(straw.getDirection());
    }
    return retval;
  }

  bool StrawDigiMC::isCrossTalk(StrawDigi::TDCChannel itdc) const {
    bool retval(false);
    size_t jtdc = static_cast<size_t>(itdc);
    if(!_stepMC[jtdc].isNull()){
      retval = _strawIndex == _stepMC[jtdc]->strawIndex();
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


  double StrawDigiMC::triggerEnergySum(StrawDigi::TDCChannel itdc) const {
    double esum(0.0);
    if(!_stepMC[(size_t)itdc].isNull()){
      for(auto imcs = _stepMCs.begin(); imcs!= _stepMCs.end(); ++ imcs){
	// if the simParticle for this step is the same as the one which fired the discrim, add the energy
	if( (*imcs)->simParticle() == _stepMC[(size_t)itdc]->simParticle() ||
	    (*imcs)->simParticle()->parent() == _stepMC[(size_t)itdc]->simParticle() ||
	    (*imcs)->simParticle() == _stepMC[(size_t)itdc]->simParticle()->parent() )
	  esum += (*imcs)->eDep();
      }
    }
    return esum;
  }

  // Print the information found in this hit.
  void StrawDigiMC::print( ostream& ost, bool doEndl ) const {

    ost << "Straw Digi MC Truth:"
      << " cluster times : "      << _cpos[0].t() << " : " << _cpos[1].t()
      << " drift distance: "      << driftDistance(StrawDigi::zero) << " : " << driftDistance(StrawDigi::one)
      << " distance to wire center: "     << distanceToMid(StrawDigi::zero) << " : " << distanceToMid(StrawDigi::one)
      << " Energy: " << energySum();

    if ( doEndl ){
      ost << endl;
    }

  }

}

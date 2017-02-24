// G4 simulations may be run with proton time==const, then a pulse
// shape can be applied at the digitization stage using this class.
//
// Andrei Gaponenko, 2012


#include "ExtinctionMonitorFNAL/Digitization/inc/ProtonPulseShape.hh"

#include <iostream>
#include <cmath>
#include <map>
#include <utility>

#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindOne.h"

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/MARSInfo.hh"

#include "Mu2eUtilities/inc/SimParticleParentGetter.hh"

//#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    ProtonPulseShape::ProtonPulseShape(const fhicl::ParameterSet& pset,
                                       art::RandomNumberGenerator::base_engine_t& rng)
      : marsMode_(pset.get<bool>("marsMode"))
      , marsInfoModuleLabel_(marsMode_ ? pset.get<std::string>("marsInfoModuleLabel") : std::string(""))
      , marsInfoInstanceName_(pset.get<std::string>("marsInfoInstanceName", ""))
      , flat_(rng)
      , pulseHalfWidth_(pset.get<double>("pulseHalfWidth"))
      , messagePrinted_(false)
    {}

    //================================================================
    void ProtonPulseShape::apply(PixelChargeCollection *inout, const art::Event& event) {
      tmm_.clear();
      tgm_.clear();

      SimParticleParentGetter pg(event);

      if(!messagePrinted_) {
        messagePrinted_ = true;
        std::cout<<"Applying ProtonPulseShape in pixel digitization. marsMode="<<marsMode_
                 <<", pulseHalfWidth="<<pulseHalfWidth_<<std::endl;
      }

      for(PixelChargeCollection::iterator i=inout->begin(); i!=inout->end(); ++i) {
        apply(&i->second, pg, event);
      }
    }

    //================================================================
    void ProtonPulseShape::apply(PixelChargeHistory *inout,
                                 const SimParticleParentGetter& pg,
                                 const art::Event& event) {

      PixelChargeHistory& in(*inout);
      PixelChargeHistory out;
      while(!in.empty()) {
        PixelTimedChargeDeposit dep = in.top();

        dep.time += getTimeShiftForPrimary(dep.particle, pg, event);

        out.push(dep);
        in.pop();
      }

      *inout = out;
    }

    //================================================================
    double ProtonPulseShape::getTimeShiftForPrimary(const art::Ptr<SimParticle>& particle,
                                                    const SimParticleParentGetter& pg,
                                                    const art::Event& event
                                                    ) {
      double dt=0;

      art::Ptr<SimParticle> g4primary = getG4Primary(particle, pg);
      if(marsMode_) {
        ProtonPathMARSId id= getPrimaryMARSId(g4primary, event);
        TimeMapMARS::const_iterator it = tmm_.find(id);
        if(it == tmm_.end()) {
          // new primary. generate its offset
          dt = pulseHalfWidth_ * (2*flat_.fire() - 1);
          tmm_.insert(std::make_pair(id, dt));
        }
        else {
          dt = it->second;
        }
      } //----------------
      else { // Not marsMode_

        TimeMapG4::const_iterator it = tgm_.find(g4primary);
        if(it == tgm_.end()) {
          // new primary. generate its offset
          dt = pulseHalfWidth_ * (2*flat_.fire() - 1);
          tgm_.insert(std::make_pair(g4primary, dt));
        }
        else {
          dt = it->second;
        }
      }

      return dt;
    }

    //================================================================
    art::Ptr<SimParticle>
    ProtonPulseShape::getG4Primary(const art::Ptr<SimParticle>& p, const SimParticleParentGetter& pg) {
      art::Ptr<SimParticle> current(p);
      art::Ptr<SimParticle> next(pg.parent(current));
      while(next) {
        current = next;
        next = pg.parent(current);
      }
      return current;
    }

    //================================================================
    ProtonPulseShape::ProtonPathMARSId
    ProtonPulseShape::getPrimaryMARSId(const art::Ptr<SimParticle>& g4primary, const art::Event& event) {

      std::vector<art::Ptr<SimParticle> > particles;
      particles.push_back(g4primary);

      AGDEBUG("here");
      art::FindOne<MARSInfo> mif(particles, event, art::InputTag(marsInfoModuleLabel_, marsInfoInstanceName_));
      AGDEBUG("here");
      MARSInfo minfo = mif.at(0).ref();
      AGDEBUG("here");

      // FIXME: need to save g4s1info. But g4s2 jobs guarantee no
      // partly correlated particles in one event, so ignoring
      // g4s1info here should be OK.
      IO::G4JobInfo g4s1info;

      return ProtonPathMARSId(minfo, g4s1info);
    }

    //================================================================
  } // namespace ExtMonFNAL
} // namespace mu2

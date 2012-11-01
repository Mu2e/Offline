// G4 simulations may be run with proton time==const, then a pulse
// shape can be applied at the digitization stage using this class.
//
// Andrei Gaponenko, 2012


#include "ExtinctionMonitorFNAL/Digitization/inc/ProtonPulseShape.hh"

#include <iostream>
#include <cmath>
#include <map>
#include <utility>

#include "cetlib/exception.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Event.h"
#include "art/Utilities/InputTag.h"

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
    {
      if(!marsMode_) {
        throw cet::exception("NOTIMPLEMENTED")
          <<"ExtMonFNAL pixel digitization: ProtonPulseShape marsMode=false not imlemented\n";
      }
    }

    //================================================================
    void ProtonPulseShape::apply(PixelChargeCollection *inout, const art::Event& event) {
      tmm_.clear();

      SimParticleParentGetter pg(event);

      art::Handle<SimParticleCollection> simh;
      event.getByLabel(marsInfoModuleLabel_, marsInfoInstanceName_, simh);
      art::FindOne<MARSInfo> mif(simh, event, art::InputTag(marsInfoModuleLabel_, marsInfoInstanceName_));

      if(!messagePrinted_) {
        messagePrinted_ = true;
        std::cout<<"Applying ProtonPulseShape in pixel digitization. marsMode="<<marsMode_
                 <<", pulseHalfWidth="<<pulseHalfWidth_<<std::endl;
      }

      for(PixelChargeCollection::iterator i=inout->begin(); i!=inout->end(); ++i) {
        apply(&i->second, pg, mif);
      }
    }

    //================================================================
    void ProtonPulseShape::apply(PixelChargeHistory *inout,
                                 const SimParticleParentGetter& pg,
                                 const art::FindOne<MARSInfo>& mif
                                 ) {

      PixelChargeHistory& in(*inout);
      PixelChargeHistory out;
      while(!in.empty()) {
        PixelTimedChargeDeposit dep = in.top();

        dep.time += getTimeShiftForPrimary(dep.particle, pg, mif);

        out.push(dep);
        in.pop();
      }

      *inout = out;
    }

    //================================================================
    double ProtonPulseShape::getTimeShiftForPrimary(const art::Ptr<SimParticle>& particle,
                                                    const SimParticleParentGetter& pg,
                                                    const art::FindOne<MARSInfo>& mif
                                                    ) {
      double dt=0;

      PrimaryMARSId id= getPrimaryMARSId(particle, pg, mif);

      TimeMapMARS::const_iterator it = tmm_.find(id);
      if(it == tmm_.end()) {
        // new primary. generate its offset
        dt = pulseHalfWidth_ * (2*flat_.fire() - 1);
        tmm_.insert(std::make_pair(id, dt));
        AGDEBUG("new primary: dt = "<<dt);
      }
      else {
        dt = it->second;
        AGDEBUG("existing primary: dt = "<<dt);
      }

      return dt;
    }

    //================================================================
    ProtonPulseShape::PrimaryMARSId ProtonPulseShape::getPrimaryMARSId(const art::Ptr<SimParticle>& p,
                                                                       const SimParticleParentGetter& pg,
                                                                       const art::FindOne<MARSInfo>& mif
                                                                       ) {

      art::Ptr<SimParticle> current(p);
      art::Ptr<SimParticle> next(pg.parent(current));

      AGDEBUG("primary sim="<<current->id()<<", pdgId="<<current->pdgId());
      while(next) {
        current = next;
        next = pg.parent(current);
      }
      AGDEBUG("primary sim="<<current->id()<<", pdgId="<<current->pdgId());

      MARSInfo info = mif.at(current.key()).ref();

      return PrimaryMARSId(info.runNumber(), info.subRunNumber(), info.protonNumber());
    }

    //================================================================
  } // namespace ExtMonFNAL
} // namespace mu2e

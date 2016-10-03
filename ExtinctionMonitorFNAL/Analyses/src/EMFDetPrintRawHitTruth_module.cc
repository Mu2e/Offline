// Printout ExtMonFNAL raw hits and associated truth
//
// Andrei Gaponenko, 2012

#include "RecoDataProducts/inc/ExtMonFNALRawHit.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawHitCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/ExtMonFNALHitTruthAssn.hh"
#include "MCDataProducts/inc/MARSInfo.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"

#include <iostream>
#include <string>
#include <vector>

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    class EMFDetPrintRawHitTruth : public art::EDAnalyzer {
      std::string hitsModuleLabel_;
      std::string truthModuleLabel_;
      std::string simParticleMarsLabel_;
      std::string geomModuleLabel_;
      const ExtMon *extmon_;

      std::ostream& printParticle(std::ostream& os, const SimParticle& p);

      std::ostream& printParticleHistory(std::ostream& os,
                                         const art::Event& event,
                                         const art::Ptr<SimParticle>& p);
    public:
      explicit EMFDetPrintRawHitTruth(const fhicl::ParameterSet& pset);
      virtual void beginRun(const art::Run& run);
      virtual void analyze(const art::Event& event);
    };

    //================================================================
    EMFDetPrintRawHitTruth::EMFDetPrintRawHitTruth(const fhicl::ParameterSet& pset)
      : art::EDAnalyzer(pset)
      , hitsModuleLabel_(pset.get<std::string>("hitsModuleLabel"))
      , truthModuleLabel_(pset.get<std::string>("truthModuleLabel"))
      , simParticleMarsLabel_(pset.get<std::string>("simParticleMarsLabel"))
      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , extmon_()
    {}

    //================================================================
    void EMFDetPrintRawHitTruth::beginRun(const art::Run& run) {
      if(!geomModuleLabel_.empty()) {
        art::Handle<ExtMon> emf;
        run.getByLabel(geomModuleLabel_, emf);
        extmon_ = &*emf;
      }
      else {
        GeomHandle<ExtMon> emf;
        extmon_ = &*emf;
      }
    }

    //================================================================
    void EMFDetPrintRawHitTruth::analyze(const art::Event& event) {

      art::Handle<ExtMonFNALRawHitCollection> ih;
      event.getByLabel(hitsModuleLabel_, ih);

      const ExtMonFNALRawHitCollection& inputs(*ih);

      art::FindManyP<SimParticle,ExtMonFNALHitTruthBits> r2t(ih, event, truthModuleLabel_);

      std::cout<<"EMFDetPrintRawHitTruth: hitsModuleLabel = "<<hitsModuleLabel_<<", truthModuleLabel = "<<truthModuleLabel_<<std::endl;

      for(ExtMonFNALRawHitCollection::const_iterator i=inputs.begin(); i!=inputs.end(); ++i) {
        std::cout<<"event "<<event.id()<<", hit "<<*i<<std::endl;

        //if(i->clock() == 43) {
        if(true) {
          std::vector<art::Ptr<SimParticle> > particles;
          std::vector<const ExtMonFNALHitTruthBits*> charges;

          r2t.get((i-inputs.begin()), particles, charges);
          std::cout<<"got truth: num particles = "<<particles.size()<<", num charges "<<charges.size()<<std::endl;;
          for(unsigned ip=0; ip<particles.size(); ++ip) {
            printParticleHistory(std::cout, event, particles[ip]);
          }
        }
      }
    }

    //================================================================
    std::ostream& EMFDetPrintRawHitTruth::printParticleHistory(std::ostream& os,
                                                               const art::Event& event,
                                                               const art::Ptr<SimParticle>& p)
    {
      art::Ptr<SimParticle> pp = p;

      printParticle(os, *pp); os<<"\n";

      while(pp->parent()) {
        pp = pp->parent();
        os<<"    parent: "; printParticle(os, *pp); os<<"\n";
      }

      art::Ptr<GenParticle> gp = pp->genParticle();
      if(gp.isAvailable()) {
        os<<"    GenParticle: "<<*gp<<"\n";
      }
      else {
        os<<"    GenParticle not available\n";
      }

      if(!simParticleMarsLabel_.empty()) {
        std::vector<art::Ptr<SimParticle> > particles;
        particles.push_back(pp);
        art::FindOneP<MARSInfo> mif(particles, event, simParticleMarsLabel_);
        art::Ptr<MARSInfo> pm = mif.at(0);
        if(pm) {
          os<<"    MARSInfo for SimParticle:"<<*pm<<"\n";
        }
        else {
          os<<"    MARSInfo for SimParticle not available\n";
        }
      }

      return os;
    }

    //================================================================
    std::ostream& EMFDetPrintRawHitTruth::printParticle(std::ostream& os, const SimParticle& p) {
      os<<"SimParticle(id="<<p.id()
        <<", pdgId="<<p.pdgId()
        <<", startMomExtMon="<<extmon_->mu2eToExtMon_momentum(p.startMomentum())
        <<", startPosExtMon="<<extmon_->mu2eToExtMon_position(p.startPosition())
        <<", startGlobalTime="<<p.startGlobalTime()
        <<", endPosExtMon="<<extmon_->mu2eToExtMon_position(p.endPosition())
        <<", endGlobalTime="<<p.endGlobalTime()
        <<")";

      return os;
    }

    //================================================================
  } // namespace ExtMonFNAL
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::EMFDetPrintRawHitTruth);

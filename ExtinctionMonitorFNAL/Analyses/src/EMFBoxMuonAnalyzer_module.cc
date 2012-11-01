// $Id: EMFBoxMuonAnalyzer_module.cc,v 1.1 2012/11/01 23:35:06 gandr Exp $
// $Author: gandr $
// $Date: 2012/11/01 23:35:06 $
//
// Original author Andrei Gaponenko, 2012

#include <iostream>
#include <string>
#include <cmath>

#include "cetlib/exception.h"

#include "CLHEP/Units/PhysicalConstants.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindOne.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "GeometryService/inc/GeomHandle.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/MARSInfo.hh"
#include "MCDataProducts/inc/MARSInfoCollection.hh"
#include "MCDataProducts/inc/GenParticleMARSAssns.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    namespace {
      struct StoppedMuon {
        double emx;
        double emy;
        double emz;
        double time;
        double endek;

        int pdgId;

        unsigned endG4Status;
        unsigned stoppingCode;

        StoppedMuon() : emx(), emy(), emz(), time(), endek(), pdgId(), endG4Status(), stoppingCode() {}

        StoppedMuon(const ExtMon& extmon, const SimParticle& sp)
          : emx()
          , emy()
          , emz()
          , time(sp.endGlobalTime())
          , endek(sp.endMomentum().e() - sp.endMomentum().m())
          , pdgId(sp.pdgId())
          , endG4Status(sp.endG4Status())
          , stoppingCode(sp.stoppingCode().id())
        {
          const CLHEP::Hep3Vector pos = extmon.mu2eToExtMon_position(sp.endPosition());
          emx = pos.x();
          emy = pos.y();
          emz = pos.z();
        }
      };

    } // namespace {}

    //================================================================
    class EMFBoxMuonAnalyzer : public art::EDAnalyzer {
      std::string g4ModuleLabel_;
      std::string g4InstanceName_;

      std::string generatorModuleLabel_;
      std::string generatorInstanceName_;

      std::string geomModuleLabel_;
      std::string geomInstanceName_;

      double cutEKineAtStop_;
      std::vector<double> cutExtMonPos_;

      const ExtMon *extmon_;

      StoppedMuon sm_;
      MARSInfo minfo_;
      TTree *nt_;

      bool isAccepted(const StoppedMuon& sm);

    public:
      explicit EMFBoxMuonAnalyzer(const fhicl::ParameterSet& pset);
      virtual void analyze(const art::Event& event);
      virtual void beginJob();
      virtual void beginRun(const art::Run& run);
    };

    //================================================================
    EMFBoxMuonAnalyzer::EMFBoxMuonAnalyzer(const fhicl::ParameterSet& pset)
      : g4ModuleLabel_(pset.get<std::string>("g4ModuleLabel"))
      , g4InstanceName_(pset.get<std::string>("g4InstanceName", ""))

      , generatorModuleLabel_(pset.get<std::string>("generatorModuleLabel"))
      , generatorInstanceName_(pset.get<std::string>("generatorInstanceName", ""))

      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))

      , cutEKineAtStop_(pset.get<double>("cutEKineAtStop"))
      , cutExtMonPos_(pset.get<std::vector<double> >("cutExtMonPos"))

      , extmon_()
      , nt_()
    {}

    //================================================================
    void EMFBoxMuonAnalyzer::beginJob() {
      art::ServiceHandle<art::TFileService> tfs;
      nt_ = tfs->make<TTree>( "sm", "Stopped muons ntuple");
      nt_->Branch("particle", &sm_, "emx/D:emy/D:emz/D:time/D:endek/D:pdgId/I:endG4Status/i:stoppingCode/i");
      nt_->Branch("minfo", &minfo_, "weight/D:protonNumber/I:subRunNumber/I");
    }

    //================================================================
    void EMFBoxMuonAnalyzer::beginRun(const art::Run& run) {
      if(!geomModuleLabel_.empty()) {
        art::Handle<ExtMon> extmon;
        run.getByLabel(geomModuleLabel_, geomInstanceName_, extmon);
        extmon_ = &*extmon;
      }
      else {
        GeomHandle<ExtMon> extmon;
        extmon_ = &*extmon;
      }
    }

    //================================================================
    void EMFBoxMuonAnalyzer::analyze(const art::Event& event) {

      art::Handle<SimParticleCollection> ih;
      event.getByLabel(g4ModuleLabel_, g4InstanceName_, ih);
      const SimParticleCollection& particles(*ih);

      art::Handle<GenParticleCollection> genh;
      event.getByLabel(generatorModuleLabel_, generatorInstanceName_, genh);
      art::FindOne<MARSInfo>
        mFinder(genh, event, art::InputTag(generatorModuleLabel_, generatorInstanceName_));

      for(SimParticleCollection::const_iterator i=particles.begin(); i!=particles.end(); ++i) {
        const SimParticle& sp = i->second;
        if(std::abs(sp.pdgId())==13) {
          const double eKine = sp.endMomentum().e() - sp.endMomentum().m();
          if(eKine < cutEKineAtStop_) {

            // compute particle kinematics
            sm_ = StoppedMuon(*extmon_, sp);

            // Only record stops in the vicinity of ExtMonFNAL
            if(isAccepted(sm_)) {

              // compute MARS data
              art::Ptr<GenParticle> gen = sp.genParticle();
              if(sp.parent()) {
                art::Ptr<SimParticle> top(sp.parent());
                while(top->parent()) {
                  top = top->parent();
                }
                gen = top->genParticle();
              }
              if(!gen) {
                throw cet::exception("BADINPUTS")<<"ERROR: no GenParticle for SimParticle "
                                                 <<sp.id()<<" in event "<<event.id()<<"\n";
              }
              minfo_ = mFinder.at(gen.key()).ref();

              // write
              nt_->Fill();
            }

          }
        }
      } // for()
    } // analyze()

    //================================================================
    bool EMFBoxMuonAnalyzer::isAccepted(const StoppedMuon& sm) {
      return
        (std::abs(sm.emx) < cutExtMonPos_[0]) &&
        (std::abs(sm.emy) < cutExtMonPos_[1]) &&
        (std::abs(sm.emz) < cutExtMonPos_[2]);
    }

    //================================================================
  } // namespace ExtMonFNAL
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::EMFBoxMuonAnalyzer);

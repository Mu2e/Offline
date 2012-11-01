// $Id: EMFBoxMuonAnalyzer_module.cc,v 1.4 2012/11/01 23:41:17 gandr Exp $
// $Author: gandr $
// $Date: 2012/11/01 23:41:17 $
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

#include "ExtinctionMonitorFNAL/Utilities/inc/EMFBoxIO.hh"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    using IO::StoppedMuon;

    namespace {

      StoppedMuon makeStoppedMuon(const ExtMon& extmon, const SimParticle& sp) {

        StoppedMuon res;

        const CLHEP::Hep3Vector pos = extmon.mu2eToExtMon_position(sp.endPosition());
        res.emx = pos.x();
        res.emy = pos.y();
        res.emz = pos.z();

        res.time = sp.endGlobalTime();
        res.endek = sp.endMomentum().e() - sp.endMomentum().m();
        res.pdgId = sp.pdgId();

        res.endG4Status = sp.endG4Status();
        res.stoppingCode = sp.stoppingCode().id();

        return res;
      }

    } // namespace {}

    //================================================================
    class EMFBoxMuonAnalyzer : public art::EDAnalyzer {
      std::string particlesModuleLabel_;
      std::string particlesInstanceName_;

      std::string marsInfoModuleLabel_;
      std::string marsInfoInstanceName_;

      std::string geomModuleLabel_;
      std::string geomInstanceName_;

      double cutEKineAtStop_;
      std::vector<double> cutExtMonPos_;

      const ExtMon *extmon_;

      StoppedMuon sm_;
      IO::MARSInfo minfo_;
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
      : particlesModuleLabel_(pset.get<std::string>("particlesModuleLabel"))
      , particlesInstanceName_(pset.get<std::string>("particlesInstanceName", ""))

      , marsInfoModuleLabel_(pset.get<std::string>("marsInfoModuleLabel"))
      , marsInfoInstanceName_(pset.get<std::string>("marsInfoInstanceName", ""))

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
      nt_->Branch("particle", &sm_, sm_.branchDescription());
      nt_->Branch("minfo", &minfo_, minfo_.branchDescription());
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
      event.getByLabel(particlesModuleLabel_, particlesInstanceName_, ih);
      const SimParticleCollection& particles(*ih);

      for(SimParticleCollection::const_iterator i=particles.begin(); i!=particles.end(); ++i) {
        const SimParticle& sp = i->second;
        if(std::abs(sp.pdgId())==13) {
          const double eKine = sp.endMomentum().e() - sp.endMomentum().m();
          if(eKine < cutEKineAtStop_) {

            // compute particle kinematics
            sm_ = makeStoppedMuon(*extmon_, sp);

            // Only record stops in the vicinity of ExtMonFNAL
            if(isAccepted(sm_)) {

              art::Ptr<SimParticle> top(ih, i->first.asUint());
              while(top->parent()) {
                top = top->parent();
              }

              // Found primary SimParticle.  It should have associated MARSInfo
              std::vector<art::Ptr<SimParticle> > vsp;
              vsp.push_back(top);

              art::FindOne<MARSInfo>
                mFinder(vsp, event,
                        art::InputTag(marsInfoModuleLabel_, marsInfoInstanceName_));

              minfo_.info = mFinder.at(0).ref();

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

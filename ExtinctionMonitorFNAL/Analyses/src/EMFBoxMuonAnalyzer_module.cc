//
// Original author Andrei Gaponenko, 2012

#include <iostream>
#include <string>
#include <cmath>

#include "cetlib_except/exception.h"

#include "CLHEP/Units/PhysicalConstants.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

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

        const ExtMonFNALPlaneStack& stack = (pos.z() < 0) ? extmon.dn() : extmon.up();
        const CLHEP::Hep3Vector stackPos = stack.mu2eToStack_position(sp.endPosition());

        // Note: we don't cut on the z position.  A muon stopped
        // within the "sensor size" from the axis will be smeared
        // in the transverse direction, even if it stopped in air
        // instead of silicon.
        res.stoppedInSensor =
          (std::abs(stackPos.x()) < stack.planes()[0].halfSize()[0]) &&
          (std::abs(stackPos.y()) < stack.planes()[0].halfSize()[1])
          ;

        return res;
      }

      //----------------------------------------------------------------
      struct BufferEntry {
        StoppedMuon sm;
        MARSInfo  minfo;
        BufferEntry(const StoppedMuon& s, const MARSInfo& m)
          : sm(s), minfo(m)
        {}
      };
      struct ProtonSort {
        bool operator()(const BufferEntry& a, const BufferEntry& b) {
          // We don't have to look at g4s1info because this module
          // processes one event at a time.
          CmpProtonIdAndSimPath cm;
          return cm(a.minfo,b.minfo);
        }
      };

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
      IO::G4JobInfo g4s1info_;
      TTree *nt_;

      unsigned numMuonStops_;
      unsigned numProtonPaths_;
      std::set<MARSInfo, CmpProtonId> uniqProtons_;

      bool isAccepted(const StoppedMuon& sm);

    public:
      explicit EMFBoxMuonAnalyzer(const fhicl::ParameterSet& pset);
      virtual void analyze(const art::Event& event);
      virtual void beginJob();
      virtual void beginRun(const art::Run& run);
      virtual void endJob();
    };

    //================================================================
    EMFBoxMuonAnalyzer::EMFBoxMuonAnalyzer(const fhicl::ParameterSet& pset)
      : art::EDAnalyzer(pset)
      , particlesModuleLabel_(pset.get<std::string>("particlesModuleLabel"))
      , particlesInstanceName_(pset.get<std::string>("particlesInstanceName", ""))

      , marsInfoModuleLabel_(pset.get<std::string>("marsInfoModuleLabel"))
      , marsInfoInstanceName_(pset.get<std::string>("marsInfoInstanceName", ""))

      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))

      , cutEKineAtStop_(pset.get<double>("cutEKineAtStop"))
      , cutExtMonPos_(pset.get<std::vector<double> >("cutExtMonPos"))

      , extmon_()
      , nt_()

      , numMuonStops_(0)
      , numProtonPaths_(0)
    {}

    //================================================================
    void EMFBoxMuonAnalyzer::beginJob() {
      art::ServiceHandle<art::TFileService> tfs;
      nt_ = tfs->make<TTree>( "sm", "Stopped muons ntuple");
      nt_->Branch("particle", &sm_, sm_.branchDescription());
      nt_->Branch("minfo", &minfo_, minfo_.branchDescription());
      nt_->Branch("g4s1info", &g4s1info_, g4s1info_.branchDescription());
    }

    //================================================================
    void EMFBoxMuonAnalyzer::endJob() {
      std::cout<<"EMFBoxMuonAnalyzer: numProtonPaths = "<<numProtonPaths_
               <<", numUniqProtons = "<<uniqProtons_.size()
               <<", numMuonStops = "<<numMuonStops_
               <<std::endl;
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

      // We want to make sure particles from a single proton are written
      // out as a contiguous chunck.  Need to accumulate all inputs from
      // one event, sort by proton, then write out in order.
      std::vector<BufferEntry> buf;

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

              // will write this out
              buf.push_back(BufferEntry(sm_, minfo_.info));
            } // isAccepted()

          }
        }
      } // for()

      //----------------
      // Got all the muons in the buffer.  Need to sort by proton.
      ProtonSort ps;
      std::sort(buf.begin(), buf.end(), ps);

      // write out, and count protons in the process
      numMuonStops_ += buf.size();
      if(!buf.empty()) {

        MARSInfo current = buf[0].minfo;
        ++numProtonPaths_;
        uniqProtons_.insert(current);

        for(unsigned i=0; i<buf.size(); ++i) {

          sm_ = buf[i].sm;

          minfo_.info = buf[i].minfo;

          g4s1info_.run = event.id().run();
          g4s1info_.subrun = event.id().subRun();
          g4s1info_.event = event.id().event();

          nt_->Fill();

          if(!sameProtonAndSimPath(current, buf[i].minfo)) {
            current = buf[i].minfo;
            ++numProtonPaths_;
            uniqProtons_.insert(current);
          }
        }
      }

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

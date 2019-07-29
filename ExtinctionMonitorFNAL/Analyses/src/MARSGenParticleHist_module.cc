// Andrei Gaponenko, 2012

#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>

#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    class MARSGenParticleHist : public art::EDAnalyzer {
      // Geometry for coordinate conversion
      std::string geomModuleLabel_;

      // GenParticleCollection to operate on
      std::string inputModuleLabel_;
      std::string inputInstanceName_;

      // Coordinate conversion
      const ProtonBeamDump *dump_;

      // Output histograms
      TH2* hzx_;
      TH2* hxy_;
      TH2* hyz_;

    public:
      explicit MARSGenParticleHist(const fhicl::ParameterSet& pset);
      virtual void beginJob();
      virtual void beginRun(const art::Run& run);
      virtual void analyze(const art::Event& event);
    };

    //================================================================
    MARSGenParticleHist::MARSGenParticleHist(const fhicl::ParameterSet& pset)
      : art::EDAnalyzer(pset)
      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel", ""))
      , inputModuleLabel_(pset.get<std::string>("inputModuleLabel"))
      , inputInstanceName_(pset.get<std::string>("inputInstanceName"))
      , dump_(0)
      , hzx_()
      , hxy_()
      , hyz_()
    {}

    //================================================================
    void MARSGenParticleHist::beginJob() {
      art::ServiceHandle<art::TFileService> tfs;
      hzx_ = tfs->make<TH2D>("dumpzx", "dump x vs z", 40, -10400., -7000., 40, -2200., 3200.);
      hxy_ = tfs->make<TH2D>("dumpxy", "dump y vs x", 40, -2200., 3200., 40, 1500., 4700.);
      hyz_ = tfs->make<TH2D>("dumpyz", "dump y vs z",  40, -10400., -7000., 40, 1500., 4700.);
    }

    //================================================================
    void MARSGenParticleHist::beginRun(const art::Run& run) {
      if(geomModuleLabel_.empty()) {
        dump_ = &*GeomHandle<ProtonBeamDump>();
      }
      else {
        // Here we retrieve geometry objects from the
        // input file, as opposed to using GeometryService that might have
        // a different version than used to produce our input data.
        art::Handle<ProtonBeamDump> dumph;
        run.getByLabel(geomModuleLabel_, "", dumph);
        dump_ = &*dumph;
      }
    }

    //================================================================
    void MARSGenParticleHist::analyze(const art::Event& event) {

      art::Handle<GenParticleCollection> particlesh;
      event.getByLabel(inputModuleLabel_, inputInstanceName_, particlesh);

      const GenParticleCollection& particles(*particlesh);
      for(unsigned i=0; i < particles.size(); ++i) {

        CLHEP::Hep3Vector posDump = dump_->mu2eToBeamDump_position(particles[i].position());
        hzx_->Fill(posDump.z(), posDump.x());
        hxy_->Fill(posDump.x(), posDump.y());
        hyz_->Fill(posDump.z(), posDump.y());

      } // for()
    } // analyze(event)

    //================================================================

  } // namespace ExtMonFNAL
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::MARSGenParticleHist);

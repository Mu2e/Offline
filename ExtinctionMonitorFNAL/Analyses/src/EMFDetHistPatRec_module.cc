// ExtMonFNAL PatRec efficiency/fake rate analysis
//
// Andrei Gaponenko, 2012

#include <iostream>
#include <string>
#include <vector>
#include <set>

#include "CLHEP/GenericFunctions/CumulativeChiSquare.hh"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "RecoDataProducts/inc/ExtMonFNALTrkParam.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFitQuality.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkClusterResiduals.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFit.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFitCollection.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TH1D.h"
#include "TH2D.h"


#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    class EMFDetHistPatRec : public art::EDAnalyzer {
      std::string patRecModuleLabel_;
      std::string patRecInstanceName_;
      std::string geomModuleLabel_; // emtpy to take info from Run
      std::string geomInstanceName_;

      const ExtMon *extmon_;

      //----------------
      TH1D *hMultiplicity_;
      TH1D *hChi2_;
      TH1D *hProb_;

      //      std::vector<TH2D*> hResiduals;

    public:
      explicit EMFDetHistPatRec(const fhicl::ParameterSet& pset);
      virtual void beginRun(const art::Run& run);
      virtual void analyze(const art::Event& event);
    };

    //================================================================
    EMFDetHistPatRec::EMFDetHistPatRec(const fhicl::ParameterSet& pset)
      : patRecModuleLabel_(pset.get<std::string>("patRecModuleLabel"))
      , patRecInstanceName_(pset.get<std::string>("patRecInstanceName", ""))
      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))

      , extmon_()

      , hMultiplicity_()
      , hChi2_()
      , hProb_()
    {}

    //================================================================
    void EMFDetHistPatRec::beginRun(const art::Run& run) {
      // This is a workaround for geometry not being available at beginJob()
      if(!extmon_) {
        if(!geomModuleLabel_.empty()) {
          art::Handle<ExtMonFNAL::ExtMon> emf;
          run.getByLabel(geomModuleLabel_, geomInstanceName_, emf);
          extmon_ = &*emf;
        }
        else {
          GeomHandle<ExtMonFNAL::ExtMon> emf;
          extmon_ = &*emf;
        }

        //----------------------------------------------------------------
        art::ServiceHandle<art::TFileService> tfs;

        hMultiplicity_ = tfs->make<TH1D>("multiplicity", "Num PatRec tracks", 200, -0.5, 199.5);
        hMultiplicity_->StatOverflows();
        hMultiplicity_->GetXaxis()->SetTitle("num PatRec tracks");

        hChi2_ = tfs->make<TH1D>("chi2", "PatRec chi2", 500, 0., 100.);
        hChi2_->StatOverflows();

        hProb_ = tfs->make<TH1D>("prob", "PatRec prob(chi2,ndf)", 500, 0., 1.);
        hProb_->StatOverflows();
      }
    }

    //================================================================
    void EMFDetHistPatRec::analyze(const art::Event& event) {

      art::Handle<ExtMonFNALTrkFitCollection> tracksh;
      event.getByLabel(patRecModuleLabel_, patRecInstanceName_, tracksh);
      const ExtMonFNALTrkFitCollection& tracks(*tracksh);

      hMultiplicity_->Fill(tracks.size());

      for(unsigned i=0; i<tracks.size(); ++i) {

        hChi2_->Fill(tracks[i].quality().chi2());

        Genfun::CumulativeChiSquare pf(tracks[i].quality().ndf());
        const double prob = 1. - pf(tracks[i].quality().chi2());
        hProb_->Fill(prob);

      } // for(i)

    } // analyze()

    //================================================================
  } // namespace ExtMonFNAL
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::EMFDetHistPatRec);

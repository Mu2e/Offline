// ExtMonFNAL PatRec efficiency/fake rate analysis
//
// Andrei Gaponenko, 2012

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "CLHEP/GenericFunctions/CumulativeChiSquare.hh"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/InputTag.h"

#include "RecoDataProducts/inc/ExtMonFNALTrkParam.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFitQuality.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkClusterResiduals.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFit.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFitCollection.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TH1D.h"
#include "TH2D.h"


#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    class EMFDetHistPatRec : public art::EDAnalyzer {
      art::InputTag patRecInTag_;
      std::string geomModuleLabel_; // emtpy to take info from Run
      std::string geomInstanceName_;

      const ExtMon *extmon_;

      //----------------
      const double maxMomentum_; // for histogramming

      TH1D *hMultiplicity_;
      TH1D *hChi2_;
      TH1D *hProb_;
      TH1D *hMomentum_;
      TH1D *hpxpz_;
      TH1D *hpypz_;
      TH1D *hAngle_;

      //      std::vector<TH2D*> hResiduals;

    public:
      explicit EMFDetHistPatRec(const fhicl::ParameterSet& pset);
      virtual void beginRun(const art::Run& run);
      virtual void analyze(const art::Event& event);
    };

    //================================================================
    EMFDetHistPatRec::EMFDetHistPatRec(const fhicl::ParameterSet& pset)
      : art::EDAnalyzer(pset)
      , patRecInTag_(pset.get<std::string>("patRecInTag"))
      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))

      , extmon_()

      , maxMomentum_(8000.) /* MeV/c */

      , hMultiplicity_()
      , hChi2_()
      , hProb_()
      , hMomentum_()
      , hpxpz_()
      , hpypz_()
      , hAngle_()
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

        hMomentum_ = tfs->make<TH1D>("momentum", "Reco track momentum", 160, 0., maxMomentum_);
        hMomentum_->GetXaxis()->SetTitle("p, MeV/c");

        hpxpz_ = tfs->make<TH1D>("pxpz", "px/pz", 200, -0.05, 0.05);
        hpxpz_->GetXaxis()->SetTitle("px/pz");

        hpypz_ = tfs->make<TH1D>("pypz", "py/pz", 200, -0.05, 0.05);
        hpypz_->GetXaxis()->SetTitle("py/pz");

        hAngle_ = tfs->make<TH1D>("angle", "Track angle w.r.t. the axis", 100, 0, 0.025);
        hAngle_->GetXaxis()->SetTitle("theta, rad");
      }
    }

    //================================================================
    void EMFDetHistPatRec::analyze(const art::Event& event) {

      auto tracksh = event.getValidHandle<ExtMonFNALTrkFitCollection>(patRecInTag_);

      hMultiplicity_->Fill(tracksh->size());

      for(const auto& track: *tracksh) {
        hChi2_->Fill(track.quality().chi2());

        Genfun::CumulativeChiSquare pf(track.quality().ndf());
        const double prob = 1. - pf(track.quality().chi2());
        hProb_->Fill(prob);

        const double pinv = extmon_->spectrometerMagnet().trackPinvFromRinv(track.pars().rinv());
        hMomentum_->Fill( (pinv > 1./maxMomentum_) ? 1./pinv : 1./maxMomentum_);

        hpxpz_->Fill(track.pars().slopex());
        hpypz_->Fill(track.pars().slopey());

        const double perp =
          std::sqrt(std::pow(track.pars().slopex(), 2) + std::pow(track.pars().slopey(), 2));
        const double theta = std::atan(perp);
        hAngle_->Fill(theta);

      } // for(i)

    } // analyze()

    //================================================================
  } // namespace ExtMonFNAL
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::EMFDetHistPatRec);

// Andrei Gaponenko, 2013

// C++ includes
#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

// cetlib includes
#include "cetlib_except/exception.h"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

// Mu2e includes
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/MuonCaptureSpectrum.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Mu2eUtilities/inc/Table.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"

namespace mu2e {

  //================================================================
  class RMCGun : public art::EDProducer {
    fhicl::ParameterSet psphys_;

    BinnedSpectrum spectrum_;
    static BinnedSpectrum parseSpectrumShape(const fhicl::ParameterSet& psphys,
                                             double *elow,
                                             double *ehi);

    int                 verbosityLevel_;
    int                 generateInternalConversion_;

    double              czmin_;
    double              czmax_;
    double              phimin_;
    double              phimax_;

    art::RandomNumberGenerator::base_engine_t& eng_;

    CLHEP::RandGeneral  randSpectrum_;
    CLHEP::RandFlat     randomFlat_;
    RandomUnitSphere    randomUnitSphere_;
    MuonCaptureSpectrum muonCaptureSpectrum_;

    RootTreeSampler<IO::StoppedParticleF> stops_;

    bool doHistograms_;

    double generateEnergy();

    TH1F* _hmomentum;
    TH1F* _hElecMom {nullptr};
    TH1F* _hPosiMom {nullptr};
    TH1F* _hTotMom {nullptr};
    TH1F* _hMee;
    TH2F* _hMeeVsE;
    TH1F* _hMeeOverE;                   // M(ee)/E(gamma)
    TH1F* _hy;				// splitting function

  public:
    explicit RMCGun(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event);
  };

  //================================================================
  RMCGun::RMCGun(const fhicl::ParameterSet& pset)
    : EDProducer{pset}
    , psphys_(pset.get<fhicl::ParameterSet>("physics"))
    , spectrum_                  (BinnedSpectrum(psphys_))
    , verbosityLevel_            (pset.get<int>   ("verbosityLevel", 0))
    , generateInternalConversion_{psphys_.get<int>("generateIntConversion", 0)}
    , czmin_                     (pset.get<double>("czmin" , -1.0))
    , czmax_                     (pset.get<double>("czmax" ,  1.0))
    , phimin_                    (pset.get<double>("phimin",  0. ))
    , phimax_                    (pset.get<double>("phimax", CLHEP::twopi ))
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randSpectrum_       (eng_, spectrum_.getPDF(), spectrum_.getNbins())
    , randomFlat_         (eng_)
    , randomUnitSphere_   (eng_, czmin_,czmax_,phimin_,phimax_)
    , muonCaptureSpectrum_(&randomFlat_,&randomUnitSphere_)
      //    , randomUnitSphere_(eng_)
    , stops_(eng_, pset.get<fhicl::ParameterSet>("muonStops"))
    , doHistograms_( pset.get<bool>("doHistograms",true ) )
  {
    produces<mu2e::GenParticleCollection>();

    if(verbosityLevel_ > 0) {
      std::cout<<"RMCGun: using = "
               <<stops_.numRecords()
               <<" stopped particles"
               <<std::endl;

      std::cout<<"RMCGun: producing photon " << std::endl;
    }

    if ( doHistograms_ ) {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "RMCGun" );

      _hmomentum     = tfdir.make<TH1F>( "hmomentum", "Produced photon momentum", 100,  40.,  140.  );

      if(generateInternalConversion_){
        _hElecMom  = tfdir.make<TH1F>("hElecMom" , "Produced electron momentum", 140,  0. , 140.);
        _hPosiMom  = tfdir.make<TH1F>("hPosiMom" , "Produced positron momentum", 140,  0. , 140.);
        _hTotMom   = tfdir.make<TH1F>("hTotMom" , "Produced total momentum", 100,  40. , 140.);
        _hMee      = tfdir.make<TH1F>("hMee"     , "M(e+e-) "           , 200,0.,200.);
        _hMeeVsE   = tfdir.make<TH2F>("hMeeVsE"  , "M(e+e-) vs E"       , 200,0.,200.,200,0,200);
        _hMeeOverE = tfdir.make<TH1F>("hMeeOverE", "M(e+e-)/E "         , 200, 0.,1);
        _hy        = tfdir.make<TH1F>("hy"       , "y = (ee-ep)/|pe+pp|", 200,-1.,1.);
      }
    }

  }

  //================================================================
  BinnedSpectrum
  RMCGun::parseSpectrumShape(const fhicl::ParameterSet& psphys,
                                                 double *elow,
                                                 double *ehi)
  {
    BinnedSpectrum res;

    const std::string spectrumShape(psphys.get<std::string>("spectrumShape"));
    if (spectrumShape == "RMC") {
      *elow = psphys.get<double>("elow");
      *ehi = psphys.get<double>("ehi");
      res.initialize<MuonCaptureSpectrum>( *elow, *ehi, psphys.get<double>("spectrumResolution") );
    }
    else if (spectrumShape == "flat") {
      *elow = psphys.get<double>("elow");
      *ehi = psphys.get<double>("ehi");
      res.initialize<SimpleSpectrum>(*elow, *ehi, *ehi-*elow, SimpleSpectrum::Spectrum::Flat );
    }
    else {
      throw cet::exception("BADCONFIG")
        << "StoppedParticleMuonGun: unknown spectrum shape "<<spectrumShape<<"\n";
    }

    return res;
  }

  //================================================================
  void RMCGun::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    const auto& stop = stops_.fire();

    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

    const double energy = generateEnergy();

    if(!generateInternalConversion_){
      output->emplace_back( PDGCode::gamma,
                          GenId::ExternalRMC,
                          pos,
                          CLHEP::HepLorentzVector( randomUnitSphere_.fire(energy), energy),
                          stop.t );

      event.put(std::move(output));
    } else {
      CLHEP::HepLorentzVector mome, momp;
      muonCaptureSpectrum_.getElecPosiVectors(energy,mome,momp);
      // Add particles to list
      auto output = std::make_unique<GenParticleCollection>(); //GenID = 42
      output->emplace_back(PDGCode::e_minus, GenId::InternalRMC,pos,mome,stop.t);
      output->emplace_back(PDGCode::e_plus , GenId::InternalRMC,pos,momp,stop.t);
      event.put(move(output));

      if(doHistograms_){
        _hElecMom ->Fill(mome.vect().mag());
        _hPosiMom ->Fill(momp.vect().mag());
        _hTotMom ->Fill(mome.vect().mag()+momp.vect().mag());

        double mee = (mome+momp).m();
        _hMee->Fill(mee);
        _hMeeVsE->Fill(energy,mee);
        _hMeeOverE->Fill(mee/energy);

        CLHEP::Hep3Vector p = mome.vect()+momp.vect();
        double y = (mome.e()-momp.e())/p.mag();

        _hy->Fill(y);
      }
    }

    if ( !doHistograms_ ) return;

    _hmomentum->Fill(energy);

  }

  //================================================================
  double RMCGun::generateEnergy() {
    return spectrum_.sample(randSpectrum_.fire());
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::RMCGun);

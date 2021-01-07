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

#include "fhiclcpp/types/OptionalAtom.h"

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
#include "TF1.h"

namespace mu2e {

  //================================================================
  class RMCGun : public art::EDProducer {
  public:
    typedef RootTreeSampler<IO::StoppedParticleF> RTS;
    
    struct PhysConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      //      fhicl::Atom<int> generateIntConversion{Name("generateIntConversion"), Comment("Generate internal conversion (>=0)"), 0};
      fhicl::Atom<std::string> genId{Name("genId"), Comment("gen process ID: InternalRMC/ExternalRMC")};
      fhicl::Atom<double> elow{Name("elow"), Comment("Minimum energy value (MeV)")};
      fhicl::Atom<double> ehi{Name("ehi"), Comment("Maximum energy value (MeV)")};
      fhicl::Atom<std::string> spectrumShape{Name("spectrumShape"), Comment("Spectrum shape (e.g. RMC or flat)")};
      fhicl::Atom<double> spectrumResolution{Name("spectrumResolution"), Comment("Energy resolution for binned spectrum (MeV)")};
      fhicl::OptionalAtom<double> kMaxUser{Name("kMaxUser"), Comment("kMax value to use (MeV)")};
      fhicl::OptionalAtom<bool> kMaxUserSet{Name("kMaxUserSet"), Comment("Whether or not to use a user kMax value (true/false)")};
      //to allow for the other binned spectrum information that could be provided
      fhicl::OptionalAtom<bool> fixMax{Name("FixMax"), Comment("Fix max (true/false)")};      
      fhicl::OptionalAtom<PDGCode::type> pdgId{Name("pdgId"), Comment("PDG ID")};      
      fhicl::OptionalAtom<unsigned> nbins{Name("nbins"), Comment("Number of bins")};      
      fhicl::OptionalAtom<std::string> spectrumFileName{Name("spectrumFileName"), Comment("Name of spectrum file")};      
      fhicl::OptionalAtom<bool> binCenter{Name("binCenter"), Comment("If spectrum uses bin centers")};      
    };
    
      
    //ignore assumed keys for a EDProducer
    struct KeysToIgnore {
      std::set<std::string> operator()()
      {
	return {"module_type"};
      }
    };

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Table<PhysConfig> psphys{Name("physics")}; //, Comment("Physics parameter set") 
      fhicl::Atom<int> verbosityLevel{Name("verbosityLevel"), Comment("verbosity level (>=0)"), 0};
      fhicl::Atom<double> czmin{Name("czmin"), Comment("Minimum cos(theta) value"), -1.};
      fhicl::Atom<double> czmax{Name("czmax"), Comment("Maximum cos(theta) value"),  1.};
      fhicl::Atom<double> phimin{Name("phimin"), Comment("Minimum phi value (0 - 2 pi))"), 0.};
      fhicl::Atom<double> phimax{Name("phimax"), Comment("Maximum phi value (0 - 2 pi)"),  CLHEP::twopi};
      fhicl::Table<RTS::Config> stops{Name("muonStops"), Comment("Muon stops parameter set")};
      fhicl::Atom<bool> doHistograms{Name("doHistograms"), Comment("Whether or not to make generation histograms (true/false)"), false};
      fhicl::Atom<bool> doCosWeights{Name("doCosWeights"), Comment("Whether or not to use a user cos(theta) generation function (true/false)"), false};
      fhicl::Atom<bool> doEnergyWeights{Name("doEnergyWeights"), Comment("Whether or not to use a user energy generation function (true/false)"), false};
      fhicl::Atom<std::string> energyFuncString{Name("energyFuncString"), Comment("User energy generation function (string)"), ""};
      fhicl::Atom<std::string> czFuncString{Name("czFuncString"), Comment("User cos(theta) generation function (string)"), ""};
    };
    typedef art::EDProducer::Table<Config> Parameters;

    fhicl::ParameterSet psphys_;

    static BinnedSpectrum parseSpectrumShape(const fhicl::ParameterSet& psphys,
                                             double *elow,
                                             double *ehi);

    int                 verbosityLevel_;
    GenId               genId_         ;    // GenId::InternalRMC or GenId::ExternalRMC

    //    int                 generateInternalConversion_;

    double              czmin_;
    double              czmax_;
    double              phimin_;
    double              phimax_;

    art::RandomNumberGenerator::base_engine_t& eng_;

    CLHEP::RandFlat     randomFlat_;
    RandomUnitSphere    randomUnitSphere_;
    MuonCaptureSpectrum muonCaptureSpectrum_;

    RootTreeSampler<IO::StoppedParticleF> stops_;

    bool doHistograms_;

    // generating weighted cos(theta) and energy spectra
    bool doCosWeights_;
    bool doEnergyWeights_;
    std::string energyFuncString_;
    std::string czFuncString_;

    BinnedSpectrum spectrum_;
    CLHEP::RandGeneral*  randSpectrum_;

    TF1* cosThetaFunc_;
    TF1* energyFunc_;

    double ehi_;
    double elow_;
    double enweightnorm_; //for weights
    double czweightnorm_;
    
    double generateEnergy();

    //only used for cos generation weights
    double generateTheta();

    // for weights to produce flat spectra from generation spectra
    double getEnergyWeight(double energy);
    double getCosWeight(double energy);


    TH1F* _hmomentum;
    TH1F* _hcos;
    TH1F* _hmomentumwt;
    TH1F* _hcoswt;
    TH1F* _hElecMom {nullptr};
    TH1F* _hPosiMom {nullptr};
    TH1F* _hTotMom {nullptr};
    TH1F* _hMee;
    TH2F* _hMeeVsE;
    TH1F* _hMeeOverE;                   // M(ee)/E(gamma)
    TH1F* _hy;				// splitting function

  public:
    explicit RMCGun(const Parameters& pset);
    virtual void produce(art::Event& event);
  };

  //================================================================
  RMCGun::RMCGun(const Parameters& pset)
    : EDProducer{pset}
    , psphys_()
    , verbosityLevel_            (pset().verbosityLevel())
    //    , generateInternalConversion_{pset().psphys().generateIntConversion()}
    , czmin_                     (pset().czmin())
    , czmax_                     (pset().czmax())
    , phimin_                    (pset().phimin())
    , phimax_                    (pset().phimax())
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randomFlat_         (eng_)
    , randomUnitSphere_   (eng_, czmin_,czmax_,phimin_,phimax_)
    , muonCaptureSpectrum_(&randomFlat_,&randomUnitSphere_)
    , stops_(eng_, pset().stops())
    , doHistograms_( pset().doHistograms() )
    , doCosWeights_( pset().doCosWeights() )
    , doEnergyWeights_( pset().doEnergyWeights() )
    , energyFuncString_( pset().energyFuncString())
    , czFuncString_    ( pset().czFuncString())
    , spectrum_ ()
    , randSpectrum_(0)
  {
    produces<mu2e::GenParticleCollection>();
    produces<mu2e::EventWeight>();

    //Binned spectrum still takes a fhicl parameter set
    const fhicl::ParameterSet physps = pset().psphys.get_PSet();
    psphys_ = fhicl::ParameterSet(physps);

    genId_  =  GenId::findByName(psphys_.get<std::string>("genId"));
//-----------------------------------------------------------------------------
// make sure the requested genID makes sense
//-----------------------------------------------------------------------------
    if ((genId_ != GenId::InternalRMC) && (genId_ != GenId::ExternalRMC)) {
      throw cet::exception("BADCONFIG") << "RMCGun: wrong process ID: " 
					<< genId_.name() << "BAIL OUT!\n";
    }

    if(verbosityLevel_ > 1) {
      std::cout << "RMCGun: Physics parameter set: \n" << physps.to_indented_string().c_str() << std::endl;
    }
    spectrum_ = BinnedSpectrum(psphys_);
    randSpectrum_ = new CLHEP::RandGeneral(eng_, spectrum_.getPDF(), spectrum_.getNbins());

    if(verbosityLevel_ > 0) {
      std::cout<<"RMCGun: using = "
               <<stops_.numRecords()
               <<" stopped particles"
               <<std::endl;

      std::cout << "RMCGun: producing " 
		<< ((genId_ == GenId::InternalRMC) ? "internal conversion" : "photon")
		<< std::endl;
    }

    if ( doHistograms_ ) {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "RMCGun" );

      _hmomentum     = tfdir.make<TH1F>( "hmomentum", "Produced photon momentum", 1000,  0.,  200.  );
      _hcos          = tfdir.make<TH1F>( "hcos"     , "Produced photon cos(#theta)", 1000,  -1.,  1.  );
      _hmomentumwt   = tfdir.make<TH1F>( "hmomentumwt", "Weighted Produced photon momentum", 1000,  0.,  200.  );
      _hcoswt        = tfdir.make<TH1F>( "hcoswt"     , "Weighted Produced photon cos(#theta)", 1000,  -1.,  1.  );

      if (genId_ == GenId::InternalRMC) {
        _hElecMom  = tfdir.make<TH1F>("hElecMom" , "Produced electron momentum", 1000,  0. , 200.);
        _hPosiMom  = tfdir.make<TH1F>("hPosiMom" , "Produced positron momentum", 1000,  0. , 200.);
        _hTotMom   = tfdir.make<TH1F>("hTotMom" , "Produced total momentum", 200,  0. , 200.);
        _hMee      = tfdir.make<TH1F>("hMee"     , "M(e+e-) "           , 200,0.,200.);
        _hMeeVsE   = tfdir.make<TH2F>("hMeeVsE"  , "M(e+e-) vs E"       , 200,0.,200.,200,0,200);
        _hMeeOverE = tfdir.make<TH1F>("hMeeOverE", "M(e+e-)/E "         , 200, 0.,1);
        _hy        = tfdir.make<TH1F>("hy"       , "y = (ee-ep)/|pe+pp|", 200,-1.,1.);
      }
    }

    //initialize TF1 for cos/energy generation if needed
    if(doCosWeights_) {
      if(genId_ == GenId::InternalRMC)
	throw cet::exception("BADCONFIG")
	  << "RMCGun: Cos generation weights not defined for internal conversions!\n";
      else {
	cosThetaFunc_ = new TF1("cosThetaFunc", czFuncString_.c_str());
	if(cosThetaFunc_->GetMinimum(czmin_,czmax_) < 0.)
	  throw cet::exception("BADCONFIG")
	    << "RMCGun: Cos generation weight function not positive for the entire spectrum!\n";
	  
	czweightnorm_ = cosThetaFunc_->Integral(czmin_,czmax_)/(czmax_-czmin_);
	if(verbosityLevel_ > 0) 
	  printf("RMCGun: Cos generation using function string %s with %.2f < cz < %.2f\n",
		 czFuncString_.c_str(), czmin_, czmax_);
      }
    }
    if(doEnergyWeights_) {
      ehi_ = psphys_.get<double>("ehi");
      elow_ = psphys_.get<double>("elow");
      if(psphys_.get<std::string>("spectrumShape") != "flat")
	throw cet::exception("BADCONFIG")
	  << "RMCGun: Energy generation weights and physics weights both defined\n";
      else {
	energyFunc_ = new TF1("energyFunc", energyFuncString_.c_str());
	if(energyFunc_->GetMinimum(elow_,ehi_) < 0.)
	  throw cet::exception("BADCONFIG")
	    << "RMCGun: Energy generation weight function not positive for the entire spectrum!\n";
	enweightnorm_ = energyFunc_->Integral(elow_,ehi_)/(ehi_-elow_);
	if(verbosityLevel_ > 0) 
	  printf("RMCGun: Energy generation using function string %s with %.2f MeV < E < %.2f MeV\n",
		 energyFuncString_.c_str(), elow_, ehi_);

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
        << "RMCGun: unknown spectrum shape "<<spectrumShape<<"\n";
    }

    return res;
  }

  //================================================================
  void RMCGun::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    const auto& stop = stops_.fire();

    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

    const double energy = generateEnergy();
    double weight = getEnergyWeight(energy);

    if (genId_ != GenId::InternalRMC) {
      if(!doCosWeights_) {
	CLHEP::HepLorentzVector mom( randomUnitSphere_.fire(energy), energy);
	double cosTheta = mom.cosTheta();
	if(doHistograms_) _hcos->Fill(cosTheta);
	if(doHistograms_) _hcoswt->Fill(cosTheta); //no cos(theta) weights

	output->emplace_back( PDGCode::gamma,
			      GenId::ExternalRMC,
			      pos,
			      mom,
			      stop.t );

	event.put(std::move(output));
	std::unique_ptr<mu2e::EventWeight> evtwt ( new EventWeight(weight) );
	event.put(std::move(evtwt));
      } else {

	CLHEP::Hep3Vector mom(0.,0.,0.);
	double theta = generateTheta();
	mom.setRThetaPhi(energy, theta, (phimax_-phimin_)*randomFlat_.fire());
	double coswt = getCosWeight(cos(theta));
	if(doHistograms_) _hcos->Fill(cos(theta));
	if(doHistograms_) _hcoswt->Fill(cos(theta), coswt);

	output->emplace_back( PDGCode::gamma,
			      GenId::ExternalRMC,
			      pos,
			      CLHEP::HepLorentzVector( mom, energy),
			      stop.t );

	event.put(std::move(output));
	weight *= coswt;
	std::unique_ptr<mu2e::EventWeight> evtwt ( new EventWeight(weight) );
	event.put(std::move(evtwt));

      }
    } else {
//-----------------------------------------------------------------------------
// generate internal RMC (off-shell photons)
//-----------------------------------------------------------------------------
      CLHEP::HepLorentzVector mome, momp;
      muonCaptureSpectrum_.getElecPosiVectors(energy,mome,momp);
      // Add particles to list
      auto output = std::make_unique<GenParticleCollection>(); //GenID = 42
      output->emplace_back(PDGCode::e_minus, GenId::InternalRMC,pos,mome,stop.t);
      output->emplace_back(PDGCode::e_plus , GenId::InternalRMC,pos,momp,stop.t);
      event.put(move(output));
      std::unique_ptr<mu2e::EventWeight> evtwt ( new EventWeight(weight) );
      event.put(std::move(evtwt));

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
	double cz = p.cosTheta();
	_hcos->Fill(cz);
	_hcoswt->Fill(cz); //no cos(theta) weights
	
        _hy->Fill(y);
      }
    }

    if ( !doHistograms_ ) return;

    _hmomentum->Fill(energy);
    _hmomentumwt->Fill(energy, getEnergyWeight(energy));

  }

  //================================================================
  double RMCGun::generateEnergy() {
    if(!doEnergyWeights_)
      return spectrum_.sample(randSpectrum_->fire());
    double maxE = energyFunc_->GetMaximum(elow_, ehi_);
    double energy = elow_ + (ehi_-elow_)*randomFlat_.fire();
    double r = randomFlat_.fire();
    while(energyFunc_->Eval(energy) < maxE*r) {
      energy = elow_ + (ehi_-elow_)*randomFlat_.fire();
      r = randomFlat_.fire();
    }
    return energy;
  }

  //================================================================  
  double RMCGun::generateTheta() {
    
    double maxC = cosThetaFunc_->GetMaximum(czmin_, czmax_);
    double r = randomFlat_.fire();
    double cz = czmin_ + (czmax_-czmin_)*randomFlat_.fire();
    while(cosThetaFunc_->Eval(cz) < maxC*r) {
      r = randomFlat_.fire();
      cz = czmin_ + (czmax_-czmin_)*randomFlat_.fire();      
    }
    double theta = acos(cz);
    return theta;
    
  }

  //================================================================
  double RMCGun::getEnergyWeight(double energy) {
    if(!doEnergyWeights_)
      return 1.;
    return enweightnorm_/energyFunc_->Eval(energy); //weight = pdf(flat,E)/pdf(gen,E)
  }

  //================================================================
  double RMCGun::getCosWeight(double cz) {
    if(!doCosWeights_)
      return 1.;
    return czweightnorm_/cosThetaFunc_->Eval(cz); //weight = pdf(flat,cz)/pdf(gen,cz)
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::RMCGun);

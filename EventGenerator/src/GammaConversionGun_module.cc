// Michael MacKenzie, 2019
// Samples photon conversions stored and outputs e+e- pair conversions
// Assumes G4 is part of the trigger path, will cause an exception otherwise

// C++ includes
#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>
#include <map>

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
#include "art/Framework/Principal/SubRun.h"
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
#include "MCDataProducts/inc/GenEventCount.hh"
#include "Mu2eUtilities/inc/Table.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "Mu2eUtilities/inc/GammaPairConversionSpectrum.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"
#include "TF1.h"

namespace mu2e {

  //================================================================
  class GammaConversionGun : public art::EDProducer {
  public:

    typedef RootTreeSampler<IO::ConversionPointF> RTS;
    
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<int> verbosityLevel{Name("verbosityLevel"), Comment("verbosity level (>=0)"), 0};
      fhicl::Table<RTS::Config> stops{Name("gammaStops"), Comment("Gamma stops parameter set")};
      fhicl::Atom<bool> doHistograms{Name("doHistograms"), Comment("Whether or not to make generation histograms (true/false)"), false};
      fhicl::Atom<double> xMin{Name("xMin"), Comment("Stop minimum x value (mm)"), -1.e9};
      fhicl::Atom<double> xMax{Name("xMax"), Comment("Stop maximum x value (mm)"),  1.e9};
      fhicl::Atom<double> yMin{Name("yMin"), Comment("Stop minimum y value (mm)"), -1.e9};
      fhicl::Atom<double> yMax{Name("yMax"), Comment("Stop maximum y value (mm)"),  1.e9};
      fhicl::Atom<double> zMin{Name("zMin"), Comment("Stop minimum z value (mm)"), -1.e9};
      fhicl::Atom<double> zMax{Name("zMax"), Comment("Stop maximum z value (mm)"),  1.e9};
      fhicl::Atom<double> rMin{Name("rMin"), Comment("Stop minimum radius in the DS (mm)"), -1.};
      fhicl::Atom<double> rMax{Name("rMax"), Comment("Stop maximum radius in the DS (mm)"),  1.e9};
      fhicl::Atom<double> czMin{Name("czMin"), Comment("Stop minimum cos(theta) value"), -1.};
      fhicl::Atom<double> czMax{Name("czMax"), Comment("Stop maximum cos(theta) value"),  1.};
      fhicl::Atom<double> pMin{Name("pMin"), Comment("Minimum photon daughter momentum (MeV/c)"), -1.};
      fhicl::Atom<double> pMax{Name("pMax"), Comment("Maximum photon generated energy (MeV/c) ( < 0 to ignore)"),  -1.};
      fhicl::Atom<int>    defaultZ{Name("defaultMaterialZ"), Comment("Override ntuple material with a given material Z"),  -1};
      fhicl::Atom<double> testE{Name("testE"), Comment("Test photon energy to override ntuple energy with (MeV/c) ( < 0 to ignore)"), -1.};
      fhicl::Atom<int>    requireCharge{Name("requireCharge"), Comment("Require a specific photon daughter to pass the cuts"), 0};
      fhicl::Atom<bool>   useCorrelatedAngleOverKE{Name("useCorrelatedAngleOverKE"), Comment("Flag to use correlated e+e- cos/KE"), false};
      fhicl::Atom<double> xOffset{Name("solenoidXOffset"), Comment("X coordinate offset for radius calculations (mm)"), -3904.};
    };
    typedef art::EDProducer::Table<Config> Parameters;


    explicit GammaConversionGun(const Parameters& pset);
    virtual void produce(art::Event& event) override;
    virtual void endSubRun(art::SubRun& sr) override;
    int                 verbosityLevel_;

    art::RandomNumberGenerator::base_engine_t& eng_;

    CLHEP::RandFlat     randomFlat_;
    RTS stops_;
    // PairProduction pairProd_;

    bool doHistograms_;
    //for restricted space generation
    double xMin_;
    double xMax_;
    double yMin_;
    double yMax_;
    double zMin_;
    double zMax_;
    double rMin_; //defined from the ds axis
    double rMax_; 
    //minimum momentum of a daughter
    double pMin_;
    //maximum momentum of a photon
    double pMax_;
    //photon cos(theta) restriction
    double czMin_;
    double czMax_;

    int defaultZ_; //for testing the pair production spectrum for a given material
    double testE_; //for testing the pair production spectrum for a given Energy
    int requireCharge_; //require specific charge daughter to pass the cuts

    bool   useCorrelatedAngleOverKE_; //correlate the cos/ke for the e+ and e-
    GammaPairConversionSpectrum* spectrum_; //pair production spectrum
    double xOffset_; //ds axis x offset

    GenEventCount::count_t genEvents_; //for normalization
    GenEventCount::count_t passedEvents_; //produced statistics

    TH1D* _hgencuts; //records number of events attempted
    TH1F* _hmomentum;
    TH1F* _hcos;
    TH1F* _hcosWt;
    TH1F* _hr; //x-y r from (x,y) = (-3904, 0)
    TH1F* _hz; 
    TH1F* _hElecMom  {nullptr};
    TH1F* _hPosiMom  {nullptr};
    TH1F* _hTotMom   {nullptr};
    TH1F* _hTotMomWt {nullptr};
    TH1F* _hTotEnergy   {nullptr};
    TH1F* _hTotEnergyWt {nullptr};
    TH1F* _hMee;
    TH2F* _hMeeVsE;
    TH1F* _hMeeOverE;                   // M(ee)/E(gamma)
    TH1F* _hy;				// splitting function
    TH2F* _hcosEvsP;                    //E(e-)/M(e-)*theta(e-) vs E(e+)/M(e+)*theta(e+) wrt photon direction
    TH1F* _henergyChange;
    TH1F* _hcosChange;
    TH1F* _hRecoilMsq; //Mass^2 of recoiling particle 
  };

  //================================================================
  GammaConversionGun::GammaConversionGun(const Parameters& pset)
    : EDProducer(pset)
    , verbosityLevel_          (pset().verbosityLevel())
    , eng_(createEngine        (art::ServiceHandle<SeedService>()->getSeed()))
    , randomFlat_              (eng_)
    , stops_                   (eng_, pset().stops())
    , doHistograms_            (pset().doHistograms())
    , xMin_                    (pset().xMin())
    , xMax_                    (pset().xMax())
    , yMin_                    (pset().yMin())
    , yMax_                    (pset().yMax())
    , zMin_                    (pset().zMin())
    , zMax_                    (pset().zMax())
    , rMin_                    (pset().rMin())
    , rMax_                    (pset().rMax())
    , pMin_                    (pset().pMin())
    , pMax_                    (pset().pMax())
    , czMin_                   (pset().czMin())
    , czMax_                   (pset().czMax())
    , defaultZ_                (pset().defaultZ())
    , testE_                   (pset().testE())
    , requireCharge_           (pset().requireCharge())
    , useCorrelatedAngleOverKE_(pset().useCorrelatedAngleOverKE())
    , spectrum_                (new GammaPairConversionSpectrum(&randomFlat_, useCorrelatedAngleOverKE_))
    , xOffset_                 (pset().xOffset())
    , genEvents_       (0)
    , passedEvents_    (0)
  {
    produces<mu2e::GenParticleCollection>();
    produces<mu2e::EventWeight>();
    produces<mu2e::GenParticleCollection>("photon"); //store photon generation energy for RMC weights
    produces<mu2e::GenEventCount, art::InSubRun>("genEvents"); //for normalization


    if(verbosityLevel_ > 0) {
      std::cout<<"GammaConversionGun: using = "
               <<stops_.numRecords()
               <<" stopped gammas"
               <<std::endl;

      std::cout<<"GammaConversionGun: producing photon " << std::endl;
      std::cout << "GammaConversionGun: useCorrelatedAngleOverKE = " << useCorrelatedAngleOverKE_ << std::endl;
    }
    if(xMin_ >= xMax_ || yMin_ >= yMax_ || zMin_ >= zMax_ || rMin_ >= rMax_)
      throw cet::exception("BADCONFIG")
	<< "GammaConversionGun: Stop (x,y,z,r) restriction error! "
	<< xMin_ << " < x < " << xMax_ << ", "
	<< yMin_ << " < y < " << yMax_ << ", "
	<< zMin_ << " < z < " << zMax_ << ", "
	<< rMin_ << " < r < " << rMax_ << "\n";

    if(defaultZ_ > 0)
      std::cout << "GammaConversionGun: Overriding Ntuple defined material and instead using Z = "
		<< defaultZ_ << std::endl;
    if(testE_ > 0) {
      std::cout << "GammaConversionGun: Overriding Ntuple defined energy and instead using "
		<< testE_ << std::endl;
      if(pMin_ > testE_)
	throw cet::exception("BADCONFIG")
	  << "GammaConversionGun: Test energy set to value below minimum electron/positron momentum!\n";
      if(pMax_ > 0. && pMax_ < testE_)
	throw cet::exception("BADCONFIG")
	  << "GammaConversionGun: Test energy set to value above maximum photon momentum!\n";

    }

    //essentially a weight histogram, so always store
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir( "GammaConversionGun" );
    _hgencuts  = tfdir.make<TH1D>( "hgencuts", "Attempts to generate vs cut number", 10,  0.5,  10.5  );
    if ( doHistograms_ ) {

      _hmomentum     = tfdir.make<TH1F>("hmomentum", "Given photon momentum", 1000,  0.,  200.  );
      _hcos          = tfdir.make<TH1F>("hcos"     , "Given photon cos(#theta)", 1000,  -1.,  1.  );
      _hcosWt        = tfdir.make<TH1F>("hcosWt"   , "Given photon cos(#theta) weighted", 1000,  -1.,  1.  );
      _hr            = tfdir.make<TH1F>("hr"       , "Given photon radius from the DS axis", 1000,  0.,  2000.  );
      _hz            = tfdir.make<TH1F>("hz"       , "Given photon Z", 1000,  0., 15000.  );
      _hElecMom      = tfdir.make<TH1F>("hElecMom" , "Produced electron momentum", 2000,  0. , 200.);
      _hPosiMom      = tfdir.make<TH1F>("hPosiMom" , "Produced positron momentum", 2000,  0. , 200.);
      _hTotMom       = tfdir.make<TH1F>("hTotMom"   , "Produced total momentum", 2000.,  0. , 200.);
      _hTotMomWt     = tfdir.make<TH1F>("hTotMomWt" , "Produced total momentum weighted", 2000.,  0. , 200.);
      _hTotEnergy    = tfdir.make<TH1F>("hEnergyMom"   , "Produced total energy", 1000.,  0. , 200.);
      _hTotEnergyWt  = tfdir.make<TH1F>("hEnergyMomWt" , "Produced total energy weighted", 1000.,  0. , 200.);
      _hMee          = tfdir.make<TH1F>("hMee"     , "M(e+e-) "           , 200,0.,200.);
      _hMeeVsE       = tfdir.make<TH2F>("hMeeVsE"  , "M(e+e-) vs E"       , 200,0.,200.,200,0,200);
      _hMeeOverE     = tfdir.make<TH1F>("hMeeOverE", "M(e+e-)/E "         , 200, 0.,1);
      _hy            = tfdir.make<TH1F>("hy"       , "y = (ee-ep)/|pe+pp|", 200,-1.,1.);
      _hcosEvsP      = tfdir.make<TH2F>("hcosEvsP" , "p*theta(e-) vs p*theta(e+) wrt photon", 150, 0., 15., 150, 0., 15.);
      _hcosChange    = tfdir.make<TH1F>("hcosChange", "cos(#theta) Change", 200,  -1.,  1.  );
      _henergyChange = tfdir.make<TH1F>("henergyChange", "E(photon) - E(e+e-)", 200,-1.,1.);
      _hRecoilMsq    = tfdir.make<TH1F>("hRecoilMsq", "Mass^2(p(photon) - p(e+e-))", 240,-10.,50.);

    }

  }

  //================================================================
  void GammaConversionGun::produce(art::Event& event) {

    
    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
    
    //fields in the stop ntuple 
    double x, y, z, t, px, py, pz;
    double weight;
    double gen_energy;
    int matN;
    int maxElem = IO::kMaxConversionMaterialElements;
    int matZ[maxElem];
    double matZeff[maxElem], matFrac[maxElem];
    

    bool passed = false;
    CLHEP::Hep3Vector pos(0.,0.,0.);
    CLHEP::HepLorentzVector mome, momp, momg;

    //loop through generations until an event passes the generation cuts
    do {
      const auto& stop = stops_.fire();
      x  = stop.x;  y  = stop.y;  z  = stop.z;  t = stop.time;
      px = stop.px; py = stop.py; pz = stop.pz; 
      weight = stop.weight; gen_energy = stop.genEnergy;
      if(defaultZ_ > 0) {matN = 1; matZ[0] = defaultZ_; matZeff[0] = defaultZ_; matFrac[0] = 1.;}
      else {
	matN = stop.matN;
	for(int index = 0; index < matN; ++index) {
	  matZ   [index] = stop.matZ   [index];
	  matZeff[index] = stop.matZeff[index];
	  matFrac[index] = stop.matFrac[index];
	}
      }
      
      if(verbosityLevel_ > 2) {
	std::cout << "Next Stop Attempt: (x,y,z,t) = (" << x << "," << y
		  << "," << z << "," << t << "), "
		  << "(px,py,pz,gen_energy) = (" << px << "," << py << ","
		  << pz << "," << gen_energy << ") " << std::endl
		  << "material: N(Elements) = " << matN << " (Z, Zeff, Fraction) -->";
	for(int index = 0; index < matN; ++index) {
	  std::cout << " (" << matZ[index] << ", " << matZeff[index] 
		    << ", " << matFrac[index] << ")";
	}
	std::cout << std::endl;
      }
      _hgencuts->Fill(1); //all generations
      ++genEvents_;
      passed = !(x < xMin_ || x > xMax_
		 || y < yMin_ || y > yMax_
		 || z < zMin_ || z > zMax_
		 || sqrt((x-xOffset_)*(x-xOffset_) + y*y) < rMin_ 
		 || sqrt((x-xOffset_)*(x-xOffset_) + y*y) > rMax_ );

      if(!passed) continue;
      _hgencuts->Fill(2); //passed spacial cut
      if(verbosityLevel_ > 2) 
	std::cout << "Passed spacial cut\n";

      CLHEP::Hep3Vector mom(px,py,pz);
      if(testE_ > 0.) mom.setMag(testE_);
      momg.setVect(mom);
      momg.setE(mom.mag());

      passed = passed && (pMax_ < 0. || pMax_ > gen_energy);
      if(!passed) continue;
      _hgencuts->Fill(3); //passed maximum gen photon momentum cut
      if(verbosityLevel_ > 2) 
	std::cout << "Passed gen energy cut\n";

      double cz = mom.cosTheta();
      passed = passed && cz >= czMin_ && cz <= czMax_;
      if(!passed) continue;
      _hgencuts->Fill(4); //Cos(theta) cut

      //can't make a daughter of pMin if energy below pmin + electron mass already
      //use slightly less than electron mass here to be safe
      double photonE = momg.e();
      if(verbosityLevel_ > 2) 
	std::cout << "Passed cos theta cut\nPhoton energy = " << photonE << std::endl;
      passed = passed && (photonE - 0.500) > pMin_;
      if(!passed) continue;
      _hgencuts->Fill(5); //passed initial minimum momentum cut
      if(verbosityLevel_ > 2) 
	std::cout << "Passed min photon energy cut\n";

      //create a corresponding material
      GammaPairConversionSpectrum::materialData material;
      for(int index = 0; index < matN; ++index) {
	material.elements.push_back(spectrum_->_elementMap[matZ[index]]);
	material.elementFractions.push_back(matFrac[index]);
      }
      //sample the spectrum
      spectrum_->fire(momg, material, mome, momp);

      if(requireCharge_ == 0) 
	passed = passed && (mome.vect().mag() > pMin_ || momp.vect().mag() > pMin_);
      else
	//charge specific cut
	passed = passed && ((requireCharge_ > 0 && momp.vect().mag() > pMin_) || (requireCharge_ < 0 && mome.vect().mag() > pMin_));
      if(!passed) continue;
      _hgencuts->Fill(6); //Final momentum cut
      if(verbosityLevel_ > 2) 
	std::cout << "Passed min daughter energy cut\n";

      pos.setX(x); pos.setY(y); pos.setZ(z);

    } while(!passed);
    _hgencuts->Fill(10); //passing all cuts
    ++passedEvents_;
    output->emplace_back(PDGCode::e_minus, GenId::gammaPairProduction, pos, mome, t);
    output->emplace_back(PDGCode::e_plus , GenId::gammaPairProduction, pos, momp, t);

    event.put(move(output));
    //add event weight and gen photon energy to the output
    std::unique_ptr<EventWeight> evtwt ( new EventWeight(weight) );
    event.put(move(evtwt));
    std::unique_ptr<GenParticleCollection> output_photon(new GenParticleCollection);
    output_photon->emplace_back(PDGCode::gamma, GenId::ExternalRMC,CLHEP::Hep3Vector(0.,0.,0.),CLHEP::HepLorentzVector(0.,0.,gen_energy,gen_energy),0.);
    event.put(move(output_photon),"photon"); 

    if ( !doHistograms_ ) return;

    _hcos->Fill((mome+momp).vect().cosTheta());
    _hcosWt->Fill((mome+momp).vect().cosTheta(),weight);
    _hr->Fill(sqrt((pos.x()+3904.)*(pos.x()+3904.)+pos.y()*pos.y()));
    _hz->Fill(pos.z());
    _hElecMom ->Fill(mome.vect().mag());
    _hPosiMom ->Fill(momp.vect().mag());
    CLHEP::Hep3Vector p = mome.vect()+momp.vect();
    double momentum = p.mag();
    _hTotMom ->Fill(momentum);
    _hTotMomWt ->Fill(momentum, weight);

    double mee = (mome+momp).m();
    _hMee->Fill(mee);
    double energy = (mome+momp).e();
    _hTotEnergy  ->Fill(energy);
    _hTotEnergyWt->Fill(energy, weight);
    _hMeeVsE->Fill(energy,mee);
    _hMeeOverE->Fill(mee/energy);

    double lepy = (mome.e()-momp.e())/energy;

    _hy->Fill(lepy);

    _hmomentum->Fill(energy);
    _henergyChange->Fill(momg.e()-energy);
    _hcosChange->Fill(momg.vect().cosTheta()-(mome+momp).vect().cosTheta());

    CLHEP::Hep3Vector p_e = mome.vect(), p_p = momp.vect(), p_g = momg.vect();
    _hcosEvsP->Fill((p_p.angle(p_g))*momp.e()/0.511, (p_e.angle(p_g))*mome.e()/0.511);

    CLHEP::HepLorentzVector recoil = momg - (mome+momp);
    double recoil_msq = recoil.m2();
    _hRecoilMsq->Fill(recoil_msq);
    
  }

  //================================================================
  void GammaConversionGun::endSubRun(art::SubRun& sr) {

    mf::LogInfo("Summary")<<"Creating GenEventCount record in the GammaConversionGun: "<< genEvents_
                          <<" generated events, where"
			  << passedEvents_ << " events were produced for "<<sr.id()<<"\n";

    sr.put(std::unique_ptr<GenEventCount>(new GenEventCount(genEvents_)), "genEvents");
  }


  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::GammaConversionGun);

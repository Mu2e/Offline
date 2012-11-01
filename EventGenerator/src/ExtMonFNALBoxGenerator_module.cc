// $Id: ExtMonFNALBoxGenerator_module.cc,v 1.3 2012/11/01 23:41:46 gandr Exp $
// $Author: gandr $
// $Date: 2012/11/01 23:41:46 $
//
// Create particle flux in the ExtMonFNAL box by randomizing
// kinematic of input particles read from a file.
// The output is in the Mu2e coordinate system.
//
//
// Original author Andrei Gaponenko, 2012

#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>
#include <iterator>
#include <map>

#include "cetlib/exception.h"

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/MassCache.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/MARSInfo.hh"
#include "MCDataProducts/inc/MARSInfoCollection.hh"
#include "MCDataProducts/inc/GenParticleMARSAssns.hh"

#include "ExtinctionMonitorFNAL/Utilities/inc/EMFBoxIO.hh"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)


namespace mu2e {
  namespace ExtMonFNAL {

    namespace {
      struct InputHit {
        IO::EMFBoxHit particle;
        double energy;
        MARSInfo minfo;
        IO::G4JobInfo g4s1info;
        IO::ParticleRandomization pr;

        InputHit(const IO::EMFBoxHit& p,
                 double e,
                 const MARSInfo& mi,
                 const IO::G4JobInfo& gi,
                 const IO::ParticleRandomization& rr)
          : particle(p), energy(e), minfo(mi), g4s1info(gi), pr(rr)
        {}
      };

      typedef std::vector<InputHit> InputHits;

      //================================================================
      struct InputStop {
        IO::StoppedMuon muon;
        MARSInfo minfo;
        IO::G4JobInfo g4s1info;
        InputStop(const IO::StoppedMuon& m, const MARSInfo& i, const IO::G4JobInfo& g4s1)
          : muon(m), minfo(i), g4s1info(g4s1)
        {}
      };

      typedef std::vector<InputStop> InputStops;

      //================================================================
    } // namespace {}

    class ExtMonFNALBoxGenerator : public art::EDProducer {
      int verbosityLevel_;
      std::string geomModuleLabel_;
      std::string geomInstanceName_;

      std::vector<std::string> inputFiles_;

      double microbunchScaling_; // the fraction of all inputs to be used to generate a microbunch

      // w==1 particles are accepted wiht 1/weightMax_ probability
      // w>max weights are truncated
      double weightMax_;

      unsigned microbunchVDHitsChunkSize_; // computed from inputs size and the two numbers above
      unsigned nextVDHit_; // index

      unsigned microbunchStoppedMuonsChunkSize_;
      unsigned nextStoppedMuon_;

      double bunchTimeHalfWidth_;
      double cutTimeMin_;
      double cutTimeMax_;
      std::vector<double> keepInBox_;

      art::RandomNumberGenerator::base_engine_t& eng_;
      CLHEP::RandFlat randFlat_;
      CLHEP::RandGaussQ randGauss_;

      double deBuncherPeriod_;
      const ExtMon *extmon_;

      MassCache mc_;

      InputHits vdhits_;
      InputStops mustops_;

      void loadInputFiles();

      void addVDHits(TFile* infile);
      void addStoppedMuons(TFile* infile);

      void generateFromVDHits(const art::Event& event,
                              GenParticleCollection *output,
                              MARSInfoCollection *info,
                              GenParticleMARSAssns *assns);
      GenParticle createOutputParticle(double randomizedTime, const InputHit& hit);
      bool inRange(int vdId, const CLHEP::Hep3Vector& posExtMon);

      void generateFromStoppedMuons(const art::Event& event,
                                    GenParticleCollection *output,
                                    MARSInfoCollection *info,
                                    GenParticleMARSAssns *assns);

      GenParticle createOutputMuon(double muonTime, const InputStop& ms);

    public:
      explicit ExtMonFNALBoxGenerator(const fhicl::ParameterSet& pset);
      virtual void produce(art::Event& event);
      virtual void beginRun(art::Run& run);
    };

    //================================================================
    ExtMonFNALBoxGenerator::ExtMonFNALBoxGenerator(const fhicl::ParameterSet& pset)
      : verbosityLevel_(pset.get<int>("verbosityLevel"))
      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))

      , inputFiles_(pset.get<std::vector<std::string> >("inputFiles"))
      , microbunchScaling_(pset.get<double>("microbunchScaling"))
      , weightMax_(pset.get<double>("weightMax"))
      , microbunchVDHitsChunkSize_()
      , nextVDHit_()
      , microbunchStoppedMuonsChunkSize_()
      , nextStoppedMuon_()

      , bunchTimeHalfWidth_(pset.get<double>("bunchTimeHalfWidth"))
      , cutTimeMin_(pset.get<double>("cutTimeMin"))
      , cutTimeMax_(pset.get<double>("cutTimeMax"))
      , keepInBox_(pset.get<std::vector<double> >("keepInBox"))

      , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
      , randFlat_(eng_)
      , randGauss_(eng_)

      , deBuncherPeriod_()
      , extmon_()
    {
      produces<mu2e::GenParticleCollection>();
      produces<mu2e::MARSInfoCollection>();
      produces<mu2e::GenParticleMARSAssns>();

      if(inputFiles_.empty()) {
        throw cet::exception("BADCONFIG")<<"Error: no inputFiles";
      }

    }

    //================================================================
    void ExtMonFNALBoxGenerator::beginRun(art::Run& run) {
      ConditionsHandle<AcceleratorParams> ch("ignored");
      deBuncherPeriod_ = ch->deBuncherPeriod;

      if(verbosityLevel_ > 0) {
        std::cout<<"ExtMonFNALBoxGenerator: using deBuncherPeriod = "<<deBuncherPeriod_<<std::endl;
      }

      if(!geomModuleLabel_.empty()) {
        art::Handle<ExtMon> extmon;
        run.getByLabel(geomModuleLabel_, geomInstanceName_, extmon);
        extmon_ = &*extmon;
      }
      else {
        GeomHandle<ExtMon> extmon;
        extmon_ = &*extmon;
      }

      if(verbosityLevel_ > 0) {
        std::cout<<"ExtMonFNALBoxGenerator: detectorCenterInMu2e = "<<extmon_->detectorCenterInMu2e()<<std::endl;
      }

      //----------------------------------------------------------------
      loadInputFiles();

      microbunchVDHitsChunkSize_ = vdhits_.size() * microbunchScaling_ * weightMax_;
      microbunchStoppedMuonsChunkSize_ = mustops_.size() * microbunchScaling_ * weightMax_;

      if(verbosityLevel_ > 0) {
        std::cout<<"ExtMonFNALBoxGenerator inputs: num hits = "<<vdhits_.size()
                 <<", num stopped muons = "<<mustops_.size()
                 <<std::endl;

        std::cout<<"ExtMonFNALBoxGenerator: microbunchVDHitsChunkSize = "<<microbunchVDHitsChunkSize_
                 <<", microbunchStoppedMuonsChunkSize = "<<microbunchStoppedMuonsChunkSize_
                 <<std::endl;
      }

      if(microbunchVDHitsChunkSize_ < 1) {
        throw cet::exception("BADCONFIG")<<"ERROR: Computed microbunchVDHitsChunkSize < 1! Increase weightMax?\n";
      }
      if(microbunchStoppedMuonsChunkSize_ < 1) {
        throw cet::exception("BADCONFIG")<<"ERROR: Computed microbunchStoppedMuonsChunkSize_ < 1! Increase weightMax?\n";
      }

    }

    //================================================================
    void ExtMonFNALBoxGenerator::loadInputFiles() {
      art::ServiceHandle<art::TFileService> tfs;

      for(unsigned i=0; i<inputFiles_.size(); ++i) {

        const std::string resolvedFileName = ConfigFileLookupPolicy()(inputFiles_[i]);
        if(verbosityLevel_ > 0) {
          std::cout<<"ExtMonFNALBoxGenerator: reading input file "<<resolvedFileName<<std::endl;
        }

        TFile *infile = tfs->make<TFile>(resolvedFileName.c_str(), "READ");
        addVDHits(infile);
        addStoppedMuons(infile);

      } // for(files)
    } // loadInputFiles()

    //================================================================
    void ExtMonFNALBoxGenerator::addVDHits(TFile* infile) {
      const std::string hitTreeName_("EMFBoxFluxAnalyzer/vdhits");

      TTree *nt = dynamic_cast<TTree*>(infile->Get(hitTreeName_.c_str()));
      if(!nt) {
        throw cet::exception("BADINPUT")<<"Could not get tree \""<<hitTreeName_
                                        <<"\" from file \""<<infile->GetName()
                                        <<"\"\n";
      }

      const Long64_t nTreeEntries = nt->GetEntries();
      if(verbosityLevel_>0) {
        std::cout<<"ExtMonFNALBoxGenerator: nTreeEntries for VD hits = "<<nTreeEntries<<std::endl;
      }

      IO::EMFBoxHit particle;
      TBranch *bhit = nt->GetBranch("particle");
      bhit->SetAddress(&particle);

      MARSInfo minfo;
      TBranch *bminf = nt->GetBranch("minfo");
      bminf->SetAddress(&minfo);

      IO::G4JobInfo g4s1info;
      TBranch *bginf = nt->GetBranch("g4s1info");
      bginf->SetAddress(&g4s1info);

      IO::ParticleRandomization pr;
      TBranch *brand = nt->GetBranch("randomization");
      brand->SetAddress(&pr);

      for(Long64_t i=0; i<nTreeEntries; ++i) {
        bhit->GetEntry(i);
        bminf->GetEntry(i);
        bginf->GetEntry(i);
        brand->GetEntry(i);

        const double mass = mc_.mass(PDGCode::type(particle.pdgId));
        const double energy = sqrt(std::pow(mass, 2) +
                                   std::pow(particle.mu2epx, 2) +
                                   std::pow(particle.mu2epy, 2) +
                                   std::pow(particle.mu2epz, 2)
                                   );

        vdhits_.push_back(InputHit(particle, energy, minfo, g4s1info, pr));
      }
    }

    //================================================================
    void ExtMonFNALBoxGenerator::addStoppedMuons(TFile* infile) {
      const std::string treeName("StoppedMuons/sm");

      TTree *nt = dynamic_cast<TTree*>(infile->Get(treeName.c_str()));
      if(!nt) {
        throw cet::exception("BADINPUT")<<"Could not get tree \""<<treeName
                                        <<"\" from file \""<<infile->GetName()
                                        <<"\"\n";
      }

      const Long64_t nTreeEntries = nt->GetEntries();
      if(verbosityLevel_>0) {
        std::cout<<"ExtMonFNALBoxGenerator: nTreeEntries for stopped muons = "<<nTreeEntries<<std::endl;
      }

      IO::StoppedMuon mu;
      TBranch *bmu = nt->GetBranch("particle");
      bmu->SetAddress(&mu);

      MARSInfo minfo;
      TBranch *bminf = nt->GetBranch("minfo");
      bminf->SetAddress(&minfo);

      IO::G4JobInfo g4s1info;
      TBranch *bginf = nt->GetBranch("g4s1info");
      bginf->SetAddress(&g4s1info);

      for(Long64_t i=0; i<nTreeEntries; ++i) {
        bmu->GetEntry(i);
        bminf->GetEntry(i);
        bginf->GetEntry(i);
        mustops_.push_back(InputStop(mu, minfo, g4s1info));
      }
    }

    //================================================================
    void ExtMonFNALBoxGenerator::produce(art::Event& event) {

      std::auto_ptr<GenParticleCollection> output(new GenParticleCollection);
      std::auto_ptr<MARSInfoCollection> info(new MARSInfoCollection());
      std::auto_ptr<GenParticleMARSAssns> assns(new GenParticleMARSAssns());

      generateFromVDHits(event, output.get(), info.get(), assns.get());
      generateFromStoppedMuons(event, output.get(), info.get(), assns.get());

      event.put(output);
      event.put(info);
      event.put(assns);
    }

    //================================================================
    void ExtMonFNALBoxGenerator::generateFromStoppedMuons(const art::Event& event,
                                                          GenParticleCollection *output,
                                                          MARSInfoCollection *info,
                                                          GenParticleMARSAssns *assns) {

      const art::ProductID particlesPID = getProductID<GenParticleCollection>(event);
      const art::EDProductGetter *particlesGetter = event.productGetter(particlesPID);

      const art::ProductID marsPID = getProductID<MARSInfoCollection>(event);
      const art::EDProductGetter *marsGetter = event.productGetter(marsPID);

      for(unsigned count=0; count < microbunchStoppedMuonsChunkSize_; ++count) {
        const InputStop& ms = mustops_[nextStoppedMuon_];

        const double randomizedStopTime =
          fmod(ms.muon.time, deBuncherPeriod_)
          + (2*randFlat_.fire() - 1.)*bunchTimeHalfWidth_;

        const double generatedMuonTime = std::max(randomizedStopTime, cutTimeMin_);

        static const double tauMuMinus = 864.; //ns, Al is conservative for Fe
        static const double tauMuPlus = 2197.; //ns, free muon

        const double tau = (ms.muon.pdgId > 0) ? tauMuMinus : tauMuPlus;
        const double weight = exp((generatedMuonTime - cutTimeMin_)/tau) * ms.minfo.weight();

        if(randFlat_.fire() * weightMax_ <= weight) {
          output->push_back(createOutputMuon(generatedMuonTime, ms));
          info->push_back(ms.minfo);
          assns->addSingle(art::Ptr<GenParticle>(particlesPID, output->size()-1, particlesGetter),
                           art::Ptr<MARSInfo>(marsPID, info->size()-1, marsGetter));
        }

        ++nextStoppedMuon_ %= mustops_.size();
      }
    }

    //================================================================
    void ExtMonFNALBoxGenerator::generateFromVDHits(const art::Event& event,
                                                    GenParticleCollection *output,
                                                    MARSInfoCollection *info,
                                                    GenParticleMARSAssns *assns) {

      const art::ProductID particlesPID = getProductID<GenParticleCollection>(event);
      const art::EDProductGetter *particlesGetter = event.productGetter(particlesPID);

      const art::ProductID marsPID = getProductID<MARSInfoCollection>(event);
      const art::EDProductGetter *marsGetter = event.productGetter(marsPID);

      for(unsigned count=0; count < microbunchVDHitsChunkSize_; ++count) {
        const InputHit& hit = vdhits_[nextVDHit_];

        const bool acceptHitWeight(randFlat_.fire() * weightMax_ <= hit.minfo.weight());
        if(acceptHitWeight) {

          const double randomizedTime =
            fmod(hit.particle.time, deBuncherPeriod_)
            + (2*randFlat_.fire() - 1.)*bunchTimeHalfWidth_;

          if((cutTimeMin_ < randomizedTime) && (randomizedTime < cutTimeMax_)) {
            output->push_back(createOutputParticle(randomizedTime, hit));
            info->push_back(hit.minfo);
            assns->addSingle(art::Ptr<GenParticle>(particlesPID, output->size()-1, particlesGetter),
                             art::Ptr<MARSInfo>(marsPID, info->size()-1, marsGetter));
          }

        }
        ++nextVDHit_ %= vdhits_.size();
      }
    }

    //================================================================
    GenParticle ExtMonFNALBoxGenerator::createOutputParticle(double randomizedTime, const InputHit& hit) {
      using CLHEP::Hep3Vector;

      Hep3Vector posExtMon(hit.particle.emx, hit.particle.emy, hit.particle.emz);
      do {
        if((hit.particle.vdId != VirtualDetectorId::EMFBoxSW) &&
           (hit.particle.vdId != VirtualDetectorId::EMFBoxNE)) {
          posExtMon.setX(hit.particle.emx +  hit.pr.sigmax * randGauss_.fire());
        }

        if((hit.particle.vdId != VirtualDetectorId::EMFBoxBottom)&&
           (hit.particle.vdId != VirtualDetectorId::EMFBoxTop)) {
          posExtMon.setY(hit.particle.emy + hit.pr.sigmay * randGauss_.fire());
        }

        if((hit.particle.vdId != VirtualDetectorId::EMFBoxFront)&&
           (hit.particle.vdId != VirtualDetectorId::EMFBoxBack)) {
          posExtMon.setZ(hit.particle.emz + hit.pr.sigmaz * randGauss_.fire());
        }
      } while(!inRange(hit.particle.vdId, posExtMon));

      const Hep3Vector posMu2e(extmon_->extMonToMu2e_position(posExtMon));


      // Draw dtheta from the Rayleigh distribution
      const double dtheta = hit.pr.rSigmaML * sqrt(-2*log(randFlat_.fire()));
      const double dphi = 2*M_PI*randFlat_.fire();

      const Hep3Vector orig(hit.particle.mu2epx, hit.particle.mu2epy, hit.particle.mu2epz);

      // Find a vector that is not collinear with the original particle direction
      const Hep3Vector n1 = (std::abs(orig.x()) < std::abs(orig.y())) ?
        ((std::abs(orig.x()) < std::abs(orig.z())) ? Hep3Vector(1,0,0) : Hep3Vector(0,0,1)) :
        ((std::abs(orig.x()) < std::abs(orig.y())) ? Hep3Vector(1,0,0) : Hep3Vector(0,1,0));

      // Construct a vector perpendicular to the original momentum
      const Hep3Vector perp = orig.cross(n1);

      // Randomize the original direction
      Hep3Vector randomized3mom(orig);
      randomized3mom.rotate(perp, dtheta);
      randomized3mom.rotate(orig, dphi);

      const CLHEP::HepLorentzVector momMu2e(randomized3mom, hit.energy);

      return GenParticle(PDGCode::type(hit.particle.pdgId),
                         GenId::MARS,
                         posMu2e,
                         momMu2e,
                         randomizedTime);
    }

    //================================================================
    bool ExtMonFNALBoxGenerator::inRange(int vdId, const CLHEP::Hep3Vector& posExtMon) {
      switch(vdId) {
      default: break;

      case VirtualDetectorId::EMFBoxFront: case VirtualDetectorId::EMFBoxBack:
        return
          (std::abs(posExtMon.x()) <= keepInBox_[0]) &&
          (std::abs(posExtMon.y()) <= keepInBox_[1]);

      case VirtualDetectorId::EMFBoxSW: case VirtualDetectorId::EMFBoxNE:
        return
          (std::abs(posExtMon.y()) <= keepInBox_[1]) &&
          (std::abs(posExtMon.z()) <= keepInBox_[2]);

      case VirtualDetectorId::EMFBoxBottom: case VirtualDetectorId::EMFBoxTop:
        return
          (std::abs(posExtMon.x()) <= keepInBox_[0]) &&
          (std::abs(posExtMon.z()) <= keepInBox_[2]);
      }

      assert(false);
      return false;
    } // inRange()

    //================================================================
    GenParticle ExtMonFNALBoxGenerator::createOutputMuon(double muonTime, const InputStop& ms) {
      using CLHEP::Hep3Vector;

      Hep3Vector posExtMon(ms.muon.emx, ms.muon.emy, ms.muon.emz);
      const Hep3Vector posMu2e(extmon_->extMonToMu2e_position(posExtMon));

      static const double muonMass = mc_.mass(PDGCode::type(13));
      const CLHEP::HepLorentzVector momMu2e(Hep3Vector(), muonMass);

      return GenParticle(PDGCode::type(ms.muon.pdgId),
                         GenId::MARS,
                         posMu2e,
                         momMu2e,
                         muonTime);
    }

    //================================================================

  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::ExtMonFNALBoxGenerator);

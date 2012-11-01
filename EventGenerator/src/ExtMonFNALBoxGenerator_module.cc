// $Id: ExtMonFNALBoxGenerator_module.cc,v 1.1 2012/11/01 23:35:24 gandr Exp $
// $Author: gandr $
// $Date: 2012/11/01 23:35:24 $
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
        MARSInfo info;
        IO::ParticleRandomization pr;

        InputHit(const IO::EMFBoxHit& p, double e, const MARSInfo& mi, const IO::ParticleRandomization& rr)
          : particle(p), energy(e), info(mi), pr(rr)
        {}
      };

      typedef std::vector<InputHit> InputHits;

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

      void loadInputFiles();
      void addVDHits(TFile* infile);
      void addStoppedMuons(TFile* infile);

      GenParticle createOutputParticle(double randomizedTime, const InputHit& hit);
      bool inRange(int vdId, const CLHEP::Hep3Vector& posExtMon);

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

      if(verbosityLevel_ > 0) {
        std::cout<<"ExtMonFNALBoxGenerator: input hits size = "<<vdhits_.size()<<std::endl;
        std::cout<<"ExtMonFNALBoxGenerator: microbunchVDHitsChunkSize = "<<microbunchVDHitsChunkSize_<<std::endl;
      }

      if(microbunchVDHitsChunkSize_ < 1) {
        throw cet::exception("BADCONFIG")<<"ERROR: Computed microbunchVDHitsChunkSize < 1! Increase weightMax?\n";
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
        //addStoppedMuons(infile);

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
        std::cout<<"ExtMonFNALBoxGenerator: nTreeEntries = "<<nTreeEntries<<std::endl;
      }

      IO::EMFBoxHit particle;
      TBranch *bhit = nt->GetBranch("particle");
      bhit->SetAddress(&particle);

      MARSInfo minfo;
      TBranch *binf = nt->GetBranch("minfo");
      binf->SetAddress(&minfo);

      IO::ParticleRandomization pr;
      TBranch *brand = nt->GetBranch("randomization");
      brand->SetAddress(&pr);

      for(Long64_t i=0; i<nTreeEntries; ++i) {
        bhit->GetEntry(i);
        binf->GetEntry(i);
        brand->GetEntry(i);

        const double mass = mc_.mass(PDGCode::type(particle.pdgId));
        const double energy = sqrt(std::pow(mass, 2) +
                                   std::pow(particle.mu2epx, 2) +
                                   std::pow(particle.mu2epy, 2) +
                                   std::pow(particle.mu2epz, 2)
                                   );

        vdhits_.push_back(InputHit(particle, energy, minfo, pr));
      }
    }

    //================================================================
    void ExtMonFNALBoxGenerator::produce(art::Event& event) {

      std::auto_ptr<GenParticleCollection> output(new GenParticleCollection);
      std::auto_ptr<MARSInfoCollection> info(new MARSInfoCollection());
      std::auto_ptr<GenParticleMARSAssns> assns(new GenParticleMARSAssns());

      const art::ProductID particlesPID = getProductID<GenParticleCollection>(event);
      const art::EDProductGetter *particlesGetter = event.productGetter(particlesPID);

      const art::ProductID marsPID = getProductID<MARSInfoCollection>(event);
      const art::EDProductGetter *marsGetter = event.productGetter(marsPID);

      for(unsigned count=0; count < microbunchVDHitsChunkSize_; ++count) {
        const InputHit& hit = vdhits_[nextVDHit_];

        const bool acceptHitWeight(randFlat_.fire() * weightMax_ <= hit.info.weight());
        if(acceptHitWeight) {

          const double randomizedTime =
            fmod(hit.particle.time, deBuncherPeriod_)
            + (2*randFlat_.fire() - 1.)*bunchTimeHalfWidth_;

          if((cutTimeMin_ < randomizedTime) && (randomizedTime < cutTimeMax_)) {
            output->push_back(createOutputParticle(randomizedTime, hit));
            info->push_back(vdhits_[nextVDHit_].info);
            assns->addSingle(art::Ptr<GenParticle>(particlesPID, output->size()-1, particlesGetter),
                             art::Ptr<MARSInfo>(marsPID, info->size()-1, marsGetter));
          }

        }
        ++nextVDHit_ %= vdhits_.size();
      }

      event.put(output);
      event.put(info);
      event.put(assns);
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

  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::ExtMonFNALBoxGenerator);

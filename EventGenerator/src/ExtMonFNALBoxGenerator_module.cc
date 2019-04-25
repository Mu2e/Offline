// $Id: ExtMonFNALBoxGenerator_module.cc,v 1.12 2013/07/30 18:45:00 wieschie Exp $
// $Author: wieschie $
// $Date: 2013/07/30 18:45:00 $
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

#include "cetlib_except/exception.h"

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GlobalConstantsService/inc/MassCache.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
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
      typedef std::set<MARSInfo, CmpProtonId> UniqProtons;

      //================================================================
    } // namespace {}

    class ExtMonFNALBoxGenerator : public art::EDProducer {
      typedef std::vector<unsigned> InputProtons;

      int verbosityLevel_;
      std::string geomModuleLabel_;
      std::string geomInstanceName_;

      std::vector<std::string> inputFiles_;

      unsigned numPrimaryProtonsPerMicrobunch_;
      double   primaryProtonHitEfficiency_;
      double   primaryProtonStopEfficiency_;

      double cutMuonTimeMin_;
      std::vector<double> keepInBox_;

      bool randomizeMomentumDirection_;

      unsigned failCountLimit_;

      // Max MARS weigh limits for accept/reject, computed from inputs
      double hitsWeightMax_;
      double stopsWeightMax_;

      art::RandomNumberGenerator::base_engine_t& eng_;
      CLHEP::RandFlat randFlat_;
      CLHEP::RandPoissonQ randPoisson_;

      const ExtMon *extmon_;

      MassCache mc_;

      InputHits  vdhits_;
      InputProtons hitProtons_;  // for each proton*simpath in vdhits, index of its first particle in vdhits_

      InputStops mustops_;
      InputProtons stopProtons_;  // for each proton*simpath in mustops, index of its first particle in vdhits_

      void loadInputFiles();

      void addVDHits(TFile* infile);
      void assembleHitProtons();
      void computeHitsWeightMax();

      void addStoppedMuons(TFile* infile);
      void assembleStopProtons();
      void computeStopsWeightMax();

      void generateFromVDHits(const art::Event& event,
                              GenParticleCollection *output,
                              MARSInfoCollection *info,
                              GenParticleMARSAssns *assns);
      GenParticle createOutputParticle(const InputHit& hit);
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
      : EDProducer{pset}
      , verbosityLevel_(pset.get<int>("verbosityLevel"))
      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))

      , inputFiles_(pset.get<std::vector<std::string> >("inputFiles"))

      , numPrimaryProtonsPerMicrobunch_(pset.get<unsigned>("numPrimaryProtonsPerMicrobunch"))
      , primaryProtonHitEfficiency_(pset.get<double>("primaryProtonHitEfficiency"))
      , primaryProtonStopEfficiency_(pset.get<double>("primaryProtonStopEfficiency"))

      , cutMuonTimeMin_(pset.get<double>("cutMuonTimeMin"))
      , keepInBox_(pset.get<std::vector<double> >("keepInBox"))

      , randomizeMomentumDirection_(pset.get<bool>("randomizeMomentumDirection"))
      , failCountLimit_(pset.get<unsigned>("failCountLimit", 1000))

      , hitsWeightMax_(0)
      , stopsWeightMax_(0)

      , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
      , randFlat_(eng_)
      , randPoisson_(eng_)

      , extmon_()
    {
      produces<mu2e::GenParticleCollection>();
      produces<mu2e::MARSInfoCollection>();
      produces<mu2e::GenParticleMARSAssns>();

      if(inputFiles_.empty()) {
        throw cet::exception("BADCONFIG")<<"Error: no inputFiles";
      }
      if(verbosityLevel_ > 0) {
        std::cout<<"ExtMonFNALBoxGenerator: numPrimaryProtonsPerMicrobunch = "<<numPrimaryProtonsPerMicrobunch_
                 <<", will use "<<numPrimaryProtonsPerMicrobunch_ * primaryProtonHitEfficiency_<<" \"hit\" protons "
                 <<" and "<<numPrimaryProtonsPerMicrobunch_*primaryProtonStopEfficiency_<<" \"stop\" protons"
                 <<std::endl;
      }
    }

    //================================================================
    void ExtMonFNALBoxGenerator::beginRun(art::Run& run) {
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
      assembleHitProtons();
      assembleStopProtons();
      computeHitsWeightMax();
      computeStopsWeightMax();

      if(verbosityLevel_ > 0) {
        std::cout<<"ExtMonFNALBoxGenerator inputs: num hits = "<<vdhits_.size()
                 <<", hitProtons = "<<hitProtons_.size()
                 <<", num stopped muons = "<<mustops_.size()
                 <<", stopProtons = "<<stopProtons_.size()
                 <<", hitsWeightMax = "<<hitsWeightMax_
                 <<", stopsWeightMax = "<<stopsWeightMax_
                 <<std::endl;

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
    void ExtMonFNALBoxGenerator::assembleHitProtons() {
      // We can simply look at contiguous particles here
      // because particles from one proton are contiguous in
      // the MARS outputs, and this is preserved by generator
      // and analyzer/dumper up in the chain.

      hitProtons_.push_back(0);
      MARSInfo mcurrent = vdhits_[0].minfo;
      IO::G4JobInfo gcurrent = vdhits_[0].g4s1info;
      for(unsigned i=1; i < vdhits_.size(); ++i) {
        if(! (sameProtonAndSimPath(vdhits_[i].minfo, mcurrent) &&
              (gcurrent == vdhits_[i].g4s1info)
              )
           )
          {
            hitProtons_.push_back(i);
            mcurrent = vdhits_[i].minfo;
            gcurrent = vdhits_[i].g4s1info;
          }
      }
    }

    //================================================================
    void ExtMonFNALBoxGenerator::computeHitsWeightMax() {
      for(InputHits::const_iterator i=vdhits_.begin(); i!=vdhits_.end(); ++i) {
        if(hitsWeightMax_ < i->minfo.weight()) {
          hitsWeightMax_ = i->minfo.weight();
        }
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
    void ExtMonFNALBoxGenerator::assembleStopProtons() {
      // We can simply look at contiguous particles here
      // because particles from one proton are contiguous in
      // the MARS outputs, and this is preserved by generator
      // and analyzer/dumper up in the chain.

      stopProtons_.push_back(0);
      MARSInfo mcurrent = mustops_[0].minfo;
      IO::G4JobInfo gcurrent = mustops_[0].g4s1info;
      for(unsigned i=1; i < mustops_.size(); ++i) {
        if(! (sameProtonAndSimPath(mustops_[i].minfo, mcurrent) &&
              (gcurrent == mustops_[i].g4s1info)
              )
           )
          {
            stopProtons_.push_back(i);
            mcurrent = mustops_[i].minfo;
            gcurrent = mustops_[i].g4s1info;
          }
      }
    }

    //================================================================
    void ExtMonFNALBoxGenerator::computeStopsWeightMax() {
      for(InputStops::const_iterator i=mustops_.begin(); i!=mustops_.end(); ++i) {
        if(stopsWeightMax_ < i->minfo.weight()) {
          stopsWeightMax_ = i->minfo.weight();
        }
      }
    }

    //================================================================
    void ExtMonFNALBoxGenerator::produce(art::Event& event) {

      std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
      std::unique_ptr<MARSInfoCollection> info(new MARSInfoCollection());
      std::unique_ptr<GenParticleMARSAssns> assns(new GenParticleMARSAssns());

      generateFromVDHits(event, output.get(), info.get(), assns.get());
      generateFromStoppedMuons(event, output.get(), info.get(), assns.get());

      event.put(std::move(output));
      event.put(std::move(info));
      event.put(std::move(assns));
    }

    //================================================================
    void ExtMonFNALBoxGenerator::generateFromStoppedMuons(const art::Event& event,
                                                          GenParticleCollection *output,
                                                          MARSInfoCollection *info,
                                                          GenParticleMARSAssns *assns) {

      const art::ProductID particlesPID = event.getProductID<GenParticleCollection>();
      const art::EDProductGetter *particlesGetter = event.productGetter(particlesPID);

      const art::ProductID marsPID = event.getProductID<MARSInfoCollection>();
      const art::EDProductGetter *marsGetter = event.productGetter(marsPID);

      const int numStopProtons = randPoisson_.fire(numPrimaryProtonsPerMicrobunch_ * primaryProtonStopEfficiency_);
      if(verbosityLevel_ > 1) {
        std::cout<<"Generating stopped muons for "<<numStopProtons<<" input protons in this event."<<std::endl;
      }

      UniqProtons seen;
      unsigned failcount = 0;
      int count = 0;

      while(count < numStopProtons) {

        const unsigned iproton = randFlat_.fireInt(stopProtons_.size());

        // NB: cases when different particles from one proton have different
        // weights are not well defined.  Use the weight for the first
        // particle corresponding to the proton
        const bool acceptStopWeight(randFlat_.fire() * stopsWeightMax_ <= mustops_[stopProtons_[iproton]].minfo.weight());
        if(acceptStopWeight) {

          // Don't put partly correlated particles in one event
          if(seen.insert(mustops_[stopProtons_[iproton]].minfo).second) {
            ++count;  // count protons *before* applying muon decay weight
            failcount = 0;

            const unsigned start = stopProtons_[iproton];
            const unsigned end = (1 + iproton == stopProtons_.size()) ?
              mustops_.size() : stopProtons_[1+iproton];

            for(unsigned istop = start; istop < end; ++istop) {
              const InputStop& mustop = mustops_[istop];

              const double generatedMuonTime = std::max(mustop.muon.time, cutMuonTimeMin_);
              static const double tauMuPlus = 2197.; //ns, free muon
              const double decayWeight = exp((generatedMuonTime - cutMuonTimeMin_)/tauMuPlus);

              if(randFlat_.fire() <= decayWeight) {
                output->push_back(createOutputMuon(generatedMuonTime, mustop));
                info->push_back(mustop.minfo);
                assns->addSingle(art::Ptr<GenParticle>(particlesPID, output->size()-1, particlesGetter),
                                 art::Ptr<MARSInfo>(marsPID, info->size()-1, marsGetter));

              }
            } // for(muons from this proton)
          }
          else {
            if(++failcount > failCountLimit_) {
              throw cet::exception("BADINPUTS")<<"Error: failed to find an uncorrelated particle after "
                                               <<failCountLimit_<<" trials\n";
            }
          } // else-skip correlated
        } // acceptStopWeight
      } // while(count)
    }// generateFromStoppedMuons()

    //================================================================
    void ExtMonFNALBoxGenerator::generateFromVDHits(const art::Event& event,
                                                    GenParticleCollection *output,
                                                    MARSInfoCollection *info,
                                                    GenParticleMARSAssns *assns) {

      const art::ProductID particlesPID = event.getProductID<GenParticleCollection>();
      const art::EDProductGetter *particlesGetter = event.productGetter(particlesPID);

      const art::ProductID marsPID = event.getProductID<MARSInfoCollection>();
      const art::EDProductGetter *marsGetter = event.productGetter(marsPID);

      const int numHitProtons = randPoisson_.fire(numPrimaryProtonsPerMicrobunch_ * primaryProtonHitEfficiency_);
      if(verbosityLevel_ > 1) {
        std::cout<<"Generating box particles for "<<numHitProtons<<" input protons in this event."<<std::endl;
      }

      UniqProtons seen;
      unsigned failcount = 0;
      int count = 0;

      while(count < numHitProtons) {

        const unsigned iproton = randFlat_.fireInt(hitProtons_.size());

        // NB: cases when different particles from one proton have different
        // weights are not well defined.  Use the weight for the first
        // particle corresponding to the proton
        const bool acceptHitWeight(randFlat_.fire() * hitsWeightMax_ <= vdhits_[hitProtons_[iproton]].minfo.weight());
        if(acceptHitWeight) {

          // Don't put partly correlated particles in one event
          if(seen.insert(vdhits_[hitProtons_[iproton]].minfo).second) {
            ++count;
            failcount = 0;

            const unsigned start = hitProtons_[iproton];
            const unsigned end = (1 + iproton == hitProtons_.size()) ?
              vdhits_.size() : hitProtons_[1+iproton];

            for(unsigned ihit = start; ihit < end; ++ihit) {
              const InputHit& hit = vdhits_[ihit];
              output->push_back(createOutputParticle(hit));
              info->push_back(hit.minfo);
              assns->addSingle(art::Ptr<GenParticle>(particlesPID, output->size()-1, particlesGetter),
                               art::Ptr<MARSInfo>(marsPID, info->size()-1, marsGetter));

            }
          }
          else {
            if(++failcount > failCountLimit_) {
              throw cet::exception("BADINPUTS")<<"Error: failed to find an uncorrelated particle after "
                                               <<failCountLimit_<<" trials\n";
            }
          } // else-skip correlated

        } // acceptHitWeight
      } // while(count)
    } // generateFromVDHits()

    //================================================================
    GenParticle ExtMonFNALBoxGenerator::createOutputParticle(const InputHit& hit) {
      using CLHEP::Hep3Vector;

      // We want to randomize particle using flat distribution in a
      // range.  The previous step computed RMS sigma; reinterpret it
      // as a width of the flat distribution.  The choices here are
      // arbitrary. Use the customary factor:
      const double sigmaScaleFactor = sqrt(12.);

      Hep3Vector posExtMon(hit.particle.emx, hit.particle.emy, hit.particle.emz);
      do {
        if((hit.particle.vdId != VirtualDetectorId::EMFBoxSW) &&
           (hit.particle.vdId != VirtualDetectorId::EMFBoxNE)) {
          posExtMon.setX(hit.particle.emx +  hit.pr.sigmax * sigmaScaleFactor * (randFlat_.fire() - 0.5));
        }

        if((hit.particle.vdId != VirtualDetectorId::EMFBoxBottom)&&
           (hit.particle.vdId != VirtualDetectorId::EMFBoxTop)) {
          posExtMon.setY(hit.particle.emy + hit.pr.sigmay * sigmaScaleFactor * (randFlat_.fire() - 0.5));
        }

        if((hit.particle.vdId != VirtualDetectorId::EMFBoxFront)&&
           (hit.particle.vdId != VirtualDetectorId::EMFBoxBack)) {
          posExtMon.setZ(hit.particle.emz + hit.pr.sigmaz * sigmaScaleFactor * (randFlat_.fire() - 0.5));
        }
      } while(!inRange(hit.particle.vdId, posExtMon));

      const Hep3Vector posMu2e(extmon_->extMonToMu2e_position(posExtMon));

      // Randomize the momentum direction if requested.
      const Hep3Vector orig(hit.particle.mu2epx, hit.particle.mu2epy, hit.particle.mu2epz);
      Hep3Vector randomized3mom(orig);
      if(randomizeMomentumDirection_) {
        // Draw dtheta from the Rayleigh distribution
        const double dtheta = hit.pr.rSigmaML * sqrt(-2*log(randFlat_.fire()));
        const double dphi = 2*M_PI*randFlat_.fire();

        // Find a vector that is not collinear with the original particle direction
        const Hep3Vector n1 = (std::abs(orig.x()) < std::abs(orig.y())) ?
          ((std::abs(orig.x()) < std::abs(orig.z())) ? Hep3Vector(1,0,0) : Hep3Vector(0,0,1)) :
          ((std::abs(orig.x()) < std::abs(orig.y())) ? Hep3Vector(1,0,0) : Hep3Vector(0,1,0));

        // Construct a vector perpendicular to the original momentum
        const Hep3Vector perp = orig.cross(n1);

        randomized3mom.rotate(perp, dtheta);
        randomized3mom.rotate(orig, dphi);
      }

      const CLHEP::HepLorentzVector momMu2e(randomized3mom, hit.energy);

      return GenParticle(PDGCode::type(hit.particle.pdgId),
                         GenId::MARS,
                         posMu2e,
                         momMu2e,
                         hit.particle.time);
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
       if(ms.muon.stoppedInSensor) {
         posExtMon.setX( 2 * extmon_->up().planes()[0].halfSize()[0] * (randFlat_.fire() - 0.5));
         posExtMon.setY( 2 * extmon_->up().planes()[0].halfSize()[1] * (randFlat_.fire() - 0.5));
       }

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

//
// Create particle flux in the ExtMonFNAL room by randomizing
// kinematic of input particles read from a file.
// The output is in the Mu2e coordinate system.
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
#include "CLHEP/Random/RandGaussQ.h"
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
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Utilities/inc/EMFBoxIO.hh"
#include "ExtinctionMonitorFNAL/Utilities/inc/EMFRandomizationParticleDefs.hh"
#include "ExtinctionMonitorFNAL/Utilities/inc/EMFRandomizationSourceDefs.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/MassCache.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/MARSInfo.hh"
#include "MCDataProducts/inc/MARSInfoCollection.hh"
#include "MCDataProducts/inc/GenParticleMARSAssns.hh"

#include "TTree.h"
#include "TFile.h"

namespace mu2e {
  namespace ExtMonFNAL {

    using namespace Randomization;

    namespace {
      //================================================================
      struct InputParticle {
        IO::EMFRoomHit particle;
        double energy;
        MARSInfo info;
        IO::ParticleRandomization pr;

        InputParticle(const IO::EMFRoomHit& p, double e, const MARSInfo& mi, const IO::ParticleRandomization& rr)
          : particle(p), energy(e), info(mi), pr(rr)
        {}
      };

      typedef std::vector<InputParticle> InputParticles;

      //================================================================
    } // namespace {}

    class ExtMonFNALRoomGenerator : public art::EDProducer {
      int verbosityLevel_;
      std::string geomModuleLabel_;
      std::string geomInstanceName_;

      std::vector<std::string> inputFiles_;

      unsigned numProtonsPerEvent_;
      unsigned failCountLimit_;
      bool randomizeMomentumDirection_;

      // Fire flat random and pick the protons from the whole input array.
      // Make sure not to pick up partly correlated particles
      // into a single event.  We should get all particles from
      // a proton generated in one MARS run, but not particles from
      // the same proton in different runs  (as rand1 and rand2).

      InputParticles particles_;

      typedef std::vector<unsigned> InputProtons;
      InputProtons protonStart_;  // for each stat. independent proton, index of its first particle in particles_

      art::RandomNumberGenerator::base_engine_t& eng_;
      CLHEP::RandFlat randFlat_;
      CLHEP::RandGaussQ randGauss_;

      const ProtonBeamDump *dump_;
      const ExtMon *extmon_;
      SourcePlaneGeom srcGeom_;

      MassCache mc_;

      void loadInputFiles();

      void addInputParticles(TFile* infile);

      void assembleProtons();

      GenParticle createOutputParticle(InputParticles::size_type ip);

    public:
      explicit ExtMonFNALRoomGenerator(const fhicl::ParameterSet& pset);
      virtual void produce(art::Event& event);
      virtual void beginRun(art::Run& run);
    };

    //================================================================
    ExtMonFNALRoomGenerator::ExtMonFNALRoomGenerator(const fhicl::ParameterSet& pset)
      : EDProducer{pset}
      , verbosityLevel_(pset.get<int>("verbosityLevel"))
      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))
      , inputFiles_(pset.get<std::vector<std::string> >("inputFiles"))

      , numProtonsPerEvent_(pset.get<unsigned>("numProtonsPerEvent"))
      , failCountLimit_(pset.get<unsigned>("failCountLimit", 1000))
      , randomizeMomentumDirection_(pset.get<bool>("randomizeMomentumDirection"))

      , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
      , randFlat_(eng_)
      , randGauss_(eng_)

      , dump_()
      , extmon_()

      , srcGeom_(pset.get<fhicl::ParameterSet>("srcGeom"))
    {
      produces<mu2e::GenParticleCollection>();
      produces<mu2e::MARSInfoCollection>();
      produces<mu2e::GenParticleMARSAssns>();

      if(inputFiles_.empty()) {
        throw cet::exception("BADCONFIG")<<"Error: no inputFiles\n";
      }
    }

    //================================================================
    void ExtMonFNALRoomGenerator::beginRun(art::Run& run) {

      if(!geomModuleLabel_.empty()) {
        art::Handle<ProtonBeamDump> dump;
        run.getByLabel(geomModuleLabel_, geomInstanceName_, dump);
        dump_ = &*dump;

        art::Handle<ExtMon> extmon;
        run.getByLabel(geomModuleLabel_, geomInstanceName_, extmon);
        extmon_ = &*extmon;
      }
      else {
        GeomHandle<ProtonBeamDump> dump;
        dump_ = &*dump;

        GeomHandle<ExtMon> extmon;
        extmon_ = &*extmon;
      }

      loadInputFiles();

      if(particles_.empty()) {
        throw cet::exception("BADINPUTs")<<"Error: no particles loaded.\n";
      }

      assembleProtons();

      std::cout<<"ExtMonFNALRoomGenerator: loaded "
               <<particles_.size()<<" particles for "
               <<protonStart_.size()<<" protons. "
               <<" randomizeMomentumDirection = "<<randomizeMomentumDirection_
               <<std::endl;
    }

    //================================================================
    void ExtMonFNALRoomGenerator::loadInputFiles() {
      art::ServiceHandle<art::TFileService> tfs;

      for(unsigned i=0; i<inputFiles_.size(); ++i) {

        const std::string resolvedFileName = ConfigFileLookupPolicy()(inputFiles_[i]);
        if(verbosityLevel_ > 0) {
          std::cout<<"ExtMonFNALRoomGenerator: reading input file "<<resolvedFileName<<std::endl;
        }

        TFile *infile = tfs->make<TFile>(resolvedFileName.c_str(), "READ");
        addInputParticles(infile);

      }
    } // loadInputFiles()

    //================================================================
    void ExtMonFNALRoomGenerator::addInputParticles(TFile* infile) {
      const std::string hitTreeName_("EMFRoomFluxAnalyzer/roomhits");

      TTree *nt = dynamic_cast<TTree*>(infile->Get(hitTreeName_.c_str()));
      if(!nt) {
        throw cet::exception("BADINPUT")<<"Could not get tree \""<<hitTreeName_
                                        <<"\" from file \""<<infile->GetName()
                                        <<"\"\n";
      }

      const Long64_t nTreeEntries = nt->GetEntries();

      IO::EMFRoomHit particle;
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

        particles_.push_back(InputParticle(particle, energy, minfo, pr));
      }
    }

    //================================================================
    void ExtMonFNALRoomGenerator::assembleProtons() {
      // We can simply look at contiguous particles here
      // because particles from one proton are contiguous in
      // the MARS outputs, and this is preserved by generator
      // and analyzer/dumper up in the chain.

      protonStart_.push_back(0);
      MARSInfo current = particles_[0].info;
      for(unsigned i=1; i < particles_.size(); ++i) {
        if(!sameProtonAndSimPath(particles_[i].info, current)) {
          protonStart_.push_back(i);
          current = particles_[i].info;
        }
      }
    }

    //================================================================
    void ExtMonFNALRoomGenerator::produce(art::Event& event) {

      std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
      std::unique_ptr<MARSInfoCollection> info(new MARSInfoCollection());
      std::unique_ptr<GenParticleMARSAssns> assns(new GenParticleMARSAssns());

      const art::ProductID particlesPID = event.getProductID<GenParticleCollection>();
      const art::EDProductGetter *particlesGetter = event.productGetter(particlesPID);

      const art::ProductID marsPID = event.getProductID<MARSInfoCollection>();
      const art::EDProductGetter *marsGetter = event.productGetter(marsPID);

      typedef std::set<MARSInfo, CmpProtonId> UniqProtons;

      UniqProtons seen;

      // Process requested number of protons
      unsigned failcount = 0;
      unsigned count = 0;
      while(count < numProtonsPerEvent_) {

        const unsigned iproton = randFlat_.fireInt(protonStart_.size());

        // Don't put partly correlated particles in one event
        if(seen.insert(particles_[protonStart_[iproton]].info).second) {
          ++count;
          failcount = 0;

          const unsigned start = protonStart_[iproton];
          const unsigned end = (1 + iproton == protonStart_.size()) ?
            particles_.size() : protonStart_[1+iproton];

          for(unsigned iparticle = start; iparticle < end; ++iparticle) {
            output->push_back(createOutputParticle(iparticle));
            info->push_back(particles_[iparticle].info);
            assns->addSingle(art::Ptr<GenParticle>(particlesPID, output->size()-1, particlesGetter),
                             art::Ptr<MARSInfo>(marsPID, info->size()-1, marsGetter));
          }
        }
        else {
          if(++failcount > failCountLimit_) {
            throw cet::exception("BADINPUTS")<<"Error: failed to find an uncorrelated particle after "
                                             <<failCountLimit_<<" trials\n";
          }
        }

      } // while(doing requested number of protons)

      event.put(std::move(output));
      event.put(std::move(info));
      event.put(std::move(assns));

    } // produce()

    //================================================================
    GenParticle ExtMonFNALRoomGenerator::createOutputParticle(InputParticles::size_type ip) {
      using CLHEP::Hep3Vector;

      const IO::EMFRoomHit& pp = particles_[ip].particle;
      const IO::ParticleRandomization& rr = particles_[ip].pr;

      const double sigmax = (pp.srcType == SourceSouthWest)||(pp.srcType == SourceNorthEast) ?
        0 : rr.sigmax;

      const double sigmay = (pp.srcType == SourceFloor)||(pp.srcType == SourceCeiling) ?
        0 : rr.sigmay;

      const double sigmaz = (pp.srcType == SourceFront)||(pp.srcType == SourceBack)||(pp.srcType == SourceSignal) ?
        0 : rr.sigmaz;

      // We want to randomize particle using flat distribution in a
      // range.  The previous step computed RMS sigma; reinterpret it
      // as a width of the flat distribution.  The choices here are
      // arbitrary. Use the customary factor:
      const double sigmaScaleFactor = sqrt(12.);

      Hep3Vector posDump(pp.dumpx, pp.dumpy, pp.dumpz);
      do{
        posDump.setX(pp.dumpx + sigmax * sigmaScaleFactor * (randFlat_.fire() - 0.5));
        posDump.setY(pp.dumpy + sigmay * sigmaScaleFactor * (randFlat_.fire() - 0.5));
        posDump.setZ(pp.dumpz + sigmaz * sigmaScaleFactor * (randFlat_.fire() - 0.5));
      } while(!inRange(SourceType(pp.srcType), posDump, srcGeom_, *dump_, *extmon_));

      const Hep3Vector posMu2e(dump_->beamDumpToMu2e_position(posDump));

      // Randomize the momentum direction if requested.
      const Hep3Vector orig(pp.mu2epx, pp.mu2epy, pp.mu2epz);
      Hep3Vector randomized3mom(orig);
      if(randomizeMomentumDirection_) {
        // Draw dtheta from the Rayleigh distribution
        const double dtheta = rr.rSigmaML * sqrt(-2*log(randFlat_.fire()));
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

      const CLHEP::HepLorentzVector momMu2e(randomized3mom, particles_[ip].energy);

      return GenParticle(PDGCode::type(pp.pdgId), GenId::MARS, posMu2e, momMu2e, pp.time);
    }

    //================================================================
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::ExtMonFNALRoomGenerator);

// $Id: ExtMonFNALRoomGenerator_module.cc,v 1.4 2012/11/01 23:39:49 gandr Exp $
// $Author: gandr $
// $Date: 2012/11/01 23:39:49 $
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

#include "cetlib/exception.h"

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "GeneralUtilities/inc/KNearestNeighbors.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/MassCache.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/MARSInfo.hh"
#include "MCDataProducts/inc/MARSInfoCollection.hh"
#include "MCDataProducts/inc/GenParticleMARSAssns.hh"

// FIXME: move utils into an appropriate package
#include "Sources/inc/ExtMonFNALMARSUtils.hh"

#include "TH1.h"
#include "TH2.h"

//#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
#define AGDEBUG(stuff)


namespace mu2e {
  namespace ExtMonFNAL {

    namespace {
      //================================================================
      // src % 3 == 0: don't smear dumpz
      // src % 3 == 1: don't smear dumpx
      // src % 3 == 2: don't smear dumpy
      enum SourceType {
        SourceFront = 0,
        SourceSouthWest,
        SourceFloor,
        SourceBack,
        SourceNorthEast,
        SourceCeiling,

        SourceSignal, // SourceFront particles which are also signal candidates

        NUM_SOURCES
      };

      //================================================================
      // These are in beam dump coordinates
      // The numbers are used to
      //
      //    1) assign input MARS particles to a source plane
      //    2) establish boundaries which (non-signal) smearing should not cross

      struct SourcePlaneGeom {
        double zFront;
        double xSW;
        double yFloor;
        double zBack;
        double xNE;
        double yCeiling;
        SourcePlaneGeom(double x1, double x2, double y1, double y2, double z1, double z2)
          : zFront(z1), xSW(x1), yFloor(y1), zBack(z2), xNE(x2), yCeiling(y2)
        {}
      };

      std::ostream& operator<<(std::ostream&os, const SourcePlaneGeom& gm) {
        return
          os<<"SourcePlaneGeom(zFront="<<gm.zFront
            <<", xSW="<<gm.xSW
            <<", zBack="<<gm.zBack
            <<", yFloor="<<gm.yFloor
            <<", xNE="<<gm.xNE
            <<", yCeiling="<<gm.yCeiling
            <<" )"
          ;
      }

      //================================================================
      enum ParticleType {
        ELECTRON, MUON, PROTON, OTHER_CHARGED,
        NEUTRON, GAMMA, OTHER_NEUTRAL,
        NUM_PARTICLE_TYPES
      };

      //================================================================
      struct InputParticle {
        PDGCode::type pdgId;
        CLHEP::HepLorentzVector momMu2e;
        CLHEP::Hep3Vector posDump;
        double time;
        MARSInfo info;
        SourceType srcType;
        InputParticle() : pdgId(), time(), srcType(NUM_SOURCES) {}
      };

      typedef std::vector<InputParticle> InputParticles;

      std::ostream& operator<<(std::ostream& os, const InputParticle& part) {
        return os<<"Particle(posDump="<<part.posDump<<", momMu2e="<<part.momMu2e<<")";
      }

      //================================================================
      struct ParticleRandomization {
        // We have some redundant data here for convenience of
        // studying the randomization approach

        // Variance of (X_neighbor - X)
        CLHEP::HepSymMatrix vpos;

        double sigmax;
        double sigmay;
        double sigmaz;
        double correlationCoefficient;  // between the two coordinates relevant for the source plane

        // "Gaussian" distribution of neighbor directions (nx', ny')
        // w.r.t this direction leads to Rayleigh distribution of
        // r = sin(theta'):  pdf ~ r/sigma^2 exp(-r^2/(2*sigma^2))
        //
        // Rayleigh mean = sigma * sqrt(pi/2)
        // Rayleigh variance = sigma^2 * (4-pi)/2
        //
        // So we can estimate sigma both ways and compare.

        double rSigmaFromMean;
        double rSigmaFromVariance;

        // ML estimate of sigma is sqrt(\sum_1^N x^2_i /(2N))
        double rSigmaML;

        ParticleRandomization()
          : sigmax(), sigmay(), sigmaz(), correlationCoefficient()
          , rSigmaFromMean(), rSigmaFromVariance(), rSigmaML()
        {}
      };

      typedef std::vector<ParticleRandomization> ParticleRandomizations;

      //================================================================
      class HistRandomization {
        TH2 *angleSigmasVM;
        TH2 *angleSigmasVL;
        //TH1 * sigmax;  ..; TH1 *corrxy; ..

        TH1 *sigmax;
        TH1 *sigmay;
        TH1 *sigmaz;
        TH1 *correlationCoefficient;

      public:
        HistRandomization(SourceType st, ParticleType pt);
        void fill(const ParticleRandomization& pr);
      };

      HistRandomization::HistRandomization(SourceType st, ParticleType pt)
        : angleSigmasVM(), angleSigmasVL()
      {
        art::ServiceHandle<art::TFileService> tfs;

        std::ostringstream os;
        os<<"st"<<st<<"_pt"<<pt;
        art::TFileDirectory tfdir = tfs->mkdir(os.str());

        angleSigmasVM = tfdir.make<TH2D>("angleSigmasVM", "Rayleigh sigma from variance vs from mean",
                                         100, 0., 1., 100, 0., 1.);

        angleSigmasVM->SetOption("colz");

        angleSigmasVL = tfdir.make<TH2D>("angleSigmasVL", "Rayleigh sigma from variance vs ML estimate",
                                         100, 0., 1., 100, 0., 1.);
        angleSigmasVL->SetOption("colz");

        sigmax = tfdir.make<TH1D>("sigmax", "sigma x", 1000, 0., 500.);
        sigmay = tfdir.make<TH1D>("sigmay", "sigma y", 1000, 0., 500.);
        sigmaz = tfdir.make<TH1D>("sigmaz", "sigma z", 1000, 0., 500.);
        correlationCoefficient = tfdir.make<TH1D>("correlation", "correlation", 101, -1., 1.);

      }

      void HistRandomization::fill(const ParticleRandomization& pr) {
        angleSigmasVM->Fill(pr.rSigmaFromMean, pr.rSigmaFromVariance);
        angleSigmasVL->Fill(pr.rSigmaML, pr.rSigmaFromVariance);

        sigmax->Fill(pr.sigmax);
        sigmay->Fill(pr.sigmay);
        sigmaz->Fill(pr.sigmaz);
        correlationCoefficient->Fill(pr.correlationCoefficient);
      }

      //================================================================
      template<class Point>
      std::ostream& operator<<(std::ostream& os, const KNearestNeighbors<Point>& knn) {
        os<<"KNearestNeighbors: {\n";
        for(unsigned i=0; i<knn.size(); ++i) {
          os<<"    particle "<<i<<": ";
          for(unsigned j=0; j<knn[i].size(); ++j) {
            const typename KNearestNeighbors<Point>::Entry ee = knn[i][j];
            os<<*ee.point<<" "<<ee.distance<<", ";
            //os<<knn[i][j]<<", ";
          }
          os<<"\n";
        }
        os<<"};\n";
        return os;
      }

      //================================================================
    } // namespace {}

    class ExtMonFNALRoomGenerator : public art::EDProducer {
      int verbosityLevel_;
      bool doHistograms_;
      std::vector<std::string> inputFiles_;

      unsigned numReadProtonsPerEvent_;
      double fracAcceptedProtons_;

      std::string geomModuleLabel_;
      std::string geomInstanceName_;

      art::RandomNumberGenerator::base_engine_t& eng_;
      CLHEP::RandFlat randFlat_;
      CLHEP::RandGaussQ randGauss_;

      const ProtonBeamDump *dump_;
      const ExtMon *extmon_;
      SourcePlaneGeom srcGeom_;
      double srcPosTolerance_;
      double signalHalfdx_;
      double signalHalfdy_;

      unsigned numNeighbors_;
      unsigned minSourceGroupStatistics_;

      GlobalConstantsHandle<ParticleDataTable> pdt_;
      MassCache mc_; // this looks redundant with pdt_, but we only use pdt_ to get charge of "rare" particles

      MARSMu2eConverter cnv_;

      InputParticles particles_;
      ParticleRandomizations pr_; // Index synced with particles_
      InputParticles::size_type numInputProtons_;
      InputParticles::size_type nextParticle_;

      InputParticles::size_type indexInParticles(const InputParticle *p) {
        return p - &particles_[0];
      }

      // inputs arranged by [SourceType][ParticleType][particleIndex]
      typedef std::vector<const InputParticle*> RandomizationGroup;
      typedef std::vector<std::vector<RandomizationGroup> > GroupedInputs;
      GroupedInputs grouped_;

      void loadInputFiles();

      // initGroups() needs to be a separate step from
      // loadInputFiles() because pointers into initial particles_ are
      // invalidated when more particles are loaded.
      void initGroups();

      SourceType classifySource(const CLHEP::Hep3Vector& posDump, const CLHEP::Hep3Vector& posMu2e);
      bool isSignal(const CLHEP::Hep3Vector& posMu2e);


      ParticleType classifyParticleType(PDGCode::type pdgId);
      double getCharge(PDGCode::type pdgId);

      void printSrcGroups();
      void mergeLowStatisticSrcGroups();
      void computeParticleRandomizations();

      ParticleRandomization computeParticleRandomization(const InputParticle *p,
                                                         const KNearestNeighbors<const InputParticle*>::Points& neighbors);

      typedef std::map<std::pair<SourceType,ParticleType>, HistRandomization> HistMap;
      HistMap histRandomization_;
      void fillRandomizationHistograms(SourceType st, ParticleType pt, const ParticleRandomization& pr);

      GenParticle createOutputParticle(InputParticles::size_type ip);
      bool inRange(SourceType st, const CLHEP::Hep3Vector& posDump);

    public:
      explicit ExtMonFNALRoomGenerator(const fhicl::ParameterSet& pset);
      virtual void produce(art::Event& event);
      virtual void beginRun(art::Run& run);
    };

    //================================================================
    ExtMonFNALRoomGenerator::ExtMonFNALRoomGenerator(const fhicl::ParameterSet& pset)
      : verbosityLevel_(pset.get<int>("verbosityLevel"))
      , doHistograms_(pset.get<bool>("doHistograms"))
      , inputFiles_(pset.get<std::vector<std::string> >("inputFiles"))

      , numReadProtonsPerEvent_(pset.get<unsigned>("numReadProtonsPerEvent"))
      , fracAcceptedProtons_(pset.get<double>("fracAcceptedProtons"))

      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))

      , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
      , randFlat_(eng_)
      , randGauss_(eng_)

      , dump_()
      , extmon_()

      , srcGeom_(pset.get<double>("xSW"),
                 pset.get<double>("xNE"),
                 pset.get<double>("yFloor"),
                 pset.get<double>("yCeiling"),
                 pset.get<double>("zFront"),
                 pset.get<double>("zBack")
                 )

      , srcPosTolerance_(pset.get<double>("srcPositionTolerance"))

        // Signal box in ExtMon coordinates is [-dx, +dx]*[-dy,+dy]
      , signalHalfdx_(pset.get<double>("signalHalfdx"))
      , signalHalfdy_(pset.get<double>("signalHalfdy"))

      , numNeighbors_(pset.get<unsigned>("numNeighbors"))
      , minSourceGroupStatistics_(pset.get<unsigned>("minSourceGroupStatistics"))

      , numInputProtons_()
      , nextParticle_()

      , grouped_(NUM_SOURCES, std::vector<RandomizationGroup>(NUM_PARTICLE_TYPES))
    {
      produces<mu2e::GenParticleCollection>();
      produces<mu2e::MARSInfoCollection>();
      produces<mu2e::GenParticleMARSAssns>();

      if(inputFiles_.empty()) {
        throw cet::exception("BADCONFIG")<<"Error: no inputFiles";
      }

      if( (fracAcceptedProtons_ <= 0.) || (1. < fracAcceptedProtons_) ) {
        throw cet::exception("BADCONFIG")
          <<"Error: got fracAcceptedProtons_ = "<<fracAcceptedProtons_
          <<".  The parameter must be in (0,1]\n";
      }

      if(srcGeom_.xSW >= srcGeom_.xNE) {
        throw cet::exception("BADCONFIG")
          <<"Error: srcGeom_.xSW ("<<srcGeom_.xSW<<") >= srcGeom_.xNE ("<<srcGeom_.xNE<<")";
      }

      if(srcGeom_.yFloor >= srcGeom_.yCeiling) {
        throw cet::exception("BADCONFIG")
          <<"Error: srcGeom_.yFloor ("<<srcGeom_.yFloor<<") >= srcGeom_.yCeiling ("<<srcGeom_.yCeiling<<")";
      }

      if(srcGeom_.zBack >= srcGeom_.zFront) {
        throw cet::exception("BADCONFIG")
          <<"Error: srcGeom_.zBack ("<<srcGeom_.zBack<<") >= srcGeom_.zFront ("<<srcGeom_.zFront<<")";
      }

      if(minSourceGroupStatistics_ < numNeighbors_) {
        throw cet::exception("BADCONFIG")
          <<"Error: minSourceGroupStatistics ("<<minSourceGroupStatistics_<<") < numNeighbors ("<<numNeighbors_<<")";
      }

      if(numNeighbors_ < 2) {
        throw cet::exception("BADCONFIG")
          <<"Error: numNeighbors="<<numNeighbors_<<", should be > 1";
      }

    }

    //================================================================
    void ExtMonFNALRoomGenerator::beginRun(art::Run& run) {
      const ExtMonFNALBuilding *building(0);

      if(!geomModuleLabel_.empty()) {
        art::Handle<ProtonBeamDump> dump;
        run.getByLabel(geomModuleLabel_, geomInstanceName_, dump);
        dump_ = &*dump;

        art::Handle<ExtMon> extmon;
        run.getByLabel(geomModuleLabel_, geomInstanceName_, extmon);
        extmon_ = &*extmon;

        art::Handle<ExtMonFNALBuilding> bh;
        run.getByLabel(geomModuleLabel_, geomInstanceName_, bh);
        building = &*bh;
      }
      else {
        GeomHandle<ProtonBeamDump> dump;
        dump_ = &*dump;

        GeomHandle<ExtMon> extmon;
        extmon_ = &*extmon;

        GeomHandle<ExtMonFNALBuilding> bh;
        building = &*bh;
      }

      const double yFloor = dump_->mu2eToBeamDump_position(CLHEP::Hep3Vector(0,building->roomInsideYmin(),0)).y();
      const double yCeiling = dump_->mu2eToBeamDump_position(CLHEP::Hep3Vector(0,building->roomInsideYmax(),0)).y();

      const double zFront = dump_->mu2eToBeamDump_position(building->coll2ShieldingCenterInMu2e()).z()
        - building->coll2ShieldingHalfSize()[2]
        ;

      const double xSW =  building->coll2ShieldingDumpXmin();
      const double xNE =  building->coll2ShieldingDumpXmax();

      if(verbosityLevel_ > 0) {
        std::cout<<"ExtMonFNALRoomGenerator INFO: src zFront offset into room   = "<<(zFront - srcGeom_.zFront)<<std::endl;
        std::cout<<"ExtMonFNALRoomGenerator INFO: src yFloor offset into room   = "<<(srcGeom_.yFloor - yFloor)<<std::endl;
        std::cout<<"ExtMonFNALRoomGenerator INFO: src yCeiling offset into room = "<<(yCeiling - srcGeom_.yCeiling)<<std::endl;
        std::cout<<"ExtMonFNALRoomGenerator INFO: src xSW offset into room      = "<<(srcGeom_.xSW - xSW)<<std::endl;
        std::cout<<"ExtMonFNALRoomGenerator INFO: src xNE offset into room      = "<<(xNE - srcGeom_.xNE)<<std::endl;
      }

      //----------------------------------------------------------------
      loadInputFiles();

      if(!numReadProtonsPerEvent_) {
        numReadProtonsPerEvent_ = numInputProtons_;
      }
      if(verbosityLevel_ > 0) {
        std::cout<<"Loaded "<<particles_.size()<<" input particles from "<<numInputProtons_<<" protons"<<std::endl;
        std::cout<<"Considering "<<numReadProtonsPerEvent_<<" input protons per event and accepting frac="
                 <<fracAcceptedProtons_<<" of them"
                 <<std::endl;
      }

      initGroups();


      if(verbosityLevel_ > 0) {
        std::cout<<"Summary of inputs before regrouping:"<<std::endl;
        printSrcGroups();
      }

      //----------------------------------------------------------------
      mergeLowStatisticSrcGroups();

      if(verbosityLevel_ > 0) {
        std::cout<<"Summary of inputs after regrouping:"<<std::endl;
        printSrcGroups();
      }

      //----------------------------------------------------------------
      computeParticleRandomizations();

      //----------------------------------------------------------------
    }

    //================================================================
    void ExtMonFNALRoomGenerator::printSrcGroups() {
      for(unsigned st = 0; st < NUM_SOURCES; ++st) {
        for(unsigned pt = 0; pt < NUM_PARTICLE_TYPES; ++pt) {
          std::cout<<"    source "<<st<<", particle type "<<pt
                   <<", num particles "<<grouped_[st][pt].size()
            ;
          if(false) {
            std::cout<<", indexes = { ";
            for(unsigned i=0; i<grouped_[st][pt].size(); ++i) {
              std::cout<<indexInParticles(grouped_[st][pt][i])<<" ";
            }
            std::cout<<" }";
          }
          std::cout<<std::endl;
        }
        std::cout<<std::endl;
      }
    }

    //================================================================
    void ExtMonFNALRoomGenerator::loadInputFiles() {
      int currentProtonNumber(-1);
      for(unsigned i=0; i<inputFiles_.size(); ++i) {
        const std::string resolvedFileName = ConfigFileLookupPolicy()(inputFiles_[i]);
        if(verbosityLevel_ > 0) {
          std::cout<<"ExtMonFNALRoomGenerator: reading input file "<<resolvedFileName<<std::endl;
        }

        const int marsSubRunNumber = 0;
        const int marsRunNumber = 99;

        std::ifstream infile(resolvedFileName.c_str());
        MARSParticle mp;
        while(readMARSLine(infile, mp)) {
          InputParticle part;

          part.pdgId = cnv_.marsToMu2eParticleCode(mp.pid);

          const double mass = mc_.mass(part.pdgId);
          const double energy = mass + cnv_.marsToMu2eEnergy(mp.kineticEnergy);
          const double p3mag = sqrt((energy-mass)*(energy+mass));

          part.momMu2e =  CLHEP::HepLorentzVector(mp.dcx * p3mag,
                                                  mp.dcy * p3mag,
                                                  mp.dcz * p3mag,
                                                  energy
                                                  );

          const CLHEP::Hep3Vector posMu2e(cnv_.marsToMu2ePosition(mp.x, mp.y, mp.z));
          part.posDump = dump_->mu2eToBeamDump_position(posMu2e);
          part.time = cnv_.marsToMu2eTime(mp.tof);
          part.info = MARSInfo(mp.weight, mp.protonNumber, marsSubRunNumber, marsRunNumber);
          part.srcType = classifySource(part.posDump, posMu2e);

          if(part.info.protonNumber() != currentProtonNumber) {
            currentProtonNumber = part.info.protonNumber();
            ++numInputProtons_;
          }

          particles_.push_back(part);

        } // while(read one file)
      } // for(files)

    } // loadInputFiles()

    //================================================================
    void ExtMonFNALRoomGenerator::initGroups() {
      for(unsigned i=0; i<particles_.size(); ++i) {
        const InputParticle& part = particles_[i];
        SourceType st = part.srcType;
        ParticleType pt = classifyParticleType(part.pdgId);
        grouped_[st][pt].push_back(&part);
        AGDEBUG("grouped_[st][pt].back() index="<<indexInParticles(grouped_[st][pt].back()));
      }
    }

    //================================================================
    SourceType ExtMonFNALRoomGenerator::classifySource(const CLHEP::Hep3Vector& posDump,
                                                       const CLHEP::Hep3Vector& posMu2e)
    {
      // The order of the tests affects the result only for corner
      // cases, where classification is ambiguous.

      if(std::abs(posDump.z() - srcGeom_.zFront) < srcPosTolerance_) {
        // check for Signal
        if(isSignal(posMu2e)) {
          return SourceSignal;
        }
        return SourceFront;
      }

      if(std::abs(posDump.x() - srcGeom_.xSW) < srcPosTolerance_) return SourceSouthWest;
      if(std::abs(posDump.y() - srcGeom_.yFloor) < srcPosTolerance_) return SourceFloor;
      if(std::abs(posDump.x() - srcGeom_.xNE) < srcPosTolerance_) return SourceNorthEast;
      if(std::abs(posDump.y() - srcGeom_.yCeiling) < srcPosTolerance_) return SourceCeiling;

      if(std::abs(posDump.z() - srcGeom_.zBack) < srcPosTolerance_) return SourceBack;

      throw cet::exception("BADINPUTS")
        <<"Error: failed to assign input posDump = "<<posDump<<" to an input source plane\n";
      //return NUM_SOURCES;
    }

    //================================================================
    bool ExtMonFNALRoomGenerator::isSignal(const CLHEP::Hep3Vector& posMu2e) {
      const CLHEP::Hep3Vector posExtMon = extmon_->mu2eToExtMon_position(posMu2e);
      return
        (std::abs(posExtMon.x()) < signalHalfdx_)&&
        (std::abs(posExtMon.y()) < signalHalfdy_);
    }

    //================================================================
    ParticleType ExtMonFNALRoomGenerator::classifyParticleType(PDGCode::type pdgId) {
      if(std::abs(pdgId)==11) {
        return ELECTRON;
      }
      if(std::abs(pdgId)==13) {
        return MUON;
      }
      if(pdgId==2212) {
        return PROTON;
      }
      else if(pdgId == 2112) {
        return NEUTRON;
      }
      else if(pdgId == 22) {
        return GAMMA;
      }
      else if(std::abs(getCharge(pdgId)) > 0.5) {
        return OTHER_CHARGED;
      }
      else {
        return OTHER_NEUTRAL;
      }
    }

    //================================================================
    double ExtMonFNALRoomGenerator::getCharge(PDGCode::type pdgId) {
      ParticleDataTable::maybe_ref info = pdt_->particle(pdgId);

      // Particles unknown to PDT are ions
      // Default ion charge:
      int charge(1); // deuterium

      if(!info.isValid()) {
        std::cout<<"ExtMonFNALRoomGenerator: no valid PDG info for pdgId = "<<pdgId<<", using charge = "<<charge<<std::endl;
      }
      else {
        charge = info.ref().charge();
      }

      return charge;
    }

    //================================================================
    void ExtMonFNALRoomGenerator::mergeLowStatisticSrcGroups() {
      for(unsigned st = 0; st < NUM_SOURCES; ++st) {

        // We merge particles types with only a few particles
        // into the highest statistics source on the same plane
        // Find the dst group.
        unsigned imax = 0;
        for(unsigned pt = 0; pt < NUM_PARTICLE_TYPES; ++pt) {
          if(grouped_[st][pt].size() > grouped_[st][imax].size()) {
            imax = pt;
          }
        }

        for(unsigned pt = 0; pt < NUM_PARTICLE_TYPES; ++pt) {
          if(!grouped_[st][pt].empty() &&
             (pt != imax) &&
             (grouped_[st][pt].size() < minSourceGroupStatistics_)) {

            if(verbosityLevel_ > 1) {
              std::cout<<"Moving particles for plane "<<st<<" from type "<<pt<<" to type "<<imax<<std::endl;
            }

            std::copy(grouped_[st][pt].begin(), grouped_[st][pt].end(),
                      std::back_inserter(grouped_[st][imax]));

            grouped_[st][pt].clear();
          }
        }

        // Verify that statistics are good after regrouping
        for(unsigned pt = 0; pt < NUM_PARTICLE_TYPES; ++pt) {
          if(!grouped_[st][pt].empty() &&
             (grouped_[st][pt].size() < minSourceGroupStatistics_)) {

            throw cet::exception("BADINPUTS")
              <<"Error: failed to regroup particles on source plane = "<<st
              <<" to satisfy minSourceGroupStatistics="<<minSourceGroupStatistics_
              <<": achieved stat="<<grouped_[st][pt].size()
              <<"\n";
          }
        }

      } // for(sources st)
    } // mergeLowStatisticSrcGroups()

    //================================================================
    class Metric {
    public:
      virtual ~Metric() {}
      virtual double operator()(const InputParticle *a, const InputParticle *b) const =0;
    };

    class MetricXY : virtual public Metric {
    public:
      virtual double operator()(const InputParticle *a, const InputParticle *b) const {
        const double dx = a->posDump.x() - b->posDump.x();
        const double dy = a->posDump.y() - b->posDump.y();
        return sqrt(dx*dx + dy*dy);
      }
    };

    class MetricYZ : virtual public Metric {
    public:
      virtual double operator()(const InputParticle *a, const InputParticle *b) const {
        const double dy = a->posDump.y() - b->posDump.y();
        const double dz = a->posDump.z() - b->posDump.z();
        return sqrt(dy*dy + dz*dz);
      }
    };

    class MetricZX : virtual public Metric {
    public:
      virtual double operator()(const InputParticle *a, const InputParticle *b) const {
        const double dx = a->posDump.x() - b->posDump.x();
        const double dz = a->posDump.z() - b->posDump.z();
        return sqrt(dx*dx + dz*dz);
      }
    };

    void ExtMonFNALRoomGenerator::computeParticleRandomizations() {
      static MetricXY mxy; // gcc 4.6.1 does not like "static const" here
      static MetricYZ myz; // gcc 4.6.1 does not like "static const" here
      static MetricZX mzx; // gcc 4.6.1 does not like "static const" here
      static const Metric *dist[] = { &mxy, &myz, &mzx };

      pr_.resize(particles_.size());

      for(unsigned st = 0; st < NUM_SOURCES; ++st) {

        for(unsigned pt = 0; pt < NUM_PARTICLE_TYPES; ++pt) {

          const RandomizationGroup& group = grouped_[st][pt];

          if(!group.empty()) {

            KNearestNeighbors<const InputParticle*>
              neighbors(numNeighbors_, group, *dist[st % 3]);

            if(verbosityLevel_ > 1) {
              std::cout<<"ExtMonFNALRoomGenerator: Neighbors dump begin"<<std::endl;
              std::cout<<neighbors<<std::endl;
              std::cout<<"ExtMonFNALRoomGenerator: Neighbors dump end"<<std::endl;
            }

            for(unsigned i=0; i<group.size(); ++i) {
              InputParticles::size_type ip = indexInParticles(group[i]);
              pr_[ip] = computeParticleRandomization(group[i], neighbors[i]);
              if(doHistograms_) {
                fillRandomizationHistograms(SourceType(st), ParticleType(pt), pr_[ip]);
              }
            }

          }
        } // for(particle types pt)
      } // for(sources st)
    } // computeParticleRandomizations()

    //================================================================
    ParticleRandomization
    ExtMonFNALRoomGenerator::computeParticleRandomization(const InputParticle *particle,
                                                          const KNearestNeighbors<const InputParticle*>::Points& neighbors)
    {
      double sumSinTheta(0);
      double sumSin2Theta(0);
      CLHEP::HepMatrix tmp(3,3);
      for(unsigned i=0; i<neighbors.size(); ++i) {
        CLHEP::HepVector dr(3);
        dr[0] = particle->posDump.x() - neighbors[i].point->posDump.x();
        dr[1] = particle->posDump.y() - neighbors[i].point->posDump.y();
        dr[2] = particle->posDump.z() - neighbors[i].point->posDump.z();
        tmp += dr*dr.T();

        const double sin2Theta =
          1. - static_cast<const CLHEP::Hep3Vector&>(particle->momMu2e).cos2Theta(neighbors[i].point->momMu2e);

        sumSin2Theta +=  sin2Theta;
        sumSinTheta += sqrt(sin2Theta);
      }

      const double n(neighbors.size());

      ParticleRandomization res;
      // Note  that here we divide by n rather than (n-1) because of the hypothesis mean==0
      // (we do symmetric smearing in coordinates, which does not bias the mean).
      res.vpos.assign(tmp/n);

      const double sinThetaMean = sumSinTheta/n;
      const double sinThetaVariance = (sumSin2Theta - n * sinThetaMean * sinThetaMean)/(n - 1);

      res.rSigmaFromMean = sqrt(2./M_PI) * sinThetaMean;
      res.rSigmaFromVariance = sqrt(2*sinThetaVariance/(4-M_PI));
      res.rSigmaML = sqrt(sumSin2Theta/(2*n));

      res.sigmax = sqrt(res.vpos[0][0]);
      res.sigmay = sqrt(res.vpos[1][1]);
      res.sigmaz = sqrt(res.vpos[2][2]);

      switch(particle->srcType) {
      case SourceFront: case SourceBack: case SourceSignal:
        res.correlationCoefficient = res.vpos[0][1]/(res.sigmax*res.sigmay);
        break;
      case SourceSouthWest: case SourceNorthEast:
        res.correlationCoefficient = res.vpos[2][1]/(res.sigmaz*res.sigmay);
        break;
      case SourceFloor: case SourceCeiling:
        res.correlationCoefficient = res.vpos[0][2]/(res.sigmax*res.sigmaz);
        break;
      default: assert(false);
      }

      return res;

    } // computeParticleRandomization()

    //================================================================
    void ExtMonFNALRoomGenerator::fillRandomizationHistograms(SourceType st, ParticleType pt, const ParticleRandomization& pr) {
      std::pair<SourceType, ParticleType> stpt(st,pt);
      HistMap::iterator ih = histRandomization_.find(stpt);
      if(ih == histRandomization_.end()) {
        // Book
        ih = histRandomization_.insert(std::make_pair(stpt, HistRandomization(st,pt))).first;
        assert(ih != histRandomization_.end());
      }

      ih->second.fill(pr);
    }

    //================================================================
    void ExtMonFNALRoomGenerator::produce(art::Event& event) {

      AGDEBUG("event "<<event.id());

      std::auto_ptr<GenParticleCollection> output(new GenParticleCollection);
      std::auto_ptr<MARSInfoCollection> info(new MARSInfoCollection());
      std::auto_ptr<GenParticleMARSAssns> assns(new GenParticleMARSAssns());

      const art::ProductID particlesPID = getProductID<GenParticleCollection>(event);
      const art::EDProductGetter *particlesGetter = event.productGetter(particlesPID);

      const art::ProductID marsPID = getProductID<MARSInfoCollection>(event);
      const art::EDProductGetter *marsGetter = event.productGetter(marsPID);


      // Process requested number of protons
      for(unsigned count = 0; count < numReadProtonsPerEvent_; ++count) {

        const InputParticle& pp = particles_[nextParticle_];
        const bool acceptProton(randFlat_.fire() <= fracAcceptedProtons_);

        // process (or skip) all particles from the current proton
        while(pp.info.protonNumber() == particles_[nextParticle_].info.protonNumber()) {
          if(acceptProton) {
            output->push_back(createOutputParticle(nextParticle_));
            info->push_back(particles_[nextParticle_].info);
            assns->addSingle(art::Ptr<GenParticle>(particlesPID, output->size()-1, particlesGetter),
                             art::Ptr<MARSInfo>(marsPID, info->size()-1, marsGetter));
          }
          ++nextParticle_ %= particles_.size();
        }
      }

      event.put(output);
      event.put(info);
      event.put(assns);
    }

    //================================================================
    GenParticle ExtMonFNALRoomGenerator::createOutputParticle(InputParticles::size_type ip) {
      using CLHEP::Hep3Vector;

      const InputParticle& pp = particles_[ip];
      const ParticleRandomization& rr = pr_[ip];

      const double sigmax = (pp.srcType == SourceSouthWest)||(pp.srcType == SourceNorthEast) ?
        0 : rr.sigmax;

      const double sigmay = (pp.srcType == SourceFloor)||(pp.srcType == SourceCeiling) ?
        0 : rr.sigmay;

      const double sigmaz = (pp.srcType == SourceFront)||(pp.srcType == SourceBack)||(pp.srcType == SourceSignal) ?
        0 : rr.sigmaz;

      Hep3Vector posDump(pp.posDump);
      do{
        posDump.setX(pp.posDump.x() + sigmax * randGauss_.fire());
        posDump.setY(pp.posDump.y() + sigmay * randGauss_.fire());
        posDump.setZ(pp.posDump.z() + sigmaz * randGauss_.fire());
      } while(!inRange(pp.srcType, posDump));

      const Hep3Vector posMu2e(dump_->beamDumpToMu2e_position(posDump));

      // Draw dtheta from the Rayleigh distribution
      const double dtheta = rr.rSigmaML * sqrt(-2*log(randFlat_.fire()));
      const double dphi = 2*M_PI*randFlat_.fire();

      const Hep3Vector& orig(pp.momMu2e);

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

      const CLHEP::HepLorentzVector momMu2e(randomized3mom, pp.momMu2e.e());

      return GenParticle(pp.pdgId, GenId::MARS, posMu2e, momMu2e, pp.time);
    }

    //================================================================
    bool ExtMonFNALRoomGenerator::inRange(SourceType st, const CLHEP::Hep3Vector& posDump) {

      switch(st) {

      default: assert(false);

      case SourceFront: case SourceBack:
        return
          (srcGeom_.xSW <= posDump.x()) && (posDump.x() <= srcGeom_.xNE) &&
          (srcGeom_.yFloor <= posDump.y()) && (posDump.y() <= srcGeom_.yCeiling);

      case SourceSouthWest: case SourceNorthEast:
        return
          (srcGeom_.zBack <= posDump.z()) && (posDump.z() <= srcGeom_.zFront) &&
          (srcGeom_.yFloor <= posDump.y()) && (posDump.y() <= srcGeom_.yCeiling);

      case SourceFloor: case SourceCeiling:
        return
          (srcGeom_.xSW <= posDump.x()) && (posDump.x() <= srcGeom_.xNE) &&
          (srcGeom_.zBack <= posDump.z()) && (posDump.z() <= srcGeom_.zFront);

      case SourceSignal: {
        CLHEP::Hep3Vector posExtMon =
          extmon_->mu2eToExtMon_position(dump_->beamDumpToMu2e_position(posDump));

        return
          (std::abs(posExtMon.x()) <= signalHalfdx_) &&
          (std::abs(posExtMon.y()) <= signalHalfdy_);
      }
      } // switch(st)
    } // inRange()

    //================================================================
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::ExtMonFNALRoomGenerator);

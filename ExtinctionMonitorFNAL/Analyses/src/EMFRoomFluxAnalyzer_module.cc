// Precompute particle randomization from MARS inputs for g4s1 ExtMonFNALRoom jobs.
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

#include <sys/times.h>
#include <stdio.h>
#include <unistd.h>

#include "cetlib_except/exception.h"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "GeneralUtilities/inc/KNearestNeighbors.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/MARSInfo.hh"
#include "MCDataProducts/inc/MARSInfoCollection.hh"
#include "MCDataProducts/inc/GenParticleMARSAssns.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "ExtinctionMonitorFNAL/Utilities/inc/EMFBoxIO.hh"
#include "ExtinctionMonitorFNAL/Utilities/inc/getCharge.hh"
#include "ExtinctionMonitorFNAL/Utilities/inc/EMFRandomizationParticleDefs.hh"
#include "ExtinctionMonitorFNAL/Utilities/inc/EMFRandomizationSourceDefs.hh"

#include "TTree.h"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    using namespace Randomization;

    namespace {

      //================================================================
      struct InputParticle {
        CLHEP::Hep3Vector posDump;
        CLHEP::HepLorentzVector momMu2e;
        double time;
        PDGCode::type pdgId;

        SourceType srcType;

        MARSInfo info;

        InputParticle(const ProtonBeamDump& dump,
                      SourceType src,
                      const GenParticle& gp,
                      const MARSInfo& inf)
          : posDump(dump.mu2eToBeamDump_position(gp.position()))
          , momMu2e(gp.momentum())
          , time(gp.time())
          , pdgId(gp.pdgId())
          , srcType(src)
          , info(inf)
        {}

      };

      typedef std::vector<InputParticle> InputParticles;

      std::ostream& operator<<(std::ostream& os, const InputParticle& part) {
        return os<<"Particle(posDump="<<part.posDump<<", momMu2e="<<part.momMu2e<<")";
      }

      //================================================================
      using IO::ParticleRandomization;
      typedef std::vector<ParticleRandomization> ParticleRandomizations;

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

    class EMFRoomFluxAnalyzer : public art::EDAnalyzer {
      int verbosityLevel_;

      std::string generatorModuleLabel_;
      std::string generatorInstanceName_;

      std::string geomModuleLabel_;
      std::string geomInstanceName_;

      const ProtonBeamDump *dump_;
      const ExtMon *extmon_;
      SourcePlaneGeom srcGeom_;

      unsigned numNeighbors_;
      unsigned minSourceGroupStatistics_;

      InputParticles particles_;
      ParticleRandomizations pr_; // Index synced with particles_

      InputParticles::size_type indexInParticles(const InputParticle *p) {
        return p - &particles_[0];
      }

      // inputs arranged by [SourceType][ParticleType][particleIndex]
      typedef std::vector<const InputParticle*> RandomizationGroup;
      typedef std::vector<std::vector<RandomizationGroup> > GroupedInputs;
      GroupedInputs grouped_;

      void initGroups();

      void printSrcGroups();
      void mergeLowStatisticSrcGroups();
      void computeParticleRandomizations();
      void writeParticleNtuple();

      ParticleRandomization computeParticleRandomization(const InputParticle *p,
                                                         const KNearestNeighbors<const InputParticle*>::Points& neighbors);

    public:
      explicit EMFRoomFluxAnalyzer(const fhicl::ParameterSet& pset);
      virtual void analyze(const art::Event& event);
      virtual void beginRun(const art::Run& run);
      virtual void endJob();
    };

    //================================================================
    EMFRoomFluxAnalyzer::EMFRoomFluxAnalyzer(const fhicl::ParameterSet& pset)
      : art::EDAnalyzer(pset)

      , verbosityLevel_(pset.get<int>("verbosityLevel"))

      , generatorModuleLabel_(pset.get<std::string>("generatorModuleLabel"))
      , generatorInstanceName_(pset.get<std::string>("generatorInstanceName", ""))

      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))

      , dump_()
      , extmon_()

      , srcGeom_(pset.get<fhicl::ParameterSet>("srcGeom"))

      , numNeighbors_(pset.get<unsigned>("numNeighbors"))
      , minSourceGroupStatistics_(pset.get<unsigned>("minSourceGroupStatistics"))

      , grouped_(NUM_SOURCES, std::vector<RandomizationGroup>(NUM_PARTICLE_TYPES))
    {
      if(minSourceGroupStatistics_ < 1 + numNeighbors_) {
        std::cout<<"WARNING: EMFRoomFluxAnalyzer: "
                 <<"minSourceGroupStatistics ("<<minSourceGroupStatistics_
                 <<") < 1 + numNeighbors (1 + "<<numNeighbors_<<")\n";
      }

      if(numNeighbors_ < 1) {
        throw cet::exception("BADCONFIG")
          <<"Error: numNeighbors="<<numNeighbors_<<", should be > 0";
      }
    }

    //================================================================
    void EMFRoomFluxAnalyzer::beginRun(const art::Run& run) {

//old_geom:      const ExtMonFNALBuilding *building(0);

      if(!geomModuleLabel_.empty()) {
        art::Handle<ProtonBeamDump> dump;
        run.getByLabel(geomModuleLabel_, geomInstanceName_, dump);
        dump_ = &*dump;

        art::Handle<ExtMon> extmon;
        run.getByLabel(geomModuleLabel_, geomInstanceName_, extmon);
        extmon_ = &*extmon;

//old_geom:        art::Handle<ExtMonFNALBuilding> bh;
//old_geom:        run.getByLabel(geomModuleLabel_, geomInstanceName_, bh);
//old_geom:        building = &*bh;

      }
      else {
        GeomHandle<ProtonBeamDump> dump;
        dump_ = &*dump;

        GeomHandle<ExtMon> extmon;
        extmon_ = &*extmon;

//old_geom:       GeomHandle<ExtMonFNALBuilding> bh;
//old_geom:       building = &*bh;
      }


//old_geom:      const double yFloor = dump_->mu2eToBeamDump_position(CLHEP::Hep3Vector(0,building->roomInsideYmin(),0)).y();
//old_geom:      const double yCeiling = dump_->mu2eToBeamDump_position(CLHEP::Hep3Vector(0,building->roomInsideYmax(),0)).y();
//old_geom:
//old_geom:      const double zFront = dump_->mu2eToBeamDump_position(building->coll2ShieldingCenterInMu2e()).z()
//old_geom:        - building->coll2ShieldingHalfSize()[2]
//old_geom:        ;
//old_geom:
//old_geom:      const double xSW =  building->coll2ShieldingDumpXmin();
//old_geom:      const double xNE =  building->coll2ShieldingDumpXmax();
//old_geom:
//old_geom:      if(verbosityLevel_ > 0) {
//old_geom:        std::cout<<"EMFRoomFluxAnalyzer INFO: src zFront offset into room   = "<<(zFront - srcGeom_.zFront)<<std::endl;
//old_geom:        std::cout<<"EMFRoomFluxAnalyzer INFO: src yFloor offset into room   = "<<(srcGeom_.yFloor - yFloor)<<std::endl;
//old_geom:        std::cout<<"EMFRoomFluxAnalyzer INFO: src yCeiling offset into room = "<<(yCeiling - srcGeom_.yCeiling)<<std::endl;
//old_geom:        std::cout<<"EMFRoomFluxAnalyzer INFO: src xSW offset into room      = "<<(srcGeom_.xSW - xSW)<<std::endl;
//old_geom:        std::cout<<"EMFRoomFluxAnalyzer INFO: src xNE offset into room      = "<<(xNE - srcGeom_.xNE)<<std::endl;
//old_geom:      }
    }

    //================================================================
    void EMFRoomFluxAnalyzer::printSrcGroups() {
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
    void EMFRoomFluxAnalyzer::initGroups() {
      for(unsigned i=0; i<particles_.size(); ++i) {
        const InputParticle& part = particles_[i];
        SourceType st = part.srcType;
        ParticleType pt = classifyParticleType(part.pdgId);
        grouped_[st][pt].push_back(&part);
      }
    }

    //================================================================
    void EMFRoomFluxAnalyzer::mergeLowStatisticSrcGroups() {
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

            if(verbosityLevel_ > 2) {
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

    void EMFRoomFluxAnalyzer::computeParticleRandomizations() {
      static MetricXY mxy; // gcc 4.6.1 does not like "static const" here
      static MetricYZ myz; // gcc 4.6.1 does not like "static const" here
      static MetricZX mzx; // gcc 4.6.1 does not like "static const" here
      static const Metric *dist[] = { &mxy, &myz, &mzx };

      pr_.resize(particles_.size());

      for(unsigned st = 0; st < NUM_SOURCES; ++st) {

        for(unsigned pt = 0; pt < NUM_PARTICLE_TYPES; ++pt) {

          const RandomizationGroup& group = grouped_[st][pt];

          if(!group.empty()) {


            if(verbosityLevel_ > 1) {
              std::cout<<"EMFRoomFluxAnalyzer: starting K nearest neighbors search for group size "
                       <<group.size()<<std::endl;
            }

            struct tms st_cpu;
            const clock_t st_time = times(&st_cpu);

            KNearestNeighbors<const InputParticle*>
              neighbors(numNeighbors_, group, *dist[st % 3]);

            struct tms en_cpu;
            const clock_t en_time = times(&en_cpu);

            if(verbosityLevel_ > 1) {
              static const double ticks_per_second =  sysconf(_SC_CLK_TCK);
              std::cout<<"EMFRoomFluxAnalyzer: KNN search for group size "
                       <<group.size()<<" took "
                       <<(en_time - st_time)/ticks_per_second<<" (real), "
                       <<(en_cpu.tms_utime - st_cpu.tms_utime)/ticks_per_second<<" (user) "
                       <<(en_cpu.tms_stime - st_cpu.tms_stime)/ticks_per_second<<" (system)"
                       <<" seconds"<<std::endl;
            }

            if(verbosityLevel_ > 2) {
              std::cout<<"EMFRoomFluxAnalyzer: Neighbors dump begin"<<std::endl;
              std::cout<<neighbors<<std::endl;
              std::cout<<"EMFRoomFluxAnalyzer: Neighbors dump end"<<std::endl;
            }

            for(unsigned i=0; i<group.size(); ++i) {
              InputParticles::size_type ip = indexInParticles(group[i]);
              pr_[ip] = computeParticleRandomization(group[i], neighbors[i]);
            }

          } // !group.empty()
        } // for(particle types pt)
      } // for(sources st)
    } // computeParticleRandomizations()

    //================================================================
    ParticleRandomization
    EMFRoomFluxAnalyzer::computeParticleRandomization(const InputParticle *particle,
                                                      const KNearestNeighbors<const InputParticle*>::Points& neighbors)
    {
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
      }

      const double n(neighbors.size());

      ParticleRandomization res;
      // Note  that here we divide by n rather than (n-1) because of the hypothesis mean==0
      // (we do symmetric smearing in coordinates, which does not bias the mean).
      const CLHEP::HepMatrix variancePos(tmp/n);
      res.sigmax = sqrt(variancePos[0][0]);
      res.sigmay = sqrt(variancePos[1][1]);
      res.sigmaz = sqrt(variancePos[2][2]);

      res.rSigmaML = sqrt(sumSin2Theta/(2*n));

      switch(particle->srcType) {
      case SourceFront: case SourceBack: case SourceSignal:
        res.correlationCoefficient = variancePos[0][1]/(res.sigmax*res.sigmay);
        break;
      case SourceSouthWest: case SourceNorthEast:
        res.correlationCoefficient = variancePos[2][1]/(res.sigmaz*res.sigmay);
        break;
      case SourceFloor: case SourceCeiling:
        res.correlationCoefficient = variancePos[0][2]/(res.sigmax*res.sigmaz);
        break;
      default: assert(false);
      }

      return res;

    } // computeParticleRandomization()

    //================================================================
    void EMFRoomFluxAnalyzer::analyze(const art::Event& event) {

      art::Handle<GenParticleCollection> genh;
      event.getByLabel(generatorModuleLabel_, generatorInstanceName_, genh);
      art::FindOne<MARSInfo>
        mFinder(genh, event, art::InputTag(generatorModuleLabel_, generatorInstanceName_));

      const GenParticleCollection& genparts(*genh);
      for(unsigned i=0; i < genparts.size(); ++i) {
        MARSInfo info = mFinder.at(i).ref();
        SourceType srcType = classifySource(genparts[i].position(), srcGeom_, *dump_, *extmon_);
        particles_.push_back(InputParticle(*dump_, srcType, genparts[i], info));
      }

    } // analyze()

    //================================================================
    void EMFRoomFluxAnalyzer::endJob() {

      initGroups();

      if(verbosityLevel_ > 0) {
        std::cout<<"Summary of inputs before regrouping:"<<std::endl;
        printSrcGroups();
      }

      mergeLowStatisticSrcGroups();

      if(verbosityLevel_ > 0) {
        std::cout<<"Summary of inputs after regrouping:"<<std::endl;
        printSrcGroups();
      }

      computeParticleRandomizations();

      writeParticleNtuple();
    }


    //================================================================
    void EMFRoomFluxAnalyzer::writeParticleNtuple() {
      art::ServiceHandle<art::TFileService> tfs;

      IO::EMFRoomHit particle;
      IO::MARSInfo minfo;
      ParticleRandomization pr;

      TTree *nt = tfs->make<TTree>("roomhits", "MARS input particles and randomization parameters");
      nt->Branch("particle", &particle, particle.branchDescription());
      nt->Branch("minfo", &minfo, minfo.branchDescription());
      nt->Branch("randomization", &pr, pr.branchDescription());

      for(unsigned i=0; i<particles_.size(); ++i) {

        particle.dumpx = particles_[i].posDump.x();
        particle.dumpy = particles_[i].posDump.y();
        particle.dumpz = particles_[i].posDump.z();

        particle.mu2epx = particles_[i].momMu2e.x();
        particle.mu2epy = particles_[i].momMu2e.y();
        particle.mu2epz = particles_[i].momMu2e.z();
        particle.time = particles_[i].time;
        particle.pdgId = particles_[i].pdgId;

        particle.srcType = particles_[i].srcType;

        minfo.info = particles_[i].info;

        pr = pr_[i];

        nt->Fill();
      }
    }

    //================================================================
  } // namespace ExtMonFNAL
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::EMFRoomFluxAnalyzer);

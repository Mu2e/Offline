// Read in a set of particles hitting ExtMonFNAL VD box (from g4s1 room jobs)
// compute randomization parameters, and write out as an ntuple.
//
// $Id: EMFBoxFluxAnalyzer_module.cc,v 1.13 2013/10/21 20:34:14 gandr Exp $
// $Author: gandr $
// $Date: 2013/10/21 20:34:14 $
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
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/MARSInfo.hh"
#include "MCDataProducts/inc/MARSInfoCollection.hh"
#include "MCDataProducts/inc/GenParticleMARSAssns.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "GlobalConstantsService/inc/MassCache.hh"

#include "ExtinctionMonitorFNAL/Utilities/inc/EMFBoxIO.hh"
#include "ExtinctionMonitorFNAL/Utilities/inc/getCharge.hh"
#include "ExtinctionMonitorFNAL/Utilities/inc/EMFRandomizationParticleDefs.hh"


#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    namespace {

      using namespace Randomization;

      unsigned const NUM_SOURCES = VirtualDetectorId::EMFBoxTop - VirtualDetectorId::EMFBoxFront + 1;
      unsigned sourceNumber(VirtualDetectorId vid) { return vid.id() - VirtualDetectorId::EMFBoxFront; }
      VirtualDetectorId vdFromSrcNum(unsigned st)  { return VirtualDetectorId(st + VirtualDetectorId::EMFBoxFront); }

      // FIXME: this number is known to ConditionsService, which is not available at beginJob()
      const double deBuncherPeriod = 1694.;

      //================================================================
      using IO::ParticleRandomization;

      //================================================================
      struct InputParticle {
        CLHEP::Hep3Vector posExtMon;
        CLHEP::Hep3Vector momMu2e;
        double time;
        PDGCode::type pdgId;

        VirtualDetectorId vd;
        ParticleType      origParticleType;

        ParticleRandomization pr;
        MARSInfo minfo;
        IO::G4JobInfo g4s1info;

        InputParticle(const ExtMon& extmon,
                      VirtualDetectorId vid,
                      const StepPointMC& hit,
                      const MARSInfo& inf,
                      const IO::G4JobInfo& g4s1i)
          : posExtMon(extmon.mu2eToExtMon_position(hit.position()))
          , momMu2e(hit.momentum())
          , time(hit.time())
          , pdgId(hit.simParticle()->pdgId())
          , vd(vid)
          , origParticleType(classifyParticleType(pdgId))

          , pr()

          , minfo(inf)
          , g4s1info(g4s1i)
        {}

      };

      typedef std::vector<InputParticle> InputParticles;

      std::ostream& operator<<(std::ostream& os, const InputParticle& part) {
        return os<<"Particle(posExtMon="<<part.posExtMon<<", momMu2e="<<part.momMu2e<<")";
      }

      // a quick hack to add transient-only TOF
      struct TransientParticle : public InputParticle {
        double beta;
        double dtbeta;
        TransientParticle(const InputParticle& p, double b, double cutMinTime)
          : InputParticle(p)
          , beta(b)
          , dtbeta(beta * (p.time - cutMinTime))
        {}
      };

      //----------------
      // Here we operate on an ensemble of g4s1 events
      // so need to look at identifiers from all the steps,
      // from MARS protonNumber to g4s1 run number.
      struct ProtonPathSort {
        bool operator()(const InputParticle& a, const InputParticle& b) {
          CmpProtonIdAndSimPath cm;
          IO::CmpG4JobInfo cg;
          return cm(a.minfo,b.minfo) || (!cm(b.minfo, a.minfo) &&
                                         cg(a.g4s1info, b.g4s1info));
        }
      };

      //================================================================
      class HistRandomization {
        TH1 *rsigmatheta;
        TH1 *sigmax;
        TH1 *sigmay;
        TH1 *sigmaz;
        TH1 *correlationCoefficient;

      public:
        HistRandomization(VirtualDetectorId st, ParticleType pt);
        void fill(const ParticleRandomization& pr);
      };

      HistRandomization::HistRandomization(VirtualDetectorId st, ParticleType pt)
        : rsigmatheta(), sigmax(), sigmay(), sigmaz(), correlationCoefficient()
      {
        art::ServiceHandle<art::TFileService> tfs;
        art::TFileDirectory tftop = tfs->mkdir("randomization");

        std::ostringstream os;
        os<<st.name()<<"_pt"<<pt;

        art::TFileDirectory tfdir = tftop.mkdir(os.str());

        rsigmatheta = tfdir.make<TH1D>("rsigmatheta", "Rayleigh sigma theta", 500, 0., 1.);
        sigmax = tfdir.make<TH1D>("sigmax", "sigma x", 1000, 0., 500.);
        sigmay = tfdir.make<TH1D>("sigmay", "sigma y", 1000, 0., 500.);
        sigmaz = tfdir.make<TH1D>("sigmaz", "sigma z", 1000, 0., 500.);
        correlationCoefficient = tfdir.make<TH1D>("correlation", "correlation", 101, -1., 1.);

      }

      void HistRandomization::fill(const ParticleRandomization& pr) {
        rsigmatheta->Fill(pr.rSigmaML);
        sigmax->Fill(pr.sigmax);
        sigmay->Fill(pr.sigmay);
        sigmaz->Fill(pr.sigmaz);
        correlationCoefficient->Fill(pr.correlationCoefficient);
      }

      //================================================================
      class HistInputs {
        TH1* histInputTimesEarly_;
        TH1* histInputTimesBunch_;
        TH1* histInputBeta_;
        TH1* histInputLogBeta_;
        TH1* histInputBetadt_;
        void init(const std::string& dirname);
      public:
        HistInputs(const std::string& dirname);
        HistInputs(VirtualDetectorId st, ParticleType pt);
        void fill(const TransientParticle& particle);
      };

      HistInputs::HistInputs(const std::string& dirname) {
        init(dirname);
      }

      HistInputs::HistInputs(VirtualDetectorId st, ParticleType pt)
        : histInputTimesEarly_(), histInputTimesBunch_()
        , histInputBeta_(), histInputLogBeta_()
        , histInputBetadt_()
      {
        art::ServiceHandle<art::TFileService> tfs;
        art::TFileDirectory tftop = tfs->mkdir("inputs");

        std::ostringstream os;
        os<<st.name()<<"_pt"<<pt;
        init(os.str());
      }

      void HistInputs::init(const std::string& subdir) {
        art::ServiceHandle<art::TFileService> tfs;
        art::TFileDirectory tftop = tfs->mkdir("inputs");
        art::TFileDirectory tfdir = tftop.mkdir(subdir);

        histInputTimesEarly_ = tfdir.make<TH1D>("hitTimesEarly", "Input hit time", 500, 0., 500.);
        histInputTimesBunch_ = tfdir.make<TH1D>("hitTimesBunch", "Input hit time MOD deBuncherPeriod", 170, 0., deBuncherPeriod);

        histInputBeta_ = tfdir.make<TH1D>("hitBeta", "Particle v/c", 500, 0., 1.);
        histInputLogBeta_ = tfdir.make<TH1D>("hitLogBeta", "Particle log10(v/c)", 600, -12., 0.);
        histInputBetadt_ = tfdir.make<TH1D>("hitBetadt", "beta*(t-tcut)", 1000, -500., 500.);
      }

      void HistInputs::fill(const TransientParticle& particle) {
        histInputTimesEarly_->Fill(particle.time);
        histInputTimesBunch_->Fill(fmod(particle.time, deBuncherPeriod));

        histInputBeta_->Fill(particle.beta);
        histInputLogBeta_->Fill(log10(particle.beta));

        histInputBetadt_->Fill(particle.dtbeta);
      }

      //----------------------------------------------------------------
      class HistInputsCollection {
        HistInputs total_;
        std::vector<std::vector<HistInputs> > hists_;
      public:
        HistInputsCollection();
        void fill(const TransientParticle& particle) {
          total_.fill(particle);
          hists_.at(sourceNumber(particle.vd)).at(particle.origParticleType).fill(particle);
        }
      };

      HistInputsCollection::HistInputsCollection()
        : total_("total")
        , hists_(std::vector<std::vector<HistInputs> >(NUM_SOURCES))
      {
        art::ServiceHandle<art::TFileService> tfs;
        art::TFileDirectory tftop = tfs->mkdir("inputs");

        for(unsigned st=0; st<NUM_SOURCES; ++st) {
          for(unsigned pt=0; pt<NUM_PARTICLE_TYPES; ++pt) {
            hists_[st].push_back(HistInputs(vdFromSrcNum(st), ParticleType(pt)));
          }
        }
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

    class EMFBoxFluxAnalyzer : public art::EDAnalyzer {
      int verbosityLevel_;

      std::string hitsModuleLabel_;
      std::string hitsInstanceName_;

      std::string marsInfoModuleLabel_;
      std::string marsInfoInstanceName_;

      std::string geomModuleLabel_;
      std::string geomInstanceName_;

      const ExtMon *extmon_;
      MassCache mc_;

      double   cutMinTime_;
      double   cutLengthScale_;

      bool noKNN_;
      unsigned numNeighbors_;
      unsigned minSourceGroupStatistics_;

      InputParticles particles_;
      InputParticles::size_type indexInParticles(const InputParticle *p) {
        return p - &particles_[0];
      }

      double particleBeta(const InputParticle& particle);
      bool isAccepted(const TransientParticle& particle);

      // inputs arranged by [VirtualDetectorId][ParticleType][particleIndex]
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

      typedef std::map<std::pair<VirtualDetectorId,ParticleType>, HistRandomization> HistMapRandomization;
      HistMapRandomization histRandomization_;
      void fillRandomizationHistograms(VirtualDetectorId st, ParticleType pt, const ParticleRandomization& pr);

      HistInputsCollection histInputs_;

    public:
      explicit EMFBoxFluxAnalyzer(const fhicl::ParameterSet& pset);
      virtual void analyze(const art::Event& event);
      virtual void beginRun(const art::Run& run);
      virtual void endJob();
    };

    //================================================================
    EMFBoxFluxAnalyzer::EMFBoxFluxAnalyzer(const fhicl::ParameterSet& pset)
      : art::EDAnalyzer(pset)
      , verbosityLevel_(pset.get<int>("verbosityLevel"))
      , hitsModuleLabel_(pset.get<std::string>("hitsModuleLabel"))
      , hitsInstanceName_(pset.get<std::string>("hitsInstanceName", ""))

      , marsInfoModuleLabel_(pset.get<std::string>("marsInfoModuleLabel"))
      , marsInfoInstanceName_(pset.get<std::string>("marsInfoInstanceName", ""))

      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))

      , extmon_()

      , cutMinTime_(pset.get<double>("cutMinTime"))
      , cutLengthScale_(pset.get<double>("cutLengthScale"))

      , noKNN_(pset.get<bool>("noKNN", false))
      , numNeighbors_(noKNN_ ? pset.get<unsigned>("numNeighbors", 0) : pset.get<unsigned>("numNeighbors"))
      , minSourceGroupStatistics_(noKNN_ ? pset.get<unsigned>("minSourceGroupStatistics", 0) : pset.get<unsigned>("minSourceGroupStatistics"))

      , grouped_(NUM_SOURCES, std::vector<RandomizationGroup>(NUM_PARTICLE_TYPES))
    {
      if(!noKNN_) {
        if(minSourceGroupStatistics_ < numNeighbors_) {
          throw cet::exception("BADCONFIG")
            <<"Error: minSourceGroupStatistics ("<<minSourceGroupStatistics_<<") < numNeighbors ("<<numNeighbors_<<")";
        }

        if(numNeighbors_ < 2) {
          throw cet::exception("BADCONFIG")
            <<"Error: numNeighbors="<<numNeighbors_<<", should be > 1";
        }

        std::cout<<"cutLengthScale/c = "<<(cutLengthScale_ / CLHEP::c_light)<<" ns"
                 <<", cutMinTime = "<<cutMinTime_<<std::endl;
      }
    }

    //================================================================
    void EMFBoxFluxAnalyzer::beginRun(const art::Run& run) {
      if(!geomModuleLabel_.empty()) {
        art::Handle<ExtMon> extmon;
        run.getByLabel(geomModuleLabel_, geomInstanceName_, extmon);
        extmon_ = &*extmon;
      }
      else {
        GeomHandle<ExtMon> extmon;
        extmon_ = &*extmon;
      }
    }

    //================================================================
    void EMFBoxFluxAnalyzer::printSrcGroups() {
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
    void EMFBoxFluxAnalyzer::initGroups() {
      for(unsigned i=0; i<particles_.size(); ++i) {
        const InputParticle& part = particles_[i];
        VirtualDetectorId st = part.vd;
        ParticleType pt = part.origParticleType;
        grouped_[sourceNumber(st)][pt].push_back(&part);
      }
    }

    //================================================================
    void EMFBoxFluxAnalyzer::mergeLowStatisticSrcGroups() {
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

            if(!noKNN_) {
              throw cet::exception("BADINPUTS")
                <<"Error: failed to regroup particles on source plane = "<<st
                <<" to satisfy minSourceGroupStatistics="<<minSourceGroupStatistics_
                <<": achieved stat="<<grouped_[st][pt].size()
                <<"\n";
            }
            else {
              std::cout<<"EMFBoxFluxAnalyzer: WARNING: "
                       <<"failed to regroup particles on source plane = "<<st
                       <<" to satisfy minSourceGroupStatistics="<<minSourceGroupStatistics_
                       <<": achieved stat="<<grouped_[st][pt].size()
                       <<std::endl;
            }
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
        const double dx = a->posExtMon.x() - b->posExtMon.x();
        const double dy = a->posExtMon.y() - b->posExtMon.y();
        return sqrt(dx*dx + dy*dy);
      }
    };

    class MetricYZ : virtual public Metric {
    public:
      virtual double operator()(const InputParticle *a, const InputParticle *b) const {
        const double dy = a->posExtMon.y() - b->posExtMon.y();
        const double dz = a->posExtMon.z() - b->posExtMon.z();
        return sqrt(dy*dy + dz*dz);
      }
    };

    class MetricZX : virtual public Metric {
    public:
      virtual double operator()(const InputParticle *a, const InputParticle *b) const {
        const double dx = a->posExtMon.x() - b->posExtMon.x();
        const double dz = a->posExtMon.z() - b->posExtMon.z();
        return sqrt(dx*dx + dz*dz);
      }
    };

    void EMFBoxFluxAnalyzer::computeParticleRandomizations() {
      static MetricXY mxy; // gcc 4.6.1 does not like "static const" here
      static MetricYZ myz; // gcc 4.6.1 does not like "static const" here
      static MetricZX mzx; // gcc 4.6.1 does not like "static const" here
      static const Metric *dist[] = { &mxy, &myz, &mzx };

      for(unsigned st = 0; st < NUM_SOURCES; ++st) {

        for(unsigned pt = 0; pt < NUM_PARTICLE_TYPES; ++pt) {

          const RandomizationGroup& group = grouped_[st][pt];

          if(!group.empty()) {


            if(verbosityLevel_ > 1) {
              std::cout<<"EMFBoxFluxAnalyzer: starting K nearest neighbors search for group size "
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
             std::cout<<"EMFBoxFluxAnalyzer: KNN search for group size "
                      <<group.size()<<" took "
                      <<(en_time - st_time)/ticks_per_second<<" (real), "
                      <<(en_cpu.tms_utime - st_cpu.tms_utime)/ticks_per_second<<" (user) "
                      <<(en_cpu.tms_stime - st_cpu.tms_stime)/ticks_per_second<<" (system)"
                      <<" seconds"<<std::endl;
            }

            if(verbosityLevel_ > 2) {
              std::cout<<"EMFBoxFluxAnalyzer: Neighbors dump begin"<<std::endl;
              std::cout<<neighbors<<std::endl;
              std::cout<<"EMFBoxFluxAnalyzer: Neighbors dump end"<<std::endl;
            }

            for(unsigned i=0; i<group.size(); ++i) {
              InputParticles::size_type ip = indexInParticles(group[i]);
              particles_[ip].pr = computeParticleRandomization(group[i], neighbors[i]);
              fillRandomizationHistograms(vdFromSrcNum(st), ParticleType(pt), particles_[ip].pr);
            }

          } // !group.empty()
        } // for(particle types pt)
      } // for(sources st)
    } // computeParticleRandomizations()

    //================================================================
    ParticleRandomization
    EMFBoxFluxAnalyzer::computeParticleRandomization(const InputParticle *particle,
                                                     const KNearestNeighbors<const InputParticle*>::Points& neighbors)
    {
      double sumSin2Theta(0);
      CLHEP::HepMatrix tmp(3,3);
      for(unsigned i=0; i<neighbors.size(); ++i) {
        CLHEP::HepVector dr(3);
        dr[0] = particle->posExtMon.x() - neighbors[i].point->posExtMon.x();
        dr[1] = particle->posExtMon.y() - neighbors[i].point->posExtMon.y();
        dr[2] = particle->posExtMon.z() - neighbors[i].point->posExtMon.z();
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

      switch(particle->vd.id()) {
      case VirtualDetectorId::EMFBoxFront:
      case VirtualDetectorId::EMFBoxBack:
        res.correlationCoefficient = variancePos[0][1]/(res.sigmax*res.sigmay);
        break;
      case VirtualDetectorId::EMFBoxSW:
      case VirtualDetectorId::EMFBoxNE:
        res.correlationCoefficient = variancePos[2][1]/(res.sigmaz*res.sigmay);
        break;
      case VirtualDetectorId::EMFBoxBottom:
      case VirtualDetectorId::EMFBoxTop:
        res.correlationCoefficient = variancePos[0][2]/(res.sigmax*res.sigmaz);
        break;
      default: assert(false);
      }

      return res;

    } // computeParticleRandomization()

    //================================================================
    void EMFBoxFluxAnalyzer::fillRandomizationHistograms(VirtualDetectorId st, ParticleType pt, const ParticleRandomization& pr) {
      std::pair<VirtualDetectorId, ParticleType> stpt(st,pt);
      HistMapRandomization::iterator ih = histRandomization_.find(stpt);
      if(ih == histRandomization_.end()) {
        // Book
        ih = histRandomization_.insert(std::make_pair(stpt, HistRandomization(st,pt))).first;
        assert(ih != histRandomization_.end());
      }

      ih->second.fill(pr);
    }

    //================================================================
    double EMFBoxFluxAnalyzer::particleBeta(const InputParticle& particle) {
      double beta = 1;
      const double mass = mc_.mass(particle.pdgId);
      if(mass > 0.) {
        const double pm2 = particle.momMu2e.mag2()/std::pow(mass, 2);
        beta = sqrt(pm2/(1.+pm2));
      }
      return beta;
    }

    //================================================================
    bool EMFBoxFluxAnalyzer::isAccepted(const TransientParticle& particle) {

      //  beta*(t - tcut) + L/c > 0:
      if(particle.dtbeta + cutLengthScale_/CLHEP::c_light > 0) {

        // Don't need to record particles exiting the box.
        // Check the momentum direction

        const CLHEP::Hep3Vector momExtMon = extmon_->mu2eToExtMon_momentum(particle.momMu2e);

        switch(particle.vd.id()) {

        case VirtualDetectorId::EMFBoxFront: return momExtMon.z() < 0;
        case VirtualDetectorId::EMFBoxBack:  return momExtMon.z() > 0;

        case VirtualDetectorId::EMFBoxSW: return momExtMon.x() > 0;
        case VirtualDetectorId::EMFBoxNE: return momExtMon.x() < 0;

        case VirtualDetectorId::EMFBoxBottom: return momExtMon.y() > 0;
        case VirtualDetectorId::EMFBoxTop:    return momExtMon.y() < 0;

        default: assert(false);
        } // switch()

      }// cut time

      return false;
    }

    //================================================================
    void EMFBoxFluxAnalyzer::analyze(const art::Event& event) {

      art::Handle<StepPointMCCollection> hitsh;
      event.getByLabel(hitsModuleLabel_, hitsInstanceName_, hitsh);
      const StepPointMCCollection& hits(*hitsh);

      for(unsigned i=0; i<hits.size(); ++i) {

        const StepPointMC& hit = hits[i];
        VirtualDetectorId vid(hit.volumeId());

        switch(vid.id()) {
        case VirtualDetectorId::EMFBoxFront:
        case VirtualDetectorId::EMFBoxSW:
        case VirtualDetectorId::EMFBoxBottom:
        case VirtualDetectorId::EMFBoxBack:
        case VirtualDetectorId::EMFBoxNE:
        case VirtualDetectorId::EMFBoxTop:
          {
            if(hit.simParticle()) {

              art::Ptr<SimParticle> top = hit.simParticle();
              while(top->parent()) {
                top = top->parent();
              }

              // Found primary SimParticle.  It should have associated MARSInfo
              std::vector<art::Ptr<SimParticle> > vsp;
              vsp.push_back(top);

              art::FindOne<MARSInfo>
                mFinder(vsp, event,
                        art::InputTag(marsInfoModuleLabel_, marsInfoInstanceName_));

              const MARSInfo minfo = mFinder.at(0).ref();

              InputParticle inparticle(*extmon_, vid, hit, minfo,
                                       IO::G4JobInfo(event.id().run(), event.id().subRun(), event.id().event())
                                       );

              TransientParticle particle(inparticle, particleBeta(inparticle), cutMinTime_);
              histInputs_.fill(particle);
              if(isAccepted(particle)) {
                particles_.push_back(particle);
              }
            }
            else {
              std::cout<<"WARNING: missing SimParticle, skipping hit "<<hit
                       <<" in event "<<event.id()<<std::endl;
            }
          }
          break;

        default:
          std::cout<<"EMFBoxFluxAnalyzer: ignoring hit on VD "<<vid.id()<<" = "<<vid.name()<<std::endl;
          break;
        } // switch(vid)
      } // for(hits)
    } // analyze()

    //================================================================
    void EMFBoxFluxAnalyzer::endJob() {

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

      if(!noKNN_) {
        computeParticleRandomizations();
      }
      else {
        std::cout<<"EMFBoxFluxAnalyzer: NOT computing particle randomizations per user request."<<std::endl;
      }

      writeParticleNtuple();
    }


    //================================================================
    void EMFBoxFluxAnalyzer::writeParticleNtuple() {
      art::ServiceHandle<art::TFileService> tfs;

      IO::EMFBoxHit particle;
      IO::MARSInfo minfo;
      IO::G4JobInfo g4s1info;
      ParticleRandomization pr;

      TTree *nt = tfs->make<TTree>("vdhits", "VD hits and randomization parameters");
      nt->Branch("particle", &particle, particle.branchDescription());
      if(!noKNN_) {
        nt->Branch("randomization", &pr, pr.branchDescription());
      }
      nt->Branch("minfo", &minfo, minfo.branchDescription());
      nt->Branch("g4s1info", &g4s1info, g4s1info.branchDescription());

      // Sort by proton
      ProtonPathSort ps;
      std::sort(particles_.begin(), particles_.end(), ps);

      // Write out and count protons in the process
      unsigned numProtonPaths = 0;
      std::set<MARSInfo, CmpProtonId> uniqProtons;
      if(!particles_.empty()) {

        MARSInfo mcurrent = particles_[0].minfo;
        IO::G4JobInfo gcurrent = particles_[0].g4s1info;
        ++numProtonPaths;
        uniqProtons.insert(mcurrent);

        for(unsigned i=0; i<particles_.size(); ++i) {
          particle.emx = particles_[i].posExtMon.x();
          particle.emy = particles_[i].posExtMon.y();
          particle.emz = particles_[i].posExtMon.z();
          particle.mu2epx = particles_[i].momMu2e.x();
          particle.mu2epy = particles_[i].momMu2e.y();
          particle.mu2epz = particles_[i].momMu2e.z();
          particle.time = particles_[i].time;
          particle.pdgId = particles_[i].pdgId;
          particle.vdId = particles_[i].vd.id();

          pr = particles_[i].pr;

          minfo.info = particles_[i].minfo;

          g4s1info = particles_[i].g4s1info;

          nt->Fill();

          if(! (sameProtonAndSimPath(mcurrent, particles_[i].minfo) &&
                (gcurrent == particles_[i].g4s1info)
                )
             )
            {
              mcurrent = particles_[i].minfo;
              gcurrent = particles_[i].g4s1info;
              ++numProtonPaths;
              uniqProtons.insert(mcurrent);
            }
        } // for(particle)
      } // !empty

      std::cout<<"EMFBoxFluxAnalyzer: numProtonPaths = "<<numProtonPaths
               <<", numUniqProtons = "<<uniqProtons.size()
               <<" for "<<particles_.size()<<" stored hits"
               <<std::endl;
    }

    //================================================================
  } // namespace ExtMonFNAL
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::EMFBoxFluxAnalyzer);

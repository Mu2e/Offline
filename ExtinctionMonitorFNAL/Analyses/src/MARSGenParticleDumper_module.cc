// Andrei Gaponenko, 2012

#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>

#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "art/Framework/Principal/SelectorBase.h"
#include "canvas/Persistency/Provenance/BranchDescription.h"

#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/MARSInfo.hh"
#include "MCDataProducts/inc/MARSInfoCollection.hh"
#include "MCDataProducts/inc/GenParticleMARSAssns.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Utilities/inc/getCharge.hh"
#include "ExtinctionMonitorFNAL/Utilities/inc/EMFRandomizationSourceDefs.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

#include "TDirectory.h"
#include "TH1.h"
#include "TTree.h"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff) do {} while(0)

namespace mu2e {
  namespace ExtMonFNAL {

    using namespace Randomization;

    //================================================================
    namespace {

      struct Triple {
        float x;
        float y;
        float z;
        Triple() : x(), y(), z() {}
        Triple(const CLHEP::Hep3Vector& v) : x(v.x()), y(v.y()), z(v.z()) {}
      };

      //================================================================
      struct Direction {
        float theta;
        float phi;
        Direction() : theta(), phi() {}
      };

      //================================================================
      // Ntuple variables
      struct GPRecord {
        int run;
        int marsSubRun;
        int protonNumber;
        int particleIndex;
        float weight;

        int srcPlane;
        int pdgId;
        int charge;

        float time;
        float  pmag;

        Triple posExtMon;
        Triple posDump;
        Triple posMu2e;

        Triple momDump;

        Direction dirSrc; // theta==0 is the inside normal,
        // phi==0 is (-zdump) except for SourceFront, where it is (+ydump)

        Direction dirExtMon; // theta==0 is (-zextmon), phi==0 is (+xextmon)

        GPRecord()
          : run(), marsSubRun(), protonNumber(), particleIndex(), weight()
          , srcPlane(),  pdgId(), charge(), time(), pmag()
        {}
      };

      //================================================================
      Direction computeSrcDirection(SourceType srcPlane, const CLHEP::Hep3Vector& momDump) {
        CLHEP::HepRotation rot;

        switch(srcPlane) { // rotate momDump to get its z parallel to the inside normal
        case SourceFront: case SourceSignal : rot.rotateY(M_PI) ; break;
        case SourceBack  : /*do nothing*/ break;

        case SourceFloor : rot.rotateX(+M_PI/2.) ; break;
        case SourceCeiling : rot.rotateX(-M_PI/2.);  break;

        case SourceSouthWest : rot.rotateY(-M_PI/2.); break;
        case SourceNorthEast : rot.rotateY(+M_PI/2.); break;
        default:
          throw cet::exception("BADINPUTS")<<" computeSrcDirection: unknown srcPlane "<<srcPlane<<"\n";
        }

        const CLHEP::Hep3Vector localMom(rot * momDump);

        Direction res;
        res.theta = acos(localMom.z()/localMom.mag());
        res.phi = atan2(localMom.y(), localMom.x());
        return res;
      }

      //================================================================
      Direction computeExtMonDirection(const CLHEP::Hep3Vector& momExtMon) {
        CLHEP::HepRotation rot;
        rot.rotateY(M_PI);
        const CLHEP::Hep3Vector localMom(rot * momExtMon);
        Direction res;
        res.theta = acos(localMom.z()/localMom.mag());
        res.phi = atan2(localMom.y(), localMom.x());
        return res;
      }

    } // namespace {}

    //================================================================
    class MARSGenParticleDumper : public art::EDAnalyzer {
      // Geometry for coordinate conversion
      std::string geomModuleLabel_;

      // GenParticleCollection to operate on
      std::string inputModuleLabel_;
      std::string inputInstanceName_;

      // Coordinate conversion
      const ProtonBeamDump *dump_;
      const ExtMon *extmon_;

      SourcePlaneGeom srcGeom_;

      // Members needed to write the ntuple
      TTree *nt_;

      GPRecord gpnt_;

    public:
      explicit MARSGenParticleDumper(const fhicl::ParameterSet& pset);
      virtual void beginJob();
      virtual void beginRun(const art::Run& run);
      virtual void analyze(const art::Event& event);
    };

    //================================================================
    MARSGenParticleDumper::MARSGenParticleDumper(const fhicl::ParameterSet& pset)
      : art::EDAnalyzer(pset)
      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel", ""))
      , inputModuleLabel_(pset.get<std::string>("inputModuleLabel"))
      , inputInstanceName_(pset.get<std::string>("inputInstanceName"))
      , dump_(0)
      , extmon_(0)
      , srcGeom_(pset.get<fhicl::ParameterSet>("srcGeom"))
      , nt_(0)
    {}

    //================================================================
    void MARSGenParticleDumper::beginJob() {
      art::ServiceHandle<art::TFileService> tfs;
      nt_ = tfs->make<TTree>("nt", "EMF room ntuple");
      nt_->Branch("common", &gpnt_, "run/I:marsSubRun/I:protonNumber/I:particleIndex/i:weight/F:srcPlane/I:pdgId/I:charge/I:time/F:pmag/F");
      nt_->Branch("posExtMon", &gpnt_.posExtMon, "emx/F:emy/F:emz/F");
      nt_->Branch("posDump", &gpnt_.posDump, "dumpx/F:dumpy/F:dumpz/F");
      nt_->Branch("posMu2e", &gpnt_.posMu2e, "mu2ex/F:mu2ey/F:mu2ez/F");
      nt_->Branch("momDump", &gpnt_.momDump, "dumppx/F:dumppy/F:dumppz/F");
      nt_->Branch("dirSrc", &gpnt_.dirSrc, "srctheta/F:srcphi/F");
      nt_->Branch("dirExtMon", &gpnt_.dirExtMon, "emtheta/F:emphi/F");
    }

    //================================================================
    void MARSGenParticleDumper::beginRun(const art::Run& run) {

      if(geomModuleLabel_.empty()) {
        dump_ = &*GeomHandle<ProtonBeamDump>();
        extmon_ = &*GeomHandle<ExtMon>();
      }
      else {
        // Here we retrieve geometry objects from the
        // input file, as opposed to using GeometryService that might have
        // a different version than used to produce our input data.
        art::Handle<ProtonBeamDump> dumph;
        run.getByLabel(geomModuleLabel_, "", dumph);
        dump_ = &*dumph;

        art::Handle<ExtMonFNAL::ExtMon> deth;
        run.getByLabel(geomModuleLabel_, "", deth);
        extmon_ = &*deth;
      }

      std::cout<<"EMFCDumper::beginRun(): got ProtonBeamDump: coreCenterInMu2e = "
               <<dump_->coreCenterInMu2e()
               <<", coreRotY = "<<dump_->coreRotY()/CLHEP::degree<<" degrees"
               <<std::endl;

      std::cout<<"MARSGenParticleDumper::beginRun(): got ExtMonFNAL: detectorCenterInMu2e = "
               <<extmon_->detectorCenterInMu2e()
               <<", rot = "<<extmon_->detectorRotationInMu2e()
               <<std::endl;
    }

    //================================================================
    void MARSGenParticleDumper::analyze(const art::Event& event) {

      art::Handle<GenParticleCollection> particlesh;
      event.getByLabel(inputModuleLabel_, inputInstanceName_, particlesh);
      art::FindOne<MARSInfo> marsInfoFinder(particlesh, event, inputModuleLabel_);

      const GenParticleCollection& particles(*particlesh);
      for(unsigned i=0; i < particles.size(); ++i) {

        const MARSInfo& info = marsInfoFinder.at(i).ref();

        gpnt_.run = info.runNumber();
        gpnt_.marsSubRun = info.subRunNumber();
        gpnt_.protonNumber =  info.protonNumber();
        gpnt_.particleIndex = i;
        gpnt_.weight = info.weight();

        SourceType srcPlane = classifySource(particles[i].position(), srcGeom_, *dump_, *extmon_);

        gpnt_.srcPlane = srcPlane;
        gpnt_.pdgId = particles[i].pdgId();
        gpnt_.charge = getCharge(particles[i].pdgId());
        gpnt_.time = particles[i].time();;

        gpnt_.posExtMon = extmon_->mu2eToExtMon_position(particles[i].position());
        gpnt_.posDump = dump_->mu2eToBeamDump_position(particles[i].position());
        gpnt_.posMu2e = particles[i].position();

        const CLHEP::Hep3Vector momDump(dump_->mu2eToBeamDump_momentum(particles[i].momentum()));
        gpnt_.momDump = momDump;
        gpnt_.pmag = momDump.mag();

        gpnt_.dirSrc = computeSrcDirection(srcPlane, momDump);

        const CLHEP::Hep3Vector momExtMon(extmon_->mu2eToExtMon_momentum(particles[i].momentum()));
        gpnt_.dirExtMon = computeExtMonDirection(momExtMon);

        nt_->Fill();

      } // for()
    } // analyze(event)

    //================================================================

  } // namespace ExtMonFNAL
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::MARSGenParticleDumper);

//
// Finds and accepts generated photons that pair produce
// Largely based on Filters/src/InflightPionHits and CommonMC/src/StoppedParticleFinder

// Mu2e includes.
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/StatusG4.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"

// Framework includes.
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

// ROOT includes
#include "TTree.h"

// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <iostream>

using namespace std;

namespace mu2e {

  class GammaConversionPoints : public art::EDFilter {
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> g4ModuleLabel{Name("g4ModuleLabel"), Comment("Geant module label")};
      fhicl::Atom<std::string> simCollectionLabel{Name("simCollectionLabel"), Comment("Sim particle collection module label ("" to use geant label)"),""};
      fhicl::Atom<bool> doFilter{Name("doFilter"), Comment("Whether or not to filter passing events (true/false)"), true};
      fhicl::Atom<std::string> processCode{Name("processCode"), Comment("Photon creation code of interest"),  "mu2eExternalRMC"};
      fhicl::Atom<int> verbosity{Name("verbosity"), Comment("Verbose level"),  0};
      fhicl::Atom<double> rMin{Name("rMin"), Comment("Stop minimum radius in the DS (mm)"), -1.};
      fhicl::Atom<double> rMax{Name("rMax"), Comment("Stop maximum radius in the DS (mm)"),  1.e9};
      fhicl::Atom<double> zMin{Name("zMin"), Comment("Stop minimum z value (mm)"), -1.e9};
      fhicl::Atom<double> zMax{Name("zMax"), Comment("Stop maximum z value (mm)"),  1.e9};
      fhicl::Atom<double> xOffset{Name("SolenoidXOffset"), Comment("X coordinate offset for radius calculations (mm)"), -3904.};
      fhicl::Atom<double> eMin{Name("eMin"), Comment("Converted photon minimum energy (MeV)"), 0.};
      fhicl::Atom<bool> saveTree{Name("saveTree"), Comment("Save a tree of photon conversion information"), false};
    };
    typedef art::EDFilter::Table<Config> Parameters;

    explicit GammaConversionPoints(const Parameters& pset);
    virtual ~GammaConversionPoints() { }

    bool filter( art::Event& event);

    virtual bool beginRun(art::Run &run);
    bool beginSubRun(art::SubRun& sr) override;

    virtual void endJob();

  private:

    // Module label of the module that made the StepPointMCCollections.
    const std::string g4ModuleLabel_;
    const std::string simCollectionLabel_;
    const ProcessCode process_;

    const bool   doFilter_; //whether or not to filter non-conversions
    const int    verbosity_;
    const double rMin_; //spacial cuts
    const double rMax_;
    const double zMin_;
    const double zMax_;
    const double xoffset_; // detector solenoid x offset
    const double eMin_; //minimum energy for converted photon

    const bool saveTree_; // store an ntuple of photon conversion information
    TTree* ntup_;

    // Number of events that pass the filter.
    int nPassed_;

    //struct of stop info
    struct ConversionPointF {
      float x;
      float y;
      float z;
      float time;
      float px;
      float py;
      float pz;
      float genEnergy;
    };

    ConversionPointF data_;
  };

  GammaConversionPoints::GammaConversionPoints(const Parameters& pset):
    art::EDFilter{pset},
    g4ModuleLabel_(pset().g4ModuleLabel()),
    simCollectionLabel_(pset().simCollectionLabel()),
    process_(ProcessCode::findByName(pset().processCode())),
    doFilter_(pset().doFilter()),
    verbosity_(pset().verbosity()),
    rMin_(pset().rMin()),
    rMax_(pset().rMax()),
    zMin_(pset().zMin()),
    zMax_(pset().zMax()),
    xoffset_(pset().xOffset()),
    eMin_(pset().eMin()),
    saveTree_(pset().saveTree()),
    ntup_(nullptr),
    nPassed_(0) {
    if(saveTree_) {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "gammaStops" );
      ntup_ = tfdir.make<TTree>("conversions", "Photon Conversions");
      ntup_->Branch("stops", &data_, "x/F:y/F:z/F:time/F:px/F:py/F:pz/F:genEnergy/F");
    }
  }

  bool GammaConversionPoints::beginSubRun(art::SubRun& sr) {
    return true;
  }

  bool GammaConversionPoints::beginRun(art::Run& ){
    return true;
  }

  bool GammaConversionPoints::filter(art::Event& event) {
    art::Handle<StatusG4> g4StatusHandle;
    event.getByLabel( g4ModuleLabel_, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;

    // Accept only events with good status from G4.
    if ( g4Status.status() > 1 ) {
      if(verbosity_ > 1) printf("GammaConversionFilter::%s: Skipping due to G4 status\n", __func__);
      return false;
    }

    art::Handle<SimParticleCollection> simsHandle;
    event.getByLabel((simCollectionLabel_.size() > 0) ? simCollectionLabel_ : g4ModuleLabel_,simsHandle);
    SimParticleCollection const& sims(*simsHandle);

    // Loop through sim particles looking for the simulated photon of interest
    for(auto const& sim : sims) {
      SimParticle const& parent = sim.second;

      if(verbosity_ > 2) {
        printf("GammaConversionPoints::%s: Sim particle: pdg = %5i, creation = %3i, stopping = %3i, p_start = %5.1f, p_end = %5.1f\n",
               __func__, parent.pdgId(), int(parent.creationCode()), int(parent.stoppingCode()),
               parent.startMomentum().vect().rho(), parent.endMomentum().vect().rho());
      }

      // Check if the generated particle
      if(parent.creationCode() == process_) {
        if(verbosity_ > 1) printf("GammaConversionPoints::%s: Particle with selected process code is found\n", __func__);

        // Check if it's a pair conversion
        if(parent.pdgId() == 22 && ProcessCode::conv == parent.stoppingCode()) {

          // Get end position and momentum vectors
          CLHEP::Hep3Vector const& pos(parent.endPosition());

          // Check if it's within the region of interest
          const double r = sqrt((pos.x()-xoffset_)*(pos.x()-xoffset_) + pos.y()*pos.y());
          if(r < rMin_ || r > rMax_ || pos.z() < zMin_ || pos.z() > zMax_) {
            if(verbosity_ > 1) printf("GammaConversionPoints::%s: Particle with r = %.2f fails filter\n", __func__, r);
            return !doFilter_;
          }

          CLHEP::Hep3Vector const& p(parent.endMomentum().vect());
          CLHEP::Hep3Vector const& pStart(parent.startMomentum().vect()); //if needed for momentum
          CLHEP::Hep3Vector pstart(pStart.x(), pStart.y(), pStart.z()); //to edit the vector
          bool useStart = false;

          // Check if the momentum is stored, if not try to calculate it
          if(p.mag() == 0.) { //in case the final momentum not saved, add the daughters
            std::vector<art::Ptr<SimParticle> > const& daughters = parent.daughters();
            if(daughters.size() < 1) {
              useStart = true; //no daughters saved, use the start momentum
            } else {
              CLHEP::Hep3Vector pp(0.,0.,0.);
              int nConv = 0;
              for(unsigned int i = 0; i < daughters.size(); ++i) {
                SimParticle const sp = *daughters[i];
                if(sp.creationCode() == ProcessCode::conv) {
                  ++nConv;
                  pp += sp.startMomentum().vect();
                } else { //remove other daughters as happening before the conversion in case not all conversions found
                  pstart -= sp.startMomentum().vect();
                }
              }
              if(nConv < 2) { //both conversion daughters weren't stored
                useStart = true;
              }
              else { //both daughters stored
                data_.px = pp.x();
                data_.py = pp.y();
                data_.pz = pp.z();
              }
            }
          } else { //end momentum is stored
            data_.px = p.x();
            data_.py = p.y();
            data_.pz = p.z();
          }
          if(useStart) { //if conversion daughters not stored and end momentum not stored
            data_.px = pstart.x();
            data_.py = pstart.y();
            data_.pz = pstart.z();
          }

          // Filter on the photon energy
          const float energy = sqrt(pow(data_.px, 2) + pow(data_.py,2) + pow(data_.pz, 2));
          if(energy < eMin_) {
            if(verbosity_ > 1) printf("GammaConversionPoints::%s: Particle with E = %.2f MeV fails filter\n", __func__, energy);
            return !doFilter_;
          }

          // Accept the event

          // Store information for debugging/output tree
          data_.x = pos.x();
          data_.y = pos.y();
          data_.z = pos.z();
          data_.time = parent.endGlobalTime();
          data_.genEnergy = pStart.mag();

          ++nPassed_;

          // Print information about the event if requested
          if(verbosity_ > 0) {
            const float pt = sqrt(pow(data_.px, 2) + pow(data_.py,2));
            const float sz = pt/energy;
            const float cz = sqrt(1. - sz*sz) * ((data_.py < 0.) ? -1. : 1.);
            const float r  = sqrt(pow(data_.x - xoffset_, 2) + pow(data_.y,2));
            printf("GammaConversionPoints::%s: Accepted event: p = %5.1f, pT = %5.1f, cos(theta_z) = %6.3f, r = %5.1f, z = %6.1f, t = %5.0f, E(gen) = %5.1f\n",
                   __func__, energy, pt, cz, r, data_.z, data_.time, data_.genEnergy);
          }

          // Store the event if requested
          if(saveTree_) ntup_->Fill();

          return true; //converted so done looking
        } else {
          if(verbosity_ > 1) printf("GammaConversionPoints::%s: Particle process code %i fails filter\n", __func__, int(parent.stoppingCode()));
          return !doFilter_; //not a conversion
        }
      }
    }

    if(verbosity_ > 1) printf("GammaConversionPoints::%s: Event had no observed photon conversions\n", __func__);
    return !doFilter_;

  } // end of ::filter.

  void GammaConversionPoints::endJob() {
    mf::LogInfo("Summary")
      << "GammaConversionPoints_module: Number of conversion events found: "
      << nPassed_
      << "\n";
  }

} //end mu2e namespace

using mu2e::GammaConversionPoints;
DEFINE_ART_MODULE(GammaConversionPoints);

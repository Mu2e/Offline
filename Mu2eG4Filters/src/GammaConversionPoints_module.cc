//
// Largely based on Filters/src/InflightPionHits and Analyses/src/StoppedParticleFinder
// Finds and writes out information for generated photons that pair produce
// Writes conversion point x, y, z, t, px, py, pz, material (z, zeff, and fraction) , event weight, and gen energy

// Mu2e includes.
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "Mu2eUtilities/inc/PhysicalVolumeMultiHelper.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"

// G4 includes.
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"

// Framework includes.
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

// Root includes
#include "TH1F.h"
#include "TH1D.h"
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
      fhicl::Atom<std::string> generatorLabel{Name("generatorLabel"), Comment("Generator module label for event weights ("" to ignore)"),""};
      fhicl::Atom<int> defaultZ{Name("defaultZ"), Comment("Assume stop is in this given material Z"),  -1};
      fhicl::Atom<bool> doFilter{Name("doFilter"), Comment("Whether or not to filter passing events (true/false)"), false};
      fhicl::Atom<double> rMin{Name("rMin"), Comment("Stop minimum radius in the DS (mm)"), -1.};
      fhicl::Atom<double> rMax{Name("rMax"), Comment("Stop maximum radius in the DS (mm)"),  1.e9};
      fhicl::Atom<double> zMin{Name("zMin"), Comment("Stop minimum z value (mm)"), -1.e9};
      fhicl::Atom<double> zMax{Name("zMax"), Comment("Stop maximum z value (mm)"),  1.e9};
      fhicl::Atom<double> xOffset{Name("SolenoidXOffset"), Comment("X coordinate offset for radius calculations (mm)"), -3904.};
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
    std::string g4ModuleLabel_;
    std::string simCollectionLabel_;
    std::string generatorLabel_;

    int    defaultZ_; //default material Z to use if defined
    bool   doFilter_; //whether or not to filter non-conversions
    double rMin_; //spacial cuts
    double rMax_; //
    double zMin_; //
    double zMax_; //
    double xoffset_; // detector solenoid x offset 
    
    TH1F* hStoppingcode_;
    TH1D* hStoppingMat_;
    TTree* ntup_;

    // Number of events that pass the filter.
    int nPassed_;

    //struct of stop info
    IO::ConversionPointF data_;

    const PhysicalVolumeInfoMultiCollection *vols_;

  };

  GammaConversionPoints::GammaConversionPoints(const Parameters& pset):
    art::EDFilter{pset},
    g4ModuleLabel_(pset().g4ModuleLabel()),
    simCollectionLabel_(pset().simCollectionLabel()),
    generatorLabel_(pset().generatorLabel()),
    defaultZ_(pset().defaultZ()),
    doFilter_(pset().doFilter()),
    rMin_(pset().rMin()),
    rMax_(pset().rMax()),
    zMin_(pset().zMin()),
    zMax_(pset().zMax()),
    xoffset_(pset().xOffset()),
    hStoppingcode_(0),
    hStoppingMat_(0),
    ntup_(0),
    nPassed_(0) {
      if(defaultZ_ > 0)
	printf("GammaConversionPionts: Using Z=%i material for all stops!\n", defaultZ_);
      if(generatorLabel_.size() == 0)
	printf("GammaConversionPionts: No generator label passed, assuming all events have weight 1\n");
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "gammaStops" );
      
      hStoppingcode_  = tfdir.make<TH1F>( "hStoppingcode", "Photon stopping code",   200, 0.,  200.);
      hStoppingMat_  = tfdir.make<TH1D>( "hStoppingMat", "Photon stopping material",  1, 0.,  1.);
    
      ntup_ = tfdir.make<TTree>( "conversions", "Photon Conversions");
      ntup_->Branch("stops", &data_, IO::ConversionPointF::branchDescription().c_str());
    }  

  bool GammaConversionPoints::beginSubRun(art::SubRun& sr) {
    if(defaultZ_ > 0) return true; //assuming a material so don't need the volume info
    art::Handle<PhysicalVolumeInfoMultiCollection> volh;
    sr.getByLabel(g4ModuleLabel_, volh);
    vols_ = &*volh;
    return true;
  }

  bool GammaConversionPoints::beginRun(art::Run& ){

    return true;
  }

  bool
  GammaConversionPoints::filter(art::Event& event) {
    
    art::Handle<StatusG4> g4StatusHandle;
    event.getByLabel( g4ModuleLabel_, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;

    // Accept only events with good status from G4.
    if ( g4Status.status() > 1 ) {
      return false;
    }

    art::Handle<SimParticleCollection> simsHandle;
    event.getByLabel((simCollectionLabel_.size() > 0) ? simCollectionLabel_ : g4ModuleLabel_,simsHandle);
    SimParticleCollection const& sims(*simsHandle);
    
    float weight = (generatorLabel_.size() > 0) ? event.getValidHandle<EventWeight>(generatorLabel_)->weight() : 1.;
    // printf("weight = %.4f\n", weight);
    //loop through sim particles looking for gen particle
    for(auto const& sim : sims) {
      //check if the generated particle
      if(sim.second.creationCode() == ProcessCode::mu2ePrimary) {
	SimParticle const& parent = sim.second;
	hStoppingcode_->Fill(parent.stoppingCode());
	//check if pair conversion
	if(parent.pdgId() == 22 && ProcessCode::conv == parent.stoppingCode()) {
	  //get end position and momentum vectors
	  CLHEP::Hep3Vector const&       pos(parent.endPosition());
	  //check if it's within the region of interest
	  double r = sqrt((pos.x()-xoffset_)*(pos.x()-xoffset_) + pos.y()*pos.y());
	  if(r < rMin_ || r > rMax_ || pos.z() < zMin_ || pos.z() > zMax_) return !doFilter_;

	  CLHEP::Hep3Vector const&       p(parent.endMomentum().vect());
	  CLHEP::Hep3Vector const&       pStart(parent.startMomentum().vect()); //for RMC weight energy or if needed for momentum
	  CLHEP::Hep3Vector pstart(pStart.x(), pStart.y(), pStart.z()); //to edit the vector
	  bool useStart = false;
	  
	  if(p.mag() == 0.) { //in case the final momentum not saved, add the daughters
	    std::vector<art::Ptr<SimParticle> > const& daughters = parent.daughters();
	    if(daughters.size() < 1) {
	      useStart = true; //no daughters saved
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

	  data_.x = pos.x();
	  data_.y = pos.y();
	  data_.z = pos.z();
	  data_.time = parent.endGlobalTime();
	  data_.weight = weight;
	  data_.genEnergy = pStart.mag();

	  //get material properties
	  if(defaultZ_ <= 0) { //find actual Z of the material
	    PhysicalVolumeMultiHelper vi(*vols_);
	    PhysicalVolumeInfo const& endVolParent = vi.endVolume(parent);
	    const std::string& materialName = endVolParent.materialName();
	    hStoppingMat_->Fill(materialName.c_str(), 1.);
	    const G4Material* material = findMaterialOrThrow(materialName);
	    const int nElements = material->GetNumberOfElements();
	    const G4ElementVector* elementVector = material->GetElementVector();
	    const double* fractionVector = material->GetFractionVector();
	    data_.matN = std::min(nElements, (int) IO::kMaxConversionMaterialElements);
	    for(int index = 0; index < data_.matN; ++index) {
	      float matZeff((*elementVector)[index]->GetZ()), matFrac(fractionVector[index]);
	      data_.matZ   [index] = (*elementVector )[index]->GetZasInt();
	      data_.matZeff[index] = matZeff;
	      data_.matFrac[index] = matFrac;
	    }
	  } else {
	    data_.matN = 1;
	    data_.matZ   [0] = defaultZ_;
	    data_.matZeff[0] = defaultZ_;
	    data_.matFrac[0] = 1.;
	  }

	  ntup_->Fill();
	  ++nPassed_; 
	  return true; //converted so done looking
	} else return !doFilter_; //not a conversion
      } 
    }

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



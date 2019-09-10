// Ntuple dumper for Primary Proton Energy Deposition.
//
// Andrei Gaponenko, 2013

#include <string>
#include <vector>
#include <limits>
#include <cmath>

#include "cetlib_except/exception.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "TDirectory.h"
#include "TH1.h"
#include "TTree.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "Mu2eUtilities/inc/SimParticleGetTau.hh"


#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "DataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"

// CLHEP includes.
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

namespace mu2e {

  //================================================================
  double getCharge(PDGCode::type pdgId) {
    // unlike generic conditions, MC particle data
    // should not change run-to-run, so static is safe
    // use static for efficiency
    static GlobalConstantsHandle<ParticleDataTable> pdt;

    ParticleDataTable::maybe_ref info = pdt->particle(pdgId);

    if(!info.isValid()) {
      throw cet::exception("MISSINGINFO")<<"No valid PDG info for pdgId = "<<pdgId<<"\n";
    }

    return info.ref().charge();
  }

  //================================================================
  double getKineticEnergy(const StepPointMC& hit) {
    // unlike generic conditions, MC particle data
    // should not change run-to-run, so static is safe
    // use static for efficiency
    static GlobalConstantsHandle<ParticleDataTable> pdt;

    ParticleDataTable::maybe_ref info = pdt->particle(hit.simParticle()->pdgId());

    if(!info.isValid()) {
      throw cet::exception("MISSINGINFO")<<"No valid PDG info for hit = "<<hit<<"\n";
    }

    const double mass = info.ref().mass();
    return sqrt(hit.momentum().mag2() + std::pow(mass, 2)) - mass;
  }

  //================================================================
  struct VDHit {
    float x;
    float y;
    float z;
    float time;

    float px;
    float py;
    float pz;
    float pmag;
    float ek;

    float totalEDep;
    float nonIonizingEDep;

    float charge;
    int   pdgId;
    unsigned particleId;
    unsigned volumeCopyNumber;


    VDHit() : x(std::numeric_limits<double>::quiet_NaN())
            , y(std::numeric_limits<double>::quiet_NaN())
            , z(std::numeric_limits<double>::quiet_NaN())

            , time(std::numeric_limits<double>::quiet_NaN())

            , px(std::numeric_limits<double>::quiet_NaN())
            , py(std::numeric_limits<double>::quiet_NaN())
            , pz(std::numeric_limits<double>::quiet_NaN())
            , pmag(std::numeric_limits<double>::quiet_NaN())
            , ek(std::numeric_limits<double>::quiet_NaN())
            , totalEDep(std::numeric_limits<double>::quiet_NaN())
            , nonIonizingEDep(std::numeric_limits<double>::quiet_NaN())

      , charge(std::numeric_limits<double>::quiet_NaN())
      , pdgId(0)
      , particleId(-1U)
      , volumeCopyNumber(-1U)
    {}

    //----------------------------------------------------------------
    VDHit(const SimParticleTimeOffset& toff, const StepPointMC& hit)
      : x(hit.position().x())
      , y(hit.position().y())
      , z(hit.position().z())

      , time(toff.timeWithOffsetsApplied(hit))

      , px(hit.momentum().x())
      , py(hit.momentum().y())
      , pz(hit.momentum().z())

      , pmag(hit.momentum().mag())
      , ek(getKineticEnergy(hit))

      , totalEDep(hit.totalEDep())
      , nonIonizingEDep(hit.nonIonizingEDep())

      , charge(getCharge(hit.simParticle()->pdgId()))

      , pdgId(hit.simParticle()->pdgId())
      , particleId(hit.simParticle()->id().asUint())
      , volumeCopyNumber(hit.volumeId())
    {}

  }; // struct VDHit

  //================================================================
  class PrimaryProtonEnergyDumper : public art::EDAnalyzer {
    typedef std::vector<std::string> VS;
    typedef std::vector<StepPointMCCollection> VspMC;

    art::InputTag hitsInputTag_;
    SimParticleTimeOffset toff_;

    bool writeProperTime_;
    VS tauHitCollections_;
    std::vector<int> decayOffCodes_;

    // Members needed to write the ntuple
    TTree *nt_;
    VDHit hit_;
    float tau_;

    TH1F* _hEnergyVsZ;
    TH1F* _hEnergyVsR;
    TH1F* _hHitX;
    TH1F* _hHitY;
    TH1F* _hHitZ;

    TH1F* _hHit_Fin1_X;
    TH1F* _hHit_Fin1_Y;
    TH1F* _hHit_Fin1_Z;
    TH1F* _hHit_Fin2_X;
    TH1F* _hHit_Fin2_Y;
    TH1F* _hHit_Fin2_Z;
    TH1F* _hHit_Fin3_X;
    TH1F* _hHit_Fin3_Y;
    TH1F* _hHit_Fin3_Z;
    TH1F* _hHit_Fin4_X;
    TH1F* _hHit_Fin4_Y;
    TH1F* _hHit_Fin4_Z;

    CLHEP::Hep3Vector _gunOrigin;
    CLHEP::HepRotation _gunRotation;
    double targetHalfLength;

    // workaround for geometry service not being available at job start
    bool booked{false};

  public:
    explicit PrimaryProtonEnergyDumper(const fhicl::ParameterSet& pset);
    virtual void beginJob();
    virtual void beginRun(art::Run const& run);
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  PrimaryProtonEnergyDumper::PrimaryProtonEnergyDumper(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , hitsInputTag_(pset.get<std::string>("hitsInputTag"))
    , toff_(pset.get<fhicl::ParameterSet>("TimeOffsets"))
    , writeProperTime_(pset.get<bool>("writeProperTime", false))
    , tauHitCollections_( writeProperTime_ ? pset.get<VS>("tauHitCollections") : VS() )
    , nt_(0)
    , tau_()
  {
    if(writeProperTime_) {
      decayOffCodes_ = pset.get<std::vector<int> >("decayOffPDGCodes");
      // must sort to use binary_search in SimParticleGetTau
      std::sort(decayOffCodes_.begin(), decayOffCodes_.end());
    }
  }

  //================================================================
  void PrimaryProtonEnergyDumper::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;

    static const char branchDesc[] = "x/F:y/F:z/F:time/F:px/F:py/F:pz/F:pmag/F:ek/F:totalEDep/F:nonIonizingEDep/F:charge/F:pdgId/I:particleId/i:volumeCopy/i";
    nt_ = tfs->make<TTree>( "nt", "PrimaryProtonEnergyDumper ntuple");

    nt_->Branch("hits", &hit_, branchDesc);
    if(writeProperTime_) {
      nt_->Branch("tau", &tau_, "tauNormalized/F");
    }
  }

  void PrimaryProtonEnergyDumper::beginRun(art::Run const& run){
    art::ServiceHandle<art::TFileService> tfs;

    if (!booked) // a workaround for geometry service not being available at job start
      //
      // geometry service to get initial proton info
      art::ServiceHandle<GeometryService> geom;

    _gunRotation = GeomHandle<ProductionTarget>()->protonBeamRotation();
    _gunOrigin = GeomHandle<ProductionTarget>()->position() + _gunRotation*CLHEP::Hep3Vector(0., 0., GeomHandle<ProductionTarget>()->halfLength());

    std::cout << "gun origin pieces in analysis module \n " << 
      GeomHandle<ProductionTarget>()->position()  << "\n"<<
      _gunRotation*CLHEP::Hep3Vector(0., 0., GeomHandle<ProductionTarget>()->halfLength()) << "\n" <<
      _gunOrigin << std::endl;
    _hEnergyVsZ = tfs->make<TH1F>("_hEnergyVsZ","Energy vs Z",300,-6300.,-6000.);
    _hHitX = tfs->make<TH1F>("_hHitX","Internal X Position of Hit, Energy Weighted",100,-10.,10.);
    _hHitY = tfs->make<TH1F>("_hHitY","Internal Y Position Of Hit, Energy Weighted",100,-10.,10.);
    std::cout << " half length = " << GeomHandle<ProductionTarget>()->halfLength() << std::endl;
    Int_t nbins = +2.0*GeomHandle<ProductionTarget>()->halfLength() + 0.5;
    std::cout << " nbins = " << nbins << std::endl;
    _hHitZ = tfs->make<TH1F>("_hHitZ","Internal Z Position of Hit, Energy Weighted",nbins+10,-10., +2.0*GeomHandle<ProductionTarget>()->halfLength());

    _hHit_Fin1_Z = tfs->make<TH1F>("_hHitZFin1","Internal Z Position of Hit, Energy Weighted, Fin1",nbins+10,-10., +2.0*GeomHandle<ProductionTarget>()->halfLength());

    _hHit_Fin2_Z = tfs->make<TH1F>("_hHitZFin2","Internal Z Position of Hit, Energy Weighted, Fin2",nbins+10,-10., +2.0*GeomHandle<ProductionTarget>()->halfLength());

    _hHit_Fin3_Z = tfs->make<TH1F>("_hHitZFin3","Internal Z Position of Hit, Energy Weighted, Fin 3",nbins+10,-10., +2.0*GeomHandle<ProductionTarget>()->halfLength());

    _hHit_Fin4_Z = tfs->make<TH1F>("_hHitZFin4","Internal Z Position of Hit, Energy Weighted, Fin 4",nbins+10,-10., +2.0*GeomHandle<ProductionTarget>()->halfLength());


    booked = true;
  }



  //================================================================
  void PrimaryProtonEnergyDumper::analyze(const art::Event& event) {
    toff_.updateMap(event);

    VspMC spMCColls;
    for ( const auto& iColl : tauHitCollections_ ){
      auto spColl = event.getValidHandle<StepPointMCCollection>(iColl);
      spMCColls.push_back( *spColl );
    }


    //    std::cout << "hitsInputTag_ " << hitsInputTag_ << std::endl;
    //std::cout << "hitsInputTag instance " << hitsInputTag_.instance() << std::endl;

    std::string hitInputTagInstance = hitsInputTag_.instance();
    const auto& ih = event.getValidHandle<StepPointMCCollection>(hitsInputTag_);

    for(const auto& i : *ih) {


      hit_ = VDHit(toff_, i);
      if (hit_.totalEDep > 0){_hEnergyVsZ->Fill(hit_.z,hit_.totalEDep);}
      CLHEP::Hep3Vector hitLoc(hit_.x,hit_.y,hit_.z);

      //
      // this rotation takes me from mu2e coordinates to internal  
      // at generation time  we apply _gunRotation to go from target coord to mu2e coord; now apply inverse
           CLHEP::Hep3Vector hitPositionInCore = _gunRotation.inverse()*(hitLoc - _gunOrigin);
      
      _hHitX->Fill(hitPositionInCore.x(),hit_.totalEDep);
      _hHitY->Fill(hitPositionInCore.y(),hit_.totalEDep);
      //
      // - sign since beam travels toward negative z in Mu2e coordinates.  make plot run from zero and look like the target...
      if (hitInputTagInstance == "ProductionTargetCore") {
      _hHitZ->Fill(-hitPositionInCore.z(),hit_.totalEDep);
      } else if (hitInputTagInstance == "ProductionTargetFin1"){
      _hHit_Fin1_Z->Fill(-hitPositionInCore.z(),hit_.totalEDep);
     } else if (hitInputTagInstance == "ProductionTargetFin2"){
       _hHit_Fin2_Z->Fill(-hitPositionInCore.z(),hit_.totalEDep);
     } else if (hitInputTagInstance == "ProductionTargetFin3"){
       _hHit_Fin3_Z->Fill(-hitPositionInCore.z(),hit_.totalEDep);
     } else if (hitInputTagInstance == "ProductionTargetFin4"){
        _hHit_Fin4_Z->Fill(-hitPositionInCore.z(),hit_.totalEDep);
      }

      if(writeProperTime_) {
        tau_ = SimParticleGetTau::calculate(i, spMCColls, decayOffCodes_);
      }

      //      nt_->Fill();
    }

  } // analyze(event)

    //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::PrimaryProtonEnergyDumper);

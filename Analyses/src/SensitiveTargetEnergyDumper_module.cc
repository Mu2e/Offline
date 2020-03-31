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
#include "TH2F.h"
#include "TNtuple.h"

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
  class SensitiveTargetEnergyDumper : public art::EDAnalyzer {
    typedef std::vector<std::string> VS;
    typedef std::vector<StepPointMCCollection> VspMC;

    art::InputTag hitsInputTag_;
    SimParticleTimeOffset toff_;

    bool writeProperTime_;
    VS tauHitCollections_;
    std::vector<int> decayOffCodes_;

    // Members needed to write the ntuple
    TNtuple *nt_;
    TNtuple *ntAll_;
    VDHit hit_;
    float tau_;

    TH1F* _hEnergyVsZ;
    TH1F* _hEnergyVsR;
    TH1F* _hHitX;
    TH1F* _hHitY;
    TH1F* _hHitZ;

    TH1F* _hHitXAll;
    TH1F* _hHitYAll;
    TH1F* _hEnergyVsZAll;

    Float_t ntMembers[4];
    TH1F* _hHitZCore;
    TH1F* _hHitZStartingCore;

    TH1F* _hHitZFin;
    TH1F* _hHitZStartingFin;

    CLHEP::Hep3Vector hitPositionInternal;

    TH2F* _hHitPosRing;
    TH2F* _hHitNegRing;

    CLHEP::Hep3Vector _gunOrigin;
    CLHEP::HepRotation _gunRotation;
    double targetHalfLength;

    int numberOfCalls{0};
    int nInCore{0};
    int nTotal{0};
    // workaround for geometry service not being available at job start
    bool booked{false};

  public:
    explicit SensitiveTargetEnergyDumper(const fhicl::ParameterSet& pset);
    virtual void beginJob();
    virtual void beginRun(art::Run const& run);
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  SensitiveTargetEnergyDumper::SensitiveTargetEnergyDumper(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , hitsInputTag_(pset.get<std::string>("hitsInputTag"))
    , toff_(pset.get<fhicl::ParameterSet>("TimeOffsets"))
    , writeProperTime_(pset.get<bool>("writeProperTime", false))
    , tauHitCollections_( writeProperTime_ ? pset.get<VS>("tauHitCollections") : VS() )
    , tau_()
  {
    if(writeProperTime_) {
      decayOffCodes_ = pset.get<std::vector<int> >("decayOffPDGCodes");
      // must sort to use binary_search in SimParticleGetTau
      std::sort(decayOffCodes_.begin(), decayOffCodes_.end());
    }
  }

  //================================================================
  void SensitiveTargetEnergyDumper::beginJob() {
   art::ServiceHandle<art::TFileService> tfs;

    //    static const char branchDesc[] = "x/F:y/F:z/F:time/F:totalEdep/F";
   //   nt_ = tfs->make<TNtuple>( "nt", "sensitiveTargetNtuple","x:y:z:totalEdep");

    ntAll_ = tfs->make<TNtuple>( "ntAll", "sensitiveTargetNtupleAllVolumes","x:y:z:totalEdep");

  }
    //   nt_->Branch("hits", &hit_, branchDesc);
    //if(writeProperTime_) {
    //  nt_->Branch("tau", &tau_, "tauNormalized/F");
    // }
    // }

    void SensitiveTargetEnergyDumper::beginRun(art::Run const& run){
      art::ServiceHandle<art::TFileService> tfs;

      if (!booked) // a workaround for geometry service not being available at job start
	//
	// geometry service to get initial proton info
	art::ServiceHandle<GeometryService> geom;

      _gunRotation = GeomHandle<ProductionTarget>()->protonBeamRotation();
      _gunOrigin = GeomHandle<ProductionTarget>()->haymanPosition();

      std::cout << "gun origin pieces in analysis module \n " << 
	GeomHandle<ProductionTarget>()->haymanPosition()  << "\n"<<
	_gunRotation*CLHEP::Hep3Vector(0., 0., GeomHandle<ProductionTarget>()->halfHaymanLength()) << "\n" <<
	_gunOrigin << std::endl;

      Int_t nbins = +2.0*GeomHandle<ProductionTarget>()->halfHaymanLength() + 0.5;
      std::cout << " nbins = " << nbins << std::endl;

      _hEnergyVsZ = tfs->make<TH1F>("_hEnergyVsZ","Energy vs Z",nbins+10
				    ,0.,+2.0*GeomHandle<ProductionTarget>()->halfHaymanLength());
      //300,-6300.,-6000.);
      _hHitX = tfs->make<TH1F>("_hHitX","Internal X Position of Hit, Energy Weighted",200,-20.,20.);
      _hHitY = tfs->make<TH1F>("_hHitY","Internal Y Position Of Hit, Energy Weighted",200,-20.,20.);
      std::cout << " half length = " << GeomHandle<ProductionTarget>()->halfHaymanLength() << std::endl;
 
      _hHitZCore = tfs->make<TH1F>("_hHitZCore","Internal Z Position of Hit, Core Section, Energy Weighted",nbins+10
				   ,-10.,+2.0*GeomHandle<ProductionTarget>()->halfHaymanLength());

      _hHitZStartingCore = tfs->make<TH1F>("_hHitZStartingCore","Internal Z Position of Hit, Starting Core Section, Energy Weighted",nbins+10
					   ,-10.,+2.0*GeomHandle<ProductionTarget>()->halfHaymanLength());

      _hHitZFin = tfs->make<TH1F>("_hHitZFin","Internal Z Position of Hit, Fin Section, Energy Weighted",nbins+10
				  ,-10.,+2.0*GeomHandle<ProductionTarget>()->halfHaymanLength());

      _hHitZStartingFin = tfs->make<TH1F>("_hHitZStartingFin","Internal Z Position of Hit, Starting Fin Section, Energy Weighted",nbins+10
					  ,-10.,+2.0*GeomHandle<ProductionTarget>()->halfHaymanLength());


      _hHitNegRing = tfs->make<TH2F>("_hHitNegRing","Scatter Plot for Ring at Beginning of Target",50,-25.,25.,50,-25.,25.);
      _hHitPosRing = tfs->make<TH2F>("_hHitPosRing","Scatter Plot for Ring at End of Target",50,-25.,25.,50,-25.,25.);

      _hEnergyVsZAll = tfs->make<TH1F>("_hEnergyVsZAll","Energy vs Z",nbins+10
				    ,0.,+2.0*GeomHandle<ProductionTarget>()->halfHaymanLength());
      _hHitXAll = tfs->make<TH1F>("_hHitXAll","Internal X Position of Hit, Energy Weighted",200,-20.,20.);
      _hHitYAll = tfs->make<TH1F>("_hHitYAll","Internal Y Position Of Hit, Energy Weighted",200,-20.,20.);

      booked = true;
    }
  



  //================================================================
  void SensitiveTargetEnergyDumper::analyze(const art::Event& event) {
    // toff_.updateMap(event);

    /*    VspMC spMCColls;
    for ( const auto& iColl : tauHitCollections_ ){
      auto spColl = event.getValidHandle<StepPointMCCollection>(iColl);
      spMCColls.push_back( *spColl );
    }
    */
    bool useThisInstance = false;
 
    ++numberOfCalls;
       std::cout << "number Of Calls = " << numberOfCalls << std::endl;
    std::string hitInputTagInstance = hitsInputTag_.instance();
 
          std::cout << "hitInputTagInstance " << hitInputTagInstance << " " << useThisInstance << std::endl;
 
    // 
    // do we want this instance?
    if (hitInputTagInstance.find("ProductionTarget") != std::string::npos) { 
      useThisInstance = true;
                 std::cout << "hitInputTagInstance Found " << hitInputTagInstance << " " << useThisInstance << std::endl;
    }


    const auto& ih = event.getValidHandle<StepPointMCCollection>(hitsInputTag_);

    for(const auto& i : *ih) {


      hit_ = VDHit(toff_, i);
      CLHEP::Hep3Vector hitLoc(hit_.x,hit_.y,hit_.z);

           std::cout << "hit x = " << hitLoc << std::endl;
      //     if (hit_.totalEDep > 0){_hEnergyVsZ->Fill(hit_.z,hit_.totalEDep);}
  
      //
      // this rotation takes me from mu2e coordinates to internal  
      // at generation time  we apply _gunRotation to go from target coord to mu2e coord.  A little tricky
      // since gunOrigin is defined to be the downstream end of the gun for the rotated target, that is, where the protons hit. Let's transform that away.
 
      hitPositionInternal = _gunRotation.inverse()*(hitLoc - _gunOrigin) - CLHEP::Hep3Vector(0.,0., GeomHandle<ProductionTarget>()->halfHaymanLength());
      //      std::cout << "hitloc, rotation, core = " << hitLoc << "\n" << _gunOrigin << "\n" << _gunRotation << "\n" << hitPositionInternal << std::endl;
      //	   std::cout << " x val " << hitLoc << std::endl;
      ntMembers[0] = hitPositionInternal.x();
      ntMembers[1] = hitPositionInternal.y();
      ntMembers[2] = hitPositionInternal.z();
      ntMembers[3] = hit_.totalEDep;
      if (useThisInstance){
      _hHitXAll->Fill(hitPositionInternal.x(),hit_.totalEDep);
      _hHitYAll->Fill(hitPositionInternal.y(),hit_.totalEDep);
      _hEnergyVsZAll->Fill(-hitPositionInternal.z(),hit_.totalEDep);
      ntAll_->Fill(ntMembers);
      ++nTotal;
            std::cout << "hitInputTagInstance " << hitInputTagInstance << " " << useThisInstance << " " << nTotal << std::endl;
      }
      //
      // - sign since beam travels toward negative z in Mu2e coordinates.  make plot run from zero and look like the target...
      if (hitInputTagInstance == "ProductionTargetCoreSection") {
	_hHitZCore->Fill(-hitPositionInternal.z(),hit_.totalEDep);
	_hHitX->Fill(hitPositionInternal.x(),hit_.totalEDep);
	_hHitY->Fill(hitPositionInternal.y(),hit_.totalEDep);
	_hEnergyVsZ->Fill(-hitPositionInternal.z(),hit_.totalEDep);
	++nInCore;
	//	std::cout << " in core section " << nInCore << std::endl; 
     } else if (hitInputTagInstance == "ProductionTargetPositiveEndRing"){
	//	std::cout << "in pos ring" << std::endl;
	//  	   std::cout << " x val " << hitLoc << std::endl;
	_hHitPosRing->Fill(hitPositionInternal.x(),hitPositionInternal.y(),hit_.totalEDep);
	_hHitX->Fill(hitPositionInternal.x(),hit_.totalEDep);
	_hHitY->Fill(hitPositionInternal.y(),hit_.totalEDep);
	_hEnergyVsZ->Fill(-hitPositionInternal.z(),hit_.totalEDep);
      } else if (hitInputTagInstance == "ProductionTargetNegativeEndRing"){
	//	std::cout << "in neg ring" << std::endl;
	//   	   std::cout << " x val " << hitLoc << std::endl;
	_hHitNegRing->Fill(hitPositionInternal.x(),hitPositionInternal.y(),hit_.totalEDep);
	_hHitX->Fill(hitPositionInternal.x(),hit_.totalEDep);
	_hHitY->Fill(hitPositionInternal.y(),hit_.totalEDep);
	_hEnergyVsZ->Fill(-hitPositionInternal.z(),hit_.totalEDep);
      } else if (hitInputTagInstance == "ProductionTargetStartingCoreSection"){
	//	std::cout << "in  starting core" << std::endl;
	//  	   std::cout << " x val " << hitLoc << std::endl;
	_hHitZStartingCore->Fill(-hitPositionInternal.z(),hit_.totalEDep);
	_hHitX->Fill(hitPositionInternal.x(),hit_.totalEDep);
	_hHitY->Fill(hitPositionInternal.y(),hit_.totalEDep);
	_hEnergyVsZ->Fill(-hitPositionInternal.z(),hit_.totalEDep);
      } else if (hitInputTagInstance == "ProductionTargetFinSection"){
	//	std::cout << "in  starting core" << std::endl;
	//  	   std::cout << " x val " << hitLoc << std::endl;
	_hHitZFin->Fill(-hitPositionInternal.z(),hit_.totalEDep);
	_hHitX->Fill(hitPositionInternal.x(),hit_.totalEDep);
	_hHitY->Fill(hitPositionInternal.y(),hit_.totalEDep);
	_hEnergyVsZ->Fill(-hitPositionInternal.z(),hit_.totalEDep);
      } else if (hitInputTagInstance == "ProductionTargetFinStartingSection"){
	//	std::cout << "in  starting core" << std::endl;
	//  	   std::cout << " x val " << hitLoc << std::endl;
	_hHitZStartingFin->Fill(-hitPositionInternal.z(),hit_.totalEDep);
	_hHitX->Fill(hitPositionInternal.x(),hit_.totalEDep);
	_hHitY->Fill(hitPositionInternal.y(),hit_.totalEDep);
	_hEnergyVsZ->Fill(-hitPositionInternal.z(),hit_.totalEDep);
      } else if (hitInputTagInstance == "ProductionTargetFinTopSection" || hitInputTagInstance == "ProductionTargetFinTopStartingSection"){
	//	std::cout << "in  starting core" << std::endl;
	//  	   std::cout << " x val " << hitLoc << std::endl;
	_hHitZStartingFin->Fill(-hitPositionInternal.z(),hit_.totalEDep);
	_hHitX->Fill(hitPositionInternal.x(),hit_.totalEDep);
	_hHitY->Fill(hitPositionInternal.y(),hit_.totalEDep);
	_hEnergyVsZ->Fill(-hitPositionInternal.z(),hit_.totalEDep);
      }

      //      if(writeProperTime_) {
      //	tau_ = SimParticleGetTau::calculate(i, spMCColls, decayOffCodes_);
      //      }
    }

  }
  // analyze(event)

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::SensitiveTargetEnergyDumper);

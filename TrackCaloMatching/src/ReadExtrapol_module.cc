//
//
//
// $Id: ReadExtrapol_module.cc,v 1.14 2014/03/04 14:36:18 murat Exp $
// $Author: murat $
// $Date: 2014/03/04 14:36:18 $
//
// Original author G. Pezzullo
//

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Geometry/HepPoint.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "KalmanTests/inc/KalRepCollection.hh"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//tracker includes
#include "BaBar/BaBar.hh"
#include "BaBar/Constants.hh"
#include "TrkBase/TrkRep.hh"
#include "TrkBase/HelixParams.hh"
#include "TrkBase/HelixTraj.hh"
#include "KalmanTrack/KalRep.hh"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "TrackCaloMatching/inc/TrkToCaloExtrapolCollection.hh"

#include "CaloCluster/inc/CaloClusterUtilities.hh"

// BaBar
#include "BaBar/BaBar.hh"
#include "BaBar/include/TrkBase/TrkDifTraj.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/KalRepCollection.hh"
#include "KalmanTests/inc/TrkFitDirection.hh"

//calorimeter includes
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "TrackCaloMatching/inc/TrackClusterLink.hh"

// Other includes.
#include "cetlib/exception.h"


// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/VisibleGenElTrack.hh"
#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"

//root includes
#include "TFile.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TMath.h"

// From the art tool-chain
#include <cmath>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <functional>
#include "cetlib/pow.h"


using namespace std;
using cet::square;
using cet::sum_of_squares;

namespace mu2e {
  double thetaWimpact(const CLHEP::Hep3Vector& mom, int vaneId){

    art::ServiceHandle<GeometryService> geom;
    if(! geom->hasElement<VaneCalorimeter>() ) return 0;
    GeomHandle<VaneCalorimeter> cg;
    Vane const v1 = cg->vane(vaneId);
    CLHEP::Hep3Vector dirMom_rotated = (v1.rotation())*mom.unit();

    if(std::fabs(dirMom_rotated.getX() ) < 1e-10){
      if(dirMom_rotated.getZ()>0.0) {
	return 90.0;
      }else{
	return -90.0;
      }
    }
    double thW = 0.0;
    thW = std::atan(-1.0*dirMom_rotated.getZ() / dirMom_rotated.getX() ) ;
    thW *= Constants::radToDegrees;
    return thW;
  }
  double thetaVimpact(const CLHEP::Hep3Vector& mom, int vaneId){//(FIXME)

    art::ServiceHandle<GeometryService> geom;
    if(! geom->hasElement<VaneCalorimeter>() ) return 0;
    GeomHandle<VaneCalorimeter> cg;
    Vane const v1 = cg->vane(vaneId);
    CLHEP::Hep3Vector dirMom_rotated = (v1.rotation())*mom.unit();

    if(std::fabs(dirMom_rotated.getX() ) < 1e-10){
      if(dirMom_rotated.getY()>0.0) {
	return 90.0;
      }else{
	return -90.0;
      }
    }
    double thV = 0.0;
    thV = std::atan(-1.0*dirMom_rotated.getY() / dirMom_rotated.getX() ) ;
    thV *= Constants::radToDegrees;
    return thV;
  }

  float thetaWimpactErr(float& z, float& dz, float& x, float& dx){
    float res = 0.0;
    if(x == 0.0) return res;
    double invDen = 1.0/(x*x + z*z);

    res += dx*dx * z*z * invDen*invDen;
    res += dz*dz * x*x * invDen*invDen;
    res = std::sqrt(res);
    res *= Constants::radToDegrees;
    return res;
  }

  float thetaVimpactErr(float& y, float& dy, float& x, float& dx){
    float res = 0.0;
    if(x == 0.0) return res;
    double invDen = 1.0/(x*x + y*y);

    res += dx*dx * y*y * invDen*invDen;
    res += dy*dy * x*x * invDen*invDen;
    res = std::sqrt(res);
    res *= Constants::radToDegrees;
    return res;
  }

  static int ncalls(0);

  class ReadExtrapol : public art::EDAnalyzer {
  public:
    explicit ReadExtrapol(fhicl::ParameterSet const& pset):
      art::EDAnalyzer(pset),

      _fitterModuleLabel(pset.get<string>("fitterModuleLabel")),
      _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
      _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
      _diagLevel(pset.get<int>("diagLevel",0)),
      _qualityCuts(pset.get<int>("qualityCuts",2)),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _makerModuleLabel(pset.get<std::string>("makerModuleLabel", "makeSH")),
      _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
      _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
      _elextractModuleLabel(pset.get<std::string>("elextractModuleLabel", "extractElData")),
      _extractElectronsData(pset.get<string>("elextractModuleLabel")),
      _trkToCaloExtrapolModuleLabel(pset.get<std::string>("trkToCaloExtrapolModuleLabel", "TrkExtrapol")),
      _Ntup(0),
      _application(nullptr),
      _directory(0){
      // construct the data product instance name
      _iname = _fdir.name() + _tpart.name();
    }

    virtual ~ReadExtrapol() {
    }
    void beginJob();
    void endJob() {}

    void analyze(art::Event const& e );

  private:

    void doExtrapolation(art::Event const& evt, bool skip);
    // Module label of the module that performed the fits.
    std::string _fitterModuleLabel;
    TrkParticle _tpart;
        
    TrkFitDirection _fdir;
        
    std::string _iname;
    // Diagnostic level
    int _diagLevel;

    //Select trks with Dave quality cuts
    int _qualityCuts;

    // Label of the generator.
    std::string _generatorModuleLabel;

    // Label of the G4 module
    std::string _g4ModuleLabel;

    //Label of makeSH
    std::string _makerModuleLabel;

    // Label of the calo readout hits maker
    std::string _caloReadoutModuleLabel;

    // Label of the calo crystal hists maker
    std::string _caloCrystalModuleLabel;

    string _elextractModuleLabel;
    string _extractElectronsData;

    // Label of the extrapolated impact points
    std::string _trkToCaloExtrapolModuleLabel;

    bool _skipEvent;

    TTree* _Ntup;//Ntupla which contains informations about the extrapolation starting from MC

    // The job needs exactly one instance of TApplication.  See note 1.
    unique_ptr<TApplication> _application;

    // Save directory from beginJob so that we can go there in endJob. See note 3.
    TDirectory* _directory;


    Int_t _match,//MC and Reco trajectories are compatible in a given tolerance
      _seedQC,//the electron generated satisfied the quality cuts
      _trkIdFound,//one reco trajectory with same trkId of the MC was found
      _vaneFound, //one reco trajectory with same trkId and intersected vane of the MC was found
      _caloInt,//flag used to distinguish the reconstructed trajectories which have at least one intersection with the calorimeter volume
      _recoIndex,//index of the intersection point; 0 is the first intersection point of the trj with the calo, 1 is the seond, and so on.
      _recoTrkId,
      _isMiss;//I use this flag to calculate how many trajectories don't have an intersection with the calorimeter; is 0 if doesn't intersect, 1 if the trj intersects the calorimeter


    Size_t _pointsInVaneIndex;//index of the MC point which belongs the intersection of a trajectory with a single vane
    //runs from 0, ... to (number of intersections of a trajectory with a single vane - 1)

    Float_t _evt,//event Id
      _dist,//distance between MC and reco point, by default is 0.0
      _seedPx ,//MC impact particle coordinates in the Mu2e general frame [mm]
      _seedPy ,
      _seedPz ,
      _seedPpu ,//MC impact momentum components in the local vane frame [MeV]
      _seedPpv ,
      _seedPpw ,
      _seedE ,//MC impact kinetic energy [MeV]
      _seedTime,//MC impact impact time [ns]
      _seedPpx ,//MC impact momentum components in the Mu2e general frame [MeV]
      _seedPpy ,
      _seedPpz ,
      _seedPu ,//MC impact particle coordinates in the local vane frame [mm]
      _seedPv ,
      _seedPw ,
      _seedThetaW,//MC impact impact angle projected to the plane identified by the vectors U, W (local vane frame) [deg]
      _seedThetaV,//MC impact angle projected to the plane identified by the vectors U, V (local vane frame) [deg]
      _recoPx ,//reconstructed impact particle coordinates in the Mu2e general frame [mm]
      _recoPy ,
      _recoPz ,
      _recoPpu ,//reconstructed impact particle momentum in the local vane frame [MeV]
      _recoPpuErr ,
      _recoPpv ,
      _recoPpvErr ,
      _recoPpw ,
      _recoPpwErr ,
      _recoE ,//reconstructed impact kinetic energy [MeV]
      _recoEErr,
      _recoTime,//reconstructed impact impact time [ns]
      _recoTimeErr,
      _recoPpx ,//reconstructed impact momentum components in the Mu2e general frame [MeV]
      _recoPpy ,
      _recoPpz ,
      _recoPu ,//reconstructed impact particle coordinates in the local vane frame [mm]
      _recoPv ,
      _recoPw ,
      _recoPuErr ,//reconstructed impact particle coordinates in the local vane frame [mm]
      _recoPvErr ,
      _recoPwErr ,
      _recoThetaW,//reconstructed impact impact angle projected to the plane identified by the vectors U, W (local vane frame) [deg]
      _recoThetaWErr,
      _recoThetaV,//reconstructed impact impact angle projected to the plane identified by the vectors U, V (local vane frame) [deg]
      _recoThetaVErr,
      _recoT0,
      _recoT0Err,
      _recoMomentumT0,
      _recoMomentumT0Err,
      _recoPathLenght,
      _recoPathLenghtErr;


  };



  bool FindTrkId(std::vector<unsigned int> vec, unsigned int t){
    bool res = false;

    unsigned int size = vec.size();
    if(size!=0){
      unsigned int cont = 0;
      while(/*!res ||*/ cont!=size){
	if(vec[cont] == t) {
	  res = true;
	}
	++cont;

      }
    }
    return res;
  }

  size_t MaxKey(std::map<size_t, unsigned int>& map){
    unsigned int max = 0;
    size_t key = map.begin()->first;
    max = map.begin()->second;
    for(std::map<size_t, unsigned int>::iterator it = map.begin(); it!= map.end(); ++it  ){
      if(it->second > max){
	max = it->second;
	key = it->first;
      }
    }
    return key;
  }

  double distVectors(CLHEP::Hep3Vector vec1, CLHEP::Hep3Vector vec2){
    double res = 0.0;

    res = cet::sum_of_squares(( vec1.y() - vec2.y()), (vec1.z() - vec2.z() ), ( vec1.z() - vec2.z()) );
    res = std::sqrt(res);
    return res;
  }


  CLHEP::Hep3Vector fromTrkToMu2eGeneralFrame(CLHEP::Hep3Vector vec){
    art::ServiceHandle<GeometryService> geom;
    double solenoidOffSetX = geom->config().getDouble("mu2e.solenoidOffset");
    double solenoidOffSetZ = -geom->config().getDouble("mu2e.detectorSystemZ0");

    Hep3Vector posMu2eFrame = vec;
    posMu2eFrame.setX(posMu2eFrame.x() - solenoidOffSetX);
    posMu2eFrame.setZ(posMu2eFrame.z() - solenoidOffSetZ);
    return posMu2eFrame;
  }

  bool findKalRep(art::Handle<TrkToCaloExtrapolCollection> &trkToCaloCollection, KalRep const* trk){
    bool res = false;
    if(trkToCaloCollection->size() == 0) return res;
    size_t i=0;
    while(!res && i<trkToCaloCollection->size()){
      KalRepPtr const& trkPtr = trkToCaloCollection->at(i).trk();
      const KalRep*  const &tmpTrk = *trkPtr;

      const TrkId &trkId1 = trk->id();
      const TrkId &trkId2 = tmpTrk->id();

      //if( &trep1 == &trep2 ){
      if( trkId1.operator long int() == trkId2.operator long int() ){
	res = true;
      }else{
	++i;
      }
    }

    return res;
  }



  void ReadExtrapol::beginJob( ) {

    cout << "start ReadExtrapol..."<<endl;

  }



  void ReadExtrapol::analyze(art::Event const& evt ) {

    ++ncalls;

    if (ncalls == 1) {

      art::ServiceHandle<art::TFileService> tfs;
      _Ntup        = tfs->make<TTree>("MCtrk", "MCtrkExtrapol trajectory info");

      _Ntup->Branch("evt", &_evt , "evt/F");
      _Ntup->Branch("match",     &_match , "match/I");
      _Ntup->Branch("isMiss",     &_isMiss , "isMiss/I");
      _Ntup->Branch("trkIdFound",     &_trkIdFound , "trkIdFound/I");
      _Ntup->Branch("vaneFound",     &_vaneFound , "vaneFound/I");
      _Ntup->Branch("pointsInVaneIndex",     &_pointsInVaneIndex , "pointsInVaneIndex/S");
      _Ntup->Branch("seedQC",     &_seedQC , "seedQC/I");
      _Ntup->Branch("caloInt",    &_caloInt , "caloInt/I");
      _Ntup->Branch("dist",       &_dist , "dist/F");
      _Ntup->Branch("seedPx",     &_seedPx , "seedPx/F");
      _Ntup->Branch("seedPy",     &_seedPy , "seedPy/F");
      _Ntup->Branch("seedPz",     &_seedPz , "seedPz/F");
      _Ntup->Branch("seedE",      &_seedE , "seedE/F");
      _Ntup->Branch("seedTime",   &_seedTime , "seedTime/F");
      _Ntup->Branch("seedPpx",    &_seedPpx , "seedPpx/F");
      _Ntup->Branch("seedPpy",    &_seedPpy , "seedPpy/F");
      _Ntup->Branch("seedPpz",    &_seedPpz , "seedPpz/F");
      _Ntup->Branch("seedPpu",    &_seedPpu , "seedPpu/F");
      _Ntup->Branch("seedPpv",    &_seedPpv , "seedPpv/F");
      _Ntup->Branch("seedPpw",    &_seedPpw , "seedPpw/F");
      _Ntup->Branch("seedPu",     &_seedPu , "seedPu/F");
      _Ntup->Branch("seedPv",     &_seedPv , "seedPv/F");
      _Ntup->Branch("seedPw",     &_seedPw , "seedPw/F");
      _Ntup->Branch("seedThetaW", &_seedThetaW , "seedThetaW/F");
      _Ntup->Branch("seedThetaV", &_seedThetaV , "seedThetaV/F");

      _Ntup->Branch("recoTrkId",  &_recoTrkId , "recoTrkId/I");
      _Ntup->Branch("recoPx",     &_recoPx , "recoPx/F");
      _Ntup->Branch("recoPy",     &_recoPy , "recoPy/F");
      _Ntup->Branch("recoPz",     &_recoPz , "recoPz/F");
      _Ntup->Branch("recoE",      &_recoE , "recoE/F");
      _Ntup->Branch("recoEErr",   &_recoEErr , "recoEErr/F");
      _Ntup->Branch("recoTime",   &_recoTime , "recoTime/F");
      _Ntup->Branch("recoTimeErr",   &_recoTimeErr , "recoTimeErr/F");
      _Ntup->Branch("recoPpx",    &_recoPpx , "recoPpx/F");
      _Ntup->Branch("recoPpy",    &_recoPpy , "recoPpy/F");
      _Ntup->Branch("recoPpz",    &_recoPpz , "recoPpz/F");
      _Ntup->Branch("recoPpu",    &_recoPpu , "recoPpu/F");
      _Ntup->Branch("recoPpv",    &_recoPpv , "recoPpv/F");
      _Ntup->Branch("recoPpw",    &_recoPpw , "recoPpw/F");
      _Ntup->Branch("recoPu",     &_recoPu , "recoPu/F");
      _Ntup->Branch("recoPv",     &_recoPv , "recoPv/F");
      _Ntup->Branch("recoPw",     &_recoPw , "recoPw/F");
      _Ntup->Branch("recoPuErr",  &_recoPuErr , "recoPuErr/F");
      _Ntup->Branch("recoPvErr",  &_recoPvErr , "recoPvErr/F");
      _Ntup->Branch("recoPwErr",  &_recoPwErr , "recoPwErr/F");
      _Ntup->Branch("recoThetaW", &_recoThetaW , "recoThetaW/F");
      _Ntup->Branch("recoThetaV", &_recoThetaV , "recoThetaV/F");
      _Ntup->Branch("recoThetaWErr", &_recoThetaWErr , "recoThetaWErr/F");
      _Ntup->Branch("recoThetaVErr", &_recoThetaVErr , "recoThetaVErr/F");
      _Ntup->Branch("recoIndex",  &_recoIndex , "recoIndex/I");//,
      _Ntup->Branch("recoT0", &_recoT0 , "recoT0/F");
      _Ntup->Branch("recoT0Err", &_recoT0Err , "recoT0Err/F");
      _Ntup->Branch("recoMomentumT0", &_recoMomentumT0 , "recoMomentumT0/F");
      _Ntup->Branch("recoMomentumT0Err", &_recoMomentumT0Err , "recoMomentumT0Err/F");
      _Ntup->Branch("recoPathLenght", &_recoPathLenght , "recoPathLenght/F");//_recoPathLenghtErr
      _Ntup->Branch("recoPathLenghtErr", &_recoPathLenghtErr , "recoPathLenghtErr/F");//_recoPathLenghtErr
    }


    doExtrapolation(evt, _skipEvent);

  } // end of analyze


  void ReadExtrapol::doExtrapolation(art::Event const& evt, bool skip){

    //Get handle to calorimeter
    art::ServiceHandle<GeometryService> geom;
    if(! geom->hasElement<VaneCalorimeter>() ) return;
    GeomHandle<VaneCalorimeter> cg;

    // Get handles to calorimeter collections
    art::Handle<CaloHitCollection> caloHits;
    evt.getByLabel(_caloReadoutModuleLabel, caloHits);

    art::Handle<CaloCrystalHitCollection>  caloCrystalHits;
    evt.getByLabel(_caloCrystalModuleLabel, caloCrystalHits);

    // Get the persistent data about pointers to StepPointMCs
    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    evt.getByLabel(_caloReadoutModuleLabel,"CaloHitMCCrystalPtr",mcptrHandle);
    PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();


    art::Handle<PtrStepPointMCVectorCollection> mcptrHandleStraw;
    evt.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandleStraw);
    PtrStepPointMCVectorCollection const* hits_mcptrStraw = mcptrHandleStraw.product();

    if (!( caloHits.isValid())) {
      return;
    }

    if (!caloCrystalHits.isValid()) {
      cout << "NO CaloCrystalHits" << endl;
      return;
    }

    // Get handles to the generated and simulated particles.
    art::Handle<GenParticleCollection> genParticles;
    evt.getByLabel(_generatorModuleLabel, genParticles);

    art::Handle<SimParticleCollection> simParticles;
    evt.getByLabel(_g4ModuleLabel, simParticles);

    art::Handle<TrkToCaloExtrapolCollection>  trjExtrapols;
    evt.getByLabel(_trkToCaloExtrapolModuleLabel, trjExtrapols);


    art::Handle<VisibleGenElTrackCollection> genEltrksHandle;
    evt.getByLabel(_extractElectronsData,genEltrksHandle);

    VisibleGenElTrackCollection const* genEltrks = genEltrksHandle.product();
    std::vector<mu2e::VisibleGenElTrack>::const_iterator genEltrk_it;

    art::Handle<StrawHitCollection> strawHandle;
    evt.getByLabel(_makerModuleLabel,strawHandle);
    StrawHitCollection const& strawHits = *strawHandle;


    if (!( hits_mcptrStraw->size() == strawHits.size() ) ) {
      mf::LogError("RANGE")
	<< "Strawhits: " << strawHits.size()
	<< " MCPtr: " << hits_mcptrStraw->size();
      return;
    }

    trkIdVector QCvec;

    if(_qualityCuts > 1){
      double trkMomCut = 100.0;//MeV
      int NtrkCut =0;
      int NtrkTot = 0;

      //mapping &counting the electrons with quality cuts in the TRK
      for ( genEltrk_it = genEltrks->begin(); genEltrk_it!= genEltrks->end(); ++genEltrk_it ){
	++NtrkTot;

	VisibleGenElTrack &iEltrk = const_cast<VisibleGenElTrack &>(*genEltrk_it);
	GenElHitData& hdil = iEltrk.getithLoopHit(0);
	GenElHitData& ldil = iEltrk.getHit((int)(iEltrk.getNumOfHit() - 1) );



	double cosTheta = iEltrk.getTrkLrntzVec().cosTheta() ;
	double cosPitch = hdil._hitMomentum.cosTheta();
	double lcosPitch = ldil._hitMomentum.cosTheta();

	if(_diagLevel>2){
	  cout << "cosThetafirst = "<<cosPitch <<", "<<"costhetaLast = "<<lcosPitch<<endl;
	  cout<< "ldil._hitMomentum.mag() = "<<ldil._hitMomentum.mag()<<"hdil._hitMomentum.mag() = "<< hdil._hitMomentum.mag()<<endl;
	}

	bool condition = true;
	//the following condition are the same used by Dave Brown for TTracker studies
	condition &= ( iEltrk.getNumOfHit() >= 20 );
	condition &= ( hdil._hitMomentum.mag() >= trkMomCut );
	condition &= ( cosTheta >= -0.5 );
	condition &= ( cosTheta <=  0.5 );
	condition &= ( cosPitch > 0.5 );
	condition &= ( cosPitch < 0.70710678118655 );// 1 / sqrt(2)

	if( condition ){
	  NtrkCut++;// # of electrons at the entrance of the TRK which have trkMomCut MeV of kinetic energy and which will light at least 20 straws
	  if( !QCvec.find(iEltrk.getTrkID().asUint()) ){
	    if(_diagLevel>2){
	      cout<<"!QCvec.find(iEltrk.getTrkID().asUint())"<<
		", trkId = "<< iEltrk.getTrkID().asUint()<<endl;
	    }
	    QCvec.push_back(iEltrk.getTrkID().asUint());
	  }


	}//end if(condition)
      }//end loop TRK mapping
      if (NtrkTot==0) return;
    }

    ElecMap recoMap, caloMap, extrMap;

    int iVane = -1;
    if (!( strawHandle.isValid())) {
      return;
    }


    //        StrawHit const* it0 = &(*strawHits.begin());
    //StrawHit const* it0 = &strawHits.front();
    // fill the map of the reconstructed trajectories
    int count = 0;

    for(size_t i=0; i<trjExtrapols->size(); ++i){
      elecData tmpElec;
      KalRepPtr const& trkPtr = trjExtrapols->at(i).trk();
      const KalRep *  const &trk = *trkPtr;
      if(i>0){
	KalRepPtr const& tmpTrkPtr = trjExtrapols->at(i-1).trk();
	const KalRep *  const &tmpTrk = *tmpTrkPtr;
	if(trk == tmpTrk){
	  ++count;
	}else{
	  count = 0;
	}
      }
      TrkHotList const* hots = /*const_cast<TrkHotList*>*/( trk->hotList() );

      //Map of track id as key, and number of occurrences as value
      map<size_t , unsigned int > StrawTracksMap;

      for(TrkHotList::hot_iterator ihot=hots->begin();ihot != hots->end();++ihot){

	TrkStrawHit const* hit = dynamic_cast<TrkStrawHit const*>(ihot.get());

	if(_diagLevel>3){
	  hit->printAll(cout);
	  cout<<"strawHits.size() : "<<strawHits.size()<<endl;
	  cout<<"hits_mcptrStraw->size() : "<<hits_mcptrStraw->size()<<endl;
	  cout<<"hit->index() = "<< hit->index()<<endl;
	}
	size_t index = hit->index();//&(hit->strawHit() ) - it0;
	if(index >= hits_mcptrStraw->size() ) continue;
	PtrStepPointMCVector const& mcptr(hits_mcptrStraw->at(index ) );

	if(mcptr.size()==0) continue;

	for (size_t j = 0; j < 1/*mcptr.size()*/; ++j) {//I use only the first stepPoint inside the hit because
	  //it generated the drift distance that it is fitted
	  //When the hit maker algorithm will change it will be usefull
	  //to loop on all the step points. (FIXME)
	  StepPointMC const& mchit = *mcptr[j];
	  art::Ptr<SimParticle> const& simptr = mchit.simParticle();
	  //art::Ptr<SimParticle>::key_type simKey(simptr.key());
	  SimParticle const& sim = *simptr;


	  // The simulated particle that made this hit.
	  SimParticleCollection::key_type trackId(mchit.trackId());


	  if(sim.isPrimary()){

	    //Find in the map if the track is already stored
	    map<size_t , unsigned int >::iterator it;
	    it = StrawTracksMap.find(trackId.asUint());
	    unsigned int temp = 0;
	    //if the contributing track id does not exist in the map
	    //add an element to the map itself, energy to the list and pdgId and genId to the vectors
	    if (it!=StrawTracksMap.end()) {
	      temp = StrawTracksMap[it->first];

	      StrawTracksMap[it->first] = temp + 1;
	    }else {
	      //insert track id in the trackId vector
	      StrawTracksMap[trackId.asUint()] = 1;
	    }
	    if(_diagLevel>3){

	      cout<<"trackId() = "<< it->first<<
		", tmp value = "<<temp<<
		", #occorrenze = "<< StrawTracksMap[it->first]<<endl;
	    }
	    break;
	  }

	}
      }//end loop on the straws

      if(StrawTracksMap.size()==0) continue;
      size_t trkId = MaxKey(StrawTracksMap);
      //if(StrawTracksMap[trkId] < 1 ) continue;
      if(_diagLevel>2){
	cout<<"filling tmpElec..."<<endl;
      }
      tmpElec._vaneId        = trjExtrapols->at(i).vaneId();
      tmpElec._index         = count;
      tmpElec._impTime       = trjExtrapols->at(i).time();
      tmpElec._impTimeErr    = trjExtrapols->at(i).timeErr();
      tmpElec._impMom3Vec    = trjExtrapols->at(i).momentum();
      tmpElec._impMom3VecErr = trjExtrapols->at(i).momentumErr();
      tmpElec._impPos        = Hep3Vector(trjExtrapols->at(i).entrancePosition().x(), trjExtrapols->at(i).entrancePosition().y(), trjExtrapols->at(i).entrancePosition().z());
      tmpElec._impPosErr     = trjExtrapols->at(i).entrancePositionErr();
      tmpElec._t0            = trjExtrapols->at(i).t0();
      tmpElec._t0Err         = trjExtrapols->at(i).t0Err();
      tmpElec._t0Momentum    = trjExtrapols->at(i).t0Momentum();
      tmpElec._t0MomentumErr = trjExtrapols->at(i).t0MomentumErr();
      tmpElec._pathLenght    = trjExtrapols->at(i).pathLengthEntrance();
      tmpElec._pathLenghtErr    = trjExtrapols->at(i).pathLenghtEntranceErr();
      recoMap[trkId].push_back(tmpElec);
      if(_diagLevel>2){
	cout<<"recoMap filled..."<<endl;
      }

    }

    if(_diagLevel>2){
      cout<<"end filling recoMap..."<<endl<<
	"start caloCrystalHits loop..."<<endl;
    }

    // fill the map of MC tracks
    double tollCaloHits = cg->crystalHalfTrans()*15.0;
    //    double tollCaloHits = cg->crystalHalfSize()*15.0;

    for(size_t i=0; i<caloCrystalHits->size(); ++i){
      elecData tmpElec;
      bool foundGen =false;
      CaloCrystalHit const& hit = (caloCrystalHits->at(i));


      std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
      if(ROIds.size()<1 ) continue;

      CaloHit const& thehit = *ROIds.at(0);
      if(_diagLevel>3){
	cout<<"step "<< i <<" : "<< caloCrystalHits->size()<<endl;
	thehit.print();
	cout<<endl;
      }


      size_t collectionPosition = ROIds.at(0).key();

      iVane = cg->vaneByRO(thehit.id());

      PtrStepPointMCVector const & mcptr(hits_mcptr->at(collectionPosition));
      if(mcptr.size() <= 0) continue;

      size_t nHitsPerCrystal = mcptr.size();
      unsigned int trkId = 0;

      for (size_t j2=0; j2<nHitsPerCrystal; ++j2) {

	StepPointMC const& mchit = *mcptr[j2];

	art::Ptr<SimParticle> const& simptr = mchit.simParticle();


	// The simulated particle that made this hit.
	SimParticleCollection::key_type trackId(mchit.trackId());
	SimParticle const& sim = *simptr;

	if(sim.fromGenerator() ){
	  if(tmpElec._impTime > mchit.time() ){
	    foundGen               = true;
	    trkId                  = trackId.asUint();
	    tmpElec._vaneId        = iVane;
	    tmpElec._impTime       = mchit.time();
	    tmpElec._impMom3Vec    = mchit.momentum();
	    tmpElec._impPos        = mchit.position();

	  }//end  if(caloMap[trackId.asUint()][iVane]._impTime > mchit.time() ){

	}else{//end if(sim.fromGenerator)
	  if(_diagLevel>2){
	    cout<<"crystal hit NOT from generator..."<<endl;
	  }
	}

      }//end loop on nHitsPerCrystal
      if(foundGen){
	if(_diagLevel>2){
	  cout<<"foundGen!!"<<endl;
	}

	if(caloMap.find(trkId) != caloMap.end()){
	  if(_diagLevel>3){
	    cout<<"trkId = "<<trkId<<
	      ", caloMap.find(trkId)!!"<<endl;
	  }
	  size_t index=0;
	  if(!caloMap[trkId].findElecData(tmpElec, index) ){
	    if(_diagLevel>3){
	      cout<<"!caloMap[trkId].findElecData(tmpElec, index)..."<<
		", index = "<<index<<
		", _vec.size() = "<< caloMap[trkId].size()<<endl;
	    }
	    if(caloMap[trkId].size() > 0){
	      bool trovato = false, condition1 = false, condition2 = false;
	      size_t c = 0;
	      double tmpDist = 0.0;
	      while(!trovato && (c != caloMap[trkId].size()) ){
		tmpDist = distVectors(tmpElec._impPos, caloMap[trkId].impPos(c));

		//condition to update the information a stored seed particle
		condition1 = (caloMap[trkId]._vec.at(c) > tmpElec);
		condition1 &= (tmpDist < tollCaloHits);

		//condition for filling a new element which represents a new intersection point of a trajectory
		condition2 = (caloMap[trkId]._vec.at(c) < tmpElec);
		condition2 &= (tmpDist > tollCaloHits);

		if(condition1){
		  trovato = true;
		  if(_diagLevel>2){
		    cout<<" step = "<<c<<
		      ", caloMap[trkId]._vec.at(c) > tmpElec"<<endl<<
		      ", distance = "<< tmpDist<<endl;
		    tmpElec.print(cout);
		    caloMap[trkId]._vec.at(c).print(cout);

		  }
		}else if(condition2){
		  trovato = true;
		  if(_diagLevel>2){
		    cout<<" step = "<<c<<
		      ", caloMap[trkId]._vec.at(c) < tmpElec"<<endl<<
		      ", distance = "<< tmpDist<<endl;
		    tmpElec.print(cout);
		    caloMap[trkId]._vec.at(c).print(cout);
		  }
		}else{
		  ++c;
		}
	      }
	      if(trovato){
		if(condition1){
		  caloMap[trkId].remove(c);
		  caloMap[trkId].push_back(tmpElec);
		}else if(condition2){
		  caloMap[trkId].push_back(tmpElec);
		}
	      }
	    }
	  }
	}else{
	  caloMap[trkId].push_back(tmpElec);
	}
      }
    }//end loop on caloCrystalHits

    //-----------------------------------------------
    double solenoidOffSetX = geom->config().getDouble("mu2e.solenoidOffset");//3904.;//[mm]
    double solenoidOffSetZ = -geom->config().getDouble("mu2e.detectorSystemZ0");//-10200.;

    double thetaWimpact = 0.0, thetaVimpact = 0.0;

    if(_diagLevel>2){
      cout<<"caloMap size = "<< caloMap.size()<<endl;

      for(ElecMap::iterator it = caloMap.begin(); it != caloMap.end(); ++it){
	cout<<"caloMap, trkId() = "<<it->first<<endl;
	for(size_t i=0; i < it->second.size(); ++i){
	  cout<<"caloMap print : "<<endl;
	  it->second._vec.at(i).print(cout);
	}
      }

      cout<<"//-----------------------------------------------//"<<endl<<
	"recoMap size = "<< recoMap.size()<<endl;
      for(ElecMap::iterator it = recoMap.begin(); it != recoMap.end(); ++it){
	cout<<"recoMap, trkId() = "<<it->first<<endl;
	for(size_t i=0; i < it->second.size(); ++i){
	  cout<<"recoMap print : "<<endl;
	  it->second._vec.at(i).print(cout);
	}
      }
    }

    unsigned int tmpVane = 0, tmpTrkId = 0;
    double toll = 4.0*cg->crystalHalfTrans();
    // double toll = 4.0*cg->crystalHalfSize();
    if(_diagLevel>2){
      cout<<"caloMap.size() = "<<caloMap.size()<<endl;
    }

    double tmpEnergyErr = 0.0, tmpErrPos = 0.0;
    Hep3Vector tmpPosVaneFrameErr, tmpPerr;

    for(ElecMap::iterator it = caloMap.begin(); it != caloMap.end(); ++it){

      _evt = evt.id().event();
      _match = 0;
      _isMiss = 0;
      _dist  = 0.0;
      _trkIdFound = 0;
      _vaneFound = 0;
      _pointsInVaneIndex = 0;
      _caloInt = 0;
      _seedPx = 0.0;
      _seedPy = 0.0;
      _seedPz = 0.0;
      _seedPpu = 0.0;
      _seedPpv = 0.0;
      _seedPpw = 0.0;
      _seedE = 0.0;
      _seedTime= 0.0;
      _seedPpx = 0.0;
      _seedPpy = 0.0;
      _seedPpz = 0.0;
      _seedPu = 0.0;
      _seedPv = 0.0;
      _seedPw = 0.0;
      _seedThetaW= 0.0;
      _seedThetaV= 0.0;
      _seedQC   = 0;
      _recoPx = 0.0;
      _recoPy = 0.0;
      _recoPz = 0.0;
      _recoPpu = 0.0;
      _recoPpuErr = 0.0;
      _recoPpv = 0.0;
      _recoPpvErr = 0.0;
      _recoPpw = 0.0;
      _recoPpwErr = 0.0;
      _recoE = 0.0;
      _recoEErr = 0.0;
      _recoTime= 0.0;
      _recoTimeErr = 0.0;
      _recoPpx = 0.0;
      _recoPpy = 0.0;
      _recoPpz = 0.0;
      _recoPu = 0.0;
      _recoPv = 0.0;
      _recoPw = 0.0;
      _recoPuErr = 0.0;
      _recoPvErr = 0.0;
      _recoPwErr = 0.0;
      _recoThetaW= 0.0;
      _recoThetaV = 0.0;
      _recoThetaWErr= 0.0;
      _recoThetaVErr = 0.0;
      _recoIndex = 0;
      _recoTrkId = 0;
      _recoT0 = 0.0;
      _recoT0Err= 0.0;
      _recoMomentumT0= 0.0;
      _recoMomentumT0Err= 0.0;
      _recoPathLenght= 0.0;
      _recoPathLenghtErr= 0.0;
      if(_diagLevel>2){
	cout<<"caloTrk.Id() = "<<it->first<<endl;
      }
      tmpTrkId = it->first;
      if( recoMap.find(tmpTrkId) != recoMap.end() ){//search the tmpTrkId in the RecoMap
	_trkIdFound = 1;
	_caloInt = 1;
	_recoTrkId = tmpTrkId;
	trkIdVector MCVaneList = caloMap[tmpTrkId].vaneList();//vector of the intersected vane indexes from MC
	trkIdVector RecoVaneList = recoMap[tmpTrkId].vaneList();//vector of the intersected vane indexes from extrapolation

	if(MCVaneList.size() == 0) continue;//MC doesn't show any intersected vane for tmpTrkId, it is a redundant check

	if(_diagLevel>2){
	  cout<<"MCVaneList.size() = "<<MCVaneList.size()<<endl;
	}

	for(size_t i=0; i< MCVaneList.size(); ++i){

	  tmpVane = MCVaneList.at(i);

	  if(_diagLevel>2){
	    cout<<"MCVaneList.at("<<i<<")"<<
	      "tmpVane = "<<tmpVane<<endl;
	  }

	  if( !RecoVaneList.find(tmpVane) ) {//search of tmpVane in the RecoMap[tmpTrkId](FIXME)
	    _Ntup->Fill();
	    continue;
	  }

	  _vaneFound = 1;
	  //vector of the index of caloMap[tmpTrkId]._vec which belong to elements with vaneId = tmpVane
	  trkIdVector CaloPointInVane = caloMap[tmpTrkId].findVane( tmpVane );

	  //vector of the index of recoMap[tmpTrkId]._vec which belong to elements with vaneId = tmpVane
	  trkIdVector RecoPointInVane = recoMap[tmpTrkId].findVane( tmpVane );


	  //DistVector distVec;
	  trkIdVector indexesUsed;
	  double tmpDist = 0.0;
	  size_t tmpIndex = 0;

	  for(size_t j=0; j<CaloPointInVane.size(); ++j){
	    if(_diagLevel>2){
	      cout<<"CaloPointInVane at "<<j<<
		", RecoPointInVane.size() = "<< RecoPointInVane.size()<<endl;
	    }
	    DistVector distVec;
	    if(RecoPointInVane.size() == 0) break;
	    for(size_t i=0; i<RecoPointInVane.size(); ++i){
	      if(_diagLevel>2){
		cout<<"RecoPointInVane at "<<i<<endl;
	      }
	      CLHEP::Hep3Vector tmpV = fromTrkToMu2eGeneralFrame(recoMap[tmpTrkId].impPos(i) );
	      tmpDist = distVectors(caloMap[tmpTrkId].impPos(j), tmpV );
	      if(_diagLevel>2){
		cout<<"tmpDist = "<<tmpDist<<endl;
	      }
	      distVec.push_back(tmpDist);
	    }
	    tmpIndex = distVec.indexMinCont();



	    _dist = distVec.at(tmpIndex);
	    _pointsInVaneIndex = j;

	    if(_diagLevel>2){
	      cout<<"found matching!"<<endl;
	    }
	    if(QCvec.find(tmpTrkId) ){
	      _seedQC = 1;
	    }else{
	      _seedQC = 0;
	    }
	    Vane const v1 = cg->vane(tmpVane );
	    Hep3Vector momentumRotUnit = (v1.rotation())*(it->second.impMom3Vec( CaloPointInVane.at(j) ).unit());

	    thetaWimpact = std::atan(-1.0*momentumRotUnit.getZ() / momentumRotUnit.getX() ) ;
	    _seedThetaW = thetaWimpact*Constants::radToDegrees;

	    thetaVimpact = std::atan(momentumRotUnit.getY() /  momentumRotUnit.getX() );
	    _seedThetaV = thetaVimpact*Constants::radToDegrees;

	    _seedE  = it->second.impMom3Vec(CaloPointInVane.at(j)).mag();
	    _seedTime = it->second.impTime(CaloPointInVane.at(j));
	    _seedPx = it->second.impPos(CaloPointInVane.at(j)).x();
	    _seedPy = it->second.impPos(CaloPointInVane.at(j)).y();
	    _seedPz = it->second.impPos(CaloPointInVane.at(j)).z();
	    Hep3Vector vaneFrame = cg->toVaneFrame(it->second.vaneId(CaloPointInVane.at(j)), it->second.impPos(CaloPointInVane.at(j)) );
	    _seedPu = vaneFrame.x();
	    _seedPv = vaneFrame.y();
	    _seedPw = vaneFrame.z();
	    _seedPpx = it->second.impMom3Vec(CaloPointInVane.at(j) ).x();
	    _seedPpy = it->second.impMom3Vec(CaloPointInVane.at(j) ).y();
	    _seedPpz = it->second.impMom3Vec(CaloPointInVane.at(j) ).z();
	    if(_diagLevel>2){
	      cout<<"momentum : "<<it->second.impMom3Vec(CaloPointInVane.at(j) )<<endl;
	    }
	    _seedPpu = it->second.impMom3Vec(CaloPointInVane.at(j) ).mag()*momentumRotUnit.x();
	    _seedPpv = it->second.impMom3Vec(CaloPointInVane.at(j) ).mag()*momentumRotUnit.y();
	    _seedPpw = it->second.impMom3Vec(CaloPointInVane.at(j) ).mag()*momentumRotUnit.z();

	    Hep3Vector momdir;
	    HepVector momvec(3);
	    momdir = recoMap[tmpTrkId].impMom3Vec(RecoPointInVane.at(tmpIndex)).unit();
	    for(int icor=0;icor<3;icor++){
	      momvec[icor] = momdir[icor];
	    }

	    tmpEnergyErr = sqrt(recoMap[tmpTrkId].impMom3VecErr(RecoPointInVane.at(tmpIndex)).covMatrix().similarity(momvec));


	    momentumRotUnit = (v1.rotation())*(recoMap[tmpTrkId].impMom3Vec(RecoPointInVane.at(tmpIndex)).unit());
	    thetaWimpact = std::atan(-1.0*momentumRotUnit.getZ() / momentumRotUnit.getX() ) ;
	    _recoThetaW = thetaWimpact*Constants::radToDegrees;

	    thetaVimpact = std::atan(momentumRotUnit.getY() /  momentumRotUnit.getX() );
	    _recoThetaV = thetaVimpact*Constants::radToDegrees;

	    _recoE = recoMap[tmpTrkId].impMom3Vec(RecoPointInVane.at(tmpIndex)).mag();
	    _recoEErr = tmpEnergyErr;

	    _recoT0 = recoMap[tmpTrkId].t0(RecoPointInVane.at(tmpIndex));
	    _recoT0Err= recoMap[tmpTrkId].t0Err(RecoPointInVane.at(tmpIndex));
	    _recoMomentumT0= recoMap[tmpTrkId].t0Momentum(RecoPointInVane.at(tmpIndex)).mag();

	    tmpEnergyErr = sqrt(recoMap[tmpTrkId].t0MomentumErr(RecoPointInVane.at(tmpIndex)).covMatrix().similarity(momvec));
	    _recoMomentumT0Err= tmpEnergyErr;
	    _recoPathLenght= recoMap[tmpTrkId].pathLenght(RecoPointInVane.at(tmpIndex));
	    _recoPathLenghtErr= recoMap[tmpTrkId].pathLenghtErr(RecoPointInVane.at(tmpIndex));

	    _recoTime = recoMap[tmpTrkId].impTime(RecoPointInVane.at(tmpIndex));
	    _recoTimeErr = recoMap[tmpTrkId].impTimeErr(RecoPointInVane.at(tmpIndex));

	    Hep3Vector posMu2eFrame(recoMap[tmpTrkId].impPos(RecoPointInVane.at(tmpIndex)) );
	    posMu2eFrame.setX(posMu2eFrame.x() - solenoidOffSetX);
	    posMu2eFrame.setZ(posMu2eFrame.z() - solenoidOffSetZ);
	    if(_diagLevel>2){
	      cout<<"posMu2eFrame : "<< posMu2eFrame<<endl;
	    }
	    _recoPx = posMu2eFrame.x();
	    _recoPy = posMu2eFrame.y();
	    _recoPz = posMu2eFrame.z();

	    vaneFrame = cg->toVaneFrame(tmpVane, posMu2eFrame);
	    _recoPu = vaneFrame.x();
	    _recoPv = vaneFrame.y();
	    _recoPw = vaneFrame.z();

	    //define the local axes V and W for the local vane frame
	    Hep3Vector Vaxes(0.0, 1.0, 0.0), Waxes(0.0, 0.0, 1.0);

	    //move these axes into the Mu2e general frame
	    Vaxes = (v1.rotation())*(Vaxes);
	    Waxes = (v1.rotation())*(Waxes);

	    double scaleErrW = 1.0/fabs( cos(thetaWimpact) );
	    double scaleErrV = 1.0/fabs( cos(thetaVimpact) );
	    momvec[0] = Vaxes.x();
	    momvec[1] = Vaxes.y();
	    momvec[2] = Vaxes.z();

	    tmpErrPos = sqrt(recoMap[tmpTrkId].impPosErr(RecoPointInVane.at(tmpIndex)).covMatrix().similarity(momvec));
	    _recoPvErr = tmpErrPos*scaleErrV;

	    momvec[0] = Waxes.x();
	    momvec[1] = Waxes.y();
	    momvec[2] = Waxes.z();

	    tmpErrPos = sqrt(recoMap[tmpTrkId].impPosErr(RecoPointInVane.at(tmpIndex)).covMatrix().similarity(momvec));
	    _recoPwErr = tmpErrPos*scaleErrW;


	    CLHEP::Hep3Vector tmpPerr = (v1.rotation())*tmpPosVaneFrameErr;//cg->toVaneFrame(tmpVane, tmpPosVaneFrameErr);
	    _recoPuErr = fabs(tmpPerr.x());

	    if(_diagLevel>2){
	      cout<<"Reco momentum : "<< recoMap[tmpTrkId].impMom3Vec(RecoPointInVane.at(tmpIndex)) <<endl;
	    }
	    _recoPpx = recoMap[tmpTrkId].impMom3Vec(RecoPointInVane.at(tmpIndex)).x();
	    _recoPpy = recoMap[tmpTrkId].impMom3Vec(RecoPointInVane.at(tmpIndex)).y();
	    _recoPpz = recoMap[tmpTrkId].impMom3Vec(RecoPointInVane.at(tmpIndex)).z();

	    _recoPpu = recoMap[tmpTrkId].impMom3Vec(RecoPointInVane.at(tmpIndex)).mag()*momentumRotUnit.x();
	    _recoPpv = recoMap[tmpTrkId].impMom3Vec(RecoPointInVane.at(tmpIndex)).mag()*momentumRotUnit.y();
	    _recoPpw = recoMap[tmpTrkId].impMom3Vec(RecoPointInVane.at(tmpIndex)).mag()*momentumRotUnit.z();

	    momvec[0] = 1.0;
	    momvec[1] = 0.0;
	    momvec[2] = 0.0;
	    tmpErrPos = sqrt(recoMap[tmpTrkId].impMom3VecErr(RecoPointInVane.at(tmpIndex)).covMatrix().similarity(momvec));
	    tmpPosVaneFrameErr.setX(tmpErrPos);

	    momvec[0] = 0.0;
	    momvec[1] = 1.0;
	    tmpErrPos = sqrt(recoMap[tmpTrkId].impMom3VecErr(RecoPointInVane.at(tmpIndex)).covMatrix().similarity(momvec));
	    tmpPosVaneFrameErr.setY(tmpErrPos);

	    momvec[2] = 1.0;
	    momvec[1] = 0.0;
	    tmpErrPos = sqrt(recoMap[tmpTrkId].impMom3VecErr(RecoPointInVane.at(tmpIndex)).covMatrix().similarity(momvec));
	    tmpPosVaneFrameErr.setZ(tmpErrPos);

	    tmpPerr = (v1.rotation())*tmpPosVaneFrameErr;
	    _recoPpuErr = tmpPerr.x();
	    _recoPpvErr = tmpPerr.y();
	    _recoPpwErr = tmpPerr.z();

	    _recoThetaWErr = thetaWimpactErr(_recoPpw, _recoPpwErr, _recoPpu, _recoPpuErr);
	    _recoThetaVErr = thetaVimpactErr(_recoPpv, _recoPpvErr, _recoPpu, _recoPpuErr);


	    _recoIndex =  recoMap[tmpTrkId].index(RecoPointInVane.at(tmpIndex));
	    if(distVec.at(tmpIndex) > toll ){
	      _match = 0;
	    }else{
	      _match = 1;
	      RecoPointInVane.remove(tmpIndex);
	    }

	    _Ntup->Fill();

	    if(_diagLevel>2){
	      cout<<"_Ntup filled..."<<endl;
	    }
	  }



	}

	if(_diagLevel>2){
	  cout<<"MCVaneList loop ends..."<<endl;
	}

      } else {//end if(same track Id)
	if(_diagLevel>2){
	  cout<<"not found TrkId..."<<endl;
	}
	_Ntup->Fill();
      }

    }//end loop on caloMap



    //---------------------------------------------------------------------------//
    //   This second part studies how many times the extrapolated trajectories   //
    //   doesn't have any intersection with the calorimeter.                     //
    //---------------------------------------------------------------------------//
    art::Handle<KalRepCollection> trksHandle;
    evt.getByLabel(_fitterModuleLabel,_iname,trksHandle);
    KalRepCollection const& trks = *trksHandle;

    if(_diagLevel>2){
      cout<<"//--------------------------//"<<endl<<
	"//  Start second part...    //"<<endl<<
	"//  Event Number : "<< evt.event()<<"       //"<< endl<<
	"//  trks.size() = "<< trks.size() <<"        //"<<endl<<
	"//--------------------------//"<<endl;
    }

    TrkToCaloExtrapolCollection missExtrapolatedTracks;

    for ( size_t i=0; i< trks.size(); ++i ){

      KalRep const* trep = trks[i];
      if ( !trep ) continue;

      if(!findKalRep(trjExtrapols, trep)){
	int tmpVane = -1;
	HelixTraj trkHel(trep->helix(0).params(),trep->helix(0).covariance());

	double lowrange = trkHel.zFlight(1740), highrange = trkHel.zFlight(3500);

	KalRepPtr tmpRecTrk(trksHandle, i);
	missExtrapolatedTracks.push_back(
					 TrkToCaloExtrapol( tmpVane,i,tmpRecTrk,lowrange, highrange) );

	if(_diagLevel>2){
	  cout<<"Trajectory reconstructed without intersections with the calorimeter volume..."<<endl;
	}

      }
    }


    for(size_t i=0; i<missExtrapolatedTracks.size(); ++i){
      elecData tmpElec;
      KalRepPtr const& trkPtr = missExtrapolatedTracks.at(i).trk();
      const KalRep *  const &trk = *trkPtr;
      TrkHotList const* hots = /*const_cast<TrkHotList*>*/( trk->hotList() );

      //Map of track id as key, and number of occurrences as value
      map<size_t , unsigned int > StrawTracksMap;

      for(TrkHotList::hot_iterator ihot=hots->begin();ihot != hots->end();++ihot){

	TrkStrawHit const* hit = dynamic_cast<TrkStrawHit const*>(ihot.get());

	if(_diagLevel>3){
	  hit->printAll(cout);
	  cout<<"strawHits.size() : "<<strawHits.size()<<endl;
	  cout<<"hits_mcptrStraw->size() : "<<hits_mcptrStraw->size()<<endl;
	  cout<<"hit->index() = "<< hit->index()<<endl;
	}
	size_t index = hit->index();
	if(index >= hits_mcptrStraw->size() ) continue;
	PtrStepPointMCVector const& mcptr(hits_mcptrStraw->at(index ) );

	if(mcptr.size() == 0) continue;
	for (size_t j = 0; j < 1/*mcptr.size()*/; ++j) {//I use only the first stepPoint inside the hit because
	  //it generated the drift distance that it is fitted
	  //When the hit maker algorithm will change it will be usefull
	  //to loop on all the step points. (FIXME)
	  StepPointMC const& mchit = *mcptr[j];
	  art::Ptr<SimParticle> const& simptr = mchit.simParticle();
	  SimParticle const& sim = *simptr;


	  // The simulated particle that made this hit.
	  SimParticleCollection::key_type trackId(mchit.trackId());


	  if(sim.isPrimary()){

	    //Find in the map if the track is already stored
	    map<size_t , unsigned int >::iterator it;
	    it = StrawTracksMap.find(trackId.asUint());
	    unsigned int temp = 0;
	    //if the contributing track id does not exist in the map
	    //add an element to the map itself, energy to the list and pdgId and genId to the vectors
	    if (it!=StrawTracksMap.end()) {
	      temp = StrawTracksMap[it->first];

	      StrawTracksMap[it->first] = temp + 1;
	    }else {
	      //insert track id in the trackId vector
	      StrawTracksMap[trackId.asUint()] = 1;
	    }
	    if(_diagLevel>3){

	      cout<<"trackId() = "<< it->first<<
		", tmp value = "<<temp<<
		", #occorrenze = "<< StrawTracksMap[it->first]<<endl;
	    }
	    break;
	  }

	}
      }//end loop on the straws
      if(_diagLevel>2){
	cout<<"loop in the straws ended..."<<endl;
      }

      if(StrawTracksMap.size()==0) continue;
      size_t trkId = MaxKey(StrawTracksMap);

      tmpElec._vaneId        = missExtrapolatedTracks.at(i).vaneId();
      tmpElec._impTime       = missExtrapolatedTracks.at(i).time();
      tmpElec._impMom3Vec    = missExtrapolatedTracks.at(i).momentum();
      tmpElec._impMom3VecErr = missExtrapolatedTracks.at(i).momentumErr();
      tmpElec._impPos        = Hep3Vector(missExtrapolatedTracks.at(i).entrancePosition().x(), missExtrapolatedTracks.at(i).entrancePosition().y(), missExtrapolatedTracks.at(i).entrancePosition().z());
      tmpElec._impPosErr     = missExtrapolatedTracks.at(i).entrancePositionErr();
      extrMap[trkId].push_back(tmpElec);
      if(_diagLevel>2){
	cout<<"extrMap[trkId].push_back(tmpElec) done..."<<endl;
      }
    }



    for(ElecMap::iterator it = extrMap.begin(); it != extrMap.end(); ++it){
      _evt = evt.id().event();
      _match = 0;
      _isMiss = 0;
      _dist  = 0.0;
      _trkIdFound = 0;
      _vaneFound = 0;
      _pointsInVaneIndex = 0;
      _caloInt = 0;
      _seedPx = 0.0;
      _seedPy = 0.0;
      _seedPz = 0.0;
      _seedPpu = 0.0;
      _seedPpv = 0.0;
      _seedPpw = 0.0;
      _seedE = 0.0;
      _seedTime= 0.0;
      _seedPpx = 0.0;
      _seedPpy = 0.0;
      _seedPpz = 0.0;
      _seedPu = 0.0;
      _seedPv = 0.0;
      _seedPw = 0.0;
      _seedThetaW= 0.0;
      _seedThetaV= 0.0;
      _seedQC   = 0;

      _recoPx = 0.0;
      _recoPy = 0.0;
      _recoPz = 0.0;
      _recoPpu = 0.0;
      _recoPpv = 0.0;
      _recoPpw = 0.0;
      _recoPpuErr = 0.0;
      _recoPpvErr = 0.0;
      _recoPpwErr = 0.0;
      _recoE = 0.0;
      _recoTime= 0.0;
      _recoEErr = 0.0;
      _recoTimeErr= 0.0;
      _recoPpx = 0.0;
      _recoPpy = 0.0;
      _recoPpz = 0.0;
      _recoPu = 0.0;
      _recoPv = 0.0;
      _recoPw = 0.0;
      _recoThetaW= 0.0;
      _recoThetaV = 0.0;
      _recoPuErr = 0.0;
      _recoPvErr = 0.0;
      _recoPwErr = 0.0;
      _recoThetaWErr= 0.0;
      _recoThetaVErr = 0.0;
      _recoTrkId = 0;
      _recoIndex = 0;
      if(_diagLevel>2){
	cout<<"caloTrk.Id() = "<<it->first<<endl;
      }
      tmpTrkId = it->first;
      if( caloMap.find(tmpTrkId) != caloMap.end() ){//search the tmpTrkId in the caloMap
	_trkIdFound = 1;
	_caloInt = 0;
	_recoTrkId = tmpTrkId;
	if(caloMap[tmpTrkId].size() == 0) continue;

	for(size_t i=0; i< caloMap[tmpTrkId].size(); ++i){
	  tmpVane = caloMap[tmpTrkId].vaneId(i);

	  _isMiss = 1;

	  _pointsInVaneIndex = i;

	  if(_diagLevel>2){
	    cout<<"found matching!"<<endl;
	  }
	  if(QCvec.find(tmpTrkId) ){
	    _seedQC = 1;
	  }else{
	    _seedQC = 0;
	  }
	  Vane const v1 = cg->vane(tmpVane );
	  Hep3Vector momentumRotUnit = (v1.rotation())*(caloMap[tmpTrkId].impMom3Vec( i ).unit());

	  thetaWimpact = std::atan(-1.0*momentumRotUnit.getZ() / momentumRotUnit.getX() ) ;
	  _seedThetaW = thetaWimpact*Constants::radToDegrees;

	  thetaVimpact = std::atan(momentumRotUnit.getY() /  momentumRotUnit.getX() );
	  _seedThetaV = thetaVimpact*Constants::radToDegrees;

	  _seedE  = caloMap[tmpTrkId].impMom3Vec(i).mag();
	  _seedTime = caloMap[tmpTrkId].impTime(i);
	  _seedPx = caloMap[tmpTrkId].impPos(i).x();
	  _seedPy = caloMap[tmpTrkId].impPos(i).y();
	  _seedPz = caloMap[tmpTrkId].impPos(i).z();
	  Hep3Vector vaneFrame = cg->toVaneFrame(caloMap[tmpTrkId].vaneId(i), caloMap[tmpTrkId].impPos(i) );
	  _seedPu = vaneFrame.x();
	  _seedPv = vaneFrame.y();
	  _seedPw = vaneFrame.z();
	  _seedPpx = caloMap[tmpTrkId].impMom3Vec(i ).x();
	  _seedPpy = caloMap[tmpTrkId].impMom3Vec(i ).y();
	  _seedPpz = caloMap[tmpTrkId].impMom3Vec(i ).z();
	  if(_diagLevel>2){
	    cout<<"momentum : "<<caloMap[tmpTrkId].impMom3Vec(i )<<endl;
	  }
	  _seedPpu = caloMap[tmpTrkId].impMom3Vec(i ).mag()*momentumRotUnit.x();
	  _seedPpv = caloMap[tmpTrkId].impMom3Vec(i ).mag()*momentumRotUnit.y();
	  _seedPpw = caloMap[tmpTrkId].impMom3Vec(i ).mag()*momentumRotUnit.z();


	  _Ntup->Fill();

	  if(_diagLevel>2){
	    cout<<"_Ntup filled..."<<endl;
	  }
	}



      } else {//end if(same track Id)
	if(_diagLevel>2){
	  cout<<"not found TrkId..."<<endl;
	}
	_Ntup->Fill();
      }

    }


    if(evt.id().event() % 100 == 0){
      cout << "Event "<<evt.id().event()<<" ReadExtrapol done..."<<endl;
    }
  }



}



using mu2e::ReadExtrapol;
DEFINE_ART_MODULE(ReadExtrapol);

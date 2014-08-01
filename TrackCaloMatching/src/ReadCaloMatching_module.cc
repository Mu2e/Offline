//
//
//
// $Id: ReadCaloMatching_module.cc,v 1.16 2014/08/01 20:57:45 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:45 $
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
#include "KalmanTests/inc/KalRepCollection.hh"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//tracker includes
//#include "TrkBase/TrkRep.hh"
//#include "KalmanTrack/KalRep.hh"
#include "KalmanTests/inc/KalFitMC.hh"
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

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Geometry/HepPoint.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "TrackCaloMatching/inc/TrkToCaloExtrapolCollection.hh"


//calorimeter includes
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "TrackCaloMatching/inc/CaloVolumeElem.hh"
#include "TrackCaloMatching/inc/CaloVolumeType.hh"
#include "TrackCaloMatching/inc/CaloSurface.hh"
#include "CaloCluster/inc/CaloClusterer.hh"
#include "CaloCluster/inc/CaloClusterTools.hh"
#include "TrackCaloMatching/inc/Calorimeter4VanesGeom.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
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

float thetaWimpact(const CLHEP::Hep3Vector& mom, int vaneId){

        GeomHandle<VaneCalorimeter> cg;
        Vane const &vane = cg->vane(vaneId);
        CLHEP::Hep3Vector dirMom_rotated = (vane.rotation())*mom.unit();

        float thW = 0.0;
        thW = std::atan(-1.0*dirMom_rotated.getZ() / dirMom_rotated.getX() ) ;
        thW *= Constants::radToDegrees;
        return thW;
}
float thetaVimpact(const CLHEP::Hep3Vector& mom, int vaneId){//(FIXME)

        GeomHandle<VaneCalorimeter> cg;
        Vane const &vane = cg->vane(vaneId);
        CLHEP::Hep3Vector dirMom_rotated = (vane.rotation())*mom.unit();

        float thV = 0.0;
        thV = std::atan(dirMom_rotated.getY() / dirMom_rotated.getX() ) ;
        thV *= Constants::radToDegrees;
        return thV;
}

float thetaWimpactErr(float& z, float& dz, float& x, float& dx){
        float res = 0.0;
        if(x == 0.0) return res;
        res += dz*fabs(z /( pow(x,2)*(1.0 + pow(z,2)/(pow(x,2) ) ) ));
        res += dx*fabs(1.0/(x*(1.0 + pow(z, 2)/( pow(x, 2) ) ) ) );
        res *= Constants::radToDegrees;
        return res;
}

float thetaVimpactErr(float& y, float& dy, float& x, float& dx){
        float res = 0.0;
        if(x == 0.0) return res;
        res += dy*fabs(y /( pow(x,2)*(1.0 + pow(y,2)/(pow(x,2) ) ) ));
        res += dx*fabs(1.0/(x*(1.0 + pow(y, 2)/( pow(x, 2) ) ) ) );
        res *= Constants::radToDegrees;
        return res;
}

struct simData{
        double entranceEnergy, exitEnergy;
        double cosTheta, cosPitch, lcosPitch, entranceTime, exitTime;
        size_t nHits;
        int isQC;
        simData(){};
};

typedef std::map<unsigned int, simData> SimMap;

static int ncalls(0);

class ReadCaloMatching : public art::EDAnalyzer {
public:
        explicit ReadCaloMatching(fhicl::ParameterSet const& pset):
          art::EDAnalyzer(pset),
        _fitterModuleLabel(pset.get<string>("fitterModuleLabel")),
        _diagLevel(pset.get<int>("diagLevel",0)),
        _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
        _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
        _makerModuleLabel(pset.get<std::string>("makerModuleLabel", "makeSH")),
        _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
        _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
        _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel", "makeCaloCluster")),
        _caloClusterAlgorithm(pset.get<std::string>("caloClusterAlgorithm", "closest")),
        _caloClusterSeeding(pset.get<std::string>("caloClusterSeeding", "energy")),
        _producerName("Algo"+mu2e::TOUpper(_caloClusterAlgorithm)+"SeededBy"+mu2e::TOUpper(_caloClusterSeeding)),
        _elextractModuleLabel(pset.get<std::string>("elextractModuleLabel", "extractElData")),
        _extractElectronsData(pset.get<string>("elextractModuleLabel")),
        _trkToCaloExtrapolModuleLabel( pset.get<std::string>("trkToCaloExtrapolModuleLabel", "TrkExtrapol")),
        _trkCaloMatchingModuleLabel( pset.get<std::string>("trkCaloMatchingModuleLabel", "CaloMatching")),
        _Ntup(0),//_trkCaloMatchingModuleLabel
        _application(nullptr),
        _directory(0){
        }

        virtual ~ReadCaloMatching() {
        }
        void beginJob();
        void endJob() {}

        void analyze(art::Event const& e );

private:

        double chiSquare(double& exV, double& exW, double& clV, double& clW, double& sigmaV2,
                        double& sigmaW2, double& extrT , double& clT, double& sigmaT2,
                        double& exEnergy, double& clEnergy, double& sigmaE2);

        void doExtrapolation(art::Event const& evt, bool skip);
        // Module label of the module that performed the fits.
        std::string _fitterModuleLabel;

        // Diagnostic level
        int _diagLevel;

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

        // Label of the calo clusters  maker
        std::string _caloClusterModuleLabel;

        string _caloClusterAlgorithm;
        string _caloClusterSeeding;
        string _producerName;

        string _elextractModuleLabel;
        string _extractElectronsData;

        // Label of the extrapolated impact points
        std::string _trkToCaloExtrapolModuleLabel;

        // Label of the matched pairs of  <TrkToCaloExtrapol, CaloCluster>
        std::string _trkCaloMatchingModuleLabel;

        bool _skipEvent;

        //const string _producerName;

        TTree* _Ntup;//Ntupla which contains informations about the extrapolation starting from MC

        // The job needs exactly one instance of TApplication.  See note 1.
        unique_ptr<TApplication> _application;

        // Save directory from beginJob so that we can go there in endJob. See note 3.
        TDirectory* _directory;
        //        bool _firstEvent;


        Int_t _match,//MC and Reco trajectories are compatible in a given tolerance
        _recoIndex,
        _recoTrkId,
        _seedQC,//the electron generated satisfied the quality cuts
        _clSize,
        _seedPdgId ,
        _seedIsGen,
	  _seedIsConv,
        _clVane,
        _mcNHits;

        Float_t _evt,//event Id
        _mcTrkEntrE,
        _mcTrkExitE,
        _mcTrkEntrT0,
        _mcTrkExitT0,
        _mcCosTh,
        _mcCosPitch,
        _recoFitCons,
        _clE ,
        _clEErr,
        _clT ,
        _clTErr,
        _clCOGx ,
        _clCOGy ,
        _clCOGz ,
        _clCOGu ,
        _clCOGv ,
        _clCOGw ,
        _clCOGuErr ,
        _clCOGvErr ,
        _clCOGwErr ,
        _clCOGrow,
        _clCOGcolumn,
        _clRows[10000],
        _clColumns[10000],
        _clShowerDir,
        _clErrShowerDir,
        _clCryEnergyMaxRow,
        _clCryEnergyMaxColumn,
        _clCryMaxEdep,
        _clWsize,
        _clVsize,
        _cryEdep[10000],
        _cryEdepTot,
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
	  _recoTrkEntrE,
	  _recoTrkExitE,
        _recoTime,//reconstructed impact impact time [ns]
	  _recoTimeErr,
	  _recoT0TrkEntr,
        _recoT0TrkEntrErr, 
	_recoT0TrkExit,
        _recoT0TrkExitErr,
	  _recoTorigin,
	  _recoToriginErr,
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
        _recoPathLenght,
        _recoPathLenghtErr,
        _recoEt0,
        _recoEt0Err,
	_recoExitPathLenght,
	_recoEntrPathLenght,
	  _recoTrkEntrPpx,  
	  _recoTrkEntrPpy,  
	  _recoTrkEntrPpz, 
	  _recoTrkEntrPx,  
	  _recoTrkEntrPy,  
	  _recoTrkEntrPz, 
	  _recoTrkExitPpx,  
	  _recoTrkExitPpy,  
	  _recoTrkExitPpz,  
	  _recoTrkExitPx,  
	  _recoTrkExitPy,  
	  _recoTrkExitPz,  
	  _recoTrk0FlightPx,  
	  _recoTrk0FlightPy,  
	  _recoTrk0FlightPz,
        _pseudoChiSquare,
        _timeChiSquare,
        _energyChiSquare,
        _posVChiSquare,
        _posWChiSquare;
};


double ReadCaloMatching::chiSquare(double& exV, double& exW, double& clV, double& clW, double& sigmaV2,
                double& sigmaW2, double& extrT , double& clT, double& sigmaT2,
                double& exEnergy, double& clEnergy, double& sigmaE2){

        double res = 0.0;
        _posVChiSquare = pow(2.35*(exV -clV), 2)/sigmaV2;
        _posWChiSquare = pow(2.35*(exW -clW), 2)/sigmaW2;
        _timeChiSquare = pow(2.35*(extrT -clT), 2)/sigmaT2;
        _energyChiSquare = pow(2.35*(exEnergy -clEnergy), 2)/sigmaE2;

        res = _posVChiSquare + _posWChiSquare /*+ _energyChiSquare */+ _timeChiSquare;
        return res;
}

bool findTrkId(std::vector<unsigned int> vec, unsigned int t){
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
        //        if(map.size()==0) return key;
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


CLHEP::Hep3Vector fromTrkToMu2eFrame(CLHEP::Hep3Vector vec){
        art::ServiceHandle<GeometryService> geom;
        double solenoidOffSetX = geom->config().getDouble("mu2e.solenoidOffset");//3904.;//[mm]
        double solenoidOffSetZ = -geom->config().getDouble("mu2e.detectorSystemZ0");//-10200.;

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



void ReadCaloMatching::beginJob( ) {
        cout << "start ReadCaloMatching..."<<endl;
}



void ReadCaloMatching::analyze(art::Event const& evt ) {

        ++ncalls;

        if (ncalls == 1) {

                art::ServiceHandle<art::TFileService> tfs;
                _Ntup        = tfs->make<TTree>("trkCalo", "trk-calo matching info");

                _Ntup->Branch("evt", &_evt , "evt/F");
                _Ntup->Branch("match",     &_match , "match/I");
                _Ntup->Branch("mcNHits",     &_mcNHits , "mcNHits/I");
                _Ntup->Branch("mcTrkEntrE",     &_mcTrkEntrE , "mcTrkEntrE/F");
                _Ntup->Branch("mcTrkExitE",     &_mcTrkExitE , "mcTrkExitE/F");
                _Ntup->Branch("mcTrkEntrT0",     &_mcTrkEntrT0 , "mcTrkEntrT0/F");
                _Ntup->Branch("mcTrkExitT0",     &_mcTrkExitT0 , "mcTrkExitT0/F");
                _Ntup->Branch("mcCosTh",     &_mcCosTh , "mcCosTh/F");
                _Ntup->Branch("mcCosPitch",     &_mcCosPitch , "mcCosPitch/F");
                _Ntup->Branch("recoFitCons",     &_recoFitCons , "recoFitCons/F");
                _Ntup->Branch("seedQC",     &_seedQC , "seedQC/I");
		_Ntup->Branch("seedIsConv",     &_seedIsConv , "seedIsConv/I");
		_Ntup->Branch("seedPdgId",     &_seedPdgId , "seedPdgId/I");
		_Ntup->Branch("seedIsGen",     &_seedIsGen , "seedIsGen/I");
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

                _Ntup->Branch("clSize", &_clSize , "clSize/I");
                _Ntup->Branch("clVane", &_clVane , "clVane/I");
                _Ntup->Branch("clE"   , &_clE    , "clE/F");
                _Ntup->Branch("clEErr"   , &_clEErr    , "clEErr/F");
                _Ntup->Branch("clT"   , &_clT   , "clT/F");
                _Ntup->Branch("clTErr"   , &_clTErr   , "clTErr/F");
                _Ntup->Branch("clCOGx"   , &_clCOGx   , "clCOGx/F");
                _Ntup->Branch("clCOGy"   , &_clCOGy   , "clCOGy/F");
                _Ntup->Branch("clCOGz"   , &_clCOGz   , "clCOGz/F");
                _Ntup->Branch("clCOGu"   , &_clCOGu   , "clCOGu/F");
                _Ntup->Branch("clCOGv"   , &_clCOGv   , "clCOGv/F");
                _Ntup->Branch("clCOGw"   , &_clCOGw   , "clCOGw/F");
                _Ntup->Branch("clCOGuErr"   , &_clCOGuErr   , "clCOGuErr/F");
                _Ntup->Branch("clCOGvErr"   , &_clCOGvErr   , "clCOGvErr/F");
                _Ntup->Branch("clCOGwErr"   , &_clCOGwErr   , "clCOGwErr/F");
                _Ntup->Branch("clCOGrow"   , &_clCOGrow   , "clCOGrow/F");
                _Ntup->Branch("clCOGcolumn"   , &_clCOGcolumn   , "clCOGcolumn/F");
                _Ntup->Branch("clRows[clSize]",     _clRows , "clRows[clSize]/F");
                _Ntup->Branch("clColumns[clSize]",  _clColumns , "clColumns[clSize]/F");
                _Ntup->Branch("clShowerDir"   , &_clShowerDir   , "clShowerDir/F");
                _Ntup->Branch("clErrShowerDir"   , &_clErrShowerDir   , "clErrShowerDir/F");
                _Ntup->Branch("clCryEnergyMaxRow"   , &_clCryEnergyMaxRow   , "clCryEnergyMaxRow/F");
                _Ntup->Branch("clCryEnergyMaxColumn"   , &_clCryEnergyMaxColumn   , "clCryEnergyMaxColumn/F");
                _Ntup->Branch("clCryMaxEdep"   , &_clCryMaxEdep   , "clCryMaxEdep/F");
                _Ntup->Branch("clWsize"   , &_clWsize   , "clWsize/F");
                _Ntup->Branch("clVsize"   , &_clVsize   , "clVsize/F");
                _Ntup->Branch("cryEdepTot"   , &_cryEdepTot   , "cryEdepTot/F");
                _Ntup->Branch("cryEdep[clSize]",  _cryEdep , "cryEdep[clSize]/F");

                _Ntup->Branch("recoTrkId",     &_recoTrkId , "recoTrkId/I");
                _Ntup->Branch("recoIndex",     &_recoIndex , "recoIndex/I");
                _Ntup->Branch("recoT0TrkEntr",      &_recoT0TrkEntr , "recoT0TrkEntr/F");
                _Ntup->Branch("recoT0TrkEntrErr",   &_recoT0TrkEntrErr , "recoT0TrkEntrErr/F");
                _Ntup->Branch("recoT0TrkExit",      &_recoT0TrkExit , "recoT0TrkExit/F");
                _Ntup->Branch("recoT0TrkExitErr",   &_recoT0TrkExitErr , "recoT0TrkExitErr/F");

                _Ntup->Branch("recoTorigin",     &_recoTorigin , "recoTorigin/F");
                _Ntup->Branch("recoToriginErr",     &_recoToriginErr , "recoToriginErr/F");

                _Ntup->Branch("recoPx",     &_recoPx , "recoPx/F");
                _Ntup->Branch("recoPy",     &_recoPy , "recoPy/F");
                _Ntup->Branch("recoPz",     &_recoPz , "recoPz/F");
                _Ntup->Branch("recoE",      &_recoE , "recoE/F");
                _Ntup->Branch("recoEErr",   &_recoEErr , "recoEErr/F");
                _Ntup->Branch("recoTrkEntrE", &_recoTrkEntrE , "recoTrkEntrE/F");
                _Ntup->Branch("recoTrkExitE", &_recoTrkExitE , "recoTrkExitE/F");
                _Ntup->Branch("recoEt0",      &_recoEt0 , "recoEt0/F");
                _Ntup->Branch("recoEt0Err",   &_recoEt0Err , "recoEt0Err/F");
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
		_Ntup->Branch("recoEntrPathLenght", &_recoEntrPathLenght , "recoEntrPathLenght/F");
                _Ntup->Branch("recoExitPathLenght", &_recoExitPathLenght , "recoExitPathLenght/F");
                _Ntup->Branch("recoPathLenght", &_recoPathLenght , "recoPathLenght/F");
                _Ntup->Branch("recoPathLenghtErr", &_recoPathLenghtErr , "recoPathLenghtErr/F");
		_Ntup->Branch("recoTrkEntrPpx", &_recoTrkEntrPpx , "recoTrkEntrPpx/F");
		_Ntup->Branch("recoTrkEntrPpy", &_recoTrkEntrPpy , "recoTrkEntrPpy/F");
		_Ntup->Branch("recoTrkEntrPpz", &_recoTrkEntrPpz , "recoTrkEntrPpz/F");
		_Ntup->Branch("recoTrkExitPx", &_recoTrkExitPx , "recoTrkExitPx/F");
		_Ntup->Branch("recoTrkExitPy", &_recoTrkExitPy , "recoTrkExitPy/F");
		_Ntup->Branch("recoTrkExitPz", &_recoTrkExitPz , "recoTrkExitPz/F");
		_Ntup->Branch("recoTrkExitPpx", &_recoTrkExitPpx , "recoTrkExitPpx/F");
		_Ntup->Branch("recoTrkExitPpy", &_recoTrkExitPpy , "recoTrkExitPpy/F");
		_Ntup->Branch("recoTrkExitPpz", &_recoTrkExitPpz , "recoTrkExitPpz/F");
		_Ntup->Branch("recoTrkExitPx", &_recoTrkExitPx , "recoTrkExitPx/F");
		_Ntup->Branch("recoTrkExitPy", &_recoTrkExitPy , "recoTrkExitPy/F");
		_Ntup->Branch("recoTrkExitPz", &_recoTrkExitPz , "recoTrkExitPz/F");
		_Ntup->Branch("recoTrk0FlightPx", &_recoTrk0FlightPx , "recoTrk0FlightPx/F");
		_Ntup->Branch("recoTrk0FlightPy", &_recoTrk0FlightPy , "recoTrk0FlightPy/F");
		_Ntup->Branch("recoTrk0FlightPz", &_recoTrk0FlightPz , "recoTrk0FlightPz/F");
		
                _Ntup->Branch("pseudoChiSquare",  &_pseudoChiSquare , "pseudoChiSquare/F");
                _Ntup->Branch("timeChiSquare",    &_timeChiSquare , "timeChiSquare/F");
                _Ntup->Branch("energyChiSquare",  &_energyChiSquare , "energyChiSquare/F");
                _Ntup->Branch("posVChiSquare",    &_posVChiSquare , "posVChiSquare/F");
                _Ntup->Branch("posWChiSquare",    &_posWChiSquare , "posWChiSquare/F");
        }


        doExtrapolation(evt, _skipEvent);

} // end of analyze


void ReadCaloMatching::doExtrapolation(art::Event const& evt, bool skip){

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


        art::Handle<TrackClusterLink>  trjCaloMatchings;
        evt.getByLabel(_trkCaloMatchingModuleLabel, trjCaloMatchings);

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
        SimMap simMap;

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
                simData sData;
                sData.nHits = iEltrk.getNumOfHit() ;
                sData.cosTheta = cosTheta ;
                sData.cosPitch = cosPitch;
                sData.lcosPitch = lcosPitch;
                sData.entranceEnergy = hdil._hitMomentum.mag();
                sData.exitEnergy = ldil._hitMomentum.mag();
                sData.isQC = 0;
                sData.entranceTime = hdil._mcHitTime;
                sData.exitTime = ldil._mcHitTime;
		

                //the following condition are the same used by Dave Brown for TTracker studies
                condition &= ( iEltrk.getNumOfHit() >= 20 );
                condition &= ( hdil._hitMomentum.mag() >= trkMomCut );
                condition &= ( cosTheta >= -0.5 );
                condition &= ( cosTheta <=  0.5 );
                condition &= ( cosPitch > 0.5 );
                condition &= ( cosPitch < 0.70710678118655 );// 1 / sqrt(2)



                if( condition ){
                        sData.isQC = 1;
                        NtrkCut++;// # of electrons at the entrance of the TRK which have trkMomCut MeV of kinetic energy and which will light at least 20 straws
                        if( !QCvec.find(iEltrk.getTrkID().asUint()) ){
                                if(_diagLevel>2){
                                        cout<<"!QCvec.find(iEltrk.getTrkID().asUint())"<<
                                                        ", trkId = "<< iEltrk.getTrkID().asUint()<<endl;
                                }
                                QCvec.push_back(iEltrk.getTrkID().asUint());
                        }


                }//end if(condition)

                simMap[iEltrk.getTrkID().asUint()] =  sData ;

        }//end loop TRK mapping
        if (NtrkTot==0) return;


        ElecMap recoMap, caloMap, extrMap;

        int iVane = -1;
        if (!( strawHandle.isValid())) {
                return;
        }

        Hep3Vector  tmpPerr, vaneFrame;

        int tmpVane = -1, tmpCluCryMaxEdepColumn = -1, tmpCluCryMaxEdepRow = -1;
        int tmpWsize = -1;
        double thV = 0.0;

        double tmpCOGv = 0.0, tmpCOGw = 0.0, tmpCOGvErr = 0.0, tmpCOGwErr = 0.0, tmpCluTime = 0.0;
        double tmpPv = 0.0, tmpPw = 0.0,tmpPtime = 0.0, tmpPtimeErr = 0.0;
        CLHEP::Hep3Vector cogVaneFrame, tmpPosVaneFrame;
        CLHEP::Hep3Vector cogVaneFrameErr, tmpPosVaneFrameErr;
        double tmpEnergy = 0.0, tmpEnergyErr = 0.0;
        double tmpErrPos = 0.0;
        std::vector<int> trjsTrkId, clustersTrkId;

        int count = 0;

        for(size_t i =0; i<trjCaloMatchings->size(); ++i){
                elecData tmpElec;
                KalRepPtr const& trkPtr = trjCaloMatchings->at(i).first->trk();
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

                TrkHotList const* hots = trk->hotList();

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

                        if(mcptr.size() == 0) continue;

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

                if(StrawTracksMap.size()==0) {
                        trjsTrkId.push_back(-1);
                        continue;
                }
                size_t trkId = MaxKey(StrawTracksMap);
                trjsTrkId.push_back((int)trkId);

                if(_diagLevel>2){
                        cout<<"filling tmpElec..."<<endl;
                }
                tmpElec._vaneId        = trjExtrapols->at(i).vaneId();
                tmpElec._index         = count;
                tmpElec._fitConsistency= trjExtrapols->at(i).fitConsistency();
                tmpElec._impTime       = trjExtrapols->at(i).time();
                tmpElec._impTimeErr    = trjExtrapols->at(i).timeErr();
                tmpElec._genTime       = trjExtrapols->at(i).tOrigin();
                tmpElec._genTimeErr    = trjExtrapols->at(i).tOriginErr();
                tmpElec._impMom3Vec    = trjExtrapols->at(i).momentum();
                tmpElec._impMom3VecErr = trjExtrapols->at(i).momentumErr();
                tmpElec._impPos        = Hep3Vector(trjExtrapols->at(i).entrancePosition().x(), trjExtrapols->at(i).entrancePosition().y(), trjExtrapols->at(i).entrancePosition().z());
                tmpElec._impPosErr     = trjExtrapols->at(i).entrancePositionErr();

                recoMap[trkId].push_back(tmpElec);
                if(_diagLevel>2){
                        cout<<"recoMap filled..."<<endl;
                }

        }
//-----------------------------------------------------------------------------
// 2013-05-17 P.Murat: commented the -r 1.10 line , uncomment it in the near future
//-----------------------------------------------------------------------------
	double tollCaloHits = cg->caloGeomInfo().crystalHalfTrans()*15.0;
	//	double tollCaloHits = cg->crystalHalfSize()*15.0;

        for(size_t i =0; i<trjCaloMatchings->size(); ++i){

                CaloCluster const& clu = *( (*trjCaloMatchings).at(i).second.get());

                elecData tmpElec;
                bool foundGen =false;
                CaloCrystalHitPtrVector caloClusterHits = clu.caloCrystalHitsPtrVector();

                unsigned int trkId = 0;

                for(size_t i=0; i<caloClusterHits.size(); ++i){
                        CaloCrystalHit const& hit = *(caloClusterHits.at(i));

                        std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
                        if(ROIds.size()<1 ) continue;

                        CaloHit const& thehit = *ROIds.at(0);
                        size_t collectionPosition = ROIds.at(0).key();

                        PtrStepPointMCVector const & mcptr(hits_mcptr->at(collectionPosition));
                        if(mcptr.size() <=0) continue;

                        if(_diagLevel>3){
                                cout<<"step "<< i <<" : "<< caloCrystalHits->size()<<endl;
                                thehit.print();
                                cout<<endl;
                        }


                        iVane = cg->vaneByRO(thehit.id());

                        size_t nHitsPerCrystal = mcptr.size();


                        for (size_t j2=0; j2<nHitsPerCrystal; ++j2) {

                                StepPointMC const& mchit = *mcptr[j2];

                                art::Ptr<SimParticle> const& simptr = mchit.simParticle();


                                // The simulated particle that made this hit.
                                SimParticleCollection::key_type trackId(mchit.trackId());
                                SimParticle const& sim = *simptr;

                                if(sim.fromGenerator() ){
                                        if(tmpElec._impTime > mchit.time() ){
                                                if(sim.genParticle().isNonnull()){
                                                        GenParticle const& gen = *sim.genParticle();
                                                        tmpElec._genMom = gen.momentum();
                                                }
                                                foundGen               = true;
                                                trkId                  = trackId.asUint();
                                                tmpElec._vaneId        = iVane;
						tmpElec._pdgId         = sim.pdgId();
						tmpElec._isGen         = sim.fromGenerator();
						tmpElec._impTime       = mchit.time();
                                                tmpElec._impMom3Vec    = mchit.momentum();
                                                tmpElec._impPos        = mchit.position();
						if(sim.fromGenerator() ){
						  GenParticle const& gen = *sim.genParticle();
						  GenId genId = gen.generatorId();
						  if(genId==GenId::conversionGun){
						    tmpElec._isConv        =  1;
						  }
						}
                                        }//end  if(caloMap[trackId.asUint()][iVane]._impTime > mchit.time() ){

                                }else{//end if(sim.fromGenerator)
                                }

                        }//end loop on nHitsPerCrystal
                }
                if(foundGen){
                        clustersTrkId.push_back(trkId);

                        if(_diagLevel>2){
                                cout<<"foundGen!!"<<endl;
                        }

                        if(caloMap.size()!= 0 ){
                                if(_diagLevel>2){
                                        cout<<"caloMap.size()!= 0..."<<endl;
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
                                }
                        }else{
                                if(_diagLevel>2){
                                        cout<<" caloMap.size() = 0 "<<
                                                        ", trkId = "<<trkId<<endl;
                                        tmpElec.print(cout);
                                }
                                caloMap[trkId].push_back(tmpElec);
                                if(_diagLevel>2){
                                        cout<<"  caloMap[trkId].push_back(tmpElec) done... "<<endl;
                                }
                        }
                }else{
                        clustersTrkId.push_back(-1);
                }

        }



        if(_diagLevel > 2){
                cout<<"trjCaloMatchings.size =  "<<trjCaloMatchings->size()<<
                                ", clustersTrkId.size() = "<< clustersTrkId.size()<<
                                ", trjsTrkId.size() = "<< trjsTrkId.size()<<endl;
        }

        count = 0;
        float thetaV = 0.0, thetaW = 0.0;
        double sigmaV2 = 0.0, sigmaW2 = 0.0, sigmaT2 = 0.0, sigmaE2 = 0.0, tmpCluEnergy = 0.0;
        float chiQ = 0.0;
        for(size_t i =0; i<trjCaloMatchings->size(); ++i){

                //initialyze all the branches
                _match = 0;
                _seedQC = 0;
                _clSize = 0;
                _seedPdgId = 0;
                _seedIsGen = 0;
                _clVane = 0;

                _evt =  0.0;
                _clE =  0.0;
                _clEErr =  0.0;
                _clT =  0.0;
                _clTErr =  0.0;
                _clCOGx =  0.0;
                _clCOGy =  0.0;
                _clCOGz =  0.0;
                _clCOGu =  0.0;
                _clCOGv =  0.0;
                _clCOGw =  0.0;
                _clCOGuErr =  0.0;
                _clCOGvErr =  0.0;
                _clCOGwErr =  0.0;
                _clCOGrow =  0.0;
                _clCOGcolumn =  0.0;

                _clShowerDir =  0.0;
                _clErrShowerDir =  0.0;
                _clCryEnergyMaxRow=  0.0;
                _clCryEnergyMaxColumn =  0.0;
                _clCryMaxEdep =  0.0;
                _clWsize =  0.0;
                _clVsize =  0.0;
                _cryEdepTot =  0.0;

                _seedQC = 0;
                _seedPx =  0.0;
                _seedPy =  0.0;
                _seedPz =  0.0;
                _seedPpu =  0.0;
                _seedPpv =  0.0;
                _seedPpw =  0.0;
                _seedE =  0.0;
                _seedTime =  0.0;
                _seedPpx =  0.0;
                _seedPpy =  0.0;
                _seedPpz =  0.0;
                _seedPu =  0.0;
                _seedPv =  0.0;
                _seedPw =  0.0;
                _seedThetaW =  0.0;
                _seedThetaV =  0.0;
		_seedIsConv = 0;
		_seedPdgId = 0;
		_seedIsGen = 0;

		_recoTorigin = 0.0,
		  _recoToriginErr = 0.0,
		_recoT0TrkEntr = 0.0;
		_recoT0TrkEntrErr = 0.0;
                _recoTrkId = 0;
                _recoPx =  0.0;
                _recoPy =  0.0;
                _recoPz =  0.0;
                _recoPpu =  0.0;
                _recoPpuErr =  0.0;
                _recoPpv =  0.0;
                _recoPpvErr =  0.0;
                _recoPpw =  0.0;
                _recoPpwErr =  0.0;
                _recoE =  0.0;
                _recoEErr =  0.0;
                _recoTime =  0.0;
                _recoTimeErr=  0.0;
                _recoPpx =  0.0;
                _recoPpy =  0.0;
                _recoPpz =  0.0;
                _recoPu =  0.0;
                _recoPv =  0.0;
                _recoPw =  0.0;
                _recoPuErr =  0.0;
                _recoPvErr =  0.0;
                _recoPwErr =  0.0;
                _recoThetaW =  0.0;
                _recoThetaWErr =  0.0;
                _recoThetaV =  0.0;
                _recoThetaVErr =  0.0;
                _recoFitCons = 0.0;
                _recoEt0 = 0.0;
                _recoEt0Err = 0.0;
                _recoPathLenght = 0.0;
                _recoPathLenghtErr = 0.0;

                _mcCosPitch     = 0.0;
                _mcCosTh        = 0.0;
                _mcNHits        = 0;
                _mcTrkEntrE     = 0.0;
                _mcTrkExitE     = 0.0;
                _mcTrkEntrT0    = 0.0;
                _mcTrkExitT0    = 0.0;
                //end

                CaloCluster const& clu = *( (*trjCaloMatchings).at(i).second.get());
                TrkToCaloExtrapol const& trkToCaloExtrapol = *( (*trjCaloMatchings).at(i).first.get() );
                _evt = evt.id().event();
                Hep3Vector momentumRotUnit;
                if(_diagLevel > 2){
                        cout<<"trjCaloMatchings at "<<i<<" step"<<endl;
                }

                _recoTrkId = trjsTrkId[i];

                if(clustersTrkId[i] == trjsTrkId[i]){
                        _match = 1;
                }

                if(QCvec.find( clustersTrkId[i] ) ){
                        _seedQC = 1;
                }

                tmpVane = clu.vaneId();


                if(simMap.find((unsigned int)trjsTrkId[i]) != simMap.end() ){
                        SimMap::iterator it = simMap.find((unsigned int)trjsTrkId[i]);
                        _mcCosPitch     = it->second.cosPitch;
                        _mcCosTh        = it->second.cosTheta;
                        _mcNHits        = it->second.nHits;
                        _mcTrkEntrE     = it->second.entranceEnergy;
                        _mcTrkExitE     = it->second.exitEnergy;
                        _mcTrkEntrT0    = it->second.entranceTime;
                        _mcTrkExitT0    = it->second.exitTime;

                }

                if(caloMap.size() != 0){
                        if(_diagLevel > 2){
                                cout<<"caloMap.size() = "<<caloMap.size()<<
                                                ", trkId = "<< clustersTrkId[i]<<
                                                ", caloMap->first = "<< caloMap.begin()->first<< endl;
                        }
                        if(caloMap.find((unsigned int)clustersTrkId[i]) != caloMap.end()){
                                if(_diagLevel > 2){
                                        cout<<"found trkId in caloMap!!"<< endl;
                                }
                                size_t ind = 0;
                                Vane const &vane1 = cg->vane(tmpVane );
                                momentumRotUnit = (vane1.rotation())*(caloMap[(unsigned int)clustersTrkId[i]].impMom3Vec(ind).unit());

                                thetaW = std::atan(-1.0*momentumRotUnit.getZ() / momentumRotUnit.getX() ) ;
                                _seedThetaW = thetaW*Constants::radToDegrees;

                                thetaV = std::atan(momentumRotUnit.getY() /  momentumRotUnit.getX() );
                                _seedThetaV = thetaV*Constants::radToDegrees;


                                _seedTime = caloMap[(unsigned int)clustersTrkId[i]].impTime(ind);
                                _seedE    = caloMap[(unsigned int)clustersTrkId[i]].impMom3Vec(ind).mag();

                                _seedIsConv = caloMap[(unsigned int)clustersTrkId[i]].isConv(ind);
                                _seedPdgId = caloMap[(unsigned int)clustersTrkId[i]].pdgId(ind);
                                _seedIsGen = caloMap[(unsigned int)clustersTrkId[i]].isGen(ind);
                                if(_diagLevel > 2){
                                        cout<<"found trkId in caloMap!! step 1"<< endl;
                                }
                                Hep3Vector impPos, impMom;
                                impPos = caloMap[(unsigned int)clustersTrkId[i]].impPos(ind);
                                impMom = caloMap[(unsigned int)clustersTrkId[i]].impMom3Vec(ind);

                                if(_diagLevel > 2){
                                        cout<<"found trkId in caloMap!! step 2"<<
                                                        ", impMom = "<< impMom<<
                                                        ", impPos = "<< impPos<<endl;
                                }
                                _seedPx = impPos.x();
                                _seedPy = impPos.y();
                                _seedPz = impPos.z();
                                if(_diagLevel > 2){
                                        cout<<"found trkId in caloMap!! step 3"<< endl;
                                }
                                _seedPpx = impMom.x();
                                _seedPpy = impMom.y();
                                _seedPpz = impMom.z();
                                if(_diagLevel > 2){
                                        cout<<"found trkId in caloMap!! step 4"<< endl;
                                }

                                Hep3Vector timpPos = cg->toVaneFrame(tmpVane, impPos);
                                Hep3Vector timpMom = cg->toVaneFrame(tmpVane,  impMom);
                                if(_diagLevel > 2){
                                        cout<<"found trkId in caloMap!! step 5"<<
                                                        ", impMom = "<< timpMom<<
                                                        ", impPos = "<< timpPos<<endl;
                                }
                                _seedPu = timpPos.x();
                                _seedPv = timpPos.y();
                                _seedPw = timpPos.z();

                                _seedPpu = timpMom.x();
                                _seedPpv = timpMom.y();
                                _seedPpw = timpMom.z();

                                if(_diagLevel > 2){
                                        cout<<"found trkId in caloMap!! step 5"<< endl;
                                }

                        }
                }
                if(_diagLevel > 2){
                        cout<<"uscito dal find in caloMap..."<<endl;
                }




                KalRepPtr const& trkPtr = trkToCaloExtrapol.trk();
                const KalRep *  const &trk = *trkPtr;
                if(i>0){
                        TrkToCaloExtrapol const& TMPtrkToCaloExtrapol = *( (*trjCaloMatchings).at(i-1).first.get() );

                        KalRepPtr const& tmpTrkPtr = TMPtrkToCaloExtrapol.trk();
                        const KalRep *  const &tmpTrk = *tmpTrkPtr;
                        if(trk == tmpTrk){
                                ++count;
                        }else{
                                count = 0;
                        }
                }
                _recoIndex = count;

                tmpVane = trkToCaloExtrapol.vaneId();
                if(_diagLevel > 2){
                        cout<<"trkVane = "<< tmpVane <<
                                        ", caloVane = "<< clu.vaneId()<<endl;
                }
                tmpEnergy = trkToCaloExtrapol.momentum().mag();
                _recoE = tmpEnergy;

                _recoEt0 = trkToCaloExtrapol.t0Momentum().mag();

                _recoFitCons = trkToCaloExtrapol.fitConsistency();

                Hep3Vector momdir;
                HepVector momvec(3);
                momdir = trkToCaloExtrapol.t0Momentum().unit();
                for(int icor=0;icor<3;icor++){
                        momvec[icor] = momdir[icor];
                }
                tmpEnergyErr = sqrt(trkToCaloExtrapol.momentumErr().covMatrix().similarity(momvec));
                _recoEt0Err = tmpEnergyErr;

                momdir = trkToCaloExtrapol.momentum().unit();
                for(int icor=0;icor<3;icor++){
                        momvec[icor] = momdir[icor];
                }
                tmpEnergyErr = sqrt(trkToCaloExtrapol.momentumErr().covMatrix().similarity(momvec));
                _recoEErr = tmpEnergyErr;




                Vane const &vane = cg->vane(tmpVane );
                if(_diagLevel > 2){
                        cout<<"momemntum trk : "<<trkToCaloExtrapol.momentum()<<endl;
                }
                momentumRotUnit = (vane.rotation())*(trkToCaloExtrapol.momentum().unit());
                if(_diagLevel > 2){
                        cout<<"momemntum trk : "<<trkToCaloExtrapol.momentum()<<endl;
                }
                thetaW = std::atan(-1.0*momentumRotUnit.getZ() / momentumRotUnit.getX() ) ;
                _recoThetaW = thetaW*Constants::radToDegrees;

                thetaV = std::atan(momentumRotUnit.getY() /  momentumRotUnit.getX() );
                _recoThetaV = thetaV*Constants::radToDegrees;

                //_recoThetaV = thetaVimpact(trkToCaloExtrapol.momentum(), trkToCaloExtrapol.vaneId());
                //_recoThetaW = thetaWimpact(trkToCaloExtrapol.momentum(), trkToCaloExtrapol.vaneId());

                tmpPtime = trkToCaloExtrapol.time();
		
		_recoTorigin = trkToCaloExtrapol.tOrigin();
		_recoToriginErr = trkToCaloExtrapol.tOriginErr();

                tmpPtimeErr = trk->t0().t0Err();

		
                const TrkStrawHit* lasthit = dynamic_cast<const TrkStrawHit*>( trk->lastHit()->kalHit()->hitOnTrack() );
                double fltlen = lasthit->fltLen();
                _recoExitPathLenght = fltlen;
		_recoTrkExitPpx = trk->momentum(fltlen).x();
		_recoTrkExitPpy = trk->momentum(fltlen).y();
		_recoTrkExitPpz = trk->momentum(fltlen).z();
		_recoTrkExitE   = trk->momentum(fltlen).mag();

		_recoTrkExitPx = trk->position(fltlen).x();
		_recoTrkExitPy = trk->position(fltlen).y();
		_recoTrkExitPz = trk->position(fltlen).z();

		_recoTrk0FlightPx = trk->position(0.0).x();
		_recoTrk0FlightPy = trk->position(0.0).y();
		_recoTrk0FlightPz = trk->position(0.0).z();

		_recoT0TrkExit = trk->arrivalTime(fltlen);
		_recoT0TrkExitErr = trk->t0().t0Err();
		
		const TrkStrawHit* firsthit = dynamic_cast<const TrkStrawHit*>( trk->firstHit()->kalHit()->hitOnTrack() );
                fltlen = firsthit->fltLen();
		_recoTrkEntrPpx = trk->momentum(fltlen).x();
		_recoTrkEntrPpy = trk->momentum(fltlen).y();
		_recoTrkEntrPpz = trk->momentum(fltlen).z();
		_recoTrkEntrE   = trk->momentum(fltlen).mag();

		_recoTrkEntrPx = trk->position(fltlen).x();
		_recoTrkEntrPy = trk->position(fltlen).y();
		_recoTrkEntrPz = trk->position(fltlen).z();

		_recoEntrPathLenght = fltlen;

		_recoT0TrkEntr = trk->arrivalTime(fltlen);// trkToCaloExtrapol.t0();
		_recoT0TrkEntrErr = trk->t0().t0Err();// trkToCaloExtrapol.t0Err();

		_recoTime = tmpPtime;
                _recoTimeErr = tmpPtimeErr;

                _recoPathLenght = trkToCaloExtrapol.pathLengthEntrance();
                _recoPathLenghtErr = trkToCaloExtrapol.pathLenghtEntranceErr();

                if(_diagLevel > 2){
                        cout<<"tmpVane = "<<tmpVane<<endl;
                        cout<<"recoThetaW = "<< _recoThetaW<<endl<<
                                        ", recoThetaV = "<< _recoThetaV<<endl;
                        cout<<"tmpEnergy = "<<tmpEnergy<<" [MeV]"<<
                                        ", tmpEnergyErr = "<< tmpEnergyErr <<" [MeV]"<<endl;
                        cout<<"thetaVimpact = "<< thV<<" [deg]"<<endl;
                        cout<<"tmpPtime = "<<tmpPtime<<" [ns]"<<
                                        "tmpPtimeErr = "<< tmpPtimeErr<<" [ns]"<<endl;
                }

                tmpPosVaneFrame.setX(trkToCaloExtrapol.entrancePosition().x());
                tmpPosVaneFrame.setY(trkToCaloExtrapol.entrancePosition().y());
                tmpPosVaneFrame.setZ(trkToCaloExtrapol.entrancePosition().z());
                tmpPosVaneFrame = fromTrkToMu2eFrame(tmpPosVaneFrame);

                if(_diagLevel > 2){
                        cout<<"Mu2e general frame: tmpPosVaneFrame = "<<tmpPosVaneFrame<<" [mm]"<<endl;
                }

                _recoPx = tmpPosVaneFrame.x();
                _recoPy = tmpPosVaneFrame.y();
                _recoPz = tmpPosVaneFrame.z();

                vaneFrame = cg->toVaneFrame(tmpVane, tmpPosVaneFrame);
		tmpPv           = vaneFrame.y();
		tmpPw           = vaneFrame.z();
                _recoPu = vaneFrame.x();
                _recoPv = vaneFrame.y();
                _recoPw = vaneFrame.z();

                //define the local axes V and W for the local vane frame
                Hep3Vector Vaxes(0.0, 1.0, 0.0), Waxes(0.0, 0.0, 1.0);

                //move these axes into the Mu2e general frame
                Vaxes = (vane.rotation())*(Vaxes);
                Waxes = (vane.rotation())*(Waxes);
                //Vaxes = cg->fromVaneFrame(tmpVane, Vaxes);
                //Waxes = cg->fromVaneFrame(tmpVane, Waxes);
                double scaleErrW = 1.0/fabs( cos(thetaW) );
                double scaleErrV = 1.0/fabs( cos(thetaV) );
                momvec[0] = Vaxes.x();
                momvec[1] = Vaxes.y();
                momvec[2] = Vaxes.z();

                tmpErrPos = sqrt(trkToCaloExtrapol.entrancePositionErr().covMatrix().similarity(momvec));
                _recoPvErr = tmpErrPos*scaleErrV;

                momvec[0] = Waxes.x();
                momvec[1] = Waxes.y();
                momvec[2] = Waxes.z();

                tmpErrPos = sqrt(trkToCaloExtrapol.entrancePositionErr().covMatrix().similarity(momvec));
                _recoPwErr = tmpErrPos*scaleErrW;

                //                momvec[0] = 1.0;
                //                momvec[1] = 0.0;
                //                momvec[2] = 0.0;
                //                tmpErrPos = sqrt(trkToCaloExtrapol.entrancePositionErr().covMatrix().similarity(momvec));
                //                tmpPosVaneFrameErr.setX(tmpErrPos);//p.x());
                //
                //                momvec[0] = 0.0;
                //                momvec[1] = 1.0;
                //                tmpErrPos = sqrt(trkToCaloExtrapol.entrancePositionErr().covMatrix().similarity(momvec));
                //                tmpPosVaneFrameErr.setY(tmpErrPos);//.y());
                //
                //                momvec[2] = 1.0;
                //                momvec[1] = 0.0;
                //                tmpErrPos = sqrt(trkToCaloExtrapol.entrancePositionErr().covMatrix().similarity(momvec));
                //                tmpPosVaneFrameErr.setZ(tmpErrPos);//trjExtrapols->at(i).entrancePositionErr().z());
                //
                //                CLHEP::Hep3Vector tmpPerr = *(vane.rotation())*tmpPosVaneFrameErr;//cg->toVaneFrame(tmpVane, tmpPosVaneFrameErr);
                _recoPuErr = 0.0;// fabs(tmpPerr.x());
                // _recoPvErr = fabs(tmpPerr.y());
                //_recoPwErr = fabs(tmpPerr.z());

                tmpPerr = trkToCaloExtrapol.momentum();
                _recoPpx = tmpPerr.x();
                _recoPpy = tmpPerr.y();
                _recoPpz = tmpPerr.z();

                _recoPpu = tmpPerr.mag()*momentumRotUnit.x();
                _recoPpv = tmpPerr.mag()*momentumRotUnit.y();
                _recoPpw = tmpPerr.mag()*momentumRotUnit.z();


                _recoThetaWErr = thetaWimpactErr(_recoPpw, _recoPpwErr, _recoPpu, _recoPpuErr);
                _recoThetaVErr = thetaVimpactErr(_recoPpv, _recoPpvErr, _recoPpu, _recoPpuErr);

                if(_diagLevel > 2){
                        cout<<"Tracker frame: tmpPosVaneFrame = "<<tmpPosVaneFrame<<" [mm]"<<
                                        "tmpPosVaneFrameErr = "<< tmpPosVaneFrameErr<<" [mm]"<<endl;
			cout<<"Track moemntum = "<< tmpPerr <<endl;
                }

                CaloClusterTools cluTool(clu);
                cogVaneFrame            = clu.cog3Vector();
		cogVaneFrame            = cg->toVaneFrame(tmpVane, cogVaneFrame);
                tmpCluCryMaxEdepRow     = cluTool.cryEnergydepMaxRow();
                tmpCluCryMaxEdepColumn  = cluTool.cryEnergydepMaxColumn();
                tmpCluTime              = clu.time();
                tmpCOGv                 = cogVaneFrame.y();
                tmpCOGw                 = cogVaneFrame.z();
                tmpCOGvErr              = clu.cog3VectorError().y();
                tmpCOGwErr              = clu.cog3VectorError().z();

		_clVane = clu.vaneId();

                tmpWsize                = cluTool.wSize();

                //tmpCluEnergy            = clu.energyDep();
                _clE = clu.energyDep();
                tmpCluEnergy = (double)_clE;
                _clEErr = cluTool.energyDepErr();//clusterEnergyErr(_clE);
                //                tmpCluEnergyErr = clusterEnergyErr(_clE);
                if(_diagLevel > 2){
                        cout<<" tmpPv = "<< tmpPv <<
                                        ", tmpPw = "<<tmpPw<<endl<<
                                        ", tmpCOGv = "<<tmpCOGv<<
                                        ", tmpCOGw = "<<tmpCOGw<<endl<<
                                        ", clE = "<< _clE<<
                                        ", clEErr = "<< _clEErr<<endl<<
                                        ", tmpPtime = "<<tmpPtime<<
                                        ", tmpCluTime = "<<tmpCluTime<<
                                        ", tmpWsize = "<<tmpWsize<<endl;
                }

                v_correction_0(_recoThetaV, tmpCOGv, tmpCOGvErr);
                w_correction_0(tmpCOGw, tmpCOGwErr, tmpCluCryMaxEdepColumn);
                w_correction_1(tmpCOGw,tmpCOGwErr, tmpWsize);
                tmpCOGvErr = 15.35 / 2.35;
                tmpCOGwErr = 15.0 / 2.35;


                if(_diagLevel > 2){
                        cout<<"after doing many corrections (0 - correction)..."<<endl;
                        cout<<" tmpPv = "<< tmpPv <<
                                        ", tmpCOGv = "<<tmpCOGv<<
                                        //         ", sigmaV2 = "<<sigmaV2<<endl<<
                                        ", tmpPw = "<<tmpPw<<
                                        ", tmpCOGw = "<<tmpCOGw<<
                                        //         ", sigmaW2 = "<<sigmaW2<<endl<<
                                        ", tmpPtime = "<<tmpPtime<<
                                        ", tmpCluTime = "<<tmpCluTime<<endl;
                        //         ", sigmaT2 = "<<sigmaT2<<endl;
                }

                //v_correction_1(tmpCluCryMaxEdepRow, tmpCOGv, tmpCOGvErr);

                if(_diagLevel > 2){
                        cout<<"after doing other corrections (1 - correction)..."<<endl;
                        cout<<" tmpPv = "<< tmpPv <<
                                        ", tmpCOGv = "<<tmpCOGv<<
                                        ", tmpCOGvErr = "<< tmpCOGvErr<<
                                        //         ", sigmaV2 = "<<sigmaV2<<endl<<
                                        ", tmpPw = "<<tmpPw<<
                                        ", tmpCOGw = "<<tmpCOGw<<
                                        ", tmpCOGwErr = "<< tmpCOGwErr<<
                                        //         ", sigmaW2 = "<<sigmaW2<<endl<<
                                        ", tmpPtime = "<<tmpPtime<<
                                        ", tmpCluTime = "<<tmpCluTime<<endl;
                        //         ", sigmaT2 = "<<sigmaT2<<endl;
                }

                _clCOGu = clu.cog3Vector().x();
                _clCOGuErr = clu.cog3VectorError().x();

                _clCOGv = tmpCOGv;
                _clCOGvErr = tmpCOGvErr;

                _clCOGw = tmpCOGw;
                _clCOGwErr = tmpCOGwErr;

                _clCryEnergyMaxRow = tmpCluCryMaxEdepRow;
                _clCryEnergyMaxColumn = tmpCluCryMaxEdepColumn;

                _clT = clu.time();
                _clTErr = cluTool.timeErr();

                sigmaE2 = std::pow(7.2, 2);//cet::sum_of_squares((float)tmpEnergyErr, _clEErr);
                sigmaV2 = std::pow(44.33, 2);//cet::sum_of_squares(_recoPvErr, (float)tmpCOGvErr);
                sigmaW2 = std::pow(22.17, 2);//cet::sum_of_squares(_recoPwErr, (float)tmpCOGwErr);
                sigmaT2 = std::pow(3.4, 2);//cet::sum_of_squares((float)tmpPtimeErr, _clTErr);

                if(_diagLevel > 2){
                        cout<<"after doing other corrections..."<<endl;
                        cout<<" tmpPv = "<< tmpPv <<
                                        ", tmpCOGv = "<<tmpCOGv<<
                                        //  ", sigmaV2 = "<<sigmaV2<<endl<<
                                        ", tmpPw = "<<tmpPw<<
                                        ", tmpCOGw = "<<tmpCOGw<<
                                        //    ", sigmaW2 = "<<sigmaW2<<endl<<
                                        ", tmpPtime = "<<tmpPtime<<
                                        ", tmpCluTime = "<<tmpCluTime<<endl;
                        //    ", sigmaT2 = "<<sigmaT2<<endl;
                }
                _clSize = clu.size();

                double tmpClVmin = cg->nCrystalR(), tmpClVmax=0., tmpClWmin = cg->nCrystalZ(), tmpClWmax = 0.;

                for(int i = 0; i<_clSize; ++i ){
                        _cryEdep[i] = clu.caloCrystalHitsPtrVector().at(i)->energyDep();//_cryEnergyDep;

                        CaloCrystalHit const& hit = *( clu.caloCrystalHitsPtrVector().at(i) );

                        std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
                        if(ROIds.size()<1 ) continue;

                        CaloHit const& thehit = *ROIds.at(0);

                        if(_cryEdep[i] > _clCryMaxEdep){
                                _clCryMaxEdep = _cryEdep[i];
                        }

                        _clRows[i] = cg->crystalRByRO(thehit.id());

                        if(_clRows[i] < tmpClVmin) {
                                tmpClVmin = _clRows[i];
                        }
                        if(_clRows[i] > tmpClVmax) {
                                tmpClVmax = _clRows[i];
                        }

                        _clColumns[i] =cg->crystalZByRO(thehit.id());

                        if(_clColumns[i] < tmpClWmin) {
                                tmpClWmin = _clColumns[i];
                        }
                        if(_clColumns[i] > tmpClWmax) {
                                tmpClWmax = _clColumns[i];
                        }

                }

                _clWsize = tmpClWmax - tmpClWmin + 1.;

                _clVsize = tmpClVmax - tmpClVmin + 1.;

                _clShowerDir = cluTool.showerDir();
                _clErrShowerDir = cluTool.errShowerDir();

                chiQ = chiSquare(tmpPv, tmpPw, tmpCOGv, tmpCOGw,
                                sigmaV2, sigmaW2, tmpPtime, tmpCluTime,
                                sigmaT2, tmpEnergy, tmpCluEnergy, sigmaE2);

                _pseudoChiSquare = chiQ;

                _Ntup->Fill();
        }
	if(evt.id().event() %1000 ==0 ){
	  cout << "Event "<<evt.id().event()<<" ReadCaloMatching done..."<<endl;
	}
}



}



using mu2e::ReadCaloMatching;
DEFINE_ART_MODULE(ReadCaloMatching);

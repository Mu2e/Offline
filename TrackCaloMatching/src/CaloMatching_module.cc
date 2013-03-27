//
//
//
// $Id: CaloMatching_module.cc,v 1.11 2013/03/27 18:37:03 murat Exp $
// $Author: murat $
// $Date: 2013/03/27 18:37:03 $
//
// Original author G. Pezzullo
//

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"

// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"

#include "KalmanTests/inc/KalRepCollection.hh"
#include "KalmanTests/inc/TrkFitDirection.hh"

#include "TrackCaloMatching/inc/TrkToCaloExtrapolCollection.hh"
#include "TrackCaloMatching/inc/TrackClusterLink.hh"

#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "CaloCluster/inc/CaloClusterTools.hh"
#include "CaloCluster/inc/CaloClusterUtilities.hh"

// Other includes.
#include "cetlib/exception.h"

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
#include "BaBar/BaBar/include/Constants.hh"


using namespace std;
using cet::square;
using cet::sum_of_squares;

namespace mu2e {

  double thetaWimpact(const CLHEP::Hep3Vector& mom, int vaneId){
    GeomHandle<VaneCalorimeter> cg;
    Vane const &vane = cg->vane(vaneId);
    CLHEP::Hep3Vector dirMom_rotated = (vane.rotation())*mom.unit();

    if(std::fabs(dirMom_rotated.getX() ) < 1e-10){
      if(dirMom_rotated.getZ()>0.0) {
	return 90.0;
      }else{
	return 0.0;
      }
    }
    double thW = 0.0;
    thW = std::atan(-1.0*dirMom_rotated.getZ() / dirMom_rotated.getX() ) ;
    thW *= Constants::radToDegrees;
    return thW;
  }
  double thetaVimpact(const CLHEP::Hep3Vector& mom, int vaneId){//(FIXME)

    GeomHandle<VaneCalorimeter> cg;
    Vane const &vane = cg->vane(vaneId);
    CLHEP::Hep3Vector dirMom_rotated = (vane.rotation())*mom.unit();

    if(std::fabs(dirMom_rotated.getX() ) < 1e-10){
      if(dirMom_rotated.getY()>0.0) {
	return 90.0;
      }else{
	return 0.0;
      }
    }
    double thV = 0.0;
    thV = std::atan(-1.0*dirMom_rotated.getY() / dirMom_rotated.getX() ) ;
    thV *= Constants::radToDegrees;
    return thV;
  }

CLHEP::Hep3Vector fromTrkToMu2eFrame(CLHEP::Hep3Vector  &vec){
        art::ServiceHandle<GeometryService> geom;
        double solenoidOffSetX = geom->config().getDouble("mu2e.solenoidOffset");
        double solenoidOffSetZ = -geom->config().getDouble("mu2e.detectorSystemZ0");
        CLHEP::Hep3Vector res;

        res.setX(vec.x() - solenoidOffSetX);
        res.setZ(vec.z() - solenoidOffSetZ);
        res.setY(vec.y());
        return res;
}


class CaloMatching : public art::EDProducer {
public:
        explicit CaloMatching(fhicl::ParameterSet const& pset):
        _fitterModuleLabel(pset.get<string>("fitterModuleLabel")),
        _diagLevel(pset.get<int>("diagLevel",0)),
        _outPutNtup(pset.get<int>("outPutNtup",0)),
        _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel", "makeCaloCluster")),
        _caloClusterAlgorithm(pset.get<std::string>("caloClusterAlgorithm", "closest")),
        _caloClusterSeeding(pset.get<std::string>("caloClusterSeeding", "energy")),
        _producerName("Algo"+mu2e::TOUpper(_caloClusterAlgorithm)+"SeededBy"+mu2e::TOUpper(_caloClusterSeeding)),
        _trkToCaloExtrapolModuleLabel(pset.get<std::string>("trkToCaloExtrapolModuleLabel", "TrkExtrapol")),
        _Ntup(0),
        _application(0),
        _directory(0),
        _firstEvent(true){
                // Tell the framework what we make.
                produces<TrackClusterLink>();
        }
        virtual ~CaloMatching() {
        }

  void beginJob();

  void endJob() {}

        void produce(art::Event & e );

private:

        double chiSquare(double& exV, double& exW, double& clV, double& clW, double& sigmaV2,
                        double& sigmaW2, double& extrT , double& clT, double& sigmaT2,
                        double& exEnergy, double& clEnergy, double& sigmaE2, int index);

        void doMatching(art::Event & evt, bool skip);
        std::string _fitterModuleLabel;

        // Diagnostic level
        int _diagLevel;

        //Ntupla for detailed information about the matching
        int _outPutNtup;

        // Label of the calo clusters  maker
        std::string _caloClusterModuleLabel;

        string _caloClusterAlgorithm;
        string _caloClusterSeeding;
        string _producerName;


        // Label of the extrapolated impact points
        std::string _trkToCaloExtrapolModuleLabel;

        std::auto_ptr<MCCaloUtilities> CaloManager;

        bool _skipEvent;

        TTree* _Ntup;

        Int_t _recoNumber,
        _recoIndex,
        _clCandidates,
        _clSize[10000],
        _clVane;

        Float_t _evt,
        _clE[10000] ,
        _clEErr[10000],
        _clT[10000] ,
        _clTErr[10000],
        _clCOGu[10000] ,
        _clCOGv[10000] ,
        _clCOGw[10000] ,
        _clCOGuErr[10000] ,
        _clCOGvErr[10000] ,
        _clCOGwErr[10000] ,
        _clCOGrow[10000],
        _clCOGcolumn[10000],
        _clShowerDir[10000],
        _clErrShowerDir[10000],
        _clCryEnergyMaxRow[10000],
        _clCryEnergyMaxColumn[10000],
        _clWsize[10000],
        _clVsize[10000],

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
        _pseudoChiSquare[10000],
        _timeChiSquare[10000],
        _energyChiSquare[10000],
        _posVChiSquare[10000],
        _posWChiSquare[10000];


        // The job needs exactly one instance of TApplication.  See note 1.
        auto_ptr<TApplication> _application;

        // Save directory from beginJob so that we can go there in endJob. See note 3.
        TDirectory* _directory;
        bool _firstEvent;


};

double CaloMatching::chiSquare(double& exV, double& exW, double& clV, double& clW, double& sigmaV2,
                double& sigmaW2, double& extrT , double& clT, double& sigmaT2,
                double& exEnergy, double& clEnergy, double& sigmaE2, int index){

        double res = 0.0;
        _posVChiSquare[index] = pow(exV -clV, 2)/sigmaV2;
        _posWChiSquare[index] = pow(exW -clW, 2)/sigmaW2;
        _timeChiSquare[index] = pow(extrT -clT, 2)/sigmaT2;
        _energyChiSquare[index] = pow(exEnergy -clEnergy, 2)/sigmaE2;

        res = _posVChiSquare[index] + _posWChiSquare[index] + _energyChiSquare[index] + _timeChiSquare[index];
        return res;
}

bool passCondition(const CaloCluster&clu, const KalRep *  const &trk){
        return true;
}




  void CaloMatching::beginJob() {

    if (_outPutNtup == 1) {

        art::ServiceHandle<art::TFileService> tfs;
        _Ntup        = tfs->make<TTree>("trkCalo", "trk-calo matching info");

        _Ntup->Branch("evt", &_evt , "evt/F");
        _Ntup->Branch("recoNumber",     &_recoNumber , "recoNumber/I");
        _Ntup->Branch("clCandidates",     &_clCandidates , "clCandidates/I");

        _Ntup->Branch("clSize[clCandidates]", _clSize , "clSize[clCandidates]/I");
        _Ntup->Branch("clVane", &_clVane , "clVane/I");
        _Ntup->Branch("clE[clCandidates]"   , _clE    , "clE[clCandidates]/F");
        _Ntup->Branch("clEErr[clCandidates]"   , _clEErr    , "clEErr[clCandidates]/F");
        _Ntup->Branch("clT[clCandidates]"   , _clT   , "clT[clCandidates]/F");
        _Ntup->Branch("clTErr[clCandidates]"   , _clTErr   , "clTErr[clCandidates]/F");
        _Ntup->Branch("clCOGu[clCandidates]"   , _clCOGu   , "clCOGu[clCandidates]/F");
        _Ntup->Branch("clCOGv[clCandidates]"   , _clCOGv   , "clCOGv[clCandidates]/F");
        _Ntup->Branch("clCOGw[clCandidates]"   , _clCOGw   , "clCOGw[clCandidates]/F");
        _Ntup->Branch("clCOGuErr[clCandidates]"   , _clCOGuErr   , "clCOGuErr[clCandidates]/F");
        _Ntup->Branch("clCOGvErr[clCandidates]"   , _clCOGvErr   , "clCOGvErr[clCandidates]/F");
        _Ntup->Branch("clCOGwErr[clCandidates]"   , _clCOGwErr   , "clCOGwErr[clCandidates]/F");//
        _Ntup->Branch("clCOGrow[clCandidates]"   ,_clCOGrow   , "clCOGrow[clCandidates]/F");
        _Ntup->Branch("clCOGcolumn[clCandidates]"   , _clCOGcolumn   , "clCOGcolumn[clCandidates]/F");
        _Ntup->Branch("clShowerDir[clCandidates]"   , _clShowerDir   , "clShowerDir[clCandidates]/F");
        _Ntup->Branch("clErrShowerDir[clCandidates]"   ,_clErrShowerDir   , "clErrShowerDir[clCandidates]/F");
        _Ntup->Branch("clCryEnergyMaxRow[clCandidates]"   , _clCryEnergyMaxRow   , "clCryEnergyMaxRow[clCandidates]/F");
        _Ntup->Branch("clCryEnergyMaxColumn[clCandidates]"   , _clCryEnergyMaxColumn   , "clCryEnergyMaxColumn[clCandidates]/F");
        _Ntup->Branch("clWsize[clCandidates]"   , _clWsize   , "clWsize[clCandidates]/F");
        _Ntup->Branch("clVsize[clCandidates]"   , _clVsize   , "clVsize[clCandidates]/F");


        _Ntup->Branch("recoIndex",     &_recoIndex , "recoIndex/I");
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

        _Ntup->Branch("pseudoChiSquare[clCandidates]",  _pseudoChiSquare , "pseudoChiSquare[clCandidates]/F");
        _Ntup->Branch("timeChiSquare[clCandidates]", _timeChiSquare , "timeChiSquare[clCandidates]/F");
        _Ntup->Branch("energyChiSquare[clCandidates]", _energyChiSquare , "energyChiSquare[clCandidates]/F");
        _Ntup->Branch("posVChiSquare[clCandidates]", _posVChiSquare , "posVChiSquare[clCandidates]/F");
        _Ntup->Branch("posWChiSquare[clCandidates]", _posWChiSquare , "posWChiSquare[clCandidates]/F");
    }
  }

void CaloMatching::produce(art::Event & evt ) {

        doMatching(evt, _skipEvent);
}


//-----------------------------------------------------------------------------
void CaloMatching::doMatching(art::Event & evt, bool skip){
  const char* oname = "CaloMatching::doMatching";
  int     nclusters, ntracks, nex, nmatches;

  //Get handle to calorimeter
  art::ServiceHandle<GeometryService> geom;
  if(! geom->hasElement<VaneCalorimeter>() ) return;
  GeomHandle<VaneCalorimeter> cg;
  
  art::Handle<KalRepCollection> trksHandle;
  evt.getByLabel(_fitterModuleLabel,"DownstreameMinus",trksHandle);
  KalRepCollection const& trks = *trksHandle;
  ntracks = trks.size();
  
  if(_diagLevel>2){
    cout<<endl<<"Event Number : "<< evt.event()<< endl;
    cout<<"\n start CaloMatching..."<<endl;
    cout<<"trks.size() = "<< trks.size() <<endl;
  }


  art::Handle<CaloClusterCollection> caloClusters;
  evt.getByLabel(_caloClusterModuleLabel,_producerName, caloClusters );
  nclusters = caloClusters->size();
  
  art::Handle<TrkToCaloExtrapolCollection>  trjExtrapols;
  evt.getByLabel(_trkToCaloExtrapolModuleLabel, trjExtrapols);
  nex = trjExtrapols->size();

  std::auto_ptr<TrackClusterLink> trackClusterLink(new TrackClusterLink);
  
  ClusterMapVector clusterMapVector;

  if (_diagLevel > 2) {
    cout<<"CaloClusters.size() = "<< nclusters << endl;
  }
//-----------------------------------------------------------------------------
// why this copy is necessary ? - doesn't seem to be
//-----------------------------------------------------------------------------
  for(size_t icl=0; icl<caloClusters->size(); ++icl) {
    CaloCluster  clu = (*caloClusters).at(icl);
    
    ClusterMap tmpCluMap(clu);

    clusterMapVector.push_back(tmpCluMap);
  }

  if(_diagLevel > 10){
    clusterMapVector.print(cout);
  }
//-----------------------------------------------------------------------------
// define more variables
//-----------------------------------------------------------------------------
  int tmpVane = -1, tmpCluCryMaxEdepColumn = -1, tmpCluCryMaxEdepRow = -1;
  int tmpWsize = -1;
  //  MatchingDataVector matchingDataVec;
  double thV = 0.0, thW = 0.0;

  double tmpCOGv(0.), tmpCOGw(0.), tmpCOGvErr(0.), tmpCOGwErr(0.), tmpCluTime(0.), tmpCluTimeErr(0.), tmpCluEnergy(0.);
  double tmpPv  (0.), tmpPw  (0.), tmpPvErr  (0.), tmpPwErr  (0.), tmpPtime  (0.), tmpPtimeErr  (0.);
  CLHEP::Hep3Vector cogVaneFrame, tmpPosVaneFrame;
  CLHEP::Hep3Vector cogVaneFrameErr, tmpPosVaneFrameErr;
  double sigmaV2 = 0.0, sigmaW2 = 0.0, sigmaT2 = 0.0, sigmaE2 = 0.0;
  double chiQ = 0.0;
  double tmpEnergy = 0.0, tmpEnergyErr = 0.0;
  double tmpErrPos = 0.0, tmpCluEnergyErr =0.0;
  Hep3Vector momdir;
  int        count = 0;
  HepVector  momvec(3);
  //define the local axes V and W for the local vane frame
  Hep3Vector Vaxes(0.0, 1.0, 0.0), Waxes(0.0, 0.0, 1.0);
  double           scaleErrV, scaleErrW, x__tmp;
  const Vane      *vane;

  
  const TrkToCaloExtrapol  *extrk;
  const KalRep             *krep, *krep0;
  const CaloCluster        *cl;
  double                   chi2_max(1.e12);
  int                      iex, icl, ltrk;
  CLHEP::Hep3Vector        tmpV;
//-----------------------------------------------------------------------------
// here the execution starts
//-----------------------------------------------------------------------------
  if(_diagLevel > 2){
    cout<<"trjExtrapols loop starts..."<< ", trjExtrapols.size() = "<< trjExtrapols->size()<<endl;
  }
//-----------------------------------------------------------------------------
// 'krep0' is used to calculate index of a track by its pointer
//-----------------------------------------------------------------------------
  if (ntracks == 0)                                         goto END;
  krep0   = trks.at(0);

  double chi2_best[100][4];
  int    iex_best [100][4];
  int    icl_best [100][4];

  for (int it=0; it<ntracks; it++) {
    for (int iv=0; iv<4; iv++) {
      chi2_best[it][iv] = chi2_max+1;
      iex_best [it][iv] = -1;
      icl_best [it][iv] = -1;
    }
  }

  for (size_t iex=0; iex<nex; ++iex){
    extrk = &trjExtrapols->at(iex);
    krep  = *extrk->trk();
//-----------------------------------------------------------------------------
// track index, important: store one, the best, intersection per track per vane
//-----------------------------------------------------------------------------
    ltrk  = krep-krep0;
    if (ltrk > 100) {
      printf("%s ERROR: more than 100 tracks, skip the rest\n", oname);
                                                        goto NEXT_INTERSECTION;
    }
//-----------------------------------------------------------------------------
// apparently, ntupling stuff mixed in
//-----------------------------------------------------------------------------
    _recoIndex  = 0;
    _recoNumber = 0;
    _evt =  0.0;
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
    
    if(_diagLevel > 2) {
      cout<<"trjExtrapols at "<<iex<<" step"<<endl;
    }
    
    if(iex>0){
      KalRepPtr const& tmpTrkPtr = trjExtrapols->at(iex-1).trk();
      const KalRep *  const &tmpTrk = *tmpTrkPtr;
      if(krep == tmpTrk){
	++count;
      }else{
	count = 0;
      }
    }

    _evt = evt.id().event();
    _recoNumber = iex;
    _recoIndex  = count;
    tmpVane     = extrk->vaneId();
    tmpEnergy   = extrk->momentum().mag();

    _clVane     = tmpVane;
    _recoE      = tmpEnergy;

    momdir      = extrk->momentum().unit();
    
    for(int icor=0;icor<3;icor++){
      momvec[icor] = momdir[icor];
    }
    tmpEnergyErr = sqrt(extrk->momentumErr().covMatrix().similarity(momvec));
    
    _recoEErr   = tmpEnergyErr;
    
    if (!clusterMapVector.searchVane(tmpVane) ) continue;
    
    thV = thetaVimpact(extrk->momentum(), extrk->vaneId());
    thW = thetaWimpact(extrk->momentum(), extrk->vaneId());
    
    _recoThetaV = thV;
    _recoThetaW = thW;
    
    tmpPtime     = extrk->time();
    tmpPtimeErr  = krep->t0().t0Err();
    _recoTime    = tmpPtime;
    _recoTimeErr = tmpPtimeErr;

    if(_diagLevel > 2){
      cout<<"tmpVane = "<<tmpVane<<endl;
      cout<<"tmpEnergy = "<<tmpEnergy<<" [MeV]"<<
	", tmpEnergyErr = "<< tmpEnergyErr <<" [MeV]"<<endl;
      cout<<"thetaVimpact = "<< thV<<" [deg]"<<endl;
      cout<<"tmpPtime = "<<tmpPtime<<" [ns]"<<
	"tmpPtimeErr = "<< tmpPtimeErr<<" [ns]"<<endl;
    }
	  
    tmpPosVaneFrame.setX(extrk->entrancePosition().x());
    tmpPosVaneFrame.setY(extrk->entrancePosition().y());
    tmpPosVaneFrame.setZ(extrk->entrancePosition().z());

    if (_diagLevel > 2){
      cout<<"Tracker frame: tmpPosVaneFrame = "<<tmpPosVaneFrame<<" [mm]"<<
	"tmpPosVaneFrameErr = "<< tmpPosVaneFrameErr<<" [mm]"<<endl;
    }

    tmpPosVaneFrame = fromTrkToMu2eFrame(tmpPosVaneFrame);
	  
    if(_diagLevel > 2){
      cout<<"Mu2e general frame: tmpPosVaneFrame = "<<tmpPosVaneFrame<<" [mm]"<<endl;
    }
    _recoPx = tmpPosVaneFrame.x();
    _recoPy = tmpPosVaneFrame.y();
    _recoPz = tmpPosVaneFrame.z();
    
    //move these axes into the Mu2e general frame
    vane = &cg->vane(tmpVane );
    
    Vaxes = vane->rotation()*(Vaxes);
    Waxes = vane->rotation()*(Waxes);
	  
    //Vaxes = cg->fromVaneFrame(tmpVane, Vaxes);
    //Waxes = cg->fromVaneFrame(tmpVane, Waxes);
	  
    thV /= Constants::radToDegrees;
    thW /= Constants::radToDegrees;
    scaleErrW = 1.0/fabs( cos(thW) );
    scaleErrV = 1.0/fabs( cos(thV) );
	  
    momvec[0] = Vaxes.x();
    momvec[1] = Vaxes.y();
    momvec[2] = Vaxes.z();
	  
    tmpErrPos = sqrt(extrk->entrancePositionErr().covMatrix().similarity(momvec));

    if(_diagLevel > 2){
      cout<<"before scaling V-pos error..."<< "posV_Err = "<< tmpErrPos<<" [mm]"<<endl;
    }
    _recoPvErr = tmpErrPos*scaleErrV;

    if(_diagLevel > 2){
      cout<<"after scaling V-pos error..."<< "posV_Err = "<< tmpErrPos<<" [mm]"<<endl;
    }

    momvec[0] = Waxes.x();
    momvec[1] = Waxes.y();
    momvec[2] = Waxes.z();
    
    tmpErrPos = sqrt(extrk->entrancePositionErr().covMatrix().similarity(momvec));

    if(_diagLevel > 2){
      cout<<"before scaling W-pos error..."<< "posW_Err = "<< tmpErrPos<<" [mm]"<<endl;
    }
    _recoPwErr = tmpErrPos*scaleErrW;

    if (_diagLevel > 2){
      cout<<"before scaling W-pos error..."<<"posW_Err = "<< tmpErrPos<<" [mm]"<<endl;
    }
    
    tmpV = cg->toVaneFrame(tmpVane, tmpPosVaneFrame);
	  
    tmpPv           = tmpV.y();
    _recoPv         = tmpPv;
    tmpPvErr        = _recoPvErr;
    tmpPw           = tmpV.z();
    _recoPw         = tmpPw;
    tmpPwErr        = _recoPwErr;
	  
    if(_diagLevel > 2){
      cout<<"tmpPosVaneFrame = "<<tmpV<<" [mm]"<<
	"tmpPosVaneFrameErr = "<< tmpPosVaneFrameErr<<" [mm]"<<endl;
    }

    _clCandidates = nclusters;
    
    //    chi2_best = chi2_max;
    for (int icl=0; icl<nclusters; icl++) {
//-----------------------------------------------------------------------------
// why it is 'indexVecCluster[j]' and not just 'j' ?
//-----------------------------------------------------------------------------
      cl  = &(*caloClusters).at(icl);
      CaloClusterTools cluTool(*cl);
      if (cl->vaneId() != tmpVane)    goto NEXT_CLUSTER;
//-----------------------------------------------------------------------------
// 2013-03-24 P.Murat: cluster center of gravity is determined in the GLOBAL 
//                     coordinate system, 
// clustering needs to be fixed to define cluster coordinates in the local frame system
// for now - transform back to the local coordinate system
//-----------------------------------------------------------------------------
      // cogVaneFrame                    = clu.cog3Vector();
      cogVaneFrame = cg->toVaneFrame(tmpVane, cl->cog3Vector());
      
      tmpCluCryMaxEdepRow      = cluTool.cryEnergydepMaxRow();
      tmpCluCryMaxEdepColumn   = cluTool.cryEnergydepMaxColumn();
      tmpCluTime               = cl->time();
      tmpCOGv                  = cogVaneFrame.y();
      tmpCOGw                  = cogVaneFrame.z();
      tmpCOGvErr               = cl->cog3VectorError().y();
      tmpCOGwErr               = cl->cog3VectorError().z();
      tmpWsize                 = cluTool.wSize();
      tmpCluEnergy             = cl->energyDep();
      tmpCluEnergyErr          = cluTool.energyDepErr();//clusterEnergyErr(tmpCluEnergy);
      tmpCluTimeErr            = cluTool.timeErr();
	    
      _clE[icl]                  = tmpCluEnergy;
      _clEErr[icl]               = tmpCluEnergyErr;
      _clT[icl]                  = tmpCluTime;
      _clTErr[icl]               = tmpCluTimeErr;
      
      _clSize[icl]               = cl->size();
      
      _clCOGrow[icl]             =  cl->cogRow();
      _clCOGcolumn[icl]          =  cl->cogColumn();
      _clShowerDir[icl]          =  cluTool.showerDir();
      _clErrShowerDir[icl]       =  cluTool.errShowerDir();
      _clCryEnergyMaxRow[icl]    =  tmpCluCryMaxEdepRow;
      _clCryEnergyMaxColumn[icl] =  tmpCluCryMaxEdepColumn;
      
      _clWsize[icl]              =  cluTool.wSize();
      _clVsize[icl]              =  cluTool.vSize();
      
      if(_diagLevel > 2){
	cout<<" tmpPv = "<< tmpPv <<
	  ", tmpCOGv = "<<  tmpCOGv<<
	  ", sigmaV2 = "<<sigmaV2<<endl<<
	  ", tmpPw = "<<tmpPw<<
	  ", tmpCOGw = "<<tmpCOGw<<
	  ", sigmaW2 = "<<tmpCOGv<<endl<<
	  ", tmpPtime = "<<tmpPtime<<
	  ", tmpCluTime = "<<tmpCluTime<<
	  ", sigmaT2 = "<<sigmaT2<<
	  ", tmpWsize = "<<tmpWsize<<endl;
      }
//-----------------------------------------------------------------------------
// 2013-03-24 : the first w-correction below is doing smth wrong.. moving 
//              the cluster in e0000001/event#1 by about 60cm.. 
//              for now, turn both W-corrections off, leaving Z unchanged
//-----------------------------------------------------------------------------
      v_correction_0(thV, tmpCOGv, tmpCOGvErr);

//                         w_correction_0(tmpCOGw, tmpCOGwErr, tmpCluCryMaxEdepColumn);
//                         w_correction_1(tmpCOGw,tmpCOGwErr, tmpWsize);

      x__tmp = tmpCOGw;
      w_correction_0(x__tmp, tmpCOGwErr, tmpCluCryMaxEdepColumn);
      w_correction_1(x__tmp,tmpCOGwErr, tmpWsize);

      tmpCOGvErr = 15.35 / 2.35;
      tmpCOGwErr = 15.0 / 2.35;
//-----------------------------------------------------------------------------
// 2013-03-13 P.Murat: use 20 MeV as the mathching resolution
//-----------------------------------------------------------------------------
      sigmaE2 = 20.; // cet::sum_of_squares(tmpEnergyErr, tmpCluEnergyErr);
      sigmaV2 = cet::sum_of_squares(tmpPvErr, tmpCOGvErr);
      sigmaW2 = cet::sum_of_squares(tmpPwErr, tmpCOGwErr);
      sigmaT2 = cet::sum_of_squares(tmpPtimeErr, tmpCluTimeErr);

      if(_diagLevel > 2){
	cout<<"after doing many corrections..."<<endl;
	cout<<" tmpPv = "<< tmpPv <<
	  ", tmpCOGv = "<<tmpCOGv<<
	  ", sigmaV2 = "<<sigmaV2<<endl<<
	  ", tmpPw = "<<tmpPw<<
	  ", tmpCOGw = "<<tmpCOGw<<
	  ", sigmaW2 = "<<sigmaW2<<endl<<
	  ", tmpPtime = "<<tmpPtime<<
	  ", tmpCluTime = "<<tmpCluTime<<
	  ", sigmaT2 = "<<sigmaT2<<endl;
      }
		  
      _clCOGu   [icl] = cogVaneFrame.x();
      _clCOGv   [icl] = tmpCOGv;
      _clCOGw   [icl] = tmpCOGw;
      _clCOGuErr[icl] = 0.0;
      _clCOGvErr[icl] = tmpCOGvErr;
      _clCOGwErr[icl] = tmpCOGwErr;

      if(_diagLevel > 2){
	cout<<"after doing other corrections..."<<endl;
	cout<<" tmpPv = "<< tmpPv <<
	  ", tmpCOGv = "<<tmpCOGv<<
	  ", sigmaV2 = "<<sigmaV2<<endl<<
	  ", tmpPw = "<<tmpPw<<
	  ", tmpCOGw = "<<tmpCOGw<<
	  ", sigmaW2 = "<<sigmaW2<<endl<<
	  ", tmpPtime = "<<tmpPtime<<
	  ", tmpCluTime = "<<tmpCluTime<<
	  ", sigmaT2 = "<<sigmaT2<<endl;
      }
		  
      chiQ = chiSquare(tmpPv, tmpPw, 
		       tmpCOGv, tmpCOGw, 
		       sigmaV2, sigmaW2, 
		       tmpPtime, tmpCluTime, sigmaT2, 
		       tmpEnergy, tmpCluEnergy, sigmaE2, icl);
	    
      if (chiQ < chi2_best[ltrk][tmpVane]) {
//-----------------------------------------------------------------------------
// new best  match
//-----------------------------------------------------------------------------
	chi2_best[ltrk][tmpVane] = chiQ;
	icl_best [ltrk][tmpVane] = icl;   // 
	iex_best [ltrk][tmpVane] = iex;
      }
    NEXT_CLUSTER:;
    }
//-----------------------------------------------------------------------------
// fill ntuple
//-----------------------------------------------------------------------------
    if (_outPutNtup > 0) {
      _Ntup->Fill();
    }
  NEXT_INTERSECTION:;
  }

//-----------------------------------------------------------------------------
// form output list of matches
//-----------------------------------------------------------------------------
  nmatches = 0;
  for (int it=0; it<ntracks; it++) {
    for (int iv=0; iv<4; iv++) {
      if (chi2_best[it][iv] < chi2_max) {
	iex = iex_best [it][iv];
	icl = icl_best [it][iv];
	trackClusterLink->addSingle(art::Ptr<TrkToCaloExtrapol>(trjExtrapols,iex), 
				    art::Ptr<CaloCluster>      (caloClusters,icl));
	nmatches += 1;
      }
    }
  }

 END:;
  evt.put(trackClusterLink);
  cout << "Event "<<evt.id().event()<<" CaloMatching done..."<<endl;

}

}

using mu2e::CaloMatching;
DEFINE_ART_MODULE(CaloMatching);

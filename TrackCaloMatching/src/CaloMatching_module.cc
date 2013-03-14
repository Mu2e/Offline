//
//
//
// $Id: CaloMatching_module.cc,v 1.9 2013/03/14 19:47:46 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/14 19:47:46 $
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


  struct matchingData{
    size_t _indexTrkToCaloExtrapol;
    std::map<size_t, double> _candidateClusterVec;

    matchingData(){};
    matchingData(size_t &index){
      _indexTrkToCaloExtrapol = index;
    };

    matchingData(size_t &indexTrk, size_t &cluster, double &value){
      _indexTrkToCaloExtrapol = indexTrk;
      _candidateClusterVec[cluster] = value;
    };

    size_t indexCluster(double &value){
      if(_candidateClusterVec.size()==0) return 0;
      size_t res = 0;
      bool trovato = false;
      std::map<size_t, double>::iterator it=_candidateClusterVec.begin();
      while(it != _candidateClusterVec.end() && !trovato){
	//for(std::map<size_t, double>::iterator it=_candidateClusterVec.begin(); it !=_candidateClusterVec.end(); ++it){
	if(it->second == value){
	  res = it->first;
	  trovato = true;
	}else{
	  ++it;
	}
      }
      if(!trovato){
	std::cout<<"index of the cluster with minimization value = "<< value <<
	  " NOT FOUND"<<endl;
      }
      return res;
    }



    size_t indexClusterMin(){
      if(_candidateClusterVec.size()==0) return 0;
      double min = _candidateClusterVec[0];
      size_t res = 0;
      for(std::map<size_t, double>::iterator it=_candidateClusterVec.begin(); it !=_candidateClusterVec.end(); ++it){
	if(it->second < min){
	  min = it->second;
	  res = it->first;
	}
      }
      return res;
    }

    double lowerMinFunc(){
      if(_candidateClusterVec.size()==0) return 0;
      double res = _candidateClusterVec[0];
      for(std::map<size_t, double>::iterator it=_candidateClusterVec.begin(); it !=_candidateClusterVec.end(); ++it){
	if(it->second < res) {
	  res = it->second;
	}
      }
      return res;
    }

    bool eraseCluster(size_t &index){
      std::map<size_t, double>::iterator it=_candidateClusterVec.find(index);
      if(it != _candidateClusterVec.end()){
	_candidateClusterVec.erase(it);
	return true;
      }
      return false;
    }


    size_t size(){
      return _candidateClusterVec.size();
    }

    //Accessors
    double &operator[](size_t n){
      return (this->_candidateClusterVec[n]);
    }

    matchingData & operator=(const matchingData other){
      _indexTrkToCaloExtrapol = other._indexTrkToCaloExtrapol;
      _candidateClusterVec = other._candidateClusterVec;
      return *this;
    }

    void print(std::ostream& os) const{
      os<<"indexTrkToCaloExtrapol = "<< this->_indexTrkToCaloExtrapol<<std::endl;
      for(std::map<size_t, double>::const_iterator it=_candidateClusterVec.begin(); it !=_candidateClusterVec.end(); ++it){
	os << "cluster index = " << it->first <<
	  ", value of the minimized function = "<< it->second<<std::endl;
      }
    }

    void addCluster(size_t &cluster, double& value){
      if(_candidateClusterVec.find(cluster) != _candidateClusterVec.end()){
	std::cout<<"This cluster is already inserted"<<endl;
	this->print(std::cout);
      }else{
	_candidateClusterVec[cluster] = value;
      }
    }

    void addCluster(unsigned int &cluster, double& value){
      if(_candidateClusterVec.find((size_t)cluster) != _candidateClusterVec.end()){
	std::cout<<"This cluster is already inserted"<<endl;
	this->print(std::cout);
      }else{
	_candidateClusterVec[cluster] = value;
      }
    }
  };

  struct MatchingDataVector{
    std::vector<matchingData> _vec;

    void push_back(matchingData &id){
      _vec.push_back(id);
    }
    size_t size(){
      return _vec.size();
    }


    bool removeCluster(size_t &clusterIndex){
      bool res = false;
      for(std::vector<matchingData>::iterator it = _vec.begin(); it != _vec.end(); ++it){
	if(it->eraseCluster(clusterIndex)){
	  _vec.erase(it);
	  res = true;
	  if(_vec.size()==0) break;
	  it = _vec.begin();

	}
      }
      return res;
    }

    bool removeTrkToCaloExtrapol(size_t &index){
      bool res = false;
      for(std::vector<matchingData>::iterator it = _vec.begin(); it != _vec.end(); ++it){
	if(it->_indexTrkToCaloExtrapol == index){
	  _vec.erase(it);
	  res = true;
	  if(_vec.size()==0) break;
	  it = _vec.begin();
	}
      }
      return res;
    }

    double minValue(){
      if(_vec.size()==0) return 0;
      double min = _vec[0].lowerMinFunc();
      for(size_t i=0; i<_vec.size(); i++){
	if(_vec[i].lowerMinFunc() < min){
	  min = _vec[i].lowerMinFunc();
	}
      }
      return min;
    }

    std::vector<size_t> findIndexMinValue(){
      std::vector<size_t> indexVec;
      if(_vec.size()==0) {
	return indexVec;
      }
      double min = _vec[0].lowerMinFunc();


      for(size_t i=0; i<_vec.size(); ++i){
	if(_vec[i].lowerMinFunc() < min){
	  min = _vec[i].lowerMinFunc();
	}
      }
      for(size_t i=0; i<_vec.size(); ++i){
	if(_vec[i].lowerMinFunc() == min){

	  indexVec.push_back(i);
	}
      }
      return indexVec;
    }

    bool hasSense(){
      if(_vec.size()==0){
	return false;
      }else {
	return true;
      }
    }

    //Accessors
    matchingData &operator[](size_t n){ return (this->_vec.at(n));}
    matchingData &at(size_t &n){ return (this->_vec.at(n));}


    matchingData &operator[](unsigned int n){ return (this->_vec.at(n));}
    matchingData &at(unsigned int &n){ return (this->_vec.at(n));}

    size_t idClusterMin(size_t &n){return _vec[n].indexClusterMin();}
    double minFuncValue(size_t &n){return _vec.at(n).lowerMinFunc();}

    size_t idClusterMin(unsigned int &n){return _vec[n].indexClusterMin();}
    double minFuncValue(unsigned int &n){return _vec.at(n).lowerMinFunc();}

    MatchingDataVector & operator=(const MatchingDataVector other){
      for(size_t i =0; i<other._vec.size(); ++i){
	_vec.push_back( other._vec.at(i) );
      }
      return *this;
    }

    MatchingDataVector(){};




  };

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


  static int ncalls(0);

  class CaloMatching : public art::EDProducer {
  public:
    explicit CaloMatching(fhicl::ParameterSet const& pset):
      _fitterModuleLabel(pset.get<string>("fitterModuleLabel")),
      _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
      _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
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
      _firstEvent(true),
      _trkdiag(0){
      // Tell the framework what we make.
      produces<TrackClusterLink>();

      // construct the data product instance name
      _iname = _fdir.name() + _tpart.name();
    }
    virtual ~CaloMatching() {
    }
    void beginJob() {}
    void endJob() {}

    void produce(art::Event & e );

  private:

    double chiSquare(double& exV, double& exW, double& clV, double& clW, double& sigmaV2,
		     double& sigmaW2, double& extrT , double& clT, double& sigmaT2,
		     double& exEnergy, double& clEnergy, double& sigmaE2, size_t& index);

    void doMatching(art::Event & evt, bool skip);
    std::string _fitterModuleLabel;

    TrkParticle _tpart;
        
    TrkFitDirection _fdir;
        
    std::string _iname;

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

       

    Int_t _sameVane,
      _recoNumber,
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
    TTree* _trkdiag;


  };

  double CaloMatching::chiSquare(double& exV, double& exW, double& clV, double& clW, double& sigmaV2,
				 double& sigmaW2, double& extrT , double& clT, double& sigmaT2,
				 double& exEnergy, double& clEnergy, double& sigmaE2, size_t& index){

    double res = 0.0;
    _posVChiSquare[index] = pow(2.35*(exV -clV), 2)/sigmaV2;
    _posWChiSquare[index] = pow(2.35*(exW -clW), 2)/sigmaW2;
    _timeChiSquare[index] = pow(2.35*(extrT -clT), 2)/sigmaT2;
    _energyChiSquare[index] = pow(2.35*(exEnergy -clEnergy), 2)/sigmaE2;

    res = _posVChiSquare[index] + _posWChiSquare[index] /*+ _energyChiSquare[index]*/ + _timeChiSquare[index];
    return res;
  }

  bool passCondition(const CaloCluster&clu, const KalRep *  const &trk){
    return true;
  }




  void CaloMatching::produce(art::Event & evt ) {

    ++ncalls;

    if (ncalls == 1 && _outPutNtup == 1) {

      art::ServiceHandle<art::TFileService> tfs;
      _Ntup        = tfs->make<TTree>("trkCalo", "trk-calo matching info");

      _Ntup->Branch("evt", &_evt , "evt/F");
      _Ntup->Branch("sameVane",     &_sameVane , "sameVane/I");
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
    doMatching(evt, _skipEvent);
  }

  void CaloMatching::doMatching(art::Event & evt, bool skip){

    //Get handle to calorimeter
    art::ServiceHandle<GeometryService> geom;
    if(! geom->hasElement<VaneCalorimeter>() ) return;
    GeomHandle<VaneCalorimeter> cg;

    art::Handle<KalRepCollection> trksHandle;
    evt.getByLabel(_fitterModuleLabel,_iname,trksHandle);
    KalRepCollection const& trks = *trksHandle;

    if(_diagLevel>2){
      cout<<endl<<"Event Number : "<< evt.event()<< endl
	  <<"start CaloMatching..."<<endl
	  <<"trks.size() = "<< trks.size() <<endl;
    }

    art::Handle<CaloClusterCollection> caloClusters;
    evt.getByLabel(_caloClusterModuleLabel,_producerName, caloClusters );

    art::Handle<TrkToCaloExtrapolCollection>  trjExtrapols;
    evt.getByLabel(_trkToCaloExtrapolModuleLabel, trjExtrapols);

    std::auto_ptr<TrackClusterLink> trackClusterLink(new TrackClusterLink);

    ClusterMapVector clusterMapVector;

    if(_diagLevel > 2){
      cout<<"CaloClusters.size() = "<<caloClusters->size()<<endl;
    }

    for(size_t icl=0; icl<caloClusters->size(); ++icl){
      CaloCluster  clu = (*caloClusters).at(icl);

      ClusterMap tmpCluMap(clu);

      clusterMapVector.push_back(tmpCluMap);
    }

    if(_diagLevel > 10){
      clusterMapVector.print(cout);
    }

    int tmpVane = -1, tmpCluCryMaxEdepColumn = -1, tmpCluCryMaxEdepRow = -1;
    int tmpWsize = -1;
    MatchingDataVector matchingDataVec;
    double thV = 0.0, thW = 0.0;

    double tmpCOGv = 0.0, tmpCOGw = 0.0, tmpCOGvErr = 0.0, tmpCOGwErr = 0.0, tmpCluTime = 0.0, tmpCluTimeErr =0.0, tmpCluEnergy = 0.0;
    double tmpPv = 0.0, tmpPw = 0.0, /*tmpPvErr = 0.0, tmpPwErr = 0.0,*/ tmpPtime = 0.0, tmpPtimeErr = 0.0;
    CLHEP::Hep3Vector cogVaneFrame, tmpPosVaneFrame;
    CLHEP::Hep3Vector cogVaneFrameErr, tmpPosVaneFrameErr;
    double sigmaV2 = 0.0, sigmaW2 = 0.0, sigmaT2 = 0.0, sigmaE2 = 0.0;
    double chiQ = 0.0;
    double tmpEnergy = 0.0, tmpEnergyErr = 0.0;
    double tmpErrPos = 0.0, tmpCluEnergyErr =0.0;

    if(_diagLevel > 2){
      cout<<"trjExtrapols loop starts..."<<
	", trjExtrapols.size() = "<< trjExtrapols->size()<<endl;

    }
    Hep3Vector momdir;
    int count = 0;
    HepVector momvec(3);
    for(size_t i=0; i<trjExtrapols->size(); ++i){

      _sameVane = 1;
      _recoIndex = 0;
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

      if(_diagLevel > 2){
	cout<<"trjExtrapols at "<<i<<" step"<<endl;
      }

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

      _evt = evt.id().event();
      _recoNumber = i;
      _recoIndex = count;
      tmpVane = trjExtrapols->at(i).vaneId();
      tmpEnergy = trjExtrapols->at(i).momentum().mag();

      _clVane = tmpVane;
      _recoE = tmpEnergy;

      momdir = trjExtrapols->at(i).momentum().unit();
      for(int icor=0;icor<3;icor++){
	momvec[icor] = momdir[icor];
      }
      tmpEnergyErr = sqrt(trjExtrapols->at(i).momentumErr().covMatrix().similarity(momvec));

      _recoEErr = tmpEnergyErr;

     

      thV = thetaVimpact(trjExtrapols->at(i).momentum(), trjExtrapols->at(i).vaneId());
      thW = thetaWimpact(trjExtrapols->at(i).momentum(), trjExtrapols->at(i).vaneId());

      _recoThetaV = thV;
      _recoThetaW = thW;

      tmpPtime = trjExtrapols->at(i).time();
      tmpPtimeErr = trk->t0().t0Err();

      _recoTime = tmpPtime;
      _recoTimeErr = tmpPtimeErr;
      if(_diagLevel > 2){
	cout<<"tmpVane = "<<tmpVane<<endl;
	cout<<"tmpEnergy = "<<tmpEnergy<<" [MeV]"<<
	  ", tmpEnergyErr = "<< tmpEnergyErr <<" [MeV]"<<endl;
	cout<<"thetaVimpact = "<< thV<<" [deg]"<<endl;
	cout<<"tmpPtime = "<<tmpPtime<<" [ns]"<<
	  "tmpPtimeErr = "<< tmpPtimeErr<<" [ns]"<<endl;
      }

      tmpPosVaneFrame.setX(trjExtrapols->at(i).entrancePosition().x());
      tmpPosVaneFrame.setY(trjExtrapols->at(i).entrancePosition().y());
      tmpPosVaneFrame.setZ(trjExtrapols->at(i).entrancePosition().z());
      if(_diagLevel > 2){
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
      
      //look the existence of a EMC cluster in the same vane
      if(!clusterMapVector.searchVane(tmpVane) ) {
	_sameVane=0;
	if(_outPutNtup>0) _Ntup->Fill();
	continue;
      }
      
      //define the local axes V and W for the local vane frame
      Hep3Vector Vaxes(0.0, 1.0, 0.0), Waxes(0.0, 0.0, 1.0);

      //move these axes into the Mu2e general frame
      Vane const &vane = cg->vane(tmpVane );

      Vaxes = (vane.rotation())*(Vaxes);
      Waxes = (vane.rotation())*(Waxes);

      //Vaxes = cg->fromVaneFrame(tmpVane, Vaxes);
      //Waxes = cg->fromVaneFrame(tmpVane, Waxes);

      thV /= Constants::radToDegrees;
      thW /= Constants::radToDegrees;
      double scaleErrW = 1.0/fabs( cos(thW) );
      double scaleErrV = 1.0/fabs( cos(thV) );

      momvec[0] = Vaxes.x();
      momvec[1] = Vaxes.y();
      momvec[2] = Vaxes.z();

      tmpErrPos = sqrt(trjExtrapols->at(i).entrancePositionErr().covMatrix().similarity(momvec));
      if(_diagLevel > 2){
	cout<<"before scaling V-pos error..."<<
	  "posV_Err = "<< tmpErrPos<<" [mm]"<<endl;
      }
      _recoPvErr = tmpErrPos*scaleErrV;
      if(_diagLevel > 2){
	cout<<"after scaling V-pos error..."<<
	  "posV_Err = "<< tmpErrPos<<" [mm]"<<endl;
      }
      momvec[0] = Waxes.x();
      momvec[1] = Waxes.y();
      momvec[2] = Waxes.z();

      tmpErrPos = sqrt(trjExtrapols->at(i).entrancePositionErr().covMatrix().similarity(momvec));
      if(_diagLevel > 2){
	cout<<"before scaling W-pos error..."<<
	  "posW_Err = "<< tmpErrPos<<" [mm]"<<endl;
      }
      _recoPwErr = tmpErrPos*scaleErrW;
      if(_diagLevel > 2){
	cout<<"before scaling W-pos error..."<<
	  "posW_Err = "<< tmpErrPos<<" [mm]"<<endl;
      }

      CLHEP::Hep3Vector tmpV = cg->toVaneFrame(tmpVane, tmpPosVaneFrame);

      tmpPv           = tmpV.y();
      _recoPv         = tmpPv;
      //tmpPvErr        = _recoPvErr;
      tmpPw           = tmpV.z();
      _recoPw         = tmpPw;
      //tmpPwErr        = _recoPwErr;

      if(_diagLevel > 2){
	cout<<"tmpPosVaneFrame = "<<tmpV<<" [mm]"<<
	  "tmpPosVaneFrameErr = "<< tmpPosVaneFrameErr<<" [mm]"<<endl;
      }

      trkIdVector indexVecCluster = clusterMapVector.findVane(tmpVane);
      matchingData tmpMatchData(i);

      if(_diagLevel > 2){
	cout<<" indexVecCluster loop starts..."<<
	  ", indexVecCluster.size() = "<< indexVecCluster.size()<<endl;
      }

      _clCandidates = indexVecCluster.size();
      for(size_t j=0; j<indexVecCluster.size(); ++j){
	if(_diagLevel > 2){
	  cout<<" indexVecCluster at "<<j<<" step"<<endl;
	}
	CaloCluster  clu = (*caloClusters).at(indexVecCluster[j]);
	CaloClusterTools cluTool(clu);

	cogVaneFrame                    = clu.cog3Vector();//caloClusters->at(indexVecCluster[j]).cog3Vector() ;
	cogVaneFrame                    = cg->toVaneFrame(tmpVane, cogVaneFrame);
	tmpCluCryMaxEdepRow             = cluTool.cryEnergydepMaxRow();
	tmpCluCryMaxEdepColumn          = cluTool.cryEnergydepMaxColumn();
	tmpCluTime                      = caloClusters->at(indexVecCluster[j]).time();
	tmpCOGv                         = cogVaneFrame.y();
	tmpCOGw                         = cogVaneFrame.z();
	tmpCOGvErr                      = caloClusters->at(indexVecCluster[j]).cog3VectorError().y();
	tmpCOGwErr                      = caloClusters->at(indexVecCluster[j]).cog3VectorError().z();
	tmpWsize                        = cluTool.wSize();
	tmpCluEnergy                    = caloClusters->at(indexVecCluster[j]).energyDep();
	tmpCluEnergyErr                 = cluTool.energyDepErr();//clusterEnergyErr(tmpCluEnergy);
	tmpCluTimeErr                   = cluTool.timeErr();

	_clE[j]                         = tmpCluEnergy;
	_clEErr[j]                      = tmpCluEnergyErr;
	_clT[j]                         = tmpCluTime;
	_clTErr[j]                      = tmpCluTimeErr;

	_clSize[j]                      = clu.size();

	_clCOGrow[j]                    =  clu.cogRow();
	_clCOGcolumn[j]                 =  clu.cogColumn();
	_clShowerDir[j]                 =  cluTool.showerDir();
	_clErrShowerDir[j]              =  cluTool.errShowerDir();
	_clCryEnergyMaxRow[j]           =  tmpCluCryMaxEdepRow;
	_clCryEnergyMaxColumn[j]        =  tmpCluCryMaxEdepColumn;

	_clWsize[j]                     =  cluTool.wSize();
	_clVsize[j]                     =  cluTool.vSize();


	if(_diagLevel > 2){
	  cout<<" tmpPv = "<< tmpPv <<
	    ", tmpCOGv = "<<tmpCOGv<<
	    ", sigmaV2 = "<<sigmaV2<<endl<<
	    ", tmpPw = "<<tmpPw<<
	    ", tmpCOGw = "<<tmpCOGw<<
	    ", sigmaW2 = "<<tmpCOGv<<endl<<
	    ", tmpPtime = "<<tmpPtime<<
	    ", tmpCluTime = "<<tmpCluTime<<
	    ", sigmaT2 = "<<sigmaT2<<
	    ", tmpWsize = "<<tmpWsize<<endl;
	}

	v_correction_0(thV, tmpCOGv, tmpCOGvErr);
	w_correction_0(tmpCOGw, tmpCOGwErr, tmpCluCryMaxEdepColumn);
	w_correction_1(tmpCOGw,tmpCOGwErr, tmpWsize);

	tmpCOGvErr = 15.35 / 2.35;
	tmpCOGwErr = 15.0 / 2.35;

	sigmaE2 = std::pow(7.2, 2);//cet::sum_of_squares(tmpEnergyErr, tmpCluEnergyErr);
	sigmaV2 = std::pow(44.33, 2);//cet::sum_of_squares(tmpPvErr, tmpCOGvErr);
	sigmaW2 = std::pow(22.17, 2);//cet::sum_of_squares(tmpPwErr, tmpCOGwErr);
	sigmaT2 = std::pow(3.4, 2);//cet::sum_of_squares(tmpPtimeErr, tmpCluTimeErr);

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

	_clCOGu[j]                      =  cogVaneFrame.x();
	_clCOGv[j]                      =  tmpCOGv;
	_clCOGw[j]                      =  tmpCOGw;
	_clCOGuErr[j]                   =  0.0;
	_clCOGvErr[j]                   = tmpCOGvErr;
	_clCOGwErr[j]                   = tmpCOGwErr;
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

	chiQ = chiSquare(tmpPv, tmpPw, tmpCOGv, tmpCOGw, sigmaV2, sigmaW2, tmpPtime, tmpCluTime, sigmaT2, tmpEnergy, tmpCluEnergy, sigmaE2, j);

	_pseudoChiSquare[j] = chiQ;

	if(passCondition(caloClusters->at(indexVecCluster[j]), trk)){
	  if(_diagLevel > 2){
	    cout<<"indexVecCluster[j] = "<<indexVecCluster[j]<<
	      "chiQ = "<<chiQ<<endl;

	  }
	  tmpMatchData.addCluster(indexVecCluster[j], chiQ);
	}
      }
      if(_outPutNtup>0){
	_Ntup->Fill();
      }
      if(indexVecCluster.size() > 0){
	matchingDataVec.push_back(tmpMatchData);
      }
    }

    if(_diagLevel > 2){
      cout<< "matchingDataVec : "<<endl;
      for(size_t i=0; i< matchingDataVec.size(); ++i){
	matchingDataVec[i].print(cout);
      }
    }



    while(matchingDataVec.hasSense() ){
      if(_diagLevel > 2){
	cout<<"'matchingDataVec.hasSense()'"<<endl;
      }
      std::vector<size_t> vecTracks = matchingDataVec.findIndexMinValue();
      size_t trkIndex = vecTracks[0];
      chiQ = matchingDataVec.minValue();

      size_t clusterIndex = matchingDataVec[trkIndex].indexCluster(chiQ);
      if(_diagLevel > 2){
	cout<<"vecTracks.size() = "<< vecTracks.size()<<
	  ", trkIndex = "<< trkIndex<<
	  ", clusterIndex = "<<clusterIndex<<
	  ", chiQ = "<< chiQ<<endl;
      }

      if(vecTracks.size() == 1) {
	chiQ = matchingDataVec.minValue();

	trackClusterLink->addSingle( art::Ptr<TrkToCaloExtrapol>(trjExtrapols, trkIndex),
				     art::Ptr<CaloCluster>(caloClusters,  clusterIndex ) );
	if(_diagLevel > 2){
	  cout<<"trackClusterLink->addSingle() done..."<<endl;
	}
      }else{
	chiQ = matchingDataVec.minValue();

	size_t indexVecTracksMin = 0;

	double clE = caloClusters->at(clusterIndex).energyDep();
	double trkE = trjExtrapols->at(trkIndex).momentum().mag();
	double tmpDeltaE = fabs(trkE - clE);

	for(size_t i=0; i<vecTracks.size(); ++i){
	  trkIndex = vecTracks[i];
	  clusterIndex = matchingDataVec[trkIndex].indexCluster(chiQ);
	  clE = caloClusters->at(clusterIndex).energyDep();
	  trkE = trjExtrapols->at(trkIndex).momentum().mag();
	  if(fabs(trkE - clE) < tmpDeltaE){
	    tmpDeltaE = fabs(trkE - clE);
	    indexVecTracksMin = i;
	  }
	}
	trkIndex = vecTracks[indexVecTracksMin];
	clusterIndex = matchingDataVec[trkIndex].indexCluster(chiQ);
	trackClusterLink->addSingle( art::Ptr<TrkToCaloExtrapol>(trjExtrapols, trkIndex),
				     art::Ptr<CaloCluster>(caloClusters, clusterIndex ) );


      }
      matchingDataVec.removeCluster(clusterIndex);
      if(_diagLevel > 2){
	cout<<"matchingDataVec.removeCluster(clusterIndex) done..."<<
	  ", clusterIndex = "<<clusterIndex<<endl;
      }
      matchingDataVec.removeTrkToCaloExtrapol(trkIndex);
      if(_diagLevel > 2){
	cout<<"matchingDataVec.removeTrkToCaloExtrapol(trkIndex) done..."<<
	  ", trkIndex = "<<trkIndex<<endl;
      }
    }

    evt.put(std::move(trackClusterLink));
    if(evt.id().event()%1000 == 0){
      cout << "Event "<<evt.id().event()<<" CaloMatching done..."<<endl;
    }
  }

}

using mu2e::CaloMatching;
DEFINE_ART_MODULE(CaloMatching);

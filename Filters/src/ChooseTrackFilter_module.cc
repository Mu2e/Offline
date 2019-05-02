//
// A Filter module that reads back the hits, track, and calorimeter clusters.  It reads in multiple collections under different track hypotheses
//
// Original author Robert Bernstein
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/unknownPDGIdName.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"

#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/HelixTraj.hh"

#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

#include "Mu2eUtilities/inc/CaloHitMCNavigator.hh"

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

#include "TDirectory.h"
#include "TNtuple.h"
#include "TTree.h"

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <vector>


#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"

using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::HepVector;
using CLHEP::keV;



namespace mu2e {


  class ChooseTrackFilter : public art::EDFilter {
     
  public:

    typedef art::Ptr<StepPointMC> StepPtr;
    typedef std::vector<StepPtr>  StepPtrs;
    typedef std::map<int,StepPtrs > HitMap;



    explicit ChooseTrackFilter(fhicl::ParameterSet const& pset);
    virtual ~ChooseTrackFilter() { }

 
    // This is called for each event.
    virtual bool filter(art::Event& e);

       



  private:
       


    int _diagLevel;
    int _nProcess;


    std::vector<std::string> _trkPatRecModuleLabel;
    std::vector<std::string> _instanceName;
    std::vector<TrkParticle> _tpart;
    std::vector<TrkFitDirection> _fdir;

    std::vector<int> _tpartTemp;
    std::vector<int> _fdirTemp;      
       

       
    int   _evt,_run;

 
    float _pmin;
    float _pmax;       

  };


  ChooseTrackFilter::ChooseTrackFilter(fhicl::ParameterSet const& pset) :
    art::EDFilter{pset},
    _diagLevel(pset.get<int>("diagLevel",0)),
    _nProcess(0),
    _trkPatRecModuleLabel(pset.get<std::vector<std::string>>("trkPatRecModuleLabel")),
    _tpartTemp(pset.get<std::vector<int>>("fitparticle")),
    _fdirTemp(pset.get<std::vector<int>>("fitdirection")),
    _pmin(pset.get<float>("pmin")),
    _pmax(pset.get<float>("pmax"))

  {

    std::cout << "constructing ChooseTrackFilter" << std::endl;

    // 
    // check these match or module is misconfigured. 
    if ( (_fdirTemp.size() != _tpartTemp.size()) || (_trkPatRecModuleLabel.size() != _tpartTemp.size())  || (_trkPatRecModuleLabel.size() != _tpartTemp.size())){
      throw cet::exception("RECO") <<"ChooseTrackFilter ctor: sizes of module label, fitparticle, and fitdirection differ" << std::endl;
    }

    for (uint i = 0; i < _fdirTemp.size(); ++i){

      //
      //kluge until this gets cleaned up...not by me...sort of a giant kluge since what I want to do is deconstruct the module label 



      if (_fdirTemp.at(i) == TrkFitDirection::downstream){
	_fdir.emplace_back(TrkFitDirection(TrkFitDirection::downstream));}
      else if (_fdirTemp.at(i) == TrkFitDirection::upstream){
	_fdir.emplace_back(TrkFitDirection(TrkFitDirection::upstream));}
    


      if (_tpartTemp.at(i) == TrkParticle::e_minus){
      	_tpart.emplace_back(TrkParticle::e_minus);}
      else if (_tpartTemp.at(i) == TrkParticle::pi_minus){
      	_tpart.emplace_back(TrkParticle(TrkParticle::pi_minus));}
      else if (_tpartTemp.at(i) == TrkParticle::e_plus){
      	_tpart.emplace_back(TrkParticle(TrkParticle::e_plus));}
      else if (_tpartTemp.at(i) == TrkParticle::pi_plus){
      	_tpart.emplace_back(TrkParticle(TrkParticle::pi_plus));}
      else throw cet::exception("RECO") << "ChooseTrackFilterParticleType misconfigured" << std::endl;
    }
 
 


    for (uint i = 0; i < _fdir.size(); ++i){
      _instanceName.emplace_back(_fdir.at(i).name() + _tpart.at(i).name());
      if (_diagLevel > 0){
	std::cout << "ChooseTrackFilter ctor: element number " << i << " instance name " << _instanceName.at(i) << std::endl;
      }
    }
  }

  bool ChooseTrackFilter::filter(art::Event& event) {

    ++_nProcess;
    if (_nProcess%10==0 && _diagLevel > 0) std::cout <<"Processing event from ChooseTrackFilter =  " <<_nProcess <<std::endl;
      
    bool cutC = false;


    // run through KalReps and see if the one we want is here
 
    for (uint ithrep = 0; ithrep < _trkPatRecModuleLabel.size();++ithrep){

      art::Handle<KalRepPtrCollection> kalRepHandle;

      if (_diagLevel > 1){
	std::cout << "ChooseTrackFilter filter : ithrep, trkPatRecModuleLabel = " << ithrep <<  "  xx" << _trkPatRecModuleLabel.at(ithrep) 
		  << "xx " << std::endl;
      }
      event.getByLabel(_trkPatRecModuleLabel.at(ithrep),_instanceName.at(ithrep),kalRepHandle);

      if (_diagLevel > 1 && !kalRepHandle.isValid() ){
	std::cout << "ChooseTrackFilter: invalid handle, loser! xx" <<  _trkPatRecModuleLabel.at(ithrep)<< "xx xx" << _instanceName.at(ithrep)
		  << "xx " << std::endl;
      }

      if (kalRepHandle.isValid()){
	if (_diagLevel > 1){
	  std::cout <<"ChooseTrackFilter: valid handle " << _trkPatRecModuleLabel.at(ithrep) << " " <<_instanceName.at(ithrep) << std::endl;
	} 
	//
	//bingo, found it.  Do stuff.
	const KalRepPtrCollection& kreps = *kalRepHandle; // the * dereferences the handle and gives ptrs; def'n of * for handles.

	for (auto const& krepPtr:kreps){ //looping over the collection and getting a ptr to each collection

	  KalRep const& krep = *krepPtr;  //derefencing the ptr

	  // Arc length from center of tracker to the most upstream point on the fitted track.
	  double s0 = krep.startValidRange();

	  // Momentum and position at s0.
	  CLHEP::Hep3Vector p0     = krep.momentum(s0);
	  HepPoint          pos0   = krep.position(s0);

	  // Some other properties at s0.
	  double loclen(0.);
	  HepVector momvec(3);
	  momvec = p0.unit();

	  const TrkSimpTraj* ltraj = krep.localTrajectory(s0,loclen);
	  BbrVectorErr momCov = krep.momentumErr(s0);
	  double fitMomErr    = sqrt(momCov.covMatrix().similarity(momvec));
	  double tanDip       = ltraj->parameters()->parameter()[HelixTraj::tanDipIndex];
	  double omega        = ltraj->parameters()->parameter()[HelixTraj::omegaIndex];
	  double d0           = ltraj->parameters()->parameter()[HelixTraj::d0Index];
	  double fitCon       = krep.chisqConsistency().significanceLevel();
	  double trkt0Err     = krep.t0().t0Err();

	  // Does this fit pass cut set C? 
	  cutC  = 
	    ( krep.fitStatus().success() >0                ) &&
	    ( fitCon             > 2.e-3                   ) &&
	    ( krep.nActive()    >= 20                      ) &&
	    ( fitMomErr          < 0.25                    ) &&
	    ( tanDip >= 0.577 && tanDip <= 1.0             ) &&
	    ( krep.t0().t0() > 700 && krep.t0().t0() < 1695.   ) &&
	    ( d0 > -80. && d0 < 105.                       ) &&
	    ( d0 + 2./omega > 450. && d0 + 2./omega < 680. ) &&
	    ( trkt0Err < 0.9                               ) &&
	    ( p0.mag() > _pmin && p0.mag() < _pmax           ) ;

	  if (_diagLevel > 1) {std:: cout << "ChooseTrackFilter t0 = " << krep.t0().t0() << " and cut set C = " << cutC << std::endl;}
	}
      }
    }
    return (cutC);
  }



}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ChooseTrackFilter);



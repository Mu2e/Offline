 //
// An EDAnalyzer module that reads back the hits created by the calorimeter and produces an ntuple
//
// $Id: CaloCheck_module.cc,v 1.4 2014/08/01 20:57:44 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:44 $
//
// Original author Bertrand Echenard
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

#include "art/Framework/Core/EDAnalyzer.h"
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
#include "TH2F.h"
#include "TH1F.h"

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <vector>




using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::keV;



namespace mu2e {


  class CaloCheck : public art::EDAnalyzer {
     
     public:

       typedef art::Ptr<StepPointMC> StepPtr;
       typedef std::vector<StepPtr>  StepPtrs;
       typedef std::map<int,StepPtrs > HitMap;



       explicit CaloCheck(fhicl::ParameterSet const& pset);
       virtual ~CaloCheck() { }

       virtual void beginJob();
       virtual void endJob();

       // This is called for each event.
       virtual void analyze(const art::Event& e);

       



     private:
       
       typedef std::vector< art::Handle<StepPointMCCollection> > HandleVector;
       typedef art::Ptr< CaloCrystalHit> CaloCrystalHitPtr;
       typedef art::Ptr<SimParticle> SimParticlePtr;


       int _diagLevel;
       int _nProcess;

       std::string _g4ModuleLabel;
       std::string _generatorModuleLabel;
       art::InputTag _simParticleTag;


       std::string _caloReadoutModuleLabel;
       std::string _caloCrystalModuleLabel;
       std::string _caloHitMCCrystalPtrLabel;
       std::string _caloClusterModuleLabel;
       std::string _caloClusterAlgorithm;
       std::string _caloClusterSeeding;
       const std::string _producerName;
       std::string _virtualDetectorLabel;
       std::string _stepPointMCLabel;
       std::string _trkPatRecModuleLabel;
       std::string _instanceName;
       TrkParticle _tpart;
       TrkFitDirection _fdir;

      
       TH2F *_crIn,*_crOut;
       TH1F *_dtime;




  };


  CaloCheck::CaloCheck(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _nProcess(0),
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _generatorModuleLabel(pset.get<string>("generatorModuleLabel")),
    _simParticleTag(pset.get<string>("simParticleTag")),
    _caloReadoutModuleLabel(pset.get<string>("caloReadoutModuleLabel")),
    _caloCrystalModuleLabel(pset.get<string>("caloCrystalModuleLabel")),
    _caloHitMCCrystalPtrLabel(pset.get<string>("calorimeterHitMCCrystalPtr")),
    _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel")),
    _virtualDetectorLabel(pset.get<string>("virtualDetectorName")),
    _stepPointMCLabel(pset.get<string>("stepPointMCLabel")),
    _trkPatRecModuleLabel(pset.get<string>("trkPatRecModuleLabel")),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _crIn(0),_crOut(0),_dtime(0)
  {
    _instanceName = _fdir.name() + _tpart.name();
  }

  void CaloCheck::beginJob(){

       art::ServiceHandle<art::TFileService> tfs;

       _crIn  = tfs->make<TH2F>("crIn","dist diff vs time diff for crystal in",100,0.,200.,100,0,100);
       _crOut = tfs->make<TH2F>("crOut","dist diff vs time diff for crystal out",100,0.,200.,100,0,100);
       _dtime = tfs->make<TH1F>("dtime","delta time in crystal",100,0.,100);

  }

	

  void CaloCheck::endJob(){
  }




  void CaloCheck::analyze(const art::Event& event) {

      ++_nProcess;
      if (_nProcess%10==0 && _diagLevel > 0) std::cout<<"Processing event from CaloCheck =  "<<_nProcess << " with instance name " << _instanceName <<std::endl;
      
 
      //Get handle to the calorimeter
      art::ServiceHandle<GeometryService> geom;
      if( ! geom->hasElement<Calorimeter>() ) return;
      Calorimeter const & cal = *(GeomHandle<Calorimeter>());
  
      //Get calo crystal hits (average from readouts)
      art::Handle<CaloCrystalHitCollection> caloCrystalHitsHandle;
      event.getByLabel(_caloCrystalModuleLabel, caloCrystalHitsHandle);
      CaloCrystalHitCollection const& caloCrystalHits(*caloCrystalHitsHandle);

      //Get calo cluster
      art::Handle<CaloClusterCollection> caloClustersHandle;
      event.getByLabel(_caloClusterModuleLabel, caloClustersHandle);
      CaloClusterCollection const& caloClusters(*caloClustersHandle);




      for (unsigned int ic=0; ic<caloCrystalHits.size();++ic) 
      {	   
	  CaloCrystalHit const& hit    = caloCrystalHits.at(ic);
	  double time0                 = hit.time();
	  CLHEP::Hep3Vector crystalPos = cal.crystal(hit.id()).position();
	  if (hit.time() < 500) continue;
	  if (hit.energyDep() < 0.9) continue;


	  int nClu(0),nCr(0);
	  std::vector<double> dmin,dtime;

	  for (auto const& clusterIt : caloClusters) 
	  {	      	      
	       for (auto const& crystalIncluster : clusterIt.caloCrystalHitsPtrVector()) 
	       {
		    if (cal.crystal(hit.id()).diskId() !=
                    cal.crystal(crystalIncluster->id()).diskId()) break;

		    if (&(*crystalIncluster) == &hit) ++nCr;
                    if (nCr > 1) std::cout<<"Warning, crystal associated to more than one cluster "<<hit.id()<<" for cluster "<<nClu<<std::endl;

		    CLHEP::Hep3Vector crystalPos2 = cal.crystal(crystalIncluster->id()).position();
		    double deltaDist = (crystalPos-crystalPos2).mag();
		    double deltaTime = abs(crystalIncluster->time()-time0);

		    dmin.push_back(deltaDist);
		    dtime.push_back(deltaTime);
               }
	       ++nClu;	         
	   }

	   if (nCr!=0) continue;
	   for (unsigned int i=0;i<dmin.size();++i) _crOut->Fill(dmin[i],dtime[i]);
       }
       
       
       for (auto const& clusterIt : caloClusters)
       {
           double dtime(0);
	   for (int i=0;i<clusterIt.size()-1;++i)
	   {
	     for (int j=i+1;j<clusterIt.size();++j)
	     {
	      double deltaTime = abs(clusterIt.caloCrystalHitsPtrVector().at(i)->time()- clusterIt.caloCrystalHitsPtrVector().at(j)->time());
	      if (deltaTime > dtime) dtime = deltaTime;
	     }	   
	   }
	   _dtime->Fill(dtime);
       }
       
       


       
   }


}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::CaloCheck);



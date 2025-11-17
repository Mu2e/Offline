 //
// An EDAnalyzer module that reads back the hits created by the calorimeter and produces an ntuple
//
//
// Original author Bertrand Echenard
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/unknownPDGIdName.hh"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"

#include "Offline/RecoDataProducts/inc/KalRepCollection.hh"
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"
#include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"


#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/TrkCaloIntersect.hh"
#include "Offline/RecoDataProducts/inc/TrkCaloMatch.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
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




using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::keV;



namespace mu2e {


  class CaloTrackMatchExample : public art::EDAnalyzer {

     public:

        explicit CaloTrackMatchExample(fhicl::ParameterSet const& pset) :
          art::EDAnalyzer(pset),
          _diagLevel(pset.get<int>("diagLevel",0)),
          _nProcess(0),
          _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel")),
          _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel")),
          _trkCaloMatchModuleLabel(pset.get<std::string>("trkCaloMatchModuleLabel")),
          _trkFitterModuleLabel(pset.get<std::string>("fitterModuleLabel")),
          _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel")),
          _virtualDetectorLabel(pset.get<std::string>("virtualDetectorName")),
          _maxChi2Match(pset.get<double>("maxChi2Match")),
          _Ntup(0)

        {
        }


        virtual ~CaloTrackMatchExample() { }

        virtual void beginJob();
        virtual void endJob();

        // This is called for each event.
        virtual void analyze(const art::Event& e);





     private:

        int findBestCluster(TrkCaloMatchCollection const& trkCaloMatches, int trkId, double maxChi2);
        int findBestTrack(  TrkCaloMatchCollection const& trkCaloMatches, int cluId, double maxChi2);

        int _diagLevel;
        int _nProcess;

        std::string      _caloCrystalModuleLabel;
        std::string      _caloClusterModuleLabel;
        std::string      _trkCaloMatchModuleLabel;
        std::string      _trkFitterModuleLabel;
        std::string      _trkfitInstanceName;
        std::string      _g4ModuleLabel;
        std::string      _virtualDetectorLabel;
        double           _maxChi2Match;





        TTree* _Ntup;

        int   _evt,_run;

        int   _nHits,_cryId[16384],_crySectionId[16384],_crySimIdx[16384],_crySimLen[16384];
        float _cryEtot,_cryTime[16384],_cryEdep[16384],_cryDose[16384],_cryPosX[16384],_cryPosY[16384],_cryPosZ[16384],_cryLeak[16384];

        int   _nCluster,_cluNcrys[16384];
        float _cluEnergy[16384],_cluTime[16384],_cluCogX[16384],_cluCogY[16384],_cluCogZ[16384];
        std::vector<std::vector<int> > _cluList;

        int   _nTrk,_trkCluIdx[8192];

        int   _nVd,_vdId[1024],_vdPdgId[1024];
        float _vdTime[1024],_vdPosX[1024],_vdPosY[1024],_vdPosZ[1024],_vdMom[1024],_vdMomX[1024],_vdMomY[1024],_vdMomZ[1024];


  };



  void CaloTrackMatchExample::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;

    _Ntup  = tfs->make<TTree>("CaloTrackMatch", "CaloTrackMatch");



    _Ntup->Branch("evt",          &_evt ,        "evt/I");
    _Ntup->Branch("run",          &_run ,        "run/I");

    _Ntup->Branch("nCry",         &_nHits ,       "nCry/I");
    _Ntup->Branch("cryId",        &_cryId ,       "cryId[nCry]/I");
    _Ntup->Branch("crySectionId", &_crySectionId, "crySectionId[nCry]/I");
    _Ntup->Branch("cryPosX",      &_cryPosX ,     "cryPosX[nCry]/F");
    _Ntup->Branch("cryPosY",      &_cryPosY ,     "cryPosY[nCry]/F");
    _Ntup->Branch("cryPosZ",      &_cryPosZ ,     "cryPosZ[nCry]/F");
    _Ntup->Branch("cryEdep",      &_cryEdep ,     "cryEdep[nCry]/F");
    _Ntup->Branch("cryTime",      &_cryTime ,     "cryTime[nCry]/F");

    _Ntup->Branch("nCluster",     &_nCluster ,    "nCluster/I");
    _Ntup->Branch("cluEnergy",    &_cluEnergy ,   "cluEnergy[nCluster]/F");
    _Ntup->Branch("cluTime",      &_cluTime ,     "cluTime[nCluster]/F");
    _Ntup->Branch("cluCogX",      &_cluCogX ,     "cluCogX[nCluster]/F");
    _Ntup->Branch("cluCogY",      &_cluCogY ,     "cluCogY[nCluster]/F");
    _Ntup->Branch("cluCogZ",      &_cluCogZ ,     "cluCogZ[nCluster]/F");
    _Ntup->Branch("cluNcrys",     &_cluNcrys ,    "cluNcrys[nCluster]/I");
    _Ntup->Branch("cluList",      &_cluList);

    _Ntup->Branch("nTrk",         &_nTrk ,        "nTrk/I");
    _Ntup->Branch("trkCluIdx",    &_trkCluIdx,    "trkCluIdx[nTrk]/I");

    _Ntup->Branch("nVd",       &_nVd ,      "nVd/I");
    _Ntup->Branch("vdId",      &_vdId ,     "vdId[nVd]/I");
    _Ntup->Branch("vdPdgId",   &_vdPdgId ,  "vdPdgId[nVd]/I");
    _Ntup->Branch("vdMom",     &_vdMom ,    "vdMom[nVd]/F");
    _Ntup->Branch("vdMomX",    &_vdMomX ,   "vdMomX[nVd]/F");
    _Ntup->Branch("vdMomY",    &_vdMomY ,   "vdMomY[nVd]/F");
    _Ntup->Branch("vdMomZ",    &_vdMomZ ,   "vdMomZ[nVd]/F");
    _Ntup->Branch("vdPosX",    &_vdPosX ,   "vdPosX[nVd]/F");
    _Ntup->Branch("vdPosY",    &_vdPosY ,   "vdPosY[nVd]/F");
    _Ntup->Branch("vdPosZ",    &_vdPosZ ,   "vdPosZ[nVd]/F");
    _Ntup->Branch("vdTime",    &_vdTime ,   "vdTime[nVd]/F");

  }



  void CaloTrackMatchExample::endJob(){
  }




  void CaloTrackMatchExample::analyze(const art::Event& event) {

      ++_nProcess;
      if (_diagLevel && _nProcess%10==0) std::cout<<"Processing event from CaloTrackMatchExample =  "<<_nProcess <<std::endl;


      // Get handle to the calorimeter
      art::ServiceHandle<GeometryService> geom;
      if( ! geom->hasElement<Calorimeter>() ) return;
      Calorimeter const & cal = *(GeomHandle<Calorimeter>());


      // Get crystal hits
      art::Handle<CaloHitCollection> CaloHitsHandle;
      event.getByLabel(_caloCrystalModuleLabel, CaloHitsHandle);
      CaloHitCollection const& CaloHits(*CaloHitsHandle);

      // Get tracks
      art::Handle<KalRepPtrCollection> trksHandle;
      event.getByLabel(_trkFitterModuleLabel,trksHandle);
      KalRepPtrCollection const& trksPtrColl(*trksHandle);

      // Get clusters
      art::Handle<CaloClusterCollection> caloClustersHandle;
      event.getByLabel(_caloClusterModuleLabel, caloClustersHandle);
      CaloClusterCollection const& caloClusters(*caloClustersHandle);

      // Get track calorimeter matches
      art::Handle<TrkCaloMatchCollection>  trkCaloMatchHandle;
      event.getByLabel(_trkCaloMatchModuleLabel, trkCaloMatchHandle);
      TrkCaloMatchCollection const& trkCaloMatches(*trkCaloMatchHandle);

      //Get virtual detector hits
      art::Handle<StepPointMCCollection> vdhits;
      event.getByLabel(_g4ModuleLabel,_virtualDetectorLabel,vdhits);








       //--------------------------  Do generated particles --------------------------------
       _evt = event.id().event();
       _run = event.run();


       //--------------------------  Do calorimeter hits --------------------------------

       _nHits = 0;
       for (unsigned int ic=0; ic<CaloHits.size();++ic)
       {
           CaloHit const& hit    = CaloHits.at(ic);
           CLHEP::Hep3Vector crystalPos = cal.crystal(hit.crystalID()).position();

           _cryTime[_nHits]      = hit.time();
           _cryEdep[_nHits]      = hit.energyDep();
           _cryPosX[_nHits]      = crystalPos.x();
           _cryPosY[_nHits]      = crystalPos.y();
           _cryPosZ[_nHits]      = crystalPos.z();
           _cryId[_nHits]        = hit.crystalID();
           _crySectionId[_nHits] = cal.crystal(hit.crystalID()).diskID();
           ++_nHits;
       }

       //--------------------------  Dump cluster info --------------------------------
       _nCluster = 0;
       _cluList.clear();
       for (CaloClusterCollection::const_iterator clusterIt = caloClusters.begin(); clusterIt != caloClusters.end(); ++clusterIt)
       {
          std::vector<int> _list;
          for (int i=0;i<clusterIt->size();++i)
             _list.push_back( int(clusterIt->caloHitsPtrVector().at(i).get()- &CaloHits.at(0)) );

          _cluEnergy[_nCluster] = clusterIt->energyDep();
          _cluTime[_nCluster]   = clusterIt->time();
          _cluNcrys[_nCluster]  = clusterIt->size();
          _cluCogX[_nCluster]   = clusterIt->cog3Vector().x();
          _cluCogY[_nCluster]   = clusterIt->cog3Vector().y();
          _cluCogZ[_nCluster]   = clusterIt->cog3Vector().z();
          _cluList.push_back(_list);

          ++_nCluster;
       }


       //--------------------------  Do tracks  --------------------------------
       // we just want the cluster that matches the track.
       // For information about the track at the entrance of the calorimeter, look at the
       // Analysis/src/ReadTrackCaloMatching_module.cc mdoule

       _nTrk = 0;
       for (unsigned int i=0;i< trksPtrColl.size(); ++i)
       {
         //KalRepPtr const& trkPtr = trksPtrColl.at(i);
         _trkCluIdx[_nTrk] = findBestCluster(trkCaloMatches,i,_maxChi2Match);
         ++_nTrk;

        }



        //--------------------------  Dump virtual detector info --------------------------------
        _nVd=0;
        for (auto iter=vdhits->begin(), ie=vdhits->end(); iter!=ie; ++iter)
        {
             const StepPointMC& hit = *iter;

             if (hit.volumeId()<73 || hit.volumeId() > 80) continue;
             if (_nVd > 999) std::cout<<"Problem, nVd = "<<_nVd <<std::endl;

             _vdId[_nVd]     = hit.volumeId();
             _vdPdgId[_nVd]  = hit.simParticle()->pdgId();
             _vdTime[_nVd]   = hit.time();
             _vdPosX[_nVd]   = hit.position().x();
             _vdPosY[_nVd]   = hit.position().y();
             _vdPosZ[_nVd]   = hit.position().z();
             _vdMom[_nVd]    = hit.momentum().mag();
             _vdMomX[_nVd]   = hit.momentum().x();
             _vdMomY[_nVd]   = hit.momentum().y();
             _vdMomZ[_nVd]   = hit.momentum().z();
             ++_nVd;
        }



          _Ntup->Fill();





  }

   int CaloTrackMatchExample::findBestCluster(TrkCaloMatchCollection const& trkCaloMatches, int trkId, double maxChi2)
   {
      double chi2Best(maxChi2), cluBest(-1);
      for (auto const& trkCaloMatch: trkCaloMatches)
      {
         if (trkCaloMatch.trkId()==trkId && trkCaloMatch.chi2() < chi2Best)
         {
           chi2Best= trkCaloMatch.chi2();
           cluBest = trkCaloMatch.cluId();
         }
      }
      return cluBest;
   }

   int CaloTrackMatchExample::findBestTrack(TrkCaloMatchCollection const& trkCaloMatches, int cluId, double maxChi2)
   {
      double chi2Best(maxChi2), trkBest(-1);
      for (auto const& trkCaloMatch: trkCaloMatches)
      {
         if (trkCaloMatch.cluId()==cluId && trkCaloMatch.chi2() < chi2Best)
         {
           chi2Best= trkCaloMatch.chi2();
           trkBest = trkCaloMatch.trkId();
         }
      }
      return trkBest;
   }


}

DEFINE_ART_MODULE(mu2e::CaloTrackMatchExample)



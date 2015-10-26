//
//
//
// $Id: ReadTrackCaloMatching_module.cc,v 1.18 2014/08/20 14:23:09 murat Exp $
// $Author: murat $
// $Date: 2014/08/20 14:23:09 $
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

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// Mu2e includes.
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "BTrk/BaBar/Constants.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "KalmanTests/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/TrkCaloIntersectCollection.hh"
#include "RecoDataProducts/inc/TrkCaloMatchCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"

#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/TrkRep.hh"
#include "BTrk/TrkBase/HelixTraj.hh"

//CLHEP includes
#include "CLHEP/Vector/TwoVector.h"
#include "BTrk/BbrGeom/HepPoint.h"


// Other includes.
#include "cetlib/exception.h"


//root includes
#include "TFile.h"
#include "TNtuple.h"


#include <cmath>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>



namespace mu2e {


  class ReadTrackCaloMatchingBis : public art::EDAnalyzer {


      public:

          explicit ReadTrackCaloMatchingBis(fhicl::ParameterSet const& pset):
            art::EDAnalyzer(pset),
            _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel")),
            _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel")),
            _trkCaloMatchModuleLabel(pset.get<std::string>("trkCaloMatchModuleLabel")),
            _trkIntersectModuleLabel(pset.get<std::string>("trkIntersectModuleLabel")),
            _trkFitterModuleLabel(pset.get<std::string>("fitterModuleLabel")),
            _tpart((TrkParticle::type)(pset.get<int>("fitparticle"))),
            _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection"))),
            _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel")),
            _virtualDetectorLabel(pset.get<std::string>("virtualDetectorName")),
            _Ntup(0)
          {
             _trkfitInstanceName = _fdir.name() + _tpart.name();
          }

          virtual ~ReadTrackCaloMatchingBis() {}
          void beginJob();
          void endJob() {}

          void analyze(art::Event const& e );



      private:

          int findBestCluster(TrkCaloMatchCollection const& trkCaloMatches, int trkId, double maxChi2);
          int findBestTrack(  TrkCaloMatchCollection const& trkCaloMatches, int cluId, double maxChi2);


          std::string      _caloCrystalModuleLabel;
          std::string      _caloClusterModuleLabel;
          std::string      _trkCaloMatchModuleLabel;
          std::string      _trkIntersectModuleLabel;
          std::string      _trkFitterModuleLabel;
          std::string      _trkfitInstanceName;
          TrkParticle      _tpart;
          TrkFitDirection  _fdir;
          std::string      _g4ModuleLabel;
          std::string      _virtualDetectorLabel;



          TTree* _Ntup;

          int    _evt;

          int    _nMatch, _mTrkId[1024],_mCluId[1024];
          float  _mChi2[1024],_mChi2Pos[1024],_mChi2Time[1024];

          int    _nTrk,_trknHit[1024],_trkStat[1024],_trkCluIdx[1024];
          float  _trkx[1024],_trky[1024],_trkz[1024],_trkFFx[1024],_trkFFy[1024],_trkFFz[1024],_trke[1024],_trkt[1024];
          float  _trkpx[1024],_trkpy[1024],_trkpz[1024],_trkprob[1024];
          float  _trkd0[1024],_trkz0[1024],_trkphi0[1024],_trkomega[1024],_trkcdip[1024],_trkdlen[1024];

          int   _nCluster,_cluNcrys[1024];
          float _cluEnergy[1024],_cluTime[1024],_cluCogX[1024],_cluCogY[1024],_cluCogZ[1024];
          float _cluE1[1024],_cluE9[1024],_cluE25[1024],_clu2Mom[1024],_cluAngle[1024];
          std::vector<std::vector<int> > _cluList;

          int   _nHits,_cryId[1024],_crySecId[1024];
          float _cryTime[1024],_cryEdep[1024],_cryPosX[1024],_cryPosY[1024],_cryPosZ[1024];

          int   _nVd,_vdId[1024],_vdPdgId[1024];
          float _vdTime[1024],_vdPosX[1024],_vdPosY[1024],_vdPosZ[1024],_vdMom[1024],_vdMomX[1024],_vdMomY[1024],_vdMomZ[1024];

  };






  void ReadTrackCaloMatchingBis::beginJob( ) {

       art::ServiceHandle<art::TFileService> tfs;
       _Ntup = tfs->make<TTree>("trkClu", "track-cluster match info");

       _Ntup->Branch("evt",      &_evt ,      "evt/I");

       _Ntup->Branch("nMatch",    &_nMatch ,   "nMatch/I");
       _Ntup->Branch("mTrkId",    &_mTrkId,    "mTrkId[nMatch]/I");
       _Ntup->Branch("mCluId",    &_mCluId,    "mCluId[nMatch]/I");
       _Ntup->Branch("mChi2",     &_mChi2,     "mChi2[nMatch]/F");
       _Ntup->Branch("mChi2Pos",  &_mChi2Pos,  "mChi2Pos[nMatch]/F");
       _Ntup->Branch("mChi2Time", &_mChi2Time, "mChi2Time[nMatch]/F");


       _Ntup->Branch("nTrk",      &_nTrk ,     "nTrk/I");
       _Ntup->Branch("trkX",      &_trkx,      "trkX[nTrk]/F");
       _Ntup->Branch("trkY",      &_trky,      "trkY[nTrk]/F");
       _Ntup->Branch("trkZ",      &_trkz,      "trkZ[nTrk]/F");
       _Ntup->Branch("trkFFX",    &_trkFFx,    "trkFFX[nTrk]/F");
       _Ntup->Branch("trkFFY",    &_trkFFy,    "trkFFY[nTrk]/F");
       _Ntup->Branch("trkFFZ",    &_trkFFz,    "trkFFZ[nTrk]/F");
       _Ntup->Branch("trkt",      &_trkt,      "trkt[nTrk]/F");
       _Ntup->Branch("trke",      &_trke,      "trke[nTrk]/F");
       _Ntup->Branch("trkpX",     &_trkpx,     "trkpX[nTrk]/F");
       _Ntup->Branch("trkpY",     &_trkpy,     "trkpY[nTrk]/F");
       _Ntup->Branch("trkpZ",     &_trkpz,     "trkpZ[nTrk]/F");
       _Ntup->Branch("trknHit",   &_trknHit,   "trknHit[nTrk]/I");
       _Ntup->Branch("trkStat",   &_trkStat,   "trkStat[nTrk]/I");
       _Ntup->Branch("trkprob",   &_trkprob,   "trkprob[nTrk]/F");
       _Ntup->Branch("trkd0",     &_trkd0,     "trkd0[nTrk]/F");
       _Ntup->Branch("trkz0",     &_trkz0,     "trkz0[nTrk]/F");
       _Ntup->Branch("trkphi0",   &_trkphi0,   "trkphi0[nTrk]/F");
       _Ntup->Branch("trkomega",  &_trkomega,  "trkomega[nTrk]/F");
       _Ntup->Branch("trkcdip",   &_trkcdip,   "trkcdip[nTrk]/F");
       _Ntup->Branch("trkdlen",   &_trkdlen,   "trkdlen[nTrk]/F");
       _Ntup->Branch("trkCluIdx", &_trkCluIdx, "trkCluIdx[nTrk]/I");

       _Ntup->Branch("nCluster",  &_nCluster ,  "nCluster/I");
       _Ntup->Branch("cluEnergy", &_cluEnergy , "cluEnergy[nCluster]/F");
       _Ntup->Branch("cluTime",   &_cluTime ,   "cluTime[nCluster]/F");
       _Ntup->Branch("cluCogX",   &_cluCogX ,   "cluCogX[nCluster]/F");
       _Ntup->Branch("cluCogY",   &_cluCogY ,   "cluCogY[nCluster]/F");
       _Ntup->Branch("cluCogZ",   &_cluCogZ ,   "cluCogZ[nCluster]/F");
       _Ntup->Branch("cluE1",     &_cluE1 ,     "cluE1[nCluster]/F");
       _Ntup->Branch("cluE9",     &_cluE1 ,     "cluE9[nCluster]/F");
       _Ntup->Branch("cluE25",    &_cluE25 ,    "cluE25[nCluster]/F");
       _Ntup->Branch("clu2Mom",   &_clu2Mom ,   "clu2Mom[nCluster]/F");
       _Ntup->Branch("cluAngle",  &_cluAngle ,  "cluAngle[nCluster]/F");
       _Ntup->Branch("cluList",   &_cluList);

       _Ntup->Branch("nCry",      &_nHits ,    "nCry/I");
       _Ntup->Branch("cryId",     &_cryId ,    "cryId[nCry]/I");
       _Ntup->Branch("crySecId",  &_crySecId,  "crySecId[nCry]/I");
       _Ntup->Branch("cryPosX",   &_cryPosX ,  "cryPosX[nCry]/F");
       _Ntup->Branch("cryPosY",   &_cryPosY ,  "cryPosY[nCry]/F");
       _Ntup->Branch("cryPosZ",   &_cryPosZ ,  "cryPosZ[nCry]/F");
       _Ntup->Branch("cryEdep",   &_cryEdep ,  "cryEdep[nCry]/F");
       _Ntup->Branch("cryTime",   &_cryTime ,  "cryTime[nCry]/F");

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



  void ReadTrackCaloMatchingBis::analyze(art::Event const& event)
  {

        Calorimeter const & cal = *(GeomHandle<Calorimeter>());

        art::Handle<KalRepPtrCollection> trksHandle;
        event.getByLabel(_trkFitterModuleLabel,_trkfitInstanceName,trksHandle);

        art::Handle<CaloClusterCollection> caloClustersHandle;
        event.getByLabel(_caloClusterModuleLabel, caloClustersHandle);
        CaloClusterCollection const& caloClusters(*caloClustersHandle);

        art::Handle<TrkCaloMatchCollection>  trkCaloMatchHandle;
        event.getByLabel(_trkCaloMatchModuleLabel, trkCaloMatchHandle);
        TrkCaloMatchCollection const& trkCaloMatches(*trkCaloMatchHandle);

        art::Handle<TrkCaloIntersectCollection>  trjIntersectHandle;
        event.getByLabel(_trkIntersectModuleLabel, trjIntersectHandle);
        TrkCaloIntersectCollection const& trkIntersect(*trjIntersectHandle);

        art::Handle<CaloCrystalHitCollection> caloCrystalHitsHandle;
        event.getByLabel(_caloCrystalModuleLabel, caloCrystalHitsHandle);
        CaloCrystalHitCollection const& caloCrystalHits(*caloCrystalHitsHandle);

        //Get virtual detector hits
        art::Handle<StepPointMCCollection> vdhits;
        event.getByLabel(_g4ModuleLabel,_virtualDetectorLabel,vdhits);



        _evt = event.id().event();


        //--------------------------  Dump matching info --------------------------------
        _nMatch=0;
        for (auto const& trkCaloMatch: trkCaloMatches)
        {
              _mTrkId[_nMatch]     = trkCaloMatch.trkId();
              _mCluId[_nMatch]     = trkCaloMatch.cluId();
              _mChi2[_nMatch]      = trkCaloMatch.chi2();
              _mChi2Pos[_nMatch]   = trkCaloMatch.chi2Pos();
              _mChi2Time[_nMatch]  = trkCaloMatch.chi2Time();

              ++_nMatch;
        }


        //--------------------------  Dump extrapolator info --------------------------------
        _nTrk=0;
        for (auto const& intersect: trkIntersect)
        {
              KalRep const& trk = *(intersect.trk());

              double pathLength             = intersect.pathLengthEntrance();
              CLHEP::Hep3Vector trkMomentum = trk.momentum(pathLength);
              double trkTime                = trk.arrivalTime(pathLength);
              double tprob                  = trk.chisqConsistency().significanceLevel();
              int    tNhits                 = trk.nActive();
              int    tStatus                = trk.fitStatus().success();


              HepPoint point = trk.position(pathLength);
              CLHEP::Hep3Vector posTrkInTracker(point.x(),point.y(),point.z());
              CLHEP::Hep3Vector posTrkInSectionFF = cal.toSectionFrameFF(intersect.sectionId(),cal.fromTrackerFrame(posTrkInTracker));

              HelixTraj trkHel(trk.helix(pathLength).params(),trk.helix(pathLength).covariance());

              _trkFFx[_nTrk]    = posTrkInSectionFF.x();
              _trkFFy[_nTrk]    = posTrkInSectionFF.y();
              _trkFFz[_nTrk]    = posTrkInSectionFF.z();
              _trkx[_nTrk]      = posTrkInTracker.x();
              _trky[_nTrk]      = posTrkInTracker.y();
              _trkz[_nTrk]      = posTrkInTracker.z();
              _trkt[_nTrk]      = trkTime;
              _trke[_nTrk]      = trkMomentum.mag();
              _trkpx[_nTrk]     = trkMomentum.x();
              _trkpy[_nTrk]     = trkMomentum.y();
              _trkpz[_nTrk]     = trkMomentum.z();
              _trknHit[_nTrk]   = tNhits;
              _trkprob[_nTrk]   = tprob;
              _trkStat[_nTrk]   = tStatus;
              _trkd0[_nTrk]     = trkHel.d0();
              _trkz0[_nTrk]     = trkHel.z0();
              _trkphi0[_nTrk]   = trkHel.phi0();
              _trkomega[_nTrk]  = trkHel.omega();
              _trkcdip[_nTrk]   = trkHel.cosDip();
              _trkdlen[_nTrk]   = intersect.pathLengthExit()-intersect.pathLengthEntrance();
              _trkCluIdx[_nTrk] = findBestCluster(trkCaloMatches,intersect.trkId(),50);

              ++_nTrk;
           }


           //--------------------------  Dump cluster info --------------------------------
           _nCluster = 0;
           _cluList.clear();
           for (CaloClusterCollection::const_iterator clusterIt = caloClusters.begin(); clusterIt != caloClusters.end(); ++clusterIt)
           {
              std::vector<int> _list;
              for (int i=0;i<clusterIt->size();++i)
                 _list.push_back( int(clusterIt->caloCrystalHitsPtrVector().at(i).get()- &caloCrystalHits.at(0)) );

              _cluEnergy[_nCluster] = clusterIt->energyDep();
              _cluTime[_nCluster]   = clusterIt->time();
              _cluNcrys[_nCluster]  = clusterIt->size();
              _cluCogX[_nCluster]   = clusterIt->cog3Vector().x();
              _cluCogY[_nCluster]   = clusterIt->cog3Vector().y();
              _cluCogZ[_nCluster]   = clusterIt->cog3Vector().z();
              _cluE1[_nCluster]     = clusterIt->e1();
              _cluE9[_nCluster]     = clusterIt->e9();
              _cluE25[_nCluster]    = clusterIt->e25();
              _clu2Mom[_nCluster]   = clusterIt->secondMoment();
              _cluAngle[_nCluster]  = clusterIt->angle();
              _cluList.push_back(_list);

              ++_nCluster;
           }


           //--------------------------  Dump crystal info --------------------------------
           _nHits=0;
           for (unsigned int ic=0; ic<caloCrystalHits.size();++ic)
           {
               CaloCrystalHit const& hit    = caloCrystalHits.at(ic);
               CLHEP::Hep3Vector crystalPos = cal.toSectionFrameFF(cal.crystal(hit.id()).sectionId(), cal.crystalOrigin(hit.id()));

               _cryTime[_nHits]      = hit.time();
               _cryEdep[_nHits]      = hit.energyDep();
               _cryPosX[_nHits]      = crystalPos.x();
               _cryPosY[_nHits]      = crystalPos.y();
               _cryPosZ[_nHits]      = crystalPos.z();
               _cryId[_nHits]        = hit.id();
               _crySecId[_nHits]     = cal.crystal(hit.id()).sectionId();

               ++_nHits;
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




   int ReadTrackCaloMatchingBis::findBestCluster(TrkCaloMatchCollection const& trkCaloMatches, int trkId, double maxChi2)
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

   int ReadTrackCaloMatchingBis::findBestTrack(TrkCaloMatchCollection const& trkCaloMatches, int cluId, double maxChi2)
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






using mu2e::ReadTrackCaloMatchingBis;
DEFINE_ART_MODULE(ReadTrackCaloMatchingBis);

//
//
//
//
// Original author G. Pezzullo
//

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
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
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/TrkCaloIntersectCollection.hh"
#include "RecoDataProducts/inc/TrkCaloMatchCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "CaloMC/inc/ClusterContentMC.hh"
#include "MCDataProducts/inc/CaloClusterMCTruthAssn.hh"

#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/TrkRep.hh"
#include "BTrk/TrkBase/HelixTraj.hh"

//CLHEP includes
#include "CLHEP/Vector/TwoVector.h"
#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/ProbTools/ChisqConsistency.hh"

// Other includes.
#include "cetlib_except/exception.h"


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


  class ReadTrackCaloMatchingMVA : public art::EDAnalyzer {


      public:

          explicit ReadTrackCaloMatchingMVA(fhicl::ParameterSet const& pset):
            art::EDAnalyzer(pset),
            _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel")),
	    _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel")),
            _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel")),
            _caloClusterTruthModuleLabel(pset.get<std::string>("caloClusterTruthModuleLabel")),
            _trkCaloMatchModuleLabel(pset.get<std::string>("trkCaloMatchModuleLabel")),
            _trkIntersectModuleLabel(pset.get<std::string>("trkIntersectModuleLabel")),
            _trkFitterModuleLabel(pset.get<std::string>("fitterModuleLabel")),
            _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel")),
            _virtualDetectorLabel(pset.get<std::string>("virtualDetectorName")),
            _Ntup(0)
          {
          }

          virtual ~ReadTrackCaloMatchingMVA() {}
          void beginJob();
          void endJob() {}

          void analyze(art::Event const& e );



      private:

          int findBestCluster(TrkCaloMatchCollection const& trkCaloMatches, int trkId, double maxChi2);
          int findBestTrack(  TrkCaloMatchCollection const& trkCaloMatches, int cluId, double maxChi2);


          std::string      _caloCrystalModuleLabel;
	  std::string      _caloReadoutModuleLabel;
          std::string      _caloClusterModuleLabel;
          std::string      _caloClusterTruthModuleLabel;
          std::string      _trkCaloMatchModuleLabel;
          std::string      _trkIntersectModuleLabel;
          std::string      _trkFitterModuleLabel;
          std::string      _trkfitInstanceName;
          std::string      _g4ModuleLabel;
          std::string      _virtualDetectorLabel;



          TTree* _Ntup;

          int    _evt;

          int    _nMatch, _mTrkId[4096],_mCluId[4096];
          float  _mChi2[4096],_mChi2Pos[4096],_mChi2Time[4096];

          int    _nTrk,_trknHit[4096],_trkStat[4096],_trkCluIdx[4096];
          float  _trkx[4096],_trky[4096],_trkz[4096],_trkFFx[4096],_trkFFy[4096],_trkFFz[4096],_trke[4096],_trkt[4096];
          float  _trkpx[4096],_trkpy[4096],_trkpz[4096],_trkprob[4096],_tkrCposX[4096],_tkrCposY[4096];
          float  _trkd0[4096],_trkz0[4096],_trkphi0[4096],_trkomega[4096],_trkcdip[4096],_trkdlen[4096];

          int   _nCluster,_nCluSim,_cluNcrys[4096];
          float _cluEnergy[4096],_cluTime[4096],_cluCogX[4096],_cluCogY[4096],_cluCogZ[4096];
          float _cluE1[4096],_cluE9[4096],_cluE25[4096],_clu2Mom[4096],_cluAngle[4096];
          int   _cluConv[16384],_cluSimIdx[16384],_cluSimLen[16384];
          std::vector<std::vector<int> > _cluList;
      
          int   _clusimId[16384],_clusimPdgId[16384],_clusimGenIdx[16384],_clusimCrCode[16384];
          float _clusimMom[16384],_clusimPosX[16384],_clusimPosY[16384],_clusimPosZ[16384],_clusimTime[16384],_clusimEdep[16384];

          int   _nHits,_cryId[4096],_crySecId[4096];
          float _cryTime[4096],_cryEdep[4096],_cryEdepErr[4096],_cryPosX[4096],_cryPosY[4096],_cryPosZ[4096];

          int   _nVd,_vdId[4096],_vdPdgId[4096],_vdenIdx[4096];
          float _vdTime[4096],_vdPosX[4096],_vdPosY[4096],_vdPosZ[4096],_vdMom[4096],_vdMomX[4096],_vdMomY[4096],_vdMomZ[4096];

  };






  void ReadTrackCaloMatchingMVA::beginJob( ) {

       art::ServiceHandle<art::TFileService> tfs;
       _Ntup = tfs->make<TTree>("trkClu", "track-cluster match info");

       _Ntup->Branch("evt",          &_evt ,      "evt/I");

       _Ntup->Branch("nMatch",       &_nMatch ,   "nMatch/I");
       _Ntup->Branch("mTrkId",       &_mTrkId,    "mTrkId[nMatch]/I");
       _Ntup->Branch("mCluId",       &_mCluId,    "mCluId[nMatch]/I");
       _Ntup->Branch("mChi2",        &_mChi2,     "mChi2[nMatch]/F");
       _Ntup->Branch("mChi2Pos",     &_mChi2Pos,  "mChi2Pos[nMatch]/F");
       _Ntup->Branch("mChi2Time",    &_mChi2Time, "mChi2Time[nMatch]/F");

       _Ntup->Branch("nTrk",         &_nTrk ,     "nTrk/I");
       _Ntup->Branch("trkX",         &_trkx,      "trkX[nTrk]/F");
       _Ntup->Branch("trkY",         &_trky,      "trkY[nTrk]/F");
       _Ntup->Branch("trkZ",         &_trkz,      "trkZ[nTrk]/F");
       _Ntup->Branch("trkFFX",       &_trkFFx,    "trkFFX[nTrk]/F");
       _Ntup->Branch("trkFFY",       &_trkFFy,    "trkFFY[nTrk]/F");
       _Ntup->Branch("trkFFZ",       &_trkFFz,    "trkFFZ[nTrk]/F");
       _Ntup->Branch("tkrCposX",     &_tkrCposX,  "tkrCposX[nTrk]/F");
       _Ntup->Branch("tkrCposY",     &_tkrCposY,  "tkrCposY[nTrk]/F");
       _Ntup->Branch("trkt",         &_trkt,      "trkt[nTrk]/F");
       _Ntup->Branch("trke",         &_trke,      "trke[nTrk]/F");
       _Ntup->Branch("trkpX",        &_trkpx,     "trkpX[nTrk]/F");
       _Ntup->Branch("trkpY",        &_trkpy,     "trkpY[nTrk]/F");
       _Ntup->Branch("trkpZ",        &_trkpz,     "trkpZ[nTrk]/F");
       _Ntup->Branch("trknHit",      &_trknHit,   "trknHit[nTrk]/I");
       _Ntup->Branch("trkStat",      &_trkStat,   "trkStat[nTrk]/I");
       _Ntup->Branch("trkprob",      &_trkprob,   "trkprob[nTrk]/F");
       _Ntup->Branch("trkd0",        &_trkd0,     "trkd0[nTrk]/F");
       _Ntup->Branch("trkz0",        &_trkz0,     "trkz0[nTrk]/F");
       _Ntup->Branch("trkphi0",      &_trkphi0,   "trkphi0[nTrk]/F");
       _Ntup->Branch("trkomega",     &_trkomega,  "trkomega[nTrk]/F");
       _Ntup->Branch("trkcdip",      &_trkcdip,   "trkcdip[nTrk]/F");
       _Ntup->Branch("trkdlen",      &_trkdlen,   "trkdlen[nTrk]/F");
       _Ntup->Branch("trkCluIdx",    &_trkCluIdx, "trkCluIdx[nTrk]/I");

       _Ntup->Branch("nCluster",     &_nCluster ,  "nCluster/I");
       _Ntup->Branch("cluEnergy",    &_cluEnergy , "cluEnergy[nCluster]/F");
       _Ntup->Branch("cluTime",      &_cluTime ,   "cluTime[nCluster]/F");
       _Ntup->Branch("cluCogX",      &_cluCogX ,   "cluCogX[nCluster]/F");
       _Ntup->Branch("cluCogY",      &_cluCogY ,   "cluCogY[nCluster]/F");
       _Ntup->Branch("cluCogZ",      &_cluCogZ ,   "cluCogZ[nCluster]/F");
       _Ntup->Branch("cluE1",        &_cluE1 ,     "cluE1[nCluster]/F");
       _Ntup->Branch("cluE9",        &_cluE1 ,     "cluE9[nCluster]/F");
       _Ntup->Branch("cluE25",       &_cluE25 ,    "cluE25[nCluster]/F");
       _Ntup->Branch("clu2Mom",      &_clu2Mom ,   "clu2Mom[nCluster]/F");
       _Ntup->Branch("cluAngle",     &_cluAngle ,  "cluAngle[nCluster]/F");
       _Ntup->Branch("cluConv",      &_cluConv ,   "cluConv[nCluster]/I");
       _Ntup->Branch("cluSimIdx",    &_cluSimIdx , "cluSimIdx[nCluster]/I");
       _Ntup->Branch("cluSimLen",    &_cluSimLen , "cluSimLen[nCluster]/I");
       _Ntup->Branch("cluList",      &_cluList);

       _Ntup->Branch("nCluSim",      &_nCluSim ,     "nCluSim/I");
       _Ntup->Branch("clusimId",     &_clusimId ,    "clusimId[nCluSim]/I");
       _Ntup->Branch("clusimPdgId",  &_clusimPdgId , "clusimPdgId[nCluSim]/I");
       _Ntup->Branch("clusimGenIdx", &_clusimGenIdx ,"clusimGenIdx[nCluSim]/I");
       _Ntup->Branch("clusimCrCode", &_clusimCrCode ,"clusimCrCode[nCluSim]/I");
       _Ntup->Branch("clusimMom",    &_clusimMom ,   "clusimMom[nCluSim]/F");
       _Ntup->Branch("clusimPosX",   &_clusimPosX ,  "clusimPosX[nCluSim]/F");
       _Ntup->Branch("clusimPosY",   &_clusimPosY ,  "clusimPosY[nCluSim]/F");
       _Ntup->Branch("clusimPosZ",   &_clusimPosZ ,  "clusimPosZ[nCluSim]/F");
       _Ntup->Branch("clusimTime",   &_clusimTime ,  "clusimTime[nCluSim]/F");
       _Ntup->Branch("clusimEdep",   &_clusimEdep ,  "clusimEdep[nCluSim]/F");

       _Ntup->Branch("nCry",         &_nHits ,    "nCry/I");
       _Ntup->Branch("cryId",        &_cryId ,    "cryId[nCry]/I");
       _Ntup->Branch("crySecId",     &_crySecId,  "crySecId[nCry]/I");
       _Ntup->Branch("cryPosX",      &_cryPosX ,  "cryPosX[nCry]/F");
       _Ntup->Branch("cryPosY",      &_cryPosY ,  "cryPosY[nCry]/F");
       _Ntup->Branch("cryPosZ",      &_cryPosZ ,  "cryPosZ[nCry]/F");
       _Ntup->Branch("cryEdep",      &_cryEdep ,  "cryEdep[nCry]/F");
       _Ntup->Branch("cryEdepErr",   &_cryEdepErr,"cryEdepErr[nCry]/F");
       _Ntup->Branch("cryTime",      &_cryTime ,  "cryTime[nCry]/F");

       _Ntup->Branch("nVd",          &_nVd ,      "nVd/I");
       _Ntup->Branch("vdId",         &_vdId ,     "vdId[nVd]/I");
       _Ntup->Branch("vdPdgId",      &_vdPdgId ,  "vdPdgId[nVd]/I");
       _Ntup->Branch("vdMom",        &_vdMom ,    "vdMom[nVd]/F");
       _Ntup->Branch("vdMomX",       &_vdMomX ,   "vdMomX[nVd]/F");
       _Ntup->Branch("vdMomY",       &_vdMomY ,   "vdMomY[nVd]/F");
       _Ntup->Branch("vdMomZ",       &_vdMomZ ,   "vdMomZ[nVd]/F");
       _Ntup->Branch("vdPosX",       &_vdPosX ,   "vdPosX[nVd]/F");
       _Ntup->Branch("vdPosY",       &_vdPosY ,   "vdPosY[nVd]/F");
       _Ntup->Branch("vdPosZ",       &_vdPosZ ,   "vdPosZ[nVd]/F");
       _Ntup->Branch("vdTime",       &_vdTime ,   "vdTime[nVd]/F");
       _Ntup->Branch("vdGenIdx",     &_vdenIdx , "vdGenIdx[nVd]/I");
  

  }



  void ReadTrackCaloMatchingMVA::analyze(art::Event const& event)
  {

        Calorimeter const & cal = *(GeomHandle<Calorimeter>());

        art::Handle<KalRepPtrCollection> trksHandle;
        event.getByLabel(_trkFitterModuleLabel,trksHandle);

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

	//Calorimeter crystal truth assignment
	art::Handle<CaloClusterMCTruthAssns> caloClusterTruthHandle;
	event.getByLabel(_caloClusterTruthModuleLabel, caloClusterTruthHandle);
	const CaloClusterMCTruthAssns& caloClusterTruth(*caloClusterTruthHandle);

        //Get virtual detector hits
        art::Handle<StepPointMCCollection> vdhits;
        event.getByLabel(_g4ModuleLabel,_virtualDetectorLabel,vdhits);

 
	std::map<art::Ptr<SimParticle>, const StepPointMC*> vdMap;
	if (vdhits.isValid())
	{
           for (auto iter=vdhits->begin(), ie=vdhits->end(); iter!=ie; ++iter)
           {
              const StepPointMC& hit = *iter;
              if (hit.volumeId()<VirtualDetectorId::EMC_Disk_0_SurfIn || hit.volumeId()>VirtualDetectorId::EMC_Disk_1_EdgeOut) continue;
	      vdMap[hit.simParticle()] = &hit;
	   }
	}


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
              CLHEP::Hep3Vector posInMu2e = cal.geomUtil().trackerToMu2e(posTrkInTracker);
              CLHEP::Hep3Vector posTrkInSectionFF = cal.geomUtil().mu2eToDiskFF(intersect.diskId(),posInMu2e);
                  
              CLHEP::Hep3Vector posTrkInCrystal = cal.geomUtil().mu2eToCrystal(cal.nearestIdxFromPosition(posInMu2e),posInMu2e);  
              
              
              
 
              HelixTraj trkHel(trk.helix(pathLength).params(),trk.helix(pathLength).covariance());

              _trkFFx[_nTrk]    = posTrkInSectionFF.x();
              _trkFFy[_nTrk]    = posTrkInSectionFF.y();
              _trkFFz[_nTrk]    = posTrkInSectionFF.z();
              _tkrCposX[_nTrk]  = posTrkInCrystal.x();
              _tkrCposY[_nTrk]  = posTrkInCrystal.y();
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
           _nCluster = _nCluSim = 0;
           _cluList.clear();
	   for (CaloClusterCollection::const_iterator clusterIt = caloClusters.begin(); clusterIt != caloClusters.end(); ++clusterIt)
	   {

               ClusterContentMC contentMC(cal, caloClusterTruth, *clusterIt);

               std::vector<int> _list;
               for (int i=0;i<clusterIt->size();++i)
               {
        	   int idx = int(clusterIt->caloCrystalHitsPtrVector().at(i).get()- &caloCrystalHits.at(0));
        	   _list.push_back(idx);
               }

               _cluEnergy[_nCluster] = clusterIt->energyDep();
               _cluTime[_nCluster]   = clusterIt->time();
               _cluNcrys[_nCluster]  = clusterIt->size();
               _cluCogX[_nCluster]   = clusterIt->cog3Vector().x(); //in disk FF frame
               _cluCogY[_nCluster]   = clusterIt->cog3Vector().y();
               _cluCogZ[_nCluster]   = clusterIt->cog3Vector().z();
               _cluConv[_nCluster]   = (contentMC.hasConversion() ? 1 : 0);
               _cluList.push_back(_list);

               _cluSimIdx[_nCluster] = _nCluSim;
               _cluSimLen[_nCluster] = contentMC.simContentMap().size();

               for (const auto& contentMap : contentMC.simContentMap() )
	       {	       
		   art::Ptr<SimParticle> sim = contentMap.first;
		   CaloContentSim       data = contentMap.second;

		   art::Ptr<SimParticle> smother(sim);
        	   while (smother->hasParent()) smother = smother->parent();
        	   int genIdx=-1;
        	   if (smother->genParticle()) genIdx = smother->genParticle()->generatorId().id();

		   //double simMom(-1);
		   CLHEP::Hep3Vector simPos(0,0,0);
		   auto vdMapEntry = vdMap.find(sim);
		   if (vdMapEntry != vdMap.end())
		   {
	              //simMom = vdMapEntry->second->momentum().mag();
		      CLHEP::Hep3Vector simPos = cal.geomUtil().mu2eToDiskFF(clusterIt->diskId(), vdMapEntry->second->position());		  
		   } 

        	   _clusimId[_nCluSim]     = sim->id().asInt();
        	   _clusimPdgId[_nCluSim]  = sim->pdgId();
        	   _clusimGenIdx[_nCluSim] = genIdx;
        	   _clusimCrCode[_nCluSim] = sim->creationCode();
        	   _clusimTime[_nCluSim]   = data.time();
        	   _clusimEdep[_nCluSim]   = data.edep();
        	   _clusimMom[_nCluSim]    = data.mom();//simMom;
        	   _clusimPosX[_nCluSim]   = simPos.x(); // in disk FF frame
        	   _clusimPosY[_nCluSim]   = simPos.y();
        	   _clusimPosZ[_nCluSim]   = simPos.z();  

        	   ++_nCluSim;
        	}

               ++_nCluster;
	   }


           //--------------------------  Dump crystal info --------------------------------
           _nHits=0;
           for (unsigned int ic=0; ic<caloCrystalHits.size();++ic)
           {
               CaloCrystalHit const& hit    = caloCrystalHits.at(ic);
               CLHEP::Hep3Vector crystalPos = cal.geomUtil().mu2eToDiskFF(cal.crystal(hit.id()).diskId(), cal.crystal(hit.id()).position());

               _cryTime[_nHits]      = hit.time();
               _cryEdep[_nHits]      = hit.energyDep();
               _cryEdepErr[_nHits]   = hit.energyDepErr();
               _cryPosX[_nHits]      = crystalPos.x();
               _cryPosY[_nHits]      = crystalPos.y();
               _cryPosZ[_nHits]      = crystalPos.z();
               _cryId[_nHits]        = hit.id();
               _crySecId[_nHits]     = cal.crystal(hit.id()).diskId();

               ++_nHits;
           }


           //--------------------------  Dump virtual detector info --------------------------------
           _nVd=0;
           for (auto iter=vdhits->begin(), ie=vdhits->end(); iter!=ie; ++iter)
           {
                const StepPointMC& hit = *iter;

                if (hit.volumeId()<73 || hit.volumeId() > 80) continue;
                if (_nVd > 999) {std::cout<<"Problem, nVd = "<<_nVd <<std::endl; continue;}

                _vdId[_nVd]     = hit.volumeId();
                _vdPdgId[_nVd]  = hit.simParticle()->pdgId();
                _vdTime[_nVd]   = hit.time();
                _vdPosX[_nVd]   = hit.position().x()+3904;
                _vdPosY[_nVd]   = hit.position().y();
                _vdPosZ[_nVd]   = hit.position().z();
                _vdMom[_nVd]    = hit.momentum().mag();
                _vdMomX[_nVd]   = hit.momentum().x();
                _vdMomY[_nVd]   = hit.momentum().y();
                _vdMomZ[_nVd]   = hit.momentum().z();
		_vdenIdx[_nVd]  = hit.simParticle()->generatorIndex();
                ++_nVd;
           }

           _Ntup->Fill();

   }




   //----------------------------------------------------------
   int ReadTrackCaloMatchingMVA::findBestCluster(TrkCaloMatchCollection const& trkCaloMatches, int trkId, double maxChi2)
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

   //----------------------------------------------------------
   int ReadTrackCaloMatchingMVA::findBestTrack(TrkCaloMatchCollection const& trkCaloMatches, int cluId, double maxChi2)
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






using mu2e::ReadTrackCaloMatchingMVA;
DEFINE_ART_MODULE(ReadTrackCaloMatchingMVA);

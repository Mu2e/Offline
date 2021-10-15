#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eDataInterface.h"

using namespace mu2e;
namespace mu2e{

  constexpr double convertToCentimeters(double val){
    return val/10;
  }
  /*------------Function delete previous event from display:-------------*/
  template <typename T, typename U> void DataLists(T data, bool Redraw, bool accumulate, std::string title, TEveElementList **List3D, TEveElementList **List2DXY = 0, TEveElementList **List2DXZ = 0 , U projection = 0){	
    if(data == 0 && Redraw){
        if (*List3D != 0){
          (*List3D)->DestroyElements();
        }

        if (*List2DXY != 0){
          (*List2DXY)->DestroyElements();

          projection->fXYMgr->ImportElements(*List2DXY, projection->fDetXYScene); 

        }if (*List2DXZ != 0){
          (*List2DXZ)->DestroyElements();

          projection->fRZMgr->ImportElements(*List2DXZ, projection->fDetRZScene);
        }
        gEve->AddElement(*List3D);
        gEve->Redraw3D(kTRUE); 
      } 
      if(data!=0){
        if (*List3D== 0) {
          *List3D = new TEveElementList((title + "3D").c_str());
          if(!accumulate){(*List3D)->DestroyElements();} if(!accumulate){(*List3D)->DestroyElements();}     
        }
        else {
          (*List3D)->DestroyElements();  
        }
        if (*List2DXY== 0) {
          *List2DXY = new TEveElementList((title + "2D").c_str());
          (*List2DXY)->IncDenyDestroy();     
        }
        else {
          if (!accumulate){(*List2DXY)->DestroyElements();}    
        }
      }
      if (*List2DXZ== 0) {
        *List2DXZ = new TEveElementList((title + "2D").c_str());
        (*List2DXZ)->IncDenyDestroy();     
      }
      else {
        if (!accumulate){(*List2DXZ)->DestroyElements();}    
      }
  }

  /*------------Function to get max and min energy in an event from the colour scheme:-------------*/
  template <typename L> void maxminE(L data, double &max, double &min){
     
      auto order = std::minmax_element(data->begin(), data->end(),
           [] (auto const& lhs, auto const& rhs) { return lhs.energyDep() < rhs.energyDep(); });

      int min_pos = order.first - data->begin();
      int max_pos = order.second - data->begin();
      min = data->at(min_pos).energyDep();
      max = data->at(max_pos).energyDep();
    }
  /*------------Function to get max and min time for time selection:-------------*/
  template <typename L> void maxminT(L data, double &max, double &min){

      auto order = std::minmax_element(data->begin(), data->end(),
         [] (auto const& lhs, auto const& rhs) { return lhs.time() < rhs.time(); });
      
      int min_pos = order.first - data->begin();
      int max_pos = order.second - data->begin();

      min = data->at(min_pos).time();
      max = data->at(max_pos).time();
            
    }
  /*------------Function to get CRV information since the words are different:-------------*/
  template <typename L> void maxminCRV(L data, double &max, double &min){
      auto order = std::minmax_element(data->begin(), data->end(),
           [] (auto const& lhs, auto const& rhs) { return lhs.GetPulseTime() < rhs.GetPulseTime(); });
      int min_pos = order.first - data->begin();
      int max_pos = order.second - data->begin();
      min = data->at(min_pos).GetPulseTime();
      max = data->at(max_pos).GetPulseTime();
  }
   
  /*------------Function to build energy list:-------------*/
  template <typename L> std::vector<double> Energies(L data, int *energylevels[]){
    std::vector<double> energies = {0, 0};
    double Max_Energy = 0;
    double Min_Energy = 1000;
    maxminE(data, Max_Energy, Min_Energy);
    
    double interval = (Max_Energy - Min_Energy)/(9);

    for(unsigned int i=0; i<data->size();i++){
      for(unsigned int n=0; n<9;n++){
        if(((*data)[i]).energyDep() >= Min_Energy + n * interval && ((*data)[i]).energyDep() <=Min_Energy + (n+1)*interval){
        (*energylevels)[i] = n;}
      }
    }
    energies.at(0) = Min_Energy;
    energies.at(1) = Max_Energy;
    
    return energies;
  }
  
  /*------------Function to build time ordering:-------------*/
  std::vector<double> TEveMu2eDataInterface::getTimeRange(bool firstloop, const ComboHitCollection *chcol, const CrvRecoPulseCollection *crvcoincol, const CaloClusterCollection *clustercol, const CaloHitCollection *cryHitcol, bool addCRVInfo, bool addHits, bool addCalo){
	  std::vector <double> time = {-1, -1};
    double max, min;
    std::vector<double> alltime;
   
    if(addCRVInfo){
      if(crvcoincol->size() !=0){
        maxminCRV(crvcoincol, max, min);
        alltime.push_back(max);
        alltime.push_back(min);
      }
    }

    if(addHits){
      if (chcol->size() != 0){
        maxminT(chcol, max, min);
        alltime.push_back(max);
        alltime.push_back(min);
      }
    }
 
    if (addCalo){
      if(clustercol->size() !=0){
        maxminT(clustercol, max, min);
        alltime.push_back(max);
        alltime.push_back(min);
      }
    }
    
    if(alltime.size() !=0){
      auto order = std::minmax_element(alltime.begin(), alltime.end(),
         [] (auto const& lhs, auto const& rhs) { return lhs < rhs; });
      int min_pos = order.first - alltime.begin();
      int max_pos = order.second - alltime.begin();
      time.at(0) = alltime.at(min_pos);
      time.at(1) = alltime.at(max_pos);
    }
 
    if (time.at(0) == -1){time.at(0) = 0;}
    if (time.at(1) == -1){time.at(1) = 100;} 
    return time;
  }
  
 /*------------Function to add CRV information to the display:-------------*/
  void TEveMu2eDataInterface::AddCRVInfo(bool firstloop, const CrvRecoPulseCollection *crvcoincol,  double min_time, double max_time, TEveMu2e2DProjection *CRV2Dproj, bool Redraw, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2){
     
    DataLists<const CrvRecoPulseCollection*, TEveMu2e2DProjection*>(crvcoincol, Redraw, accumulate,  "CRVRecoPulse", &fCrvList3D, &fCrvList2DXY,&fCrvList2DYZ, CRV2Dproj);
    
    if(crvcoincol->size() !=0){
    
      TEveElementList *CrvList2DXY = new TEveElementList("CrvData2DXY");
      TEveElementList *CrvList2DYZ = new TEveElementList("CrvData2DYZ");
      TEveElementList *CrvList3D = new TEveElementList("CrvData3D");
      GeomHandle<CosmicRayShield> CRS;
      for(unsigned int i=0; i <crvcoincol->size(); i++)
      {
        const CrvRecoPulse &crvRecoPulse = crvcoincol->at(i);
        TEveMu2eCRVEvent *teve_crv3D = new TEveMu2eCRVEvent(crvRecoPulse);
        TEveMu2eCRVEvent *teve_crv2DXY = new TEveMu2eCRVEvent(crvRecoPulse);
	TEveMu2eCRVEvent *teve_crv2DYZ = new TEveMu2eCRVEvent(crvRecoPulse);
        const CRSScintillatorBarIndex &crvBarIndex = crvRecoPulse.GetScintillatorBarIndex();
        const CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndex);
        CLHEP::Hep3Vector crvCounterPos = crvCounter.getPosition();
        CLHEP::Hep3Vector pointInMu2e = crvCounterPos;
        hep3vectorTocm(crvCounterPos);
        string pos3D = "(" + to_string((double)crvCounterPos.x()) + ", " + to_string((double)crvCounterPos.y()) + ", " + to_string((double)crvCounterPos.z()) + ")";
        if((min_time == -1 && max_time == -1) or (crvRecoPulse.GetPulseTime() > min_time and crvRecoPulse.GetPulseTime() < max_time)){
          teve_crv3D->DrawHit3D("CRVHits3D, Position = " + pos3D + ", Pulse Time = " + to_string(crvRecoPulse.GetPulseTime()) + ", Pulse Height = "+
          to_string(crvRecoPulse.GetPulseHeight()) + "Pulse Width = " + to_string(crvRecoPulse.GetPulseTime()),  i + 1, pointInMu2e, CrvList3D);
          
	  teve_crv2DXY->DrawHit2DXY("CRVHits2D, Position = " + pos3D + ", Pulse Time = " + to_string(crvRecoPulse.GetPulseTime()) + ", Pulse Height = "+ to_string(crvRecoPulse.GetPulseHeight()) + "Pulse Width = " +
          to_string(crvRecoPulse.GetPulseTime()),  i + 1, crvCounterPos, CrvList2DXY);
          teve_crv2DYZ->DrawHit2DYZ("CRVHits2D, Position = " + pos3D + ", Pulse Time = " + to_string(crvRecoPulse.GetPulseTime()) + ", Pulse Height = "+ to_string(crvRecoPulse.GetPulseHeight()) + "Pulse Width = " +
          to_string(crvRecoPulse.GetPulseTime()),  i + 1, crvCounterPos, CrvList2DYZ);
          fCrvList3D->AddElement(CrvList3D);
          fCrvList2DXY->AddElement(CrvList2DXY);
          fCrvList2DYZ->AddElement(CrvList2DYZ);
        }
      }
          
      CRV2Dproj->fXYMgr->ImportElements(fCrvList2DXY, CRV2Dproj->fEvtXYScene);
      CRV2Dproj->fRZMgr->ImportElements(fCrvList2DYZ, CRV2Dproj->fEvtRZScene);
      
      gEve->AddElement(fCrvList3D);
      gEve->Redraw3D(kTRUE);
    }
  }
            

    
  /*------------Function to add ComboHits to Tracker in 3D and 2D displays:-------------*/
  std::vector<double> TEveMu2eDataInterface::AddComboHits(bool firstloop, const ComboHitCollection *chcol, TEveMu2e2DProjection *tracker2Dproj, bool Redraw, double min_energy, double max_energy, double min_time, double max_time, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2){
   
    std::vector<double> energies = {0,0};
    DataLists<const ComboHitCollection*, TEveMu2e2DProjection*>(chcol, Redraw, accumulate, "ComboHit", &fHitsList3D, &fHitsList2DXY, &fHitsList2DXZ, tracker2Dproj);
    /*
    TXYMgr->ImportElements(fHitsList2D, scene1); 
    TRZMgr->ImportElements(fHitsList2D, scene2); */
    if(chcol!=0){
      TEveElementList *HitList2DXY = new TEveElementList("ComboHits2DXY");
      TEveElementList *HitList2DXZ = new TEveElementList("ComboHits2DXZ");
      TEveElementList *HitList3D = new TEveElementList("ComboHits3D");

      int *energylevels = new int[chcol->size()];
      energies = Energies<const ComboHitCollection*>(chcol, &energylevels);

      for(size_t i=0; i<chcol->size();i++){
        ComboHit hit = (*chcol)[i];
        TEveMu2eHit *teve_hit2DXY = new TEveMu2eHit(hit);
	TEveMu2eHit *teve_hit2DXZ = new TEveMu2eHit(hit);
        TEveMu2eHit *teve_hit3D = new TEveMu2eHit(hit);
        
        CLHEP::Hep3Vector HitPos(hit.pos().x(), hit.pos().y(), hit.pos().z());
        GeomHandle<DetectorSystem> det;
        CLHEP::Hep3Vector pointInMu2e = det->toMu2e(HitPos);
        string energy = to_string(teve_hit3D->GetEnergy());
        string pos3D = "(" + to_string((double)pointInMu2e.x()) + ", " + to_string((double)pointInMu2e.y()) + ", " + to_string((double)pointInMu2e.z()) + ")";
        string pos2D = "(" + to_string((double)hit.pos().x()) + ", " + to_string((double)hit.pos().y()) + ", " + to_string((double)hit.pos().z()) + ")";
        if (((min_time == -1 && max_time== -1) || (hit.time() > min_time && hit.time() < max_time)) && ((hit.energyDep() >= min_energy && hit.energyDep() <= max_energy) || (min_energy == -1 && max_energy == -1))){
          teve_hit3D->DrawHit3D("ComboHits3D, Position = " + pos3D + ", Energy = " + energy + ", Time = " + to_string(hit.time()) + ", ", i + 1,  pointInMu2e, energylevels[i], HitList3D);
          teve_hit2DXY->DrawHit2DXY("ComboHits2D, Position = " + pos2D + ", Energy = " + energy + ", Time = " + to_string(hit.time()) + ", ", i + 1, HitPos,energylevels[i], HitList2DXY);
          teve_hit2DXZ->DrawHit2DXZ("ComboHits2D, Position = " + pos2D + ", Energy = " + energy + ", Time = " + to_string(hit.time()) + ", ", i + 1, HitPos,energylevels[i], HitList2DXZ);

          fHitsList2DXY->AddElement(HitList2DXY);
          fHitsList2DXZ->AddElement(HitList2DXZ); 
          fHitsList3D->AddElement(HitList3D); 
        }
      }
      tracker2Dproj->fXYMgr->ImportElements(fHitsList2DXY, tracker2Dproj->fEvtXYScene); 
      tracker2Dproj->fRZMgr->ImportElements(fHitsList2DXZ, tracker2Dproj->fEvtRZScene);
      /*TXYMgr->ImportElements(fHitsList2D, scene1);
      TRZMgr->ImportElements(fHitsList2D, scene2);*/
      gEve->AddElement(HitList3D);
      gEve->Redraw3D(kTRUE); 
    }
    return energies;
  }

  /*------------Function to add Calo Clusters to 3D and 2D display:-------------*/
  std::vector<double> TEveMu2eDataInterface::AddCaloClusters(bool firstloop, const CaloClusterCollection *clustercol, TEveMu2e2DProjection *calo2Dproj, bool Redraw, double min_energy, double max_energy, double min_time, double max_time, bool accumulate, TEveProjectionManager *CfXYMgr, TEveProjectionManager *CfRZMgr, TEveScene *scene1, TEveScene *scene2) {
    vector <double> energies = {0, 0};

    DataLists<const CaloClusterCollection*, TEveMu2e2DProjection*>(clustercol, Redraw, accumulate, "CaloCluster", &fClusterList3D, &fClusterList2D_disk0, &fClusterList2D_disk0,calo2Dproj);
    DataLists<const CaloClusterCollection*, TEveMu2e2DProjection*>(clustercol, Redraw, accumulate, "CaloCluster", &fClusterList3D, &fClusterList2D_disk1, &fClusterList2D_disk1,calo2Dproj);
    //CfXYMgr->ImportElements(fClusterList2D_disk0, scene1); 
    //CfRZMgr->ImportElements(fClusterList2D_disk1, scene2); 
    
    if(clustercol->size()!=0){ 
      TEveElementList *ClusterList3D = new TEveElementList("CaloClusters3D");
      TEveElementList *ClusterList2D_disk0 = new TEveElementList("CaloClusters2D_Disk0");
      TEveElementList *ClusterList2D_disk1 = new TEveElementList("CaloClusters2D_Disk1");
      int *energylevels = new int[clustercol->size()];
      energies = Energies<const CaloClusterCollection*>(clustercol, &energylevels);
      Calorimeter const &cal = *(GeomHandle<Calorimeter>());
      std::vector<CLHEP::Hep3Vector> hits;
      bool addHits = false; //To add crystal hits 
      bool addCrystals = true; // To highlight the crystals included in a give cluster
      for(unsigned int i=0; i<clustercol->size();i++){
        CaloCluster const  &cluster= (*clustercol)[i];
        TEveMu2eCluster *teve_cluster3D = new TEveMu2eCluster(cluster);
        TEveMu2eCluster *teve_cluster2D = new TEveMu2eCluster(cluster);
        
        if(addCrystals){
          for(unsigned h =0 ; h < cluster.caloHitsPtrVector().size();h++)     {
            TEvePointSet *crystals2D = new TEvePointSet();
            art::Ptr<CaloHit>  crystalhit = cluster.caloHitsPtrVector()[h];
            int cryID = crystalhit->crystalID();
            
            int diskID = cal.crystal(crystalhit->crystalID()).diskID();
            Crystal const &crystal = cal.crystal(cryID);
            double crystalXLen = pointmmTocm(crystal.size().x());
            double crystalYLen = pointmmTocm(crystal.size().y());
            double crystalZLen = pointmmTocm(crystal.size().z());

            CLHEP::Hep3Vector crystalPos   = cal.geomUtil().mu2eToDiskFF(0,crystal.position());
            Double_t origin[3];
            crystals2D ->SetMarkerStyle(9); //TODO - use name
            crystals2D ->SetMarkerSize(1);
            crystals2D ->SetMarkerColor(kRed);
            crystalPos = PointToCalo(crystalPos, diskID);
            hep3vectorTocm(crystalPos);
            origin [0] = crystalPos.x();
            origin [1] = crystalPos.y();
            origin [2] = crystalPos.z();
            TEveGeoShape *crystalShape   = new TEveGeoShape();
            crystalShape->SetFillColor(kRed);
            crystalShape->SetShape(new TGeoBBox("cryHit", (crystalXLen/2), (crystalYLen/2), convertToCentimeters((crystalZLen/2)), origin));
            crystals2D->SetNextPoint(origin [0],origin [1],origin [2]);
            if(diskID == 0){
              fClusterList2D_disk0->AddElement(crystals2D);
              crystals2D->SetPickable(kTRUE);
              std::string info = "Crystal Id: "+to_string(cryID)+" Disk Id: "+to_string(diskID)+" Energy Dep: "+to_string(crystalhit->energyDep()) + "+/-"+to_string(crystalhit->energyDepErr())+" time: "+to_string(crystalhit->time())+" +/- "+to_string(crystalhit->timeErr());
              crystals2D->SetTitle((info).c_str());
            } 
            if(diskID == 1){
              fClusterList2D_disk1->AddElement(crystals2D);
              crystals2D->SetPickable(kTRUE);
              std::string info = "Crystal Id: "+to_string(cryID)+" Disk Id: "+to_string(diskID)+" Energy Dep: "+to_string(crystalhit->energyDep()) + "+/-"+to_string(crystalhit->energyDepErr())+" time: "+to_string(crystalhit->time())+" +/- "+to_string(crystalhit->timeErr());
              crystals2D->SetTitle((info).c_str());
            }
          }
        }
        if(addHits){
          for(unsigned h =0 ; h < cluster.caloHitsPtrVector().size();h++)     {
            art::Ptr<CaloHit>  crystalhit = cluster.caloHitsPtrVector()[h];
            int diskID = cal.crystal(crystalhit->crystalID()).diskID();
            CLHEP::Hep3Vector HitPos(cal.geomUtil().mu2eToDiskFF(diskID, cal.crystal(crystalhit->crystalID()).position()));
            CLHEP::Hep3Vector hit = PointToCalo(HitPos,diskID);
            hep3vectorTocm(hit);
            hits.push_back(hit);
          }
        }
        
        CLHEP::Hep3Vector COG(cluster.cog3Vector().x(),cluster.cog3Vector().y(), cluster.cog3Vector().z());
        
        CLHEP::Hep3Vector pointInMu2e2D = PointToCalo(COG,cluster.diskID());
        CLHEP::Hep3Vector pointInMu2e3D = PointToCaloCM(COG,cluster.diskID());
       
        string pos3D = "(" + to_string(convertToCentimeters((double)pointInMu2e3D.x())) + ", " + to_string(convertToCentimeters((double)pointInMu2e3D.y())) + ", " + to_string(convertToCentimeters((double)pointInMu2e3D.z())) + ")";
        string pos2D = "(" + to_string((double)COG.x()) + ", " + to_string((double)COG.y()) + ", " + to_string((double)COG.z()) + ")";

        if (((min_time == -1 && max_time == -1) || (cluster.time() > min_time &&  cluster.time() < max_time )) && ((cluster.energyDep() >= min_energy && cluster.energyDep() <= max_energy) || (min_energy == -1 && max_energy == -1))){
          teve_cluster3D->DrawCluster("CaloCluster3D, Cluster #" + to_string(i + 1) + ", Position =" + pos3D + ", Energy = " + to_string(cluster.energyDep()) + "+/- " + to_string(cluster.energyDepErr()) + ", Time = " + to_string(cluster.time()) + " +/- " + to_string(cluster.timeErr()),pointInMu2e3D, energylevels[i], ClusterList3D, hits, addHits);
          fClusterList3D->AddElement(ClusterList3D); 
          
          if(cluster.diskID()==0){
            teve_cluster2D->DrawCluster("CaloCluster3D, Cluster #" + to_string(i + 1) + ", Position =" + pos2D + ", Energy = " + to_string(cluster.energyDep()) + "+/- " + to_string(cluster.energyDepErr()) +  ", Time = " + to_string(cluster.time()) + " +/- " + to_string(cluster.timeErr()), pointInMu2e2D,energylevels[i], ClusterList2D_disk0, hits, addHits); 
            fClusterList2D_disk0->AddElement(ClusterList2D_disk0); 
            calo2Dproj->fXYMgr->ImportElements(fClusterList2D_disk0, calo2Dproj->fDetXYScene); 
            //CfXYMgr->ImportElements(fClusterList2D_disk0, scene1); For Multiview

          }
          if(cluster.diskID()==1){
            teve_cluster2D->DrawCluster("CaloCluster3D, Cluster #" + to_string(i + 1) + ", Position =" + pos2D + ", Energy = " + to_string(cluster.energyDep()) + "+/- " + to_string(cluster.energyDepErr()) + ", Time = " + to_string(cluster.time())+ " +/- " + to_string(cluster.timeErr()), pointInMu2e2D,energylevels[i], ClusterList2D_disk1, hits, addHits); 
            fClusterList2D_disk1->AddElement(ClusterList2D_disk1); 
            calo2Dproj->fRZMgr->ImportElements(fClusterList2D_disk1, calo2Dproj->fDetRZScene); 
            //CfRZMgr->ImportElements(fClusterList2D_disk1, scene2); For MultiView
          }
           
        }
      }
      gEve->AddElement(fClusterList3D);
      gEve->Redraw3D(kTRUE);
    }
      return energies;
  }

  /*------------Function to add Crystal hits to 2D and 3D Calo displays:-------------*/
  void TEveMu2eDataInterface::AddCrystalHits(bool firstloop, const CaloHitCollection *cryHitcol, TEveMu2e2DProjection *calo2Dproj, double min_time, double max_time, bool Redraw, bool accumulate, TEveProjectionManager *CfXYMgr, TEveProjectionManager *CfRZMgr, TEveScene *scene1, TEveScene *scene2){
    vector <double> energies = {0, 0};
    Calorimeter const &cal = *(GeomHandle<Calorimeter>());

    DataLists<const CaloHitCollection*, TEveMu2e2DProjection*>(cryHitcol, Redraw, accumulate, "CrystalHit", &fHitsList3D, &fHitsList2DXY,&fHitsList2DXZ, calo2Dproj);
    /*CfXYMgr->ImportElements(fClusterList2D_disk0, scene1); 
    CfRZMgr->ImportElements(fClusterList2D_disk1, scene2); */
    if(cryHitcol!=0){

      TEveElementList *HitList = new TEveElementList("CrystalHits");
      int *energylevels = new int[cryHitcol->size()];

      energies = Energies<const CaloHitCollection*>(cryHitcol, &energylevels);

      for(unsigned int i=0; i<cryHitcol->size();i++){
        TEveMu2eHit *teve_hit = new TEveMu2eHit();
        CaloHit const  &hit = (*cryHitcol)[i];
        int diskID = cal.crystal(hit.crystalID()).diskID();
        CLHEP::Hep3Vector HitPos(cal.geomUtil().mu2eToDiskFF(diskID, cal.crystal(hit.crystalID()).position()));
        CLHEP::Hep3Vector pointInMu2e = PointToCalo(HitPos,diskID);
        if (((min_time == -1 && max_time == -1) || (hit.time() > min_time && hit.time() < max_time ))){
          teve_hit->DrawHit3D("CrystalHits",  1, pointInMu2e, energylevels[i], HitList);
          fCrystalHitList->AddElement(HitList);    
        }
      }
    /*CfXYMgr->ImportElements(fClusterList2D_disk0, scene1); 
    CfRZMgr->ImportElements(fClusterList2D_disk1, scene2); */
    gEve->AddElement(fCrystalHitList);
    gEve->Redraw3D(kTRUE);  
    }
  }

  /*------------Function to add Kalman Reco Helix to the display:-------------*/
  void TEveMu2eDataInterface::AddHelixPieceWise3D(bool firstloop, std::tuple<std::vector<std::string>, std::vector<const KalSeedCollection*>> track_tuple, TEveMu2e2DProjection *tracker2Dproj, double min_time, double max_time, bool Redraw, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2){
  
    std::vector<const KalSeedCollection*> track_list = std::get<1>(track_tuple);
    std::vector<std::string> names = std::get<0>(track_tuple);
    std::vector<int> colour;
    for(unsigned int j=0; j< track_list.size(); j++){
      const KalSeedCollection* seedcol = track_list[j];
      colour.push_back(j+3);
      DataLists<const KalSeedCollection*, TEveMu2e2DProjection*>(seedcol, Redraw, accumulate, "HelixTrack", &fTrackList3D, &fTrackList2DXY,&fTrackList2DXZ, tracker2Dproj);

      /*TXYMgr->ImportElements(fTrackList2D, scene1); 
      TRZMgr->ImportElements(fTrackList2D, scene2); */
      if(seedcol!=0){  
        for(unsigned int k = 0; k < seedcol->size(); k = k + 20){   
          KalSeed kseed = (*seedcol)[k];
          const std::vector<mu2e::KalSegment> &segments = kseed.segments();
          size_t nSegments=segments.size();
          if(nSegments==0) continue;
          const mu2e::KalSegment &segmentFirst = kseed.segments().front();
          const mu2e::KalSegment &segmentLast = kseed.segments().back();
          double fltLMin=segmentFirst.fmin();
          double fltLMax=segmentLast.fmax();
          TEveMu2eCustomHelix *line = new TEveMu2eCustomHelix();
          TEveMu2eCustomHelix *line_twoDXY = new TEveMu2eCustomHelix();
          TEveMu2eCustomHelix *line_twoDXZ = new TEveMu2eCustomHelix();
          line->fKalSeed_ = kseed;
          line->SetSeedInfo(kseed);
  
          for(size_t m=0; m<nSegments; m++){
            const mu2e::KalSegment &segment = segments.at(m);
            fltLMin=segment.fmin();
            fltLMax=segment.fmax();
            if(m>0){
              double fltLMaxPrev=segments.at(m-1).fmax();
              fltLMin=(fltLMin+fltLMaxPrev)/2.0;
            }
            if(m+1<nSegments){
              double fltLMinNext=segments.at(m+1).fmin();
              fltLMax=(fltLMax+fltLMinNext)/2.0;
            }
            for(double fltL=fltLMin; fltL<=fltLMax; fltL+=1.0){
	            GeomHandle<DetectorSystem> det;
              XYZVectorF pos;
              segment.helix().position(fltL,pos);
              CLHEP::Hep3Vector p = GenVector::Hep3Vec(pos);
              CLHEP::Hep3Vector InMu2e = det->toMu2e(p);
              line->SetPostionAndDirectionFromKalRep((InMu2e.z()));              
              line->SetNextPoint((InMu2e.x()), (InMu2e.y()), (InMu2e.z()));
              line_twoDXY->SetNextPoint(pointmmTocm(p.x()), pointmmTocm(p.y()), pointmmTocm(p.z()));
              line_twoDXZ->SetNextPoint(pointmmTocm(p.x()), pointmmTocm(p.y())+1000, pointmmTocm(p.z()));
            }
          }
          
          line_twoDXY->SetLineColor(colour[j]);
          line_twoDXY->SetLineWidth(3);
          fTrackList2DXY->AddElement(line_twoDXY);
          
          line_twoDXZ->SetLineColor(colour[j]);
          line_twoDXZ->SetLineWidth(3);
          fTrackList2DXZ->AddElement(line_twoDXZ);

          line->SetPickable(kTRUE);
          const std::string title = "Helix: " + names[j] + "PDG Code = " + to_string(line->PDGcode_) +", Momentum = " + to_string(line->Momentum_)+ ", Time = " + to_string(line->Time_);
          line->SetTitle(Form(title.c_str()));
          line->SetLineColor(colour[j]);
          line->SetLineWidth(3);
          fTrackList3D->AddElement(line);
        }
        
        /*TXYMgr->ImportElements(fTrackList2D, scene1);
        TRZMgr->ImportElements(fTrackList2D, scene2);*/
        tracker2Dproj->fXYMgr->ImportElements(fTrackList2DXY, tracker2Dproj->fEvtXYScene);
        tracker2Dproj->fRZMgr->ImportElements(fTrackList2DXZ, tracker2Dproj->fEvtRZScene);
        gEve->AddElement(fTrackList3D);
        gEve->Redraw3D(kTRUE);
      }
      
    }
  }

  /*------------Function to add No Field cosmic track fit to 2D and 3D display:-------------*/
  void TEveMu2eDataInterface::AddCosmicTrack(bool firstloop, const CosmicTrackSeedCollection *cosmiccol, TEveMu2e2DProjection *tracker2Dproj, double min_time, double max_time, bool Redraw, bool accumulate,  TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2){
  if(cosmiccol !=0){
    DataLists<const CosmicTrackSeedCollection*, TEveMu2e2DProjection*>(cosmiccol, Redraw, accumulate,"CosmicTrack", &fTrackList3D, &fTrackList2DXY,&fTrackList2DXZ, tracker2Dproj);
    /*TXYMgr->ImportElements(fTrackList2D, scene1); 
    TRZMgr->ImportElements(fTrackList2D, scene2);*/ 
    TEveMu2eStraightTrack *line = new TEveMu2eStraightTrack();
    TEveMu2eStraightTrack *line2DXY = new TEveMu2eStraightTrack();
    TEveMu2eStraightTrack *line2DXZ = new TEveMu2eStraightTrack();
    for(unsigned int ist = 0; ist < cosmiccol->size(); ++ist){

      CosmicTrackSeed sts =(*cosmiccol)[ist];
      CosmicTrack st = sts._track;

      line->SetLineColor(kGreen);

      double endY = 0;
      double startY = 0;

      endY = sts._straw_chits[sts._straw_chits.size()-1].pos().y();
      startY = sts._straw_chits[0].pos().y(); 

      Float_t ty1 = startY;
      Float_t ty2 = endY;
      Float_t tx1 = st.MinuitParams.A0  - st.MinuitParams.A1*ty1;
      Float_t tx2 = st.MinuitParams.A0  - st.MinuitParams.A1*ty2;
      Float_t tz1 = st.MinuitParams.B0  - st.MinuitParams.B1*ty1;
      Float_t tz2 = st.MinuitParams.B0  - st.MinuitParams.B1*ty2; 	
      line->AddLine((tx1-3904), (ty1), (tz1+10171), (tx2-3904), (ty2), (tz2+10171));
      const std::string title = "CosmicTrack #" + to_string(ist + 1) + ", Parameters: A0:" + to_string(st.MinuitParams.A0) + ", A1:" + to_string(st.MinuitParams.A1) + ", B0:" + to_string(st.MinuitParams.B0) + ", B1:" + to_string(st.MinuitParams.B1);
      line->SetTitle(Form(title.c_str()));

      line2DXY->AddLine(tx1, ty1, tz1, tx2, ty2, tz2);	
      line2DXY->SetPickable(kTRUE);
      line2DXY->SetLineColor(kGreen);
      line2DXY->SetLineWidth(3);
      fTrackList2DXY->AddElement(line2DXY);
      
      line2DXZ->AddLine(tx1, ty1+1000, tz1, tx2, ty2+1000, tz2);	
      line2DXZ->SetPickable(kTRUE);
      line2DXZ->SetLineColor(kGreen);
      line2DXZ->SetLineWidth(3);
      fTrackList2DXZ->AddElement(line2DXZ);

      line->SetPickable(kTRUE);
      line->SetLineColor(kGreen);
      line->SetLineWidth(3);
      fTrackList3D->AddElement(line);

    }

    tracker2Dproj->fXYMgr->ImportElements(fTrackList2DXY, tracker2Dproj->fEvtXYScene);
    tracker2Dproj->fRZMgr->ImportElements(fTrackList2DXZ, tracker2Dproj->fEvtRZScene);

    /*TXYMgr->ImportElements(fTrackList2D, scene1);
    TRZMgr->ImportElements(fTrackList2D, scene2);*/
    gEve->AddElement(fTrackList3D);
    gEve->Redraw3D(kTRUE);
    }
  }
  
  /*TODO - Note this function is untested and awaits compatibility in current Offline. We made it for future upgrades.
  void TEveMu2eDataInterface::AddTrackExitTrajectories(bool firstloop, const TrkExtTrajCollection *trkextcol, double min_time, double max_time, bool Redraw, bool accumulate) {
  DataLists<const TrkExtTrajCollection*, TEveMu2e2DProjection*>(trkextcol, Redraw, accumulate, "TrackExitTraj", &fExtTrackList3D);
    if(trkextcol!=0){
       for(unsigned int i = 0; i < trkextcol->size(); i++){
        TrkExtTraj trkext = (*trkextcol)[i];
        TEveMu2eCustomHelix *line = new TEveMu2eCustomHelix();
        line->fTrkExtTraj_ = trkext;

        line->SetMomentumExt();
        line->SetParticleExt();

        line->SetPoint(0,trkext.front().x(),trkext.front().y(),trkext.front().z());
        for(unsigned int k = 0 ; k< trkext.size(); k+=10){
          const mu2e::TrkExtTrajPoint &trkextpoint = trkext[k];
          line->SetNextPoint(trkextpoint.x(), trkextpoint.y(), trkextpoint.z());  //TODO accurate translation
        }
        //     fExtTrackList2D->AddElement(line);
        //     tracker2Dproj->fXYMgr->ImportElements(fExtTrackList2D, tracker2Dproj->fDetXYScene);
        //     tracker2Dproj->fRZMgr->ImportElements(fExtTrackList2D, tracker2Dproj->fDetRZScene);
        line->SetPickable(kTRUE);
        const std::string title = "Helix #" + to_string(i + 1) + ", Momentum = " + to_string(line->Momentum_);
        line->SetTitle(Form(title.c_str()));
        line->SetLineColor(kRed);
        line->SetLineWidth(3);
        fExtTrackList3D->AddElement(line);

      }
        gEve->AddElement(fExtTrackList3D);
        gEve->Redraw3D(kTRUE);
    }

  }*/
}

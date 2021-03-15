#include <TApplication.h>
#include <TEvePad.h>
#include <TObject.h>
#include <TSystem.h>
#include <limits>
#include <vector>
#include <tuple>
#include <algorithm>
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eDataInterface.h"
#include "TEveEventDisplay/src/dict_classes/GeomUtils.h"


using namespace mu2e;
namespace mu2e{

  template <typename T, typename U> void DataLists(T data, bool Redraw, bool accumulate, std::string title, TEveElementList **List3D, TEveElementList **List2D = 0, U projection = 0){	
    if(data == 0 && Redraw){
      if (*List3D != 0){
        (*List3D)->DestroyElements();
      }
      if (*List2D != 0){
          (*List2D)->DestroyElements();
        }
	projection->fXYMgr->ImportElements(*List2D, projection->fDetXYScene); 
        projection->fRZMgr->ImportElements(*List2D, projection->fDetRZScene);
      
      gEve->AddElement(*List3D);
      gEve->Redraw3D(kTRUE); 
    }
    if(data!=0){
      if (*List3D==0) {
        *List3D = new TEveElementList((title + "3D").c_str());
        (*List3D)->IncDenyDestroy();     
      }
      else {
        if (!accumulate){(*List3D)->DestroyElements();}   
      }
      if (*List2D== 0) {
        *List2D = new TEveElementList((title + "2D").c_str());
        (*List2D)->IncDenyDestroy();     
      }
      else {
        if (!accumulate){(*List2D)->DestroyElements();}   
      }
    }
  }

 
template <typename L> void maxminE(L data, double &max, double &min){
    auto order = std::minmax_element(data->begin(), data->end(),
         [] (auto const& lhs, auto const& rhs) { return lhs.energyDep() < rhs.energyDep(); });
    int min_pos = order.first - data->begin();
    int max_pos = order.second - data->begin();
    min = data->at(min_pos).energyDep();
    max = data->at(max_pos).energyDep();
  }
template <typename L> void maxminT(L data, double &max, double &min){
    auto order = std::minmax_element(data->begin(), data->end(),
       [] (auto const& lhs, auto const& rhs) { return lhs.time() < rhs.time(); });
    int min_pos = order.first - data->begin();
    int max_pos = order.second - data->begin();
    min = data->at(min_pos).time();
    max = data->at(max_pos).time();
  }
template <typename L> void maxminCRV(L data, double &max, double &min){
    auto order = std::minmax_element(data->begin(), data->end(),
         [] (auto const& lhs, auto const& rhs) { return lhs.GetPulseTime() < rhs.GetPulseTime(); });
    int min_pos = order.first - data->begin();
    int max_pos = order.second - data->begin();
    min = data->at(min_pos).GetPulseTime();
    max = data->at(max_pos).GetPulseTime();
}
 
  
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

  std::vector<double> TEveMu2eDataInterface::getTimeRange(bool firstloop, const ComboHitCollection *chcol, const CrvRecoPulseCollection *crvcoincol, const CaloClusterCollection *clustercol, const CaloHitCollection *cryHitcol){
	  std::vector <double> time = {-1, -1};
    double max, min;
    std::vector<double> alltime;
    if (crvcoincol != 0){
      maxminCRV(crvcoincol, max, min);
      alltime.push_back(max);
      alltime.push_back(min);
    }
    
    if (chcol != 0){
      maxminT(chcol, max, min);
      alltime.push_back(max);
      alltime.push_back(min);
    }

    if (clustercol != 0){
      maxminT(clustercol, max, min);
      alltime.push_back(max);
      alltime.push_back(min);
    }

    if (cryHitcol != 0){
      maxminT(cryHitcol, max, min);
      alltime.push_back(max);
      alltime.push_back(min);
    }
    auto order = std::minmax_element(alltime.begin(), alltime.end(),
       [] (auto const& lhs, auto const& rhs) { return lhs < rhs; });
    int min_pos = order.first - alltime.begin();
    int max_pos = order.second - alltime.begin();
    time.at(0) = alltime.at(min_pos);
    time.at(1) = alltime.at(max_pos);
    
    if (time.at(0) == -1){time.at(0) = 0;}
    if (time.at(1) == -1){time.at(1) = 100;} 
    return time;
  }

   void TEveMu2eDataInterface::AddCRVInfo(bool firstloop, const CrvRecoPulseCollection *crvcoincol,  double min_time, double max_time, bool Redraw, bool accumulate){
      DataLists<const CrvRecoPulseCollection*, TEveMu2e2DProjection*>(crvcoincol, Redraw, accumulate,  "CRVRecoPulse", &fCrvList3D, &fCrvList2D);

        if(crvcoincol!=0){
          TEveElementList *CrvList3D = new TEveElementList("CrvData3D");
          GeomHandle<CosmicRayShield> CRS;
          for(unsigned int i=0; i <crvcoincol->size(); i++)
          {
            const CrvRecoPulse &crvRecoPulse = crvcoincol->at(i);

            TEveMu2eCRVEvent *teve_crv3D = new TEveMu2eCRVEvent(crvRecoPulse);
            const CRSScintillatorBarIndex &crvBarIndex = crvRecoPulse.GetScintillatorBarIndex();
            const CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndex);
            CLHEP::Hep3Vector crvCounterPos = crvCounter.getPosition(); 
            hep3vectorTocm(crvCounterPos);
            string pos3D = "(" + to_string((double)crvCounterPos.x()) + ", " + to_string((double)crvCounterPos.y()) + ", " + to_string((double)crvCounterPos.z()) + ")";
            if( (min_time == -1 && max_time == -1) or (crvRecoPulse.GetPulseTime() > min_time and crvRecoPulse.GetPulseTime() < max_time)){
              teve_crv3D->DrawHit3D("CRVHits3D, Position = " + pos3D + ", Pulse Time = " + to_string(crvRecoPulse.GetPulseTime()) + ", Pulse Height = "+ to_string(crvRecoPulse.GetPulseHeight()) + + "Pulse Width = " + to_string(crvRecoPulse.GetPulseTime()),  i + 1, crvCounterPos, CrvList3D);
              fCrvList3D->AddElement(CrvList3D); 
            }
          } 
          gEve->AddElement(fCrvList3D);
          gEve->Redraw3D(kTRUE); 
        }
    }

   std::vector<double> TEveMu2eDataInterface::AddComboHits(bool firstloop, const ComboHitCollection *chcol, TEveMu2e2DProjection *tracker2Dproj, bool Redraw, double min_energy, double max_energy, double min_time, double max_time, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2){
    
    std::vector<double> energies = {0,0};
    DataLists<const ComboHitCollection*, TEveMu2e2DProjection*>(chcol, Redraw, accumulate, "ComboHit", &fHitsList3D, &fHitsList2D, tracker2Dproj);
	TXYMgr->ImportElements(fHitsList2D, scene1); 
        TRZMgr->ImportElements(fHitsList2D, scene2); 
    if(chcol!=0){
      TEveElementList *HitList2D = new TEveElementList("ComboHits2D");
      TEveElementList *HitList3D = new TEveElementList("ComboHits3D");

      int *energylevels = new int[chcol->size()];
      energies = Energies<const ComboHitCollection*>(chcol, &energylevels);

      for(size_t i=0; i<chcol->size();i++){
        ComboHit hit = (*chcol)[i];
        TEveMu2eHit *teve_hit2D = new TEveMu2eHit(hit);
        TEveMu2eHit *teve_hit3D = new TEveMu2eHit(hit);
        
        CLHEP::Hep3Vector HitPos(hit.pos().x(), hit.pos().y(), hit.pos().z());
        GeomHandle<DetectorSystem> det;
        CLHEP::Hep3Vector pointInMu2e = det->toMu2e(HitPos);
        string energy = to_string(teve_hit3D->GetEnergy());
        string pos3D = "(" + to_string((double)pointInMu2e.x()) + ", " + to_string((double)pointInMu2e.y()) + ", " + to_string((double)pointInMu2e.z()) + ")";
        string pos2D = "(" + to_string((double)hit.pos().x()) + ", " + to_string((double)hit.pos().y()) + ", " + to_string((double)hit.pos().z()) + ")";
        if (((min_time == -1 && max_time== -1) || (hit.time() > min_time && hit.time() < max_time)) && ((hit.energyDep() >= min_energy && hit.energyDep() <= max_energy) || (min_energy == -1 && max_energy == -1))){
          teve_hit3D->DrawHit3D("ComboHits3D, Position = " + pos3D + ", Energy = " + energy + ", Time = " + to_string(hit.time()) + ", ", i + 1,  pointInMu2e, energylevels[i], HitList3D);
          teve_hit2D->DrawHit2D("ComboHits2D, Position = " + pos2D + ", Energy = " + energy + ", Time = " + to_string(hit.time()) + ", ", i + 1, HitPos,energylevels[i], HitList2D);

          fHitsList2D->AddElement(HitList2D); 
          fHitsList3D->AddElement(HitList3D); 
    
          //TXYMgr->ImportElements(HitList2D);
	  //TRZMgr->ImportElements(HitList3D);
        }
      }
      
      //tracker2Dproj->fXYMgr->ImportElements(HitList2D, tracker2Dproj->fDetXYScene); 
      //tracker2Dproj->fRZMgr->ImportElements(HitList2D, tracker2Dproj->fDetRZScene);
      
      TXYMgr->ImportElements(fHitsList2D, scene1);
      TRZMgr->ImportElements(fHitsList2D, scene2);
      gEve->AddElement(HitList3D);
      gEve->Redraw3D(kTRUE); 
    }
    return energies;
  }


  std::vector<double> TEveMu2eDataInterface::AddCaloClusters(bool firstloop, const CaloClusterCollection *clustercol, TEveMu2e2DProjection *calo2Dproj, bool Redraw, double min_energy, double max_energy, double min_time, double max_time, bool accumulate, TEveProjectionManager *CfXYMgr, TEveProjectionManager *CfRZMgr, TEveScene *scene1, TEveScene *scene2){
    vector <double> energies = {0, 0};
    std::cout<<"Time interval"<<min_time<<" "<<max_time<<std::endl;
  DataLists<const CaloClusterCollection*, TEveMu2e2DProjection*>(clustercol, Redraw, accumulate, "CaloCluster", &fClusterList3D, &fClusterList2D_disk0, calo2Dproj);
  DataLists<const CaloClusterCollection*, TEveMu2e2DProjection*>(clustercol, Redraw, accumulate, "CaloCluster", &fClusterList3D, &fClusterList2D_disk1, calo2Dproj);
  CfXYMgr->ImportElements(fClusterList2D_disk0, scene1); 
  CfRZMgr->ImportElements(fClusterList2D_disk1, scene2); 

    if(clustercol!=0){ 
      TEveElementList *ClusterList3D = new TEveElementList("CaloClusters3D");
      TEveElementList *ClusterList2D_disk0 = new TEveElementList("CaloClusters2D_Disk0");
      TEveElementList *ClusterList2D_disk1 = new TEveElementList("CaloClusters2D_Disk1");
      int *energylevels = new int[clustercol->size()];
      energies = Energies<const CaloClusterCollection*>(clustercol, &energylevels);
      Calorimeter const &cal = *(GeomHandle<Calorimeter>());
      std::vector<CLHEP::Hep3Vector> hits;
      bool addHits = false;
      for(unsigned int i=0; i<clustercol->size();i++){
        CaloCluster const  &cluster= (*clustercol)[i];
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
        TEveMu2eCluster *teve_cluster3D = new TEveMu2eCluster(cluster);
        TEveMu2eCluster *teve_cluster2D = new TEveMu2eCluster(cluster);
       
        CLHEP::Hep3Vector COG(cluster.cog3Vector().x(),cluster.cog3Vector().y(), cluster.cog3Vector().z());
        
        CLHEP::Hep3Vector pointInMu2e = PointToCalo(COG,cluster.diskID());
       
        string pos3D = "(" + to_string((double)pointInMu2e.x()) + ", " + to_string((double)pointInMu2e.y()) + ", " + to_string((double)pointInMu2e.z()) + ")";
        string pos2D = "(" + to_string((double)COG.x()) + ", " + to_string((double)COG.y()) + ", " + to_string((double)COG.z()) + ")";

        if (((min_time == -1 && max_time == -1) || (cluster.time() > min_time &&  cluster.time() < max_time )) && ((cluster.energyDep() >= min_energy && cluster.energyDep() <= max_energy) || (min_energy == -1 && max_energy == -1))){
          teve_cluster3D->DrawCluster("CaloCluster3D, Cluster #" + to_string(i + 1) + ", Position =" + pos3D + ", Energy = " + to_string(cluster.energyDep()) + ", Time = " + to_string(cluster.time()), pointInMu2e, energylevels[i], ClusterList3D, hits, addHits);
            fClusterList3D->AddElement(ClusterList3D); 
            
            //fClusterList2D->AddElement(CrystalList2D); 
            if(cluster.diskID()==0){
              teve_cluster2D->DrawCluster("CaloCluster3D, Cluster #" + to_string(i + 1) + ", Position =" + pos2D + ", Energy = " + to_string(cluster.energyDep()) + ", Time = " + to_string(cluster.time()), pointInMu2e,energylevels[i], ClusterList2D_disk0, hits, addHits); 
              fClusterList2D_disk0->AddElement(ClusterList2D_disk0); 
              calo2Dproj->fXYMgr->ImportElements(fClusterList2D_disk0, calo2Dproj->fDetXYScene); 
              CfXYMgr->ImportElements(fClusterList2D_disk0, scene1); 
            }
            if(cluster.diskID()==1){
             teve_cluster2D->DrawCluster("CaloCluster3D, Cluster #" + to_string(i + 1) + ", Position =" + pos2D + ", Energy = " + to_string(cluster.energyDep()) + ", Time = " + to_string(cluster.time()), pointInMu2e,energylevels[i], ClusterList2D_disk1, hits, addHits); 
             fClusterList2D_disk1->AddElement(ClusterList2D_disk1); 
             calo2Dproj->fRZMgr->ImportElements(fClusterList2D_disk1, calo2Dproj->fDetRZScene); 
             
             CfRZMgr->ImportElements(fClusterList2D_disk1, scene2); 
            }
        }
      }
      
      gEve->AddElement(fClusterList3D);
      gEve->Redraw3D(kTRUE);
    }
      return energies;
  }

  
  void TEveMu2eDataInterface::AddCrystalHits(bool firstloop, const CaloHitCollection *cryHitcol, TEveMu2e2DProjection *calo2Dproj, double min_time, double max_time, bool Redraw, bool accumulate, TEveProjectionManager *CfXYMgr, TEveProjectionManager *CfRZMgr, TEveScene *scene1, TEveScene *scene2){
    vector <double> energies = {0, 0};
    Calorimeter const &cal = *(GeomHandle<Calorimeter>());
    DataLists<const CaloHitCollection*, TEveMu2e2DProjection*>(cryHitcol, Redraw, accumulate, "CrystalHit", &fHitsList3D, &fHitsList2D, calo2Dproj);
    CfXYMgr->ImportElements(fClusterList2D_disk0, scene1); 
    CfRZMgr->ImportElements(fClusterList2D_disk1, scene2); 
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
    CfXYMgr->ImportElements(fClusterList2D_disk0, scene1); 
    CfRZMgr->ImportElements(fClusterList2D_disk1, scene2); 
    gEve->AddElement(fCrystalHitList);
    gEve->Redraw3D(kTRUE);  
    }
  }

  void TEveMu2eDataInterface::AddHelixPieceWise2D(bool firstloop, std::tuple<std::vector<std::string>, std::vector<const KalSeedCollection*>> track_tuple, TEveMu2e2DProjection *tracker2Dproj, double min_time, double max_time, bool Redraw, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2){
  
    std::vector<const KalSeedCollection*> track_list = std::get<1>(track_tuple);

    std::vector<std::string> names = std::get<0>(track_tuple);
    std::vector<int> colour;
    for(unsigned int j=0; j< track_list.size(); j++){
      const KalSeedCollection* seedcol = track_list[j];
      colour.push_back(j+3);
      DataLists<const KalSeedCollection*, TEveMu2e2DProjection*>(seedcol, Redraw, accumulate, "HelixTrack", &fTrackList3D, &fTrackList2D, tracker2Dproj);
      TXYMgr->ImportElements(fTrackList2D, scene1); 
      TRZMgr->ImportElements(fTrackList2D, scene2); 
      if(seedcol!=0){  
        for(unsigned int k = 0; k < seedcol->size(); k = k + 20){
          KalSeed kseed = (*seedcol)[k];
          TEveMu2eCustomHelix *line = new TEveMu2eCustomHelix();
          line->fKalSeed = kseed;
          line->SetSeedInfo(kseed);

          unsigned int nSteps = 50;  
          double kStepSize = 61;//nSteps/(TrackerLength())+70; //+ 70;///TODO CaloLength() +
          for(unsigned int i = 0 ; i< nSteps; i++){
            double zpos = (i*kStepSize)-TrackerLength()/2;
            if(i==0) {
              line->SetPostionAndDirectionFromKalRep(zpos);
              CLHEP::Hep3Vector Pos(line->Position.x(), line->Position.y(), zpos+line->Position.z());
              double x = (Pos.x()+line->Direction.x()*line->Momentum);
              double y = (Pos.y()+line->Direction.y()*line->Momentum);
              double z = (Pos.z());
              CLHEP::Hep3Vector Point(x,y,z);
              line->SetPoint(i,pointmmTocm(Point.x()), pointmmTocm(Point.y()), pointmmTocm(Point.z()));
            } else {
              line->SetPostionAndDirectionFromKalRep(zpos);
              CLHEP::Hep3Vector Pos(line->Position.x(), line->Position.y(), zpos+line->Position.z());
              double x = (Pos.x()+line->Direction.x()*line->Momentum);
              double y = (Pos.y()+line->Direction.y()*line->Momentum);
              double z = (Pos.z());
              CLHEP::Hep3Vector Point(x,y,z);
              line->SetNextPoint(pointmmTocm(Point.x()), pointmmTocm(Point.y()), pointmmTocm(Point.z()));
              }
          }
        
          line->SetLineColor(colour[j]);
          line->SetLineWidth(3);

          const std::string title = "Helix: " + names[j] + "PDG Code = " + to_string(line->PDGcode) +", Momentum = " + to_string(line->Momentum)+ ", Time = " + to_string(line->Time);
          line->SetTitle(Form(title.c_str()));
          fTrackList2D->AddElement(line);
      }
        
      TXYMgr->ImportElements(fTrackList2D, scene1);
      TRZMgr->ImportElements(fTrackList2D, scene2);
      gEve->AddElement(fTrackList3D);
      gEve->Redraw3D(kTRUE);
      }
      
    }
    }
    
  void TEveMu2eDataInterface::AddHelixPieceWise3D(bool firstloop, std::tuple<std::vector<std::string>, std::vector<const KalSeedCollection*>> track_tuple, TEveMu2e2DProjection *tracker2Dproj, double min_time, double max_time, bool Redraw, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2){
  
    std::vector<const KalSeedCollection*> track_list = std::get<1>(track_tuple);

    std::vector<std::string> names = std::get<0>(track_tuple);
    std::vector<int> colour;
    for(unsigned int j=0; j< track_list.size(); j++){

      const KalSeedCollection* seedcol = track_list[j];
      colour.push_back(j+3);
      DataLists<const KalSeedCollection*, TEveMu2e2DProjection*>(seedcol, Redraw, accumulate, "HelixTrack", &fTrackList3D, &fTrackList2D, tracker2Dproj);
      TXYMgr->ImportElements(fTrackList2D, scene1); 
      TRZMgr->ImportElements(fTrackList2D, scene2); 
      if(seedcol!=0){  
        for(unsigned int k = 0; k < seedcol->size(); k = k + 20){
          KalSeed kseed = (*seedcol)[k];
          TEveMu2eCustomHelix *line = new TEveMu2eCustomHelix();
          TEveMu2eCustomHelix *line_twoD = new TEveMu2eCustomHelix();
          line->fKalSeed = kseed;
          line->SetSeedInfo(kseed);

          unsigned int nSteps = 500;  
          double kStepSize = 6.1;//nSteps/(TrackerLength())+70; //+ 70;///TODO CaloLength() +
          for(unsigned int i = 0 ; i< nSteps; i++){
            double zpos = (i*kStepSize)-TrackerLength()/2;
            GeomHandle<DetectorSystem> det;
            if(i==0) {
              line->SetPostionAndDirectionFromKalRep(zpos);
              CLHEP::Hep3Vector Pos(line->Position.x(), line->Position.y(), zpos+line->Position.z());
              std::cout<<"track "<<Pos.x()<<" "<<Pos.y()<<" "<<Pos.z()<<std::endl;
              double x = (Pos.x()+line->Direction.x()*line->Momentum);
              double y = (Pos.y()+line->Direction.y()*line->Momentum);
              double z = (Pos.z());
              CLHEP::Hep3Vector Point(x,y,z);
              CLHEP::Hep3Vector InMu2e = det->toMu2e(Point);
              line->SetPoint(i,pointmmTocm(InMu2e.x()), pointmmTocm(InMu2e.y()), pointmmTocm(InMu2e.z()));
              line_twoD->SetPoint(i,pointmmTocm(Point.x()), pointmmTocm(Point.y()), pointmmTocm(Point.z()));
            } else {
              line->SetPostionAndDirectionFromKalRep(zpos);
              CLHEP::Hep3Vector Pos(line->Position.x(), line->Position.y(), zpos+line->Position.z());
              std::cout<<Pos.x()<<" "<<Pos.y()<<" "<<Pos.z()<<std::endl;
              double x = (Pos.x()+line->Direction.x()*line->Momentum);
              double y = (Pos.y()+line->Direction.y()*line->Momentum);
              double z = (Pos.z());
              CLHEP::Hep3Vector Point(x,y,z);
              CLHEP::Hep3Vector InMu2e = det->toMu2e(Point);
              line->SetNextPoint(pointmmTocm(InMu2e.x()), pointmmTocm(InMu2e.y()), pointmmTocm(InMu2e.z()));
              line_twoD->SetNextPoint(pointmmTocm(Point.x()), pointmmTocm(Point.y()), pointmmTocm(Point.z()));
              }
          }
        
          line_twoD->SetLineColor(colour[j]);
          line_twoD->SetLineWidth(3);
          fTrackList2D->AddElement(line_twoD);

        line->SetPickable(kTRUE);
        const std::string title = "Helix: " + names[j] + "PDG Code = " + to_string(line->PDGcode) +", Momentum = " + to_string(line->Momentum)+ ", Time = " + to_string(line->Time);
        line->SetTitle(Form(title.c_str()));
        line->SetLineColor(colour[j]);
        line->SetLineWidth(3);
        fTrackList3D->AddElement(line);
      }
        
      TXYMgr->ImportElements(fTrackList2D, scene1);
      TRZMgr->ImportElements(fTrackList2D, scene2);
      gEve->AddElement(fTrackList3D);
      gEve->Redraw3D(kTRUE);
      }
      
    }
  }

  void TEveMu2eDataInterface::AddCosmicTrack(bool firstloop, const CosmicTrackSeedCollection *cosmiccol, TEveMu2e2DProjection *tracker2Dproj, double min_time, double max_time, bool Redraw, bool accumulate,  TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2){
     if(cosmiccol !=0){
      DataLists<const CosmicTrackSeedCollection*, TEveMu2e2DProjection*>(cosmiccol, Redraw, accumulate,"CosmicTrack", &fTrackList3D, &fTrackList2D, tracker2Dproj);
	TXYMgr->ImportElements(fTrackList2D, scene1); 
        TRZMgr->ImportElements(fTrackList2D, scene2); 
      TEveMu2eStraightTrack *line = new TEveMu2eStraightTrack();
      TEveMu2eStraightTrack *line2D = new TEveMu2eStraightTrack();
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
        line->AddLine(pointmmTocm(tx1-3904), pointmmTocm(ty1), pointmmTocm(tz1+10171), pointmmTocm(tx2-3904), pointmmTocm(ty2), pointmmTocm(tz2+10171));
        const std::string title = "CosmicTrack #" + to_string(ist + 1) + ", Parameters: A0:" + to_string(st.MinuitParams.A0) + ", A1:" + to_string(st.MinuitParams.A1) + ", B0:" + to_string(st.MinuitParams.B0) + ", B1:" + to_string(st.MinuitParams.B1);
        line->SetTitle(Form(title.c_str()));
        
          line2D->AddLine(tx1, ty1, tz1, tx2, ty2, tz2);	
          line2D->SetPickable(kTRUE);
          line2D->SetLineColor(kGreen);
          line2D->SetLineWidth(3);
          fTrackList2D->AddElement(line2D);
        
        line->SetPickable(kTRUE);

        line->SetLineColor(kGreen);
        line->SetLineWidth(3);
        fTrackList3D->AddElement(line);

        }
        
          tracker2Dproj->fXYMgr->ImportElements(fTrackList2D, tracker2Dproj->fDetXYScene);
          tracker2Dproj->fRZMgr->ImportElements(fTrackList2D, tracker2Dproj->fDetRZScene);
        
        TXYMgr->ImportElements(fTrackList2D, scene1);
      TRZMgr->ImportElements(fTrackList2D, scene2);
        gEve->AddElement(fTrackList3D);
        gEve->Redraw3D(kTRUE);
    }
  }
  
  //TODO - Note this function is untested and awaits compatibility in current Offline. We made it for future upgrades.
  void TEveMu2eDataInterface::AddTrackExitTrajectories(bool firstloop, const TrkExtTrajCollection *trkextcol, double min_time, double max_time, bool Redraw, bool accumulate) {
  DataLists<const TrkExtTrajCollection*, TEveMu2e2DProjection*>(trkextcol, Redraw, accumulate, "TrackExitTraj", &fExtTrackList3D);
    if(trkextcol!=0){
       for(unsigned int i = 0; i < trkextcol->size(); i++){
        TrkExtTraj trkext = (*trkextcol)[i];
        TEveMu2eCustomHelix *line = new TEveMu2eCustomHelix();
        line->fTrkExtTraj = trkext;

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
        const std::string title = "Helix #" + to_string(i + 1) + ", Momentum = " + to_string(line->Momentum);
        line->SetTitle(Form(title.c_str()));
        line->SetLineColor(kRed);
        line->SetLineWidth(3);
        fExtTrackList3D->AddElement(line);

      }
        gEve->AddElement(fExtTrackList3D);
        gEve->Redraw3D(kTRUE);
    }

}
  

}

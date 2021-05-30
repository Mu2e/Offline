
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eMCInterface.h"

using namespace mu2e;
namespace mu2e{

  /*------------Function to clear lists for new events:-------------*/
  template <typename T, typename U> void DataLists(T data, bool Redraw, bool accumulate, string title, TEveElementList **List3D, TEveElementList **List2D = 0, U projection = 0){	
      if(data == 0 && Redraw){
        if (*List3D != 0){
          (*List3D)->DestroyElements();
        }
        
          if (*List2D != 0){
            (*List2D)->DestroyElements();
          
          projection->fXYMgr->ImportElements(*List2D, projection->fDetXYScene); 
          projection->fRZMgr->ImportElements(*List2D, projection->fDetRZScene);
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
        if (*List2D== 0) {
          *List2D = new TEveElementList((title + "2D").c_str());
          (*List2D)->IncDenyDestroy();     
        }
        else {
          if (!accumulate){(*List2D)->DestroyElements();}    
        }
    }
  }

  /*------------Function to add straight line MC Trajectory i.e. for Comsics in No field:-------------*/
  void TEveMu2eMCInterface::AddSimpleMCTrajectory(bool firstloop, const MCTrajectoryCollection *trajcol, TEveMu2e2DProjection *tracker2Dproj, bool Redraw, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2){    
        DataLists<const MCTrajectoryCollection*, TEveMu2e2DProjection*>(trajcol, Redraw, accumulate, "MC Trajectory", &fTrackList3D, &fTrackList2D, tracker2Dproj);
    if(trajcol!=0){
      TEveElementList *HitList3D = new TEveElementList("MCtraj3D");
      TEveElementList *HitList2D = new TEveElementList("MCtraj2D");
      std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator trajectoryIter;
      for(trajectoryIter=trajcol->begin(); trajectoryIter!=trajcol->end(); trajectoryIter++)
      {
        const std::vector<MCTrajectoryPoint> &points = trajectoryIter->second.points();
        string pdgId= to_string(trajectoryIter->first->pdgId());
        CLHEP::Hep3Vector StartHitPos(points[0].x(), points[0].y(), points[0].z());
        CLHEP::Hep3Vector EndHitPos(points[points.size()-1].x(), points[points.size()-1].y(), points[points.size()-1].z());
        TEveMu2eMCTraj *teve_hit3D = new TEveMu2eMCTraj();
        string energy = to_string(points[0].kineticEnergy());
        teve_hit3D->DrawSimpleLine("MCTraj PDG " + pdgId + "Energy = " + energy  + ", ",  StartHitPos, EndHitPos, HitList3D);

        fTrackList3D->AddElement(HitList3D);

        GeomHandle<DetectorSystem> det;
        StartHitPos = det->toMu2e(StartHitPos);
        EndHitPos = det->toMu2e(EndHitPos);
        TEveMu2eMCTraj *teve_hit2D = new TEveMu2eMCTraj();
        teve_hit2D->DrawSimpleLine("MCTraj PDG " + pdgId + "Energy = " + energy + ", ", StartHitPos, EndHitPos, HitList2D);
        fTrackList2D->AddElement(HitList2D); 

      }
      
        tracker2Dproj->fXYMgr->ImportElements(fTrackList2D, tracker2Dproj->fDetXYScene); 
        tracker2Dproj->fRZMgr->ImportElements(fTrackList2D, tracker2Dproj->fDetRZScene);
      
      gEve->AddElement(fTrackList3D);
      gEve->Redraw3D(kTRUE);
      }
  
  }
  
  /*------------Function to help user select a list of PDG codes to display:-------------*/
  int TEveMu2eMCInterface::Contains(std::vector<int> v, int x)
  {
    return std::count(v.begin(), v.end(), x);
  }
  
  /*------------
  Function to display MCTracjories of any shape, these are made up of a series of TEveLines, in the same way as Reco Helices:      
  -------------*/
    /*void TEveMu2eMCInterface::AddFullMCTrajectory(bool firstloop, const MCTrajectoryCollection *trajcol, TEveMu2e2DProjection *tracker2Dproj, bool Redraw, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2, std::vector<int> particleIds){
	DataLists<const MCTrajectoryCollection*, TEveMu2e2DProjection*>(trajcol, Redraw, accumulate, "MC Trajectory", &fTrackList3D, &fTrackList2D, tracker2Dproj);
	      std::cout<<"Full 1 "<<std::endl;
	      
        //TXYMgr->ImportElements(fTrackList2D, scene1); 
        //TRZMgr->ImportElements(fTrackList2D, scene2); 
        if(trajcol!=0){
          std::cout<<"Full 2 "<<std::endl;
          TEveMu2eCustomHelix *line = new TEveMu2eCustomHelix();
          TEveMu2eCustomHelix *line_twoD = new TEveMu2eCustomHelix();
          
          std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator trajectoryIter;
          unsigned int g = 0;
          std::cout<<"Full 3 "<<std::endl;
          for(trajectoryIter=trajcol->begin(); trajectoryIter!=trajcol->end(); trajectoryIter++)
          {
            //if(!Contains(particleIds_,abs(trajectoryIter->first->pdgId()))) { continue;}
            const std::vector<MCTrajectoryPoint> &points = trajectoryIter->second.points();
            for(unsigned int i=0; i<points.size();i++){

              CLHEP::Hep3Vector Pos(points[i].x(), points[i].y(), points[i].z());
              GeomHandle<DetectorSystem> det;
              CLHEP::Hep3Vector HitPos2D = det->toDetector(Pos);

              if(i==0) {
                    line->SetPoint(i,pointmmTocm(Pos.x()), pointmmTocm(Pos.y()),pointmmTocm(Pos.z()));
                    line_twoD->SetPoint(i,pointmmTocm(HitPos2D.x()), pointmmTocm(HitPos2D.y()),pointmmTocm(HitPos2D.z()));
                } else {
                    line->SetNextPoint(pointmmTocm(Pos.x()), pointmmTocm(Pos.y()),pointmmTocm(Pos.z()));
                    line_twoD->SetNextPoint(pointmmTocm(HitPos2D.x()), pointmmTocm(HitPos2D.y()),pointmmTocm(HitPos2D.z()));
              }
            }

            string energy = to_string(points[0].kineticEnergy());
            //string pdgId= to_string(trajectoryIter->first->pdgId());
            line_twoD->SetLineColor(kRed);
            line_twoD->SetLineWidth(3);
            fTrackList2D->AddElement(line_twoD);

            line->SetPickable(kTRUE);
            const std::string title = " MCTrajectory # "+to_string(g)+" Energy = "+ energy;//+" PDG Code " + pdgId;
            line->SetTitle(Form(title.c_str()));
            if(abs(trajectoryIter->first->pdgId()) == 11) line->SetLineColor(kRed); //electrons
            else if(abs(trajectoryIter->first->pdgId()) == 13) line->SetLineColor(kGreen); //muons
            else if(abs(trajectoryIter->first->pdgId()) == 2112) line->SetLineColor(kYellow); //neutron  
            else if(abs(trajectoryIter->first->pdgId()) == 2212) line->SetLineColor(kRed); //proton
            else line->SetLineColor(kCyan); //nuceli
            line->SetLineWidth(3);
            fTrackList3D->AddElement(line);

            g++;
          }
          //TXYMgr->ImportElements(fTrackList2D, scene1);
          //TRZMgr->ImportElements(fTrackList2D, scene2);
          tracker2Dproj->fXYMgr->ImportElements(fTrackList2D, tracker2Dproj->fDetXYScene); 
          tracker2Dproj->fRZMgr->ImportElements(fTrackList2D, tracker2Dproj->fDetRZScene);

          gEve->AddElement(fTrackList3D);
          gEve->Redraw3D(kTRUE);
               
      }
    }*/
    
      
    void TEveMu2eMCInterface::AddFullMCTrajectory(bool firstloop, const MCTrajectoryCollection *trajcol, TEveMu2e2DProjection *tracker2Dproj, bool Redraw, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2, std::vector<int> particleIds){
	DataLists<const MCTrajectoryCollection*, TEveMu2e2DProjection*>(trajcol, Redraw, accumulate, "MC Trajectory", &fTrackList3D, &fTrackList2D, tracker2Dproj);

        if(trajcol!=0){
          TEveMu2eCustomHelix *line_twoD = new TEveMu2eCustomHelix();
      
          std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator trajectoryIter;
          for(trajectoryIter=trajcol->begin(); trajectoryIter!=trajcol->end(); trajectoryIter++)
          { 
            TEveMu2eCustomHelix *line = new TEveMu2eCustomHelix();
            
            int x = Contains(particleIds,trajectoryIter->first->pdgId()); 
            if(x == 1){
              const std::vector<MCTrajectoryPoint> &points = trajectoryIter->second.points();
              for(unsigned int i=0; i<points.size();i++){
                CLHEP::Hep3Vector Pos(points[i].x(), points[i].y(), points[i].z());
                GeomHandle<DetectorSystem> det;
                CLHEP::Hep3Vector HitPos2D = det->toDetector(Pos);

                if(i==0) {
                      line->SetPoint(i,pointmmTocm(Pos.x()), pointmmTocm(Pos.y()),pointmmTocm(Pos.z()));
                      line_twoD->SetPoint(i,pointmmTocm(HitPos2D.x()), pointmmTocm(HitPos2D.y()),pointmmTocm(HitPos2D.z()));
                } else {
                    line->SetNextPoint(pointmmTocm(Pos.x()), pointmmTocm(Pos.y()),pointmmTocm(Pos.z()));
                    line_twoD->SetNextPoint(pointmmTocm(HitPos2D.x()), pointmmTocm(HitPos2D.y()),pointmmTocm(HitPos2D.z()));
                }
              }
              
              string energy = to_string(points[0].kineticEnergy());
              //string pdgId= to_string(trajectoryIter->first->pdgId());
              const char* pid = "pid";
              auto t = new TEveText(pid);
              t->SetFontSize(15);
              line_twoD->SetLineColor(kRed);
              line_twoD->SetLineWidth(3);
              fTrackList2D->AddElement(line_twoD);

              const std::string title = " MCTrajectory "+ energy;// + "PDG Code = " + pdgId;
              line->SetTitle(Form(title.c_str()));
              if(trajectoryIter->first->pdgId() == 11){
                      line->SetLineColor(kRed);//electrons
                      pid = "electron -";
                      t->SetText(pid);
                      t->SetMainColor(kRed);
                      t->RefMainTrans().SetPos(0.0,400.0,3000.0);
              }
              else if(trajectoryIter->first->pdgId() == -11){
                           line->SetLineColor(kYellow);//positrons
                           pid = "positron +";
                           t->SetText(pid);
                           t->SetMainColor(kYellow);
                           t->RefMainTrans().SetPos(0.0,550.0,3000.0);
              }
              else if(trajectoryIter->first->pdgId() == 13){
                           line->SetLineColor(kGreen);//muons
                           pid = "muon - ";
                           t->SetText(pid);
                           t->SetMainColor(kGreen);
                           t->RefMainTrans().SetPos(0.0,650.0,3000.0);
              }
              else if(trajectoryIter->first->pdgId() == -211){
                           line->SetLineColor(kMagenta);//pion
                           pid = "pion -";
                           t->SetText(pid);
                           t->SetMainColor(kMagenta);
                           t->RefMainTrans().SetPos(0.0,750.0,3000.0);
              }
              else if(trajectoryIter->first->pdgId() == 2212){
                           line->SetLineColor(kBlue); //proton
                           pid = "proton";
                           t->SetText(pid);
                           t->SetMainColor(kBlue);
                           t->RefMainTrans().SetPos(0.0,850.0,3000.0);
              }
              else line->SetLineColor(kCyan); //nuceli
              line->SetLineWidth(3);
              fTrackList3D->AddElement(line);
              fTrackList3D->AddElement(t);
              }
             else std::cout<<"Not a user selected particle"<<std::endl;
            }
            tracker2Dproj->fXYMgr->ImportElements(fTrackList2D, tracker2Dproj->fDetXYScene); 
            tracker2Dproj->fRZMgr->ImportElements(fTrackList2D, tracker2Dproj->fDetRZScene);
            gEve->AddElement(fTrackList3D);
            gEve->Redraw3D(kTRUE);
        }
    }
}

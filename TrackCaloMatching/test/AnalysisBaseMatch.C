// Analysis match dumps the data to train the regression to predict the cluster position
// There are three modes available (uncomment the one you want below)
//   - the position of the 6 most energetic crystals (w.r.t. most energetic one) in XY coordinates
//   - the position of the 6 most energetic crystals (w.r.t. most energetic one) in UV coordinates
//   - the energy of the neighboring crystals (2 layers), using -1 to indicate if there is no crystal
//
// Note: the XY and UV distributions are also attainable by performing the coordinate transform in the TMVAReadRegress.C matrix
//   

#define AnalysisBaseMatch_cxx
#include "AnalysisBaseMatch.h"
#include <iostream>
#include <iomanip>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TH2Poly.h>
#include <TVector3.h>
#include <TAxis.h>
#include <TTree.h>
#include <TBranch.h>
#include "Src/Disk.hh"

#include <vector>
#include <cmath>
#include <set>
#include <list>
#include <fstream>
#include <algorithm>


int Nplot=0;

TH1F *_sener,*_distOut,*_dpos,*_phi,*_theta,*_output,*_cryDist,*_dClu,*_dphi,*_dtheta,*_dtrk,*_hx,*_hy,*_ha,*_div;
TH2F *_targ,*_duv,*_view[10];


TTree* dTuple;
TBranch *b_e0,*b_e1,*b_e2,*b_e3,*b_e4,*b_e5,*b_e6,*b_e7,*b_e8,*b_e9;
TBranch *b_e10,*b_e11,*b_e12,*b_e13,*b_e14,*b_e15,*b_e16,*b_e17,*b_e18,*b_e19;
TBranch *b_e20,*b_e21,*b_e22,*b_e23,*b_e24;
TBranch *b_t0,*b_t1,*b_r0,*b_z0,*b_ip0,*b_ip1,*b_c0,*b_c1,*b_c2,*b_c3,*b_c4,*b_c5,*b_ffx,*b_ffy,*b_ffz,*b_vdx,*b_vdy,*b_vdz;
TBranch *b_targetX,*b_targetY,*b_targetA;

Float_t e0_,e1_,e2_,e3_,e4_,e5_,e6_,e7_,e8_,e9_;
Float_t e10_,e11_,e12_,e13_,e14_,e15_,e16_,e17_,e18_,e19_;
Float_t e20_,e21_,e22_,e23_,e24_;
Float_t t0_,t1_,r0_,z0_,ip0_,ip1_,c0_,c1_,c2_,c3_,c4_,c5_,ffx_,ffy_,ffz_,vdx_,vdy_,vdz_;
Float_t targetX_,targetY_,targetA_;




double _cellSize = 34.3;
bool is_train    = false;
int mode         = 3;





void AnalysisBaseMatch::Begin(TTree * /*tree*/)
{

   TString option = GetOption();
   if (option.Sizeof()>=2) histname=option;   
   fHistFile = new TFile(histname,"RECREATE");
   cout<<endl<<"Histname "<<histname<<endl;
 
   _sener    = Book1DF("sener",       "Rcog sig",          100,      0,    1, "E",      "MeV");
   _distOut  = Book1DF("distOut",     "Dist out",          100,      0, 1000, "d",      "mm");
   _dpos     = Book1DF("dpos",        "Delta pos",         100,   -200,  200, "dx",     "mm");
   _phi      = Book1DF("phi",         "phi",               100,      0,  6.3, "phi",    "");
   _theta    = Book1DF("theta",       "theta",             100,      0,  3.2, "theta",  "");
   _cryDist  = Book1DF("cryDist",     "cryDist",           100,      0,  100, "dist",   "mm");
   _dClu     = Book1DF("dCLu",        "cryDist",           100,      0,  500, "dist",   "mm");
   _dphi     = Book1DF("dphi",        "cryDist",           100,   -0.5,  0.5, "#Delta#phi",  "");
   _dtheta   = Book1DF("dtheta",      "cryDist",           100,   -0.5,  0.5, "#Delta#theta","");
   _dtrk     = Book1DF("dtrk",        "cryDist",           100,    -30,   30, "#DeltaX","");
   _hx       = Book1DF("hx",          "X residual",         50,    -30,   30, "dx",     "mm");
   _hy       = Book1DF("hy",          "Y residual",         50,    -30,   30, "dy",     "mm");
   _ha       = Book1DF("ha",          "alpha",             200,      0,  200, "alpha",  "mm");
   _div      = Book1DF("div",         "div",                50,      0,  200, "div",    "mm");
   _targ     = new TH2F("targ",       "targ",              100,     -2,    2, 100, -2, 2);
   _duv      = new TH2F("duv",        "targ uv",           100,     -2,    2, 100, -2, 2);
    
   for (int i=0;i<10;++i) _view[i] = new TH2F(Form("view_%i",i),"view", 50, -5, 5, 50, -5, 5);

   disk = new Disk(1,374,660,34.300); 


   dTuple = new TTree("TreeReg","TreeReg");
   b_e0 = dTuple->Branch("e0", &e0_, "e0/F");  
   b_e1 = dTuple->Branch("e1", &e1_, "e1/F");  
   b_e2 = dTuple->Branch("e2", &e2_, "e2/F");  
   b_e3 = dTuple->Branch("e3", &e3_, "e3/F");  
   b_e4 = dTuple->Branch("e4", &e4_, "e4/F");  
   b_e5 = dTuple->Branch("e5", &e5_, "e5/F");  
   b_e6 = dTuple->Branch("e6", &e6_, "e6/F");  
   b_e7 = dTuple->Branch("e7", &e7_, "e7/F");  
   b_e8 = dTuple->Branch("e8", &e8_, "e8/F");  
   b_e9 = dTuple->Branch("e9", &e9_, "e9/F");  

   b_e10 = dTuple->Branch("e10", &e10_, "e10/F");  
   b_e11 = dTuple->Branch("e11", &e11_, "e11/F");  
   b_e12 = dTuple->Branch("e12", &e12_, "e12/F");  
   b_e13 = dTuple->Branch("e13", &e13_, "e13/F");  
   b_e14 = dTuple->Branch("e14", &e14_, "e14/F");  
   b_e15 = dTuple->Branch("e15", &e15_, "e15/F");  
   b_e16 = dTuple->Branch("e16", &e16_, "e16/F");  
   b_e17 = dTuple->Branch("e17", &e17_, "e17/F");  
   b_e18 = dTuple->Branch("e18", &e18_, "e18/F");  
   b_e19 = dTuple->Branch("e19", &e19_, "e19/F");  

   b_e20 = dTuple->Branch("e20", &e20_, "e20/F");  
   b_e21 = dTuple->Branch("e21", &e21_, "e21/F");  
   b_e22 = dTuple->Branch("e22", &e22_, "e22/F");  
   b_e23 = dTuple->Branch("e23", &e23_, "e23/F");  
   b_e24 = dTuple->Branch("e24", &e24_, "e24/F");  


   b_t0 = dTuple->Branch("t0", &t0_, "t0/F");  
   b_t1 = dTuple->Branch("t1", &t1_, "t1/F");  

   b_r0 = dTuple->Branch("r0", &r0_, "r0/F");  
   b_z0 = dTuple->Branch("z0", &z0_, "z0/F");  
   b_ip0 = dTuple->Branch("ip0", &ip0_, "ip0/F");  
   b_ip1 = dTuple->Branch("ip1", &ip1_, "ip1/F");  

   b_c0 = dTuple->Branch("c0", &c0_, "c0/F");  
   b_c1 = dTuple->Branch("c1", &c1_, "c1/F");  
   b_c2 = dTuple->Branch("c2", &c2_, "c2/F");  
   b_c3 = dTuple->Branch("c3", &c3_, "c3/F");  
   b_c4 = dTuple->Branch("c4", &c4_, "c4/F");  
   b_c5 = dTuple->Branch("c5", &c5_, "c5/F");  

   b_ffx = dTuple->Branch("ffx", &ffx_, "ffx/F");  
   b_ffy = dTuple->Branch("ffy", &ffy_, "ffy/F");  
   b_ffz = dTuple->Branch("ffz", &ffz_, "ffz/F");  
   b_vdx = dTuple->Branch("vdx", &vdx_, "vdx/F");  
   b_vdy = dTuple->Branch("vdy", &vdy_, "vdy/F");  
   b_vdz = dTuple->Branch("vdz", &vdz_, "vdz/F");  

   b_targetX = dTuple->Branch("targetX", &targetX_, "targetX/F");  
   b_targetY = dTuple->Branch("targetY", &targetY_, "targetY/F");  
   b_targetA = dTuple->Branch("targetA", &targetA_, "targetA/F");  
   
   is_train=false;
   if (histname.Contains("Train")) is_train = true;
   
   Nplot = 0;

}

      
      
      

      
      

Bool_t AnalysisBaseMatch::Process(Long64_t entry)
{
     
     fChain->GetTree()->GetEntry(entry);
          
     double vdPhi(0),vdTheta(0),vdmom(-1);
     
     int VDTrk = vdFinder(); 
     if (VDTrk > -1)
     {
         vdmom   = sqrt(vdMomX[VDTrk]*vdMomX[VDTrk]+vdMomY[VDTrk]*vdMomY[VDTrk]+vdMomZ[VDTrk]*vdMomZ[VDTrk]);
         vdPhi   = (atan2(vdMomY[VDTrk],vdMomX[VDTrk]));
         vdTheta = acos(vdMomZ[VDTrk]/vdmom);
         _phi->Fill(vdPhi);
         _theta->Fill(vdTheta);
     }


     //find a track that is matched to a conversion cluster
     int itrkMatch(-1);
     for (int im=0;im<nMatch;++im)
     {
          int ic = mCluId[im];
          if (cluConv[ic]!=1) continue; 
          itrkMatch = mTrkId[im];
     } 


     if (itrkMatch == -1)  return kTRUE; 
     if (trkprob[itrkMatch] < 1e-3) return kTRUE;
     if (trkStat[itrkMatch] !=1) return kTRUE;
     
     double trkp = sqrt(trkpX[itrkMatch]*trkpX[itrkMatch]+trkpY[itrkMatch]*trkpY[itrkMatch]+trkpZ[itrkMatch]*trkpZ[itrkMatch]);

     double trkPhi   = (atan2(trkpY[itrkMatch],trkpX[itrkMatch]));
     double trkTheta = acos(trkpZ[itrkMatch]/trkp);
          
     if (VDTrk > -1) _dphi->Fill(vdPhi-trkPhi);
     if (VDTrk > -1) _dtheta->Fill(vdTheta-trkTheta);
     


     
     
     for (int ic=0;ic<nCluster;++ic) 
     {     
	  if (cluConv[ic]!=1)     continue;	
          
          double dx0 = cluCogX[ic]-trkFFX[itrkMatch];
          double dy0 = cluCogY[ic]-trkFFY[itrkMatch];
          if (sqrt(dx0*dx0+dy0*dy0)>200) continue; 
          
	  	  	  
	  std::vector<int> crystalL = (*cluList)[ic];
	  int seedId     = cryId[crystalL[0]];
	  double centerX = cryPosX[crystalL[0]];
	  double centerY = cryPosY[crystalL[0]];
	  if (seedId >= disk->nCrystal()) seedId -= disk->nCrystal();


          int ivm(-1);
          for (int iv=0;iv<nVd;++iv) if (vdMom[iv] > 80 && (vdPdgId[iv]==11 ||vdPdgId[iv]==13) ) ivm=iv;
          
          int ism(-1);
          for (int is=cluSimIdx[ic];is<cluSimIdx[ic]+cluSimLen[ic];++is)if (clusimMom[is]>90) ism = is;


          //if (sqrt(centerX*centerX+centerY*centerY)>450) continue;          
          //if (trkFFZ[itrkMatch]>1) continue; 
          //std::cout<<entry<<" "<<nCluster<<" "<<ic<<" "<<centerX<<" "<<centerY<<" "<<endl;
	  	  

	  double matrixX[10]={0};
	  double matrixY[10]={0};
	  double matrixE[10]={0};
	  double etot25b(0);
	  for (unsigned int ii=0; ii<crystalL.size() && ii <10; ++ii)
          {	  
	      int cryIdx  = crystalL[ii];
	      matrixX[ii] = (cryPosX[cryIdx]-centerX)/_cellSize;
	      matrixY[ii] = (cryPosY[cryIdx]-centerY)/_cellSize;
	      matrixE[ii] = cryEdep[cryIdx]/105.0;
	      etot25b += cryEdep[cryIdx];
              //cout<<"- "<<cryPosX[cryIdx]<<" "<<cryPosY[cryIdx]<<" "<<cryEdep[cryIdx]<<endl;
	  }
	  _sener->Fill(etot25b/cluEnergy[ic]);
	  
          
          TVector3 cogNew = calcCog(ic);  

	  
 	  double targetX =(trkFFX[itrkMatch]-centerX)/_cellSize;
 	  double targetY =(trkFFY[itrkMatch]-centerY)/_cellSize;
	  
          TVector2 dxy(targetX,targetY); 
          TVector2 du(cos(trkPhi),sin(trkPhi));
          TVector2 dv(-sin(trkPhi),cos(trkPhi));
          TVector2 deltaUV0(dxy*du,dxy*dv);
          
          
          _targ->Fill(targetX,targetY);
          _duv->Fill(deltaUV0.X(),deltaUV0.Y());
          _dtrk->Fill(trkFFX[itrkMatch]-clusimPosX[ic]);
          _dtrk->Fill(trkFFY[itrkMatch]-clusimPosY[ic]);

          int miv(-1);
          for (int iv=0;iv<nVd;iv++) if (vdMom[iv]> 90 && (vdPdgId[iv]==11 || vdPdgId[iv]==13 )) miv=iv;
          if (ivm==-1) return kTRUE;
          if (miv > -1) 
          {
             double divx = vdPosX[miv]-trkFFX[itrkMatch];
             double divy = vdPosY[miv]-trkFFY[itrkMatch];
             double div  = sqrt(divx*divx+divy*divy);

             
             _div->Fill(div);
             if (div > 50) return kTRUE;
          }



	  //this fills the position of the 7 most energetic crystals (w.r.t. most energetic one) in XY coordinates  
	  if (mode==1)
          {
             e0_  = matrixX[0];
	     e1_  = matrixY[0];
	     e2_  = matrixE[0];
	     e3_  = matrixX[1];
	     e4_  = matrixY[1];
	     e5_  = matrixE[1];
	     e6_  = matrixX[2];
	     e7_  = matrixY[2];
	     e8_  = matrixE[2];
	     e9_  = matrixX[3];
	     e10_ = matrixY[3];
	     e11_ = matrixE[3];
	     e12_ = matrixX[4];
	     e13_ = matrixY[4];
	     e14_ = matrixE[4];
	     e15_ = matrixX[5];
	     e16_ = matrixY[5];
	     e17_ = matrixE[5];
	     e18_ = matrixX[6];
	     e19_ = matrixY[6];
	     e20_ = matrixE[6];

             targetX_ = targetX;
             targetY_ = targetY;
             targetA_ = targetExtend(ic,itrkMatch);

	     t0_ = trkPhi;
	     t1_ = trkTheta;
             r0_ = sqrt(centerX*centerX+centerY*centerY)/1000;
             z0_ = trkFFZ[itrkMatch]/200;
             ip0_ = trkCposX[itrkMatch]/10;
             ip1_ = trkCposY[itrkMatch]/10;

	     c0_  = cluCogX[ic];
	     c1_  = cluCogY[ic]; 	     
             c2_  = cogNew.X();
	     c3_  = cogNew.Y();
             c4_  = clusimPosX[ic];
	     c5_  = clusimPosY[ic];
             ffx_ = trkFFX[itrkMatch];
             ffy_ = trkFFY[itrkMatch];
             ffz_ = trkFFZ[itrkMatch];
             vdx_ = (ivm>-1) ? vdPosX[ivm] : -999;
             vdy_ = (ivm>-1) ? vdPosY[ivm] : -999;
             vdz_ = (ivm>-1) ? vdPosZ[ivm] : -999;
             
          }


	  //this fills the position of the 7 most energetic crystals (w.r.t. most energetic one) in UV coordinates  
	  if (mode==2)
          {
              e0_  = matrixX[0];
	      e1_  = matrixY[0];
	      e2_  = matrixE[0];
	      e3_  = cos(trkPhi)*matrixX[1]+sin(trkPhi)*matrixY[1];
	      e4_  = -sin(trkPhi)*matrixX[1]+cos(trkPhi)*matrixY[1];
	      e5_  = matrixE[1];
	      e6_  = cos(trkPhi)*matrixX[2]+sin(trkPhi)*matrixY[2];
	      e7_  = -sin(trkPhi)*matrixX[2]+cos(trkPhi)*matrixY[2];
	      e8_  = matrixE[2];
	      e9_  = cos(trkPhi)*matrixX[3]+sin(trkPhi)*matrixY[3];
	      e10_ = -sin(trkPhi)*matrixX[3]+cos(trkPhi)*matrixY[3];
	      e11_ = matrixE[3];
	      e12_ = cos(trkPhi)*matrixX[4]+sin(trkPhi)*matrixY[4];
	      e13_ = -sin(trkPhi)*matrixX[4]+cos(trkPhi)*matrixY[4];
	      e14_ = matrixE[4];
	      e15_ = cos(trkPhi)*matrixX[5]+sin(trkPhi)*matrixY[5];
	      e16_ = -sin(trkPhi)*matrixX[5]+cos(trkPhi)*matrixY[5];
	      e17_ = matrixE[5];
	      e18_ = cos(trkPhi)*matrixX[6]+sin(trkPhi)*matrixY[6];
	      e19_ = -sin(trkPhi)*matrixX[6]+cos(trkPhi)*matrixY[6];
	      e20_ = matrixE[6];

              targetX_ = deltaUV0.X();
              targetY_ = deltaUV0.Y();
              targetA_ = targetExtend(ic,itrkMatch);

	      t0_ = trkPhi;
	      t1_ = trkTheta;
              r0_ = sqrt(centerX*centerX+centerY*centerY)/1000;
              z0_ = trkFFZ[itrkMatch]/200;
              ip0_ = trkCposX[itrkMatch]/10;
              ip1_ = trkCposY[itrkMatch]/10;

	      c0_ =  cos(trkPhi)*cluCogX[ic] + sin(trkPhi)*cluCogY[ic];
	      c1_ = -sin(trkPhi)*cluCogX[ic] + cos(trkPhi)*cluCogY[ic];
              c2_ =  cos(trkPhi)*cogNew.X() + sin(trkPhi)*cogNew.Y();
	      c3_ = -sin(trkPhi)*cogNew.X() + cos(trkPhi)*cogNew.Y();
 	      c4_ =  cos(trkPhi)*clusimPosX[ic] + sin(trkPhi)*clusimPosY[ic];
	      c5_ = -sin(trkPhi)*clusimPosX[ic] + cos(trkPhi)*clusimPosY[ic];

              ffx_ = trkFFX[itrkMatch];
              ffy_ = trkFFY[itrkMatch];
              ffz_ = trkFFZ[itrkMatch];
              vdx_ = (ivm>-1) ? vdPosX[ivm] : -999;
              vdy_ = (ivm>-1) ? vdPosY[ivm] : -999;
              vdz_ = (ivm>-1) ? vdPosZ[ivm] : -999;
          }
          
          
          //here we take the two rings around the most energetic         
	  if (mode==3)
          {
              int offset(0);
              if (cryId[crystalL[0]]>=disk->nCrystal()) offset=disk->nCrystal();

              std::vector<int> neighbors1 = disk->findLocalNeighbors(cryId[crystalL[0]]-offset,1);
              std::vector<int> neighbors2 = disk->findLocalNeighbors(cryId[crystalL[0]]-offset,2);
              neighbors1.insert(neighbors1.end(), neighbors2.begin(), neighbors2.end());

              double etot19(cryEdep[crystalL[0]]);
              std::vector<double> evec;
              evec.push_back(cryEdep[crystalL[0]]);

              for (auto in : neighbors1)
              {
                 if (in == -1){evec.push_back(-1);continue;}

                 double et(0);
                 for (auto icr : crystalL) if (in==cryId[icr]-offset) et=cryEdep[icr];
                 evec.push_back(et);
                 etot19+=et;
              }

	      //this fills the energy of the central / first layer / second layer of crystals
	      e0_  = evec[0]/100;
	      e1_  = evec[1]/50;
	      e2_  = evec[2]/50;
	      e3_  = evec[3]/50;
	      e4_  = evec[4]/50;
	      e5_  = evec[5]/50;
	      e6_  = evec[6]/50;
	      e7_  = evec[7];
	      e8_  = evec[8];
	      e9_  = evec[9];
	      e10_ = evec[10];
	      e11_ = evec[11];
	      e12_ = evec[12];
	      e13_ = evec[13];
	      e14_ = evec[14];
	      e15_ = evec[15];
	      e16_ = evec[16];
	      e17_ = evec[17];
	      e18_ = evec[18];
	      e19_ = (cluEnergy[ic]-etot19)/100;
	      e20_ = 0;

              targetX_ = targetX;
              targetY_ = targetY;
              //targetX_ = deltaUV0.X();
              //targetY_ = deltaUV0.Y();
              targetA_ = targetExtend(ic,itrkMatch)/100;

	      t0_ = trkPhi;
	      t1_ = trkTheta;
              r0_ = sqrt(centerX*centerX+centerY*centerY)/1000;
              z0_ = trkFFZ[itrkMatch]/200;
              ip0_ = trkCposX[itrkMatch]/10;
              ip1_ = trkCposY[itrkMatch]/10;
              //ip0_ = sqrt(trkCposX[itrkMatch]*trkCposX[itrkMatch]+trkCposY[itrkMatch]*trkCposY[itrkMatch])/10;
              //ip0_ = std::max(abs(trkCposX[itrkMatch]),abs(trkCposY[itrkMatch]))/10;
              //ip0_ = std::max(abs(trkCposX[itrkMatch]),abs(trkCposY[itrkMatch]))> 10 ? 1 : 0;
              

	      c0_  = cluCogX[ic];
	      c1_  = cluCogY[ic];
              c2_  = cogNew.X();
	      c3_  = cogNew.Y();
 	      c4_  = clusimPosX[ic];
	      c5_  = clusimPosY[ic];              
              ffx_ = trkFFX[itrkMatch];
              ffy_ = trkFFY[itrkMatch];
              ffz_ = trkFFZ[itrkMatch]; 
              vdx_ = (ivm>-1) ? vdPosX[ivm] : -999;
              vdy_ = (ivm>-1) ? vdPosY[ivm] : -999;
              vdz_ = (ivm>-1) ? vdPosZ[ivm] : -999;
          }

          
          if (is_train)	  
	     {if (abs(targetX_)<2 && abs(targetY_)<2) dTuple->Fill();}           
	  else 
 	     {dTuple->Fill();}

	  
	  if ( (mode==1 || mode==2) && Nplot<10)
          {
             _view[Nplot]->Fill(e0_,e1_,e2_);
             _view[Nplot]->Fill(e3_,e4_,e5_);
             _view[Nplot]->Fill(e6_,e7_,e8_);
             _view[Nplot]->Fill(e9_,e10_,e11_);
             _view[Nplot]->Fill(e12_,e13_,e14_);
             _view[Nplot]->Fill(e15_,e16_,e17_);
             _view[Nplot]->Fill(e18_,e19_,e20_);
             _view[Nplot]->Fill(e21_,e22_,e23_);
             ++Nplot;          
          }
	  
	  
        
     }
     


     return kTRUE;
}




double AnalysisBaseMatch::targetExtend(int ic, int it)
{
     double trkp  = sqrt(trkpX[it]*trkpX[it] + trkpY[it]*trkpY[it]+ trkpZ[it]*trkpZ[it]);
     double trkpt = sqrt(trkpX[it]*trkpX[it]+trkpY[it]*trkpY[it]);
     
     TVector3 trkPos0(trkFFX[it],trkFFY[it],trkFFZ[it]);
     TVector3 trkDir(trkpX[it]/trkp,trkpY[it]/trkp,trkpZ[it]/trkp);

     TVector3 cog = calcCog(ic);
     
     double dlenT(0);
     TVector3 diff = trkPos0-cog ;
     double distMin = sqrt(diff.X()*diff.X()+diff.Y()*diff.Y());
     for (double dlent=0;dlent<400;dlent+=1)
     {
	  //TVector3 trkPos = trkPos0+dlent*trkDir;
	  TVector3 trkPos = getPosTrk(it,dlent);
                    
	  diff = trkPos-cog;
          double dist = sqrt(diff.X()*diff.X()+diff.Y()*diff.Y());
          if (dist < distMin) {distMin=dist;dlenT=dlent;}
          if (trkPos.Z() > 200) break; 
     }
     
     TVector3 TrkPosMin =  trkPos0+dlenT*trkDir;
     double dx = cog.X() - TrkPosMin.X();
     double dy = cog.Y() - TrkPosMin.Y();
         
     _hx->Fill(dx);
     _hy->Fill(dy);
     _ha->Fill(dlenT);
     
     return dlenT;
}


TVector3 AnalysisBaseMatch::getPosTrk(int it, double dlen)
{
    double sinDip  = sqrt(1-trkcdip[it]*trkcdip[it]);
    double length  = (trkZ[it] - trkz0[it])/sinDip + dlen;
    double posXTrk =  sin(trkphi0[it] + trkomega[it]*trkcdip[it]*length)/trkomega[it] - (trkd0[it] + 1.0/trkomega[it])*sin(trkphi0[it]);
    double posYTrk = -cos(trkphi0[it] + trkomega[it]*trkcdip[it]*length)/trkomega[it] + (trkd0[it] + 1.0/trkomega[it])*cos(trkphi0[it]);
    double posZtrk = dlen/sinDip+trkFFZ[it];

    return TVector3(posXTrk,posYTrk,posZtrk);
}


TVector3 AnalysisBaseMatch::calcCog(int ic)
{
      
      TVector3 aVector(0,0,0);
      double sumWeights(0);    

      std::vector<int> crystalL = (*cluList)[ic];
      for (unsigned int il=0;il<crystalL.size();++il) 
      {         

	int icry = crystalL[il];
	//double weight = -5.45 + 2.63*log(cryEdep[icry]);
        double weight = -4.934 + cryEdep[icry]; 
	if (weight < 0) weight = 0;

	aVector[0] += (cryPosX[icry])*weight;
	aVector[1] += cryPosY[icry]*weight;
	sumWeights += weight;
     }

     if (sumWeights>1e-3) { aVector[0] /= sumWeights; aVector[1] /= sumWeights;}
     else                 { aVector[0] = aVector[1] = 0;}

     return aVector;   
}




int AnalysisBaseMatch::vdFinder()
{
     int nv(0),vIdxMax(-1);
     double dmax(0);
     
     for (int i=0;i<nVd;++i)
     {       
       double vdmom = sqrt(vdMomX[i]*vdMomX[i]+vdMomY[i]*vdMomY[i]+vdMomZ[i]+vdMomZ[i]);       
       if (vdmom < 60  || vdGenIdx[i]!=0) continue;
       double dist =sqrt(vdPosX[i]*vdPosX[i]+vdPosY[i]*vdPosY[i]);	 
       if (dist > dmax) {dmax = dist; vIdxMax=i;} 
       ++nv;
     }
     if (nv==1 || dmax > 370) return vIdxMax;     
     return -1;
}

double AnalysisBaseMatch::cluDistance(int ic1, int ic2) 
{
    std::vector<int> crystalL1 = (*cluList)[ic1];
    std::vector<int> crystalL2 = (*cluList)[ic2];
    
    double dist(1e6);
    for (unsigned int i1 = 0;i1 < crystalL1.size();++i1)
    {
      double x1 = cryPosX[crystalL1[i1]];
      double y1 = cryPosY[crystalL1[i1]];
      for (unsigned int i2 = 0;i2 < crystalL2.size();++i2)
      {
        double x2 = cryPosX[crystalL2[i2]];
        double y2 = cryPosY[crystalL2[i2]];
        double dd = sqrt( (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
	if (dd <  dist) dist = dd;
      }    
    }
    return dist;
}    









void AnalysisBaseMatch::Terminate()
{
    fHistFile->cd();  
    fHistFile->Write();
    fHistFile->Close();  
}


void AnalysisBaseMatch::SetStyleGr(TGraph *gr)
{
    gr->SetTitle(" ");
    gr->GetXaxis()->SetTitle("Sig efficiency");
    gr->GetYaxis()->SetTitle("Event fraction");
    gr->GetXaxis()->SetTitleSize(0.04);
    gr->GetXaxis()->SetLabelSize(0.04);
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->GetYaxis()->SetLabelSize(0.04);
    gr->GetYaxis()->SetTitleOffset(1.2);
    gr->GetYaxis()->SetRangeUser(0,0.4);
}








TH1F* AnalysisBaseMatch::Book1DF(TString name, TString title, Int_t nbins, Double_t low, Double_t high, TString xaxis, TString units)
{
   
   TH1F* histo = new TH1F(name,title,nbins,low,high);
   SetStyle(histo,xaxis,units);
   return histo;
}

TH1I* AnalysisBaseMatch::Book1DI(TString name, TString title, Int_t nbins, Double_t low, Double_t high, TString xaxis, TString units)
{   
   TH1I* histo = new TH1I(name,title,nbins,low,high);
   SetStyle(histo,xaxis,units);
   return histo;
}

TH2F* AnalysisBaseMatch::Book2DF(TString name, TString title, Int_t nbinsx, Double_t xlow, Double_t xhigh, Int_t nbinsy, Double_t ylow, Double_t yhigh, TString xaxis, TString yaxis)
{
   TH2F* histo = new TH2F(name,title,nbinsx,xlow,xhigh,nbinsy,ylow,yhigh);
   histo->GetYaxis()->SetTitleOffset(1.2);
   histo->GetXaxis()->SetTitleSize(0.045);
   histo->GetYaxis()->SetTitleSize(0.045);
   histo->GetXaxis()->SetLabelSize(0.04);
   histo->GetYaxis()->SetLabelSize(0.04);
   histo->SetTitle("");
   histo->GetXaxis()->SetTitle(xaxis);     
   histo->GetYaxis()->SetTitle(yaxis);
   
   return histo;
}

void AnalysisBaseMatch::SetStyle(TH1* histo,TString xaxis, TString units)
{
   
   histo->SetMarkerStyle(20);
   histo->SetMarkerSize(0.6);
   histo->SetLineWidth(2);
   histo->GetYaxis()->SetTitleOffset(1.2);
   histo->GetXaxis()->SetTitleSize(0.045);
   histo->GetYaxis()->SetTitleSize(0.045);
   histo->GetXaxis()->SetLabelSize(0.04);
   histo->GetYaxis()->SetLabelSize(0.04);
   histo->SetTitle("");

   Char_t   ylabel[30];
   Double_t binsize = histo->GetBinWidth(1);
   sprintf(ylabel,"Entries / %g ",binsize);

   histo->GetXaxis()->SetTitle(xaxis+" ("+units+")"); 
   if (units=="")histo->GetXaxis()->SetTitle(xaxis);

   histo->GetYaxis()->SetTitle(ylabel+units);
   if (units=="") histo->GetYaxis()->SetTitle("Entries");

   return;
}
/*
	  double output = (clusimPosY[ic]-centerY)/_cellSize;
	  //double output = ((cluCogX[ic]-clusimPosX[ic]-3904)*cos(vdPhi)+(cluCogY[ic]-clusimPosY[ic])*sin(vdPhi))/_cellSize;	  	  
	  _output->Fill(output);
	  
	  vector<double> res;
	  
	  res.push_back(vdPhi);
	  res.push_back(vdTheta);
	  res.push_back(cluCogX[ic]);
	  res.push_back(cluCogY[ic]);
	  res.push_back(clusimPosX[ic]);
	  res.push_back(clusimPosY[ic]);
	  
	  res.push_back(matrixE[0]);
	  for (int i=1;i<10;++i) {res.push_back(matrixX[i]);res.push_back(matrixY[i]);res.push_back(matrixE[i]);}	  
  
	  res.push_back(vdPhi);
	  res.push_back(vdTheta);
	  res.push_back(output);	  
	  _results.push_back(res);
*/

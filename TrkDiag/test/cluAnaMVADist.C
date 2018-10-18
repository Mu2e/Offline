// Plot the MVA inpyt variables for 
//  bkg failing MVA - black
//  signal clusters passing MVA - red
//  signal clusters failing MVA - green
//
// In you analyzer, you must you must run
//     BD: @local::BD
//    TCD : @local::TCD
//and set 
// physics.producers.FlagBkgHits.SaveBkgClusters :true
// physics.analyzers.TCD.diagLevel :3
//
// you also need to add CaloTrig in the Analyzer module list

#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <vector>
#include <iostream>

int mode=1;


void cluAnaMVADist()
{
      
   gStyle->SetOptStat(0);
   
   
   const int Nvar = 12;
   TString htitle[Nvar]= {"HitRho","HitRhoSpread","ClusterRho","TimeSPread","Zmin","Zmax","ZGap","Nplanes","NexpectedPlanes","PlaneFrac","NplanesHit","nHits"};
   double  hmin[Nvar]  = {  0,  0,350,-2,-1500,-1500,   0, 0, 0,0,   0,  0};
   double  hmax[Nvar]  = {100,100,700,30, 1500, 1500,3000,40,40,1.1,15,100};
   
   TH1F *_var[Nvar][3];
   for (int i=0;i<Nvar;++i){   
    _var[i][0]  = new TH1F(Form("var%i_0",i),htitle[i],100,hmin[i],hmax[i]);
    _var[i][1]  = new TH1F(Form("var%i_1",i),htitle[i],100,hmin[i],hmax[i]);
    _var[i][2]  = new TH1F(Form("var%i_2",i),htitle[i],100,hmin[i],hmax[i]);
    _var[i][1]->SetLineColor(2);
    _var[i][2]->SetLineColor(8);
    _var[i][0]->SetMinimum(0.5);
    _var[i][1]->SetLineWidth(2);
    _var[i][2]->SetLineWidth(2);
   }
   
   TH1F* _mva0 = new TH1F("mvaout0","mvaout0",100,0,1);
   TH1F* _mva1 = new TH1F("mvaout1","mvaout1",100,0,1);
   _mva1->SetLineColor(2);
   _mva0->SetMinimum(0.5);

   TTreeReader fReader;  
   TTreeReaderValue<Int_t> iev      = {fReader, "iev"};
   TTreeReaderValue<Int_t> cluIdx   = {fReader, "cluIdx"};
   TTreeReaderValue<Int_t> nchits   = {fReader, "nchits"};
   TTreeReaderValue<Int_t> nshits   = {fReader, "nshits"};
   TTreeReaderValue<Float_t> mvaout = {fReader, "mvaout"};
   TTreeReaderValue<Float_t> pmom   = {fReader, "pmom"};
   TTreeReaderValue<Int_t> ppid     = {fReader, "ppid"};
   TTreeReaderValue<Int_t> ppdg     = {fReader, "ppdg"};
   TTreeReaderValue<Int_t> pgen     = {fReader, "pgen"};
   TTreeReaderValue<Int_t> pproc    = {fReader, "pproc"};
   TTreeReaderValue<Int_t> nprimary = {fReader, "nprimary"};
   TTreeReaderValue<Int_t> nconv    = {fReader, "nconv"};
   TTreeReaderValue<Int_t> ndelta   = {fReader, "ndelta"};
   TTreeReaderValue<Int_t> ncompt   = {fReader, "ncompt"};
   TTreeReaderValue<Int_t> ngconv   = {fReader, "ngconv"};
   TTreeReaderValue<Int_t> nebkg    = {fReader, "nebkg"};
   TTreeReaderValue<Int_t> nprot    = {fReader, "nprot"};
   TTreeReaderValue<Int_t> ncontrib = {fReader, "ncontrib"};
   TTreeReaderArray<Int_t> icontrib = {fReader, "icontrib"};
   TTreeReaderValue<Int_t> nindex   = {fReader, "nindex"};
   TTreeReaderArray<Int_t> hitidx   = {fReader, "hitidx"};
   TTreeReaderValue<Float_t> HitRho          = {fReader, "HitRho"};
   TTreeReaderValue<Float_t> HitRhoSpread    = {fReader, "HitRhoSpread"};
   TTreeReaderValue<Float_t> ClusterRho      = {fReader, "ClusterRho"};
   TTreeReaderValue<Float_t> TimeSpread      = {fReader, "TimeSpread"};
   TTreeReaderValue<Float_t> ZMin            = {fReader, "ZMin"};
   TTreeReaderValue<Float_t> ZMax            = {fReader, "ZMax"};
   TTreeReaderValue<Float_t> ZGap            = {fReader, "ZGap"};
   TTreeReaderValue<Float_t> NPlanes         = {fReader, "NPlanes"};
   TTreeReaderValue<Float_t> NExpectedPlanes = {fReader, "NExpectedPlanes"};
   TTreeReaderValue<Float_t> PlaneFraction   = {fReader, "PlaneFraction"};
   TTreeReaderValue<Float_t> NPlaneHits      = {fReader, "NPlaneHits"};
   TTreeReaderValue<Float_t> NHits           = {fReader, "NHits"};

   TTreeReader fReader2;  
   TTreeReaderValue<Int_t>   iev2      = {fReader2, "iev"};
   TTreeReaderValue<Int_t>   nhits     = {fReader2, "nhits"};
   TTreeReaderArray<Int_t>   hitNcombo = {fReader2, "hitNcombo"};
   TTreeReaderArray<Float_t> hitRad    = {fReader2, "hitRad"};
   TTreeReaderArray<Float_t> hitPhi    = {fReader2, "hitPhi"};
   TTreeReaderArray<Float_t> hitTime   = {fReader2, "hitTime"};
   TTreeReaderArray<Int_t>   hitPdg    = {fReader2, "hitPdg"};
   TTreeReaderArray<Int_t>   hitCrCode = {fReader2, "hitCrCode"};
   TTreeReaderArray<Int_t>   hitGen    = {fReader2, "hitGen"};
   
   TTreeReader fReader3;  
   TTreeReaderValue<Int_t>   iev3      = {fReader3, "evt"};
   TTreeReaderValue<Int_t>   nCluster  = {fReader3, "nCluster"};
   TTreeReaderArray<Float_t> cluEnergy = {fReader3, "cluEnergy"};
   TTreeReaderArray<Float_t> cluTime   = {fReader3, "cluTime"};
   TTreeReaderArray<Float_t> cluCogX   = {fReader3, "cluCogX"};
   TTreeReaderArray<Float_t> cluCogY   = {fReader3, "cluCogY"};
 
   TFile f("trkDiag.root");
   
   TTree* tree = (TTree*) f.Get("BD/bkgdiag");   
   TTree* tree2 = (TTree*) f.Get("BD/bkgdiag2");   
   TTree* tree3 = (TTree*) f.Get("CaloTrig/Calo");   

   fReader.SetTree(tree);      
   fReader2.SetTree(tree2);      
   fReader3.SetTree(tree3);      
   
   
   int nTot(0),nPass(0),nPass2(0),nConvTot(0),nConvPass(0),nConvPass2(0);
   for (int entry = 1; entry < tree2->GetEntriesFast();++entry){
   
       //set fReader2 to correct entry
       int in2(0);
       for (in2=0;in2<tree2->GetEntriesFast();++in2){
           fReader2.SetEntry(in2);       
           if (*iev2 == entry) break;
       }

       if (*iev2 != entry) continue;

       fReader3.SetEntry(in2); 
       vector<double> caloCluTime,caloCluPhi;  
       for (size_t i=0;i<*nCluster;++i) {
         if (cluEnergy[i]>60) {caloCluTime.push_back(cluTime[i]);caloCluPhi.push_back(atan2(cluCogY[i],cluCogX[i]));}
       }
       
       vector<double> badCluTime,badCluPhi;  
       for (int in=0;in<tree->GetEntriesFast();++in){
          fReader.SetEntry(in);
          if (*iev != entry) continue;
          bool sigSel = (*NHits <25) && (*ZGap < 50 || *ZGap >800) && (*PlaneFraction< 0.41 || *PlaneFraction >0.95);
          if (sigSel && *mvaout >0.001) {
          badCluTime.push_back(hitTime[hitidx[0]]);
          badCluPhi.push_back(hitPhi[hitidx[0]]);
          }
       }
       
       

       int nC(0),nCP(0);
       bool show(false);
       for (int in=0;in<tree->GetEntriesFast();++in){
         fReader.SetEntry(in);
         if (*iev != entry) continue;
 
         bool isSel(false);
         bool sigSel = (*NHits <25) && (*ZGap < 50 || *ZGap >800) && (*PlaneFraction< 0.41 || *PlaneFraction >0.95);
         
         
         if (*nshits<3 || *NPlanes<2) {
           //nTot     += *nchits;
           //nConvTot += *nconv;
           //nC       += *nconv;
           //nConvPass += *nconv;
           //nPass += *nchits;
           //nCP += *nconv;         
           //nConvPass2 += *nconv;
           //nPass2 += *nchits;          
           continue; 
         }


         nTot     += *nchits;
         nConvTot += *nconv;
         nC       += *nconv;
         if (*mvaout < 0.5) nConvPass += *nconv;
         if (*mvaout < 0.5) nPass += *nchits;
         if (*mvaout < 0.5) nCP += *nconv;         

         //regular selection
         if (sigSel){ 
           nConvPass2 += *nconv;
           nPass2 += *nchits;
           isSel = true;
         }
         
         //recover hits based on track cluster 
         
         if (!sigSel && *NHits <20){
             int nCluOk(0);
             for (size_t i=0;i<badCluTime.size();++i){ 
               double dPhi = badCluPhi[i]-hitPhi[hitidx[0]];
               while (dPhi < -3.14159)dPhi +=6.283185;
               while (dPhi > 3.14159) dPhi -=6.283185;           
               if ((std::abs(badCluTime[i]-hitTime[hitidx[0]]) < 30) && abs(dPhi) < 1.2) ++nCluOk;
             }
             if (nCluOk>3){
                nConvPass2 += *nconv;
                nPass2 += *nchits;
                isSel = true;
            }          
         }
         
                 
         //recover hits based on calo cluster timing
         /*
         if (!sigSel && *NHits <20){
             int nCluOk(0);
             for (size_t i=0;i<caloCluPhi.size();++i){ 
               double dPhi = badCluPhi[i]-hitPhi[hitidx[0]];
               while (dPhi < -3.14159)dPhi +=6.283185;
               while (dPhi > 3.14159) dPhi -=6.283185;           
               if ((std::abs(caloCluTime[i]-hitTime[hitidx[0]]) < 40) && abs(dPhi) < 2.0) ++nCluOk;
             }
            
             if (nCluOk>0){
                nConvPass2 += *nconv;
                nPass2 += *nchits;
                isSel = true;
            }          
         }
         */
         
         
         
         
         if (*nconv>0 && !isSel) show = true;
         if (*nconv>0 && !isSel) cout<<"Failure for entry="<<entry<<"  "<<*mvaout<<"  "<<*nconv<<"/"<<*nshits<<endl;
         
         
         //uncomment to see effect of cuts on the plots
         //if ((*NHits >25) || (*ZGap > 50 && *ZGap <900) || (*PlaneFraction> 0.4 && *PlaneFraction <0.95)) continue;
         //if ((*NHits >25) ||  (*PlaneFraction> 0.4 && *PlaneFraction <0.95)) continue;
         
         if (*nconv==0) _mva0->Fill(*mvaout);
         if (*nconv>0)  _mva1->Fill(*mvaout);
         
         if (*nconv==0 && *mvaout > 0.5){
           _var[0][0]->Fill(*HitRho);
           _var[1][0]->Fill(*HitRhoSpread);
           _var[2][0]->Fill(*ClusterRho);
           _var[3][0]->Fill(*TimeSpread);
           _var[4][0]->Fill(*ZMin);
           _var[5][0]->Fill(*ZMax);
           _var[6][0]->Fill(*ZGap);
           _var[7][0]->Fill(*NPlanes);
           _var[8][0]->Fill(*NExpectedPlanes);
           _var[9][0]->Fill(*PlaneFraction);
           _var[10][0]->Fill(*NPlaneHits);
           _var[11][0]->Fill(*NHits);         
         } 
         
         if (*nconv>0 && *mvaout < 0.5) {
           _var[0][1]->Fill(*HitRho);
           _var[1][1]->Fill(*HitRhoSpread);
           _var[2][1]->Fill(*ClusterRho);
           _var[3][1]->Fill(*TimeSpread);
           _var[4][1]->Fill(*ZMin);
           _var[5][1]->Fill(*ZMax);
           _var[6][1]->Fill(*ZGap);
           _var[7][1]->Fill(*NPlanes);
           _var[8][1]->Fill(*NExpectedPlanes);
           _var[9][1]->Fill(*PlaneFraction);
           _var[10][1]->Fill(*NPlaneHits);
           _var[11][1]->Fill(*NHits);         
         }        
         
         if (*nconv>0 && *mvaout >= 0.5){
           _var[0][2]->Fill(*HitRho);
           _var[1][2]->Fill(*HitRhoSpread);
           _var[2][2]->Fill(*ClusterRho);
           _var[3][2]->Fill(*TimeSpread);
           _var[4][2]->Fill(*ZMin);
           _var[5][2]->Fill(*ZMax);
           _var[6][2]->Fill(*ZGap);
           _var[7][2]->Fill(*NPlanes);
           _var[8][2]->Fill(*NExpectedPlanes);
           _var[9][2]->Fill(*PlaneFraction);
           _var[10][2]->Fill(*NPlaneHits);
           _var[11][2]->Fill(*NHits);                    
         } 
         
        }
        if (show) std::cout<<"  Fraction of selected signal hits= "<<nCP<<"/"<<nC<<endl;
   }
    
   
   
   cout<<endl<<endl;
   cout<<"Fraction of signal hit selected by MVA="<<float(nConvPass)/float(nConvTot)<<endl;
   cout<<"Fraction of signal hit selected by CUT="<<float(nConvPass2)/float(nConvTot)<<endl;
   cout<<"Fraction of hits rejected by MVA="<<1.0-float(nPass)/float(nTot)<<endl;
   cout<<"Fraction of hits rejected by CUT="<<1.0-float(nPass2)/float(nTot)<<endl;
   cout<<endl;
     
   TCanvas *c1 = new TCanvas("c1","c1",1000,1000);   
   c1->SetLogy(1);
   c1->Divide(3,2);
   
   c1->cd(1)->SetLogy(1);   
   _var[0][0]->Draw();
   _var[0][1]->Draw("same");
   _var[0][2]->Draw("same");
   c1->cd(2)->SetLogy(1);;   
   _var[1][0]->Draw();
   _var[1][1]->Draw("same");
   _var[1][2]->Draw("same");
   c1->cd(3)->SetLogy(1);   
   _var[2][0]->Draw();
   _var[2][1]->Draw("same");
   _var[2][2]->Draw("same");
   c1->cd(4)->SetLogy(1);   
   _var[3][0]->Draw();
   _var[3][1]->Draw("same");
   _var[3][2]->Draw("same");
   c1->cd(5)->SetLogy(1);   
   _var[4][0]->Draw();
   _var[4][1]->Draw("same");
   _var[4][2]->Draw("same");
   c1->cd(6)->SetLogy(1);   
   _var[5][0]->Draw();
   _var[5][1]->Draw("same");
   _var[5][2]->Draw("same");
   c1->SaveAs("cluMvaAna1.pdf");
   
   c1->cd(1);   
   _var[6][0]->Draw();
   _var[6][1]->Draw("same");
   _var[6][2]->Draw("same");
   c1->cd(2);   
   _var[7][0]->Draw();
   _var[7][1]->Draw("same");
   _var[7][2]->Draw("same");
   c1->cd(3);   
   _var[8][0]->Draw();
   _var[8][1]->Draw("same");
   _var[8][2]->Draw("same");
   c1->cd(4);   
   _var[9][0]->Draw();
   _var[9][1]->Draw("same");
   _var[9][2]->Draw("same");
   c1->cd(5);   
   _var[10][0]->Draw();
   _var[10][1]->Draw("same");
   _var[10][2]->Draw("same");
   c1->cd(6);   
   _var[11][0]->Draw();
   _var[11][1]->Draw("same");
   _var[11][2]->Draw("same");
   c1->SaveAs("cluMvaAna2.pdf");

   TCanvas *c2 = new TCanvas("c2","c2",1000,1000);   
   c2->SetLogy();
   _mva0->Draw();
   _mva1->Draw("same");
   c2->SaveAs("cluMvaDist.pdf");

   f.Close();
}

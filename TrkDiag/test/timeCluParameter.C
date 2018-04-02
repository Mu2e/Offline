// Plot the variables used for the MVA in timeClusterFinder. 
// Must run without the cluster refinment: TrkPatRec/fcl/prolog.fcl -->RefineClusters : false
//
// In you analyzer, you must you must run
//     BD: @local::BD
//    TCD : @local::TCD
//and set 
// physics.producers.FlagBkgHits.SaveBkgClusters :true
// physics.analyzers.TCD.diagLevel :3
//

#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <vector>
#include <iostream>




void timeCluParameter()
{

   TH1F *_time[3],*_phi[3],*_rad[3];
   TH2F *_view[3];
   for (int i=0;i<3;++i)_time[i]  = new TH1F(Form("time%i",i),"delta time",100,1,51);
   for (int i=0;i<3;++i)_rad[i]   = new TH1F(Form("rad%i",i), "rad",100,350,700);
   for (int i=0;i<3;++i)_phi[i]   = new TH1F(Form("phi%i",i), "delta phi",100,0.,6.3);
   for (int i=0;i<3;++i)_view[i]  = new TH2F(Form("view%i",i), "delta phi vs rad",100,0.,6.3,100,350,700);

   for (int i=0;i<3;++i){   
     _time[i]->SetLineColor(i+1);
     _rad[i]->SetLineColor(i+1);
     _phi[i]->SetLineColor(i+1);
   }


   TTreeReader fReader;  
   TTreeReaderValue<Int_t>   iev        = {fReader, "iev"};
   TTreeReaderValue<Int_t>   tcnHits    = {fReader, "tcnHits"};
   TTreeReaderArray<Int_t>   tcHitIdx   = {fReader, "tcHitIdx"};
   TTreeReaderValue<Int_t>   tcover     = {fReader, "tcover"};
   TTreeReaderValue<Int_t>   tcnce      = {fReader, "tcnce"};
   TTreeReaderValue<Float_t> tctime     = {fReader, "tctime"};
   TTreeReaderValue<Float_t> tcphi      = {fReader, "tcphi"};
   TTreeReaderValue<Float_t> tcminhtime = {fReader, "tcminhtime"};
   TTreeReaderValue<Float_t> tcmaxhtime = {fReader, "tcmaxhtime"};

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
  

   TFile f("trkDiag.root");
   
   TTree* tree = (TTree*) f.Get("TCD/tcdiag2");   
   TTree* tree2 = (TTree*) f.Get("BD/bkgdiag2");   
   fReader.SetTree(tree);      
   fReader2.SetTree(tree2);      
   


   for (int entry=1;entry <600 ; ++entry){         

      for (int in2=0;in2<tree2->GetEntriesFast();++in2){
          fReader2.SetEntry(in2);       
          if (*iev2 == entry) break;
      }

      if (*iev2 != entry) continue;

      for (int in=0;in<tree->GetEntriesFast();++in) {
         fReader.SetEntry(in);
         if (*iev != entry) continue;

         for (int id=0;id<*tcnHits;++id){
            int ih = tcHitIdx[id];
            int ipass = (hitGen[ih]==2) ? 1 : 0;
            double dphi = hitPhi[ih]-*tcphi;
            if (dphi > 6.283185) dphi -= 6.283185;
            //if (abs(dphi) > 1.57) continue;
            _phi[ipass]->Fill(abs(dphi));
            _time[ipass]->Fill(abs(hitTime[ih]-*tctime));
            _rad[ipass]->Fill(hitRad[ih]);
            _view[ipass]->Fill(abs(dphi),hitRad[ih]);
         }
      }
   
   }
   
   _phi[0]->SetMinimum(0);
   _rad[0]->SetMinimum(0);
   _time[0]->SetMinimum(0);
   
   
   //_time[0]->Scale(_time[1]->GetSumOfWeights()/_time[0]->GetSumOfWeights());
   //_rad[0]->Scale(_rad[1]->GetSumOfWeights()/_rad[0]->GetSumOfWeights());
   
   TCanvas *c1 = new TCanvas("c1","c1",1000,1000);   
   c1->Divide(2,2);

   c1->cd(1);   
   _phi[0]->Draw("");
   _phi[1]->Draw("same");
   c1->cd(2);   
   _time[0]->Draw("");
   _time[1]->Draw("same");
   c1->cd(3);   
   _rad[0]->Draw("");
   _rad[1]->Draw("same");
   c1->cd(4);   
   _view[0]->Draw("colz");

   c1->SaveAs("tcCluParameter.pdf");
   f.Close();
}

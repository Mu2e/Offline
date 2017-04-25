#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"


using namespace TMVA;
   
void TMVAReadRegress(TString testDS, int mode=1) 
{

   if (mode!=1 && mode !=2 && mode!=3 ) {cout<<"Mose must be 1 or 2 or 3"<<endl; return;}
   
   TMVA::Tools::Instance();

   gStyle->SetFrameBorderMode(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetPadColor(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetStatColor(0);
   gStyle->SetTitleFillColor(0);
   gStyle->SetOptStat(0);
   gStyle->SetHistLineWidth(2);
   gStyle->SetLineWidth(2);
   gStyle->SetPadLeftMargin(0.12);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetPalette(1);
   gStyle->SetFrameLineWidth(2);
      
   Float_t e0,e1,e2,e3,e4,e5,e6,e7,e8,e9;
   Float_t e10,e11,e12,e13,e14,e15,e16,e17,e18,e19;
   Float_t e20,e21,e22,e23,e24;
   Float_t r0,z0,t0,t1,ip0,ip1,c0,c1,c2,c3,c4,c5,FFX,FFY,FFZ,vdx,vdy,vdz;
   Float_t targetX,targetY;


   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    
   if (mode==3)reader->AddVariable( "e0", &e0); 
   if (mode==3)reader->AddVariable( "e1", &e1); 
   reader->AddVariable( "e2",   &e2);
   reader->AddVariable( "e3",   &e3);
   reader->AddVariable( "e4",   &e4);
   reader->AddVariable( "e5",   &e5);
   reader->AddVariable( "e6",   &e6);
   reader->AddVariable( "e7",   &e7);
   reader->AddVariable( "e8",   &e8);
   reader->AddVariable( "e9",   &e9);
   reader->AddVariable( "e10",  &e10);
   reader->AddVariable( "e11",  &e11);
   reader->AddVariable( "e12",  &e12);
   reader->AddVariable( "e13",  &e13);
   reader->AddVariable( "e14",  &e14);
   reader->AddVariable( "e15",  &e15);
   reader->AddVariable( "e16",  &e16);
   reader->AddVariable( "e17",  &e17);
   reader->AddVariable( "e18",  &e18);
   reader->AddVariable( "e19",  &e19);
   if (mode!=3) reader->AddVariable( "e20",  &e20);
   reader->AddVariable( "t0",   &t0);
   reader->AddVariable( "t1",   &t1);
   reader->AddVariable( "r0",   &r0);
   reader->AddVariable( "z0",   &z0);
   //reader->AddVariable( "ip0",  &ip0);
   //reader->AddVariable( "ip1",  &ip1);
   reader->AddSpectator( "c0",  &c0);
   reader->AddSpectator( "c1",  &c1);
   reader->AddSpectator( "c2",  &c2);
   reader->AddSpectator( "c3",  &c3); 
   reader->AddSpectator( "c4",  &c4); 
   reader->AddSpectator( "c5",  &c5); 
   reader->AddSpectator( "ffx", &FFX); 
   reader->AddSpectator( "ffy", &FFY); 
   reader->AddSpectator( "ffz", &FFZ); 
   reader->AddSpectator( "vdx", &vdx); 
   reader->AddSpectator( "vdy", &vdy); 
   reader->AddSpectator( "vdz", &vdz); 
   reader->BookMVA("BDTG method", "weightsX/weights/TMVARegression_BDTG.weights.xml" ); 
   
   TMVA::Reader *reader2 = new TMVA::Reader( "!Color:!Silent" );    
   if (mode==3) reader2->AddVariable( "e0", &e0); 
   if (mode==3) reader2->AddVariable( "e1", &e1); 
   reader2->AddVariable( "e2",   &e2);
   reader2->AddVariable( "e3",   &e3);
   reader2->AddVariable( "e4",   &e4);
   reader2->AddVariable( "e5",   &e5);
   reader2->AddVariable( "e6",   &e6);
   reader2->AddVariable( "e7",   &e7);
   reader2->AddVariable( "e8",   &e8);
   reader2->AddVariable( "e9",   &e9);
   reader2->AddVariable( "e10",  &e10);
   reader2->AddVariable( "e11",  &e11);
   reader2->AddVariable( "e12",  &e12);
   reader2->AddVariable( "e13",  &e13);
   reader2->AddVariable( "e14",  &e14);
   reader2->AddVariable( "e15",  &e15);
   reader2->AddVariable( "e16",  &e16);
   reader2->AddVariable( "e17",  &e17);
   reader2->AddVariable( "e18",  &e18);
   reader2->AddVariable( "e19",  &e19);
   if (mode!=3) reader2->AddVariable( "e20",  &e20);
   reader2->AddVariable( "t0",   &t0);
   reader2->AddVariable( "t1",   &t1);
   reader2->AddVariable( "r0",   &r0);
   reader2->AddVariable( "z0",   &z0);
   //reader2->AddVariable( "ip0",  &ip0);
   //reader2->AddVariable( "ip1",  &ip1);
   reader2->AddSpectator( "c0",  &c0);
   reader2->AddSpectator( "c1",  &c1);
   reader2->AddSpectator( "c2",  &c2);
   reader2->AddSpectator( "c3",  &c3); 
   reader2->AddSpectator( "c4",  &c4); 
   reader2->AddSpectator( "c5",  &c5); 
   reader2->AddSpectator( "ffx", &FFX); 
   reader2->AddSpectator( "ffy", &FFY); 
   reader2->AddSpectator( "ffz", &FFZ); 
   reader2->AddSpectator( "vdx", &vdx); 
   reader2->AddSpectator( "vdy", &vdy); 
   reader2->AddSpectator( "vdz", &vdz); 
   reader2->BookMVA("BDTG method", "weightsY/weights/TMVARegression_BDTG.weights.xml" ); 



   TFile input(testDS);
   TTree *regTree = (TTree*)input.Get("TreeReg");

   regTree->SetBranchAddress( "e0", &e0 );
   regTree->SetBranchAddress( "e1", &e1 );

   regTree->SetBranchAddress( "e2", &e2 );
   regTree->SetBranchAddress( "e3", &e3 );
   regTree->SetBranchAddress( "e4", &e4 );
   regTree->SetBranchAddress( "e5", &e5 );
   regTree->SetBranchAddress( "e6", &e6 );
   regTree->SetBranchAddress( "e7", &e7 );
   regTree->SetBranchAddress( "e8", &e8 );
   regTree->SetBranchAddress( "e9", &e9 );
   regTree->SetBranchAddress( "e10", &e10 );
   regTree->SetBranchAddress( "e11", &e11 );
   regTree->SetBranchAddress( "e12", &e12 );
   regTree->SetBranchAddress( "e13", &e13 );
   regTree->SetBranchAddress( "e14", &e14 );
   regTree->SetBranchAddress( "e15", &e15 );
   regTree->SetBranchAddress( "e16", &e16 );
   regTree->SetBranchAddress( "e17", &e17 );
   regTree->SetBranchAddress( "e18", &e18 );
   regTree->SetBranchAddress( "e19", &e19 );
   regTree->SetBranchAddress( "e20", &e20 );
   regTree->SetBranchAddress( "t0", &t0 );
   regTree->SetBranchAddress( "t1", &t1 );
   regTree->SetBranchAddress( "r0", &r0 );
   regTree->SetBranchAddress( "z0", &z0 );
   regTree->SetBranchAddress( "ip0", &ip0 );
   regTree->SetBranchAddress( "ip1", &ip1 );
   regTree->SetBranchAddress( "c0", &c0 );
   regTree->SetBranchAddress( "c1", &c1 );
   regTree->SetBranchAddress( "c2", &c2 );
   regTree->SetBranchAddress( "c3", &c3 );
   regTree->SetBranchAddress( "c4", &c4 );
   regTree->SetBranchAddress( "c5", &c5 );
   regTree->SetBranchAddress( "ffx", &FFX );
   regTree->SetBranchAddress( "ffy", &FFY );
   regTree->SetBranchAddress( "ffz", &FFZ );
   regTree->SetBranchAddress( "vdx", &vdx );
   regTree->SetBranchAddress( "vdy", &vdy );
   regTree->SetBranchAddress( "vdz", &vdz );
   regTree->SetBranchAddress( "targetX", &targetX );
   regTree->SetBranchAddress( "targetY", &targetY );


   TH1F hh1("hh1",  "",100,-50,50);
   TH1F hh2("hh2","",100,-50,50);
   TH2F hh3("hh3","hh4",200,-50,50,200,-50,50);

   TH1F hh4("hh4",  "",100,-50,50);
   TH1F hh5("hh5","",100,-50,50);
   TH2F hh6("hh6","hh8",200,-50,50,200,-50,50);




   for (Long64_t ievt=0; ievt<regTree->GetEntries();ievt++){
      regTree->GetEntry(ievt);
      
      Float_t val1 = (reader->EvaluateRegression("BDTG method"))[0];
      Float_t val2 = (reader2->EvaluateRegression("BDTG method"))[0];

      //cout<<e2<<" "<<e3<<" "<<e4<<" "<<e5<<" "<<e6<<" "<<e7<<" "<<e8<<" "<<e9<<" "
      //    <<e10<<" "<<e11<<" "<<e12<<" "<<e13<<" "<<e14<<" "<<e15<<" "<<e16<<" "<<e17<<" "
      //    <<e18<<" "<<e19<<" "<<e20<<" "<<e21<<" "<<e22<<" "<<t0<<" "<<t1<<" "<<r0<<endl;

      //cout<<val1<<" "<<val2<<endl;    
      
      hh1.Fill((val1-targetX)*34.3);      
      hh2.Fill((val2-targetY)*34.3);
      hh3.Fill((val1-targetX)*34.3,(val2-targetY)*34.3);
      
      TVector2 dxy((val1-targetX)*34.3,(val2-targetY)*34.3); 
      TVector2 du(cos(t0),sin(t0));
      TVector2 dv(-sin(t0),cos(t0));
      TVector2 deltaUV0(dxy*du,dxy*dv);
      hh4.Fill(deltaUV0.X());
      hh5.Fill(deltaUV0.Y());
      hh6.Fill(deltaUV0.X(),deltaUV0.Y());
   }


   TLatex text;
   text.SetTextSize(0.045);
    

   hh1.SetLineWidth(2);
   hh1.GetYaxis()->SetTitleOffset(1.2);
   hh1.GetXaxis()->SetTitleSize(0.045);
   hh1.GetYaxis()->SetTitleSize(0.045);
   hh1.GetXaxis()->SetLabelSize(0.04);
   hh1.GetYaxis()->SetLabelSize(0.04);
 
   hh2.SetLineWidth(2);
   hh2.GetYaxis()->SetTitleOffset(1.2);
   hh2.GetXaxis()->SetTitleSize(0.045);
   hh2.GetYaxis()->SetTitleSize(0.045);
   hh2.GetXaxis()->SetLabelSize(0.04);
   hh2.GetYaxis()->SetLabelSize(0.04);

   TF1 *myfit = new TF1("myfit","gaus(0)+gaus(3)", -40,40);
   myfit->SetParameter(0,1000);
   myfit->SetParameter(1,0);
   myfit->SetParameter(2,8);
   myfit->SetParameter(3,100);
   myfit->SetParameter(4,0);
   myfit->SetParameter(5,30);
 
   TCanvas *ca1 = new TCanvas("ca1","ca1");
   hh1.Draw();
   hh1.SetLineWidth(2);
   hh1.Fit("myfit","","",-40,40);   
   hh1.GetXaxis()->SetTitle("X_{rec}-X_{gen} (mm)");
   hh1.GetYaxis()->SetTitle("Entries / 1 mm");
   text.DrawLatexNDC(0.17,0.82,Form("#sigma_{core} =%4.2g #pm %2.1g mm",myfit->GetParameter(2),myfit->GetParError(2)));
   text.DrawLatexNDC(0.17,0.75,Form("#sigma_{tail} =%4.2g #pm %2.1g mm",myfit->GetParameter(5),myfit->GetParError(5)));
   text.DrawLatexNDC(0.17,0.68,Form("frac = %4.3g ",myfit->GetParameter(3)/myfit->GetParameter(0)));
   ca1->SaveAs("plots/cluX.pdf");
   cout<<endl;
   cout<<"Results sigma core / tail = "<<myfit->GetParameter(2)<<"  /  "<<myfit->GetParameter(5)<<endl;
   cout<<"Results RMS = "<<hh1.GetRMS()<<endl;
   cout<<endl;
   cout<<endl;
   

   hh2.Draw();
   hh2.SetLineWidth(2);
   hh2.Fit("myfit","","",-40,40);
   hh2.GetXaxis()->SetTitle("Y_{rec}-Y_{gen} (mm)");
   hh2.GetYaxis()->SetTitle("Entries / 1 mm");
   text.DrawLatexNDC(0.17,0.82,Form("#sigma_{core} =%4.2g #pm %2.1g mm",myfit->GetParameter(2),myfit->GetParError(2)));
   text.DrawLatexNDC(0.17,0.75,Form("#sigma_{tail} =%4.2g #pm %2.1g mm",myfit->GetParameter(5),myfit->GetParError(5)));
   text.DrawLatexNDC(0.17,0.68,Form("frac = %4.3g ",myfit->GetParameter(3)/myfit->GetParameter(0)));
   ca1->SaveAs("plots/cluY.pdf");
   cout<<endl;
   cout<<"Results sigma core / tail = "<<myfit->GetParameter(2)<<"  /  "<<myfit->GetParameter(5)<<endl;
   cout<<"Results RMS = "<<hh2.GetRMS()<<endl;
   cout<<endl;
   cout<<endl;

   hh4.Draw();
   hh4.SetLineWidth(2);
   hh4.Fit("myfit","","",-40,40);
   hh4.GetXaxis()->SetTitle("U_{rec}-U_{gen} (mm)");
   hh4.GetYaxis()->SetTitle("Entries / 1 mm");
   text.DrawLatexNDC(0.17,0.82,Form("#sigma_{core} =%4.2g #pm %2.1g mm",myfit->GetParameter(2),myfit->GetParError(2)));
   text.DrawLatexNDC(0.17,0.75,Form("#sigma_{tail} =%4.2g #pm %2.1g mm",myfit->GetParameter(5),myfit->GetParError(5)));
   text.DrawLatexNDC(0.17,0.68,Form("frac = %4.3g ",myfit->GetParameter(3)/myfit->GetParameter(0)));
   ca1->SaveAs("plots/cluU.pdf");
   cout<<endl;

   hh5.Draw();
   hh5.Fit("myfit","","",-40,40);
   hh5.SetLineWidth(2);
   hh5.GetXaxis()->SetTitle("V_{rec}-V_{gen} (mm)");
   hh5.GetYaxis()->SetTitle("Entries / 1 mm");
   text.DrawLatexNDC(0.17,0.82,Form("#sigma_{core} =%4.2g #pm %2.1g mm",myfit->GetParameter(2),myfit->GetParError(2)));
   text.DrawLatexNDC(0.17,0.75,Form("#sigma_{tail} =%4.2g #pm %2.1g mm",myfit->GetParameter(5),myfit->GetParError(5)));
   text.DrawLatexNDC(0.17,0.68,Form("frac = %4.3g ",myfit->GetParameter(3)/myfit->GetParameter(0)));
   ca1->SaveAs("plots/cluV.pdf");
   cout<<endl;
   
   
   hh4.Fit("myfit","q","",-40,40);   
   hh5.Fit("myfit","q","",-40,40);   
   
   
   
   
   TFile res(Form("dataMu2e/results_%i.root",mode),"RECREATE");
   res.Add(&hh1);
   res.Add(&hh2);
   res.Add(&hh3);
   res.Add(&hh4);
   res.Add(&hh5);
   res.Add(&hh6);
   res.Write();
   res.Close();   
   
}

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
#include "TVector3.h"


using namespace TMVA;
   
void TMVAReadExtend(TString testDS, int mode=1) 
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
   Float_t r0,z0,t0,t1,c0,c1,c2,c3,c4,c5,FFX,FFY,FFZ;
   Float_t targetX,targetY,targetA;


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
   reader->AddVariable( "e20",  &e20);
   reader->AddVariable( "t0",   &t0);
   reader->AddVariable( "t1",   &t1);
   reader->AddVariable( "r0",   &r0);
   reader->AddVariable( "z0",   &z0);
   reader->AddSpectator( "c0",  &c0);
   reader->AddSpectator( "c1",  &c1);
   reader->AddSpectator( "c2",  &c2);
   reader->AddSpectator( "c3",  &c3); 
   reader->AddSpectator( "c4",  &c4); 
   reader->AddSpectator( "c5",  &c5); 
   reader->AddSpectator( "ffx", &FFX); 
   reader->AddSpectator( "ffy", &FFY); 
   reader->AddSpectator( "ffz", &FFZ); 
   reader->BookMVA("BDTG method", "weightsA/weights/TMVARegression_BDTG.weights.xml" ); 
   


   TFile input(testDS);
   TTree *regTree = (TTree*)input.Get("TreeReg");

   regTree->SetBranchAddress( "e0",  &e0 );
   regTree->SetBranchAddress( "e1",  &e1 );
   regTree->SetBranchAddress( "e2",  &e2 );
   regTree->SetBranchAddress( "e3",  &e3 );
   regTree->SetBranchAddress( "e4",  &e4 );
   regTree->SetBranchAddress( "e5",  &e5 );
   regTree->SetBranchAddress( "e6",  &e6 );
   regTree->SetBranchAddress( "e7",  &e7 );
   regTree->SetBranchAddress( "e8",  &e8 );
   regTree->SetBranchAddress( "e9",  &e9 );
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
   regTree->SetBranchAddress( "t0",  &t0 );
   regTree->SetBranchAddress( "t1",  &t1 );
   regTree->SetBranchAddress( "r0",  &r0 );
   regTree->SetBranchAddress( "z0",  &z0 );
   regTree->SetBranchAddress( "c0",  &c0 );
   regTree->SetBranchAddress( "c1",  &c1 );
   regTree->SetBranchAddress( "c2",  &c2 );
   regTree->SetBranchAddress( "c3",  &c3 );
   regTree->SetBranchAddress( "c4",  &c4 );
   regTree->SetBranchAddress( "c5",  &c5 );
   regTree->SetBranchAddress( "ffx", &FFX );
   regTree->SetBranchAddress( "ffy", &FFY );
   regTree->SetBranchAddress( "ffz", &FFZ );
   regTree->SetBranchAddress( "targetX", &targetX );
   regTree->SetBranchAddress( "targetY", &targetY );
   regTree->SetBranchAddress( "targetA", &targetA );


   TH1F hh1("hh1",  "",100,-50,50);
   TH1F hh2("hh2",  "",100,-50,50);
   TH1F hh3("hh3",  "",100,-50,50);


   for (Long64_t ievt=0; ievt<regTree->GetEntries();ievt++)
   {
      regTree->GetEntry(ievt);
      
      Float_t val1 = 100*((reader->EvaluateRegression("BDTG method"))[0]);
      
      
      TVector3 trkDir(sin(t1)*cos(t0),sin(t1)*sin(t0),cos(t1));
      TVector3 pos0(FFX,FFY,FFZ);
      
      TVector3 pos = pos0+val1*trkDir;
      double dx = c2-pos.X();
      double dy = c3-pos.Y();      
      
      hh1.Fill(dx);      
      hh2.Fill(dy);
      hh3.Fill(val1-targetA);
      
      //need to add stuff to calculate X ad Y position
      
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
   myfit->SetParameter(5,20);
 
   TCanvas *ca1 = new TCanvas("ca1","ca1");
   hh1.Draw();
   hh1.SetLineWidth(2);
   hh1.Fit("myfit","","",-40,40);   
   hh1.GetXaxis()->SetTitle("X_{rec}-X_{gen} (mm)");
   hh1.GetYaxis()->SetTitle("Entries / 2 mm");
   text.DrawLatexNDC(0.17,0.82,Form("#sigma_{core} =%4.2g #pm %2.1g mm",myfit->GetParameter(2),myfit->GetParError(2)));
   text.DrawLatexNDC(0.17,0.75,Form("#sigma_{tail} =%4.2g #pm %2.1g mm",myfit->GetParameter(5),myfit->GetParError(5)));
   ca1->SaveAs("cluAX.pdf");
   cout<<endl;
   cout<<"Results sigma core / tail = "<<myfit->GetParameter(2)<<"  /  "<<myfit->GetParameter(5)<<endl;
   cout<<"Results RMS = "<<hh1.GetRMS()<<endl;
   cout<<endl;
   cout<<endl;
   

   hh2.Draw();
   hh2.SetLineWidth(2);
   hh2.Fit("myfit","","",-40,40);
   hh2.GetXaxis()->SetTitle("Y_{rec}-Y_{gen} (mm)");
   hh2.GetYaxis()->SetTitle("Entries / 2 mm");
   text.DrawLatexNDC(0.17,0.82,Form("#sigma_{core} =%4.2g #pm %2.1g mm",myfit->GetParameter(2),myfit->GetParError(2)));
   text.DrawLatexNDC(0.17,0.75,Form("#sigma_{tail} =%4.2g #pm %2.1g mm",myfit->GetParameter(5),myfit->GetParError(5)));
   ca1->SaveAs("cluAY.pdf");
   cout<<endl;
   cout<<"Results sigma core / tail = "<<myfit->GetParameter(2)<<"  /  "<<myfit->GetParameter(5)<<endl;
   cout<<"Results RMS = "<<hh2.GetRMS()<<endl;
   cout<<endl;
   cout<<endl;
    
   TFile res(Form("dataMu2e/results_%i.root",mode),"RECREATE");
   res.Add(&hh1);
   res.Add(&hh2);
   res.Add(&hh3);
   res.Write();
   res.Close();   
   
}

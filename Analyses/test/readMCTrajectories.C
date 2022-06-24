//
// Look at some events made from ReadMCTrajectories.
//
//
// Original author Rob Kutschke
//

{

  gROOT->Reset();

  gROOT->SetStyle("Plain");

  gStyle->SetOptStat("emruoi");

  TFile* file = new TFile( "readMCTrajectories.root");

  TNtuple* nt;    file->GetObject("readMCTraj/ntup",nt);
  TH1F* nTraj;    file->GetObject("readMCTraj/nTraj",nTraj);
  TH1F* nPoints1; file->GetObject("readMCTraj/nPoints1",nPoints1);
  TH1F* nPoints2; file->GetObject("readMCTraj/nPoints2",nPoints2);

  TCanvas *canvas = new TCanvas("c", "MC Trajectoreis", 900, 900 );

  canvas->Divide(2,3);

  canvas->cd(1); nTraj->Draw();
  canvas->cd(2); nPoints1->Draw();
  canvas->cd(3); nPoints2->Draw();
  gPad->SetLogy();

  //gStyle->SetOptStat("emruoi");
  canvas->cd(4); nt->Draw( "x:y:z","evt==1","");
  canvas->cd(5); nt->Draw( "x:y:z","evt==2","");
  canvas->cd(6); nt->Draw( "x:y:z","evt==3","");


}

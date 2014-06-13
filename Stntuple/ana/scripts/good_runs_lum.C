//-----------------------------------------------------------------------------
// this macro plots STNMAKER (dev_224) memory consumption (in SGI2 pages) as a 
// function of the number of processed events
// the plots show mmory leak in TauFinder module, purify points to the 
// following place:
//
//-----------------------------------------------------------------------------
good_runs_lum() {
   c1 = new TCanvas("good_runs_lum","Luminosity in the good runs",0,0,500,700);
   //   c1->SetFillColor(42);

   c1->Divide(1,3);
   c1->cd(1);

   c1->SetGrid();

      // create a 2-d histogram to define the range

   TH2F *hr = new TH2F("grl","Luminosity in the good runs, pb^-1",
		       2,0.,20,2,0,5);
   hr->SetXTitle("cut on run lumi(nb^-1))");
   hr->SetYTitle("Integrated lumi");

   hr->GetYaxis()->SetLabelSize  (0.025);
   hr->GetYaxis()->SetTitleOffset(1.5);
   hr->GetYaxis()->SetTitleSize  (0.025);

   hr->Draw();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(250);

      // create first graph
   Int_t n1 = 13;
   Float_t x1[]  = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 15., 20.};
   Float_t y1[]  = {4.160, 4.107, 4.003, 3.854, 3.736, 3.562, 
                    3.497, 3.328, 3.172, 3.046, 2.886, 2.014, 1.418};
   gr1 = new TGraph(n1,x1,y1);
   gr1->SetMarkerColor(kBlue);
   gr1->SetMarkerStyle(21);
   gr1->Draw("LP");

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
   c1->cd(2);

   TH2F *h1 = new TH2F("grl","Fraction of Luminosity in the good runs, pb^-1",
		       2,0.,20,2,0,1.1);
   h1->SetXTitle("cut on run lumi(nb^-1))");
   h1->SetYTitle("Fraction of luminosity in the good runs");

   h1->Draw();

   float y10[13];

   for (int i=0; i<13; i++) y10[i] = y1[i]/y1[0];
   gr1 = new TGraph(n1,x1,y10);
   gr1->SetMarkerColor(kBlue);
   gr1->SetMarkerStyle(21);
   gr1->Draw("LP");
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
   c1->cd(3);
      // create second graph
   TH2F *h2 = new TH2F("grl","number of runs to process",
		       2,0.,20,2,0,600);

   h2->Draw();

   Int_t n2 = 13;
   Float_t y2[]  = {575, 443, 373, 313, 279, 240, 228, 202, 181, 166, 149, 79, 44};
   gr2 = new TGraph(n2,x1,y2);
    gr2->SetMarkerColor(kRed);
    gr2->SetMarkerStyle(20);
    gr2->Draw("LP");

//     TText* text1 = new TText(2300,7200,"TauFinder OFF");
//     text1->Draw();

//     TText* text2 = new TText(1800,9700,"TauFinder ON");
//     text2->Draw();
}

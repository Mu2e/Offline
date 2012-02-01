void mu2e(double ndio=100000, double nconv=100000) {
  double nstopped(7.52e17);
  double capfrac(0.609); 
  double decayfrac = 1.0 - capfrac;
  double ndecay = nstopped*decayfrac;
  double conprob(1e-15);

  unsigned nbins(25);
  double mmin(101);
  double mmax(106);
  double mevperbin = mmax-mmin/nbins;

  TChain* dio = new TChain("ReadKalFits/trkdiag");
  char diofile[100];
  for(unsigned idio=0;idio<9;idio++){
    snprintf(diofile,100,"/data/7432_%i/flatDIO.root",idio);
    dio->Add(diofile);
  }

  TChain* con = new TChain("ReadKalFits/trkdiag");
  char confile[100];
  for(unsigned icon=0;icon<9;icon++){
    snprintf(confile,100,"/data/7436_%i/mixConversion.root",icon);
    con->Add(confile);
  }

  TH1F* timeshift = new TH1F("timeshift","T0 - muon conversion time;nsec",100,0,100);
  con->Project("timeshift","t0-mct0","fitstatus>0");

  TProfile* truedio = new TProfile("truedio","Weighted MC true DIO spectrum;MeV;d#Gamma/#Gamma (MeV^{-1})",100,100,105.5,0,1);
  dio->Project("truedio","diowt:mcmom");
  truedio->SetMaximum(1e-12);
  truedio->SetMinimum(1e-24);

  TCanvas* valid = new TCanvas("valid","validation",1200,800);
  valid->Divide(1,2);
  valid->cd(1);
  gPad->SetLogy();
  truedio->Draw();
  valid->cd(2);
  timeshift->Draw();

  TH1F* diospec[4];
  TH1F* conspec[4];
// basic cuts
  TCut reco("fitstatus>0");
  TCut pitch("td>0.577 && td <1.0");
  TCut livegate("t0>800");
// cuts for different tightness of selection
  TCut ncuts[4], t0cuts[4], momcuts[4], fitcuts[4];
  ncuts[0] = "nactive>=20";
  ncuts[1] = "nactive>=20";
  ncuts[2] = "nactive>=25";
  ncuts[3] = "nactive>=30";
  t0cuts[0] = "";
  t0cuts[1] = "t0err<1.5";
  t0cuts[2] = "t0err<1.0";
  t0cuts[3] = "t0err<0.9";
  momcuts[0] = "";
  momcuts[1] = "fitmomerr<0.2";
  momcuts[2] = "fitmomerr<0.18";
  momcuts[3] = "fitmomerr<0.15";
  fitcuts[0] = "";
  fitcuts[1] = "fitcon>1e-4";
  fitcuts[2] = "fitcon>1e-3";
  fitcuts[3] = "fitcon>1e-2";


  for(unsigned ires=0;ires<4;ires++){
    char dioname[50];
    snprintf(dioname,50,"diospec%i",ires);
    char conname[50];
    snprintf(conname,50,"conspec%i",ires);
 
    diospec[ires] = new TH1F(dioname,"Reconstructed momentum;MeV;Events per experiment",nbins,mmin,mmax);
    diospec[ires]->SetStats(0);
    diospec[ires]->SetLineColor(kBlue);

    conspec[ires] = new TH1F(conname,"Reconstructed momentum;MeV;Events per experiment",nbins,mmin,mmax);
    conspec[ires]->SetStats(0);
    conspec[ires]->SetLineColor(kRed);

    TCut quality = ncuts[ires] && t0cuts[ires] && momcuts[ires] && fitcuts[ires];

    dio->Project(dioname,"fitmom","diowt"*(reco+pitch+livegate+quality));
    diospec[ires]->Scale(ndecay*mevperbin/ndio);

    double conscale = ndecay*conprob/nconv;
    con->Project(conname,"fitmom",reco+pitch+livegate+quality);
    conspec[ires]->Scale(conscale);

  }

  TLegend* leg = new TLegend(0.7,0.8,0.9,0.9);
  leg->AddEntry(diospec[0],"DIO","L");
  leg->AddEntry(conspec[0],"Conversion","L");

  TPaveText* info = new TPaveText(0.4,0.8,0.7,0.9,"NDC");
  char text[80];
  snprintf(text,80,"%g stopped muons",nstopped);
  TString snstop(text);
  info->AddText(snstop);
  snprintf(text,80,"%g Conversion Rate",conprob);
  TString sconprob(text);
  info->AddText(sconprob);
  info->SetBorderSize(0);
 
// log scale
  TCanvas* logcan = new TCanvas("mu2elog","mu2e log scale",1200,800);
  logcan->Clear();
  logcan->Divide(2,2);

  for(unsigned ires=0;ires<4;ires++){
    logcan->cd(ires+1);
    gPad->SetLogy();
    diospec[ires]->Draw();
    conspec[ires]->Draw("same");

    TPaveText* logtext = new TPaveText(0.1,0.5,0.4,0.8,"NDC");  
    char line[40];
    sprintf(line,"%s",ncuts[ires].GetTitle());
    logtext->AddText(line);
    sprintf(line,"%s",t0cuts[ires].GetTitle());
    logtext->AddText(line);
    sprintf(line,"%s",momcuts[ires].GetTitle());
    logtext->AddText(line);
    sprintf(line,"%s",fitcuts[ires].GetTitle());
    logtext->AddText(line);
    logtext->Draw();

  }
  logcan->cd(2);
  leg->Draw();
  info->Draw();
 

// linear scale
  TCanvas* lincan = new TCanvas("mu2elin","mu2e linear scale",1200,800);
  lincan->Clear();
  lincan->Divide(2,2);
   for(unsigned ires=0;ires<4;ires++){
    lincan->cd(ires+1);
    TH1F* diocopy = diospec[ires]->DrawCopy();
    diocopy->SetMaximum(6);
    conspec[ires]->Draw("same");

    TPaveText* lintext = new TPaveText(0.1,0.5,0.4,0.8,"NDC");  
    char line[40];
    sprintf(line,"%s",ncuts[ires].GetTitle());
    lintext->AddText(line);
    sprintf(line,"%s",t0cuts[ires].GetTitle());
    lintext->AddText(line);
    sprintf(line,"%s",momcuts[ires].GetTitle());
    lintext->AddText(line);
    sprintf(line,"%s",fitcuts[ires].GetTitle());
    lintext->AddText(line);
    lintext->Draw();

  }
  lincan->cd(2);
  leg->Draw();
  info->Draw();
 
}

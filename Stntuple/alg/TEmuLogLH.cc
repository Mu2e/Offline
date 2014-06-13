//-----------------------------------------------------------------------------

#include "Stntuple/alg/TEmuLogLH.hh"
#include "TMath.h"
// #include "murat/plot/smooth_new.hh"

ClassImp(TEmuLogLH)
//-----------------------------------------------------------------------------
// default constructor
//-----------------------------------------------------------------------------
TEmuLogLH::TEmuLogLH() {
  fEleDtHist = NULL;
  fEleEpHist = NULL;
  fEleXsHist = NULL;
  fEleDeHist = NULL;

  fMuoDtHist = NULL;
  fMuoEpHist = NULL;
  fMuoXsHist = NULL;
  fMuoDeHist = NULL;
}

//-----------------------------------------------------------------------------
// real costructor
//-----------------------------------------------------------------------------
TEmuLogLH::TEmuLogLH(const char* FnEle,  const char* FnMuo) {
  //  TH1  *h1;

  printf(" >>> ERROR : TEmuLogLH::TEmuLogLH(const char* FnEle,  const char* FnMuo) called by mistake\n");
  
  // h1 = get_mu2e_histogram(FnEle,"TCalm002/Hist/trk_1/dt");
  // fEleDtHist = (TH1*) h1->Clone("ele_dt");

  // // rebin if necessary
  // fEleDtHist->Rebin(4);
  // fEleDtFunc = new smooth_new(fEleDtHist);

  // h1 = get_mu2e_histogram(FnEle,"TCalm002/Hist/trk_1/ep");
  // fEleEpHist = (TH1*) h1->Clone("ele_ep");

  // // rebin if necessary
  // fEleEpHist->Rebin(1);
  // fEleEpFunc = new smooth_new(fEleEpHist);

//-----------------------------------------------------------------------------
// muon part
//-----------------------------------------------------------------------------
  // h1 = get_mu2e_histogram(FnMuo,"TCalm002/Hist/trk_1/dt");
  // fMuoDtHist = (TH1*) h1->Clone("muo_dt");

  // // rebin if necessary
  // fMuoDtHist->Rebin(4);
  // fMuoDtFunc = new smooth_new(fMuoDtHist);

  // h1 = get_mu2e_histogram(FnMuo,"TCalm002/Hist/trk_1/ep");
  // fMuoEpHist = (TH1*) h1->Clone("muo_ep");

  // // rebin if necessary
  // fMuoEpHist->Rebin(1);
  // fMuoEpFunc = new smooth_new(fMuoEpHist);
}


//-----------------------------------------------------------------------------
TEmuLogLH::~TEmuLogLH() {
  if (fEleDtHist) delete fEleDtHist;
  if (fEleEpHist) delete fEleEpHist;
  if (fEleXsHist) delete fEleXsHist;
  if (fEleDeHist) delete fEleDeHist;

  if (fMuoDtHist) delete fMuoDtHist;
  if (fMuoEpHist) delete fMuoEpHist;
  if (fMuoXsHist) delete fMuoXsHist;
  if (fMuoDeHist) delete fMuoDeHist;
}

//-----------------------------------------------------------------------------
// for the time being, assume that there is no under/overflows
//-----------------------------------------------------------------------------
void  TEmuLogLH::SetEleDtHist(TH1* Hist) { 	
  if (fEleDtHist) delete fEleDtHist;
  fEleDtHist         = (TH1*) Hist->Clone(Form("%s_ele_dt_hist_clone",Hist->GetName())); 
  fEleDtHist->Scale(1./Hist->Integral());
}

//-----------------------------------------------------------------------------
void  TEmuLogLH::SetEleEpHist(TH1* Hist) { 
  if (fEleEpHist) delete fEleEpHist;
  fEleEpHist         = (TH1*) Hist->Clone(Form("%s_EleEpHist_clone",Hist->GetName())); 
  fEleEpHist->Scale(1./Hist->Integral());
}

//-----------------------------------------------------------------------------
void  TEmuLogLH::SetEleXsHist(TH1* Hist) { 
  if (fEleXsHist) delete fEleXsHist;
  fEleXsHist         = (TH1*) Hist->Clone(Form("%s_EleXsHist_clone",Hist->GetName())); 
  fEleXsHist->Scale(1./Hist->Integral());
}

// //-----------------------------------------------------------------------------
// void  TEmuLogLH::SetEleDeHist(TH1* Hist) { 
//   if (fEleDeHist) delete fEleDeHist;
//   fEleDeHist         = (TH1*) Hist->Clone(Form("%s_EleDeHist_clone",Hist->GetName())); 
//   fEleDeHist->Scale(1./Hist->Integral());
// }

//-----------------------------------------------------------------------------
void  TEmuLogLH::SetMuoDtHist(TH1* Hist) { 
  if (fMuoDtHist) delete fMuoDtHist;
  fMuoDtHist         = (TH1*) Hist->Clone(Form("%s_MuoDtHist_clone",Hist->GetName())); 
  fMuoDtHist->Scale(1./Hist->Integral());
}

//-----------------------------------------------------------------------------
void  TEmuLogLH::SetMuoEpHist(TH1* Hist) { 
  if (fMuoEpHist) delete fMuoEpHist;
  fMuoEpHist         = (TH1*) Hist->Clone(Form("%s_MuoEpHist_clone",Hist->GetName())); 
  fMuoEpHist->Scale(1./Hist->Integral());
}

//-----------------------------------------------------------------------------
void  TEmuLogLH::SetMuoXsHist(TH1* Hist) { 
  if (fMuoXsHist) delete fMuoXsHist;
  fMuoXsHist         = (TH1*) Hist->Clone(Form("%s_MuoXsHist_clone",Hist->GetName())); 
  fMuoXsHist->Scale(1./Hist->Integral());
}


// //-----------------------------------------------------------------------------
// void  TEmuLogLH::SetMuoDeHist(TH1* Hist) { 
//   if (fMuoDeHist) delete fMuoDeHist;
//   fMuoDeHist         = (TH1*) Hist->Clone(Form("%s_MuoDeHist_clone",Hist->GetName())); 
//   fMuoDeHist->Scale(1./Hist->Integral());
// }


//-----------------------------------------------------------------------------
double TEmuLogLH::LogLHDt(double Dt, int PdgCode) {

  double log_lh, p1;
  int    i1;

  TH1    *h1(0);   

  if      (abs(PdgCode) == 11) {
//-----------------------------------------------------------------------------
// electron branch , assume constant bin
//-----------------------------------------------------------------------------
    h1 = fEleDtHist;
  }
  else if (abs(PdgCode) == 13) {
//-----------------------------------------------------------------------------
// muon branch 
//-----------------------------------------------------------------------------
    h1 = fMuoDtHist;
  }

  i1 = (Dt-h1->GetXaxis()->GetXmin())/h1->GetBinWidth(1)+1;
  p1 = h1->GetBinContent(i1);
  if (p1 < 1.e-15) p1 = 1.e-15;

  log_lh = TMath::Log(p1);
  
  return log_lh;
}


//-----------------------------------------------------------------------------
double TEmuLogLH::LogLHEp(double Ep, int PdgCode) {

  double log_lh, ep, p1;
  int    i1;

  TH1    *h1(0);   

  if      (abs(PdgCode) == 11) {
//-----------------------------------------------------------------------------
// electron branch , assume constant bin
//-----------------------------------------------------------------------------
    h1 = fEleEpHist;
  }
  else if (abs(PdgCode) == 13) {
//-----------------------------------------------------------------------------
// muon branch 
//-----------------------------------------------------------------------------
    h1 = fMuoEpHist;
  }

  i1 = (Ep-h1->GetXaxis()->GetXmin())/h1->GetBinWidth(1)+1;
  p1 = h1->GetBinContent(i1);
  if (p1 < 1.e-15) p1 = 1.e-15;

  log_lh = TMath::Log(p1);
  
  return log_lh;
}


//-----------------------------------------------------------------------------
double TEmuLogLH::LogLHXs(double Xs, int PdgCode) {

  double log_lh, p1;
  int    i1;

  TH1    *h1(0);   

  if      (abs(PdgCode) == 11) {
//-----------------------------------------------------------------------------
// electron branch , assume constant bin
//-----------------------------------------------------------------------------
    h1 = fEleXsHist;
  }
  else if (abs(PdgCode) == 13) {
//-----------------------------------------------------------------------------
// muon branch 
//-----------------------------------------------------------------------------
    h1 = fMuoXsHist;
  }

  i1 = (Xs-h1->GetXaxis()->GetXmin())/h1->GetBinWidth(1)+1;
  p1 = h1->GetBinContent(i1);
  if (p1 < 1.e-15) p1 = 1.e-15;

  log_lh = TMath::Log(p1);
  
  return log_lh;
}


// //-----------------------------------------------------------------------------
// double TEmuLogLH::LogLHDe(double De, int PdgCode) {

//   double llh, ep, p1;
//   int    i1;

//   TH1    *h1(0);   

//   if      (abs(PdgCode) == 11) {
// //-----------------------------------------------------------------------------
// // electron branch , assume constant bin
// //-----------------------------------------------------------------------------
//     h1 = fEleDeHist;
//   }
//   else if (abs(PdgCode) == 13) {
// //-----------------------------------------------------------------------------
// // muon branch 
// //-----------------------------------------------------------------------------
//     h1 = fMuoDeHist;
//   }

//   i1 = (De-h1->GetXaxis()->GetXmin())/h1->GetBinWidth(1)+1;
//   p1 = h1->GetBinContent(i1);
//   if (p1 < 1.e-15) p1 = 1.e-15;

//   llh = TMath::Log(p1);
  
//   return llh;
// }


//-----------------------------------------------------------------------------
double TEmuLogLH::LogLHCal(CalData_t* Data, int PdgCode) {

  double llh, dt, ep;
  //  int    i1, i2;

  //  TH1    *h1(0), *h2(0);   

  dt = Data->fDt;
  ep = Data->fEp;

  llh = LogLHDt(dt,PdgCode)+LogLHEp(ep,PdgCode);

  return llh;
}

// //-----------------------------------------------------------------------------
// double TEmuLogLH::LogLHTrk(TrkData_t* Data, int PdgCode) {

//   double llh, dt, ep, p1, p2;
//   int    i1, i2;

//   TH1    *h1(0), *h2(0);   

//   xs = Data->fXs;
//   de = Data->fDe;

//   llh = LogLHXs(xs,PdgCode)+LogLHDe(de,PdgCode);

//   return llh;
// }

//-----------------------------------------------------------------------------
// 11 and 13 - electron and muon PDG codes
//-----------------------------------------------------------------------------
double TEmuLogLH::LogLHRDt(double Dt) {
  double llhr;
  llhr = LogLHDt(Dt,11)-LogLHDt(Dt,13);
  return llhr;
}

//-----------------------------------------------------------------------------
double TEmuLogLH::LogLHREp(double Ep) {
  double llhr;
  llhr = LogLHEp(Ep,11)-LogLHEp(Ep,13);
  return llhr;
}

//-----------------------------------------------------------------------------
double TEmuLogLH::LogLHRXs(double Xs) {
  double llhr;
  llhr = LogLHXs(Xs,11)-LogLHXs(Xs,13);
  return llhr;
}

// //-----------------------------------------------------------------------------
// double TEmuLogLH::LogLHRDe(double De) {
//   double llhr;
//   llhr = LogLHDe(De,11)-LogLHDe(De,13);
//   return llhr;
// }


//-----------------------------------------------------------------------------
// 11 and 13 - electron and muon PDG codes
//-----------------------------------------------------------------------------
double TEmuLogLH::LogLHRCal(CalData_t* Data) {
  double llhr = LogLHRDt(Data->fDt)+LogLHREp(Data->fEp);
  return llhr;
}

// double TEmuLogLH::LogLHRTrk(TrkData_t* Data) {
//   double llhr = LogLHRXs(Data->fXs)+LogLHRDe(Data->fDe);
//   return llhr;
// }

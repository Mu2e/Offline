Double_t crystalball (Double_t *x, Double_t *par) {
  // par[0] : norm
  // par[1] : x0
  // par[2] : sigma
  // par[3] : n
  // par[4] : alpha
  // par[5] : fraction of 2nd Gaussian
  // par[6] : tail gaussian sigma

  double DeltaX = (x[0]- par[1]);
  double absSigma = fabs(par[2]);
  double fval=0.0;
  if ( DeltaX/absSigma > -1.*par[4]) {
    double g = par[0]*TMath::Gaus(x[0], par[1], par[2]);
//    double g2 = par[5]*par[0]*TMath::Gaus(x[0], par[1], par[6]);
//    return g1+g2;
    double e = 0.0;//par[0]*par[5]*exp(-DeltaX/par[6])/par[6];
    //return g+e;
    fval = g+e;
  }
  else {
	double absAlpha = fabs(par[4]);
    double A = pow(par[3]/absAlpha, par[3])*exp(-0.5*par[4]*par[4]);
    double B = par[3]/absAlpha - absAlpha;
    //return par[0]*A*pow(B-DeltaX/absSigma, -1.*par[3]);
    fval = par[0]*A*pow(B-DeltaX/absSigma, -1.*par[3]);
  }
  return fval;
}

Double_t doubleCrystalball (Double_t *x, Double_t *par) {
	  // par[0] : norm
	  // par[1] : x0
	  // par[2] : sigma
	  // par[3] : n
	  // par[4] : alpha
	  // par[5] : fraction core
	  // par[6] : tail gaussian sigma
	  // par[7] : tail n
	  // par[8] : tail alpha
	  
	  const double invSqrt2pi=0.398942280401432702863;

	  double DeltaX = (x[0]- par[1]);
	  double absSigma = fabs(par[2]);
	  double fval = par[0]*invSqrt2pi/par[2];
	  if ( DeltaX/absSigma > -1.*par[4]) {
	    fval *= TMath::Gaus(x[0], par[1], par[2]);
	  }
	  else {
	    double absAlpha = fabs(par[4]);
	    double A = pow(par[3]/absAlpha, par[3])*exp(-0.5*par[4]*par[4]);
	    double B = par[3]/absAlpha - absAlpha;
	    fval *= A*pow(B-DeltaX/absSigma, -1.*par[3]);
	  }

	  double tailFval=par[0]*invSqrt2pi/par[6];
	  double tailAbsSigma = fabs(par[6]);
	  DeltaX*=-1.0;
	  if ( DeltaX/tailAbsSigma > -1.*par[8]) {
	    tailFval *= TMath::Gaus(x[0], par[1], par[6]);
	  }
	  else {
	    double tailAbsAlpha = fabs(par[8]);
	    double tailA = pow(par[7]/tailAbsAlpha, par[7])*exp(-0.5*par[8]*par[8]);
	    double tailB = par[7]/tailAbsAlpha - tailAbsAlpha;
	    tailFval *= tailA*pow(tailB-DeltaX/tailAbsSigma, -1.*par[7]);
	  }

	  return fval*par[5] + (1.0-par[5])*tailFval;
}


void MomRes_1(TFile *d, int nGenEv=-1, TString fitOpt="L", TString fold="kalmanFit"){

  if (nGenEv<=0) nGenEv=((TH1F*)d->Get("g4run/totalMultiplicity"))->GetEntries();
  TString treeName = fold+"/trfit";
  TTree* t = (TTree *) d->Get(treeName.Data());
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat("oumr");

  TCut enterSim("0.299792458/startpar.paramVec.m[2]*TMath::Sqrt(1.0+pow(startpar.paramVec.m[4],2))>100&&fitinfo.nhits>50");
  TCut fitcut("TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))>1.0e-4");
  TCut reco("fitinfo.fit>0");

  float norm = t->Draw("startpar.paramVec.m[2]",enterSim);

  TCut final = (reco+fitcut+enterSim);

  TH1F *hMomRes = new TH1F("hMomRes","momentum resolution;Mev/c;Entries / (25 keV/c)",200,-2.5,2.5);

  float integral = t->Draw("0.299792458/recopar.paramVec.m[2]*TMath::Sqrt(1.0+pow(recopar.paramVec.m[4],2))-0.299792458/startpar.paramVec.m[2]*TMath::Sqrt(1.0+pow(startpar.paramVec.m[4],2))>>hMomRes",final);

  /*
  TF1 *dgaus = new TF1("dgaus","gaus(0)+gaus(3)");
  dgaus->SetParName(0,"Norm");
  dgaus->SetParName(1,"#mu_{core}");
  dgaus->SetParName(2,"#sigma_{core}");
  dgaus->SetParName(3,"Norm_{tail}");
  dgaus->SetParName(4,"#mu_{tail}");
  dgaus->SetParName(5,"#sigma_{tail}");
  dgaus->SetParameters(integral,hMomRes->GetMean(1)*0.5,hMomRes->GetRMS(1)*0.5,integral/10,hMomRes->GetMean(1)*2.0,hMomRes->GetRMS(1)*1.5);
  dgaus->SetParLimits(0,1,integral*2);
  dgaus->SetParLimits(1,-0.2,0.4);
  dgaus->SetParLimits(2,0,hMomRes->GetRMS(1));
  dgaus->SetParLimits(3,1,integral*0.5);
  dgaus->SetParLimits(4,-2.5,2.5);
  dgaus->SetParLimits(5,0,hMomRes->GetRMS(1)*3.0);

  //dgaus->SetRange(-0.75,2.5);
  dgaus->SetNpx(1000);
  hMomRes->Fit("dgaus","","",-0.75,2.5);
  */

  /*
  TF1* cball = new TF1("cball",crystalball,-2.0,1.5,7);
  cball->SetParName(0,"Norm");
  cball->SetParName(1,"x0");
  cball->SetParName(2,"sigma");
  cball->SetParName(3,"n");
  cball->SetParName(4,"alpha");
  cball->SetParName(5,"tailfrac");
  cball->SetParName(6,"taillambda");
  cball->SetParameters(integral,hMomRes->GetMean(1)*0.5,hMomRes->GetRMS(1)*0.5,1.0,1.0,0.05,0.5);
  cball->SetParLimits(0,1,integral*2);
  cball->SetParLimits(1,-0.2,0.4);
  cball->SetParLimits(2,0,hMomRes->GetRMS(1));
  cball->SetParLimits(3,0,100.);
  cball->SetParLimits(4,-10,10);
  cball->SetParLimits(5,0,1.0);
  cball->SetParLimits(6,0,10.0);
  cball->FixParameter(5,0.0);
  cball->FixParameter(6,1.0);

  cball->SetNpx(1000);
  hMomRes->Fit("cball","","",-2.5,2.5);
  */

  TF1* dcball = new TF1("dcball",doubleCrystalball,-2.5,2.5,9);
  dcball->SetParName(0,"Norm");
  dcball->SetParName(1,"x0");
  dcball->SetParName(2,"sigma");
  dcball->SetParName(3,"n");
  dcball->SetParName(4,"alpha");
  dcball->SetParName(5,"frac_{core}");
  dcball->SetParName(6,"sigma_{tail}");
  dcball->SetParName(7,"n_{tail}");
  dcball->SetParName(8,"alpha_{tail}");

  dcball->SetParameters(integral,hMomRes->GetMean(1)*0.5,hMomRes->GetRMS(1)*0.25,1.0,1.0
		  ,0.8,hMomRes->GetRMS(1),1.0,1.0);
  dcball->SetParLimits(0,1,integral*2);
  dcball->SetParLimits(1,-0.2,0.4);
  dcball->SetParLimits(2,0.04,hMomRes->GetRMS(1)*0.7);
  dcball->SetParLimits(3,0,20.);
  dcball->SetParLimits(4,0,5);
  dcball->SetParLimits(5,0,1.0);
  dcball->SetParLimits(6,0.1,hMomRes->GetRMS(1)*2.0);
  dcball->SetParLimits(7,0,20.);
  dcball->SetParLimits(8,0,5);

  dcball->SetNpx(1000);
  hMomRes->Fit("dcball","","",-2.5,2.5);
  hMomRes->Fit("dcball","","",-2.5,2.5);
  hMomRes->Fit("dcball",fitOpt.Data(),"",-2.5,2.5);

  double covMat[9][9];
  gMinuit->mnemat(&covMat[0][0],9);
  double corrF = covMat[2][6]/(sqrt(covMat[2][2]*covMat[6][6]));
  double wSigmaCore = dcball->GetParameter(2)*dcball->GetParameter(5);
  double wSigmaTail = dcball->GetParameter(6)*(1.0-dcball->GetParameter(5));
  if (fitOpt.Contains("L")) {corrF = sqrt(2.0)*fabs(corrF);}
  double meanSigma = sqrt(wSigmaCore*wSigmaCore + wSigmaTail*wSigmaTail + 2.0*corrF*wSigmaCore*wSigmaTail);
  cout<<"sigmas correlation "<<corrF<<endl;
  //TVirtualFitter *fitter0 = TVirtualFitter::Fitter(dcball);
  //fitter0->g

  int intFirstBin = 1 + (int)(((dcball->GetParameter(1)-3.0*meanSigma) - hMomRes->GetXaxis()->GetXmin() )/hMomRes->GetBinWidth(1));
  int intLastBin  = 1 + (int)(((dcball->GetParameter(1)+3.0*meanSigma) - hMomRes->GetXaxis()->GetXmin() )/hMomRes->GetBinWidth(1));

  float nEvtIn3sigma = hMomRes->Integral(intFirstBin,intLastBin);
  cout<<"meanSigma "<<meanSigma<<" intFirstBin "<<intFirstBin<<" intLastBin "<<intLastBin<<" nEvtIn3sigma "<<nEvtIn3sigma<<endl;

  TPaveText* ttext = new TPaveText(0.1,0.70,0.4,0.9,"NDC");
  ttext->SetTextSizePixels(18);
  ttext->AddText("Truth Cuts");
  ttext->AddText("momin>100");
  ttext->AddText("nHits>50");
  ttext->AddText(Form("accept = %.2f %%",(norm/((float)nGenEv))*100.0));
  ttext->Draw();
  TPaveText* rtext = new TPaveText(0.1,0.45,0.4,0.70,"NDC");
  rtext->SetTextSizePixels(18);
  rtext->AddText("Reco Cuts");
  rtext->AddText("fitStatus>0");
  rtext->AddText("fitProb>1.0e-4");
  rtext->AddText(Form("eff = %.2f %%",(integral/norm)*100.0));
  rtext->Draw();
  TPaveText* vtext = new TPaveText(0.1,0.30,0.4,0.45,"NDC");
  vtext->SetTextSizePixels(18);
  vtext->AddText(Form("#bar{#sigma} = %.3f",meanSigma));
  vtext->AddText(Form("eff (in x0#pm 3#bar{#sigma}) = %.2f %%",(nEvtIn3sigma/norm)*100.0));
  vtext->Draw();

}

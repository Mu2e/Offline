#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TCut.h"
#include "TPaveText.h"
#include "TLine.h"
#include <iostream>


// the following approximation is from Czarnecki etal, 'Muon decay in orbit:spectrum of high-energy electrons',
// for E>85 MeV
Double_t DIOCZ(Double_t *x, Double_t *par) {
  double ee = x[0];
  double norm = par[0];
  double mal(25133);
//    double mmu(105.654);
  double emu(105.194);
//    double emue(104.973);
//    double me(0.511);
  double a5(8.6434e-17);
  double a6(1.16874e-17);
  double a7(-1.87828e-19);
  double a8(9.16327e-20);
  double delta = emu - ee - ee*ee/(2*mal);
  return norm*(a5*pow(delta,5) + a6*pow(delta,6) + a7*pow(delta,7) + a8*pow(delta,8));
}

void mu2e_IT(TTree* dio, TTree* con, double diogenrange, double ndio, double ncon, TString addCut="",bool weightdio=true, int useExtapMom=0, double dioeff=1.0, double coneff=1.0, Long64_t NEntries=-1, Long64_t skipentries=0) {

  Long64_t nentries=1000000000;
  if (NEntries>0 && NEntries<1000000000) {nentries = NEntries;}
  if (skipentries>1000000000) {skipentries=0;}


// diogenrange is the momentum range over which the DIO events were generated
  double nstopped(7.56e17);
  double capfrac(0.609); 
  double decayfrac = 1.0 - capfrac;
  double ndecay = nstopped*decayfrac;
  double ncap = nstopped*capfrac;
  double conprob(1e-15);
  double momlow(103.3);
  double momhigh(104.7);
  double trueconvmom(104.973);

  unsigned nbins(100);
  double mmin(101);
  double mmax(106);
  double mevperbin = (mmax-mmin)/nbins;

  int minnhits(50);
  double momErrScale(1.0/*0.8*//*0.61*/);

  bool plotDiO = dio!=0x0;

  double conscale = ncap*conprob/ncon;
  cout << "Conversion scale factor =" << conscale << endl;
// dio spectrum
  TF1* diocz_f = new TF1("diocz_f",DIOCZ,85.0,105,1);
  diocz_f->SetLineColor(kGreen);
  diocz_f->SetParameter(0,1.0);
// integrate the DIO spectrum over the range specified.  This is relative to the free decay rate
  //double dioSpcNorm = diocz_f->Integral(85,trueconvmom);
  double dioint = diocz_f->Integral(trueconvmom-diogenrange,trueconvmom);//dioSpcNorm;
  double dioscale(1.0);
  if (plotDiO) {
          if(weightdio){
                  dioscale =ndecay*diogenrange/ndio;
          } else {
                  dioscale = /*dioint**/ndecay/ndio;
                  //dioscale = 1.0/dioscale;
          }
  }
  //cout<<" Total Cz int "<<dioSpcNorm<<" dioint% "<<dioint<<endl;
  cout<<" Prob simul DIO range "<<dioint<<endl;
  cout << "DIO scale factor = " << dioscale << endl;

  TH1F* diospec[4];
  TH1F* conspec[4];
  TH1F* conspec_1[4];
  TH1F* conspec_1_1[4];
  TH1F* conspec_2[4];
// basic cuts
  TCut reco("fitinfo.fit==1&&recopar.parvec.m[2]>0");
  TCut onlyFrstSD("fitinfo.iseed==0");
  reco+=onlyFrstSD;
  TCut addcut(addCut.Data());
  double tdlow(0.57735027);
  double tdhigh(1.0);
  double t0min(710);
  char ctext[80];
  snprintf(ctext,80,"recopar.parvec.m[4]>%f&&recopar.parvec.m[4]<%f",tdlow,tdhigh);
  TCut pitch(ctext);
  snprintf(ctext,80,"fitinfo.t0fit>%f",t0min);
  TCut livegate(ctext);
  TCut goodMatCorr("(fitinfo.fitmombeam-fitinfo.fitmom)>0&&(fitinfo.fitmombeam-fitinfo.fitmom)<1.0");
  //snprintf(ctext,80,"fitinfo.fitmom>%f&&fitinfo.fitmom<%f",momlow,momhigh);
  //TCut momwin(ctext);
  double momlowTmp[4] = {momlow, momlow, momlow, momlow};
  double momhighTmp[4] = {momhigh, momhigh, momhigh, momhigh};

// cuts for different tightness of selection
  TCut ncuts[4], t0cuts[4], momcuts[4], fitcuts[4];
  ncuts[0] = Form("fitinfo.nhits>=%i",minnhits+10);
  ncuts[1] = Form("fitinfo.nhits>=%i",minnhits+10);
  ncuts[2] = Form("fitinfo.nhits>=%i",minnhits+20);
  ncuts[3] = Form("fitinfo.nhits>=%i",minnhits+20);//+30);
  t0cuts[0] = "fitinfo.errt0<3";
  t0cuts[1] = "fitinfo.errt0<2";
  t0cuts[2] = "fitinfo.errt0<1.2";//1
  t0cuts[3] = "fitinfo.errt0<1.2";//1
  momcuts[0] = Form("(%1.2f*fitinfo.fitmomerr)<0.4",momErrScale);
  momcuts[1] = Form("(%1.2f*fitinfo.fitmomerr)<0.3",momErrScale);
  momcuts[2] = Form("(%1.2f*fitinfo.fitmomerr)<0.2",momErrScale);
  momcuts[3] = Form("(%1.2f*fitinfo.fitmomerr)<0.2",momErrScale);
  //fitcuts[0] = "TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))>1e-6";
  //fitcuts[1] = "TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))>1e-4";
  //fitcuts[2] = "TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))>1e-3";
  //fitcuts[3] = "TMath::Prob(fitinfo.chi2,(fitinfo.nhits-5))>1e-2";
  fitcuts[0] = "fitinfo.fitcon>1e-6";//1e-4;
  fitcuts[1] = "fitinfo.fitcon>1e-5";//1e-3;
  fitcuts[2] = "fitinfo.fitcon>1e-4";//5e-3;
  fitcuts[3] = "fitinfo.fitcon>5e-3";//1e-2;

  TCut ncutsText[4], t0cutsText[4], momcutsText[4], fitcutsText[4];
  ncutsText[0] = Form("nactive>=%i",minnhits+10);
  ncutsText[1] = Form("nactive>=%i",minnhits+10);
  ncutsText[2] = Form("nactive>=%i",minnhits+20);
  ncutsText[3] = Form("nactive>=%i",minnhits+20);//+30);
  t0cutsText[0] = "t0err<3";
  t0cutsText[1] = "t0err<2";
  t0cutsText[2] = "t0err<1.2";//1
  t0cutsText[3] = "t0err<1.2";//1
  if (momErrScale>0.99) {
          momcutsText[0] = "fitmomerr<0.4";
          momcutsText[1] = "fitmomerr<0.3";
          momcutsText[2] = "fitmomerr<0.2";
          momcutsText[3] = "fitmomerr<0.2";
  } else {
          momcutsText[0] = Form("(%1.2f*fitmomerr)<0.4",momErrScale);
          momcutsText[1] = Form("(%1.2f*fitmomerr)<0.3",momErrScale);
          momcutsText[2] = Form("(%1.2f*fitmomerr)<0.2",momErrScale);
          momcutsText[3] = Form("(%1.2f*fitmomerr)<0.2",momErrScale);
  }

  fitcutsText[0] = "fitcon>1e-6";//1e-4;
  fitcutsText[1] = "fitcon>1e-5";//1e-3;
  fitcutsText[2] = "fitcon>1e-4";//5e-3;
  fitcutsText[3] = "fitcon>5e-3";//1e-2;


  for(unsigned ires=0;ires<4;ires++){
    char dioname[50];
    snprintf(dioname,50,"diospec%i",ires);
    char conname[50];
    snprintf(conname,50,"conspec%i",ires);
    char conname_1[50];
    snprintf(conname_1,50,"conspec_1_%i",ires);
    char conname_1_1[50];
    snprintf(conname_1_1,50,"conspec_1_1_%i",ires);
    char conname_2[50];
    snprintf(conname_2,50,"conspec_2_%i",ires);
 
    diospec[ires] = new TH1F(dioname,"Reconstructed momentum;MeV/c;Events per experiment",nbins,mmin,mmax);
    diospec[ires]->SetStats(0);
    diospec[ires]->SetLineColor(kBlue);
    diospec[ires]->Sumw2();

    conspec[ires] = new TH1F(conname,"Reconstructed momentum;MeV/c;Events per experiment",nbins,mmin,mmax);
    conspec[ires]->SetStats(0);
    conspec[ires]->SetLineColor(kRed);
    conspec[ires]->Sumw2();
    conspec_1[ires] = new TH1F(conname_1,"Reconstructed momentum;MeV/c;Events per experiment",nbins,mmin,mmax);
    conspec_1[ires]->SetStats(0);
    conspec_1[ires]->SetLineColor(kGreen);
    conspec_1[ires]->Sumw2();
    conspec_1_1[ires] = new TH1F(conname_1_1,"Reconstructed momentum;MeV/c;Events per experiment",nbins,mmin,mmax);
    conspec_1_1[ires]->SetStats(0);
    conspec_1_1[ires]->SetLineColor(kYellow-3);
    conspec_1_1[ires]->Sumw2();
    conspec_2[ires] = new TH1F(conname_2,"Reconstructed momentum;MeV/c;Events per experiment",nbins,mmin,mmax);
    conspec_2[ires]->SetStats(0);
    conspec_2[ires]->SetLineColor(kViolet);
    conspec_2[ires]->Sumw2();

    TCut quality = ncuts[ires] && t0cuts[ires] && momcuts[ires] && fitcuts[ires];
    TCut final = (reco+pitch+livegate+quality+addcut);
    //TCut fromIW ("fitinfo.trkfrstHitPosZ>-990");
    TCut fromIW ("fitinfo.frstHitPosZ>-990");

    if (useExtapMom==1) {
            final+=goodMatCorr;
    }

    if (plotDiO) {
    if(weightdio){
            if (useExtapMom==2) {
                    dio->Project(dioname,"fitinfo.fitmombeam","fitinfo.wgt"*(final+goodMatCorr),"",nentries,skipentries);
                    dio->Project(Form("+%s",dioname),"fitinfo.fitmom","fitinfo.wgt"*(final+!goodMatCorr),"",nentries,skipentries);
            } else if (useExtapMom==1) {
                    dio->Project(dioname,"fitinfo.fitmombeam","fitinfo.wgt"*final,"",nentries,skipentries);
            } else {
                    dio->Project(dioname,"fitinfo.fitmom","fitinfo.wgt"*final,"",nentries,skipentries);
            }
    } else {
            if (useExtapMom==2) {
                    dio->Project(dioname,"fitinfo.fitmombeam","fitinfo.wgt"*(final+goodMatCorr),"",nentries,skipentries);
                    dio->Project(Form("+%s",dioname),"fitinfo.fitmom","fitinfo.wgt"*(final+!goodMatCorr),"",nentries,skipentries);
            } else if (useExtapMom==1) {
                    dio->Project(dioname,"fitinfo.fitmombeam","fitinfo.wgt"*final,"",nentries,skipentries);
            } else {
                    dio->Project(dioname,"fitinfo.fitmom","fitinfo.wgt"*final,"",nentries,skipentries);
            }
            //diospec[ires]->Scale(ndecay);
    }
    diospec[ires]->Scale(dioscale*dioeff);
    }

    if (useExtapMom==2) {
            con->Project(conname_1,"fitinfo.fitmombeam",final+goodMatCorr,"",nentries,skipentries);
            con->Project(conname_1_1,"fitinfo.fitmombeam",final+goodMatCorr+fromIW,"",nentries,skipentries);
            con->Project(conname_2,"fitinfo.fitmom",final+!goodMatCorr,"",nentries,skipentries);
            conspec[ires]->Add(conspec_1[ires]);
            conspec[ires]->Add(conspec_2[ires]);
            conspec_1[ires]->Scale(conscale*coneff);
            conspec_1_1[ires]->Scale(conscale*coneff);
            conspec_2[ires]->Scale(conscale*coneff);
    } else if (useExtapMom==1) {
            con->Project(conname,"fitinfo.fitmombeam",final,"",nentries,skipentries);
    } else {
            con->Project(conname,"fitinfo.fitmom",final,"",nentries,skipentries);
    }
    conspec[ires]->Scale(conscale*coneff);

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

  TLegend* legCE = new TLegend(0.65,0.30,0.9,0.5);
  //legCE->SetTextSize(10);
  if (useExtapMom==2) {
          legCE->AddEntry(conspec_1[0],"CE coor for mat","L");
          legCE->AddEntry(conspec_2[0],"CE no coor","L");
          legCE->AddEntry(conspec_1_1[0],"CE coor for mat from IW","L");
  }

// linear scale
  TCanvas* lincan = new TCanvas("mu2elin","mu2e linear scale",1200,800);
  lincan->Clear();
  lincan->Divide(2,2);
  for(unsigned ires=0;ires<4;ires++){
    lincan->cd(ires+1);
    TH1* diocopy = diospec[ires]->DrawCopy();
    diocopy->SetMinimum(-0.2);
    diocopy->SetMaximum(5);
    conspec[ires]->Draw("same");
    if (useExtapMom==2) {
            conspec_1[ires]->Draw("same");
            conspec_1_1[ires]->Draw("same");
            conspec_2[ires]->Draw("same");
            legCE->Draw();
    }

    int istart = diospec[ires]->FindFixBin(momlow+0.5*mevperbin);
    int istop = diospec[ires]->FindFixBin(momhigh-0.5*mevperbin);
//    cout << "Integration low edge " << diospec[ires]->GetBinLowEdge(istart) << " for cut at " << momlow << endl;
//    cout << "Integration high edge " << diospec[ires]->GetBinLowEdge(istop)+mevperbin << " for cut at " << momhigh << endl;
    double dint_err, cint_err;
    double dint = diospec[ires]->IntegralAndError(istart,istop,dint_err);
    double cint = conspec[ires]->IntegralAndError(istart,istop,cint_err);

    TPaveText* inttext = new TPaveText(0.5,0.65,0.9,0.8,"NDC");
    char itext[50];
    snprintf(itext,50,"%4.2f MeV/c < P < %4.2f MeV/c",momlow,momhigh);
    inttext->AddText(itext);
    snprintf(itext,50,"DIO integral = %5.3f #pm %4.3f",dint,dint_err);
    inttext->AddText(itext);
    snprintf(itext,50,"Conv. integral = %5.3f #pm %4.3f",cint,cint_err);
    inttext->AddText(itext);
    inttext->Draw();

    if (plotDiO && dint<0.24) {
          while (/*(dint-2.0*dint_err)<0.25*/dint<0.24 && istart>1) {
                --istart;
                dint = diospec[ires]->Integral(istart,istop);
          }
          ++istart;
          dint = diospec[ires]->IntegralAndError(istart,istop,dint_err);
          cint = conspec[ires]->IntegralAndError(istart,istop,cint_err);

          double momlowNew(diospec[ires]->GetBinLowEdge(istart));
          TPaveText* inttextNew = new TPaveText(0.5,0.50,0.9,0.65,"NDC");
          inttextNew->AddText(Form("%4.2f MeV/c < P < %4.2f MeV/c",momlowNew,momhigh));
          inttextNew->AddText(Form("DIO integral = %5.3f #pm %4.3f",dint,dint_err));
          inttextNew->AddText(Form("Conv. integral = %5.3f #pm %4.3f",cint,cint_err));
          inttextNew->Draw();

          TLine* momlowl_new = new TLine(momlowNew,0.0,momlowNew,1.5*conspec[ires]->GetBinContent(conspec[ires]->GetMaximumBin()));
          momlowl_new->SetLineColor(kBlack);
          momlowl_new->SetLineStyle(2);
          momlowl_new->SetLineWidth(2);
          momlowl_new->Draw();

          momlowTmp[ires] = momlowNew;
    }

    TPaveText* lintext = new TPaveText(0.1,0.5,0.4,0.8,"NDC");  
    char line[40];
    snprintf(line,80,"%4.3f<tan(#lambda)<%4.3f",tdlow,tdhigh);
    lintext->AddText(line);
    snprintf(line,80,"t0>%5.1f nsec",t0min);
    lintext->AddText(line);
    sprintf(line,"%s",ncutsText[ires].GetTitle());
    lintext->AddText(line);
    sprintf(line,"%s",t0cutsText[ires].GetTitle());
    lintext->AddText(line);
    sprintf(line,"%s",momcutsText[ires].GetTitle());
    lintext->AddText(line);
    sprintf(line,"%s",fitcutsText[ires].GetTitle());
    lintext->AddText(line);
    if (!addCut.IsNull()) {
            TString addCutText = addCut;
            addCutText.ReplaceAll("fitinfo.","");
            addCutText.ReplaceAll(">","\>");
            addCutText.ReplaceAll("<","\<");
            int iSubs=0;
            TString subAddCut[10];
            while (addCutText.Contains("&") && iSubs<10) {
                    Ssiz_t andpos = addCutText.First("&");
                    TSubString splitAddCut = addCutText(0,andpos);
                    subAddCut[iSubs]=splitAddCut;
                    cout<<"andpos "<< andpos<<" "<<splitAddCut<<endl;
                    if (andpos>0) {
                            lintext->AddText(subAddCut[iSubs].Data());
                            addCutText=addCutText.Remove(0,andpos+2);
                    }
                    ++iSubs;
            }
            cout<<"addCutText "<<addCutText<<endl;
      lintext->AddText(addCutText.Data());
    }
    lintext->Draw();

    TLine* momlowl = new TLine(momlow,0.0,momlow,1.5*conspec[ires]->GetBinContent(conspec[ires]->GetMaximumBin()));
    momlowl->SetLineColor(kBlack);
    momlowl->SetLineStyle(2);
    momlowl->SetLineWidth(2);
    momlowl->Draw();

    TLine* momhighl = new TLine(momhigh,0.0,momhigh,1.5*conspec[ires]->GetBinContent(conspec[ires]->GetMaximumBin()));
    momhighl->SetLineColor(kBlack);
    momhighl->SetLineStyle(2);
    momhighl->SetLineWidth(2);
    momhighl->Draw();
    leg->Draw();
    info->Draw();
 
  }
  lincan->cd(0);
  lincan->SaveAs("mu2e_lin.png");

  if (plotDiO) {
  TCanvas* dioc = new TCanvas("dio","dio",1200,800);
  dioc->Divide(2,2);

  Double_t dmhi = trueconvmom;
  Double_t dmlow = trueconvmom - diogenrange;
  TH1F* diogen = new TH1F("diogen","True DIO momentum;MeV",nbins,dmlow,dmhi);
  TH1F* diowt = new TH1F("diowt","True DIO momentum;MeV",nbins,dmlow,dmhi);
//  diowt->Sumw2();
  dio->Project("diogen","fitinfo.genmom","","",nentries,skipentries);
  double diofscale(1.0);
  if(weightdio){
          dio->Project("diowt","fitinfo.genmom","fitinfo.wgt","",nentries,skipentries);
          // dead-reconing on spectrum, accounting for bins
          diofscale = ndecay*(dmhi-dmlow)/nbins;
  } else {
          dio->Project("diowt","fitinfo.genmom","fitinfo.wgt","",nentries,skipentries);
          //diowt->Scale(ndecay);
          //diowt->Scale(1.0/dioscale);
          diofscale = ndecay;//dioint*dioscale;
          //diofscale*=dioscale;
  }
  diowt->Scale(dioscale);
  diowt->SetLineColor(kBlue);
  diogen->SetLineColor(kRed);
  diowt->SetStats(0);
  diogen->SetStats(0);
  dioc->cd(1);
  gPad->SetLogy();
  diocz_f->SetParameter(0,diofscale);
  diowt->Draw();
  diocz_f->Draw("same");
  diogen->Draw("same");
  TLegend* dioleg = new TLegend(.2,.4,.6,.6);
  dioleg->AddEntry(diogen,"Generated","l");
  dioleg->AddEntry(diowt,"Weighted","l");
  dioleg->AddEntry(diocz_f,"Czarnecki etal","l");
  dioleg->Draw();

  dioc->cd(3)->SetLogy(kTRUE);
//  gPad->SetLogy();
  Int_t colors[4] = {kRed,kBlue,kGreen,kBlack};
  TH1F* diogenwin[4] = {0,0,0,0};
  const char* dopt[4] = {"","same","same","same"};
  const char* cutset[4] = {"Cutset A","Cutset B","Cutset C","Cutset D"};
  TLegend* dgenwinleg = new TLegend(.5,.6,.7,.9);
  TLegend* dgenwinPROMleg = new TLegend(.1,.6,.45,.9);
  diogenwin[0] = new TH1F("diogenwin_0","True momentum of DIO in signal box;MeV",100,dmlow,dmhi);
  for(unsigned ires=0;ires<4;ires++){
    char dioname[50];
    snprintf(dioname,50,"diogenwin%i",ires);
    diogenwin[ires] = new TH1F(dioname,"True momentum of DIO in signal box;MeV",100,dmlow,dmhi);
    diogenwin[ires]->SetStats(0);
//   TH1F* diogood[ires] = new TH1F("diogood","True DIO momentum",100,dmlow,dmhi);
//    dio->Project("diogoodwt","fitinfo.genmom",goodfit);

    TCut quality = ncuts[ires] && t0cuts[ires] && momcuts[ires] && fitcuts[ires];
    TCut final = (reco+pitch+livegate+quality+addcut);
    TCut momwin(Form("fitinfo.fitmom>=%f&&fitinfo.fitmom<=%f",momlowTmp[ires],momhighTmp[ires]));
    TCut momwinCorr(Form("fitinfo.fitmombeam>=%f&&fitinfo.fitmombeam<=%f",momlowTmp[ires],momhighTmp[ires]));
    //final.Print();
    if (useExtapMom==2) {
            dio->Project(dioname,"fitinfo.genmom",final+goodMatCorr+momwinCorr,"",nentries,skipentries);
            dio->Project(Form("+%s",dioname),"fitinfo.genmom",final+!goodMatCorr+momwin,"",nentries,skipentries);
    } else if (useExtapMom==1) {
            dio->Project(dioname,"fitinfo.genmom",final+momwinCorr,"",nentries,skipentries);
    } else {
            dio->Project(dioname,"fitinfo.genmom",final+momwin,"",nentries,skipentries);
    }
    diogenwin[ires]->SetFillColor(colors[ires]);
    dgenwinleg->AddEntry(diogenwin[ires],cutset[ires],"f");
    diogenwin[ires]->Draw(dopt[ires]);

    int binmin = (int)((momlowTmp[ires]-diogenwin[ires]->GetBinLowEdge(1))/diogenwin[ires]->GetBinWidth(1));
    //int binmax = (int)((momhighTmp[ires]-diogenwin[ires]->GetBinLowEdge(1))/diogenwin[ires]->GetBinWidth(1));
    double nPromtotedDIO = diogenwin[ires]->Integral(1,binmin);
    double promotedDIOfrac = nPromtotedDIO/diogenwin[ires]->GetEntries();
    //cout<<"for cut "<<ires<<" momlow = "<<momlowTmp[ires]<<" momhigh = "<<momhighTmp[ires]<<endl;
    //std::cout<<"DIO fraction with gen momentum out of reconstruction window (cut "<<ires<<") = "<<promotedDIOfrac<<std::endl;
    dgenwinPROMleg->AddEntry(diogenwin[ires],Form("DIO promoted %.2f%%",(promotedDIOfrac*100.0)),"f");
    //std::cout<<"check : entries = "<<diogenwin[ires]->GetEntries()<<" bin to minp = "<<binmin<<" integral from 1 to minp = "<<nPromtotedDIO<<" % tot integ = "<<diogenwin[ires]->Integral(1,100)/diogenwin[ires]->GetEntries()<<endl;
  }
  dgenwinleg->Draw();
  dgenwinPROMleg->Draw();

  dioc->cd(2);
  for(unsigned ires=0;ires<4;ires++){
    diospec[ires]->SetFillColor(colors[ires]);
    if(ires==0)
      diospec[ires]->Draw("Hist");
    else
      diospec[ires]->Draw("Histsame");
  }
  dgenwinleg->Draw();
  TLine* momlowl = new TLine(momlow,0.0,momlow,1.5*diospec[0]->GetBinContent(diospec[0]->GetMaximumBin()));
  momlowl->SetLineColor(kBlack);
  momlowl->SetLineStyle(2);
  momlowl->SetLineWidth(2);
  momlowl->Draw();

  TLine* momhighl = new TLine(momhigh,0.0,momhigh,1.5*diospec[0]->GetBinContent(diospec[0]->GetMaximumBin()));
  momhighl->SetLineColor(kBlack);
  momhighl->SetLineStyle(2);
  momhighl->SetLineWidth(2);
  momhighl->Draw();

  dioc->SaveAs("diocan.png");
  }
}

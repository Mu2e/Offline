#include <iostream>
#include <iomanip>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TH2Poly.h>
#include <TVector3.h>
#include <TVector2.h>
#include <TAxis.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>
#include <TMinuit.h>
#include <TFitter.h>
#include <TFitter.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TEllipse.h>
#include <TArrow.h>

#include <vector>
#include <cmath>
#include <set>
#include <list>


void     fillFinal(double* par);
void     fitFinal(TString title, TH1F *h, int ichoice);

void     fcn(int& npar, double* deriv, double& f, double par[], int flag);
double   getChi2(double* par);
TVector3 calcCog(int ic, double *par);
TVector3 getPosTrk(int it, double dlen);
bool     checkVD(bool wantFace=false, bool want400=false, bool wantSide=false);
double   calcLen(int ic);

TH2Poly* buildDisk(TString hname,double _radiusIn, double _radiusOut, double _cellSize);
void     getXYPosition_disk(int index, double&x, double& y) ;
bool     isInsideDisk(double x, double y, double _cellSize, double _radiusIn, double _radiusOut) ;
double   calcDistToSide(TVector2& a, TVector2& b) ;

void     setTree();
void     makeHist();
void     histHelper(TH1F* hist, TString title);
void     visualize();

double   Tail(double *x, double *par);
double   PartGau(double *x, double *par);
double   CBall(double *x, double *par);
double   fitFun(double *x, double *par);
double   FHWM(TF1 *func, double xmin, double xmax);


TTree *tree;

Int_t         evt;

Int_t         nMatch;
Int_t         mTrkId[1024];
Int_t         mCluId[1024];
Float_t       mChi2[1024];
Float_t       mChi2Pos[1024];
Float_t       mChi2Time[1024];

Int_t         nTrk;
Float_t       trkX[1024];
Float_t       trkY[1024];
Float_t       trkZ[1024];
Float_t       trkFFX[1024];
Float_t       trkFFY[1024];
Float_t       trkFFZ[1024];
Float_t       trkt[1024];
Float_t       trke[1024];
Float_t       trkpX[1024];
Float_t       trkpY[1024];
Float_t       trkpZ[1024];
Int_t         trknHit[1024];
Int_t         trkStat[1024];
Float_t       trkprob[1024];
Float_t       trkd0[1024];
Float_t       trkz0[1024];
Float_t       trkphi0[1024];
Float_t       trkomega[1024];
Float_t       trkcdip[1024];
Float_t       trkdlen[1024];

Int_t         nCluster;
Float_t       cluEnergy[1024];
Float_t       cluTime[1024];
Float_t       cluCogX[1024];
Float_t       cluCogY[1024];
Float_t       cluCogZ[1024];
Float_t       clu2Mom[1024];
Float_t       cluE25[1024];
vector<vector<int> > *cluList;

Int_t         nCry;
Int_t         cryId[1024];   
Int_t         crySecId[1024];   
Float_t       cryPosX[1024];   
Float_t       cryPosY[1024];   
Float_t       cryPosZ[1024];  
Float_t       cryEdep[1024];  
Float_t       cryTime[1024];   
Int_t         nVd;
Int_t         vdId[1024];   
Int_t         vdPdgId[1024];   
Float_t       vdMom[1024];  
Float_t       vdMomX[1024];   
Float_t       vdMomY[1024];   
Float_t       vdMomZ[1024];   
Float_t       vdPosX[1024];   
Float_t       vdPosY[1024];   
Float_t       vdPosZ[1024];   
Float_t       vdTime[1024];   
Int_t         vdSimId[1024];   




TH2Poly *_dedep[30];
int numPoly=0;
double visX[100],visY[100];
double visPX[100],visPY[100];
double visCX[100],visCY[100];
int visN[100];
double visLX[100][100];
double visLY[100][100];

TH1F *_hx,*_hy,*_hxtrk,*_hytrk,*_hxy,*_hxytrk,*_hu,*_hutrk,*_hv,*_hvtrk,*_hr,*_huc,*_hvc,*_hu0,*_hv0,*_ha,*_hncr,*_prob1,*_prob2,*_chi2prob;
TH2F *_hxy2,*_hu2nd,*_hv2nd,*_huE25,*_huDist,*_huRad,*_hucDip,*_prob2d;
int _dispCount = 0;
TGraphErrors *gr;



int    WeightMode = 0;
double RatioEcut  = 0.0;
int    NevtFit    = 10000;
bool   checkFace  = false;
bool   checkSide  = false;

TString modelName[5]={"linear","sqrt","log"};


double duOffset = 0;
double dvOffset = 0;



void FitCog(int wmode = 0)
{
     WeightMode = wmode;

     TFile file("dataMu2e/ElecTrain/tt.root");
     tree = (TTree *)file.Get("ReadTrackCaloMatching/trkClu");
     tree->SetMakeClass(1);
     tree->LoadBaskets();
     cluList = 0;

     setTree();     
     makeHist();
    
    
     double param[10];
     int npar(5);

     TFitter minuit(99);
     minuit.SetFCN(fcn);

     switch (WeightMode) {

       //linear
       case 0:
	 minuit.SetParameter(0, "par0",  0., 1.,  0.0, 0.0);
	 //minuit.SetParameter(1, "par1",  1., 1.,  0.0, 0.0);
	 npar=1;
	 break;

       //sqrt(x)
       case 1:
	 minuit.SetParameter(1, "par1",  1., 1., -20.0, 20.0);
	 npar=1;
	 break;

       //log(x)
       case 2:
	 minuit.SetParameter(0, "par0",  0.1, 0.01, 0,0);
	 minuit.SetParameter(1, "par1",  1, 0.01, 0.001, 100.0);
	 npar=2;
	 break;

       default :
         cout<<"No wieght "<< WeightMode<<" available"<<endl;
	 return;
     }


     double argF[2] = {2000,0.01};
     minuit.ExecuteCommand("MIGRAD",argF,2);
     for (int i=0; i<npar; i++)  param[i] = minuit.GetParameter(i);
     fillFinal(param);



//double param2[2]={-10.01,0};
//double param2[2]={0,0};
//fillFinal(param2);
     
     cout<<endl<<endl;
     fitFinal("Fit U0  ",_hu0,3);
     fitFinal("Fit V0  ",_hv0,3);
     fitFinal("Fit R   ",_hr,4);
     
     cout<<endl;
     fitFinal("Fit U   ",_hu,3);
     fitFinal("Fit V   ",_hv,3);


     cout<<endl;
     fitFinal("Fit X    ",_hx,1);
     fitFinal("Fit Y    ",_hy,1);
     fitFinal("Fit X+Y  ",_hxy,1);
     cout<<endl<<endl;





     TCanvas *c1 = new TCanvas("c1","c1",800,600);
     c1->Divide(2);
     c1->cd(1); _hx->Draw("e");
     c1->cd(2); _hy->Draw("e");
     c1->SaveAs(Form("plots/FitCogNew%i1.eps",wmode));

     TCanvas *c2 = new TCanvas("c2","c2",800,600);
     _hxy2->Draw("colz");
     c2->SaveAs(Form("plots/FitCogNew%i2.eps",wmode));

 
     TCanvas *c4 = new TCanvas("c4","c4",800,600);
     c4->Divide(2);
     c4->cd(1); _hu->Draw("e");
     c4->cd(2); _hv->Draw("e");
     c4->SaveAs(Form("plots/FitCogNew%i4.eps",wmode));


     TCanvas *c5 = new TCanvas("c5","c5",800,600);
     _hxy->Draw("e");
     c5->SaveAs(Form("plots/FitCogNew%i4a.eps",wmode));
     _hv->Draw("e");
     c5->SaveAs(Form("plots/FitCogNew%i4b.eps",wmode));



//visualize();
     
     TFile fres("hres.root","RECREATE");
     fres.Add(_hx);
     fres.Add(_hy);
     fres.Add(_hxy);
     fres.Add(_hxy2);
     fres.Add(_hu0);
     fres.Add(_hv0);
     fres.Add(_hu);
     fres.Add(_hr);
     fres.Add(_hv);
     fres.Add(_huc);
     fres.Add(_hvc);
     fres.Add(_hu2nd);
     fres.Add(_hv2nd);
     fres.Add(_huE25);
     fres.Add(_huDist);
     fres.Add(_hucDip);
     fres.Add(_huRad);
     fres.Add(_ha);
     fres.Add(_hncr);
     fres.Add(_prob1);
     fres.Add(_prob2);
     fres.Add(_prob2d);
     fres.Add(_chi2prob);
     for (int i=0;i<19;i++)fres.Add(_dedep[i]);
     fres.Write();
     fres.Close(); 
 
}














void fillFinal(double* par)
{
     
    int nevent = tree->GetEntries();

    for (int i=0;i<nevent;i++)
    {
        tree->GetEntry(i);
        for (int im=0;im<nMatch;im++)
        {
	   if (mChi2[im]>100) continue;
	   int it = mTrkId[im];
	   int ic = mCluId[im];


	   if (cluEnergy[ic]<50) continue;
	   if (trkStat[it]!=1 || trkprob[it]<0.01 || trknHit[it] < 15) continue;

           //if (checkFace  && trkFFZ[it] >1 ) continue;
           //if (checkSide  && trkFFZ[it] <1 ) continue;



	   double distShowerMax = calcLen(ic);

	   double trkp  = sqrt(trkpX[it]*trkpX[it] + trkpY[it]*trkpY[it]+ trkpZ[it]*trkpZ[it]);
	   double trkpt = sqrt(trkpX[it]*trkpX[it]+trkpY[it]*trkpY[it]);
	   TVector2 du(trkpX[it]/trkpt,trkpY[it]/trkpt);
	   TVector2 dv(-trkpY[it]/trkpt,trkpX[it]/trkpt);
	   TVector2 dx_inuv(trkpX[it]/trkpt,-trkpY[it]/trkpt);
	   TVector2 dy_inuv(trkpY[it]/trkpt,trkpX[it]/trkpt);



           //TVector3 trkPos0 = getPosTrk(it,0);
           TVector3 trkPos0(trkFFX[it],trkFFY[it],0);
           TVector3 trkDir(trkpX[it]/trkp,trkpY[it]/trkp,trkpZ[it]/trkp);
	   
	   TVector3 cog = calcCog(ic, par);
	   TVector3 cluDir(0,0,1);

           TVector3 initDiff = cog-trkPos0;
	   if ( sqrt(initDiff.X()*initDiff.X()+initDiff.Y()*initDiff.Y()) > 200) continue;

	   
	 

	   //min distan0 ce in xy cordinate
	   double dlenT(0);
           TVector3 diff = trkPos0-cog ;
           double distMin = sqrt(diff.X()*diff.X()+diff.Y()*diff.Y());
	   for (double dlent=0;dlent<200;dlent+=1.0)
	   {
	        //TVector3 trkPos = getPosTrk(it,dlent);
	        TVector3 trkPos = trkPos0+dlent*trkDir;
	        diff = trkPos-cog;
                double dist = sqrt(diff.X()*diff.X()+diff.Y()*diff.Y());
                if (dist < distMin) {distMin=dist;dlenT=dlent;}
	   }
           
         
           TVector3 TrkPosMin =  trkPos0+dlenT*trkDir;
	   double dx = cog.X() - TrkPosMin.X();
	   double dy = cog.Y() - TrkPosMin.Y();
           TVector2 xy(dx,dy);
	   
           TVector2 deltaUV0(xy*du,xy*dv);
	   TVector2 deltaUV( xy*du-duOffset,xy*dv-dvOffset);
	   TVector2 deltaXY(dx,dy);
            _ha->Fill(dlenT);
	   
/*         
	   
	   double dlenC(0),dlenT(0);
           double distMinGen = (trkPos0-cog).Mag();
	   for (double dlent=0;dlent<100;dlent+=1.0)
	   {
	        TVector3 trkPos = getPosTrk(it,dlent);
	        //TVector3 trkPos = trkPos0+dlent*trkDir;
	        for (double dlenc=0;dlenc<100;dlenc+=1.0)
	        {
	            TVector3 cluPos = cog + (dlenc*cluDir);
	            double dist = (cluPos-trkPos).Mag();
		    if (dist < distMinGen) {distMinGen=dist;;dlenT=dlent;;dlenC=dlenc;}
	        }
	   }

	   TVector3 TrkPosMin =  trkPos0+dlenT*trkDir;
	   double dx = cog.X() - TrkPosMin.X();
	   double dy = cog.Y() - TrkPosMin.Y();
           TVector2 xy(dx,dy);
	   
           TVector2 deltaUV0(xy*du,xy*dv);
	   TVector2 deltaUV( xy*du-duOffset,xy*dv-dvOffset);
	   TVector2 deltaXY(dx,dy);
*/
	   
	    

	   _hx->Fill(deltaXY.X());
	   _hy->Fill(deltaXY.Y());
	   _hxy->Fill(deltaXY.X());
	   _hxy->Fill(deltaXY.Y());
	   _hxy2->Fill(deltaXY.X(),deltaXY.Y());
	   _hr->Fill(deltaXY.Mod2());
	   _hu0->Fill(deltaUV0.X());
	   _hv0->Fill(deltaUV0.Y());
	   _hu->Fill(deltaUV.X());
	   _hv->Fill(deltaUV.Y());
          


	   _hu2nd->Fill(clu2Mom[ic]/1000,deltaUV0.X());
	   _hv2nd->Fill(clu2Mom[ic]/1000,deltaUV0.Y());
	   _huE25->Fill(cluE25[ic]/cluEnergy[ic],deltaUV0.X());
	   _huDist->Fill(distShowerMax,deltaUV0.X());
	   _hucDip->Fill(trkcdip[it],deltaUV0.X());
	   _huRad->Fill(sqrt(cog.X()*cog.X()+cog.Y()*cog.Y()),deltaUV0.X());
	  



           int Ncr(0);
	   std::vector<int> crystalL = (*cluList)[ic];               
	   for (unsigned int il=0;il<crystalL.size();++il) if (cryEdep[crystalL.at(il)] > 5) ++Ncr;
	   _hncr->Fill(Ncr);



	   if (deltaUV.X() >25 && numPoly<30)
	   {
		visX[numPoly]=TrkPosMin.X();
		visY[numPoly]=TrkPosMin.Y();
		visPX[numPoly]=trkpX[it]/trkpt;
		visPY[numPoly]=trkpY[it]/trkpt;	       
		visCX[numPoly]=cog.X();
		visCY[numPoly]=cog.Y();

		for (int il=0;il<nCry;++il) if (cryEdep[il] > 0.9)  _dedep[numPoly]->Fill(cryPosX[il],cryPosY[il],cryEdep[il]);

		std::vector<int> crystalList = (*cluList)[ic];               
		visN[numPoly]=crystalList.size();              
		for (unsigned int il=0;il<crystalList.size();++il)
		{	      
	      	   int icry = crystalList.at(il);
		   visLX[numPoly][il]=cryPosX[icry];
	           visLY[numPoly][il]=cryPosY[icry];
		} 

		++numPoly;	
	   }  	   
	  
        }
    }
         
    return;
}



void fitFinal(TString title, TH1F *h, int ichoice)
{
     cout<<title;
     if (ichoice==0)
     {
       h->Fit("gaus","q","",-20,20);
       cout<<"Results "<<setprecision(3)<<h->GetFunction("gaus")->GetParameter(1)<<" +- "<<h->GetFunction("gaus")->GetParError(1)<<"     "
           <<h->GetFunction("gaus")->GetParameter(2)<<" +- "<<h->GetFunction("gaus")->GetParError(2)<<endl;
     } 

     if (ichoice==1)
     {       
       TF1 func("f","gaus(0)+gaus(3)");
       func.SetParameters(0.5*h->GetMaximum()/2.5,0,6.0,0.5*h->GetMaximum()/5.5,10,1);
       h->Fit("f","q");
       double fwhm = FHWM(&func,-30,30);
      
       cout<<"Results "<<setprecision(3)<<func.GetParameter(2)<<" +- "<<func.GetParError(2)<<"  /  "
                       <<func.GetParameter(5)<<" +- "<<func.GetParError(5)<<"     FWHM/2.34   ="<<fwhm<<" RMS "<<h->GetRMS()<<endl;


 
     }

     if (ichoice==2)
     {
       h->Fit("gaus","q");
       cout<<"Results "<<setprecision(3)<<h->GetFunction("gaus")->GetParameter(2)<<" +- "<<h->GetFunction("gaus")->GetParError(2)<<"  /  "
                       <<h->GetFunction("gaus")->GetParameter(1)<<" +- "<<h->GetFunction("gaus")->GetParError(1)<<endl;
     } 

     if (ichoice==3)
     {
       TF1 func("fitFun",fitFun,h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),5);
       func.SetParameters(h->GetSumOfWeights()/6/h->GetRMS(),h->GetBinCenter(h->GetMaximumBin()),h->GetRMS(),5,1);
       func.SetParLimits(3,0,50);
       h->Fit("fitFun","q");
       double fwhm = FHWM(&func,-30,30);
       cout<<"Results "<<setprecision(3)<<func.GetParameter(2)<<" +- "<<func.GetParError(2)<<"  /  "
            <<func.GetParameter(1)<<" +- "<<func.GetParError(1)<<"     FWHM/2.34   ="<<fwhm<<" RMS "<<h->GetRMS()<<endl;

 
     } 

     if (ichoice==4)
     {       
       TF1 func("f","expo(0)+expo(2)");
       func.SetParameters(8,-0.015,5.7,-0.004);
       h->Fit("f","q");
             double fwhm = FHWM(&func,-30,30);

       cout<<"Results "<<setprecision(5)<<func.GetParameter(0)<<"  "<<func.GetParameter(1)<<"  /  "
                       <<func.GetParameter(2)<<"  "<<func.GetParameter(3)<<"     FWHM/2.34   ="<<fwhm<<" RMS "<<h->GetRMS()<<endl;
;
     }
}



void visualize()
{
     
    TCanvas *ce = new TCanvas("ce","ce",1200,1200);

    for (int i=0;i<20;i++)
    {
       _dedep[i]->Draw("colz");
       TEllipse *el = new TEllipse(visX[i],visY[i],10,10);
       el->SetFillStyle(3000);
       el->SetFillColor(1);
       el->Draw();

       TEllipse *el2 = new TEllipse(visCX[i],visCY[i],10,10);
       el2->SetFillStyle(3000);
       el2->SetFillColor(5);
       el2->Draw();

       for (int ic=0;ic<visN[i];++ic) {
	 TEllipse *el3 = new TEllipse(visLX[i][ic],visLY[i][ic],2,2);
         el3->SetFillStyle(3000);
         el3->SetFillColor(1);
         el3->Draw();
       }

       TArrow *ar = new TArrow(visX[i],visY[i],visX[i]+50*visPX[i],visY[i]+50*visPY[i],0.01);
       ar->Draw();

       ce->SaveAs(Form("plots/vis%i.eps",i));
    }

}


double FHWM(TF1 *func, double xmin, double xmax)
{

  double step(0.01);
  
  double vmax(-1),xvmax(0);
  for (double x=xmin;x<xmax;x+=step)  {
    double val = func->Eval(x);
    if (val > vmax) {vmax = val ; xvmax=x;}
  }
      
  double vlim1(0),vlim2(0);
  for (double x=xmin;x<xmax;x+=step) { 
    double val = func->Eval(x);
    if (val > (vmax/2.0)) {vlim1 = x ; break;}
  }
  
  for (double x=xmax;x>xmin;x-=step) { 
    double val = func->Eval(x);
    if (val > (vmax/2.0)) {vlim2 = x+step ; break;}
  }

  return (vlim2-vlim1)/2.34;
}







void fcn(int& npar, double* deriv, double& f, double par[], int flag)
{
  f=getChi2(par);
  //cout<<par[0]<<" "<<" "<<f<<endl;
}                        



double getChi2(double* par)
{
   
    double chi2(0);

    int nevent = tree->GetEntries();
    if (NevtFit > 0 && NevtFit < nevent) nevent = NevtFit;

    for (int i=0;i<nevent;i++)
    {
       tree->GetEntry(i);
       for (int im=0;im<nMatch;im++) {

	 if (mChi2[im]>100) continue;
	 int it = mTrkId[im];
	 int ic = mCluId[im];


	 if (cluEnergy[ic]<60) continue;
	 if (trkStat[it]!=1 || trkprob[it]<0.001) continue;

         if (checkFace  && trkFFZ[it] >1 ) continue;
         if (checkSide  && trkFFZ[it] <1 ) continue;


	 double trkp = sqrt(trkpX[it]*trkpX[it] + trkpY[it]*trkpY[it]+ trkpZ[it]*trkpZ[it]);
	 TVector3 trkDir(trkpX[it]/trkp,trkpY[it]/trkp,trkpZ[it]/trkp);
         TVector3 trkPos0 = getPosTrk(it,0);
         
	 TVector3 cog = calcCog(ic, par);
	 TVector3 cluDir(0,0,1);



	 //min distance in xy cordinate
         TVector3 diff = trkPos0-cog ;
         double distMin = sqrt(diff.X()*diff.X()+diff.Y()*diff.Y());
	 for (double dlent=0;dlent<100;dlent+=1.0)
	 {
	      //TVector3 trkPos = getPosTrk(it,dlent);
	      TVector3 trkPos = trkPos0+dlent*trkDir;
	      diff = trkPos-cog;
              double dist = sqrt(diff.X()*diff.X()+diff.Y()*diff.Y());
              if (dist < distMin) distMin=dist;
	 }

	 
         
/*       
         // min distance in 3d
         double distMin = (trkPos0-cog).Mag();
	 for (double dlent=0;dlent<100;dlent+=1.0)
	 {
	      TVector3 trkPos = getPosTrk(it,dlent);
	      //TVector3 trkPos = trkPos0+dlent*trkDir;
	      for (double dlenc=0;dlenc<100;dlenc+=1.0)
	      {
	          TVector3 cluPos = cog + (dlenc*cluDir);
	          double dist = (cluPos-trkPos).Mag();
		  if (dist < distMin) distMin=dist;
	      }
	 }
*/         


	 chi2+=distMin/100;

       }
    } 

    return chi2;
}

TVector3 getPosTrk(int it, double dlen)
{
	 double sinDip  = sqrt(1-trkcdip[it]*trkcdip[it]);
	 double length  = (trkZ[it] - trkz0[it])/sinDip + dlen;
	 double posXTrk =  sin(trkphi0[it] + trkomega[it]*trkcdip[it]*length)/trkomega[it] - (trkd0[it] + 1.0/trkomega[it])*sin(trkphi0[it]);
	 double posYTrk = -cos(trkphi0[it] + trkomega[it]*trkcdip[it]*length)/trkomega[it] + (trkd0[it] + 1.0/trkomega[it])*cos(trkphi0[it]);
	 double posZtrk = dlen/sinDip;
         
	 return TVector3(posXTrk,posYTrk,posZtrk);
}


TVector3 calcCog(int ic, double *par)
{

      TVector3 aVector(0,0,0);
      double sumWeights(0);    

      std::vector<int> crystalList = (*cluList)[ic];

      for (unsigned int il=0;il<crystalList.size();++il) 
      {         
	int icry = crystalList.at(il);

	double energy    = cryEdep[icry];
	if (energy < 1e-6) continue;
//	if (energy < 3) continue;

	//weight function
	double weight(0),arg(0);
	double eRatio = energy/cluEnergy[ic];
	

	switch(WeightMode)
	{
           case 0: 
	       if ( eRatio> RatioEcut) weight = par[0]+energy;
	       break;

           case 1: 
	       if ( eRatio > RatioEcut && (par[0]+energy)>1e-5) weight = sqrt(par[0]+energy); 
	       break;

           case 2: 
	       if ( eRatio > RatioEcut) weight = par[0]+par[1]*log(energy); 
	       break;

	}
	if (weight < 0) weight=0;
        

	aVector[0] += cryPosX[icry]*weight;
	aVector[1] += cryPosY[icry]*weight;
	sumWeights += weight;
     }

     if (sumWeights>1e-3) { aVector[0] /= sumWeights; aVector[1] /= sumWeights;}
     else                 { aVector[0]=aVector[1]=0;}

     return aVector;   
}

double calcLen(int ic)
{

      double distMax(0);
      double Emin = 1;
      
      std::vector<int> crystalList = (*cluList)[ic];

      for (unsigned int i=1;i<crystalList.size();++i) 
      {         
	int icry = crystalList.at(i);
	if (cryEdep[icry]<Emin) continue;
       
        for (unsigned int j=0;j<i;++j)
	{
 	   int jcry = crystalList.at(j);
   	   if (cryEdep[jcry]<Emin) continue;

	   double dx = cryPosX[icry]-cryPosX[jcry];
	   double dy = cryPosY[icry]-cryPosY[jcry];
	   double dist = sqrt(dx*dx+dy*dy);
	   if (dist > distMax) distMax = dist;
	}
      }	
     return distMax;   
}

bool checkVD(bool wantFace, bool want400, bool wantSide)
{
     int iVD(-1);
     for (int iv=0;iv<nVd;++iv) if (vdPdgId[iv]==11 && vdMom[iv] > 90) {iVD = iv; break;}
     if (iVD==-1) return true;

     bool   isFace = vdId[iVD]==73 || vdId[iVD]==75;
     bool isInside = vdId[iVD]==77 || vdId[iVD]==79;
     bool  isDisk0 = vdId[iVD]==73 || vdId[iVD]==77;

     if (wantFace && !isFace) return false;
     if (want400 && sqrt(vdPosX[iVD]*vdPosX[iVD]+vdPosY[iVD]*vdPosY[iVD])<400) return false;
     if (wantSide && !isInside) return false;
     return true;
}




double fitFun(double *x, double *par) {return par[0]*CBall(x,par);}


double CBall(double *x, double *par) {

  double arg    = (x[0]-par[1])/fabs(par[2]);
  double lim    = par[4];

  if (arg < lim) return PartGau(x,par);
  return Tail(x,par);

}

double PartGau(double *x, double *par) {
  double arg    = (x[0]-par[1])/par[2];
  double Gauss  = TMath::Exp(-0.5*arg*arg);
  return Gauss;
}

double Tail(double *x, double *par) {

  double arg    = (x[0]-par[1])/fabs(par[2]);
  double power  = 1.0-par[3];
  double abso   = par[3]/fabs(par[4]);
  double A      = pow(abso,par[3])*TMath::Exp(-0.5*par[4]*par[4]);
  double B      = abso-fabs(par[4]);
  double tail   = B+arg;
  double resu   = A*pow(tail,-par[3]);
  return resu;
}








void setTree()
{
     tree->SetBranchAddress("nMatch",    &nMatch);
     tree->SetBranchAddress("mTrkId",    mTrkId);
     tree->SetBranchAddress("mCluId",    mCluId);
     tree->SetBranchAddress("mChi2",     mChi2);
     tree->SetBranchAddress("mChi2Pos",  mChi2Pos);
     tree->SetBranchAddress("mChi2Time", mChi2Time);

     tree->SetBranchAddress("nTrk",    &nTrk);
     tree->SetBranchAddress("trkX",    trkX);
     tree->SetBranchAddress("trkY",    trkY);
     tree->SetBranchAddress("trkZ",    trkZ);
     tree->SetBranchAddress("trkFFX",  trkFFX);
     tree->SetBranchAddress("trkFFY",  trkFFY);
     tree->SetBranchAddress("trkFFZ",  trkFFZ);
     tree->SetBranchAddress("trkt",    trkt);
     tree->SetBranchAddress("trke",    trke);
     tree->SetBranchAddress("trkpX",   trkpX);
     tree->SetBranchAddress("trkpY",   trkpY);
     tree->SetBranchAddress("trkpZ",   trkpZ);
     tree->SetBranchAddress("trknHit", trknHit);
     tree->SetBranchAddress("trkStat", trkStat);
     tree->SetBranchAddress("trkprob", trkprob);
     tree->SetBranchAddress("trkd0",   trkd0);
     tree->SetBranchAddress("trkz0",   trkz0);
     tree->SetBranchAddress("trkphi0", trkphi0);
     tree->SetBranchAddress("trkomega",trkomega);
     tree->SetBranchAddress("trkcdip", trkcdip);
     tree->SetBranchAddress("trkdlen", trkdlen);

     tree->SetBranchAddress("nCluster",  &nCluster);
     tree->SetBranchAddress("cluEnergy", cluEnergy);
     tree->SetBranchAddress("cluTime",   cluTime);
     tree->SetBranchAddress("cluCogX",   cluCogX);
     tree->SetBranchAddress("cluCogY",   cluCogY);
     tree->SetBranchAddress("cluCogZ",   cluCogZ);
     tree->SetBranchAddress("clu2Mom",   clu2Mom);
     tree->SetBranchAddress("cluE25",    cluE25);
     tree->SetBranchAddress("cluList",   &cluList);

     tree->SetBranchAddress("nCry",     &nCry);
     tree->SetBranchAddress("cryId",    cryId);
     tree->SetBranchAddress("crySecId", crySecId);
     tree->SetBranchAddress("cryPosX",  cryPosX);
     tree->SetBranchAddress("cryPosY",  cryPosY);
     tree->SetBranchAddress("cryPosZ",  cryPosZ);
     tree->SetBranchAddress("cryEdep",  cryEdep);
     tree->SetBranchAddress("cryTime",  cryTime);
     tree->SetBranchAddress("nVd",      &nVd);
     tree->SetBranchAddress("vdId",     vdId);
     tree->SetBranchAddress("vdPdgId",  vdPdgId);
     tree->SetBranchAddress("vdMom",    vdMom);
     tree->SetBranchAddress("vdMomX",   vdMomX);
     tree->SetBranchAddress("vdMomY",   vdMomY);
     tree->SetBranchAddress("vdMomZ",   vdMomZ);
     tree->SetBranchAddress("vdPosX",   vdPosX);
     tree->SetBranchAddress("vdPosY",   vdPosY);
     tree->SetBranchAddress("vdPosZ",   vdPosZ);
     tree->SetBranchAddress("vdTime",   vdTime);
}







void makeHist()
{
   
     gStyle->SetFrameBorderMode(0);
     gStyle->SetCanvasBorderMode(0);
     gStyle->SetPadBorderMode(0);
     gStyle->SetPadColor(0);
     gStyle->SetCanvasColor(0);
     gStyle->SetStatColor(0);
     gStyle->SetTitleFillColor(0);
     gStyle->SetOptStat(0);
     gStyle->SetHistLineWidth(3);
     gStyle->SetLineWidth(2);
     gStyle->SetPadLeftMargin(0.12);
     gStyle->SetPalette(1);
     gStyle->SetFrameLineWidth(2);
     gStyle->SetLegendBorderSize(0);
     gStyle->SetLegendFillColor(0);


     _hx     = new TH1F("hx","X residual",50,-30,30);
     _hy     = new TH1F("hy","Y residual",50,-30,30);
     _hr     = new TH1F("hr","R residual",200,0,2000);
     _hxy    = new TH1F("hxy","hxy",100,-30,30);
     _hxy2   = new TH2F("hxy2","",100,-50.,50.,100,-50.,50.);
     _hu2nd  = new TH2F("hu2nd","",100,0,200.,100,-50.,50.);
     _hv2nd  = new TH2F("hv2nd","",100,0,200.,100,-50.,50.);
     _huE25  = new TH2F("huE25","",100,0,1.,100,-50.,50.);
     _huDist = new TH2F("huDist","",100,0,500.,100,-50.,50.);
     _hucDip = new TH2F("hucDip","",100,0,1.,100,-50.,50.);
     _huRad  = new TH2F("huRad","",100,300,700.,100,-50.,50.);
     _prob2d  = new TH2F("prob2d","",50,0,1.,50,0.,1.);

     _hu0      = new TH1F("hu0","U residual init",150,-50,50);
     _hv0      = new TH1F("hv0","V residual init",150,-50,50);
     _hu       = new TH1F("hu","U residual",80,-30,50);
     _hv       = new TH1F("hv","V residual",100,-50,50);
     _huc      = new TH1F("huc","U residual",150,-50,50);
     _hvc      = new TH1F("hvc","V residual",150,-50,50);
     _ha       = new TH1F("ha","Depth",100,0,100);
     _hncr     = new TH1F("hncr","Ncr",20,0,20);
     _prob1    = new TH1F("prob1","Prob",50,0,1);
     _prob2    = new TH1F("prob2","Prob",50,0,1);
     _chi2prob = new TH1F("chi2prob","Prob",50,0,1);

 
      histHelper(_hx,"#DeltaX (mm)");
      histHelper(_hy,"#DeltaY (mm)");
      histHelper(_hr,"R (mm)");
      histHelper(_hu0,"#DeltaU (mm)");
      histHelper(_hv0,"#DeltaV (mm)");
      histHelper(_hu,"#DeltaU (mm)");
      histHelper(_hv,"#DeltaV (mm)");
      

     _hxy2->GetXaxis()->SetTitle("X residual (mm)");
     _hxy2->GetYaxis()->SetTitle("Y residual (mm)");
     _hxy2->GetXaxis()->SetTitleSize(0.045);
     _hxy2->GetYaxis()->SetTitleSize(0.045);
     

     for (int i=0;i<30;++i)_dedep[i] = buildDisk(Form("Edep%i_a",i),351,660,33.065); 

}

void histHelper(TH1F* hist,TString title)
{
     hist->GetXaxis()->SetTitle(title);
//     hist->GetYaxis()->SetTitle(Form("Entries / %f mm",hist->GetBinWidth(1)));
//     hist->GetYaxis()->SetTitleOffset(2);
     hist->GetYaxis()->SetLabelSize(0.045);
     hist->GetXaxis()->SetLabelSize(0.045);
     hist->GetXaxis()->SetTitleSize(0.055);
     hist->SetLineWidth(2);
     hist->SetLineColor(1);

}


TH2Poly* buildDisk(TString hname, double _radiusIn, double _radiusOut, double _cellSize)
{

       TH2Poly* _h2p = new TH2Poly(hname,"",-_radiusOut-_cellSize,_radiusOut+_cellSize,-_radiusOut-_cellSize,_radiusOut+_cellSize);
       _h2p->SetContour(50);

       int nRings    = int(1.5*_radiusOut/_cellSize)+1;
       int nHexagons = 1 + 3*nRings*(nRings-1);

       for (int i=0;i<nHexagons;++i) {
	   double x0(0),y0(0);
	   getXYPosition_disk(i,x0,y0);
	   x0 *= _cellSize;
	   y0 *= _cellSize;

	   if ( !isInsideDisk(x0,y0,_cellSize,_radiusIn,_radiusOut) ) continue;

	   double x[6]={0},y[6]={0};
           x[0]= x0-0.288675135*_cellSize; y[0]= y0-0.5*_cellSize;
           x[1]= x0+0.288675135*_cellSize; y[1]= y0-0.5*_cellSize;
           x[2]= x0+0.577350269*_cellSize; y[2]= y0;
           x[3]= x0+0.288675135*_cellSize; y[3]= y0+0.5*_cellSize;
           x[4]= x0-0.288675135*_cellSize; y[4]= y0+0.5*_cellSize;
           x[5]= x0-0.577350269*_cellSize; y[5]= y0;
           _h2p->AddBin(6, x, y);
       }

       return _h2p;
}

void getXYPosition_disk(int index, double&x, double& y) 
{        

   if (index==0) {x=0;y=0;return;} 

   int _step_l[6]={0,-1,-1,0,1,1};
   int _step_k[6]={1,1,0,-1,-1,0};

   int nRing = int(0.5+sqrt(0.25+(float(index)-1.0)/3.0));
   int nSeg  = (index -1 -3*nRing*(nRing-1))/nRing;
   int nPos  = (index -1 -3*nRing*(nRing-1))%nRing;

   int l =  nRing+(nPos+1)*_step_l[nSeg];
   int k = -nRing+(nPos+1)*_step_k[nSeg];

   for (int i=0;i<nSeg;++i) {
      l += _step_l[i]*nRing;
      k += _step_k[i]*nRing;
   }

   x = (l+k)*sqrt(3.0)/2.0;
   y = (l-k)/2.0;
   return;
}

bool isInsideDisk(double x, double y, double _cellSize, double _radiusIn, double _radiusOut) 
{          

   double apexX[7]={-0.2886751,+0.2886751,+0.5773502,+0.2886751,-0.2886751,-0.5773502,-0.2886751};
   double apexY[7]={-0.5,-0.5,0,0.5,0.5,0,-0.5};

   for (int ip=0;ip<6;++ip){  
       TVector2 p1(x+_cellSize*apexX[ip],   y+_cellSize*apexY[ip]);
       TVector2 p2(x+_cellSize*apexX[ip+1], y+_cellSize*apexY[ip+1]);
       if (calcDistToSide(p1,p2) < _radiusIn) return false;
       if (p1.Mod() > _radiusOut)             return false;      
   }
   return true;
}

double calcDistToSide(TVector2& a, TVector2& b) 
{
    double t = -3.0*(b-a)*a;
    if (t<0) return a.Mod();
    if (t>1) return b.Mod();
    return (a+t*(b-a)).Mod();
}





























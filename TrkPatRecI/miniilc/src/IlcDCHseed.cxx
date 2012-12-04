/**************************************************************************
 * Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 *                                                                        *
// Author: The ILC Off-line Project. 
 // Part of the code has been developed by Alice Off-line Project. 
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/




//-----------------------------------------------------------------
//           Implementation of the DCH seed class
//        This class is used by the IlcDCHtrackerMI class
//      Origin: Marian Ivanov, CERN, Marian.Ivanov@cern.ch
//-----------------------------------------------------------------
#include "TClonesArray.h"
#include "IlcDCHseed.h"
#include "TClass.h"

//ClassImp(IlcDCHseed)



IlcDCHseed::IlcDCHseed():
  IlcDCHtrack(),
//  fEsd(0x0),
  fClusterOwner(kFALSE),
  fPoints(0x0),
  fEPoints(0x0),
  fLayer(0),
  fErrorY2(1e10),
  fErrorZ2(1e10),
  fCurrentCluster(0x0),
  fCurrentClusterIndex(-1),
  fInDead(kFALSE),
  fIsSeeding(kFALSE),
  fSort(0),
  fBSigned(kFALSE),
  fSeedType(0),
  fSeed1(-1),
  fSeed2(-1),
  fMAngular(0),
  fCircular(0)
{
  //
  for (Int_t i=0;i<kMaxLayer;i++) for (Int_t j=0;j<kMaxInLayer;j++) SetClusterIndex(i,-3,j);
  for (Int_t i=0;i<kMaxLayer;i++) for (Int_t j=0;j<kMaxInLayer;j++) fClusterPointer[i][j]=0;
  for (Int_t i=0;i<kMaxLayer;i++) for (Int_t j=0;j<kMaxInLayer;j++) fClusterChi2[i][j]=-1;
  for (Int_t i=0;i<3;i++)   fKinkIndexes[i]=0;
  for (Int_t i=0;i<IlcPID::kSPECIES;i++)   fDCHr[i]=0.2;
  for (Int_t i=0;i<4;i++) {
    fDEDX[i] = 0.;
    fSDEDX[i] = 1e10;
    fNCDEDX[i] = 0;
  }
  for (Int_t i=0;i<12;i++) fOverlapLabels[i] = -1;
}

IlcDCHseed::IlcDCHseed(const IlcDCHseed &s, Bool_t clusterOwner):
  IlcDCHtrack(s),
  //  fEsd(0x0),
  fClusterOwner(clusterOwner),
  fPoints(0x0),
  fEPoints(0x0),
  fLayer(0),
  fErrorY2(1e10),
  fErrorZ2(1e10),
  fCurrentCluster(0x0),
  fCurrentClusterIndex(-1),
  fInDead(kFALSE),
  fIsSeeding(kFALSE),
  fSort(0),
  fBSigned(kFALSE),
  fSeedType(0),
  fSeed1(-1),
  fSeed2(-1),
  fMAngular(0),
  fCircular(0)
{
  //---------------------
  // dummy copy constructor
  //-------------------------
  for (Int_t i=0;i<kMaxLayer;i++) 
  for (Int_t j=0;j<kMaxInLayer;j++) {
    fClusterPointer[i][j]=0;
    if (fClusterOwner){
      if (s.fClusterPointer[i][j])
	fClusterPointer[i][j] = new IlcDCHcluster(*(s.fClusterPointer[i][j]));
    }else{
      fClusterPointer[i][j] = s.fClusterPointer[i][j];
    }
    fTrackPoints[i][j] = s.fTrackPoints[i][j];
  }
  for (Int_t i=0;i<kMaxLayer;i++) for (Int_t j=0;j<kMaxInLayer;j++) fClusterChi2[i][j]=s.fClusterChi2[i][j];
  for (Int_t i=0;i<kMaxLayer;i++) for (Int_t j=0;j<kMaxInLayer;j++)  fIndex[i][j] = s.fIndex[i][j];
  for (Int_t i=0;i<kMaxLayer;i++) fNInLayer[i] = s.fNInLayer[i];
  for (Int_t i=0;i<IlcPID::kSPECIES;i++)   fDCHr[i]=s.fDCHr[i];
  for (Int_t i=0;i<4;i++) {
    fDEDX[i] = s.fDEDX[i];
    fSDEDX[i] = s.fSDEDX[i];
    fNCDEDX[i] = s.fNCDEDX[i];
  }
  for (Int_t i=0;i<12;i++) fOverlapLabels[i] = s.fOverlapLabels[i];
}


IlcDCHseed::IlcDCHseed(const IlcDCHtrack &t):
  IlcDCHtrack(t),
  //  fEsd(0x0),
  fClusterOwner(kFALSE),
  fPoints(0x0),
  fEPoints(0x0),
  fLayer(0),
  fErrorY2(1e10),
  fErrorZ2(1e10),
  fCurrentCluster(0x0),
  fCurrentClusterIndex(-1),
  fInDead(kFALSE),
  fIsSeeding(kFALSE),
  fSort(0),
  fBSigned(kFALSE),
  fSeedType(0),
  fSeed1(-1),
  fSeed2(-1),
  fMAngular(0),
  fCircular(0)
{
  //
  // Constructor from IlcDCHtrack
  //
  fFirstPoint =0;
  for (Int_t i=0;i<5;i++)   fDCHr[i]=0.2;
  for (Int_t i=0;i<kMaxLayer;i++) 
  for (Int_t j=0;j<kMaxInLayer;j++){
    fClusterPointer[i][j] = 0;
    Int_t index = t.GetClusterIndex(i,j);
    if (index>=-1){ 
      SetClusterIndex(i,index,j);
    }
    else{
      SetClusterIndex(i,-3,j); 
    }    
  }
  for (Int_t i=0;i<kMaxLayer;i++) for (Int_t j=0;j<kMaxInLayer;j++) fClusterChi2[i][j]=-1;
  for (Int_t i=0;i<4;i++) {
    fDEDX[i] = 0.;
    fSDEDX[i] = 1e10;
    fNCDEDX[i] = 0;
  }
  for (Int_t i=0;i<12;i++) fOverlapLabels[i] = -1;
}

IlcDCHseed::IlcDCHseed(Double_t xr, Double_t alpha, const Double_t xx[5],
		       const Double_t cc[15], Int_t index):      
  IlcDCHtrack(xr, alpha, xx, cc, index),
  //  fEsd(0x0),
  fClusterOwner(kFALSE),
  fPoints(0x0),
  fEPoints(0x0),
  fLayer(0),
  fErrorY2(1e10),
  fErrorZ2(1e10),
  fCurrentCluster(0x0),
  fCurrentClusterIndex(-1),
  fInDead(kFALSE),
  fIsSeeding(kFALSE),
  fSort(0),
  fBSigned(kFALSE),
  fSeedType(0),
  fSeed1(-1),
  fSeed2(-1),
  fMAngular(0),
  fCircular(0)
{
  //
  // Constructor
  //
  fFirstPoint =0;
  for (Int_t i=0;i<kMaxLayer;i++) for (Int_t j=0;j<kMaxInLayer;j++) SetClusterIndex(i,-3,j);
  for (Int_t i=0;i<kMaxLayer;i++) for (Int_t j=0;j<kMaxInLayer;j++) fClusterPointer[i][j]=0;
  for (Int_t i=0;i<kMaxLayer;i++) for (Int_t j=0;j<kMaxInLayer;j++) fClusterChi2[i][j]=-1;
  for (Int_t i=0;i<5;i++)   fDCHr[i]=0.2;
  for (Int_t i=0;i<4;i++) {
    fDEDX[i] = 0.;
    fSDEDX[i] = 1e10;
    fNCDEDX[i] = 0;
  }
  for (Int_t i=0;i<12;i++) fOverlapLabels[i] = -1;
}

IlcDCHseed::~IlcDCHseed(){
  //
  // destructor
  if (fPoints) delete fPoints;
  fPoints =0;
  if (fEPoints) delete fEPoints;
  fEPoints = 0;
  if (fClusterOwner){
    for (Int_t icluster=0; icluster<kMaxLayer; icluster++)
      for (Int_t j=0;j<kMaxInLayer;j++){
	delete fClusterPointer[icluster][j];
      }
  }
}

IlcDCHTrackerPoint * IlcDCHseed::GetTrackPoint(Int_t i,Int_t j)
{
  //
  // 
  return &fTrackPoints[i][j];
}

void IlcDCHseed::RebuildSeed()
{
  //
  // rebuild seed to be ready for storing
  IlcDCHcluster cldummy;
  cldummy.SetQ(0);
  IlcDCHTrackPoint pdummy;
  pdummy.GetTPoint().SetShared(10);
  for (Int_t i=0;i<kMaxLayer;i++)for (Int_t j=0;j<kMaxInLayer;j++){
    IlcDCHcluster * cl0 = fClusterPointer[i][j];
    IlcDCHTrackPoint *trpoint = (IlcDCHTrackPoint*)fPoints->UncheckedAt(i*kMaxInLayer+j);     
    if (cl0){
      trpoint->GetTPoint() = *(GetTrackPoint(i,j));
      trpoint->GetCPoint() = *cl0;
      trpoint->GetCPoint().SetQ(TMath::Abs(cl0->GetQ()));
    }
    else{
      *trpoint = pdummy;
      trpoint->GetCPoint()= cldummy;
    }
    
  }

}


Double_t IlcDCHseed::GetDensityFirst(Int_t n)
{
  //
  //
  // return cluster for n layers bellow first point
  Int_t nfoundable = 1;
  Int_t nfound      = 1;
  for (Int_t i=fLastPoint-1;i>0&&nfoundable<n; i--)for(int j=0;j<kMaxInLayer;j++){
    Int_t index = GetClusterIndex(i,j);
    if (index!=-1&&j==0) nfoundable++;
    if (index>0) nfound++;
  }
  if (double(nfoundable)<double(n)*0.4) return 0;
  return Double_t(nfound)/Double_t(nfoundable);

}


void IlcDCHseed::GetClusterStatistic(Int_t first, Int_t last, Int_t &found, Int_t &foundable, Int_t &shared, Bool_t plus2)
{
  // get cluster stat.  on given region
  //
  found       = 0;
  foundable   = 0;
  shared      =0;
  for (Int_t i=first;i<=last; i++)for (Int_t j=0;j<kMaxInLayer;j++){
    Int_t index = GetClusterIndex(i,j);
    if (index!=-1&&j==0) foundable++;
    if (fClusterPointer[i][j]) {
      found++;
    }
    else 
      continue;

    if (fClusterPointer[i][j]->IsUsed(10)) {
      shared++;
      continue;
    }
    if (!plus2) continue; //take also neighborhoud
    //
    if ( (i>0) && fClusterPointer[i-1][0]){
      if (fClusterPointer[i-1][0]->IsUsed(10)) {
	shared++;
	continue;
      }
    }
    if ( fClusterPointer[i+1][0]){
      if (fClusterPointer[i+1][0]->IsUsed(10)) {
	shared++;
	continue;
      }
    }
    
  }

  IlcDebug(5,Form("from first layer=%i to last=%i, found =%i, foundable=%i,shared=%i,removal=%i",first,last,found,foundable,shared,GetRemoval()));
  //if (shared>found){
    //Error("IlcDCHseed::GetClusterStatistic","problem\n");
  //}
}





void IlcDCHseed::Reset(Bool_t all)
{
  //
  //
  SetNumberOfClusters(0);
  fNFoundable = 0;
  SetChi2(0);
  ResetCovariance(10.);
  /*
  if (fTrackPoints){
    for (Int_t i=0;i<8;i++){
      delete [] fTrackPoints[i];
    }
    delete fTrackPoints;
    fTrackPoints =0;
  }
  */
  
  for(Int_t i=0;i<kMaxLayer;i++) for (Int_t j=0;j<kMaxInLayer;j++) fClusterPointer[i][j]=0;

  if (all){   
    for (Int_t i=0;i<kMaxLayer;i++)for (Int_t j=0;j<kMaxInLayer;j++)SetClusterIndex(i,-3,j);
  }

}


void IlcDCHseed::Modify(Double_t factor)
{

  //------------------------------------------------------------------
  //This function makes a track forget its history :)  
  //------------------------------------------------------------------
  if (factor<=0) {
    ResetCovariance(10.);
    return;
  }
  ResetCovariance(factor);

  SetNumberOfClusters(0);
  fNFoundable =0;
  SetChi2(0);
  fRemoval = 0;
  //fFirstPoint = kMaxLayer;
  //fLastPoint  = 0;
}


//_________________________________________________________________________________________


Int_t IlcDCHseed::Compare(const TObject *o) const {
  //-----------------------------------------------------------------
  // This function compares tracks according to the sector - for given sector according z
  //-----------------------------------------------------------------
  IlcDCHseed *t=(IlcDCHseed*)o;

  if (fSort == 0){
    Double_t z2 = t->GetZ();
    Double_t z1 = GetZ();
    if (z2>z1) return 1;
    if (z2<z1) return -1;
    return 0;
  }
  else {
    Float_t f2 =1;
    f2 = 1-20*TMath::Sqrt(t->GetSigma1Pt2())/(TMath::Abs(t->Get1Pt())+0.0066);
    if (t->fBConstrain) f2=1.2;

    Float_t f1 =1;
    f1 = 1-20*TMath::Sqrt(GetSigma1Pt2())/(TMath::Abs(Get1Pt())+0.0066);

    if (fBConstrain)   f1=1.2;
 
    if (t->GetNumberOfClusters()*f2 <GetNumberOfClusters()*f1) return -1;
    else return +1;
  }
}


//_____________________________________________________________________________
Double_t IlcDCHseed::GetPredictedChi2(const IlcCluster *c) const 
{
  //-----------------------------------------------------------------
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  Error("IlcDCHseed::Update","Plane chi2 doesn't have sence forr DCH");
  Double_t p[2]={c->GetY(), c->GetZ()};
  Double_t cov[3]={fErrorY2, 0., fErrorZ2};
  return IlcExternalTrackParam::GetPredictedChi2(p,cov);
}

#include <iostream>
using namespace std; 
//_____________________________________________________________________________
Double_t IlcDCHseed::GetPredictedChi2(const IlcCluster *c2,
				      double rwire,double anglewire,double tanstereo,bool withZ) const 
{
  //-----------------------------------------------------------------
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  const IlcDCHcluster *c=dynamic_cast<const IlcDCHcluster*>(c2);
  if(!c) return kFALSE;

  Double_t p[2]={c->GetImP(), c->GetZLocal()};
  Double_t cov[3]={c->GetSigmaImP2(), 0., c->GetSigmaZ2()};
  Double_t chi2= IlcExternalTrackParam::GetPredictedChi2(rwire,anglewire,tanstereo,GetBz(),
							 p,cov,withZ,kTRUE);
  if(IlcDebugLevelClass()>10){
    double dd[2],ddz[10];
    GetDistance2Wire(rwire,anglewire,tanstereo,GetBz(),ddz,dd,kTRUE);
    cout<<"IlcDCHseed::Chi2: chi2="<<chi2<<" clErY="<<sqrt(fErrorY2)<<" "<<sqrt(fErrorZ2)<<endl;
    cout<<"dist to wire "<<dd[0]<<" "<<dd[1]<<" mes "<<p[0]<<" "<<p[1]
	<<" phi tr="<<Phi_0_2pi(GetPhi()-TMath::ATan(GetZ()/rwire*tanstereo))
	<<" wire0="<<Phi_0_2pi(anglewire)
	<<" mes="<<Phi_0_2pi(TMath::ATan2(c->GetY(),c->GetX())-TMath::ATan(c->GetZ()/rwire*tanstereo))
	<<" trackid="<<c->GetLabel(0)<<" "<<c->GetLabel(1)<<endl;  
  }
  return chi2;
}



//_____________________________________________________________________________
Bool_t IlcDCHseed::Update(const IlcCluster */*c*/, Double_t /*chisq*/, Int_t /*index*/)
{
  //-----------------------------------------------------------------
  // This function associates a cluster with this track.
  //-----------------------------------------------------------------
  Error("IlcDCHseed::Update","Plane update doesn't have sence forr DCH");
  return kFALSE;
}

Bool_t IlcDCHseed::Update(const IlcCluster* c2, Double_t chi2, Int_t /*i*/,
			  double rwire,double anglewire,double tanstereo,bool withZ,int signImP,bool calcChi2){
  const IlcDCHcluster *c=dynamic_cast<const IlcDCHcluster*>(c2);
  if(!c) return kFALSE;
  
  Double_t p[2]={(signImP!=0?signImP:1.)*c->GetImP(), c->GetZLocal()};
  Double_t cov[3]={fErrorY2, 0., fErrorZ2};

  if(signImP==0)
    if(fabs(c->GetImP())<3*sqrt(c->GetSigmaImP2())) {
      cov[0]=fErrorY2+c->GetImP()*c->GetImP();
      p[0]=0;
    }

  if(IlcDebugLevelClass()>10){
    double dd[2],ddz[10];
    GetDistance2Wire(rwire,anglewire,tanstereo,GetBz(),ddz,dd,kTRUE);
    Double_t p2[2]={c->GetImP(), c->GetZLocal()};
    Double_t cov2[3]={c->GetSigmaImP2(), 0., c->GetSigmaZ2()};
    double chi2=IlcExternalTrackParam::GetPredictedChi2(dd,ddz,p2,cov2,withZ);
    cout<<"IlcDCHseed::Update: chi2="<<chi2<<" clErY="<<sqrt(cov[0])<<" "<<sqrt(fErrorZ2)
	<<" tr "<<sqrt(GetSigmaY2())<<" "<<sqrt(GetSigmaZ2())<<" totChi2 "<<GetChi2()<<endl;
    Print("param");
    double phit=Phi_0_2pi(GetPhi()-TMath::ATan(GetZ()/rwire*tanstereo));
    double phiw=Phi_0_2pi(anglewire);
    double phim=Phi_0_2pi(TMath::ATan2(c->GetY(),c->GetX())-TMath::ATan(c->GetZ()/rwire*tanstereo));
    cout<<"ncls="<<GetNumberOfClusters()<<" lay="<<GetLayer()
	<<" dist to wire "<<dd[0]<<" "<<dd[1]<<" mes "<<p[0]<<" "<<p[1]
	<<" lr "<<(Phi_mpi_pi(phit-phiw)*Phi_mpi_pi(phim-phiw)>=0)
	<<" phi tr="<<phit<<" wire0="<<phiw<<" mes="<<phim
	<<" trackid="<<c->GetLabel(0)<<" "<<c->GetLabel(1)<<endl;
  }

  if(calcChi2){
    double deriv[10],rztrack[2];
    GetDistance2Wire(rwire,anglewire,tanstereo,GetBz(),deriv,rztrack,kTRUE);
    Double_t p2[2]={c->GetImP(), c->GetZLocal()};
    Double_t cov2[3]={c->GetSigmaImP2(), 0., c->GetSigmaZ2()};
    chi2=IlcExternalTrackParam::GetPredictedChi2(rztrack,deriv,p2,cov2,withZ);
    //    std::cout<<chi2<<" "<<signImP<<" "<<p[0]<<" "<<p[1]<<" "<<p2[0]<<" "<<p2[1]<<endl;
    if(!UpdateWithWire(rztrack,deriv,p,cov,withZ)) return kFALSE;
  }else{
    if (!IlcExternalTrackParam::UpdateWithWire(rwire,anglewire,tanstereo,GetBz(),
					       p,cov,withZ,kTRUE)) return kFALSE;
  }

  Int_t n=GetNumberOfClusters();
  //  fIndex[n]=index;
  SetNumberOfClusters(n+1);
  SetChi2(GetChi2()+chi2);

  return kTRUE;
}


//_____________________________________________________________________________
Float_t IlcDCHseed::CookdEdx(Double_t low, Double_t up,Int_t i1, Int_t i2, Bool_t onlyused) {
  //-----------------------------------------------------------------
  // This funtion calculates dE/dX within the "low" and "up" cuts.
  //-----------------------------------------------------------------

  Float_t amp[kMaxLayer];
  Float_t angular[kMaxLayer];
  Float_t weight[kMaxLayer];
  Int_t index[kMaxLayer];
  //Int_t nc = 0;
  //  TClonesArray & arr = *fPoints; 
  Float_t meanlog = 100.;
  
  Float_t mean[4]  = {0,0,0,0};
  Float_t sigma[4] = {1000,1000,1000,1000};
  Int_t nc[4]      = {0,0,0,0};
  Float_t norm[4]    = {1000,1000,1000,1000};
  //
  //
  fNShared =0;

  for (Int_t of =0; of<4; of++){    
    for (Int_t i=of+i1;i<i2;i+=4)
      for(int j=0;j<kMaxInLayer;j++)
      {
	Int_t index = fIndex[i][j];
	if (index<0||index&ClusterIndex::NotAccept) continue;

	IlcDCHTrackerPoint * point = GetTrackPoint(i,j);

	if (point==0) continue;
	IlcDCHcluster * cl = fClusterPointer[i][j];
	if (cl==0) continue;	
	if (onlyused && (!cl->IsUsed(10))) continue;
	if (cl->IsUsed(11)) {
	  fNShared++;
	  continue;
	}

	Float_t angley = point->GetAngleY();
	Float_t anglez = point->GetAngleZ();

	Float_t rsigmay2 =  point->GetSigmaY();
	Float_t rsigmaz2 =  point->GetSigmaZ();

	Float_t rsigma = TMath::Sqrt(rsigmay2*rsigmaz2);

	Float_t ampc   = 0;     // normilczation to the number of electrons
	
	ampc = 1.*cl->GetQ();

	ampc *= 2.0;     // put mean value to channel 50
	//ampc *= 0.58;     // put mean value to channel 50
	Float_t w      =  1.;
	//	if (type>0)  w =  1./(type/2.-0.5); 
	//	Float_t z = TMath::Abs(cl->GetZ());
	if (i<64) {
	  ampc /= 0.6;
	  //ampc /= (1+0.0008*z);
	} else
	  if (i>128){
	    ampc /=1.5;
	    //ampc /= (1+0.0008*z);
	  }else{
	    //ampc /= (1+0.0008*z);
	  }
	
	if (rsigma>1.5) ampc/=1.3;  // if big backround
	amp[nc[of]]        = ampc;
	angular[nc[of]]    = TMath::Sqrt(1.+angley*angley+anglez*anglez);
	weight[nc[of]]     = w;
	nc[of]++;
      }
    
    TMath::Sort(nc[of],amp,index,kFALSE);
    Float_t sumamp=0;
    Float_t sumamp2=0;
    Float_t sumw=0;
    //meanlog = amp[index[Int_t(nc[of]*0.33)]];
    meanlog = 50;
    for (Int_t i=int(nc[of]*low+0.5);i<int(nc[of]*up+0.5);i++){
      Float_t ampl      = amp[index[i]]/angular[index[i]];
      ampl              = meanlog*TMath::Log(1.+ampl/meanlog);
      //
      sumw    += weight[index[i]]; 
      sumamp  += weight[index[i]]*ampl;
      sumamp2 += weight[index[i]]*ampl*ampl;
      norm[of]    += angular[index[i]]*weight[index[i]];
    }
    if (sumw<1){ 
      SetdEdx(0);  
    }
    else {
      norm[of] /= sumw;
      mean[of]  = sumamp/sumw;
      sigma[of] = sumamp2/sumw-mean[of]*mean[of];
      if (sigma[of]>0.1) 
	sigma[of] = TMath::Sqrt(sigma[of]);
      else
	sigma[of] = 1000;
      
    mean[of] = (TMath::Exp(mean[of]/meanlog)-1)*meanlog;
    //mean  *=(1-0.02*(sigma/(mean*0.17)-1.));
    //mean *=(1-0.1*(norm-1.));
    }
  }

  Float_t dedx =0;
  fSdEdx =0;
  fMAngular =0;
  //  mean[0]*= (1-0.05*(sigma[0]/(0.01+mean[1]*0.18)-1));
  //  mean[1]*= (1-0.05*(sigma[1]/(0.01+mean[0]*0.18)-1));

  
  //  dedx = (mean[0]* TMath::Sqrt((1.+nc[0]))+ mean[1]* TMath::Sqrt((1.+nc[1])) )/ 
  //  (  TMath::Sqrt((1.+nc[0]))+TMath::Sqrt((1.+nc[1])));

  Int_t norm2 = 0;
  Int_t norm3 = 0;
  for (Int_t i =0;i<4;i++){
    if (nc[i]>2&&nc[i]<1000){
      dedx      += mean[i] *nc[i];
      fSdEdx    += sigma[i]*(nc[i]-2);
      fMAngular += norm[i] *nc[i];    
      norm2     += nc[i];
      norm3     += nc[i]-2;
    }
    fDEDX[i]  = mean[i];             
    fSDEDX[i] = sigma[i];            
    fNCDEDX[i]= nc[i]; 
  }

  if (norm3>0){
    dedx   /=norm2;
    fSdEdx /=norm3;
    fMAngular/=norm2;
  }
  else{
    SetdEdx(0);
    return 0;
  }

  
  SetdEdx(dedx);
  return dedx;
}
Double_t IlcDCHseed::Bethe(Double_t bg){
  //
  // This is the Bethe-Bloch function normilcsed to 1 at the minimum
  //
  Double_t bg2=bg*bg;
  Double_t bethe;
  if (bg<3.5e1) 
    bethe=(1.+ bg2)/bg2*(log(5940*bg2) - bg2/(1.+ bg2));
  else // Density effect ( approximately :) 
    bethe=1.15*(1.+ bg2)/bg2*(log(3.5*5940*bg) - bg2/(1.+ bg2));
  return bethe/11.091;
}

void IlcDCHseed::CookPID()
{
  //
  // cook PID information according dEdx
  //
  Double_t fRange = 10.;
  Double_t fRes   = 0.1;
  Double_t fMIP   = 47.;
  //
  Int_t ns=IlcPID::kSPECIES;
  Double_t sumr =0;
  for (Int_t j=0; j<ns; j++) {
    Double_t mass=IlcPID::ParticleMass(j);
    Double_t mom=GetP();
    Double_t dedx=fdEdx/fMIP;
    Double_t bethe=Bethe(mom/mass); 
    Double_t sigma=fRes*bethe;
    if (sigma>0.001){
      if (TMath::Abs(dedx-bethe) > fRange*sigma) {
	fDCHr[j]=TMath::Exp(-0.5*fRange*fRange)/sigma;
	sumr+=fDCHr[j];
	continue;
      }
      fDCHr[j]=TMath::Exp(-0.5*(dedx-bethe)*(dedx-bethe)/(sigma*sigma))/sigma;
      sumr+=fDCHr[j];
    }
    else{
      fDCHr[j]=1.;
      sumr+=fDCHr[j];
    }
  }
  for (Int_t j=0; j<ns; j++) {
    fDCHr[j]/=sumr;           //normilcze
  }
}


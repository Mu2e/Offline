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

/* $Id: IlcDCHtrack.cxx,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */

//-----------------------------------------------------------------
//           Implementation of the DCH track class
//        This class is used by the IlcDCHtracker class
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include <Riostream.h>

#include "IlcDCHtrack.h"
#include "IlcCluster.h"
#include "IlcDCHcluster.h"
#include "IlcTracker.h"
//#include "IlcESDtrack.h"

//ClassImp(IlcDCHtrack)

//_________________________________________________________________________
IlcDCHtrack::IlcDCHtrack(): 
  IlcKalmanTrack(),
  fdEdx(0),
  fSdEdx(1e10),
  fNFoundable(0),
  fBConstrain(kFALSE),
  fLastPoint(-1),
  fFirstPoint(-1),
  fRemoval(0),
  fTrackType(0),
  fLab2(-1),
  fNShared(0),
  fReference()
{
  //-------------------------------------------------
  // default constructor
  //-------------------------------------------------
  for (Int_t i=0; i<kMaxLayer;i++) for (Int_t j=0; j<kMaxInLayer;j++) fIndex[i][j]=-2;
  for (Int_t i=0; i<kMaxLayer;i++) fNInLayer[i]=0;
  for (Int_t i=0; i<4;i++) fPoints[i]=0.;
  for (Int_t i=0; i<12;i++) fKinkPoint[i]=0.;
  for (Int_t i=0; i<3;i++) fKinkIndexes[i]=0;
  for (Int_t i=0; i<3;i++) fV0Indexes[i]=0;
}

//_________________________________________________________________________



IlcDCHtrack::IlcDCHtrack(Double_t x, Double_t alpha, const Double_t p[5],
		         const Double_t cov[15], Int_t index) :
  IlcKalmanTrack(),
  fdEdx(0),
  fSdEdx(1e10),
  fNFoundable(0),
  fBConstrain(kFALSE),
  fLastPoint(0),
  fFirstPoint(0),
  fRemoval(0),
  fTrackType(0),
  fLab2(0),
  fNShared(0),
  fReference()
{
  //-----------------------------------------------------------------
  // This is the main track constructor.
  //-----------------------------------------------------------------
  Double_t cnv=1./(IlcTracker::GetBz()*kB2C);

  Double_t pp[5]={
    p[0],
    p[1],
    x*p[4] - p[2],
    p[3],
    p[4]*cnv
  };

  Double_t c22 = x*x*cov[14] - 2*x*cov[12] + cov[5];
  Double_t c32 = x*cov[13] - cov[8];
  Double_t c20 = x*cov[10] - cov[3], 
           c21 = x*cov[11] - cov[4], c42 = x*cov[14] - cov[12];

  Double_t cc[15]={
    cov[0 ],
    cov[1 ],     cov[2 ],
    c20,         c21,         c22,
    cov[6 ],     cov[7 ],     c32,     cov[9 ],
    cov[10]*cnv, cov[11]*cnv, c42*cnv, cov[13]*cnv, cov[14]*cnv*cnv
  };

//   Double_t p0=TMath::Sign(1/kMostProbablePt,pp[4]);
//   Double_t w0=cc[14]/(cc[14] + p0*p0), w1=p0*p0/(cc[14] + p0*p0);
//   pp[4] = w0*p0 + w1*pp[4]; 
//   cc[10]*=w1; cc[11]*=w1; cc[12]*=w1; cc[13]*=w1; cc[14]*=w1;

  Set(x,alpha,pp,cc);

  SetNumberOfClusters(0);
  
  fIndex[0][0]=index;
  for (Int_t i=1; i<kMaxLayer;i++)for (Int_t j=0; j<kMaxInLayer;j++) fIndex[i][j]=-2;
  for (Int_t i=1; i<kMaxLayer;i++) fNInLayer[i]=0;
  for (Int_t i=0; i<4;i++) fPoints[i]=0.;
  for (Int_t i=0; i<12;i++) fKinkPoint[i]=0.;
  for (Int_t i=0; i<3;i++) fKinkIndexes[i]=0;
  for (Int_t i=0; i<3;i++) fV0Indexes[i]=0;
}
/*
//_____________________________________________________________________________
IlcDCHtrack::IlcDCHtrack(const IlcESDtrack& t) :
  IlcKalmanTrack(),
  fdEdx(t.GetTPCsignal()),
  fSdEdx(1e10),
  fNFoundable(0),
  fBConstrain(kFALSE),
  fLastPoint(0),
  fFirstPoint(0),
  fRemoval(0),
  fTrackType(0),
  fLab2(0),
  fNShared(0),
  fReference()
{
  //-----------------------------------------------------------------
  // Conversion IlcESDtrack -> IlcDCHtrack.
  //-----------------------------------------------------------------
  for (Int_t i=0; i<kMaxLayer;i++) fIndex[i]=-2;
  SetNumberOfClusters(t.GetTPCclusters(fIndex,kMaxLayer));
  SetLabel(t.GetLabel());
  SetMass(t.GetMass());
  for (Int_t i=0; i<4;i++) fPoints[i]=0.;
  for (Int_t i=0; i<12;i++) fKinkPoint[i]=0.;
  for (Int_t i=0; i<3;i++) fKinkIndexes[i]=0;
  for (Int_t i=0; i<3;i++) fV0Indexes[i]=0;

  //  Set(t.GetX(),t.GetAlpha(),t.GetParameter(),t.GetCovariance());
  Set(&t);

  if ((t.GetStatus()&IlcESDtrack::kTIME) == 0) return;
  StartTimeIntegral();
  Double_t times[10]; t.GetIntegratedTimes(times); SetIntegratedTimes(times);
  SetIntegratedLength(t.GetIntegratedLength());
}
*/
//_____________________________________________________________________________
IlcDCHtrack::IlcDCHtrack(const IlcDCHtrack& t) :
  IlcKalmanTrack(t),
  fdEdx(t.fdEdx),
  fSdEdx(t.fSdEdx),
  fNFoundable(t.fNFoundable),
  fBConstrain(t.fBConstrain),
  fLastPoint(t.fLastPoint),
  fFirstPoint(t.fFirstPoint),
  fRemoval(t.fRemoval),
  fTrackType(t.fTrackType),
  fLab2(t.fLab2),
  fNShared(t.fNShared),
  fReference(t.fReference)

{
  //-----------------------------------------------------------------
  // This is a track copy constructor.
  //-----------------------------------------------------------------
  //Set(t.GetX(),t.GetAlpha(),t.GetParameter(),t.GetCovariance());
  Set(&t);

  for (Int_t i=0; i<kMaxLayer; i++)for (Int_t j=0; j<kMaxInLayer; j++) fIndex[i][j]=t.fIndex[i][j];
  for (Int_t i=0; i<kMaxLayer; i++) fNInLayer[i]=t.fNInLayer[i];
  for (Int_t i=0; i<4;i++) fPoints[i]=t.fPoints[i];
  for (Int_t i=0; i<12;i++) fKinkPoint[i]=t.fKinkPoint[i];
  for (Int_t i=0; i<3;i++) fKinkIndexes[i]=t.fKinkIndexes[i];
  for (Int_t i=0; i<3;i++) fV0Indexes[i]=t.fV0Indexes[i];
}

//_____________________________________________________________________________
Int_t IlcDCHtrack::Compare(const TObject *o) const {
  //-----------------------------------------------------------------
  // This function compares tracks according to the their curvature
  //-----------------------------------------------------------------
  IlcDCHtrack *t=(IlcDCHtrack*)o;
  Double_t co=t->GetSigmaY2()*t->GetSigmaZ2();
  Double_t c =GetSigmaY2()*GetSigmaZ2();
  if (c>co) return 1;
  else if (c<co) return -1;
  return 0;
}

//_____________________________________________________________________________
Double_t IlcDCHtrack::GetPredictedChi2(const IlcCluster *c) const 
{
  //-----------------------------------------------------------------
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  Error("IlcDCHseed::Update","Plane chi2 doesn't have sence forr DCH");
  Double_t p[2]={c->GetY(), c->GetZ()};
  Double_t cov[3]={c->GetSigmaY2(), 0., c->GetSigmaZ2()};
  return IlcExternalTrackParam::GetPredictedChi2(p,cov);
}

//_____________________________________________________________________________
Double_t IlcDCHtrack::GetPredictedChi2(const IlcCluster *c2,
				      double rwire,double anglewire,double tanstereo,bool withZ) const 
{
  //-----------------------------------------------------------------
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  const IlcDCHcluster *c=dynamic_cast<const IlcDCHcluster*>(c2);
  if(!c) return kFALSE;

  Double_t p[2]={c->GetImP(), c->GetZ()};
  Double_t cov[3]={c->GetSigmaImP2(), 0., c->GetSigmaZ2()};
  return IlcExternalTrackParam::GetPredictedChi2(rwire,anglewire,tanstereo,GetBz(),
						 p,cov,withZ);
}


//_____________________________________________________________________________
Bool_t IlcDCHtrack::PropagateTo(Double_t xk,Double_t rho,Double_t x0) {
  //-----------------------------------------------------------------
  //  This function propagates a track to a reference plane x=xk.
  //  rho - density of the crossed matrial (g/cm^3)
  //  x0  - radiation length of the crossed material (g/cm^2) 
  //-----------------------------------------------------------------
  Double_t oldX=GetX(), oldY=GetY(), oldZ=GetZ();
  Double_t oldR=TMath::Hypot(oldX,oldY);

  Double_t bz=GetBz();
  if (!IlcExternalTrackParam::PropagateToR(xk,bz)) return kFALSE;

  Double_t d = TMath::Sqrt((GetX()-oldX)*(GetX()-oldX) + 
                           (GetY()-oldY)*(GetY()-oldY) + 
                           (GetZ()-oldZ)*(GetZ()-oldZ));
  if (IsStartedTimeIntegral() && GetX()>oldR) AddTimeStep(d);
  
  if ((xk-oldR)*GetDir()>0) d = -d;
  if (!IlcExternalTrackParam::CorrectForMaterial(d*rho/x0,x0,GetMass())) 
    return kFALSE;
  
  return kTRUE;
}
/*
//_____________________________________________________________________________
Bool_t 
IlcDCHtrack::PropagateToVertex(const IlcESDVertex *v,Double_t d,Double_t x0) 
{
  //-----------------------------------------------------------------
  // This function propagates tracks to the vertex.
  //-----------------------------------------------------------------
  Double_t bz=GetBz();
  if (PropagateToDCA(v,bz,kVeryBig))
   if (IlcExternalTrackParam::CorrectForMaterial(d,x0,GetMass())) return kTRUE;
  return kFALSE;
}
*/
//_____________________________________________________________________________
Bool_t IlcDCHtrack::Update(const IlcCluster */*c*/, Double_t /*chisq*/, Int_t /*index*/)
{
  //-----------------------------------------------------------------
  // This function associates a cluster with this track.
  //-----------------------------------------------------------------
  Error("IlcDCHseed::Update","Plane update doesn't have sence forr DCH");
  return kFALSE;
}

Bool_t IlcDCHtrack::Update(const IlcCluster* c2, Double_t chi2, Int_t index,
			   double rwire,double anglewire,double tanstereo,bool withZ){
  Error("IlcDCHtrack::Update","not correct forr use");
  const IlcDCHcluster *c=dynamic_cast<const IlcDCHcluster*>(c2);
  if(!c) return kFALSE;
  
  Double_t p[2]={c->GetImP(), c->GetZ()};
  Double_t cov[3]={c->GetSigmaImP2(), 0., c->GetSigmaZ2()};

  if (!IlcExternalTrackParam::UpdateWithWire(rwire,anglewire,tanstereo,GetBz(),
					     p,cov,withZ)) return kFALSE;

  Int_t n=GetNumberOfClusters();
  fIndex[n][0]=index;
  SetNumberOfClusters(n+1);
  SetChi2(GetChi2()+chi2);

  return kTRUE;
}

////////////////////////////////////////////////////////////////////////
// MI ADDITION

Float_t IlcDCHtrack::Density(Int_t layer0, Int_t layer1)
{
  //
  // calculate cluster density
  Int_t good  = 0;
  Int_t found = 0;
  //if (layer0<fFirstPoint) layer0 = fFirstPoint;
  if (layer1>fLastPoint) layer1 = fLastPoint;

  
  for (Int_t i=layer0;i<=layer1;i++){ 
  for (Int_t j=0;j<=kMaxInLayer;j++){ 
    //    Int_t index = fClusterIndex[i];
    Int_t index = fIndex[i][j];
    if (index!=-1)  good++;
    if (index>=0)    found++;
  }}
  Float_t density=0;
  if (good>0) density = Float_t(found)/Float_t(good);
  return density;
}


Float_t IlcDCHtrack::Density2(Int_t layer0, Int_t layer1)
{
  //
  // calculate cluster density
  Int_t good  = 0;
  Int_t found = 0;
  //  
  for (Int_t i=layer0;i<=layer1;i++){     
  for (Int_t j=0;j<=kMaxInLayer;j++){ 
    Int_t index = fIndex[i][j];
    if (index!=-1)  good++;
    if (index>=0)    found++;
  }}
  Float_t density=0;
  if (good>0) density = Float_t(found)/Float_t(good);
  return density;
}

void  IlcDCHtrack::UpdatePoints()
{
  //--------------------------------------------------
  //calculates first ,amx dens and last points
  //--------------------------------------------------
  Float_t density[kMaxLayer];
  for (Int_t i=0;i<kMaxLayer;i++) density[i]=-1.;
  fPoints[0]= kMaxLayer;
  fPoints[1] = -1;
  //
  Int_t ngood=0;
  Int_t undeff=0;
  Int_t nall =0;
  Int_t range=20;
  for (Int_t i=0;i<kMaxLayer;i++){
    Int_t last = i-range;
    if (nall<range) nall++;
    if (last>=0){
      if (fIndex[last][0]>=0&& (fIndex[last][0]&ClusterIndex::NotAccept)==0) ngood--;
      if (fIndex[last][0]==-1) undeff--;
    }
    if (fIndex[i][0]>=0&& (fIndex[i][0]&ClusterIndex::NotAccept)==0)   ngood++;
    if (fIndex[i][0]==-1) undeff++;
    if (nall==range &&undeff<range/2) {
      density[i-range/2] = Float_t(ngood)/Float_t(nall-undeff);
//       cout<<i-range/2<<" "<<density[i-range/2]<<" "<<ngood<<" "<<nall<<" "<<undeff<<endl;
    }
  }
  Float_t maxdens=0;
  Int_t indexmax =0;
  for (Int_t i=0;i<kMaxLayer;i++){
    if (density[i]<0) continue;
    if (density[i]>maxdens){
      maxdens=density[i];
      indexmax=i;
    }
  }
  //
  //max dens point
  fPoints[3] = maxdens;
  fPoints[1] = indexmax;
  //
  // last point
  for (Int_t i=indexmax;i<kMaxLayer;i++){
    if (density[i]<0) continue;
    if (density[i]<maxdens/2.5) {
      break;
    }
    fPoints[2]=i;
  }
  //
  // first point
  for (Int_t i=indexmax;i>0;i--){
    if (density[i]<0) continue;
    if (density[i]<maxdens/2.) {
      break;
    }
    fPoints[0]=i;
  }
  //
}

Double_t IlcDCHtrack::GetBz() const {
  //
  // returns Bz component of the magnetic field (kG)
  //
  if (IlcTracker::UniformField()) return IlcTracker::GetBz();
  Double_t r[3]; GetXYZ(r);
  return IlcTracker::GetBz(r);
}

Double_t IlcDCHtrack::GetYat(Double_t xk) const {
//-----------------------------------------------------------------
// This function calculates the Y-coordinate of a track at the plane x=xk.
//-----------------------------------------------------------------
  double y;
  if(!GetYAt(xk,GetBz(),y)) return 0;
  return y;
}



Int_t  IlcDCHtrack::GetProlongation(Double_t xk, Double_t &y, Double_t & z) const
{
  //-----------------------------------------------------------------
  // This function find proloncation of a track to a reference plane x=xk.
  // doesn't change internal state of the track
  //-----------------------------------------------------------------
  double xyz[3];
  bool localCoor=kTRUE;
  if(!GetXYZAt(xk,GetBz(),xyz,localCoor)) return 0;
  y=xyz[1];
  z=xyz[2];
  return 1;

}



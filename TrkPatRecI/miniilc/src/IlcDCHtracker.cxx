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


//-------------------------------------------------------
//          Implementation of the DCH tracker
//
//   Origin: Marian Ivanov   Marian.Ivanov@cern.ch
// 
//  IlcDCH parallel tracker
//-------------------------------------------------------


/* $Id: IlcDCHtracker.cxx,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */

#include "Riostream.h"
#include <TClonesArray.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TTree.h>
//#include "IlcLog.h"

#include "IlcComplexCluster.h"
//#include "IlcESD.h"
//#include "IlcKink.h"
//#include "IlcV0.h"
#include "IlcHelix.h"
//#include "IlcRunLoader.h"
#include "IlcDCHReconstructor.h"
#include "IlcDCHpolyTrack.h"
#include "IlcDCHreco.h"
#include "IlcDCHseed.h" 
#include "IlcDCHtracker.h"
#include "TStopwatch.h"
#include "IlcDCHReconstructor.h"
#include "IlcPID.h"
#include "TTreeStream.h"
//#include "IlcAlignObj.h"
#include "IlcTrackPointArray.h"
#include "IlcDCHrecoParam.h"
#include "IlcDCHParam.h"
#include "IlcDCHwireposition.h"

#include <map>
//
//ClassImp(IlcDCHtracker)


class IlcDCHFastMath {
public:
  IlcDCHFastMath();  
  static Double_t FastAsin(Double_t x);   
 private: 
  static Double_t fgFastAsin[20000];  //lookup table for fast asin computation
};

Double_t IlcDCHFastMath::fgFastAsin[20000];
IlcDCHFastMath gIlcDCHFastMath; // needed to fill the LUT

IlcDCHFastMath::IlcDCHFastMath(){
  //
  // initialized lookup table;
  for (Int_t i=0;i<10000;i++){
    fgFastAsin[2*i] = TMath::ASin(i/10000.);
    fgFastAsin[2*i+1] = (TMath::ASin((i+1)/10000.)-fgFastAsin[2*i]);
  }
}

Double_t IlcDCHFastMath::FastAsin(Double_t x){
  //
  // return asin using lookup table
  if (x>0){
    Int_t index = int(x*10000);
    return fgFastAsin[2*index]+(x*10000.-index)*fgFastAsin[2*index+1];
  }
  x*=-1;
  Int_t index = int(x*10000);
  return -(fgFastAsin[2*index]+(x*10000.-index)*fgFastAsin[2*index+1]);
}




Int_t IlcDCHtracker::UpdateTrack(IlcDCHseed * track, Int_t accept){
  //
  //update track information using current cluster - track->fCurrentCluster


  IlcDCHcluster* c =track->GetCurrentCluster();
  if (accept>0) track->SetCurrentCluster(c,track->GetCurrentClusterIndex() | ClusterIndex::NotAccept);  //sign not accepted clusters

  Int_t dchindex = track->GetCurrentClusterIndex();

  int layer=(dchindex&ClusterIndex::Layer)>>ClusterIndex::LayerOffset;
  track->SetLayer(layer);

  int ninlayer=track->GetClusterPointerI(layer,c);
  track->SetClusterIndex(track->GetLayer(), dchindex,ninlayer);  
  track->SetClusterPointer(track->GetLayer(),c,ninlayer);  

  if (track->GetFirstPoint()>track->GetLayer()) 
    track->SetFirstPoint(track->GetLayer());
  if (track->GetLastPoint()<track->GetLayer()) 
    track->SetLastPoint(track->GetLayer());
  

  //

  Double_t angle2 = track->GetSnp()*track->GetSnp();
  //
  //SET NEW Track Point
  //
  if (angle2<1) //PH sometimes angle2 is very big. To be investigated...
  {
    angle2 = TMath::Sqrt(angle2/(1-angle2)); 
    IlcDCHTrackerPoint   &point =*(track->GetTrackPoint(track->GetLayer(),ninlayer));
    //
    point.SetSigmaY(c->GetSigmaImP2());
    point.SetSigmaZ(c->GetSigmaZ2());
    point.SetErrY(sqrt(track->GetErrorY2()));
    point.SetErrZ(sqrt(track->GetErrorZ2()));
    //
    double xyz[3];
    track->GetXYZ(xyz);
    point.SetX(xyz[0]);
    point.SetY(xyz[1]);
    point.SetZ(xyz[2]);
    point.SetAngleY(angle2);
    point.SetAngleZ(track->GetTgl());
    if (point.IsShared()){
      track->SetErrorY2(track->GetErrorY2()*4);
      track->SetErrorZ2(track->GetErrorZ2()*4);
    }
  }  
  
  if (accept>0) return 0;

  //if(c->IsUsed(10)) track->SetNShared(track->GetNShared()+1);
 
  IlcDCHLayer& lay=fSector[layer];
  double chi2p=track->GetChi2();
  // Double_t chi2 = //track->GetChi2();
        // track->GetPredictedChi2(c,
	// 			lay.GetR(),lay.GetWirePhi(c->GetW()),lay.GetTStereoAngle(),
	// 			kUseZCoordinate);

  bool calcChi2=true;
  bool status=track->Update(c,chi2p,dchindex,
			    lay.GetR(),lay.GetWirePhi(c->GetW()),lay.GetTStereoAngle(),
			    kUseZCoordinate,0,calcChi2);
  
  //  track->SetClusterChi2(track->GetLayer(),/*track->GetChi2()-chi2p*/chi2,ninlayer);
  track->SetClusterChi2(track->GetLayer(),track->GetChi2()-chi2p,ninlayer);
  //std::cout<<"chi22calc "<<chi2<<" "<<track->GetChi2()-chi2p<<" "<<"prev "<<chi2p<<" "<<track->GetChi2()<<std::endl;
  return status;
}



Int_t IlcDCHtracker::AcceptCluster(IlcDCHseed * seed, IlcDCHcluster * cluster, Float_t factor, 
                                      Float_t cory, Float_t corz)
{
  //
  // decide according desired precision to accept given 
  // cluster for tracking
  
  seed->SetErrorY2(cluster->GetSigmaImP2()*cory);
  seed->SetErrorZ2(cluster->GetSigmaZ2()*corz);
 
  IlcDCHLayer &lay=fSector[GetLayerFromId(cluster->GetId())];
  Double_t chi2 = seed->GetPredictedChi2(cluster,
					 lay.GetR(),lay.GetWirePhi(cluster->GetW()),
					 lay.GetTStereoAngle(),
					 kUseZCoordinate);

  if(fDebug>15) cout<<"Try to accept cluster chi2="<<chi2
		    <<" trackid "<<cluster->GetLabel(0)<<" "<<cluster->GetLabel(1)<<endl;
  
  if(chi2/(kUseZCoordinate?2:1)>kMaxCHI2cluster) return 3;
  if(chi2/(kUseZCoordinate?2:1)*factor>kMaxCHI2cluster) return 2;
  return 0;

}




//_____________________________________________________________________________
IlcDCHtracker::IlcDCHtracker(const IlcDCHParam *par,IlcDCHwireposition *wpos): 
IlcTracker()
{
  //---------------------------------------------------------------------
  // The main DCH tracker constructor
  //---------------------------------------------------------------------
  if(wpos){
    fSector.Setup(par,wpos); 
  }else{
    IlcDCHwireposition wirepos;
    fSector.Setup(par,&wirepos); 
  }

  fSeeds=0;
  fNtracks = 0;
  kUseZCoordinate=(fabs(par->GetDrop())<1e-4);
  fParam = par;  
  Int_t nlayers = fSector.GetNLayers();

  if (fDebug) { IlcInfo(Form("Number of layers in DCH = %i",nlayers)); }
  
  //
  fRecPar=new IlcDCHrecoParam;
  fDebug     = IlcDebugLevelClass()/10;
  //  fEvent     =0;
  fDebugStreamer = new TTreeSRedirector("DCHdebug.root");
}
//________________________________________________________________________
IlcDCHtracker::IlcDCHtracker(const IlcDCHtracker &t):
  IlcTracker(t)
{
  //------------------------------------
  // dummy copy constructor
  //------------------------------------------------------------------
}
IlcDCHtracker & IlcDCHtracker::operator=(const IlcDCHtracker& /*r*/){
  //------------------------------
  // dummy 
  //--------------------------------------------------------------
  return *this;
}
//_____________________________________________________________________________
IlcDCHtracker::~IlcDCHtracker() {
  //------------------------------------------------------------------
  // DCH tracker destructor
  //------------------------------------------------------------------
  if (fSeeds) {
    fSeeds->Delete(); 
    delete fSeeds;
  }
  if (fDebugStreamer) delete fDebugStreamer;
}

Int_t IlcDCHtracker::GetLayerFromId(ULong_t clid) const{
  return (clid/100000-1)*(fParam->GetRingNum()-1)+(clid%100000)/1000;
}
/*
void IlcDCHtracker::FillESD(TObjArray* arr)
{
  //
  //
  //fill esds using updated tracks
  if (fEvent){
    // write tracks to the event
    // store index of the track
    Int_t nseed=arr->GetEntriesFast();
    //FindKinks(arr,fEvent);
    int nlayers=fSector.GetNLayers();
    for (Int_t i=0; i<nseed; i++) {
      IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i);    
      if (!pt) continue; 
      pt->UpdatePoints();
      //      pt->PropagateTo(fParam->GetInnerRadius());
      if (pt->GetKinkIndex(0)<=0){  //don't propagate daughter tracks 
	pt->PropagateTo(fParam->GetInnerRadius());
      }
 
      if (( pt->GetPoints()[2]- pt->GetPoints()[0])>5 && pt->GetPoints()[3]>0.8){
	IlcESDtrack iotrack;
	iotrack.UpdateTrackParams(pt,IlcESDtrack::kTPCin);
	iotrack.SetTPCPoints(pt->GetPoints());
	iotrack.SetKinkIndexes(pt->GetKinkIndexes());
	iotrack.SetV0Indexes(pt->GetV0Indexes());
	//	iotrack.SetDCHpid(pt->fDCHr);
	//iotrack.SetDCHindex(i);
	fEvent->AddTrack(&iotrack);
	continue;
      }
       
      if ( (pt->GetNumberOfClusters()>0.44*nlayers)&& (Float_t(pt->GetNumberOfClusters())/Float_t(pt->GetNFoundable()))>fRecPar->GetMinDensity()) {
	IlcESDtrack iotrack;
	iotrack.UpdateTrackParams(pt,IlcESDtrack::kTPCin);
	iotrack.SetTPCPoints(pt->GetPoints());
	//iotrack.SetDCHindex(i);
	iotrack.SetKinkIndexes(pt->GetKinkIndexes());
	iotrack.SetV0Indexes(pt->GetV0Indexes());
	//	iotrack.SetDCHpid(pt->fDCHr);
	fEvent->AddTrack(&iotrack);
	continue;
      } 
      //
      // short tracks  - maybe decays

      if ( (pt->GetNumberOfClusters()>0.18*nlayers) && (Float_t(pt->GetNumberOfClusters())/Float_t(pt->GetNFoundable()))>0.70) {
	Int_t found,foundable,shared;
	pt->GetClusterStatistic(0,int(0.4*nlayers),found, foundable,shared,kFALSE);
	if ( (found>20) && (pt->GetNShared()/float(pt->GetNumberOfClusters())<0.2)){
	  IlcESDtrack iotrack;
	  iotrack.UpdateTrackParams(pt,IlcESDtrack::kTPCin);	
	  //iotrack.SetDCHindex(i);
	  iotrack.SetTPCPoints(pt->GetPoints());
	  iotrack.SetKinkIndexes(pt->GetKinkIndexes());
	  iotrack.SetV0Indexes(pt->GetV0Indexes());
	  //iotrack.SetDCHpid(pt->fDCHr);
	  fEvent->AddTrack(&iotrack);
	  continue;
	}
      }       
      
      if ( (pt->GetNumberOfClusters()>0.12*nlayers) && (Float_t(pt->GetNumberOfClusters())/Float_t(pt->GetNFoundable()))>0.8) {
	Int_t found,foundable,shared;
	pt->GetClusterStatistic(0,int(0.4*nlayers),found, foundable,shared,kFALSE);
	if (found<0.12*nlayers) continue;
	if (pt->GetNShared()/float(pt->GetNumberOfClusters())>0.2) continue;
	//
	IlcESDtrack iotrack;
	iotrack.UpdateTrackParams(pt,IlcESDtrack::kTPCin);	
	iotrack.SetTPCPoints(pt->GetPoints());
	iotrack.SetKinkIndexes(pt->GetKinkIndexes());
	iotrack.SetV0Indexes(pt->GetV0Indexes());
	//iotrack.SetDCHpid(pt->fDCHr);
	//iotrack.SetDCHindex(i);
	fEvent->AddTrack(&iotrack);
	continue;
      }   
      // short tracks  - secondaties
      //
      if ( (pt->GetNumberOfClusters()>0.19*nlayers) ) {
	Int_t found,foundable,shared;
	pt->GetClusterStatistic(int(nlayers-1-0.19*nlayers),nlayers-1,found, foundable,shared,kFALSE);
	if ( (found>0.12*nlayers) && (pt->GetNShared()/float(pt->GetNumberOfClusters())<0.2) &&float(found)/float(foundable)>0.8){
	  IlcESDtrack iotrack;
	  iotrack.UpdateTrackParams(pt,IlcESDtrack::kTPCin);	
	  iotrack.SetTPCPoints(pt->GetPoints());
	  iotrack.SetKinkIndexes(pt->GetKinkIndexes());
	  iotrack.SetV0Indexes(pt->GetV0Indexes());
	  //iotrack.SetDCHpid(pt->fDCHr);	
	  //iotrack.SetDCHindex(i);
	  fEvent->AddTrack(&iotrack);
	  continue;
	}
      }       
      
      if ( (pt->GetNumberOfClusters()>=fRecPar->GetMinNOfClusetrs())) {
	if (pt->GetNShared()/float(pt->GetNumberOfClusters())>0.2) continue;
	IlcESDtrack iotrack;
	iotrack.UpdateTrackParams(pt,IlcESDtrack::kTPCin);	
	iotrack.SetTPCPoints(pt->GetPoints());
	iotrack.SetKinkIndexes(pt->GetKinkIndexes());
	iotrack.SetV0Indexes(pt->GetV0Indexes());
	//	iotrack.SetDCHpid(pt->fDCHr);
	//iotrack.SetDCHindex(i);
	fEvent->AddTrack(&iotrack);
	continue;
      } 
      if (pt->GetKinkIndex(0)!=0){
	IlcESDtrack iotrack;
	iotrack.UpdateTrackParams(pt,IlcESDtrack::kTPCin);	
	iotrack.SetTPCPoints(pt->GetPoints());
	iotrack.SetKinkIndexes(pt->GetKinkIndexes());
	iotrack.SetV0Indexes(pt->GetV0Indexes());
	fEvent->AddTrack(&iotrack);
	continue;
      }

    }
  }
  printf("Number of filled ESDs-\t%d\n",fEvent->GetNumberOfTracks());
}
*/

//_____________________________________________________________________________
Double_t IlcDCHtracker::F1old(Double_t x1,Double_t y1,
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature
  //-----------------------------------------------------------------
  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);
  if ( xr*xr+yr*yr<=0.00000000000001) return 100;
  return -xr*yr/sqrt(xr*xr+yr*yr); 
}



//_____________________________________________________________________________
Double_t IlcDCHtracker::F1(Double_t x1,Double_t y1,
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature
  //-----------------------------------------------------------------
  x3 -=x1;
  x2 -=x1;
  y3 -=y1;
  y2 -=y1;
  //  
  Double_t det = x3*y2-x2*y3;
  if (det==0) {
    return 100;
  }
  //
  Double_t u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
  Double_t x0 = x3*0.5-y3*u;
  Double_t y0 = y3*0.5+x3*u;
  Double_t c2 = 1/TMath::Sqrt(x0*x0+y0*y0);
  if (det<0) c2*=-1;
  return c2;
}


Double_t IlcDCHtracker::F2(Double_t x1,Double_t y1,
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature
  //-----------------------------------------------------------------
  x3 -=x1;
  x2 -=x1;
  y3 -=y1;
  y2 -=y1;
  //  
  Double_t det = x3*y2-x2*y3;
  if (det==0) {
    return 100;
  }
  //
  Double_t u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
  Double_t x0 = x3*0.5-y3*u; 
  Double_t y0 = y3*0.5+x3*u;
  Double_t c2 = 1/TMath::Sqrt(x0*x0+y0*y0);
  if (det<0) c2*=-1;
  x0+=x1;
  x0*=c2;  
  return x0;
}



//_____________________________________________________________________________
Double_t IlcDCHtracker::F2old(Double_t x1,Double_t y1,
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature times center of curvature
  //-----------------------------------------------------------------
  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);
  
  return -a/(d*y1-b)*xr/sqrt(xr*xr+yr*yr);
}

//_____________________________________________________________________________
Double_t IlcDCHtracker::F3(Double_t x1,Double_t y1, 
                   Double_t x2,Double_t y2,
                   Double_t z1,Double_t z2) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the tangent of the track dip angle
  //-----------------------------------------------------------------
  return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}


Double_t IlcDCHtracker::F3n(Double_t x1,Double_t y1, 
                   Double_t x2,Double_t y2,
                   Double_t z1,Double_t z2, Double_t c) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the tangent of the track dip angle
  //-----------------------------------------------------------------

  //  Double_t angle1;
  
  //angle1    =  (z1-z2)*c/(TMath::ASin(c*x1-ni)-TMath::ASin(c*x2-ni));
  //
  Double_t d  =  TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
  if (TMath::Abs(d*c*0.5)>1) return 0;
  //  Double_t   angle2    =  TMath::ASin(d*c*0.5);
  //  Double_t   angle2    =  IlcDCHFastMath::FastAsin(d*c*0.5);
  Double_t   angle2    = (d*c*0.5>0.1)? TMath::ASin(d*c*0.5): IlcDCHFastMath::FastAsin(d*c*0.5);

  if(fabs(angle2)<kAlmost0) return 0;

  angle2  = (z1-z2)*c/(angle2*2.);
  return angle2;
}

Bool_t   IlcDCHtracker::GetProlongation(Double_t x1, Double_t x2, Double_t x[5], Double_t &y, Double_t &z)
{//-----------------------------------------------------------------
  // This function find proloncation of a track to a reference plane x=x2.
  //-----------------------------------------------------------------
  
  Double_t dx=x2-x1;

  if (TMath::Abs(x[4]*x1 - x[2]) >= 0.999) {   
    return kFALSE;
  }

  Double_t c1=x[4]*x1 - x[2], r1=sqrt(1.- c1*c1);
  Double_t c2=x[4]*x2 - x[2], r2=sqrt(1.- c2*c2);  
  y = x[0];
  z = x[1];
  
  Double_t dy = dx*(c1+c2)/(r1+r2);
  Double_t dz = 0;
  //
  Double_t delta = x[4]*dx*(c1+c2)/(c1*r2 + c2*r1);
  
  if (TMath::Abs(delta)>0.01){
    dz = x[3]*TMath::ASin(delta)/x[4];
  }else{
    dz = x[3]*IlcDCHFastMath::FastAsin(delta)/x[4];
  }
  
  //dz = x[3]*IlcDCHFastMath::FastAsin(delta)/x[4];

  y+=dy;
  z+=dz;
  
  return kTRUE;  
}

//Int_t  IlcDCHtracker::LoadClusters (TObjArray* fClusters)
Int_t  IlcDCHtracker::LoadClusters (TTree *tree)
{
  //
  //
  TObjArray* fClusters=new TObjArray;
  if (ReadClusters(fClusters,tree)) {
     IlcError("Problem with reading the clusters !");
     return 1;
  }
  //int ncl = fClusters->GetEntriesFast();
  LoadClusters(fClusters);
  
  fClusters->Delete();
  delete fClusters;
  return 0;
}

Int_t  IlcDCHtracker::LoadClusters (TObjArray* fClusters){
  //
  // load clusters to the memory
  //
  fClusters->Sort();
  int ncl = fClusters->GetEntriesFast();

  int ibegin=0,iend=1;
  while(ibegin<ncl) {
    int layer=GetLayerFromId(((IlcDCHcluster*)fClusters->At(ibegin))->GetId());
    while(iend<ncl&&
	  layer==GetLayerFromId(((IlcDCHcluster*)fClusters->At(iend))->GetId())){
      iend++;
    }
    if(layer>=fSector.GetNLayers()){
      IlcError(Form("Layer=%i from cluster id=%i bigger than total numberr of layers=%i",
		    layer,(Int_t)((IlcDCHcluster*)fClusters->At(ibegin))->GetId(),fSector.GetNLayers()));
      ibegin=iend;
      continue;
    }
    IlcDCHLayer * dchlayer = &(fSector[layer]);    
    dchlayer->SetN(iend-ibegin);
    dchlayer->SetClusters(new IlcDCHcluster[dchlayer->GetN()]);
    if(fDebug>=10) 
      if(iend-ibegin)
	cout<<"Load clusters layer "<<layer<<" N "<<dchlayer->GetN();
    for (Int_t i=0;i<dchlayer->GetN();i++) {
      IlcDCHcluster* cl=(IlcDCHcluster*)fClusters->At(i+ibegin);
      cl->SetX(dchlayer->GetR()*TMath::Cos(dchlayer->GetWirePhi(cl->GetW())));
      cl->SetY(dchlayer->GetR()*TMath::Sin(dchlayer->GetWirePhi(cl->GetW())));
//         cl->SetZ(0.);
//         cl->SetZLocal(0.);
      dchlayer->SetCluster(i, *(IlcDCHcluster*)(fClusters->At(i+ibegin)));
      if(fDebug>=10) cout<<" ("<<dchlayer->GetCluster(i)->GetX()<<","
			 <<dchlayer->GetCluster(i)->GetY()<<","
			 <<dchlayer->GetCluster(i)->GetZ()<<","
			 <<Phi_0_2pi(TMath::ATan2(dchlayer->GetCluster(i)->GetY(),dchlayer->GetCluster(i)->GetX()))<<","
			 <<dchlayer->GetCluster(i)->GetImP()<<")";
      dchlayer->GetCluster(i)->SetImP(fabs(dchlayer->GetCluster(i)->GetImP()));
    };
    if(fDebug>=10) if(dchlayer->GetN())cout<<endl;
    ibegin=iend;
  }

  Int_t nlayers = fSector.GetNLayers();
  UInt_t index=0;
  for (Int_t layer = 0;layer<nlayers;layer++){
    IlcDCHLayer*  dchlayer = &(fSector[layer]);  
    Int_t ncl = dchlayer->GetN();
    while (ncl--) {
      IlcDCHcluster *c= (dchlayer->GetCluster(ncl));
      index=(layer<<ClusterIndex::LayerOffset)+(ncl<<ClusterIndex::ClusterOffset);
      dchlayer->InsertCluster(c,index);
    }
    dchlayer->SetCoeffConv(1./dchlayer->GetMaxZ()/1.02);
    //
    // write indexes for fast acces
    //
    for (Int_t i=0;i<IlcDCHLayer::NFAST;i++)
      dchlayer->SetFastClusterZ(i,-1);
    for (Int_t i=0;i<dchlayer->GetN();i++){
      Int_t zi = Int_t(((*dchlayer)[i]->GetZ()*dchlayer->GetCoeffConv()+1.)*IlcDCHLayer::NFAST/2);
      if(zi<0||zi>=IlcDCHLayer::NFAST) 
	Warning("LoadOuterSectors","Cluster Z out of detector size cl->z:%f, Z of detector:%f",
		(*dchlayer)[i]->GetZ(),dchlayer->GetMaxZ());
      dchlayer->SetFastClusterZ(zi,i);  // write index
    }
    Int_t last = 0;
    for (Int_t i=0;i<IlcDCHLayer::NFAST;i++){
      if (dchlayer->GetFastClusterZ(i)<0)
	dchlayer->SetFastClusterZ(i,last);
      else
	last = dchlayer->GetFastClusterZ(i);
    }
  }
  if (fDebug) { IlcInfo(Form("Number of loaded clusters: %i",ncl)); }

  return 0;
}

//_____________________________________________________________________________
Int_t IlcDCHtracker::ReadClusters(TObjArray *array, TTree *clusterTree) const
{
  //
  // Reads IlcDCHcluster) 
  // from the file. The names of the cluster tree and branches 
  // should match the ones used in IlcDCHclusterizer::WriteClusters()
  //
  if(!clusterTree){
    IlcError("Can't get the branch !");
    return 1;    
  }

  TClonesArray *clusterArray = new TClonesArray("IlcDCHcluster"); 
  
  TBranch *branch = clusterTree->GetBranch("DCHcluster");
  if (!branch) {
    IlcError("Can't get the branch !");
    return 1;
  }
  branch->SetAddress(&clusterArray); 
  
  // Loop through all entries in the tree
  Int_t nEntries   = (Int_t) clusterTree->GetEntries();
  Int_t nbytes     = 0;
  IlcDCHcluster *c = 0;
  int nclusters=0;
  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {    
    
    // Import the tree
    nbytes += clusterTree->GetEvent(iEntry);  
    
    // Get the number of points in the detector
    Int_t nCluster = clusterArray->GetEntriesFast();  
    
    // Loop through all DCH digits
    for (Int_t iCluster = 0; iCluster < nCluster; iCluster++) { 
      c = (IlcDCHcluster *) clusterArray->UncheckedAt(iCluster);
      //      cout<<iCluster<<" "<<c->GetImP()<<" "<<c->IsUsed(100)<<" from "<<nCluster<<endl;
      if(c->IsUsed(100)) continue;
      IlcDCHcluster *co = new IlcDCHcluster(*c);
      array->AddLast(co);
      clusterArray->RemoveAt(iCluster); 
      nclusters++;
    }

  }
  //  IlcInfo(Form("Number of loaded clusters: %i",nclusters));
  delete clusterArray;

  return 0;

}

void IlcDCHtracker::UnloadClusters()
{
  //
  // unload clusters from the memory
  //
  Int_t nlayers = fSector.GetNLayers();
  for (Int_t layer = 0;layer<nlayers;layer++){
    IlcDCHLayer&  dchlayer = fSector[layer];
    dchlayer.ResetClusters();
  }
}
  
void   IlcDCHtracker::Transform(IlcCluster * /*cluster*/){
  //
  // 
  //
  return;
}

//_________________________________________________________________________
IlcDCHcluster *IlcDCHtracker::GetCluster(Int_t index) const {
  //--------------------------------------------------------------------
  //       Return pointer to a given cluster
  //--------------------------------------------------------------------
  if (index<0) return 0; // no cluster
  Int_t layer=(index&ClusterIndex::Layer)>>ClusterIndex::LayerOffset; 
  Int_t ncl=(index&ClusterIndex::Cluster)>>ClusterIndex::ClusterOffset;

  const IlcDCHLayer * dchlayer=0;
  IlcDCHcluster * cllayer =0;


  dchlayer = &(fSector[layer]);
  if (dchlayer==0) return 0;
  if (dchlayer->GetN()<=ncl) return 0;
  cllayer = dchlayer->GetClusters();
  return &(cllayer[ncl]);      
  
}



Int_t IlcDCHtracker::FollowToNext(IlcDCHseed& t, Int_t nr) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation to next pad layer
  //-----------------------------------------------------------------
  //
  Double_t  x= GetRlayer(nr,t.GetZ());
  IlcDCHcluster *cl=0;
  Int_t dchindex= t.GetClusterIndex(nr,0);
  bool cannotProp=false;
  //
  if (fIteration>0 && dchindex>=0){  //if we have already clusters 
    if (fDebug>15) cout<<"Follow to next layer "<<nr<<" "<<dchindex<<endl;
    if (!t.PropagateTo(x)){
      cannotProp=true;
      if((dchindex&ClusterIndex::NotAccept)) return -1;
    }
    t.SetLayer(nr);
    //
    for(int j=0;j<kMaxInLayer;j++){
      dchindex= t.GetClusterIndex(nr,j);
      if(dchindex<0) break;
      cl = t.GetClusterPointer(nr,j);
      if (!cl) cl = GetCluster(dchindex);
     
      t.SetCurrentCluster(cl,dchindex); 
      if (fDebug>15) cout<<"Try Accept cluster layer "<<nr<<" inlay "<<j<<" yz "<<cl->GetY()<<" "<<cl->GetZ()<<endl;
      Int_t accept = AcceptCluster(&t,t.GetCurrentCluster(),1.);
      if ((dchindex&ClusterIndex::NotAccept)==0) accept =0;
      if (accept<3) { 
	//if founded cluster is acceptible
	if (cl->IsUsed(11)) {  // id cluster is shared inrease uncertainty
	  t.SetErrorY2(t.GetErrorY2()*3);
	  t.SetErrorZ2(t.GetErrorZ2()*3); 
	}
	t.SetNFoundable(t.GetNFoundable()+1);
	UpdateTrack(&t,accept);
      }
    }  
  }
  if(cannotProp) return -1;
  if (TMath::Abs(t.GetSnp())>IlcDCHReconstructor::GetMaxSnpTracker()) return 0;  // cut on angle
  if (fIteration>1){
    // not look for new cluster during refitting    
    if(dchindex!=-1)
      t.SetNFoundable(t.GetNFoundable()+1);
    return 0;
  }
  //
//  UInt_t index=0;
  if (fDebug>15) {
    cout<<"Current Position ";t.Print("params");
    cout<<"Propagate to layer="<<nr<<" r="<<x<<endl;
  }
  //
  if (!t.PropagateTo(x)) {
    if (fIteration==0) t.SetRemoval(10);
    return -1;
  }
  Double_t y=t.GetPhi(); 
  Double_t z=t.GetZ();
  //

  if (!IsActive(nr)) {
    t.SetInDead(kTRUE);
    t.SetClusterIndex(nr,-1,0); 
    return 0;
  }

  const IlcDCHLayer &klayer=GetLayer(nr);
  if ( (t.GetSigmaY2()<0) || t.GetSigmaZ2()<0) return 0;

  //adjust position forr hyperboloid
  if(klayer.IsStereo()){
    x=klayer.GetRAtZ(t.GetZ());
    if (!t.PropagateTo(x)) {
      if (fIteration==0) t.SetRemoval(10);
      return -1;
    }
    
    y=t.GetPhi(); 
    z=t.GetZ();
  }

  Double_t  roady  =fRecPar->GetRoadY()+fabs((GetRlayer(nr)-GetRlayer(nr>0?nr-1:nr+1))*1*t.GetSnp());
  Double_t  roadz = fRecPar->GetRoadZ()+fabs((GetRlayer(nr)-GetRlayer(nr>0?nr-1:nr+1))*2*t.GetTgl());
  if (fDebug>15) cout<<"Follow to next: Try Find new cluster yz "<<y<<" "<<z<<" "<<roady<<" "<<roadz<<endl; 
  //
  if(fDebug>15) 
    cout<<"Good Rigion x="<<x<<"  z="<<z<<
      " max="<<IlcDCHReconstructor::GetCtgRange()*x+fRecPar->GetSafetyZ()<<
      " zlength="<<klayer.GetMaxZ()<<
      " Track snp="<<t.GetSnp()<<" max="<<IlcDCHReconstructor::GetMaxSnpTracker()<<" fndb "<<t.GetNFoundable()<<endl;
  if (!(TMath::Abs(z)<(IlcDCHReconstructor::GetCtgRange()*x+fRecPar->GetSafetyZ()) && 
	TMath::Abs(z)<klayer.GetMaxZ()+fRecPar->GetSafetyZ() && 
	(TMath::Abs(t.GetSnp())<IlcDCHReconstructor::GetMaxSnpTracker())))
    return 0;

  if(TMath::Abs(z)<=klayer.GetMaxZ()){
    if(klayer.IsCrossedCell(y,z,t.GetSnp(),t.GetTgl())){
      t.SetNFoundable(t.GetNFoundable()+1);
    }else{
      t.SetClusterIndex(nr,-1,0);
    }
  }
  int found=FindAndUpdate(t,nr,roady,roadz);  //calculate 
  if(!found){  
    if ( fIteration==0 && t.GetNFoundable()*fRecPar->GetMinDensity() > t.GetNumberOfClusters()) {
      if(fDebug>15) 
	cout<<"Set stopped track "<<t.GetNumberOfClusters()<<" "<<t.GetNFoundable()
	    <<" "<<fRecPar->GetMinDensity()<<endl;
      t.SetRemoval(10);
    }
  }
  return 1;
}

Int_t IlcDCHtracker::FindAndUpdate(IlcDCHseed& t, Int_t nr,double roady,double roadz){
  const IlcDCHLayer &klayer=GetLayer(nr);
  IlcDCHcluster *cl=0;
  
  if (klayer) {
    double y=t.GetPhi(); 
    double z=t.GetZ();

    IlcDCHcluster* clar[5]={0,0,0,0,0};
    UInt_t indar[5]={0,0,0,0,0};
    klayer.FindNearest2(y,t.GetSnp(),z,roady,roadz,4,clar,indar);
    
    std::vector<std::pair<double,int> > chi2;
    for(int i=0;i<5;i++){
      if(clar[i]&&!t.IsAlready(nr,clar[i])){
	if(fDebug>15) cout<<i<<" - ";
	double chi2n=t.GetPredictedChi2(clar[i],klayer.GetR(),klayer.GetWirePhi(clar[i]->GetW()),klayer.GetTStereoAngle(),kUseZCoordinate);
	chi2.push_back(std::pair<double,int>(chi2n,i));
      }
    }
    std::sort(chi2.begin(),chi2.end());
    
    for(int i=0;i<(int)chi2.size();i++){
      int ii=chi2[i].second;
      if(clar[ii]&&!t.IsAlready(nr,clar[ii])){
	cl=clar[ii];
	UInt_t index=indar[ii];
	t.SetCurrentCluster(cl,klayer.GetIndex(index)); 
	t.SetLayer(nr);
	if (fIteration==2&&cl->IsUsed(10)&&i==0) return 0; 
	Int_t accept = AcceptCluster(&t,t.GetCurrentCluster(),1.);
	if(fDebug>15) cout<<"Accept cluster yz "<<cl->GetY()<<" "<<cl->GetZ()<<" accept "<<accept<<endl;
	if (fIteration==2&&cl->IsUsed(11)) {
	  t.SetErrorY2(t.GetErrorY2()*3);
	  t.SetErrorZ2(t.GetErrorZ2()*3); 
	}
	if (accept<3) UpdateTrack(&t,accept);
      }
    }
  }
  return cl?1:0;

}


//_________________________________________________________________________
Bool_t IlcDCHtracker::GetTrackPoint(Int_t index, IlcTrackPoint &p ) const
{
  // Get track space point by index
  // return false in case the cluster doesn't exist
  IlcDCHcluster *cl = GetCluster(index);
  if (!cl) return kFALSE;
  Float_t xyz[3];
  xyz[0] = cl->GetX();
  xyz[1] = cl->GetY();
  xyz[2] = cl->GetZ();
  Float_t cov[6];
  Float_t sigmaY2 = cl->GetSigmaY2();
  Float_t sigmaZ2 = cl->GetSigmaZ2();
  cov[0] = sigmaY2;
  cov[1] = 0.;
  cov[2] = sigmaY2;
  cov[3] = 0.;
  cov[4] = 0.;
  cov[5] = sigmaZ2;
  p.SetXYZ(xyz[0],xyz[1],xyz[2],cov);
  //IlcAlignObj::ELayerID iLayer;
  UShort_t volid = 0; //IlcAlignObj::LayerToVolUID(iLayer,idet);
  p.SetVolumeID(volid);
  return kTRUE;
}



Int_t IlcDCHtracker::UpdateClusters(IlcDCHseed& t,  Int_t nr) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation to next pad layer
  //-----------------------------------------------------------------
  t.SetCurrentCluster(0,-1);
  
  Double_t xt= t.GetR();
  Int_t     layer = GetLayerNumber(xt,t.GetZ())-1; 

  if (layer < nr) return 1; // don't prolongate if not information until now -

  Double_t x= GetRlayer(nr,t.GetZ());
  Double_t y,z;

  if (!t.PropagateTo(x)){
    return 0;
  }
  //
  y=t.GetPhi();
  z=t.GetZ();

  //
  if (TMath::Abs(t.GetSnp())>IlcDCHReconstructor::GetMaxSnpTracker()) return 0;

  if (!IsActive(nr)) {
    t.SetInDead(kTRUE);
    t.SetClusterIndex(nr,-1,0); 
    return 0;
  }
  //IlcInfo(Form("A - Sector%d phi %f - alpha %f", t.fRelativeSector,y/x, t.GetAlpha()));

  IlcDCHLayer &klayer=GetLayer(nr);

  //adjust position forr hyperboloid
  if(klayer.IsStereo()){
    x=klayer.GetRAtZ(t.GetZ());
    if (!t.PropagateTo(x)) return 0;
    
    y=t.GetPhi(); 
    z=t.GetZ();
  }

  if (!(TMath::Abs(t.GetZ())<(IlcDCHReconstructor::GetCtgRange()*t.GetX()+fRecPar->GetSafetyZ())&&
	(TMath::Abs(t.GetSnp())<IlcDCHReconstructor::GetMaxSnpTracker())))
    return 0;

  if(TMath::Abs(t.GetZ())<klayer.GetMaxZ()){
    if(klayer.IsCrossedCell(y,z,t.GetSnp(),t.GetTgl())){
      t.SetNFoundable(t.GetNFoundable()+1);
    }else{
      t.SetClusterIndex(nr,-1,0);
    }
  }
  //
  Double_t roady = fRecPar->GetRoadY();
  Double_t roadz = fRecPar->GetRoadZ();
  //

  IlcDCHcluster *cl=0;
  Int_t index = t.GetClusterIndex(nr,0);    

  if ( (index>=0) && (index&ClusterIndex::NotAccept)==0){
    cl = t.GetClusterPointer(nr,0);
    if ( (cl==0) && (index>=0)) cl = GetCluster(index);
    if (cl) {
      t.SetCurrentCluster(cl,index);
      return 1;
    }
  }
  
  UInt_t uindex = 0;

  if (klayer.GetN()>0) {    
    IlcDCHcluster* clar[5]={0,0,0,0,0};
    UInt_t indar[5]={0,0,0,0,0};
    klayer.FindNearest2(y,t.GetSnp(),z,roady,roadz,4,clar,indar);

    if(clar[0]){
      std::vector<std::pair<double,int> > chi2;
      for(int i=0;i<5;i++){
	if(clar[i]){
	  double chi2n=t.GetPredictedChi2(clar[i],klayer.GetR(),klayer.GetWirePhi(clar[i]->GetW()),klayer.GetTStereoAngle(),kUseZCoordinate);
	  chi2.push_back(std::pair<double,int>(chi2n,i));
	}
      }
      std::sort(chi2.begin(),chi2.end());
      cl=clar[chi2[0].second];
      uindex=indar[chi2[0].second];
    }
  }

  if (cl) {
    t.SetCurrentCluster(cl,klayer.GetIndex(uindex));
  }

  return 1;
}


Int_t IlcDCHtracker::FollowToNextCluster(IlcDCHseed & t, Int_t nr) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation to next pad layer
  //-----------------------------------------------------------------

  //update error according neighborhoud

  if (t.GetCurrentCluster()) {
    t.SetLayer(nr); 
    Int_t accept = AcceptCluster(&t,t.GetCurrentCluster(),1.);
    
    if (t.GetCurrentCluster()->IsUsed(10)){
      //
      //
      t.SetNShared(t.GetNShared()+1);
      if (t.GetNShared()>0.7*t.GetNumberOfClusters()) {
	t.SetRemoval(10);
	return 0;
      }
    }   
    if (fIteration>0) accept = 0;
    if (accept<3)  UpdateTrack(&t,accept);  
    if(accept<3) {
      Double_t roady = fRecPar->GetRoadY();
      Double_t roadz = fRecPar->GetRoadZ();
      FindAndUpdate(t,nr,roady,roadz);
    }
  } else {
    if (fIteration==0){
      if ( ( (t.GetSigmaY2()+t.GetSigmaZ2())>fRecPar->GetCutTrackSigma2())&& 
	   t.GetNumberOfClusters()>=fRecPar->GetMinNOfClusetrs()) t.SetRemoval(10);      
      if (  t.GetNumberOfClusters()>=fRecPar->GetMinNOfClusetrs()&&
	    t.GetChi2()/((kUseZCoordinate?2:1)*t.GetNumberOfClusters()-5)>kMaxCHI2track ) t.SetRemoval(10);      

    }
    if((t.GetSigmaY2()+t.GetSigmaZ2())>1e20) t.SetRemoval(10);
  }
  return 1;
}



//_____________________________________________________________________________
Int_t IlcDCHtracker::FollowProlongation(IlcDCHseed& t, Int_t rf, Int_t step,int contin) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation.
  //-----------------------------------------------------------------
  Double_t xt=t.GetR();
  int nlayers=fSector.GetNLayers();
  //

  Int_t first = GetLayerNumber(xt,t.GetZ())-(t.GetNumberOfClusters()>0?1:0);
  if(first<nlayers-1&&!contin){
    if((t.GetClusterIndex(first+1,0)&ClusterIndex::NotAccept)==0)
      first++;
  }
  if(contin) first=TMath::Min(first,t.GetLayer()-1);

  // std::cout<<"FollowProlongation r "<<xt<<" z "<<t.GetZ()<<" current layer "<<t.GetLayer()
  // 	   <<" layer position "<<GetLayerNumber(xt,t.GetZ())<<" will be first "<<first<<std::endl;

  
  for (Int_t nr= first; nr>=rf; nr-=step) {
    // update kink info
    if (t.GetKinkIndexes()[0]>0){
      for (Int_t i=0;i<3;i++){
	Int_t index = t.GetKinkIndexes()[i];
	if (index==0) break;
	if (index<0) continue;
	///
        /*
	IlcKink * kink = 0;//(IlcKink*)fEvent->GetKink(index-1);
	if (!kink){
	  printf("PROBLEM\n");
	}
	else{
	  Int_t kinklayer = kink->GetTPCRow0()+2+Int_t(0.5/(0.05+kink->GetAngle(2)));
	  if (kinklayer==nr){
	    IlcExternalTrackParam paramd(t);
	    kink->SetDaughter(paramd);
	    kink->SetStatus(2,5);
	    kink->Update();
	  }
	}*/
      }
    }

    if (nr==int(nlayers*0.5)) t.UpdateReference();
    if (FollowToNext(t,nr)<=0) 
      if (!t.IsActive()) 
	return 0;
    
  }   
  return 1;
}

Int_t IlcDCHtracker::FollowBackProlongation(IlcDCHseed& t, Int_t rf) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation.
  //-----------------------------------------------------------------
  //
  //std::cout<<"from first to last "<<t.GetFirstPoint()<<" "<<t.GetLastPoint()<<endl;
  Double_t xt=t.GetR();  
    
  Int_t first = t.GetFirstPoint();
  int lay=GetLayerNumber(xt,t.GetZ());
  if (first<lay) first = lay;
  //
  if (first<0) first=0;
  for (Int_t nr=first; nr<=rf; nr++) {
    //    if ( (TMath::Abs(t.GetSnp())>0.95)) break;//patch 28 fev 06
    if (t.GetKinkIndexes()[0]<0){
      for (Int_t i=0;i<3;i++){
	/*Int_t index = t.GetKinkIndexes()[i];
	if (index==0) break;
	if (index>0) continue;
	index = TMath::Abs(index);
	IlcKink * kink = 0;//(IlcKink*)fEvent->GetKink(index-1);
	if (!kink){
	  printf("PROBLEM\n");
	}
	else{
	  Int_t kinklayer = kink->GetTPCRow0()-2-Int_t(0.5/(0.05+kink->GetAngle(2)));
	  if (kinklayer==nr){
	    IlcExternalTrackParam paramm(t);
	    kink->SetMother(paramm);
	    kink->SetStatus(2,1);
	    kink->Update();
	  }
	  }*/
      }      
    }
    //

    if(FollowToNext(t,nr)<0) break;    
    if(TMath::Abs(t.GetZ())>fSector[nr].GetMaxZ()+fRecPar->GetSafetyZ()) break;
  }   
  return 1;
}




   
Float_t IlcDCHtracker::OverlapFactor(IlcDCHseed * s1, IlcDCHseed * s2, Int_t &sum1, Int_t & sum2)
{
  //
  //
  sum1=0;
  sum2=0;
  Int_t sum=0;
  //
  //  Int_t offset =0;
  Int_t firstpoint = TMath::Min(s1->GetFirstPoint(),s2->GetFirstPoint());
  Int_t lastpoint = TMath::Max(s1->GetLastPoint(),s2->GetLastPoint());
  if (lastpoint>kMaxLayer) 
    lastpoint =kMaxLayer;
  if (firstpoint<0) 
    firstpoint = 0;
  if (firstpoint>lastpoint) {
    firstpoint =lastpoint;
    //    lastpoint  =kMaxLayer;
  }
    
  
  for (Int_t i=firstpoint-1;i<lastpoint+1;i++){
    for(int j=0;j<kMaxInLayer;j++)
      if (s1->GetClusterIndex(i,j)>0) sum1++;
    for(int j=0;j<kMaxInLayer;j++)
      if (s2->GetClusterIndex(i,j)>0) sum2++;
    for(int j1=0;j1<kMaxInLayer;j1++)
      for(int j2=0;j2<kMaxInLayer;j2++)
	if (s1->GetClusterIndex(i,j1)==s2->GetClusterIndex(i,j2) && s1->GetClusterIndex(i,j1)>=0) {
	  sum++;
	}
  }
  if (sum<5) return 0;

  Float_t summin = TMath::Min(sum1+1,sum2+1);
  Float_t ratio = (sum+1)/Float_t(summin);
  return ratio;
}

void  IlcDCHtracker::SignShared(IlcDCHseed * s1, IlcDCHseed * s2)
{
  //
  //
  if (TMath::Abs(s1->GetC()-s2->GetC())>fRecPar->GetNearestByC()) return;
  if (TMath::Abs(s1->GetTgl()-s2->GetTgl())>fRecPar->GetNearestByTgl()) return;

  Float_t dz2 =(s1->GetZ() - s2->GetZ());
  dz2*=dz2;
  Float_t dy2 =(s1->GetY() - s2->GetY());
  dy2*=dy2;
  Float_t distance = dz2+dy2;
  if (distance>fRecPar->GetNearestByDistance2()) return ; // if there are far away  - not overlap - to reduce combinatorics
  
  //
  Int_t sumshared=0;
  //
  Int_t firstpoint = TMath::Max(s1->GetFirstPoint(),s2->GetFirstPoint());
  Int_t lastpoint = TMath::Min(s1->GetLastPoint(),s2->GetLastPoint());
  //
  if (firstpoint>=lastpoint-5) return;;

  for (Int_t i=firstpoint;i<lastpoint;i++){
    for(int j1=0;j1<kMaxInLayer;j1++)
      for(int j2=0;j2<kMaxInLayer;j2++)
	if (s1->GetClusterIndex(i,j1)==s2->GetClusterIndex(i,j2) && s1->GetClusterIndex(i,j1)>=0) {
	  sumshared++;
	}
  }
  if (sumshared>4&&s1->IsActive()&&s2->IsActive()){
    // sign clusters
    //
    for (Int_t i=firstpoint;i<lastpoint;i++){
      for(int j1=0;j1<kMaxInLayer;j1++)
	for(int j2=0;j2<kMaxInLayer;j2++)
	  if ( (s1->GetClusterIndex(i,j1))==(s2->GetClusterIndex(i,j2)) && s1->GetClusterIndex(i,j1)>=0) {
	    IlcDCHTrackerPoint *p1  = s1->GetTrackPoint(i,j1);
	    IlcDCHTrackerPoint *p2  = s2->GetTrackPoint(i,j2);; 
	    p1->SetShared(kTRUE);
	    p2->SetShared(kTRUE); 
	  }
    }
  }
  //  
  if (sumshared>10){
    for (Int_t i=0;i<4;i++){
      if (s1->GetOverlapLabel(3*i)==0){
	s1->SetOverlapLabel(3*i,  s2->GetLabel());
	s1->SetOverlapLabel(3*i+1,sumshared);
	s1->SetOverlapLabel(3*i+2,s2->GetUniqueID());
	break;
      }	
    }
    for (Int_t i=0;i<4;i++){
      if (s2->GetOverlapLabel(3*i)==0){
	s2->SetOverlapLabel(3*i,  s1->GetLabel());
	s2->SetOverlapLabel(3*i+1,sumshared);
	s2->SetOverlapLabel(3*i+2,s1->GetUniqueID());
	break;
      }	
    }    
  }
  
}

void  IlcDCHtracker::SignShared(TObjArray * arr)
{
  //
  //sort trackss according sectors
  //  
  for (Int_t i=0; i<arr->GetEntriesFast(); i++) {
    IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i);    
    if (!pt) continue;
    //if (pt) RotateToLocal(pt);
    pt->SetSort(0);
  }
  arr->UnSort();
  arr->Sort();  // sorting according z
  arr->Expand(arr->GetEntries());
  //
  //
  Int_t nseed=arr->GetEntriesFast();
  for (Int_t i=0; i<nseed; i++) {
    IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i);    
    if (!pt) continue;
    for (Int_t j=0;j<12;j++){
      pt->SetOverlapLabel(j,0);
    }
  }
  for (Int_t i=0; i<nseed; i++) {
    IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i);    
    if (!pt) continue;
    if (pt->GetRemoval()>10) continue;
    for (Int_t j=i+1; j<nseed; j++){
      IlcDCHseed *pt2=(IlcDCHseed*)arr->UncheckedAt(j);
      //      if (pt2){
      if (pt2->GetRemoval()<=10) {
	SignShared(pt,pt2);
      }
    }  
  }
}

void  IlcDCHtracker::RemoveDouble(TObjArray * arr, Float_t factor1, Float_t factor2,  Int_t removilcndex)
{
  //
  //sort trackss according sectors
  //
  if (fDebug&1) {
    Info("RemoveDouble","Number of tracks before double removal- %d\n",arr->GetEntries());
  }
  //
  for (Int_t i=0; i<arr->GetEntriesFast(); i++) {
    IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i);    
    if (!pt) continue;
    pt->SetSort(0);
  }
  arr->UnSort();
  arr->Sort();  // sorting according z
  arr->Expand(arr->GetEntries());
  //
  //reset overlap labels
  //
  Int_t nseed=arr->GetEntriesFast();
  for (Int_t i=0; i<nseed; i++) {
    IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i);    
    if (!pt) continue;
    pt->SetUniqueID(i);
    for (Int_t j=0;j<12;j++){
      pt->SetOverlapLabel(j,0);
    }
  }
  //
  //sign shared tracks
  for (Int_t i=0; i<nseed; i++) {
    IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i);    
    if (!pt) continue;
    if (pt->GetRemoval()>10) continue;
    Float_t deltac = pt->GetC()*0.1;
    for (Int_t j=i+1; j<nseed; j++){
      IlcDCHseed *pt2=(IlcDCHseed*)arr->UncheckedAt(j);
      //      if (pt2){
      if (pt2->GetRemoval()<=10) {
	if (TMath::Abs(pt->GetC()  -pt2->GetC())>deltac) continue;
	//if (TMath::Abs(pt->GetTgl()-pt2->GetTgl())>0.05) continue;
	//Checking in SignShared
	SignShared(pt,pt2);
      }
    }
  }
  //
  // remove highly shared tracks
  for (Int_t i=0; i<nseed; i++) {
    IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i);    
    if (!pt) continue;
    if (pt->GetRemoval()>10) continue;
    //
    Int_t sumshared =0;
    for (Int_t j=0;j<4;j++){
      sumshared = pt->GetOverlapLabel(j*3+1);      
    }
    Float_t factor = factor1;
    if (pt->GetRemoval()>0) factor = factor2;
    if (sumshared/pt->GetNumberOfClusters()>factor){
      for (Int_t j=0;j<4;j++){
	if (pt->GetOverlapLabel(3*j)==0) continue;
	if (pt->GetOverlapLabel(3*j+1)<5) continue; 
	if (pt->GetRemoval()==removilcndex) continue;      
	IlcDCHseed * pt2 = (IlcDCHseed*)arr->UncheckedAt(pt->GetOverlapLabel(3*j+2));
	if (!pt2) continue;
	if (pt2->GetSigma2C()<pt->GetSigma2C()){
	  //	  pt->fRemoval = removilcndex;
	  delete arr->RemoveAt(i);	  
	  break;
	}
      }      
    }
  }
  arr->Compress();
  if (fDebug&1) {
    Info("RemoveDouble","Number of tracks after double removal- %d\n",arr->GetEntries());
  }
}





void IlcDCHtracker::SortTracks(TObjArray * arr, Int_t mode) const
{
  //
  //sort tracks in array according mode criteria
  Int_t nseed = arr->GetEntriesFast();    
  for (Int_t i=0; i<nseed; i++) {
    IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }
    pt->SetSort(mode);
  }
  arr->UnSort();
  arr->Sort();
}

void IlcDCHtracker::RemoveUsed(TObjArray * arr, Float_t factor1,  Float_t factor2, Int_t removilcndex)
{

  //Loop over all tracks and remove "overlaps"
  //
  //
  Int_t nseed = arr->GetEntriesFast();  
  Int_t good =0;
  int nlayers=fSector.GetNLayers();
  if (fDebug>0){
    Info("RemoveUsed","\n*****\nNumber of good tracks before shared removal\t%d\n",nseed);
  }

  for (Int_t i=0; i<nseed; i++) {
    IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i);    
    if (!pt) {
      delete arr->RemoveAt(i);
    }
    else{
      pt->SetSort(1);
      pt->SetBSigned(kFALSE);
    }
  }
  arr->Compress();
  nseed = arr->GetEntriesFast();
  arr->UnSort();
  arr->Sort();
  //
  //unsign used
  UnsignClusters();
  //
  for (Int_t i=0; i<nseed; i++) {
    IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }    
    Int_t found,foundable,shared;
    if (pt->IsActive()) 
      pt->GetClusterStatistic(0,nlayers-1,found, foundable,shared,kFALSE);
    else
      pt->GetClusterStatistic(0,nlayers-1,found, foundable,shared,kTRUE); 
    //
    Double_t factor = factor2;
    if (pt->GetBConstrain()) factor = factor1;

    if ((Float_t(shared)/Float_t(found))>factor){
      pt->Desactivate(removilcndex);
      continue;
    }

    good++;
    for (Int_t i=0; i<nlayers; i++) 
      for(int j=0;j<kMaxInLayer;j++){
	Int_t dchindex=pt->GetClusterIndex(i,j);
	if (dchindex<0 || dchindex&ClusterIndex::NotAccept ) continue;
	IlcDCHcluster *c= pt->GetClusterPointer(i,j);        
	if (!c) continue;
	c->Use(10);  
      }
  }
  fNtracks = good;
  if (fDebug>0){
    Info("RemoveUsed","\n*****\nNumber of good tracks after shared removal\t%d\n",fNtracks);
  }
}


void IlcDCHtracker::RemoveUsed2(TObjArray * arr, Float_t factor1,  Float_t factor2, Int_t minimal)
{

  //Loop over all tracks and remove "overlaps"
  //
  //
  UnsignClusters();
  //
  Int_t nseed = arr->GetEntriesFast();  
  double * quilcty = new double[nseed];
  Int_t   * indexes = new Int_t[nseed];
  Int_t good =0;
  //
  if (IlcDebugLevelClass()>0){
    Info("RemoveUsed2","\n*****\nNumber of good tracks before shared removal\t%d\n",nseed);
  }
  //
  for (Int_t i=0; i<nseed; i++) {
    IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i);    
    if (!pt){
      quilcty[i]=-1;
      continue;
    }
    pt->UpdatePoints();    //select first last max dens points
    Float_t * points = pt->GetPoints();
    if (points[3]<0.8) quilcty[i] =-1;
    //
    quilcty[i] = (points[2]-points[0]+1)*points[3]+
      pt->GetNumberOfClusters()*(1.+0.3*TMath::Prob(pt->GetChi2(),(kUseZCoordinate?2:1)*pt->GetNumberOfClusters()-5));
      //-pt->GetChi2()/((kUseZCoordinate?2:1)*pt->GetNumberOfClusters()-5)*1e-4;
    if(IlcDebugLevelClass()>4){
      CookLabel(pt,0.1);
      cout<<i<<" "<<quilcty[i]<<" "<<pt->GetNumberOfClusters()<<" densed from "<<points[2]<<" to "<<points[0]
	  <<" maxdens="<<points[3]<<" index="<<points[1]<<" lbl="<<pt->GetLabel()<<" fake="<<pt->GetFakeRatio()
	  <<" chi2="<<pt->GetChi2()<<" prob="<<TMath::Prob(pt->GetChi2(),(kUseZCoordinate?2:1)*pt->GetNumberOfClusters()-5)<<endl;
    }
  }
  TMath::Sort(nseed,quilcty,indexes);
  //
  //
  for (Int_t itrack=0; itrack<nseed; itrack++) {
    Int_t trackindex = indexes[itrack];
    IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(trackindex);    
    if (quilcty[trackindex]<0){
      if (pt) {
	delete arr->RemoveAt(trackindex);
      }
      else{
	arr->RemoveAt(trackindex);
      }
      continue;
    }
    //
//     Int_t first = Int_t(pt->GetPoints()[0]);
//     Int_t last  = Int_t(pt->GetPoints()[2]);
    Int_t first = pt->GetFirstPoint();
    Int_t last  = pt->GetLastPoint();
    Double_t factor = (pt->GetBConstrain()) ? factor1: factor2;
    //
    Int_t found,foundable,shared;
    if(fDebug>15) cout<<"seedid "<<trackindex<<" ";
    pt->GetClusterStatistic(first,last, found, foundable,shared,kFALSE);
    Float_t sharedfactor = Float_t(shared+1)/Float_t(found+1);
    Bool_t itsgold =kFALSE;
    //   if (pt->GetESD()){
    //  if (pt->GetESD()->GetVXDclusters(0)>4) itsgold= kTRUE;
    // }
    if (!itsgold){
      //
      if (Float_t(shared+1)/Float_t(found+1)>factor){
	if (pt->GetKinkIndexes()[0]!=0) continue;  //don't remove tracks  - part of the kinks
	delete arr->RemoveAt(trackindex);
	continue;
      }      
      if (pt->GetNumberOfClusters()<fRecPar->GetMinNOfClusetrs()&&(found-0.5*shared)<minimal){  //remove short tracks
	if (pt->GetKinkIndexes()[0]!=0) continue;  //don't remove tracks  - part of the kinks
	delete arr->RemoveAt(trackindex);
	continue;
      }
    }

    good++;
    if (sharedfactor>0.4) continue;
    if (pt->GetKinkIndexes()[0]>0) continue;
    if(fDebug>15) cout<<"Take It "<<trackindex<<" "<<endl;
    for (Int_t i=first; i<last; i++)
      for(int j=0;j<kMaxInLayer;j++){
	Int_t dchindex=pt->GetClusterIndex(i,j);
	if (dchindex<0 || dchindex&ClusterIndex::NotAccept ) continue;
	IlcDCHcluster *c= pt->GetClusterPointer(i,j);        
	if (!c) continue;
	c->Use(10);  
      }    
  }
  fNtracks = good;
  if (IlcDebugLevelClass()>0){
    Info("RemoveUsed2","\n*****\nNumber of good tracks after shared removal\t%d\n",fNtracks);
  }

  delete []quilcty;
  delete []indexes;
}

void IlcDCHtracker::UnsignClusters() 
{
  //
  // loop over all clusters and unsign them
  //

  for(Int_t layer=0;layer<fSector.GetNLayers();layer++){
    IlcDCHcluster *cl = fSector[layer].GetClusters();
    for (Int_t icl =0;icl< fSector[layer].GetN();icl++)
      cl[icl].Use(-1);
  }
	
  
}



void IlcDCHtracker::SignClusters(TObjArray * arr, Float_t fnumber, Float_t fdensity)
{
  //
  //sign clusters to be "used"
  //
  // snumber and sdensity sign number of sigmas - bellow mean value to be accepted
  // loop over "primaries"
  int nlayers=fSector.GetNLayers();
  
  Float_t sumdens=0;
  Float_t sumdens2=0;
  Float_t sumn   =0;
  Float_t sumn2  =0;
  Float_t sumchi =0;
  Float_t sumchi2 =0;

  Float_t sum    =0;

  Int_t nseed = arr->GetEntriesFast();
  for (Int_t i=0; i<nseed; i++) {
    IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }    
    if (!(pt->IsActive())) continue;
    Float_t dens = pt->GetNumberOfClusters()/Float_t(pt->GetNFoundable());
    Float_t chi2 = pt->GetChi2()/((kUseZCoordinate?2:1)*pt->GetNumberOfClusters()-5);
    if ( (dens>0.1) && (pt->GetNumberOfClusters()>fSector.GetNLayers()*0.4-30)&&chi2<2.5){
      sumdens += dens;
      sumdens2+= dens*dens;
      sumn    += pt->GetNumberOfClusters();
      sumn2   += pt->GetNumberOfClusters()*pt->GetNumberOfClusters();
//       if (chi2>2.5) chi2=2.5;
      sumchi  +=chi2;
      sumchi2 +=chi2*chi2;
      sum++;
    }
  }

  Float_t mdensity = 0.9;
  Float_t meann    = fSector.GetNLayers()*0.4;
  Float_t meanchi  = 1;
  Float_t sdensity = 0.1;
  Float_t smeann    = 10;
  Float_t smeanchi  =0.4;
  
  if (sum>20){
    mdensity = sumdens/sum;
    meann    = sumn/sum;
    meanchi  = sumchi/sum;
    //
    sdensity = sumdens2/sum-mdensity*mdensity;
    if(sdensity<=0) IlcWarning(Form("Negative sdensity:%f, den2:%f, mden:%f",sdensity,
				    sumdens2/sum,mdensity));
    sdensity = TMath::Sqrt(fabs(sdensity));
    //
    smeann   = sumn2/sum-meann*meann;
    if(smeann<=0) IlcWarning(Form("Negative smeann:%f, sumn2:%f, meann:%f",smeann,
				    sumn2/sum,meann));
    smeann   = TMath::Sqrt(fabs(smeann));
    //
    smeanchi = sumchi2/sum - meanchi*meanchi;
    if(smeanchi<=0) IlcWarning(Form("Negative smeanchi:%f, sumn2:%f, meanchi:%f",smeanchi,
				    sumchi2/sum,meanchi));
    smeanchi = TMath::Sqrt(fabs(smeanchi));
  }

  if(IlcDebugLevelClass()>=5){
    IlcDebug(5,Form("meann %f ,smeann %f, meanchi %f,smeanchi %f, meandens %f, smeandens %f",meann,smeann,meanchi,smeanchi,mdensity,sdensity));
  }

  //REMOVE  SHORT DELTAS or tracks going out of sensitive volume of DCH
  //
  for (Int_t i=0; i<nseed; i++) {
    IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }
    if (pt->GetBSigned()) continue;
    if (pt->GetBConstrain()) continue;    

    Bool_t isok =kFALSE;
    if ( (pt->GetNShared()/pt->GetNumberOfClusters()<0.5) &&pt->GetNumberOfClusters()>=fRecPar->GetMinNOfClusetrs())
      isok = kTRUE;
    if ((TMath::Abs(1/pt->GetC())<100.) && (pt->GetNShared()/pt->GetNumberOfClusters()<0.7))
      isok =kTRUE;
    if  (TMath::Abs(pt->GetZ()/pt->GetR())>1.1)
      isok =kTRUE;
    if ( (TMath::Abs(pt->GetSnp())>0.7 && pt->GetD(0,0)>60.))
      isok =kTRUE;
    
    if (isok)     
      for (Int_t i=0; i<nlayers; i++) 
	for(int j=0;j<kMaxInLayer;j++){
	  Int_t index=pt->GetClusterIndex(i,j);
	  if (index<0) continue;
	  IlcDCHcluster *c= pt->GetClusterPointer(i,j);
	  if (!c) continue;
	  //	if(pt->GetClusterChi2(i)<9)
	  c->Use(10);  
	}
  }
  
  
  //
  Float_t bfactor=1;
  Double_t maxchi  = meanchi+2.*smeanchi;
  Double_t mindens = TMath::Max(double(mdensity-sdensity*fdensity*bfactor),0.33*fRecPar->GetMinDensity());
  Double_t minn    = TMath::Max(Int_t(meann-fnumber*smeann*bfactor),fRecPar->GetMinNOfClusetrs());
  
  if(IlcDebugLevelClass()>=5)
    IlcDebug(5,Form("minimum acceptable dens (forr sign clusters) %f, nrecpoints %f, maxchi2=%f",mindens,minn,maxchi));
  
  int nsigned=0;
  int nsignedcl=0;
  int nsignedclunic=0;

  for (Int_t i=0; i<nseed; i++) {
    IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }    
    //if (!(pt->IsActive())) continue;
    if (pt->GetBSigned()) continue;
    if(IlcDebugLevelClass()>=5){
      IlcDebug(5,Form("Signed %i, seed chi2 %f, NofClusters %i, first point %i , last point %i ",
		      i,pt->GetChi2(),pt->GetNumberOfClusters(),pt->GetFirstPoint(),pt->GetLastPoint()));
      IlcDebug(5,Form("dens=%f, fisrtdensity=%f",pt->GetNumberOfClusters()/Float_t(pt->GetNFoundable()),pt->GetDensityFirst(2*fRecPar->GetMinNOfClusetrs())));
    }
    Double_t chi     = pt->GetChi2()/((kUseZCoordinate?2:1)*pt->GetNumberOfClusters()-5);
    if (chi>maxchi) continue;

    Float_t dens = pt->GetNumberOfClusters()/Float_t(pt->GetNFoundable());
   
    //sign only tracks with enoug big density at the beginning
    
    if ((pt->GetDensityFirst(2*fRecPar->GetMinNOfClusetrs())<fRecPar->GetMinDensity()) && pt->GetNumberOfClusters()<fRecPar->GetMinNOfClusetrs()) continue; 
    
    
    
    //    if (pt->fBConstrain) mindens = TMath::Max(mdensity-sdensity*fdensity*bfactor,0.65);
    if ( (pt->GetRemoval()==10) && (pt->GetSnp()>0.8)&&(dens>mindens))
      minn=0;

    if ((dens>mindens && pt->GetNumberOfClusters()>=minn) && chi<maxchi ){
      if(IlcDebugLevelClass()>=5)
	IlcDebug(5,Form("sign cluster forr track"));
      //Int_t noc=pt->GetNumberOfClusters();
      pt->SetBSigned(kTRUE);
      for (Int_t i=0; i<nlayers; i++) 
	for(int j=0;j<kMaxInLayer;j++){
	  Int_t index=pt->GetClusterIndex(i,j);
	  if (index<0) continue;
	  IlcDCHcluster *c= pt->GetClusterPointer(i,j);
	  if (!c) continue;
	  if(!c->IsUsed(10)) nsignedclunic++;
	  c->Use(10);  
	  nsignedcl++;
	}
      if (IlcDCHReconstructor::StreamLevel()>0){
	CookLabel(pt,0.1);
	TTreeSRedirector &cstream = *fDebugStreamer;
	cstream<<"Sign"
	       <<"seed.="<<pt<<"\n"; 
      }
      nsigned++;
    }
  }
  //  gLastCheck = nseed;
  //  arr->Compress();
  if (IlcDebugLevelClass()>3){
    Info("SignClusters",Form("Number of signed tracks %i clusters: %i unic clusters: %i",nsigned,nsignedcl,nsignedclunic));
  }
}


void  IlcDCHtracker::StopNotActive(TObjArray * arr, Int_t layer0, Float_t th0, Float_t th1, Float_t th2) const
{
  // stop not active tracks
  // take th1 as threshold for number of founded to number of foundable on last 10 active layers
  // take th2 as threshold for number of founded to number of foundable on last 20 active layers 
  Int_t nseed = arr->GetEntriesFast();  
  //
  for (Int_t i=0; i<nseed; i++) {
    IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }
    if (!(pt->IsActive())) continue;
    StopNotActive(pt,layer0,th0, th1,th2);
  }
}



void  IlcDCHtracker::StopNotActive(IlcDCHseed * seed, Int_t layer0, Float_t th0, Float_t th1,
 Float_t th2) const
{
  // stop not active tracks
  // take th1 as threshold for number of founded to number of foundable on last 10 active layers
  // take th2 as threshold for number of founded to number of foundable on last 20 active layers 
  Int_t sumgood1  = 0;
  Int_t sumgood2  = 0;
  Int_t foundable = 0;
  Int_t maxindex = seed->GetLastPoint();  //last foundable layer
  if (seed->GetNFoundable()*th0 > seed->GetNumberOfClusters()) {
    seed->Desactivate(10) ;
    return;
  }

  for (Int_t i=layer0; i<maxindex; i++){
    Int_t index = seed->GetClusterIndex(i,0);
    if (index!=-1) foundable++;
    //if (!c) continue;
    if (foundable<=2*fRecPar->GetMinNOfClusetrs()) sumgood1++;
    if (foundable<=3*fRecPar->GetMinNOfClusetrs()) {
      sumgood2++;
    }
    else{ 
      break;
    }        
  }
  if (foundable>=2*fRecPar->GetMinNOfClusetrs()){ 
     if (sumgood1<(th1*2*fRecPar->GetMinNOfClusetrs())) seed->Desactivate(10);
  }
  if (foundable>=3*fRecPar->GetMinNOfClusetrs())
    if (sumgood2<(th2*3*fRecPar->GetMinNOfClusetrs())) seed->Desactivate(10);
}

/*
Int_t IlcDCHtracker::RefitInward(IlcESD *event)
{
  //
  // back propagation of ESD tracks
  //
  //return 0;
  fEvent = event;
  ReadSeeds(event,2);
  fIteration=2;
  //PrepareForProlongation(fSeeds,1);
  PropagateForward2(fSeeds);

  Int_t ntracks=0;
  Int_t nseed = fSeeds->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    IlcDCHseed * seed = (IlcDCHseed*) fSeeds->UncheckedAt(i);
    if (!seed) continue;
    if (seed->GetKinkIndex(0)>0) {
      UpdateKinkQuilctyD(seed);  // update quilcty informations for kinks
    }else{
      seed->PropagateTo(fParam->GetInnerRadius());
    };
    seed->UpdatePoints();
    IlcESDtrack *esd=event->GetTrack(i);
    seed->CookdEdx(0.02,0.6);
    CookLabel(seed,0.1); //For comparison only
    //
    if (IlcDCHReconstructor::StreamLevel()>0 && seed!=0&&esd!=0) {
      TTreeSRedirector &cstream = *fDebugStreamer;
      cstream<<"Crefit"<<
	"Esd.="<<esd<<
	"Track.="<<seed<<
	"\n"; 
    }
    if (seed->GetNumberOfClusters()>=fRecPar->GetMinNOfClusetrs()){
      esd->UpdateTrackParams(seed,IlcESDtrack::kTPCrefit); 
      esd->SetTPCPoints(seed->GetPoints());
      esd->SetTPCPointsF(seed->GetNFoundable());
      Int_t ndedx   = seed->GetNCDEDX(0)+seed->GetNCDEDX(1)+seed->GetNCDEDX(2)+seed->GetNCDEDX(3);
      Float_t sdedx = (seed->GetSDEDX(0)+seed->GetSDEDX(1)+seed->GetSDEDX(2)+seed->GetSDEDX(3))*0.25;
      Float_t dedx  = seed->GetdEdx();
      esd->SetTPCsignal(dedx, sdedx, ndedx);
      //
      // add seed to the esd track in Calib level
      //
      if (IlcDCHReconstructor::StreamLevel()>0){
	IlcDCHseed * seedCopy = new IlcDCHseed(*seed, kTRUE); 
	esd->AddCalibObject(seedCopy);
      }
      ntracks++;
    }
    else{
      //printf("problem\n");
    }
  }
  //FindKinks(fSeeds,event);
  Info("RefitInward","Number of refitted tracks %d",ntracks);
  fEvent =0;
  //WriteTracks();
  return 0;
}
*/
/*
Int_t IlcDCHtracker::PropagateBack(IlcESD *event)
{
  //
  // back propagation of ESD tracks
  //

  fEvent = event;
  fIteration = 1;
  ReadSeeds(event,1);
  PropagateBack(fSeeds); 
  RemoveUsed2(fSeeds,0.4,0.4,fRecPar->GetMinNOfClusetrs());
  //
  Int_t nseed = fSeeds->GetEntriesFast();
  Int_t ntracks=0;
  for (Int_t i=0;i<nseed;i++){
    IlcDCHseed * seed = (IlcDCHseed*) fSeeds->UncheckedAt(i);
    if (!seed) continue;
    if (seed->GetKinkIndex(0)<0)  UpdateKinkQuilctyM(seed);  // update quilcty informations for kinks
    seed->UpdatePoints();
    IlcESDtrack *esd=event->GetTrack(i);
    seed->CookdEdx(0.02,0.6);
    CookLabel(seed,0.1); //For comparison only
    if (seed->GetNumberOfClusters()>=fRecPar->GetMinNOfClusetrs()){
      esd->UpdateTrackParams(seed,IlcESDtrack::kTPCout);
      esd->SetTPCPoints(seed->GetPoints());
      esd->SetTPCPointsF(seed->GetNFoundable());
      Int_t ndedx   = seed->GetNCDEDX(0)+seed->GetNCDEDX(1)+seed->GetNCDEDX(2)+seed->GetNCDEDX(3);
      Float_t sdedx = (seed->GetSDEDX(0)+seed->GetSDEDX(1)+seed->GetSDEDX(2)+seed->GetSDEDX(3))*0.25;
      Float_t dedx  = seed->GetdEdx();
      esd->SetTPCsignal(dedx, sdedx, ndedx);
      ntracks++;
      Int_t eventnumber = event->GetEventNumber();// patch 28 fev 06
      if (IlcDCHReconstructor::StreamLevel()>0){
	(*fDebugStreamer)<<"Cback"<<
	  "Tr0.="<<seed<<
	  "EventNr="<<eventnumber<<
	  "\n"; // patch 28 fev 06   
      }
    }
  }
  //FindKinks(fSeeds,event);
  Info("PropagateBack","Number of back propagated tracks %d",ntracks);
  fEvent =0;
  //WriteTracks();

  return 0;
}
*/

void IlcDCHtracker::DeleteSeeds()
{
  //
  //delete Seeds

  Int_t nseed = fSeeds->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    IlcDCHseed * seed = (IlcDCHseed*)fSeeds->At(i);
    if (seed) delete fSeeds->RemoveAt(i);
  }
  delete fSeeds;

  fSeeds =0;
}

//void IlcDCHtracker::ReadSeeds(IlcESD *event, Int_t direction)
//{
//  //
//  //read seeds from the event
//  int nlayers=fSector.GetNLayers();
//  Int_t nentr=event->GetNumberOfTracks();
//  if (fDebug>0){
//    Info("ReadSeeds", "Number of ESD tracks: %d\n", nentr);
//  }
//  if (fSeeds)
//    DeleteSeeds();
//  if (!fSeeds){
//    fSeeds = new TObjArray(nentr);
//  }
//  UnsignClusters();
//  //  Int_t ntrk=0;
//  for (Int_t i=0; i<nentr; i++) {
//    IlcESDtrack *esd=event->GetTrack(i);
//    ULong_t status=esd->GetStatus();
//    if (!(status&IlcESDtrack::kTPCin)) continue;
//    IlcDCHtrack t(*esd);
//    t.SetNumberOfClusters(0);
//    //    IlcDCHseed *seed = new IlcDCHseed(t,t.GetAlpha());
//    IlcDCHseed *seed = new IlcDCHseed(t/*,t.GetAlpha()* /);
//    for (Int_t ikink=0;ikink<3;ikink++) {
//      Int_t index = esd->GetKinkIndex(ikink);
//      seed->GetKinkIndexes()[ikink] = index;
//      if (index==0) continue;
//      index = TMath::Abs(index);
//      IlcESDkink * kink = fEvent->GetKink(index-1);
//      if (kink&&esd->GetKinkIndex(ikink)<0){
//	if ((status & IlcESDtrack::kTRDrefit) != 0) kink->SetStatus(1,2);
//	if ((status & IlcESDtrack::kVXDout) != 0)   kink->SetStatus(1,0);
//      }
//      if (kink&&esd->GetKinkIndex(ikink)>0){
//	if ((status & IlcESDtrack::kTRDrefit) != 0) kink->SetStatus(1,6);
//	if ((status & IlcESDtrack::kVXDout) != 0)   kink->SetStatus(1,4);
//      }
//
//    }
//    if (((status&IlcESDtrack::kVXDout)==0)&&(direction==1)) seed->ResetCovariance(10.);
//    if ( direction ==2 &&(status & IlcESDtrack::kTRDrefit) == 0 ) seed->ResetCovariance(10.);
//    if ( direction ==2 && ((status & IlcESDtrack::kTPCout) == 0) ) {
//      fSeeds->AddAt(0,i);
//      delete seed;
//      continue;
//    }
//    if ( direction ==2 &&(status & IlcESDtrack::kTRDrefit) > 0 )  {
//      Double_t par0[5],par1[5],alpha,x;
//      esd->GetInnerExternalParameters(alpha,x,par0);
//      esd->GetExternalParameters(x,par1);
//      Double_t delta1 = TMath::Abs(par0[4]-par1[4])/(0.000000001+TMath::Abs(par0[4]+par1[4]));
////       Double_t delta2 = TMath::Abs(par0[3]-par1[3]);
//      Double_t trdchi2=0;
//      if (esd->GetTRDncls()>0) trdchi2 = esd->GetTRDchi2()/esd->GetTRDncls();
//      //reset covariance if suspicious
//      if ( (delta1>0.1)  ||trdchi2>7.) //|| (delta2>0.006)
//	seed->ResetCovariance(10.);
//    }
//
//    if (fDebug>0){
//      Info("ReadSeeds", "Number of accepted ESD tracks: %d\n", fSeeds->GetEntries());
//    }
//
//    seed->SetESD(esd);
//    // sign clusters
//    if (esd->GetKinkIndex(0)<=0){
//      for (Int_t ilayer=0;ilayer<nlayers;ilayer++){
//	Int_t index = seed->GetClusterIndex(ilayer);
//	if (index<0) continue;
//
//	IlcDCHcluster * cl = GetCluster(index);
//	seed->SetClusterPointer(ilayer,cl);
//	if (cl){
//	  if ((index & ClusterIndex::NotAccept)==0){
//	    cl->Use(10);  // accepted cluster
//	  }else{
//	    cl->Use(6);   // close cluster not accepted
//	  }
//	}else{
//	  Info("ReadSeeds","Not found cluster");
//	}
//      }
//    }
//    fSeeds->AddAt(seed,i);
//  }
//}

//_____________________________________________________________________________
void IlcDCHtracker::MakeSeeds3Stereo(TObjArray * arr, Int_t i1, Int_t i2,  Float_t cuts[4]) {
  //-----------------------------------------------------------------
  // This function creates track seeds.
  // SEEDING WITH VERTEX CONSTRAIN 
  //-----------------------------------------------------------------
  // cuts[0]   - fP4 cut
  // cuts[1]   - tan(phi)  cut
  // cuts[2]   - zvertex cut
  // cuts[3]   - fP3 cut
  TStopwatch timer;
  timer.Start();

  Int_t nin0  = 0;
  Int_t nin1  = 0;
  Int_t nin2  = 0;
  Int_t nin   = 0;
  Int_t nout1 = 0;
  Int_t nout2 = 0;

  Double_t x[5], c[15];
  //  Int_t di = i1-i2;
  //
  IlcDCHseed * seed = new IlcDCHseed;
  
  Double_t x1 =GetRlayer(i1);

  Double_t x3g=GetX(), y3g=GetY(), z3g=GetZ();//-50.0;
  Double_t x3=x3g,y3=y3g,z3=z3g;
  Double_t phi1;

  Int_t imiddle = (i2+i1)/2;    //middle pad layer index

  int layind[3]={imiddle,imiddle-1,imiddle+1};
  const IlcDCHLayer* krm[3]={&GetLayer(imiddle),
			     &GetLayer(imiddle-1),
			     &GetLayer(imiddle+1)}; //middle pad -layer

  Int_t inext = (i2+i1)/2-abs(i2-i1);
  const IlcDCHLayer* krn[3]={0,0,0};
  if(inext>0){
    krn[0]=&GetLayer(inext);
    krn[1]=&GetLayer(inext-1);
    krn[2]=&GetLayer(inext+1);
  }
  //

  const IlcDCHLayer& kr1=GetLayer(i1);
  const IlcDCHLayer& kr2=GetLayer(i2);

  //
  // change cut on curvature if it can't reach this layer
  // maximal curvature set to reach it
  Double_t dvertexmax  = TMath::Sqrt((x1-x3)*(x1-x3));
  if (dvertexmax*0.5*cuts[0]>1.1){
    if(fDebug>15) cout<<"Rmin before "<<1/cuts[0]<<" after "<<(dvertexmax*0.5)/1.1<<endl;
    cuts[0] = 1.1/(dvertexmax*0.5);
  }
  Double_t r2min = 1/(cuts[0]*cuts[0]);  //minimal square of radius given by cut
  if(IlcDebugLevelClass()>3) {
    cout<<"Seeding3: i1layer="<<i1<<" i2layer="<<i2<<" im "<<imiddle<<
      " Minimum rmin:"<<TMath::Sqrt(r2min)<<endl;
  }

  // loop over clusters  
  for (Int_t is=0; is < kr1; is++) {
    //
    const IlcDCHcluster *kcl1 = kr1[is];
    if (kcl1->IsUsed(10)) continue;


    //
    for (Int_t js=0; js < kr2; js++) {
      const IlcDCHcluster *kcl2 = kr2[js];
      if (kcl2->IsUsed(10)) continue;	

      //curvature (radius) cut
      if(fDebug>15){ 
	cout<<" trackid "<<kcl1->GetLabel(0)<<" "<<kcl1->GetLabel(1)
	    <<" "<<kcl2->GetLabel(0)<<" "<<kcl2->GetLabel(1)<<" ij "<<is<<" "<<js<<" ll="<<i1<<" "<<i2<<endl;
      }
      //     if (r02<r2min) continue;		
      
      nin0++;	
      //
      // do we have cluster at the middle ?
      UInt_t mindex; 
      IlcDCHcluster * kcm=0;

      double rwire[3]={kr2.GetR(),0.,kr1.GetR()};
      double tangle[3]={kr2.GetTStereoAngle(),0.,kr1.GetTStereoAngle()};
      double phi[3]={kr2.GetWirePhi(kcl2->GetW()),0.,kr1.GetWirePhi(kcl1->GetW())};
      int i0=0;

      int noutfirst=nout1;
      int foundmiddleprev=-1;
      int selectSecond=0;
      double phiselect[3]={-123,-123,-123};
      int wireselect[3]={-1,-1,-1};
      for(int ks=0;ks<6;ks++){
	if(nout1!=noutfirst) break;
	if(ks==3&&foundmiddleprev<0) break;
	if(ks==3){selectSecond=1;}

	int im=ks%3;
	
	double phisamestereo=123.;
	if(tangle[2]*krm[im]->GetTStereoAngle()>0) {
	  phisamestereo=phi[2];
	  i0=1;
	}else if(tangle[0]*krm[im]->GetTStereoAngle()>0){
	  phisamestereo=phi[0];
	  i0=0;
	}else{
	  if(fDebug>15) cout<<"Different stereo angles "<<tangle[0]<<" "<<tangle[2]<<" "<<krm[im]->GetTStereoAngle()<<endl;
	  continue;
	}
	

	if(selectSecond){
	  if(wireselect[im]<0) continue;
	  double phisamestereotry=phisamestereo-Phi_mpi_pi(phiselect[im]-phisamestereo);
	  kcm = krm[im]->FindNearest2(phisamestereotry,0.,0.,fRecPar->GetRoadY()+rwire[2]-rwire[0],1e6,mindex,0);
	  if(!kcm||kcm->GetW()==wireselect[im])
	    kcm = krm[im]->FindNearest2(phisamestereo,0.,0.,fRecPar->GetRoadY()+rwire[2]-rwire[0],1e6,mindex,selectSecond);
	  else
	    phisamestereo=phisamestereotry;
	}else{
	  kcm = krm[im]->FindNearest2(phisamestereo,0.,0.,fRecPar->GetRoadY()+rwire[2]-rwire[0],1e6,mindex,selectSecond);
	}
	if(fDebug>15) cout<<"lay "<<layind[im] 
			  <<" on middle "<<phisamestereo<<" road "<<fRecPar->GetRoadY()+rwire[2]-rwire[0];
	if ((!kcm) || (kcm->IsUsed(10))) {	  
	  if(fDebug>15) cout<<" not found on middles "<<layind[im]<<endl;
	  continue;
	}
	
	if(fDebug>15) cout<<" "<<mindex<<" trackid "<<kcm->GetLabel(0)<<" "<<kcm->GetLabel(1);
	
	wireselect[im]=kcm->GetW();

	int foundmiddle=im;
	foundmiddleprev=foundmiddle;

	nin1++;              

	rwire[1]=krm[foundmiddle]->GetR();
	tangle[1]=krm[foundmiddle]->GetTStereoAngle();
	phi[1]=krm[foundmiddle]->GetWirePhi(kcm->GetW());

	phiselect[im]=phi[1];
	
	phi1=(phi[1-i0]+phi[1-i0+1])/2.+((fabs(phi[1-i0]-phi[1-i0+1])>TMath::Pi()/2)?(TMath::Pi()):0.);
	double zz=TMath::Tan(phi1-phi[1-i0])*rwire[1-i0]/tangle[1-i0];
	
	if(fDebug>15) cout<<" zz "<<zz<<" phis "<<phi[1-i0]<<" "<<phi[1-i0+1]<<endl;
	if(fabs(zz)>kr2.GetMaxZ()+fRecPar->GetSafetyZ()) continue;


	
	
	const IlcDCHcluster *clseed[3]={kcl2,kcm,kcl1};
	const UInt_t index[3]={kr2.GetIndex(js),krm[im]->GetIndex(mindex),kr1.GetIndex(is)};
	const int ilay[3]={i2,layind[im],i1};
	
	for(int i=0;i<15;i++)c[i]=0.;
	for(int i=0;i<5;i++) c[(i+1)*(i+2)/2-1]=1;
	c[0]=fRecPar->GetSeedBeamSigma2();
	c[2]=fRecPar->GetSeedBeamSigmaZ2();
	c[14]=1.e-3;
	
	
	double mes[3]={kcl2->GetImP(),kcm->GetImP(),kcl1->GetImP()};
	
	for(int lr=0;lr<8;lr++){
	  int leftr[3]={(lr%2?1:-1),((lr/2)%2?1:-1),((lr/4)%2?1:-1)};
	  if(fDebug>15)
	    cout<<"Left-Right "<<leftr[2]<<" "<<leftr[1]<<" "<<leftr[0]
		<<" trackid "<<clseed[2]->GetLabel(0)<<" "<<clseed[2]->GetLabel(1)<<","
		<<clseed[1]->GetLabel(0)<<" "<<clseed[1]->GetLabel(1)<<","
		<<clseed[0]->GetLabel(0)<<" "<<clseed[0]->GetLabel(1)<<" "<<endl;
	  
// 	  if(clseed[0]->GetLabel(0)!=clseed[1]->GetLabel(0)||
// 	     clseed[0]->GetLabel(0)!=clseed[2]->GetLabel(0)) continue;
	  
	  x[0]=y3;
	  x[1]=z3;
	  double dphi0=((phi[i0]+leftr[i0]*mes[i0]/rwire[i0]+TMath::ATan(zz*tangle[i0]/rwire[i0]))-
			(phi[i0+1]+leftr[i0+1]*mes[i0+1]/rwire[i0+1]+
			 TMath::ATan((zz+(rwire[i0+1]-rwire[i0])*x[3])*tangle[i0+1]/rwire[i0+1]))
			);
	  dphi0=Phi_mpi_pi(dphi0);
	  x[2]=-dphi0*rwire[i0]/(rwire[i0+1]-rwire[i0]);
	  //tgl to sin
	  x[2]=x[2]/TMath::Hypot(1.,x[2]);
	  x[4]=2*x[2]/rwire[i0+1];
	  double ratzw=TMath::Hypot(rwire[1-i0],zz*tangle[1-i0]);
	  //	  x[3]=(zz-z3)/TMath::Hypot(rwire[1-i0],zz*tangle[1-i0]);
	  x[3]=(zz-z3)/(TMath::ASin(x[2])/x[2]*ratzw);
	  //x[3]=-(zz-z3)/((TMath::Pi()-TMath::ASin(fabs(x[2])))/fabs(x[2])*ratzw);
	  

	  double xyz[3];			

	  IlcDCHseed *track=new(seed) IlcDCHseed(x3+1e-4, phi1, x, c, index[2]);
	  if(fDebug>15){
	    cout<<"yz "<<x[0]<<" "<<x[1]<<" tgl "<<x[3]<<" r "<<x[4]<<" sphi "<<x[2]
		<<" dphi0 "<<dphi0<<" r1-r0 "<<rwire[i0+1]<<" "<<rwire[i0]<<endl;
	    track->Print();
	  }
	  double chi2=1e30;
	  double propagateOK=true;
	  int kk;
	  for(kk=0;kk<10;kk++){
	    //	    if(!track->IlcExternalTrackParam::PropagateTo(0.001,IlcTracker::GetBz())){propagateOK=false;break;}
	    track->IlcExternalTrackParam::PropagateTo(0.001,IlcTracker::GetBz());

	    if(fDebug>15) track->Print("params");
	    track->Modify(1000);//ResetCovariance(100);
	    //	    if(kk>0) if(!track->IlcExternalTrackParam::Update(x,c)){propagateOK=false;break;}	
	    if(kk<3) {
	      	c[0]=fRecPar->GetSeedBeamSigma2();
		c[2]=fRecPar->GetSeedBeamSigma2();
	    }

	    if(kk>0) if(!track->IlcExternalTrackParam::UpdateWithWire(0.,0.,0.,IlcTracker::GetBz(),
								      x,c,kTRUE,kTRUE)){propagateOK=false;break;}		
	    chi2=0.;
	    for(int i=0;i<3;i++){
	      track->GetXYZAtStereo(rwire[i],tangle[i],IlcTracker::GetBz(),xyz);
	      if(!track->PropagateTo(TMath::Hypot(xyz[0],xyz[1]),0.,1.e10)&&i<2){propagateOK=false;break;}

	      double dphi=track->GetPhi()-phi[i]-TMath::ATan(track->GetZ()/rwire[i]*tangle[i]);
	      dphi=Phi_mpi_pi(dphi);
	      track->SetErrorY2(clseed[i]->GetSigmaImP2());
	      bool calcChi2=true;
	      double chi2=track->GetPredictedChi2(clseed[i],rwire[i],phi[i],tangle[i],kUseZCoordinate);
	      if(!track->Update(clseed[i],chi2,0,
				rwire[i],phi[i],tangle[i],kUseZCoordinate,((dphi*leftr[i])<0?-1:1),calcChi2))
		{propagateOK=false;break;}
	      if(fabs(track->GetTgl())>IlcDCHReconstructor::GetCtgRange()||
		 fabs(track->GetPt())<kSmallPt){propagateOK=false;break;}
	      //    track->Print();
	    }
	    chi2=track->GetChi2();
	    if(fDebug>15) cout<<"chi2 "<<chi2<<" pOK "<<propagateOK<<" tz "<<track->GetZ()
			      <<" "<<kr1.GetMaxZ()+fRecPar->GetSafetyZ()<<endl;
	    if(fDebug>15) track->Print("params");
	    if(!propagateOK) break;
	    if(kk>=3&&chi2<1e-5) break;
	    if(fabs(track->GetZ())>2*(kr1.GetMaxZ()+fRecPar->GetSafetyZ())) break;
	  }
	  //	  std::cout<<"seedmake "<<kk<<" "<<chi2<<" "<<track->GetZ()<<" "<<propagateOK<<std::endl;
       
	  if(!propagateOK) continue;
	  if(chi2>1.) continue;
	  if(fabs(track->GetZ())>1.5*(kr1.GetMaxZ()+fRecPar->GetSafetyZ())) continue;


	  int foundnext=-1;
	  if(krn[0]){
	    //	    double r=krn[0].GetRAtZ(track->GetZ());
	    for(int i=0;i<3;i++){
	      track->GetXYZAtStereo(krn[i]->GetR(),krn[i]->GetTStereoAngle(),IlcTracker::GetBz(),xyz);
	      double phiat=TMath::ATan2(xyz[1],xyz[0]);
	      if(fDebug>15) cout<<"on next "<<Phi_0_2pi(phiat)<<" road "<<fRecPar->GetRoadY()<<" z= "<<xyz[2]<<" ";
	      const IlcDCHcluster *kcm=krn[i]->FindNearest2(phiat,track->GetSnp(),xyz[2],fRecPar->GetRoadY(),1e6,mindex);
	      if ((!kcm) || (kcm->IsUsed(10))) {	  
		if(fDebug>15) cout<<"not found on next "<<inext+(i/2-(i%2))<<endl;
	      }else{
		double chi2=track->GetPredictedChi2(kcm,krn[i]->GetR(),krn[i]->GetWirePhi(kcm->GetW()),
						    krn[i]->GetTStereoAngle(),kUseZCoordinate);
		if(fDebug>15) cout<<" "<<mindex<<" trackid "<<kcm->GetLabel(0)<<" "<<kcm->GetLabel(1)<<" chi2 "<<chi2<<endl;
		if(chi2/(kUseZCoordinate?2:1)<16.){
		  foundnext=i;
		  break;
		}
	      }
	    }
	    if(foundnext<0) continue;
	  }
	  
	  nin2++;

	  if(track->GetTgl()<0) track->ChangeDirection();

	  track->Modify(1000.);
	  double ratz=TMath::Hypot(track->GetZ()*tangle[2],rwire[2]);
	  track->PropagateTo(ratz);
	  if(fDebug>15){
	    track->Print();
	    cout<<"double rwire[3]={"<<kr2.GetR()<<","<<krm[foundmiddle]->GetR()<<","<<kr1.GetR()<<"};"<<endl;
	    cout<<"double angle[3]={"<<kr2.GetTStereoAngle()<<","<<krm[foundmiddle]->GetTStereoAngle()<<","
		<<kr1.GetTStereoAngle()<<"};"<<endl;
	    cout<<"double phi[3]={"<<kr2.GetWirePhi(kcl2->GetW())<<","<<krm[foundmiddle]->GetWirePhi(kcm->GetW())
		<<","<<kr1.GetWirePhi(kcl1->GetW())<<"};"<<endl;
	    cout<<"double mes[3][2]={{"<<kcl2->GetImP()<<","<<kcl2->GetZLocal()<<"},"<<endl
		<<"{"<<kcm->GetImP()<<","<<kcm->GetZLocal()<<"},"<<endl
		<<"{"<<kcl1->GetImP()<<","<<kcl1->GetZLocal()<<"}};"<<endl;
	    cout<<"double covm[3][3]={{"<<kcl2->GetSigmaImP2()<<",0,"<<kcl2->GetSigmaZ2()<<"},"<<endl
		<<"{"<<kcm->GetSigmaImP2()<<",0,"<<kcm->GetSigmaZ2()<<"},"<<endl
		<<"{"<<kcl1->GetSigmaImP2()<<",0,"<<kcl1->GetSigmaZ2()<<"}};"<<endl;
	  }
	
	  track->SetIsSeeding(kTRUE);
	  track->SetSeed1(i1);
	  track->SetSeed2(i2);
	  track->SetSeedType(3);
	  track->SetNumberOfClusters(0);
	  for(int i=2;i>=0;i--){				
	    track->SetClusterIndex(ilay[i],index[i],0);
	    track->SetClusterPointer(ilay[i],(IlcDCHcluster*)clseed[i],0);
	  }
	  fIteration=1;
	  FollowProlongation(*track, imiddle,1,0);
	  fIteration=0;
	  if(fDebug>21) track->Dump();
	  Int_t foundable,found,shared;
	  track->GetClusterStatistic(imiddle,i1, found, foundable, shared, kTRUE);
	  if(fDebug>15) 
	    cout<<"llllllll "<<found<<" "<<foundable<<" "<<shared
		<<" "<<(track->GetSigmaY2()+0*track->GetSigmaZ2())
		<<" "<<fRecPar->GetCutTrackSigma2()*5<<endl;
	  
	  if ((found<fRecPar->GetMinDensity()*foundable)  || shared>0.5*found || (track->GetSigmaY2()+0*track->GetSigmaZ2())>fRecPar->GetCutTrackSigma2()*5){
	    //            cout<<"caput 1\n";
	    seed->Reset();
	    seed->~IlcDCHseed();
	    continue;
	  }

	  if(track->GetTgl()<0) track->ChangeDirection();

	  nin++;
	  fIteration=1;
	  FollowProlongation(*track, i2,1,1);
	  fIteration=0;
	  
	  
	  
	  if(fDebug>15) 
	    cout<<"mmmmmmmm "<<track->GetNumberOfClusters()
		<<" "<<i1<<" "<<i2<<" "<<fRecPar->GetMinDensity()
		<<" "<<track->GetNFoundable()
		<<" "<<track->GetNShared()<<" chi2 "<<track->GetChi2()<<" "<<chi2
		<<" l "<<clseed[2]->GetLabel(0)<<" "<<clseed[2]->GetLabel(1)<<","
		<<clseed[1]->GetLabel(0)<<" "<<clseed[1]->GetLabel(1)<<","
		<<clseed[0]->GetLabel(0)<<" "<<clseed[0]->GetLabel(1)<<" "<<endl;
	  if(track->GetNumberOfClusters()>7&&track->GetChi2()/((kUseZCoordinate?2:1)*track->GetNumberOfClusters()-5)>4){
	    if(fDebug>15) 
	      cout<<"Try to refit segment"<<endl;

	    track->Modify(1000);
	    fIteration=1;
	    FollowBackProlongation(*track, i1);
	    track->Modify(1000);
	    FollowProlongation(*track, i2,1,1);
	    fIteration=0;
	    if(fDebug>15)
	      cout<<"m2m2m2m2m2m2m2m2 "<<track->GetNumberOfClusters()
		  <<" "<<i1<<" "<<i2<<" "<<fRecPar->GetMinDensity()
		  <<" "<<track->GetNFoundable()
		  <<" "<<track->GetNShared()<<" chi2 "<<track->GetChi2()<<" "<<chi2
		  <<" l "<<clseed[2]->GetLabel(0)<<" "<<clseed[2]->GetLabel(1)<<","
		  <<clseed[1]->GetLabel(0)<<" "<<clseed[1]->GetLabel(1)<<","
		  <<clseed[0]->GetLabel(0)<<" "<<clseed[0]->GetLabel(1)<<" "<<endl;
	    
	  }
      
	  if (track->GetNumberOfClusters()<(i1-i2)*fRecPar->GetMinDensity() || 
	      track->GetNumberOfClusters() < track->GetNFoundable()*fRecPar->GetMinDensity() || 
	      track->GetNShared()>0.4*track->GetNumberOfClusters()//){ 
	      ||track->GetNumberOfClusters()<6
	      ||track->GetChi2()/((kUseZCoordinate?2:1)*track->GetNumberOfClusters()-5)>4) {
	    seed->Reset();
	    seed->~IlcDCHseed();
	    continue;
	  }
	  nout1++;

	  //Int_t rc = 1;
	  track->SetBConstrain(1);
	  track->SetLastPoint(i1);  // first cluster in track position
	  track->SetFirstPoint(track->GetLastPoint());

	  if(foundnext>=0){
	    int nclustbefore=track->GetNumberOfClusters();
	    int nfoundablebefore=track->GetNFoundable();
	    int i0=inext+(foundnext/2-(foundnext%2));
	    FollowProlongation(*track, i0,1,1);
	    if(fDebug>15) 
	      cout<<"nnnnnnnn "<<track->GetNumberOfClusters()
		  <<" "<<i1<<" "<<i2<<" "<<i0<<" "<<fRecPar->GetMinDensity()
		  <<" "<<track->GetNFoundable()
		  <<" "<<track->GetNShared()<<" chi2 "<<track->GetChi2()<<" "<<chi2
		  <<" l "<<clseed[2]->GetLabel(0)<<" "<<clseed[2]->GetLabel(1)<<","
		  <<clseed[1]->GetLabel(0)<<" "<<clseed[1]->GetLabel(1)<<","
		  <<clseed[0]->GetLabel(0)<<" "<<clseed[0]->GetLabel(1)<<" "<<endl;
	    if (track->GetNumberOfClusters()<(i1-i0)*fRecPar->GetMinDensity() || 
		track->GetNShared()>0.4*track->GetNumberOfClusters()
		||(track->GetNumberOfClusters()-nclustbefore)<=(track->GetNFoundable()-nfoundablebefore)*fRecPar->GetMinDensity()
		||(track->GetChi2()/((kUseZCoordinate?2:1)*track->GetNumberOfClusters()-5)>4)) {
	      seed->Reset();
	      seed->~IlcDCHseed();
	      continue;
	    }
	  }
	  // Z VERTEX CONDITION
	  Double_t zv, bz=GetBz();
	  if ( !track->GetZAt(0.,bz,zv) ) continue;
	  if (TMath::Abs(zv-z3)>cuts[2]) {
	    FollowProlongation(*track, TMath::Max(i2-20,0),1,1);
	    if ( !track->GetZAt(0.,bz,zv) ) continue;
	    if (TMath::Abs(zv-z3)>cuts[2]){
	      FollowProlongation(*track, TMath::Max(i2-40,0),1,1);
	      if ( !track->GetZAt(0.,bz,zv) ) continue;
	      if (TMath::Abs(zv-z3)>cuts[2] &&(track->GetNumberOfClusters() > track->GetNFoundable()*fRecPar->GetMinDensity())){
		// make seed without constrain
		IlcDCHseed * track2 = MakeSeed(track,0.2,0.5,1.);
		if(track2->GetTgl()<0) track2->ChangeDirection();
		track2->SetNumberOfClusters(0);
		FollowProlongation(*track2, i2,1);
		if(track2->GetTgl()<0) track2->ChangeDirection();
		track2->SetBConstrain(kFALSE);
		track2->SetSeedType(1);
		track2->SetRemoval(0);
		arr->AddLast(track2); 
		if(fDebug>15) 
		  cout<<"Accept seed not constrained id="<<arr->GetEntries()<<endl;
		seed->Reset();
		seed->~IlcDCHseed();
		continue;		
	      }
	      else{
		seed->Reset();
		seed->~IlcDCHseed();
		continue;
		
	      }
	    }
	  }
	  if(track->GetTgl()<0) track->ChangeDirection();

	  track->SetSeedType(0);
	  track->SetRemoval(0);
	  arr->AddLast(track); 

	  if (IlcDCHReconstructor::StreamLevel()>0){
	    CookLabel(track,0.1);
	    TTreeSRedirector &cstream = *fDebugStreamer;
	    cstream<<"Seed"
		   <<"seed.="<<track<<"\n"; 
	  }
	  if(fDebug>15) 
	    cout<<"Accept seed id="<<arr->GetEntries()-1<<endl;
	  seed = new IlcDCHseed; 	
	  nout2++;
	}
      }
    }
  }
  if (IlcDebugLevelClass()>3){
    timer.Print();
    Info("MakeSeeds3Stereo","N passed z sel&Rmin, z cut with dphi cor, middle check, density, nout1, nout2 \nSeeding statistic:\t%d\t%d\t%d\t%d\t%d\t%d",
	 nin0,nin1,nin2,nin,nout1,nout2);
  }
  delete seed;
 
}

void IlcDCHtracker::MakeSeeds3StereoLine(TObjArray * arr, Int_t i1, Int_t i2,  Float_t cuts[4]) {
  //-----------------------------------------------------------------
  // This function creates track seeds.
  // SEEDING WITH LINE and momentum predefined
  //-----------------------------------------------------------------
  // cuts[0]   - fP4 cut
  // cuts[1]   - tan(phi)  cut
  // cuts[2]   - zvertex cut
  // cuts[3]   - fP3 cut
  TStopwatch timer;
  timer.Start();

  Int_t nin0  = 0;
  Int_t nin1  = 0;
  Int_t nin2  = 0;
  Int_t nin   = 0;
  Int_t nout1 = 0;
  Int_t nout2 = 0;

  Double_t x[5], c[15];
  //  Int_t di = i1-i2;
  //
  IlcDCHseed * seed = new IlcDCHseed;
  
  Double_t x1 =GetRlayer(i1);

  Double_t x3g=GetX();//, y3g=GetY(), z3g=GetZ();//-50.0;
  Double_t x3=x3g;//,y3=y3g,z3=z3g;
  Double_t phi1,phi2;

  Int_t imiddle = (i2+i1)/2;    //middle pad layer index

  int layind[4]={0,+2,-1,1};
  const IlcDCHLayer* krm[4]={&GetLayer(imiddle),
			     &GetLayer(imiddle+2),
			     &GetLayer(imiddle-1),
			     &GetLayer(imiddle+1)}; //middle pad -layer

  Int_t inext = (i2+i1)/2-abs(i2-i1);
//  const IlcDCHLayer* krn[3]={0,0,0};
//  if(inext>0){
//    krn[0]=&GetLayer(inext);
//    krn[1]=&GetLayer(inext-1);
//    krn[2]=&GetLayer(inext+1);
//  }
  //

  const IlcDCHLayer& kr1=GetLayer(i1);
  const IlcDCHLayer& kr2=GetLayer(i2);
  const IlcDCHLayer& kr1n=GetLayer(i1-1);
  const IlcDCHLayer& kr2n=GetLayer(i2+1);

  //
  // change cut on curvature if it can't reach this layer
  // maximal curvature set to reach it
  Double_t dvertexmax  = TMath::Sqrt((x1-x3)*(x1-x3));
  if (dvertexmax*0.5*cuts[0]>1.1){
    if(fDebug>15) cout<<"Rmin before "<<1/cuts[0]<<" after "<<(dvertexmax*0.5)/1.1<<endl;
    cuts[0] = 1.1/(dvertexmax*0.5);
  }
  Double_t r2min = 1/(cuts[0]*cuts[0]);  //minimal square of radius given by cut
  if(IlcDebugLevelClass()>3) {
    cout<<"Seeding3: i1layer="<<i1<<" i2layer="<<i2<<" im "<<imiddle<<
      " Minimum rmin:"<<TMath::Sqrt(r2min)<<endl;
  }

  double rwire[6]={kr1.GetR(),kr1n.GetR(),kr2.GetR(),kr2n.GetR()};
  double tangle[6]={kr1.GetTStereoAngle(),kr1n.GetTStereoAngle(),kr2.GetTStereoAngle(),kr2n.GetTStereoAngle()};
  double phi[6];
  // loop over clusters  
  for (Int_t is=0; is < kr1; is++) {
    //
    const IlcDCHcluster *kcl1 = kr1[is];
    if (kcl1->IsUsed(10)) continue;
    phi[0]=kr1.GetWirePhi(kcl1->GetW());

    for (Int_t isn=0; isn < kr1n; isn++) {
      const IlcDCHcluster *kcl1n = kr1n[isn];
      if (kcl1n->IsUsed(10)) continue;
      phi[1]=kr1n.GetWirePhi(kcl1n->GetW());

      phi1=(phi[0]+phi[1])/2.+((fabs(phi[1]-phi[0])>TMath::Pi()/2)?(TMath::Pi()):0.);
      double zz1=TMath::Tan(phi1-phi[0])*rwire[0]/tangle[0];
	
      if(fabs(zz1)>kr1.GetMaxZ()+fRecPar->GetSafetyZ()) continue;

      
      //
      for (Int_t js=0; js < kr2; js++) {
	const IlcDCHcluster *kcl2 = kr2[js];
	if (kcl2->IsUsed(10)) continue;	
	phi[2]=kr2.GetWirePhi(kcl2->GetW());

	for (Int_t jsn=0; jsn < kr2n; jsn++) {
	  const IlcDCHcluster *kcl2n = kr2n[jsn];
	  if (kcl2n->IsUsed(10)) continue;
	  phi[3]=kr2n.GetWirePhi(kcl2n->GetW());
	  
	  phi2=(phi[2]+phi[3])/2.+((fabs(phi[2]-phi[3])>TMath::Pi()/2)?(TMath::Pi()):0.);
	  double zz2=TMath::Tan(phi1-phi[2])*rwire[2]/tangle[2];
	  
	  if(fabs(zz2)>kr2.GetMaxZ()+fRecPar->GetSafetyZ()) continue;

	  //curvature (radius) cut
	  if(fDebug>15){ 
	    cout<<" phiz1 "<<phi1<<" "<<zz1<<" phiz2 "<<phi2<<" "<<zz2
		<<" "<<" iinjjn "<<is<<" "<<isn<<" "<<js<<" "<<jsn<<" "<<" ll="<<i1<<" "<<i2
		<<" phi "<<phi[0]<<" "<<phi[1]<<" "<<phi[2]<<" "<<phi[3]<<endl;
	  }
      
	  nin0++;	

	  //
	  // do we have cluster at the middle ?

	  IlcDCHcluster * kcm=0;
	  /*
	  int i0=0;

	  int noutfirst=nout1;*/
//	  int foundmiddleprev=-1;
	  int selectSecond=0;
	  for(int ks=0;ks<1;ks++){
	    /*	    if(nout1!=noutfirst) break;
	    if(ks==1&&foundmiddleprev<0) break;
	    if(ks==2){selectSecond=1; foundmiddleprev=-1;}
	    */	    
	    //int foundmiddle=-1;
	    int nfoundm[2]={0,0};
	    UInt_t mindex;
	    for(int i=0;i<4;i++){
	      double rmid=krm[i]->GetR();
	      double phisamestereo=(i<2?(phi[3]*(rmid-rwire[0])+phi[0]*(rwire[3]-rmid))/(rwire[3]-rwire[0]):
				    (phi[2]*(rmid-rwire[1])+phi[1]*(rwire[2]-rmid))/(rwire[2]-rwire[1]));

	      if(fDebug>15) cout<<"lay "<<imiddle+layind[i] 
				<<" on middle "<<phisamestereo<<" road "<<krm[i]->GetDeltaPhi()*krm[i]->GetR();
	      
	      kcm = krm[i]->FindNearest2(phisamestereo,0.,0.,krm[i]->GetDeltaPhi()*rmid*300,1e6,mindex,selectSecond);
	      if ((!kcm) || (kcm->IsUsed(10))) {	  
		if(fDebug>15) cout<<" not found"<<endl;
	      }else{
		if(fDebug>15) cout<<" index "<<mindex<<endl;
		nfoundm[i/2]++;
	      }
	    }
	    if(nfoundm[0]==0||nfoundm[1]==0)
	      continue;

	    nin1++;              
	
	
	    UInt_t index=kr1.GetIndex(is);

	    /*
	    rwire[1]=krm[foundmiddle]->GetR();
	    tangle[1]=krm[foundmiddle]->GetTStereoAngle();
	    phi[1]=krm[foundmiddle]->GetWirePhi(kcm->GetW());
	    
	    const IlcDCHcluster *clseed[3]={kcl2,kcm,kcl1};
	    */
	    for(int i=0;i<15;i++)c[i]=0.;
	    for(int i=0;i<5;i++) c[(i+1)*(i+2)/2-1]=1;
	    c[0]=kcl1->GetSigmaImP2()*16;
	    c[2]=c[0]/tangle[0]/tangle[0];
	    c[5]=c[0]*2/(rwire[2]-rwire[0])/(rwire[2]-rwire[0]);
	    c[9]=c[2]*2/(rwire[2]-rwire[0])/(rwire[2]-rwire[0]);
	    c[14]=1.e-2/30./30.;
	
	    
	    //	    if(fabs(zz)>kr2.GetMaxZ()+fRecPar->GetSafetyZ()) continue;
	    
	    double mes[4]={kcl1->GetImP(),kcl1n->GetImP(),kcl2->GetImP(),kcl2n->GetImP()};
	
	    //over all left-rigth
	    for(int lr=0;lr<16;lr++){
	      int leftr[4]={(lr%2?1:-1),((lr>>1)%2?1:-1),((lr>>2)%2?1:-1),((lr>>3)%2?1:-1)};
	      if(fDebug>15)
		cout<<"Left-Right "<<leftr[0]<<" "<<leftr[1]<<" "<<leftr[2]<<" "<<leftr[3]<<endl;
	  
	      x[0]=0;
	      x[1]=zz1;
	      double dphi0=phi[0]+leftr[0]*mes[0]/rwire[0]-(phi[2]+leftr[2]*mes[2]/rwire[2]);
	      dphi0=Phi_mpi_pi(dphi0);
	      x[2]=-dphi0*rwire[0]/(rwire[0]-rwire[2]);
	      //tgl to sin
	      x[2]=x[2]/TMath::Hypot(1.,x[2]);
	      x[3]=(zz1-zz2)/TMath::Hypot(rwire[0]-rwire[2],dphi0*rwire[2]);
	      x[4]=TMath::Sign(1/30.*sqrt(1.+x[2]*x[2]),x[2]);
	  
//	      double xyz[3];
	      double ratz=TMath::Hypot(zz1*tangle[0],rwire[0]);
	      
	      x[2]=-x[2]+ratz*x[4];
	      IlcDCHseed *track=new(seed) IlcDCHseed(ratz, phi1, x, c, index);
	      if(fDebug>15){
		cout<<"yz "<<x[0]<<" "<<x[1]<<" tgl "<<x[3]<<" r "<<x[4]<<" sphi "<<x[2]
		    <<" dphi0 "<<dphi0<<" r0-r1 "<<rwire[0]<<" "<<rwire[2]<<endl;
		track->Print();
	      }

	      /*
	      int foundnext=-1;
	      UInt_t mindex;
	      if(krn[0]){
		//	    double r=krn[0].GetRAtZ(track->GetZ());
		for(int i=0;i<3;i++){
		  track->GetXYZAtStereo(krn[i]->GetR(),krn[i]->GetTStereoAngle(),IlcTracker::GetBz(),xyz);
		  double phiat=TMath::ATan2(xyz[1],xyz[0]);
		  if(fDebug>15) cout<<"on next "<<Phi_0_2pi(phiat)<<" road "<<fRecPar->GetRoadY()<<" z= "<<xyz[2]<<" ";
		  const IlcDCHcluster *kcm=krn[i]->FindNearest2(phiat,track->GetSnp(),xyz[2],fRecPar->GetRoadY(),1e6,mindex);
		  if ((!kcm) || (kcm->IsUsed(10))) {	  
		    if(fDebug>15) cout<<"not found on next "<<inext+(i/2-(i%2))<<endl;
		  }else{
		    double chi2=track->GetPredictedChi2(kcm,krn[i]->GetR(),krn[i]->GetWirePhi(kcm->GetW()),
							krn[i]->GetTStereoAngle(),kUseZCoordinate);
		    if(fDebug>15) cout<<" "<<mindex<<" trackid "<<kcm->GetLabel(0)<<" "<<kcm->GetLabel(1)<<" chi2 "<<chi2<<endl;
		    if(chi2/(kUseZCoordinate?2:1)<16.){
		      foundnext=i;
		      break;
		    }
		  }
		}
		if(foundnext<0) continue;
	      }
	      */
	      nin2++;
	      
	      if(track->GetTgl()<0) track->ChangeDirection();

	      //	      track->ResetCovariance(1000.);
	      //double ratz=TMath::Hypot(track->GetZ()*tangle[0],rwire[0]);
	      //track->PropagateTo(ratz);
	      if(fDebug>15){
		track->Print();
		int foundmiddle=0;
		cout<<"double rwire[3]={"<<kr2.GetR()<<","<<krm[foundmiddle]->GetR()<<","<<kr1.GetR()<<"};"<<endl;
		cout<<"double angle[3]={"<<kr2.GetTStereoAngle()<<","<<krm[foundmiddle]->GetTStereoAngle()<<","
		    <<kr1.GetTStereoAngle()<<"};"<<endl;
		cout<<"double phi[3]={"<<kr2.GetWirePhi(kcl2->GetW())
		  //		    <<","<<krm[foundmiddle]->GetWirePhi(kcm->GetW())
		    <<","<<kr1.GetWirePhi(kcl1->GetW())<<"};"<<endl;
		cout<<"double mes[3][2]={{"<<kcl2->GetImP()<<","<<kcl2->GetZLocal()<<"},"<<endl
		  //		    <<"{"<<kcm->GetImP()<<","<<kcm->GetZLocal()<<"},"<<endl
		    <<"{"<<kcl1->GetImP()<<","<<kcl1->GetZLocal()<<"}};"<<endl;
		cout<<"double covm[3][3]={{"<<kcl2->GetSigmaImP2()<<",0,"<<kcl2->GetSigmaZ2()<<"},"<<endl
		  //		    <<"{"<<kcm->GetSigmaImP2()<<",0,"<<kcm->GetSigmaZ2()<<"},"<<endl
		    <<"{"<<kcl1->GetSigmaImP2()<<",0,"<<kcl1->GetSigmaZ2()<<"}};"<<endl;
	      }
	
	      track->SetIsSeeding(kTRUE);
	      track->SetSeed1(i1);
	      track->SetSeed2(i2);
	      track->SetSeedType(3);
	      track->SetNumberOfClusters(0);
	      track->SetClusterIndex(i1,index,0);
	      track->SetClusterPointer(i1,(IlcDCHcluster*)kcl1,0);
	  
	      FollowProlongation(*track, imiddle,1,0);
	      if(fDebug>21) track->Dump();
	      Int_t foundable,found,shared;
	      track->GetClusterStatistic(imiddle,i1, found, foundable, shared, kTRUE);
	      if(fDebug>15) 
		cout<<"llllllll "<<found<<" "<<foundable<<" "<<shared
		    <<" "<<(track->GetSigmaY2()+0*track->GetSigmaZ2())
		    <<" "<<fRecPar->GetCutTrackSigma2()*5<<endl;
	      
	      if ((found<fRecPar->GetMinDensity()*foundable)  || shared>0.5*found || (track->GetSigmaY2()+0*track->GetSigmaZ2())>fRecPar->GetCutTrackSigma2()*5){
		//            cout<<"caput 1\n";
		seed->Reset();
		seed->~IlcDCHseed();
		continue;
	      }

	      if(track->GetTgl()<0) track->ChangeDirection();
	      
	      nin++;
	      FollowProlongation(*track, i2,1,1);
	  
	  
	      
	      if(fDebug>15) 
		cout<<"mmmmmmmm "<<track->GetNumberOfClusters()
		    <<" "<<i1<<" "<<i2<<" "<<fRecPar->GetMinDensity()
		    <<" "<<track->GetNFoundable()
		    <<" "<<track->GetNShared()<<" chi2 "<<track->GetChi2()<<endl;
	      //		    <<" l "<<clseed[2]->GetLabel(0)<<" "<<clseed[2]->GetLabel(1)<<","
	      //	    <<clseed[1]->GetLabel(0)<<" "<<clseed[1]->GetLabel(1)<<","
	      //	    <<clseed[0]->GetLabel(0)<<" "<<clseed[0]->GetLabel(1)<<" "<<endl;
	      
	      if (track->GetNumberOfClusters()<(i1-i2)*fRecPar->GetMinDensity() || 
		  track->GetNumberOfClusters() < track->GetNFoundable()*fRecPar->GetMinDensity() || 
		  track->GetNShared()>0.4*track->GetNumberOfClusters()//){ 
		  ||track->GetNumberOfClusters()<6
		  ||track->GetChi2()/((kUseZCoordinate?2:1)*track->GetNumberOfClusters()-5)>4) {
		seed->Reset();
		seed->~IlcDCHseed();
		continue;
	      }
	      nout1++;

	      //Int_t rc = 1;
	      track->SetBConstrain(1);
	      track->SetLastPoint(i1);  // first cluster in track position
	      track->SetFirstPoint(track->GetLastPoint());

	      int foundnext=-1;
	      if(foundnext>=0){
		int nclustbefore=track->GetNumberOfClusters();
		int nfoundablebefore=track->GetNFoundable();
		int i0=inext+(foundnext/2-(foundnext%2));
		FollowProlongation(*track, i0,1,1);
		if(fDebug>15) 
		  cout<<"nnnnnnnn "<<track->GetNumberOfClusters()
		      <<" "<<i1<<" "<<i2<<" "<<i0<<" "<<fRecPar->GetMinDensity()
		      <<" "<<track->GetNFoundable()
		      <<" "<<track->GetNShared()<<" chi2 "<<track->GetChi2()<<" "<<endl;
		//		      <<" l "<<clseed[2]->GetLabel(0)<<" "<<clseed[2]->GetLabel(1)<<","
		//      <<clseed[1]->GetLabel(0)<<" "<<clseed[1]->GetLabel(1)<<","
		//      <<clseed[0]->GetLabel(0)<<" "<<clseed[0]->GetLabel(1)<<" "<<endl;
		if (track->GetNumberOfClusters()<(i1-i0)*fRecPar->GetMinDensity() || 
		    track->GetNShared()>0.4*track->GetNumberOfClusters()
		    ||(track->GetNumberOfClusters()-nclustbefore)<=(track->GetNFoundable()-nfoundablebefore)*fRecPar->GetMinDensity()
		    ||(track->GetChi2()/((kUseZCoordinate?2:1)*track->GetNumberOfClusters()-5)>4)) {
		  seed->Reset();
		  seed->~IlcDCHseed();
		  continue;
		}
	      }
	      

	      track->SetSeedType(0);
	      track->SetRemoval(0);
	      arr->AddLast(track); 

	      if (IlcDCHReconstructor::StreamLevel()>0){
		CookLabel(track,0.1);
		TTreeSRedirector &cstream = *fDebugStreamer;
		cstream<<"Seed"
		       <<"seed.="<<track<<"\n"; 
	      }
	      if(fDebug>15) 
		cout<<"Accept seed id="<<arr->GetEntries()-1<<endl;
	      seed = new IlcDCHseed; 	
	      nout2++;
	    }
	  }
	}
      }
    }
  }
  if (IlcDebugLevelClass()>-1){
    timer.Print();
    Info("MakeSeeds3StereoLine","N passed z sel&Rmin, z cut with dphi cor, middle check, density, nout1, nout2 \nSeeding statistic:\t%d\t%d\t%d\t%d\t%d\t%d",
	 nin0,nin1,nin2,nin,nout1,nout2);
  }
  delete seed;
 
}

//_____________________________________________________________________________
void IlcDCHtracker::MakeSeeds3(TObjArray * arr, Int_t i1, Int_t i2,  Float_t cuts[4]) {
  //-----------------------------------------------------------------
  // This function creates track seeds.
  // SEEDING WITH VERTEX CONSTRAIN 
  //-----------------------------------------------------------------
  // cuts[0]   - fP4 cut
  // cuts[1]   - tan(phi)  cut
  // cuts[2]   - zvertex cut
  // cuts[3]   - fP3 cut
  if(!kUseZCoordinate) return MakeSeeds3Stereo(arr,i1,i2,cuts);

  Int_t nin0  = 0;
  Int_t nin1  = 0;
  Int_t nin2  = 0;
  Int_t nin   = 0;
  Int_t nout1 = 0;
  Int_t nout2 = 0;

  Double_t x[5], c[15];
  //  Int_t di = i1-i2;
  //
  IlcDCHseed * seed = new IlcDCHseed;
  
  Double_t x1 =GetRlayer(i1);
  Double_t xx2=GetRlayer(i2);

  Double_t x3g=GetX(), y3g=GetY(), z3g=GetZ();
  Double_t x3=x3g,y3=y3g,z3=z3g;

  Int_t imiddle = (i2+i1)/2;    //middle pad layer index
  Double_t xm[3] = {GetRlayer(imiddle),
		    GetRlayer(imiddle-1),
		    GetRlayer(imiddle+1)}; // radius of middle pad-layer

  const IlcDCHLayer* krm[3]={&GetLayer(imiddle),
			     &GetLayer(imiddle-1),
			     &GetLayer(imiddle+1)}; //middle pad -layer
  //

  const IlcDCHLayer& kr1=GetLayer(i1);

  //
  // change cut on curvature if it can't reach this layer
  // maximal curvature set to reach it
  Double_t dvertexmax  = TMath::Sqrt((x1-x3)*(x1-x3));
  if (dvertexmax*0.5*cuts[0]>1.1){
    if(fDebug>15) cout<<"Rmin before "<<1/cuts[0]<<" after "<<(dvertexmax*0.5)/1.1<<endl;
    cuts[0] = 1.1/(dvertexmax*0.5);
  }
  Double_t r2min = 1/(cuts[0]*cuts[0]);  //minimal square of radius given by cut
  if(fDebug>15) {
    cout<<"Seeding3: i1layer="<<i1<<" i2layer="<<i2<<" im "<<imiddle<<
      " Minimum rmin:"<<TMath::Sqrt(r2min)<<endl;
  }

  // loop over clusters  
  for (Int_t is=0; is < 2*kr1; is++) {
    //
    const IlcDCHcluster *kcl1 = kr1[is/2];
    if (kcl1->IsUsed(10)) continue;
    Double_t y1=(is%2?1:-1)*kcl1->GetImP(), z1=kcl1->GetZ();
    Double_t phi1=kr1.GetWirePhi(kcl1->GetW(),z1); 

    double cosphi,sinphi; sincos(phi1,&sinphi,&cosphi);
    x3= x3g*cosphi+y3g*sinphi;
    y3=-x3g*sinphi+y3g*cosphi;

    x1=GetRlayer(i1,z1);
    // find possible directions    
    Float_t anglez = (z1-z3)/(x1-x3); 
    Float_t extraz = z1 - anglez*(x1-xx2);  // extrapolated z 
    xx2=GetRlayer(i2,extraz);
    extraz=z1 - anglez*(x1-xx2);
    //
    //
    //find   rotation angles relative to line given by vertex and point 1
    Double_t dvertex2 = (x1-x3)*(x1-x3)+(y1-y3)*(y1-y3);
    Double_t dvertex  = TMath::Sqrt(dvertex2);
    Double_t angle13  = TMath::ATan((y1-y3)/(x1-x3));
    Double_t cs13     = cos(-angle13), sn13 = sin(-angle13);            
    
    
    Double_t dddz1=0;  // direction of delta inclination in z axis
    Double_t dddz2=0;
    if ( (z1-z3)>0)
      dddz1 =1;    
    else
      dddz2 =1;
    //

    IlcDCHLayer&  kr2  = GetLayer(i2);
    Int_t  index1 = TMath::Max(kr2.Find(extraz-fRecPar->GetRoadZ()-dddz1*TMath::Abs(z1)*0.05)-1,0);
    Int_t  index2 = TMath::Min(kr2.Find(extraz+fRecPar->GetRoadZ()+dddz2*TMath::Abs(z1)*0.05)+1,kr2);
    
    // rotation angles to p1-p3
    Double_t cs13r     = cos(-angle13)/dvertex, sn13r = sin(-angle13)/dvertex;            
    Double_t x2,   y2,   z2; 
    //
    
    //
    for (Int_t js=2*index1; js < 2*index2; js++) {
      const IlcDCHcluster *kcl2 = kr2[js/2];
      if (kcl2->IsUsed(10)) continue;	

      double phiw=kr2.GetWirePhi(kcl2->GetW(),extraz)+(js%2?1:-1)*kcl2->GetImP()/kr2.GetR()-phi1;
      x2 = TMath::Cos(phiw)*xx2; 
      y2 = TMath::Sin(phiw)*xx2;
      z2 = kcl2->GetZ();	  

      Double_t dxx0 =  (x2-x3)*cs13r;
      Double_t dyy0 =  (x2-x3)*sn13r;
      //
      //calcutate parameters
      //	
      Double_t yy0 =  dyy0 +(y2-y3)*cs13r;
      // stright track
      if (TMath::Abs(yy0)<0.000001) continue;
      Double_t xx0 =  dxx0 -(y2-y3)*sn13r;
      Double_t y0  =  0.5*(xx0*xx0+yy0*yy0-xx0)/yy0;
      Double_t r02 = (0.25+y0*y0)*dvertex2;//dddddddddddddddddddd	
      //curvature (radius) cut
      if(fDebug>15){ 
	cout<<" trackid "<<kcl1->GetLabel(0)<<" "<<kcl1->GetLabel(1)
	    <<" "<<kcl2->GetLabel(0)<<" "<<kcl2->GetLabel(1)<<endl;
	cout<<"r of seed="<<TMath::Sqrt(r02)
	    <<" sign1="<<(is%2?1:-1)<<" sign2="<<(js%2?1:-1)
	    <<" "<<is<<" "<<js<<" "<<index1<<" "<<index2<<endl;
      }
      if (r02<r2min) continue;		
      
      nin0++;	
      //
      Double_t c0  = 1/TMath::Sqrt(r02);
      if (yy0>0) c0*=-1.;	
      
      Double_t dfi0   = 2.*IlcDCHFastMath::FastAsin(dvertex*c0*0.5);
      Double_t dfi1   = 2.*IlcDCHFastMath::FastAsin(TMath::Sqrt(yy0*yy0+(1-xx0)*(1-xx0))*dvertex*c0*0.5);  
      if(fabs(dfi1)<kAlmost0||fabs(dfi0)<kAlmost0) continue;
      //
      //
      Double_t zzzz2    = z1-(z1-z3)*dfi1/dfi0;
      if(fDebug>15){ 
	cout<<"z estimate "<<zzzz2<<" clZ2="<<z2<<" clZ1="<<z1
	    <<" dfi2="<<dfi0<<" dfi1="<<dfi1<<" r2="<<xx2<<" r1="<<x1
	    <<" roadZ="<<fRecPar->GetRoadZ()+fabs((z1-z2)/(x1-xx2)*(xm[2]-xm[1])*0.5)<<endl;
      }
      if (TMath::Abs(zzzz2-z2)>fRecPar->GetRoadZ()+3*sqrt(fRecPar->GetSeedBeamSigmaZ2())) continue;       
      nin1++;              
      //	
      Double_t dip    = (z1-z2)*c0/dfi1;        
      Double_t x0 = (0.5*cs13+y0*sn13)*dvertex*c0;
      //
      
      x[0] = y1;
      x[1] = z1;
      x[2] = x0;
      x[3] = dip;
      x[4] = c0;
      //
      //
      // do we have cluster at the middle ?
      Double_t ym,zm;
      UInt_t dummy; 
      IlcDCHcluster * cm=0;
      bool foundmiddle=kFALSE;

      for(int i=0;i<3;i++){
	if(!GetProlongation(x1,xm[i],x,ym,zm)) continue;
	if(fabs(zm)>krm[i]->GetMaxZ()+fRecPar->GetSafetyZ()) continue;
	if(fDebug>15) 
	  cout<<"on middle "<<TMath::ATan2(ym,xm[i])+phi1<<" "<<zm<<endl;
	cm = krm[i]->FindNearest2(TMath::ATan2(ym,xm[i])+phi1,0.,zm,
				  fRecPar->GetRoadY(),
				  fRecPar->GetRoadZ()+fabs((z1-z2)/(x1-xx2)*(xm[2]-xm[1])*0.5),
				  dummy);
	if ((!cm) || (cm->IsUsed(10))) {	  
	  if(fDebug>15) 
	    cout<<"not found on middles "<<imiddle+(i/2-(i%2))<<endl;
	}else{
	  foundmiddle=kTRUE;
	  break;
	}
      }
      if(!foundmiddle)
	continue;
      
      nin2++;


      UInt_t index=kr1.GetIndex(is/2);
      x[0]=y3;
      x[1]=z3;
      for(int i=0;i<15;i++) c[i]=0.;
      c[0]=fRecPar->GetSeedBeamSigma2();
      c[2]=fRecPar->GetSeedBeamSigmaZ2();
      c[5]=c[9]=1.;c[14]=1.e-4;

      IlcDCHseed *track=new(seed) IlcDCHseed(x3+1e-4, phi1, x, c, index);
      if(fDebug>15){
	track->Print();
      }
      track->PropagateTo(xx2);
      track->SetErrorY2(kcl2->GetSigmaImP2());
      track->SetErrorZ2(kcl2->GetSigmaZ2());
      double dphi=track->GetPhi()-kr2.GetWirePhi(kcl2->GetW(),extraz);
      if(dphi<-TMath::Pi()) dphi+=TMath::Pi()*2;
      if(dphi>TMath::Pi())  dphi-=TMath::Pi()*2;

      track->Update(kcl2,0.,0,
	   	    kr2.GetR(),kr2.GetWirePhi(kcl2->GetW()),kr2.GetTStereoAngle(),
	   	    kUseZCoordinate,((dphi*(js%2?1:-1))<0?-1:1));
      if(fDebug>15)
	track->Print();
      track->PropagateTo(x1);
      track->SetErrorY2(kcl1->GetSigmaImP2());
      track->SetErrorZ2(kcl1->GetSigmaZ2());
      
      dphi=track->GetPhi()-kr1.GetWirePhi(kcl1->GetW(),z1);
      if(dphi<-TMath::Pi()) dphi+=TMath::Pi()*2;
      if(dphi>TMath::Pi())  dphi-=TMath::Pi()*2;

      track->Update(kcl1,0.,0,
	   	    kr1.GetR(),kr1.GetWirePhi(kcl1->GetW()),kr1.GetTStereoAngle(),
	   	    kUseZCoordinate,((dphi*(is%2?1:-1))<0?-1:1));
      track->ResetCovariance(100.);

      dphi=track->GetPhi()-kr1.GetWirePhi(kcl1->GetW(),z1);
      if(dphi<-TMath::Pi()) dphi+=TMath::Pi()*2;
      if(dphi>TMath::Pi())  dphi-=TMath::Pi()*2;

      track->Update(kcl1,0.,0,
 	   	    kr1.GetR(),kr1.GetWirePhi(kcl1->GetW()),kr1.GetTStereoAngle(),
 	   	    kUseZCoordinate,((dphi*(is%2?1:-1))<0?-1:1));
      
      track->PropagateTo(x1);
      if(fDebug>15){
	track->Print();
	cout<<"double rwire[2]={"<<kr2.GetR()<<","<<kr1.GetR()<<"};"<<endl;
	cout<<"double angle[2]={"<<kr2.GetTStereoAngle()<<","
	    <<kr1.GetTStereoAngle()<<"};"<<endl;
	cout<<"double phi[2]={"<<kr2.GetWirePhi(kcl2->GetW())<<","
	    <<kr1.GetWirePhi(kcl1->GetW())<<"};"<<endl;
	cout<<"double mes[2][2]={{"<<kcl2->GetImP()<<","<<kcl2->GetZ()<<"},"<<endl
	    <<"{"<<kcl1->GetImP()<<","<<kcl1->GetZ()<<"}};"<<endl;
	cout<<"double covm[2][3]={{"<<kcl2->GetSigmaImP2()<<",0,"
	    <<kcl2->GetSigmaZ2()<<"},"<<endl
	    <<"{"<<kcl1->GetSigmaImP2()<<",0,"<<kcl1->GetSigmaZ2()<<"}};"<<endl;
      }

      track->SetIsSeeding(kTRUE);
      track->SetSeed1(i1);
      track->SetSeed2(i2);
      track->SetSeedType(3);
      track->SetNumberOfClusters(0);
      track->SetClusterIndex(i1,index,0);
      track->SetClusterPointer(i1,(IlcDCHcluster*)kcl1,0);

      FollowProlongation(*track, imiddle,1);
      if(fDebug>21) track->Dump();
      Int_t foundable,found,shared;
      track->GetClusterStatistic(imiddle,i1, found, foundable, shared, kTRUE);
      if(fDebug>15) 
	cout<<"llllllll "<<found<<" "<<foundable<<" "<<shared
	    <<" "<<(track->GetSigmaY2()+track->GetSigmaZ2())
	    <<" "<<fRecPar->GetCutTrackSigma2()*5<<endl;
      
      if ((found<fRecPar->GetMinDensity()*foundable)  || shared>0.5*found || (track->GetSigmaY2()+track->GetSigmaZ2())>fRecPar->GetCutTrackSigma2()*5){
	//            cout<<"caput 1\n";
	seed->Reset();
	seed->~IlcDCHseed();
	continue;
      }
            
      nin++;
      FollowProlongation(*track, i2,1,1);
      
      
      //Int_t rc = 1;
      track->SetBConstrain(1);
      track->SetLastPoint(i1);  // first cluster in track position
      track->SetFirstPoint(track->GetLastPoint());

      if(fDebug>15) 
	cout<<"mmmmmmmm "<<track->GetNumberOfClusters()
	    <<" "<<i1<<" "<<i2<<" "<<fRecPar->GetMinDensity()
	    <<" "<<track->GetNFoundable()
	    <<" "<<track->GetNShared()<<endl;
      
      if (track->GetNumberOfClusters()<(i1-i2)*fRecPar->GetMinDensity() || 
	  track->GetNumberOfClusters() < track->GetNFoundable()*fRecPar->GetMinDensity() || 
	  track->GetNShared()>0.4*track->GetNumberOfClusters() ) {
	seed->Reset();
	seed->~IlcDCHseed();
	continue;
      }
      nout1++;
      // Z VERTEX CONDITION
      Double_t zv, bz=GetBz();
      if ( !track->GetZAt(0.,bz,zv) ) continue;
      if (TMath::Abs(zv-z3)>cuts[2]) {
	FollowProlongation(*track, TMath::Max(i2-20,0),1,1);
	if ( !track->GetZAt(0.,bz,zv) ) continue;
	if (TMath::Abs(zv-z3)>cuts[2]){
	  FollowProlongation(*track, TMath::Max(i2-40,0),1,1);
	  if ( !track->GetZAt(0.,bz,zv) ) continue;
	  if (TMath::Abs(zv-z3)>cuts[2] &&(track->GetNumberOfClusters() > track->GetNFoundable()*fRecPar->GetMinDensity())){
	    // make seed without constrain
	    IlcDCHseed * track2 = MakeSeed(track,0.2,0.5,1.);
	    FollowProlongation(*track2, i2,1);
	    track2->SetBConstrain(kFALSE);
	    track2->SetSeedType(1);
	    track2->SetRemoval(0);
	    arr->AddLast(track2); 
	    if(fDebug>15) 
	      cout<<"Accept not constrained seed id="<<arr->GetEntries()<<endl;
	    seed->Reset();
	    seed->~IlcDCHseed();
	    continue;		
	  }
	  else{
	    seed->Reset();
	    seed->~IlcDCHseed();
	    continue;
	    
	  }
	}
      }
      track->SetSeedType(0);
      track->SetRemoval(0);
      arr->AddLast(track); 
      if(fDebug>15) 
	cout<<"Accept seed id="<<arr->GetEntries()-1<<endl;
      seed = new IlcDCHseed; 	
      nout2++;

    }
  }
  if (fDebug>3){
    Info("MakeSeeds3","N passed z sel&Rmin, z cut with dphi cor, middle check, density, nout1, nout2 \nSeeding statistic:\t%d\t%d\t%d\t%d\t%d\t%d",
	 nin0,nin1,nin2,nin,nout1,nout2);
  }
  delete seed;
}


void IlcDCHtracker::MakeSeeds5(TObjArray * arr, Int_t i1, Int_t i2,  Float_t cuts[4]) {
  if(!kUseZCoordinate) return ;


  //-----------------------------------------------------------------
  // This function creates track seeds.
  //-----------------------------------------------------------------
  // cuts[0]   - fP4 cut
  // cuts[1]   - tan(phi)  cut
  // cuts[2]   - zvertex cut
  // cuts[3]   - fP3 cut


  Int_t nin0  = 0;
  Int_t nin1  = 0;
  Int_t nin2  = 0;
  Int_t nin   = 0;
  Int_t nout1 = 0;
  Int_t nout2 = 0;
  Int_t nout3 =0;
  Double_t x[5], c[15];

  if(i1-8<0||i2<0) return;
  if(i1>=fSector.GetNLayers()) return;
  //
  // make temporary seed
  IlcDCHseed * seed = new IlcDCHseed;
  //
  //

  // first 3 padlayers
  Double_t r1  = GetRlayer(i1-1);
  const    IlcDCHLayer& kr1 =GetLayer(i1-1);
  //
  Double_t r1p = GetRlayer(i1);
  const    IlcDCHLayer& kr1p=GetLayer(i1);
  //
  Double_t r1m = GetRlayer(i1-2);
  const    IlcDCHLayer& kr1m=GetLayer(i1-2);

  //
  //last 3 padlayer for seeding
  IlcDCHLayer&  kr3  = GetLayer(i1-7);
  Double_t    r3   =  GetRlayer(i1-7);
  //
  IlcDCHLayer&  kr3p  = GetLayer(i1-6);
  Double_t    r3p   = GetRlayer(i1-6);
  //
  IlcDCHLayer&  kr3m  = GetLayer(i1-8);
  Double_t    r3m   = GetRlayer(i1-8);

  //
  //
  // middle padlayer
  Int_t im = i1-4;                           //middle pad layer index
  Double_t rmm[3] = {GetRlayer(im),
		     GetRlayer(im-1),
		     GetRlayer(im+1)}; // radius of middle pad-layer

  double rm=rmm[0];

  const IlcDCHLayer* krm[3]={&GetLayer(im),
                             &GetLayer(im-1),
                             &GetLayer(im+1)}; //middle pad -layer
   //
  //
  Double_t deltax  = r1-r3;
  Double_t dymax   = deltax*cuts[1];
  Double_t dzmax   = deltax*cuts[3];
  //
  // loop over clusters  
  for (Int_t is=0; is < 2*kr1.GetN(); is++) {
    const IlcDCHcluster *kcl1 = kr1[is/2];
    //
    if (kcl1->IsUsed(10)) continue;
    Double_t x1=r1, y1=(is%2?1:-1)*kcl1->GetImP(), z1=kcl1->GetZ();    
    Double_t phi1=kr1.GetWirePhi(kcl1->GetW(),z1); 
    //
    Int_t  index1 = TMath::Max(kr3.Find(z1-dzmax)-1,0);
    Int_t  index2 = TMath::Min(kr3.Find(z1+dzmax)+1,kr3);
    //    
    Double_t x3,y3,z3;
    //
    //
    UInt_t index;
    for (Int_t js=2*index1; js < 2*index2; js++) {
      const IlcDCHcluster *kcl3 = kr3[js/2];
      if (kcl3->IsUsed(10)) continue;
      nin0++;

      double phi3=kr3.GetWirePhi(kcl3->GetW(),kcl3->GetZ())+(js%2?1:-1)*kcl3->GetImP()/kr3.GetR()-phi1;
      x3 = TMath::Cos(phi3)*r3; 
      y3 = TMath::Sin(phi3)*r3;
      z3 = kcl3->GetZ();	

      if(fDebug>15){ 
	cout<<" trackid "<<kcl1->GetLabel(0)<<" "<<kcl1->GetLabel(1)
	    <<" "<<kcl3->GetLabel(0)<<" "<<kcl3->GetLabel(1)<<endl;
	cout<<"x1,y1,z1="<<x1<<" "<<y1<<" "<<z1<<" "<<phi1
	    <<" x3,y3="<<x3<<" "<<y3<<" "<<z3<<" "<<phi3
	    <<" sign1="<<(is%2?1:-1)<<" sign2="<<(js%2?1:-1)
	    <<" dy,zmax "<<dymax<<" "<<dzmax<<endl;
      }

      // apply angular cuts
      if (TMath::Abs(y1-y3)>dymax) continue;
      if (TMath::Abs(z1-z3)>dzmax) continue;
      //
      Double_t angley = (y1-y3)/(x1-x3);
      Double_t anglez = (z1-z3)/(x1-x3);
      //
      Double_t erry = TMath::Abs(angley)*(r1-r1m)*0.5+
	3*sqrt(kcl3->GetSigmaImP2()+kcl1->GetSigmaImP2());
      Double_t errz = TMath::Abs(anglez)*(r1-r1m)*0.5+
	3*sqrt(kcl3->GetSigmaZ2()+kcl1->GetSigmaZ2());
      //
      Double_t yyym , zzzm;

      const IlcDCHcluster *kcm =0;
      int foundmiddle=-1;
      for(int i=0;i<3;i++){
	rm=rmm[i];
	yyym = angley*(rm-x1)+y1;
	zzzm = anglez*(rm-x1)+z1;

        if(fDebug>15) 
          cout<<"on middle "<<TMath::ATan2(yyym,rm)+phi1<<" "<<zzzm<<endl;
        kcm = krm[i]->FindNearest2(TMath::ATan2(yyym,rm)+phi1,0.,zzzm,
				   erry,errz,index);
        if ((!kcm) || (kcm->IsUsed(10))) {          
          if(fDebug>15) 
            cout<<"not found on middles "<<im+(i/2-(i%2))<<endl;
        }else{
          foundmiddle=i;
          break;
        }
      }
      if(foundmiddle<0)
        continue;
 
      erry = TMath::Abs(angley)*(r1-r1m)*0.4+
	3*sqrt(kcl3->GetSigmaImP2()+kcl1->GetSigmaImP2()+kcm->GetSigmaImP2());
      errz = TMath::Abs(anglez)*(r1-r1m)*0.4+
	3*sqrt(kcl3->GetSigmaZ2()+kcl1->GetSigmaZ2()+kcm->GetSigmaZ2());
      //
      //
      //
      Int_t used  =0;
      Int_t found =0;
      //
      // look around first
      yyym = angley*(r1m-x1)+y1;
      zzzm = anglez*(r1m-x1)+z1;
      const IlcDCHcluster *kc1m = kr1m.FindNearest2(TMath::ATan2(yyym,r1m)+phi1,0.,zzzm,erry,errz,index);

      //
      if (kc1m){
	found++;
	if (kc1m->IsUsed(10)) used++;
      }
      yyym = angley*(r1p-x1)+y1;
      zzzm = anglez*(r1p-x1)+z1;
      const IlcDCHcluster *kc1p = kr1p.FindNearest2(TMath::ATan2(yyym,r1p)+phi1,0.,zzzm,erry,errz,index);

      //
      if (kc1p){
	found++;
	if (kc1p->IsUsed(10)) used++;
      }
      if (used>1||found<1){
	if(fDebug>15) 
	  cout<<"not found near middles "<<i1-1<<" used "<<used<<" found "<<found
	      <<" erry "<<erry<<" errz "<<errz<<endl;
	continue;
      }

      //
      // look around last
      Double_t xxxm=(r3m-r3)/(r1-r3)*(x1-x3)+x3;
      yyym = angley*(xxxm-x3)+y3;
      zzzm = anglez*(xxxm-x3)+z3;
      const IlcDCHcluster *kc3m = kr3m.FindNearest2(TMath::ATan2(yyym,xxxm)+phi1,0.,zzzm,erry,errz,index);

      //
      if (kc3m){
	found++;
	if (kc3m->IsUsed(10)) used++;
      }

      xxxm=(r3p-r3)/(r1-r3)*(x1-x3)+x3;
      yyym = angley*(xxxm-x3)+y3;
      zzzm = anglez*(xxxm-x3)+z3;
      const IlcDCHcluster *kc3p = kr3p.FindNearest2(TMath::ATan2(yyym,xxxm)+phi1,0.,zzzm,erry,errz,index);
      //
      if (kc3p){
	found++;
	if (kc3p->IsUsed(10)) used++;
      }

      if (used>1||found<3){
	if(fDebug>15) 
	  cout<<"not found near middles "<<i1-7<<" used "<<used<<" found "<<found
	      <<" erry "<<erry<<" errz "<<errz<<endl;
	continue;
      }
      nin1++;

      //
      Double_t x2,y2,z2;
      for(int jms=0;jms<2;jms++){
	nin2++;

	double phi2=krm[foundmiddle]->GetWirePhi(kcm->GetW(),kcm->GetZ())+
	  (jms%2?1:-1)*kcm->GetImP()/krm[foundmiddle]->GetR()-phi1;
	x2 = TMath::Cos(phi2)*rm; 
	y2 = TMath::Sin(phi2)*rm;
	z2 = kcm->GetZ();
	//
	if(fDebug>15){ 
	  cout<<" trackid "<<kcl1->GetLabel(0)<<" "<<kcl1->GetLabel(1)
	      <<" "<<kcl3->GetLabel(0)<<" "<<kcl3->GetLabel(1)
	      <<" "<<kcm->GetLabel(0)<<" "<<kcm->GetLabel(1)<<endl;
	  cout<<"x1,y1,z1="<<x1<<" "<<y1<<" "<<z1<<" "<<kcl1->GetImP()<<" "<<r1<<" "<<phi1<<endl
	      <<"x2,y2,z2="<<x2<<" "<<y2<<" "<<z2<<" "<<kcm->GetImP()<<" "<<rm<<" "<<phi2<<endl
	      <<"x3,y3,z3="<<x3<<" "<<y3<<" "<<z3<<" "<<kcl3->GetImP()<<" "<<r3<<" "<<phi3<<endl
	      <<" sign1="<<(is%2?1:-1)<<" sign2="<<(js%2?1:-1)<<" sign3="<<(jms%2?1:-1)<<endl;
	}
                  	
	x[0]=y3;
	x[1]=z3;
	x[4]=F1(x1,y1,x2,y2,x3,y3);
	//if (TMath::Abs(x[4]) >= cuts[0]) continue;
	//
	x[2]=F2(x1,y1,x2,y2,x3,y3);
	//
	x[3]=F3n(x1,y1,x2,y2,z1,z2,x[4]);
	//if (TMath::Abs(x[3]) > cuts[3]) continue;
	//
	//
	for(int i=0;i<15;i++) c[i]=0.;
	c[0]=c[2]=100.;
	c[5]=c[9]=1.;c[14]=1.e-4;

	UInt_t index=kr1.GetIndex(is/2);
	IlcDCHseed *track=new(seed) IlcDCHseed(x3, phi1, x, c, index);

	if(fDebug>15)
	  track->Print();
	track->PropagateTo(r3);
	track->SetErrorY2(kcl3->GetSigmaImP2());
	track->SetErrorZ2(kcl3->GetSigmaZ2());
	track->Update(kcl3,0.,0,
		      kr3.GetR(),kr3.GetWirePhi(kcl3->GetW()),kr3.GetTStereoAngle(),
		      kUseZCoordinate);
	if(fDebug>15)
	  track->Print();
	track->PropagateTo(rm);
	track->SetErrorY2(kcm->GetSigmaImP2());
	track->SetErrorZ2(kcm->GetSigmaZ2());
	track->Update(kcm,0.,0,
		      krm[foundmiddle]->GetR(),krm[foundmiddle]->GetWirePhi(kcm->GetW()),
		      krm[foundmiddle]->GetTStereoAngle(),
		      kUseZCoordinate);
	if(fDebug>15)
	  track->Print();
	track->PropagateTo(r1);
	track->SetErrorY2(kcl1->GetSigmaImP2());
	track->SetErrorZ2(kcl1->GetSigmaZ2());
	track->Update(kcl1,0.,0,
		      kr1.GetR(),kr1.GetWirePhi(kcl1->GetW()),kr1.GetTStereoAngle(),
		      kUseZCoordinate);
	track->PropagateTo(r1);
	if(fDebug>15)
	  track->Print();

	track->SetIsSeeding(kTRUE);
	track->SetNumberOfClusters(1);
	track->SetClusterIndex(i1-1,index,0);
	track->SetClusterPointer(i1-1,(IlcDCHcluster*)kcl1,0);

      
	
	nin++;      
	FollowProlongation(*track, i1-7,1);
	nout1++;
	if (fDebug>15) 
	  cout<<"Try to find1 "<<track->GetNumberOfClusters()<<" "<<(i1-i2)*fRecPar->GetMinDensity()<<" "
	      <<track->GetNumberOfClusters()<<" "<<track->GetNFoundable()*fRecPar->GetMinDensity()<<" "
	      <<track->GetSigmaY2()+ track->GetSigmaZ2()<<" "<<fRecPar->GetCutTrackSigma2()<<endl;
	if (track->GetNumberOfClusters() < track->GetNFoundable()*fRecPar->GetMinDensity() || 
	    track->GetNShared()>0.6*track->GetNumberOfClusters() || ( track->GetSigmaY2()+ track->GetSigmaZ2())>fRecPar->GetCutTrackSigma2()){
	  seed->Reset();
	  seed->~IlcDCHseed();
	  continue;
	}
	//Int_t rc = 1;
	FollowProlongation(*track, i2,1);
	track->SetBConstrain(0);
	track->SetLastPoint(i1);  // first cluster in track position
	track->SetFirstPoint(track->GetLastPoint());
	nout2++;	
	if (fDebug>15) 
	  cout<<"Try to find2 "<<track->GetNumberOfClusters()<<" "<<(i1-i2)*fRecPar->GetMinDensity()<<" "
	      <<track->GetNumberOfClusters()<<" "<<track->GetNFoundable()*fRecPar->GetMinDensity()<<" "
	      <<track->GetSigmaY2()+ track->GetSigmaZ2()<<" "<<fRecPar->GetCutTrackSigma2()<<endl;
      
	if (track->GetNumberOfClusters()<(i1-i2)*fRecPar->GetMinDensity() || 
	    track->GetNumberOfClusters()<track->GetNFoundable()*fRecPar->GetMinDensity() || 
	    track->GetNShared()>2. || track->GetChi2()/((kUseZCoordinate?2:1)*track->GetNumberOfClusters()-5)>3 || 
	    ( track->GetSigmaY2()+ track->GetSigmaZ2())>fRecPar->GetCutTrackSigma2() ) {
	  seed->Reset();
	  seed->~IlcDCHseed();
	  continue;
	}
	if (fDebug>15) cout<<"MakeSeed5 Number of seeds "<<arr->GetEntriesFast()<<endl;
// 	{
// 	  FollowProlongation(*track, TMath::Max(i2-10,0),1);
// 	  IlcDCHseed * track2 = MakeSeed(track,0.2,0.5,0.9);
// 	  track2->Print();
	  
// 	  FollowProlongation(*track2, i2,1);
// 	  track2->SetBConstrain(kFALSE);
// 	  track2->SetSeedType(4);
// 	  arr->AddLast(track2); 
// 	  seed->Reset();
// 	  seed->~IlcDCHseed();
// 	}
	track->SetBConstrain(kFALSE);
	track->SetSeedType(4);
 	arr->AddLast(track); 
 	seed = new IlcDCHseed;    
	nout3++;

      }
    }
  }
  
  if (fDebug>3){
    Info("MakeSeeds5","\nSeeding statiistic:\t%d\t%d\t%d\t%d\t%d\t%d\t%d",nin0,nin1,nin2,nin,nout1,nout2,nout3);
  }
  delete seed;
}


//_____________________________________________________________________________
void IlcDCHtracker::MakeSeeds2(TObjArray * arr, Int_t i1, Int_t i2, Float_t */*cuts[4]*/) {
  //-----------------------------------------------------------------
  // This function creates track seeds - without vertex constraint
  //-----------------------------------------------------------------
  // cuts[0]   - fP4 cut        - not applied
  // cuts[1]   - tan(phi)  cut
  // cuts[2]   - zvertex cut    - not applied 
  // cuts[3]   - fP3 cut
  Int_t nin0=0;
  Int_t nin1=0;
  Int_t nin2=0;
  Int_t nin3=0;

  Int_t layer0 = (i1+i2)/2;
  Int_t dlayer = (i1-i2)/2;
  const IlcDCHLayer& kr0=fSector[layer0];
  IlcDCHLayer * kr=0;

  IlcDCHpolyTrack polytrack;
  Int_t nclusters=fSector[layer0];
  IlcDCHseed * seed = new IlcDCHseed;

  Int_t sumused=0;
  Int_t cused=0;
  Int_t cnused=0;
  for (Int_t is=0; is < nclusters; is++) {  //LOOP over clusters
    Int_t nfound =0;
    Int_t nfoundable =0;
    for (Int_t iter =1; iter<2; iter++){   //iterations
      const IlcDCHLayer& krm=fSector[layer0-iter];
      const IlcDCHLayer& krp=fSector[layer0+iter];      
      const IlcDCHcluster * cl= kr0[is];
      
      if (cl->IsUsed(10)) {
	cused++;
      }
      else{
	cnused++;
      }
      Double_t x = kr0.GetR();
      // Initiilczation of the polytrack
      nfound =0;
      nfoundable =0;
      polytrack.Reset();
      //
      Double_t y0= cl->GetY();
      Double_t z0= cl->GetZ();
      Float_t erry = 0;
      Float_t errz = 0;
            
      erry = (0.5)*cl->GetSigmaY2()/TMath::Sqrt(cl->GetQ())*6;	    
      errz = (0.5)*cl->GetSigmaZ2()/TMath::Sqrt(cl->GetQ())*6;      
      polytrack.AddPoint(x,y0,z0,erry, errz);

      sumused=0;
      if (cl->IsUsed(10)) sumused++;


      Float_t roady = (5*TMath::Sqrt(cl->GetSigmaY2())+2*fRecPar->GetRoadY())*iter;
      Float_t roadz = (5*TMath::Sqrt(cl->GetSigmaZ2())+2*fRecPar->GetRoadZ())*iter;
      //
      x = krm.GetR();
      IlcDCHcluster * cl1 = krm.FindNearest(y0,z0,roady,roadz);
      if (cl1) {
	erry = (0.5)*cl1->GetSigmaY2()/TMath::Sqrt(cl1->GetQ())*3;	    
	errz = (0.5)*cl1->GetSigmaZ2()/TMath::Sqrt(cl1->GetQ())*3;
	if (cl1->IsUsed(10))  sumused++;
	polytrack.AddPoint(x,cl1->GetY(),cl1->GetZ(),erry,errz);
      }
      //
      x = krp.GetR();
      IlcDCHcluster * cl2 = krp.FindNearest(y0,z0,roady,roadz);
      if (cl2) {
	erry = (0.5)*cl2->GetSigmaY2()/TMath::Sqrt(cl2->GetQ())*3;	    
	errz = (0.5)*cl2->GetSigmaZ2()/TMath::Sqrt(cl2->GetQ())*3;
	if (cl2->IsUsed(10)) sumused++;	 
	polytrack.AddPoint(x,cl2->GetY(),cl2->GetZ(),erry,errz);
      }
      //
      if (sumused>0) continue;
      nin0++;
      polytrack.UpdateParameters();
      // follow polytrack
      roadz = fRecPar->GetRoadZ()*1.2;
      roady = fRecPar->GetRoadY()*1.2;
      //
      Double_t yn,zn;
      nfoundable = polytrack.GetN();
      nfound     = nfoundable; 
      //
      for (Int_t ddlayer = iter+1; ddlayer<dlayer;ddlayer++){
	Float_t maxdist = 0.8*(1.+3./(ddlayer));
	for (Int_t delta = -1;delta<=1;delta+=2){
	  Int_t layer = layer0+ddlayer*delta;
	  kr = &(fSector[layer]);
	  Double_t xn = kr->GetR();
	  polytrack.GetFitPoint(xn,yn,zn);
	  nfoundable++;
	  IlcDCHcluster * cln = kr->FindNearest(yn,zn,roady,roadz);
	  if (cln) {
	    Float_t dist =  TMath::Sqrt(  (yn-cln->GetY())*(yn-cln->GetY())+(zn-cln->GetZ())*(zn-cln->GetZ()));
	    if (dist<maxdist){
	      erry = (dist+0.3)*cln->GetSigmaY2()/TMath::Sqrt(cln->GetQ())*(1.+1./(ddlayer));	    
	      errz = (dist+0.3)*cln->GetSigmaZ2()/TMath::Sqrt(cln->GetQ())*(1.+1./(ddlayer));
	      /*if (cln->IsUsed(10)) {
		//	printf("used\n");
		sumused++;
		erry*=2;
		errz*=2;
	      }
	      erry=0.1;
	      errz=0.1;*/
	      polytrack.AddPoint(xn,cln->GetY(),cln->GetZ(),erry, errz);
	      nfound++;
	    }
	  }
	}
	if ( (sumused>3) || (sumused>0.5*nfound) || (nfound<fRecPar->GetMinDensity()*nfoundable))  break;     
	polytrack.UpdateParameters();
      }           
    }
    if ( (sumused>3) || (sumused>0.5*nfound))  {
      //printf("sumused   %d\n",sumused);
      continue;
    }
    nin1++;
    Double_t dy,dz;
    polytrack.GetFitDerivation(kr0.GetR(),dy,dz);
    IlcDCHpolyTrack track2;
    
    polytrack.Refit(track2,0.5+TMath::Abs(dy)*0.3,0.4+TMath::Abs(dz)*0.3);
    if (track2.GetN()<fRecPar->GetMinDensity()*nfoundable) continue;
    nin2++;

    if ((nfound>fRecPar->GetMinDensity()*nfoundable) &&( nfoundable>fRecPar->GetMinDensity()*(i1-i2))) {
      //
      // test seed with and without constrain
      for (Int_t constrain=0; constrain<=0;constrain++){
	// add polytrack candidate

	Double_t x[5], c[15];
	Double_t x1,x2,x3,y1,y2,y3,z1,z2,z3;
	track2.GetBoundaries(x3,x1);	
	x2 = (x1+x3)/2.;
	track2.GetFitPoint(x1,y1,z1);
	track2.GetFitPoint(x2,y2,z2);
	track2.GetFitPoint(x3,y3,z3);
	//
	//is track pointing to the vertex ?
	Double_t x0,y0,z0;
	x0=0;
	polytrack.GetFitPoint(x0,y0,z0);

	if (constrain) {
	  x2 = x3;
	  y2 = y3;
	  z2 = z3;
	  
	  x3 = 0;
	  y3 = 0;
	  z3 = 0;
	}
	x[0]=y1;
	x[1]=z1;
	x[4]=F1(x1,y1,x2,y2,x3,y3);
		
	//	if (TMath::Abs(x[4]) >= cuts[0]) continue;  //
	x[2]=F2(x1,y1,x2,y2,x3,y3);
	
	//if (TMath::Abs(x[4]*x1-x[2]) >= cuts[1]) continue;
	//x[3]=F3(x1,y1,x2,y2,z1,z2);
	x[3]=F3n(x1,y1,x3,y3,z1,z3,x[4]);
	//if (TMath::Abs(x[3]) > cuts[3]) continue;

	
	Double_t sy =0.1, sz =0.1;
	Double_t sy1=fRecPar->GetSeedTrackSigma2(), sz1=fRecPar->GetSeedTrackSigma2();
	Double_t sy2=fRecPar->GetSeedTrackSigma2(), sz2=fRecPar->GetSeedTrackSigma2();
	Double_t sy3=fRecPar->GetSeedTrackSigma2();

	if (constrain){
	  sy3=25000*x[4]*x[4]+fRecPar->GetSeedBeamSigma2();
	}
	
	Double_t f40=(F1(x1,y1+sy,x2,y2,x3,y3)-x[4])/sy;
	Double_t f42=(F1(x1,y1,x2,y2+sy,x3,y3)-x[4])/sy;
	Double_t f43=(F1(x1,y1,x2,y2,x3,y3+sy)-x[4])/sy;
	Double_t f20=(F2(x1,y1+sy,x2,y2,x3,y3)-x[2])/sy;
	Double_t f22=(F2(x1,y1,x2,y2+sy,x3,y3)-x[2])/sy;
	Double_t f23=(F2(x1,y1,x2,y2,x3,y3+sy)-x[2])/sy;

	Double_t f30=(F3(x1,y1+sy,x3,y3,z1,z3)-x[3])/sy;
	Double_t f31=(F3(x1,y1,x3,y3,z1+sz,z3)-x[3])/sz;
	Double_t f32=(F3(x1,y1,x3,y3+sy,z1,z3)-x[3])/sy;
	Double_t f34=(F3(x1,y1,x3,y3,z1,z3+sz)-x[3])/sz;

	
	c[0]=sy1;
	c[1]=0.;       c[2]=sz1;
	c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
	c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
	c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
	c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
	c[13]=f30*sy1*f40+f32*sy2*f42;
	c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
	
	Int_t layer1 = GetLayerNumber(x1,z1);

	UInt_t index=0;
	//kr0.GetIndex(is);
	IlcDCHseed *track=new (seed) IlcDCHseed(x1,/*sec*alpha+shift*/0,x,c,index);
	track->SetIsSeeding(kTRUE);
	Int_t rc=FollowProlongation(*track, i2);	
	if (constrain) track->SetBConstrain(1);
	else
	  track->SetBConstrain(0);
	track->SetLastPoint(layer1);  // first cluster in track position
	track->SetFirstPoint(track->GetLastPoint());

	if (rc==0 || track->GetNumberOfClusters()<(i1-i2)*fRecPar->GetMinDensity() || 
	    track->GetNumberOfClusters() < track->GetNFoundable()*fRecPar->GetMinDensity() || 
	    track->GetNShared()>0.4*track->GetNumberOfClusters()) {
	  //delete track;
	  seed->Reset();
	  seed->~IlcDCHseed();
	}
	else {
	  arr->AddLast(track);
	  seed = new IlcDCHseed;
	}
	nin3++;
      }
    }  // if accepted seed
  }
  if (fDebug>3){
    Info("MakeSeeds2","\nSeeding statiistic:\t%d\t%d\t%d\t%d",nin0,nin1,nin2,nin3);
  }
  delete seed;
}


IlcDCHseed *IlcDCHtracker::MakeSeed(IlcDCHseed *track, Float_t r0, Float_t r1, Float_t r2)
{
  //
  //
  //reseed using track points
  Int_t p0 = int(r0*track->GetNumberOfClusters());     // point 0 
  Int_t p1 = int(r1*track->GetNumberOfClusters());
  Int_t p2 = int(r2*track->GetNumberOfClusters());   // last point
  Int_t pp2=0;
  Double_t  x0[3],x1[3],x2[3];
  IlcDCHcluster /**cl0=0,*cl1=0,*cl2=0,*/*cl=0;
//  int lay0=0,lay1=0,lay2=0;
  for (Int_t i=0;i<3;i++){
    x0[i]=-kVeryBig;
    x1[i]=-kVeryBig;
    x2[i]=-kVeryBig;
  }

  // find track position at given ratio of the length
  Int_t index=-1,indexlast=-1;
  int nlayers=fSector.GetNLayers();
  for (Int_t i=0;i<nlayers;i++) for (Int_t j=0;j<kMaxInLayer;j++) {
      if (track->GetClusterPointer(i,j)){
      index++;
      IlcDCHTrackerPoint   *trpoint =track->GetTrackPoint(i,j);
      double rpoint2=trpoint->GetX()*trpoint->GetX()+trpoint->GetY()*trpoint->GetY();
      if ( (index<p0) || x0[0]<-kVeryBig*0.5 ){
	if (rpoint2>1){
	  cl = track->GetClusterPointer(i,j);
	  if (cl!=0){	
	    x0[0] = trpoint->GetX();
	    x0[1] = trpoint->GetY();
	    x0[2] = trpoint->GetZ();
//	    cl0=cl;
//	    lay0=i;
	  }
	}
      }

      if ( (index<p1) && (rpoint2>1)){
	cl = track->GetClusterPointer(i,j);
	if (cl!=0){
	  x1[0] = trpoint->GetX();
	  x1[1] = trpoint->GetY();
	  x1[2] = trpoint->GetZ();
//	  cl1=cl;
//	  lay1=i;
	}
      }
      if ( (index<p2) && (rpoint2>1)){
	cl = track->GetClusterPointer(i,j);
	if (cl!=0){
	  x2[0] = trpoint->GetX();
	  x2[1] = trpoint->GetY();
	  x2[2] = trpoint->GetZ(); 
	  pp2 = i;
//	  cl2=cl;
//	  lay2=i;
	  indexlast=track->GetClusterIndex(i,j);
	}
      }
    }
  }
  
  if(fDebug>15)
    cout<<"xyz0 "<<x0[0]<<" "<<x0[1]<<" "<<x0[2]<<" "
	<<"xyz1 "<<x1[0]<<" "<<x1[1]<<" "<<x1[2]<<" "
	<<"xyz2 "<<x2[0]<<" "<<x2[1]<<" "<<x2[2]<<" "<<endl;

  double rr2=TMath::Hypot(x2[0],x2[1]);
  double phi2=TMath::ATan2(x2[1],x2[0]);
  double cphi,sphi;sincos(phi2,&sphi,&cphi);
  x2[0]=rr2;x2[1]=0;
  double xx=x1[0],yy=x1[1];
  x1[0]= cphi*xx+sphi*yy;
  x1[1]=-sphi*xx+cphi*yy;
  xx=x0[0],yy=x0[1];
  x0[0]= cphi*xx+sphi*yy;
  x0[1]=-sphi*xx+cphi*yy;

  if(fDebug>15)
    cout<<"xyz0 "<<x0[0]<<" "<<x0[1]<<" "<<x0[2]<<" "
	<<"xyz1 "<<x1[0]<<" "<<x1[1]<<" "<<x1[2]<<" "
	<<"xyz2 "<<x2[0]<<" "<<x2[1]<<" "<<x2[2]<<" "<<endl;
  //
  //
  //
  Double_t x[5],c[15];
  //
  x[0]=x2[1];
  x[1]=x2[2];
  x[4]=F1(x2[0],x2[1],x1[0],x1[1],x0[0],x0[1]);
  x[2]=F2(x2[0],x2[1],x1[0],x1[1],x0[0],x0[1]);
  x[3]=F3n(x2[0],x2[1],x0[0],x0[1],x2[2],x0[2],x[4]);
  //  
  Double_t sy =0.1,  sz =0.1; //only forr derivative
  //
  Double_t sy1=fRecPar->GetSeedTrackSigma2()+track->GetSigmaY2(), 
    sz1=fRecPar->GetSeedTrackSigma2()+track->GetSigmaZ2();
  Double_t sy2=fRecPar->GetSeedTrackSigma2()+track->GetSigmaY2(), 
    sz2=fRecPar->GetSeedTrackSigma2()+track->GetSigmaZ2();
  Double_t sy3=fRecPar->GetSeedTrackSigma2()+track->GetSigmaY2();
  //
  Double_t f40=(F1(x2[0],x2[1]+sy,x1[0],x1[1],x0[0],x0[1])-x[4])/sy;
  Double_t f42=(F1(x2[0],x2[1],x1[0],x1[1]+sy,x0[0],x0[1])-x[4])/sy;
  Double_t f43=(F1(x2[0],x2[1],x1[0],x1[1],x0[0],x0[1]+sy)-x[4])/sy;
  Double_t f20=(F2(x2[0],x2[1]+sy,x1[0],x1[1],x0[0],x0[1])-x[2])/sy;
  Double_t f22=(F2(x2[0],x2[1],x1[0],x1[1]+sy,x0[0],x0[1])-x[2])/sy;
  Double_t f23=(F2(x2[0],x2[1],x1[0],x1[1],x0[0],x0[1]+sy)-x[2])/sy;
  //
  Double_t f30=(F3(x2[0],x2[1]+sy,x0[0],x0[1],x2[2],x0[2])-x[3])/sy;
  Double_t f31=(F3(x2[0],x2[1],x0[0],x0[1],x2[2]+sz,x0[2])-x[3])/sz;
  Double_t f32=(F3(x2[0],x2[1],x0[0],x0[1]+sy,x2[2],x0[2])-x[3])/sy;
  Double_t f34=(F3(x2[0],x2[1],x0[0],x0[1],x2[2],x0[2]+sz)-x[3])/sz;
  
  
  c[0]=sy1;
  c[1]=0.;       c[2]=sz1;
  c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
  c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
  c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
  c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
  c[13]=f30*sy1*f40+f32*sy2*f42;
  c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
  
  IlcDCHseed *seed=new  IlcDCHseed(x2[0], phi2, x, c, indexlast);

  seed->SetLastPoint(pp2);
  seed->SetFirstPoint(pp2);
  

  return seed;
}


IlcDCHseed *IlcDCHtracker::ReSeed(IlcDCHseed *track, Float_t r0, Float_t r1, Float_t r2)
{
  //
  //
  //reseed using founded clusters 
  //
  // Find the number of clusters
  int nlayers=fSector.GetNLayers();
  Int_t nclusters = 0;
  for (Int_t ilayer=0;ilayer<nlayers;ilayer++){
    for(int j=0;j<kMaxInLayer;j++)
      if (track->GetClusterIndex(ilayer,j)>=0) nclusters++;
  }
  //
  Int_t ipos[3];
  ipos[0] = TMath::Max(int(r0*nclusters),0);             // point 0 cluster
  ipos[1] = TMath::Min(int(r1*nclusters),nclusters-1);   // 
  ipos[2] = TMath::Min(int(r2*nclusters),nclusters-1);   // last point
  //
  //
  Double_t  xyz[3][3],eyz[3][3];
  Int_t     layer[3],iclayer[3];
  //
  // find track layer position at given ratio of the length
  Int_t index=-1;
  for (Int_t ilayer=0;ilayer<nlayers;ilayer++) for(int j=0;j<kMaxInLayer;j++){    
    if (track->GetClusterIndex(ilayer,j)<0) continue;
    index++;
    for (Int_t ipoint=0;ipoint<3;ipoint++){
      if (index<=ipos[ipoint]) {
	layer[ipoint] = ilayer;
	iclayer[ipoint] = j;
      }
    }        
  }
  //
  //Get cluster and sector position
  for (Int_t ipoint=0;ipoint<3;ipoint++){
    Int_t clindex = track->GetClusterIndex(layer[ipoint],iclayer[ipoint]);
    IlcDCHcluster * cl = GetCluster(clindex);
    if (cl==0) return 0;
    xyz[ipoint][0]  = GetRlayer(layer[ipoint]);
    xyz[ipoint][1]  = cl->GetY();
    xyz[ipoint][2]  = cl->GetZ();
    eyz[ipoint][1]  = cl->GetSigmaY2();
    eyz[ipoint][2]  = cl->GetSigmaZ2();
  }
  //
  //
  // Calculate seed state vector and covariance matrix

  //
  //
  //
  Double_t x[5],c[15];
  //
  x[0]=xyz[2][1];
  x[1]=xyz[2][2];
  x[4]=F1(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]);
  x[2]=F2(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]);
  x[3]=F3n(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1],xyz[2][2],xyz[0][2],x[4]);
  //  
  Double_t sy =0.1,  sz =0.1; //only forr derivative
  //
  Double_t sy1=eyz[0][1], sz1=eyz[0][2];
  Double_t sy2=eyz[1][1], sz2=eyz[1][2];
  Double_t sy3=eyz[2][1];
  //
  Double_t f40=(F1(xyz[2][0],xyz[2][1]+sy,xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1])-x[4])/sy;
  Double_t f42=(F1(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1]+sy,xyz[0][0],xyz[0][1])-x[4])/sy;
  Double_t f43=(F1(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]+sy)-x[4])/sy;
  Double_t f20=(F2(xyz[2][0],xyz[2][1]+sy,xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1])-x[2])/sy;
  Double_t f22=(F2(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1]+sy,xyz[0][0],xyz[0][1])-x[2])/sy;
  Double_t f23=(F2(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]+sy)-x[2])/sy;
  //
  Double_t f30=(F3(xyz[2][0],xyz[2][1]+sy,xyz[0][0],xyz[0][1],xyz[2][2],xyz[0][2])-x[3])/sy;
  Double_t f31=(F3(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1],xyz[2][2]+sz,xyz[0][2])-x[3])/sz;
  Double_t f32=(F3(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1]+sy,xyz[2][2],xyz[0][2])-x[3])/sy;
  Double_t f34=(F3(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1],xyz[2][2],xyz[0][2]+sz)-x[3])/sz;
  
  
  c[0]=sy1;
  c[1]=0.;       c[2]=sz1;
  c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
  c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
  c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
  c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
  c[13]=f30*sy1*f40+f32*sy2*f42;
  c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
  
  IlcDCHseed *seed=new  IlcDCHseed(xyz[2][0], /*sec[2]*fSector.GetAlpha()+fSector.GetAlphaShift()*/0, x, c, 0);
  seed->SetLastPoint(layer[2]);
  seed->SetFirstPoint(layer[2]);  
  return seed;
}


IlcDCHseed *IlcDCHtracker::ReSeed(IlcDCHseed *track,Int_t r0, Bool_t forward)
{
  //
  //
  //reseed using founded clusters 
  //
  Double_t  xyz[3][3],eyz[3][3];
  Int_t     layer[3]={0,0,0};
  //
  // forward direction
  int nlayers=fSector.GetNLayers();
  if (forward){
    for (Int_t ilayer=r0;ilayer<nlayers;ilayer++){
      if (track->GetClusterIndex(ilayer,0)>=0){
	layer[0] = ilayer;
	break;
      }
    }
    for (Int_t ilayer=nlayers;ilayer>r0;ilayer--){
      if (track->GetClusterIndex(ilayer,0)>=0){
	layer[2] = ilayer;
	break;
      }
    }
    for (Int_t ilayer=layer[2]-int(0.7*fRecPar->GetMinNOfClusetrs());ilayer>layer[0];ilayer--){
      if (track->GetClusterIndex(ilayer,0)>=0){
	layer[1] = ilayer;
	break;
      }
    }
    //
  }
  if (!forward){
    for (Int_t ilayer=0;ilayer<r0;ilayer++){
      if (track->GetClusterIndex(ilayer,0)>=0){
	layer[0] = ilayer;
	break;
      }
    }
    for (Int_t ilayer=r0;ilayer>0;ilayer--){
      if (track->GetClusterIndex(ilayer,0)>=0){
	layer[2] = ilayer;
	break;
      }
    }    
    for (Int_t ilayer=layer[2]-int(0.7*fRecPar->GetMinNOfClusetrs());ilayer>layer[0];ilayer--){
      if (track->GetClusterIndex(ilayer,0)>=0){
	layer[1] = ilayer;
	break;
      }
    } 
  }
  //
  if ((layer[2]-layer[0])<fRecPar->GetMinNOfClusetrs()-1) return 0;
  if (layer[1]==0) return 0;
  //
  //
  //Get cluster and sector position
  for (Int_t ipoint=0;ipoint<3;ipoint++){
    Int_t clindex = track->GetClusterIndex(layer[ipoint],0);
    IlcDCHcluster * cl = GetCluster(clindex);
    if (cl==0) {
      return 0;
    }
    xyz[ipoint][0]  = GetRlayer(layer[ipoint]);
    IlcDCHTrackerPoint * point = track->GetTrackPoint(layer[ipoint],0);    
    if (point&&ipoint<2){
      //
       xyz[ipoint][1]  = point->GetY();
       xyz[ipoint][2]  = point->GetZ();
    }
    else{
      xyz[ipoint][1]  = cl->GetY();
      xyz[ipoint][2]  = cl->GetZ();
    }
    eyz[ipoint][1]  = cl->GetSigmaY2();
    eyz[ipoint][2]  = cl->GetSigmaZ2();
  }
  //
  //
  //
  //
  // Calculate seed state vector and covariance matrix

  Double_t x[5],c[15];
  //
  x[0]=xyz[2][1];
  x[1]=xyz[2][2];
  x[4]=F1(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]);
  x[2]=F2(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]);
  x[3]=F3n(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1],xyz[2][2],xyz[0][2],x[4]);
  //  
  Double_t sy =0.1,  sz =0.1; //only forr derivation
  //
  Double_t sy1=eyz[0][1], sz1=eyz[0][2];
  Double_t sy2=eyz[1][1], sz2=eyz[1][2];
  Double_t sy3=eyz[2][1];
  //
  Double_t f40=(F1(xyz[2][0],xyz[2][1]+sy,xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1])-x[4])/sy;
  Double_t f42=(F1(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1]+sy,xyz[0][0],xyz[0][1])-x[4])/sy;
  Double_t f43=(F1(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]+sy)-x[4])/sy;
  Double_t f20=(F2(xyz[2][0],xyz[2][1]+sy,xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1])-x[2])/sy;
  Double_t f22=(F2(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1]+sy,xyz[0][0],xyz[0][1])-x[2])/sy;
  Double_t f23=(F2(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]+sy)-x[2])/sy;
  //
  Double_t f30=(F3(xyz[2][0],xyz[2][1]+sy,xyz[0][0],xyz[0][1],xyz[2][2],xyz[0][2])-x[3])/sy;
  Double_t f31=(F3(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1],xyz[2][2]+sz,xyz[0][2])-x[3])/sz;
  Double_t f32=(F3(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1]+sy,xyz[2][2],xyz[0][2])-x[3])/sy;
  Double_t f34=(F3(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1],xyz[2][2],xyz[0][2]+sz)-x[3])/sz;
  
  
  c[0]=sy1;
  c[1]=0.;       c[2]=sz1;

  c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
  c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
  c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
  c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
  c[13]=f30*sy1*f40+f32*sy2*f42;
  c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
  
  IlcDCHseed *seed=new  IlcDCHseed(xyz[2][0], /*sec[2]*fSector.GetAlpha()+fSector.GetAlphaShift()*/0, x, c, 0);
  seed->SetLastPoint(layer[2]);
  seed->SetFirstPoint(layer[2]);  
  for (Int_t i=layer[0];i<layer[2];i++)for (Int_t j=0;j<kMaxInLayer;j++){
      seed->SetClusterIndex(i, track->GetClusterIndex(i,j),j);
  }

  return seed;
}
/*
void  IlcDCHtracker::FindKinks(TObjArray * array, IlcESD *esd)
{
  //
  //  find kinks
  //
  //
  int nlayers=fSector.GetNLayers();

  TObjArray *kinks= new TObjArray(10000);
  //  TObjArray *v0s= new TObjArray(10000);
  Int_t nentries = array->GetEntriesFast();
  IlcHelix *helixes      = new IlcHelix[nentries];
  Int_t    *sign         = new Int_t[nentries];
  Int_t    *nclusters    = new Int_t[nentries];
  Float_t  *alpha        = new Float_t[nentries];
  IlcKink  *kink         = new IlcKink();
  Int_t      * usage     = new Int_t[nentries];
  Float_t  *zm           = new Float_t[nentries];
  Float_t  *z0           = new Float_t[nentries]; 
  Float_t  *fim          = new Float_t[nentries];
  Float_t  *shared       = new Float_t[nentries];
  Bool_t   *circular     = new Bool_t[nentries];
  Float_t *dca          = new Float_t[nentries];
  //const IlcESDVertex * primvertex = esd->GetVertex();
  //
  //  nentries = array->GetEntriesFast();
  //
  
  //
  //
  for (Int_t i=0;i<nentries;i++){
    sign[i]=0;
    usage[i]=0;
    IlcDCHseed* track = (IlcDCHseed*)array->At(i);    
    if (!track) continue;
    track->SetCircular(0);
    shared[i] = kFALSE;
    track->UpdatePoints();
    if (( track->GetPoints()[2]- track->GetPoints()[0])>5 && track->GetPoints()[3]>0.8){
    }
    nclusters[i]=track->GetNumberOfClusters();
    alpha[i] = track->GetAlpha();
    new (&helixes[i]) IlcHelix(*track);
    Double_t xyz[3];
    helixes[i].Evaluate(0,xyz);
    sign[i] = (track->GetC()>0) ? -1:1;
    Double_t x,y,z;
    x=GetRlayer(int(0.5*nlayers)); //dddddddddddddddddddddd
    if (track->GetProlongation(x,y,z)){
      zm[i]  = z;
      fim[i] = alpha[i]+TMath::ATan2(y,x);
    }
    else{
      zm[i]  = track->GetZ();
      fim[i] = alpha[i];
    }   
    z0[i]=1000;
    circular[i]= kFALSE;
    if (track->GetProlongation(0,y,z))  z0[i] = z;
    dca[i] = track->GetD(0,0);    
  }
  //
  //
  TStopwatch timer;
  timer.Start();
  Int_t ncandidates =0;
  Int_t nall =0;
  Int_t ntracks=0; 
  Double_t phase[2][2],radius[2];

  //
  // Find circling track
  TTreeSRedirector &cstream = *fDebugStreamer;
  //
  for (Int_t i0=0;i0<nentries;i0++){
    IlcDCHseed * track0 = (IlcDCHseed*)array->At(i0);
    if (!track0) continue;
    if (track0->GetNumberOfClusters()<3*fRecPar->GetMinNOfClusetrs()) continue;
    if (TMath::Abs(1./track0->GetC())>fParam->GetOuterRadius()) continue;
    for (Int_t i1=i0+1;i1<nentries;i1++){
      IlcDCHseed * track1 = (IlcDCHseed*)array->At(i1);
      if (!track1) continue;
      if (track1->GetNumberOfClusters()<3*fRecPar->GetMinNOfClusetrs())                  continue;
      if ( TMath::Abs(track1->GetTgl()+track0->GetTgl())>0.1) continue;
      if (track0->GetBConstrain()&&track1->GetBConstrain()) continue;
      if (TMath::Abs(1./track1->GetC())>fParam->GetOuterRadius()) continue;
      if (track1->Get1Pt()*track0->Get1Pt()>0)      continue;
      if (track1->GetTgl()*track0->GetTgl()>0)      continue;
      if (max(TMath::Abs(1./track0->GetC()),TMath::Abs(1./track1->GetC()))>fParam->GetOuterRadius()) continue;
      if (track0->GetBConstrain()&&TMath::Abs(track1->Get1Pt())<TMath::Abs(track0->Get1Pt())) continue; //returning - lower momenta
      if (track1->GetBConstrain()&&TMath::Abs(track0->Get1Pt())<TMath::Abs(track1->Get1Pt())) continue; //returning - lower momenta
       //
      Float_t mindcar = TMath::Min(TMath::Abs(dca[i0]),TMath::Abs(dca[i1]));
      if (mindcar<fRecPar->GetMinDCAr())   continue;
      Float_t mindcaz = TMath::Min(TMath::Abs(z0[i0]-GetZ()),TMath::Abs(z0[i1]-GetZ()));
      if (mindcaz<fRecPar->GetMinDCAz()) continue;
      if (mindcar+mindcaz<2*(fRecPar->GetMinDCAr()+fRecPar->GetMinDCAz())) continue;
      //
      //
      Float_t xc0 = helixes[i0].GetHelix(6);
      Float_t yc0 = helixes[i0].GetHelix(7);
      Float_t r0  = helixes[i0].GetHelix(8);
      Float_t xc1 = helixes[i1].GetHelix(6);
      Float_t yc1 = helixes[i1].GetHelix(7);
      Float_t r1  = helixes[i1].GetHelix(8);
	
      Float_t rmean = (r0+r1)*0.5;
      Float_t delta =TMath::Sqrt((xc1-xc0)*(xc1-xc0)+(yc1-yc0)*(yc1-yc0));
      //if (delta>30) continue;
      if (delta>rmean*0.25) continue;
      if (TMath::Abs(r0-r1)/rmean>0.3) continue; 
      //
      Int_t npoints = helixes[i0].GetRPHIintersections(helixes[i1], phase, radius,10);
      if (npoints==0) continue;
      helixes[i0].GetClosestPhases(helixes[i1], phase);
      //
      Double_t xyz0[3];
      Double_t xyz1[3];
      Double_t hangles[3];
      helixes[i0].Evaluate(phase[0][0],xyz0);
      helixes[i1].Evaluate(phase[0][1],xyz1);

      helixes[i0].GetAngle(phase[0][0],helixes[i1],phase[0][1],hangles);
      Double_t deltah[2],deltabest;
      if (hangles[2]<2.8) continue;

//      cstream<<"C"<<track0->fLab<<track1->fLab<<
//	track0->fP3<<track1->fP3<<
//	track0->fP4<<track1->fP4<<
//	delta<<rmean<<npoints<<
//	hangles[0]<<hangles[2]<<
//	xyz0[2]<<xyz1[2]<<radius[0]<<"\n";

      if (npoints>0){
	Int_t ibest=0;
	helixes[i0].ParabolicDCA(helixes[i1],phase[0][0],phase[0][1],radius[0],deltah[0],2);
	if (npoints==2){
	  helixes[i0].ParabolicDCA(helixes[i1],phase[1][0],phase[1][1],radius[1],deltah[1],2);
	  if (deltah[1]<deltah[0]) ibest=1;
	}
	deltabest  = TMath::Sqrt(deltah[ibest]);
	helixes[i0].Evaluate(phase[ibest][0],xyz0);
	helixes[i1].Evaluate(phase[ibest][1],xyz1);
	helixes[i0].GetAngle(phase[ibest][0],helixes[i1],phase[ibest][1],hangles);
	Double_t radiusbest = TMath::Sqrt(radius[ibest]);
	//
	if (deltabest>6) continue;
	if (mindcar+mindcaz<4*(fRecPar->GetMinDCAr()+fRecPar->GetMinDCAz()) && (hangles[2]<3.12||deltabest>3)) continue;
	Bool_t sign =kFALSE;
	if (hangles[2]>3.06) sign =kTRUE;
	//
	if (sign){
	  circular[i0] = kTRUE;
	  circular[i1] = kTRUE;
	  if (TMath::Abs(track0->Get1Pt())<TMath::Abs(track1->Get1Pt())){
	    track0->SetCircular(track0->GetCircular()+1);
	    track1->SetCircular(track1->GetCircular()+2);
	  }
	  else{
	    track1->SetCircular(track1->GetCircular()+1);
	    track0->SetCircular(track0->GetCircular()+2);
	  }
	}		
	if (sign&&IlcDCHReconstructor::StreamLevel()>1){	  
	  //debug stream	  
	  Int_t lab0=track0->GetLabel();
	  Int_t lab1=track1->GetLabel();
	  cstream<<"Curling"<<
	    "lab0="<<lab0<<
	    "lab1="<<lab1<<   
	    "Tr0.="<<track0<<
	    "Tr1.="<<track1<<	   
	    "dca0="<<dca[i0]<<
	    "dca1="<<dca[i1]<<
	    "mindcar="<<mindcar<<
	    "mindcaz="<<mindcaz<<
	    "delta="<<delta<<
	    "rmean="<<rmean<<
	    "npoints="<<npoints<<                      
	    "hangles0="<<hangles[0]<<
	    "hangles2="<<hangles[2]<<                    
	    "xyz0="<<xyz0[2]<<
	    "xyzz1="<<xyz1[2]<<
	    "z0="<<z0[i0]<<
	    "z1="<<z0[i1]<<
	    "radius="<<radiusbest<<
	    "deltabest="<<deltabest<< 
	    "phase0="<<phase[ibest][0]<<
	    "phase1="<<phase[ibest][1]<<
	    "\n"; 	  	  
	}
      }
    }
  }
  //
  //  Finf kinks loop
  // 
  //
  for (Int_t i =0;i<nentries;i++){
    if (sign[i]==0) continue;
    IlcDCHseed * track0 = (IlcDCHseed*)array->At(i);
    ntracks++;
    //
    Double_t cradius0 = 0.5*fParam->GetInnerRadius()*0.5*fParam->GetInnerRadius();
    Double_t cradius1 = (fParam->GetOuterRadius()+0.5*fParam->GetInnerRadius())*
      (fParam->GetOuterRadius()+0.5*fParam->GetInnerRadius());
    Double_t cdist1=8.;
    Double_t cdist2=8.;
    Double_t cdist3=0.55; 
    for (Int_t j =i+1;j<nentries;j++){
      nall++;
      if (sign[j]*sign[i]<1) continue;
      if ( (nclusters[i]+nclusters[j])>200) continue;
      if ( (nclusters[i]+nclusters[j])<80) continue;
      if ( TMath::Abs(zm[i]-zm[j])>fRecPar->GetDeltaZ()) continue;
      if ( TMath::Abs(fim[i]-fim[j])>0.6 && TMath::Abs(fim[i]-fim[j])<5.7 ) continue;
      //IlcDCHseed * track1 = (IlcDCHseed*)array->At(j);  Double_t phase[2][2],radius[2];    
      Int_t npoints = helixes[i].GetRPHIintersections(helixes[j], phase, radius,20);
      if (npoints<1) continue;
      // cuts on radius      
      if (npoints==1){
	if (radius[0]<cradius0||radius[0]>cradius1) continue;
      }
      else{
	if ( (radius[0]<cradius0||radius[0]>cradius1) && (radius[1]<cradius0||radius[1]>cradius1) ) continue;
      }
      //      
      Double_t delta1=10000,delta2=10000;
      // cuts on the intersection radius
      helixes[i].LinearDCA(helixes[j],phase[0][0],phase[0][1],radius[0],delta1);
      if (radius[0]<20&&delta1<1) continue; //intersection at vertex
      if (radius[0]<10&&delta1<3) continue; //intersection at vertex
      if (npoints==2){ 
	helixes[i].LinearDCA(helixes[j],phase[1][0],phase[1][1],radius[1],delta2);
	if (radius[1]<20&&delta2<1) continue;  //intersection at vertex
	if (radius[1]<10&&delta2<3) continue;  //intersection at vertex	
      }
      //
      Double_t distance1 = TMath::Min(delta1,delta2);
      if (distance1>cdist1) continue;  // cut on DCA linear approximation
      //
      npoints = helixes[i].GetRPHIintersections(helixes[j], phase, radius,20);
      helixes[i].ParabolicDCA(helixes[j],phase[0][0],phase[0][1],radius[0],delta1);
      if (radius[0]<20&&delta1<1) continue; //intersection at vertex
      if (radius[0]<10&&delta1<3) continue; //intersection at vertex
      //
      if (npoints==2){ 
	helixes[i].ParabolicDCA(helixes[j],phase[1][0],phase[1][1],radius[1],delta2);	
	if (radius[1]<20&&delta2<1) continue;  //intersection at vertex
	if (radius[1]<10&&delta2<3) continue;  //intersection at vertex	
      }            
      distance1 = TMath::Min(delta1,delta2);
      Float_t rkink =0;
      if (delta1<delta2){
	rkink = TMath::Sqrt(radius[0]);
      }
      else{
	rkink = TMath::Sqrt(radius[1]);
      }
      if (distance1>cdist2) continue;
      //
      //
      IlcDCHseed * track1 = (IlcDCHseed*)array->At(j);
      //
      //
      Int_t layer0 = GetLayerNumber(rkink,0.); 
      if (layer0<10)  continue;
      if (layer0>nlayers-10) continue;
      //
      //
      Float_t dens00=-1,dens01=-1;
      Float_t dens10=-1,dens11=-1;
      //
      Int_t found,foundable,shared;
      track0->GetClusterStatistic(0,layer0-5, found, foundable,shared,kFALSE);
      if (foundable>5) dens00 = Float_t(found)/Float_t(foundable);
      track0->GetClusterStatistic(layer0+5,nlayers-5, found, foundable,shared,kFALSE);
      if (foundable>5) dens01 = Float_t(found)/Float_t(foundable);
      //
      track1->GetClusterStatistic(0,layer0-5, found, foundable,shared,kFALSE);
      if (foundable>10) dens10 = Float_t(found)/Float_t(foundable);
      track1->GetClusterStatistic(layer0+5,nlayers-5, found, foundable,shared,kFALSE);
      if (foundable>10) dens11 = Float_t(found)/Float_t(foundable);
      //     
      if (dens00<dens10 && dens01<dens11) continue;
      if (dens00>dens10 && dens01>dens11) continue;
      if (TMath::Max(dens00,dens10)<0.1)  continue;
      if (TMath::Max(dens01,dens11)<0.3)  continue;
      //
      if (TMath::Min(dens00,dens10)>0.6)  continue;
      if (TMath::Min(dens01,dens11)>0.6)  continue;

      //
      IlcDCHseed * ktrack0, *ktrack1;
      if (dens00>dens10){
	ktrack0 = track0;
	ktrack1 = track1;
      }
      else{
	ktrack0 = track1;
	ktrack1 = track0;
      }
      if (TMath::Abs(ktrack0->GetC())>5) continue; // cut on the curvature for mother particle
      IlcExternalTrackParam paramm(*ktrack0);
      IlcExternalTrackParam paramd(*ktrack1);
      if (layer0>60&&ktrack1->GetReference().GetR()>fParam->GetInnerRadius()) new (&paramd) IlcExternalTrackParam(ktrack1->GetReference()); 
      //
      //
      kink->SetMother(paramm);
      kink->SetDaughter(paramd);
      kink->Update();

      Float_t x[3] = { kink->GetPosition()[0],kink->GetPosition()[1],kink->GetPosition()[2]};
      layer0 = GetLayerNumber(x[0],x[2]); 

      if (kink->GetR()<fParam->GetInnerRadius()+fRecPar->GetDistFromEdge()) continue;
      if (kink->GetR()>fParam->GetOuterRadius()-fRecPar->GetDistFromEdge()) continue;
      if (kink->GetPosition()[2]/kink->GetR()>IlcDCHReconstructor::GetCtgRange()) continue;  //out of fiducial volume
      if (kink->GetDistance()>cdist3) continue;
      Float_t dird = kink->GetDaughterP()[0]*kink->GetPosition()[0]+kink->GetDaughterP()[1]*kink->GetPosition()[1];  // rough direction estimate
      if (dird<0) continue;

      Float_t dirm = kink->GetMotherP()[0]*kink->GetPosition()[0]+kink->GetMotherP()[1]*kink->GetPosition()[1];  // rough direction estimate
      if (dirm<0) continue;
      Float_t mpt = TMath::Sqrt(kink->GetMotherP()[0]*kink->GetMotherP()[0]+kink->GetMotherP()[1]*kink->GetMotherP()[1]);
      if (mpt<0.2) continue;

      if (mpt<1){
	//for high momenta momentum not defined well in first iteration
	Double_t qt   =  TMath::Sin(kink->GetAngle(2))*ktrack1->GetP();
	if (qt>0.35) continue; 
      }
      
      kink->SetLabel(CookLabel(ktrack0,0.4,0,layer0),0);
      kink->SetLabel(CookLabel(ktrack1,0.4,layer0,nlayers),1);
      if (dens00>dens10){
	kink->SetTPCDensity(dens00,0,0);
	kink->SetTPCDensity(dens01,0,1);
	kink->SetTPCDensity(dens10,1,0);
	kink->SetTPCDensity(dens11,1,1);
	kink->SetIndex(i,0);
	kink->SetIndex(j,1);
      }
      else{
	kink->SetTPCDensity(dens10,0,0);
	kink->SetTPCDensity(dens11,0,1);
	kink->SetTPCDensity(dens00,1,0);
	kink->SetTPCDensity(dens01,1,1);
	kink->SetIndex(j,0);
	kink->SetIndex(i,1);
      }

      if (mpt<1||kink->GetAngle(2)>0.1){
	//	angle and densities  not defined yet
	if (kink->GetTPCDensityFactor()<0.8) continue;
	if ((2-kink->GetTPCDensityFactor())*kink->GetDistance() >0.25) continue;
	if (kink->GetAngle(2)*ktrack0->GetP()<0.003) continue; //too small angle
	if (kink->GetAngle(2)>0.2&&kink->GetTPCDensityFactor()<1.15) continue;
	if (kink->GetAngle(2)>0.2&&kink->GetTPCDensity(0,1)>0.05) continue;

	Float_t criticalangle = track0->GetSigmaSnp2()+track0->GetSigmaTgl2();
	criticalangle+= track1->GetSigmaSnp2()+track1->GetSigmaTgl2();
	criticalangle= 3*TMath::Sqrt(criticalangle);
	if (criticalangle>0.02) criticalangle=0.02;
	if (kink->GetAngle(2)<criticalangle) continue;
      }
      //
      kink->SetShapeFactor(-1.);
      kinks->AddLast(kink);
      kink = new IlcKink;
      ncandidates++;
    }
  }
  //
  // sort the kinks according quilcty - and refit them towards vertex
  //
  Int_t       nkinks    = kinks->GetEntriesFast();
  Float_t    *quilcty   = new Float_t[nkinks];
  Int_t      *indexes   = new Int_t[nkinks];
  IlcDCHseed *mothers   = new IlcDCHseed[nkinks];
  IlcDCHseed *daughters = new IlcDCHseed[nkinks];
  //
  //
  for (Int_t i=0;i<nkinks;i++){
    quilcty[i] =100000;
    IlcKink *kink = (IlcKink*)kinks->At(i);
    //
    // refit kinks towards vertex
    // 
    Int_t index0 = kink->GetIndex(0);
    Int_t index1 = kink->GetIndex(1);
    IlcDCHseed * ktrack0 = (IlcDCHseed*)array->At(index0);
    IlcDCHseed * ktrack1 = (IlcDCHseed*)array->At(index1);
    //
    Int_t sumn=ktrack0->GetNumberOfClusters()+ktrack1->GetNumberOfClusters();
    //
    // Refit Kink under if too small angle
    //
    if (kink->GetAngle(2)<0.05){
      kink->SetTPCRow0(GetLayerNumber(kink->GetR(),kink->GetPosition()[2]));
      Int_t layer0 = kink->GetTPCRow0();
      Int_t dlayer = Int_t(2.+0.5/(0.05+kink->GetAngle(2)));
      //
      //
      Int_t last  = layer0-dlayer;
      if (last<40) last=40;
      if (last<ktrack0->GetFirstPoint()+25) last = ktrack0->GetFirstPoint()+25;
      IlcDCHseed* seed0 = ReSeed(ktrack0,last,kFALSE);
      //
      //
      Int_t first = layer0+dlayer;
      if (first>130) first=130;
      if (first>ktrack1->GetLastPoint()-25) first = TMath::Max(ktrack1->GetLastPoint()-25,30);
      IlcDCHseed* seed1 = ReSeed(ktrack1,first,kTRUE);
      //
      if (seed0 && seed1){
	kink->SetStatus(1,8);
	if (RefitKink(*seed0,*seed1,*kink)) kink->SetStatus(1,9);
	layer0 = GetLayerNumber(kink->GetR(),kink->GetPosition()[2]);
	sumn = seed0->GetNumberOfClusters()+seed1->GetNumberOfClusters();
	new (&mothers[i])   IlcDCHseed(*seed0);
	new (&daughters[i]) IlcDCHseed(*seed1); 
      }
      else{
	delete kinks->RemoveAt(i);
	if (seed0) delete seed0;
	if (seed1) delete seed1;
	continue;
      }
      if (kink->GetDistance()>0.5 || kink->GetR()<fParam->GetInnerRadius()+fRecPar->GetDistFromEdge() || kink->GetR()>fParam->GetOuterRadius()-fRecPar->GetDistFromEdge()) {
	delete kinks->RemoveAt(i);
	if (seed0) delete seed0;
	if (seed1) delete seed1;
	continue;
      }
      //
      delete seed0;
      delete seed1;            
    }
    //
    if (kink) quilcty[i] = 160*((0.1+kink->GetDistance())*(2.-kink->GetTPCDensityFactor()))/(sumn+40.);  //the longest -clossest will win
  }
  TMath::Sort(nkinks,quilcty,indexes,kFALSE);
  //
  //remove double find kinks
  //
  for (Int_t ikink0=1;ikink0<nkinks;ikink0++){
    IlcKink * kink0 = (IlcKink*) kinks->At(indexes[ikink0]);
    if (!kink0) continue;
    //
    for (Int_t ikink1=0;ikink1<ikink0;ikink1++){
      if (!kink0) continue;
      IlcKink * kink1 = (IlcKink*) kinks->At(indexes[ikink1]);
      if (!kink1) continue;
      // if not close kink continue
      if (TMath::Abs(kink1->GetPosition()[2]-kink0->GetPosition()[2])>10) continue;
      if (TMath::Abs(kink1->GetPosition()[1]-kink0->GetPosition()[1])>10) continue;
      if (TMath::Abs(kink1->GetPosition()[0]-kink0->GetPosition()[0])>10) continue;
      //
      IlcDCHseed &mother0   = mothers[indexes[ikink0]];
      IlcDCHseed &daughter0 = daughters[indexes[ikink0]];
      IlcDCHseed &mother1   = mothers[indexes[ikink1]];
      IlcDCHseed &daughter1 = daughters[indexes[ikink1]];
      Int_t layer0 = (kink0->GetTPCRow0()+kink1->GetTPCRow0())/2;
      //
      Int_t same  = 0;
      Int_t both  = 0;
      Int_t samem = 0;
      Int_t bothm = 0;
      Int_t samed = 0;
      Int_t bothd = 0;
      //
      for (Int_t i=0;i<layer0;i++){
	if (mother0.GetClusterIndex(i)>=0 && mother1.GetClusterIndex(i)>=0){
	  both++;
	  bothm++;
	  if (mother0.GetClusterIndex(i)==mother1.GetClusterIndex(i)){
	    same++;
	    samem++;
	  }
	}
      }

      for (Int_t i=layer0;i<nlayers;i++){
	if (daughter0.GetClusterIndex(i)>=0 && daughter0.GetClusterIndex(i)>=0){
	  both++;
	  bothd++;
	  if (mother0.GetClusterIndex(i)==mother1.GetClusterIndex(i)){
	    same++;
	    samed++;
	  }
	}
      }
      Float_t ratio = Float_t(same+1)/Float_t(both+1);
      Float_t ratiom = Float_t(samem+1)/Float_t(bothm+1);
      Float_t ratiod = Float_t(samed+1)/Float_t(bothd+1);
      if (ratio>0.3 && ratiom>0.5 &&ratiod>0.5) {
	Int_t sum0 = mother0.GetNumberOfClusters()+daughter0.GetNumberOfClusters();
	Int_t sum1 = mother1.GetNumberOfClusters()+daughter1.GetNumberOfClusters();
	if (sum1>sum0){
	  shared[kink0->GetIndex(0)]= kTRUE;
	  shared[kink0->GetIndex(1)]= kTRUE;	  
	  delete kinks->RemoveAt(indexes[ikink0]);
	}
	else{
	  shared[kink1->GetIndex(0)]= kTRUE;
	  shared[kink1->GetIndex(1)]= kTRUE;	  
	  delete kinks->RemoveAt(indexes[ikink1]);
	}
      }
    }
  }


  for (Int_t i=0;i<nkinks;i++){
    IlcKink * kink = (IlcKink*) kinks->At(indexes[i]);
    if (!kink) continue;
    kink->SetTPCRow0(GetLayerNumber(kink->GetR(),kink->GetPosition()[2]));
    Int_t index0 = kink->GetIndex(0);
    Int_t index1 = kink->GetIndex(1);
    if (circular[index0]||(circular[index1]&&kink->GetDistance()>0.2)) continue;
    kink->SetMultiple(usage[index0],0);
    kink->SetMultiple(usage[index1],1);
    if (kink->GetMultiple()[0]+kink->GetMultiple()[1]>2) continue;
    if (kink->GetMultiple()[0]+kink->GetMultiple()[1]>0 && quilcty[indexes[i]]>0.2) continue;
    if (kink->GetMultiple()[0]+kink->GetMultiple()[1]>0 && kink->GetDistance()>0.2) continue;
    if (circular[index0]||(circular[index1]&&kink->GetDistance()>0.1)) continue;

    IlcDCHseed * ktrack0 = (IlcDCHseed*)array->At(index0);
    IlcDCHseed * ktrack1 = (IlcDCHseed*)array->At(index1);
    if (!ktrack0 || !ktrack1) continue;
    Int_t index = esd->AddKink(kink);
    //
    //
    if ( ktrack0->GetKinkIndex(0)==0 && ktrack1->GetKinkIndex(0)==0) {  //best kink
      if (mothers[indexes[i]].GetNumberOfClusters()>20 && daughters[indexes[i]].GetNumberOfClusters()>20 && (mothers[indexes[i]].GetNumberOfClusters()+daughters[indexes[i]].GetNumberOfClusters())>100){
	new (ktrack0) IlcDCHseed(mothers[indexes[i]]);
	new (ktrack1) IlcDCHseed(daughters[indexes[i]]);
      }
    }
    //
    ktrack0->SetKinkIndex(usage[index0],-(index+1));
    ktrack1->SetKinkIndex(usage[index1], (index+1));
    usage[index0]++;
    usage[index1]++;
  }
  //
  // Remove tracks corresponding to shared kink's
  //
  for (Int_t i=0;i<nentries;i++){
    IlcDCHseed * track0 = (IlcDCHseed*)array->At(i);
    if (!track0) continue;
    if (track0->GetKinkIndex(0)!=0) continue;
    if (shared[i]) delete array->RemoveAt(i);
  }

  //
  //
  RemoveUsed2(array,0.5,0.4,fRecPar->GetMinNOfClusetrs());
  UnsignClusters();
  for (Int_t i=0;i<nentries;i++){
    IlcDCHseed * track0 = (IlcDCHseed*)array->At(i);
    if (!track0) continue;
    track0->CookdEdx(0.02,0.6);
    track0->CookPID();
  }
  //
  for (Int_t i=0;i<nentries;i++){
    IlcDCHseed * track0 = (IlcDCHseed*)array->At(i);
    if (!track0) continue;
    if (track0->GetPt()<1.4) continue;
    //remove double high momenta tracks - overlapped with kink candidates
    Int_t shared=0;
    Int_t all   =0;
    for (Int_t icl=track0->GetFirstPoint();icl<track0->GetLastPoint(); icl++){
      if (track0->GetClusterPointer(icl)!=0){
	all++;
	if (track0->GetClusterPointer(icl)->IsUsed(10)) shared++;
      }
    }
    if (Float_t(shared+1)/Float_t(all+1)>0.5) {
      Info("FindKinks",Form("Delete tracks shared points=%i, all=%i, ratio=%f",shared,all,Double_t(shared+1)/Double_t(all+1)));
      delete array->RemoveAt(i);
      continue;
    }
    //
    if (track0->GetKinkIndex(0)!=0) continue;
    if (track0->GetNumberOfClusters()<80) continue;

    IlcDCHseed *pmother = new IlcDCHseed();
    IlcDCHseed *pdaughter = new IlcDCHseed();
    IlcKink *pkink = new IlcKink;

    IlcDCHseed & mother = *pmother;
    IlcDCHseed & daughter = *pdaughter;
    IlcKink & kink = *pkink;
    if (CheckKinkPoint(track0,mother,daughter, kink)){
      if (mother.GetNumberOfClusters()<30||daughter.GetNumberOfClusters()<20) {
	delete pmother;
	delete pdaughter;
	delete pkink;
	continue;  //too short tracks
      }
      if (mother.GetPt()<1.4) {
	delete pmother;
	delete pdaughter;
	delete pkink;
	continue;
      }
      Int_t layer0= kink.GetTPCRow0();
      if (kink.GetDistance()>0.5 || kink.GetR()<fParam->GetInnerRadius()+fRecPar->GetDistFromEdge() || kink.GetR()>fParam->GetOuterRadius()-fRecPar->GetDistFromEdge()) {
//      cout<<"buttata 5\n";
	delete pmother;
	delete pdaughter;
	delete pkink;
	continue;
      }
      //
      Int_t index = esd->AddKink(&kink);      
      mother.SetKinkIndex(0,-(index+1));
      daughter.SetKinkIndex(0,index+1);
      if (mother.GetNumberOfClusters()>50) {
	delete array->RemoveAt(i);
	array->AddAt(new IlcDCHseed(mother),i);
      }
      else{
	array->AddLast(new IlcDCHseed(mother));
      }
      array->AddLast(new IlcDCHseed(daughter));      
      for (Int_t icl=0;icl<layer0;icl++) {
	if (mother.GetClusterPointer(icl)) mother.GetClusterPointer(icl)->Use(20);
      }
      //
      for (Int_t icl=layer0;icl<nlayers;icl++) {
	if (daughter.GetClusterPointer(icl)) daughter.GetClusterPointer(icl)->Use(20);
      }
      //
    }
    delete pmother;
    delete pdaughter;
    delete pkink;
  }

  delete [] daughters;
  delete [] mothers;
  //
  //
  delete [] dca;
  delete []circular;
  delete []shared;
  delete []quilcty;
  delete []indexes;
  //
  delete kink;
  delete[]fim;
  delete[] zm;
  delete[] z0;
  delete [] usage;
  delete[] alpha;
  delete[] nclusters;
  delete[] sign;
  delete[] helixes;
  kinks->Delete();
  delete kinks;

  Info("FindKinks","Ncandidates:\t nkinks=%d\t ntracks=%d\t ncanditates=%d\t nall=%d\n",esd->GetNumberOfKinks(),ntracks,ncandidates,nall);
  timer.Print();
}
*/
/*
void  IlcDCHtracker::FindV0s(TObjArray * array, IlcESD *esd)
{
  //
  //  find V0s
  //
  //
  TObjArray *dchv0s      = new TObjArray(100000);
  Int_t     nentries     = array->GetEntriesFast();
  IlcHelix *helixes      = new IlcHelix[nentries];
  Int_t    *sign         = new Int_t[nentries];
  Float_t  *alpha        = new Float_t[nentries];
  Float_t  *z0           = new Float_t[nentries]; 
  Float_t  *dca          = new Float_t[nentries];
  Float_t  *sdcar        = new Float_t[nentries];
  Float_t  *cdcar        = new Float_t[nentries];
  Float_t  *pulldcar     = new Float_t[nentries];
  Float_t  *pulldcaz     = new Float_t[nentries];
  Float_t  *pulldca      = new Float_t[nentries];
  Bool_t   *isPrim       = new Bool_t[nentries];  
  const IlcESDVertex * primvertex = esd->GetVertex();
  Double_t             zvertex = primvertex->GetZv(); 
  //
  //  nentries = array->GetEntriesFast();
  //
  for (Int_t i=0;i<nentries;i++){
    sign[i]=0;
    isPrim[i]=0;
    IlcDCHseed* track = (IlcDCHseed*)array->At(i);    
    if (!track) continue;
    track->GetV0Indexes()[0] = 0;  //rest v0 indexes
    track->GetV0Indexes()[1] = 0;  //rest v0 indexes
    track->GetV0Indexes()[2] = 0;  //rest v0 indexes
    //
    alpha[i] = track->GetAlpha();
    new (&helixes[i]) IlcHelix(*track);
    Double_t xyz[3];
    helixes[i].Evaluate(0,xyz);
    sign[i] = (track->GetC()>0) ? -1:1;
    Double_t x,y,z;
    x=160;
    z0[i]=1000;
    if (track->GetProlongation(0,y,z))  z0[i] = z;
    dca[i] = track->GetD(0,0);
    // 
    // dca error parrameterezation + pulls
    //
    sdcar[i]      = TMath::Sqrt(0.150*0.150+(100*track->GetC())*(100*track->GetC()));
    if (TMath::Abs(track->GetTgl())>1) sdcar[i]*=2.5;
    cdcar[i]      = TMath::Exp((TMath::Abs(track->GetC())-0.0106)*525.3);
    pulldcar[i]   = (dca[i]-cdcar[i])/sdcar[i];
    pulldcaz[i]   = (z0[i]-zvertex)/sdcar[i];
    pulldca[i]    = TMath::Sqrt(pulldcar[i]*pulldcar[i]+pulldcaz[i]*pulldcaz[i]);
    if (track->DCHrPID(1)+track->DCHrPID(2)+track->DCHrPID(3)>0.5) {
      if (pulldca[i]<3.) isPrim[i]=kTRUE;  //pion, muon and Kaon  3 sigma cut
    }
    if (track->DCHrPID(4)>0.5) {
      if (pulldca[i]<0.5) isPrim[i]=kTRUE;  //proton 0.5 sigma cut
    }
    if (track->DCHrPID(0)>0.4) {
      isPrim[i]=kFALSE;  //electron no  sigma cut
    }
  }
  //
  //
  TStopwatch timer;
  timer.Start();
  Int_t ncandidates =0;
  Int_t nall =0;
  Int_t ntracks=0; 
  Double_t phase[2][2],radius[2];
  //
  //  Finf V0s loop
  // 
  //
  // //  
  TTreeSRedirector &cstream = *fDebugStreamer; 
  Float_t fprimvertex[3]={GetX(),GetY(),GetZ()};
  IlcV0 vertex; 
  Double_t cradius0 = fRecPar->GetMinRofV0()*fRecPar->GetMinRofV0();
  Double_t cradius1 = fRecPar->GetMaxRofV0()*fRecPar->GetMaxRofV0();
  Double_t cdist1=fRecPar->GetCutOnDistanceInV0();
  Double_t cdist2=fRecPar->GetCutOnDistanceInV0();
  Double_t cpointAngle = 0.95;
  //
  Double_t delta[2]={10000,10000};
  for (Int_t i =0;i<nentries;i++){
    if (sign[i]==0) continue;
    IlcDCHseed * track0 = (IlcDCHseed*)array->At(i);
    if (!track0) continue;
    if (IlcDCHReconstructor::StreamLevel()>0){
      cstream<<"Tracks"<<
	"Tr0.="<<track0<<
	"dca="<<dca[i]<<
	"z0="<<z0[i]<<
	"zvertex="<<zvertex<<
	"sdcar0="<<sdcar[i]<<
	"cdcar0="<<cdcar[i]<<
	"pulldcar0="<<pulldcar[i]<<
	"pulldcaz0="<<pulldcaz[i]<<
	"pulldca0="<<pulldca[i]<<
	"isPrim="<<isPrim[i]<<
	"\n";
    }
    //
    if (track0->Get1Pt()<0) continue;
    if (track0->GetKinkIndex(0)>0||isPrim[i]) continue;   //daughter kink
    //
    if (TMath::Abs(helixes[i].GetHelix(4))<0.000000001) continue;
    ntracks++;
    // debug output
    
    
    for (Int_t j =0;j<nentries;j++){
      IlcDCHseed * track1 = (IlcDCHseed*)array->At(j);
      if (!track1) continue;
      if (track1->GetKinkIndex(0)>0 || isPrim[j]) continue; //daughter kink
      if (sign[j]*sign[i]>0) continue; 
      if (TMath::Abs(helixes[j].GetHelix(4))<0.000001) continue;
      if (track0->GetCircular()+track1->GetCircular()>1) continue;    //circling -returning track
      nall++;
      //
      // DCA to prim vertex cut
      //
      //
      delta[0]=10000;
      delta[1]=10000;

      Int_t npoints = helixes[i].GetRPHIintersections(helixes[j], phase, radius,cdist2);
      if (npoints<1) continue;
      Int_t iclosest=0;
      // cuts on radius      
      if (npoints==1){
	if (radius[0]<cradius0||radius[0]>cradius1) continue;
	helixes[i].LinearDCA(helixes[j],phase[0][0],phase[0][1],radius[0],delta[0]);
	if (delta[0]>cdist1) continue;
      }
      else{
	if (TMath::Max(radius[0],radius[1])<cradius0|| TMath::Min(radius[0],radius[1])>cradius1) continue;	
	helixes[i].LinearDCA(helixes[j],phase[0][0],phase[0][1],radius[0],delta[0]);	
	helixes[i].LinearDCA(helixes[j],phase[1][0],phase[1][1],radius[1],delta[1]);
	if (delta[1]<delta[0]) iclosest=1;
	if (delta[iclosest]>cdist1) continue;
      }
      helixes[i].ParabolicDCA(helixes[j],phase[iclosest][0],phase[iclosest][1],radius[iclosest],delta[iclosest]);
      if (radius[iclosest]<cradius0 || radius[iclosest]>cradius1 || delta[iclosest]>cdist1) continue;
      //
      Double_t pointAngle = helixes[i].GetPointAngle(helixes[j],phase[iclosest],fprimvertex);
      if (pointAngle<cpointAngle) continue;
      //
      Bool_t isGamma = kFALSE;
      vertex.SetP(*track0); //track0 - plus
      vertex.SetM(*track1); //track1 - minus
      vertex.Update(fprimvertex);
      if (track0->DCHrPID(0)>0.3&&track1->DCHrPID(0)>0.3&&vertex.GetAnglep()[2]<0.15) isGamma=kTRUE;              // gamma conversion candidate
      Double_t pointAngle2 = vertex.GetPointAngle();
      //continue;
      if (vertex.GetPointAngle()<cpointAngle && (!isGamma)) continue; // point angle cut
      if (vertex.GetDist2()>2&&(!isGamma)) continue;         // point angle cut
      Float_t sigmae     = 0.15*0.15;
      if (vertex.GetRr()<fParam->GetInnerRadius()) 
	sigmae += (sdcar[i]*sdcar[i]+sdcar[j]*sdcar[j])*(1.-vertex.GetRr()/fParam->GetInnerRadius())*(1.-vertex.GetRr()/fParam->GetInnerRadius());
      sigmae+= TMath::Sqrt(sigmae);
      if (vertex.GetDist2()/sigmae>3.&&(!isGamma)) continue; 
      Float_t densb0=0,densb1=0,densa0=0,densa1=0;
      Int_t layer0 = GetLayerNumber(vertex.GetRr(),vertex.GetXr(2));
      if (layer0>fRecPar->GetMinNOfClusetrs()){
	if (vertex.GetDist2()>0.2) continue;             
	densb0     = track0->Density2(0,layer0-5);          
	densb1     = track1->Density2(0,layer0-5);         
	if (densb0>0.3|| densb1>0.3) continue;            //clusters before vertex
	densa0     = track0->Density2(layer0+5,layer0+40);    
	densa1     = track1->Density2(layer0+5,layer0+40);    
	if ((densa0<0.4|| densa1<0.4)&&(!isGamma)) continue;            //missing clusters after vertex
      }
      else{
	densa0     = track0->Density2(0,40);  //cluster density
	densa1     = track1->Density2(0,40);  //cluster density
	if ((vertex.GetRr()<fParam->GetInnerRadius()&&densa0+densa1<1.)&&(!isGamma)) continue;
      }
      vertex.SetLab(0,track0->GetLabel());
      vertex.SetLab(1,track1->GetLabel());
      vertex.SetChi2After((densa0+densa1)*0.5);
      vertex.SetChi2Before((densb0+densb1)*0.5);
      vertex.SetIndex(0,i);
      vertex.SetIndex(1,j);
      vertex.SetStatus(1); // DCH v0 candidate
      vertex.SetRp(track0->DCHrPIDs());
      vertex.SetRm(track1->DCHrPIDs());
      dchv0s->AddLast(new IlcESDv0(vertex));      
      ncandidates++;
      {
	Int_t eventNr = esd->GetEventNumber();
	Double_t radiusm= (delta[0]<delta[1])? TMath::Sqrt(radius[0]):TMath::Sqrt(radius[1]);  
	Double_t deltam= (delta[0]<delta[1])? TMath::Sqrt(delta[0]):TMath::Sqrt(delta[1]);  
	if (IlcDCHReconstructor::StreamLevel()>0) {
	  Int_t lab0=track0->GetLabel();
	  Int_t lab1=track1->GetLabel();
	  Char_t c0=track0->GetCircular();
	  Char_t c1=track1->GetCircular();
	  cstream<<"V0"<<
	  "Event="<<eventNr<<
	  "vertex.="<<&vertex<<
	  "Tr0.="<<track0<<
	  "lab0="<<lab0<<
	  "Helix0.="<<&helixes[i]<<	
	  "Tr1.="<<track1<<
	  "lab1="<<lab1<<
	  "Helix1.="<<&helixes[j]<<
	  "pointAngle="<<pointAngle<<
	  "pointAngle2="<<pointAngle2<<
	  "dca0="<<dca[i]<<
	  "dca1="<<dca[j]<<
	  "z0="<<z0[i]<<
	  "z1="<<z0[j]<<
	  "zvertex="<<zvertex<<
	  "circular0="<<c0<<
	  "circular1="<<c1<<
	  "npoints="<<npoints<<
	  "radius0="<<radius[0]<<
	  "delta0="<<delta[0]<<
	  "radius1="<<radius[1]<<
	  "delta1="<<delta[1]<<
	  "radiusm="<<radiusm<<
	  "deltam="<<deltam<<
	  "sdcar0="<<sdcar[i]<<
	  "sdcar1="<<sdcar[j]<<
	  "cdcar0="<<cdcar[i]<<
	  "cdcar1="<<cdcar[j]<<
	  "pulldcar0="<<pulldcar[i]<<
	  "pulldcar1="<<pulldcar[j]<<
	  "pulldcaz0="<<pulldcaz[i]<<
	  "pulldcaz1="<<pulldcaz[j]<<
	  "pulldca0="<<pulldca[i]<<
	  "pulldca1="<<pulldca[j]<<
	  "densb0="<<densb0<<
	  "densb1="<<densb1<<
	  "densa0="<<densa0<<
	  "densa1="<<densa1<<
	  "sigmae="<<sigmae<<
	    "\n";}
      }
    }
  }    
  Float_t *quilcty = new Float_t[ncandidates];
  Int_t *indexes = new Int_t[ncandidates];
  Int_t naccepted =0;
  for (Int_t i=0;i<ncandidates;i++){
    quilcty[i]     = 0; 
    IlcESDv0 *v0 = (IlcESDv0*)dchv0s->At(i);
    quilcty[i]     = 1./(1.00001-v0->GetPointAngle());   //base point angle
    // quilcty[i]    /= (0.5+v0->GetDist2());  
    // quilcty[i]    *= v0->GetChi2After();               //density factor
    Double_t minpulldca = TMath::Min(2.+pulldca[v0->GetIndex(0)],(2.+pulldca[v0->GetIndex(1)]) );     //pull
    Int_t index0 = v0->GetIndex(0);
    Int_t index1 = v0->GetIndex(1);
    IlcDCHseed * track0 = (IlcDCHseed*)array->At(index0);
    IlcDCHseed * track1 = (IlcDCHseed*)array->At(index1);
    if (track0->DCHrPID(0)>0.3&&track1->DCHrPID(0)>0.3&&v0->GetAnglep()[2]<0.15) quilcty[i]+=1000000;              // gamma conversion candidate
    if (track0->DCHrPID(4)>0.9||(track1->DCHrPID(4)>0.9&&minpulldca>4)) quilcty[i]*=10;    // lambda candidate candidate
  }

  TMath::Sort(ncandidates,quilcty,indexes,kTRUE);
  //
  //
  for (Int_t i=0;i<ncandidates;i++){
    IlcESDv0 * v0 = (IlcESDv0*)dchv0s->At(indexes[i]);
    if (!v0) continue;
    Int_t index0 = v0->GetIndex(0);
    Int_t index1 = v0->GetIndex(1);
    IlcDCHseed * track0 = (IlcDCHseed*)array->At(index0);
    IlcDCHseed * track1 = (IlcDCHseed*)array->At(index1);
    if (!track0||!track1) {
      printf("Bug\n");
      continue;
    }
    Bool_t accept =kTRUE;  //default accept
    Int_t *v0indexes0 = track0->GetV0Indexes();
    Int_t *v0indexes1 = track1->GetV0Indexes();
    //
    Int_t order0 = (v0indexes0[0]!=0) ? 1:0;
    Int_t order1 = (v0indexes1[0]!=0) ? 1:0;    
    if (v0indexes0[1]!=0) order0 =2;
    if (v0indexes1[1]!=0) order1 =2;      
    //
    if (v0indexes0[2]!=0) {order0=3; accept=kFALSE;}
    if (v0indexes0[2]!=0) {order1=3; accept=kFALSE;}
    //
    IlcESDv0 * v02 = v0;
    if (accept){
      v0->SetOrder(0,order0);
      v0->SetOrder(1,order1);
      v0->SetOrder(1,order0+order1);     
      v0->SetOnFlyStatus(kTRUE);
      Int_t index = esd->AddV0(v0);
      v02 = esd->GetV0(index);
      v0indexes0[order0]=index;
      v0indexes1[order1]=index;
      naccepted++;
    }
    {
      Int_t eventNr = esd->GetEventNumber();
      if (IlcDCHReconstructor::StreamLevel()>0) {
	Int_t lab0=track0->GetLabel();
	Int_t lab1=track1->GetLabel();
	cstream<<"V02"<<
	"Event="<<eventNr<<
	"vertex.="<<v0<<	
	"vertex2.="<<v02<<
	"Tr0.="<<track0<<
	"lab0="<<lab0<<
	"Tr1.="<<track1<<
	"lab1="<<lab1<<
	"dca0="<<dca[index0]<<
	"dca1="<<dca[index1]<<
	"order0="<<order0<<
	"order1="<<order1<<
	"accept="<<accept<<
	"quilcty="<<quilcty[i]<<
	"pulldca0="<<pulldca[index0]<<
	"pulldca1="<<pulldca[index1]<<
	"index="<<i<<
	  "\n";}
    }
  }    


  //
  //
  delete []quilcty;
  delete []indexes;
//
  delete [] isPrim;
  delete [] pulldca;
  delete [] pulldcaz;
  delete [] pulldcar;
  delete [] cdcar;
  delete [] sdcar;
  delete [] dca;
  //
  delete[] z0;
  delete[] alpha;
  delete[] sign;
  delete[] helixes;
  printf("DCH V0 finder : naccepted\t%d\tncandidates\t%d\tntracks\t%d\tnall\t%d\n",naccepted,ncandidates,ntracks,nall);
  timer.Print();
}
*/
/*
Int_t IlcDCHtracker::RefitKink(IlcDCHseed &mother, IlcDCHseed &daughter, IlcESDkink &knk)
{
  //
  // refit kink towards to the vertex
  //
  //
  int nlayers=fSector.GetNLayers();
  IlcKink &kink=(IlcKink &)knk;

  Int_t layer0 = GetLayerNumber(kink.GetR(),kink.GetPosition()[2]);
  FollowProlongation(mother,0);
  mother.Reset(kFALSE);
  //
  FollowProlongation(daughter,layer0);
  daughter.Reset(kFALSE);
  FollowBackProlongation(daughter,nlayers-1);
  daughter.Reset(kFALSE);
  Int_t first = TMath::Max(layer0-fRecPar->GetMinNOfClusetrs(),fRecPar->GetMinNOfClusetrs()); 
  Int_t last  = TMath::Min(layer0+fRecPar->GetMinNOfClusetrs(),nlayers-fRecPar->GetMinNOfClusetrs());
  //
  const Int_t kNdiv =5;
  IlcDCHseed  param0[kNdiv];  // parameters along the track
  IlcDCHseed  param1[kNdiv];  // parameters along the track
  IlcKink     kinks[kNdiv];   // corresponding kink parameters
  //
  Int_t layers[kNdiv];
  for (Int_t ilayer=0; ilayer<kNdiv;ilayer++){
    layers[ilayer] = first +((last-first)*ilayer)/(kNdiv-1);
  }
  // store parameters along the track
  //
  for (Int_t ilayer=0;ilayer<kNdiv;ilayer++){
    FollowBackProlongation(mother, layers[ilayer]);
    FollowProlongation(daughter,layers[kNdiv-1-ilayer]);       
    new(&param0[ilayer])     IlcDCHseed(mother);
    new(&param1[kNdiv-1-ilayer])   IlcDCHseed(daughter);
  }
  //
  // define kinks 
  for (Int_t ilayer=0; ilayer<kNdiv-1;ilayer++){
    if (param0[ilayer].GetNumberOfClusters()<kNdiv||param1[ilayer].GetNumberOfClusters()<kNdiv) continue;
    kinks[ilayer].SetMother(param0[ilayer]);
    kinks[ilayer].SetDaughter(param1[ilayer]);
    kinks[ilayer].Update();
  }
  //
  // choose kink with best "quilcty"
  Int_t index =-1;
  Double_t mindist = 10000;
  for (Int_t ilayer=0;ilayer<kNdiv;ilayer++){
    if (param0[ilayer].GetNumberOfClusters()<fRecPar->GetMinNOfClusetrs()||param1[ilayer].GetNumberOfClusters()<fRecPar->GetMinNOfClusetrs()) continue;
    if (TMath::Abs(kinks[ilayer].GetR())>fParam->GetOuterRadius()-fRecPar->GetDistFromEdge()) continue; //240.
    if (TMath::Abs(kinks[ilayer].GetR())<fParam->GetInnerRadius()+fRecPar->GetDistFromEdge()) continue; //100.
    //
    Float_t normdist = TMath::Abs(param0[ilayer].GetX()-kinks[ilayer].GetR())*(0.1+kink.GetDistance());
    normdist/= (param0[ilayer].GetNumberOfClusters()+param1[ilayer].GetNumberOfClusters()+2*fRecPar->GetMinNOfClusetrs());
    if (normdist < mindist){
      mindist = normdist;
      index = ilayer;
    }
  }
  //
  if (index==-1) return 0;
  //
  //
  param0[index].Reset(kTRUE);
  FollowProlongation(param0[index],0);
  //
  new (&mother) IlcDCHseed(param0[index]);
  new (&daughter) IlcDCHseed(param1[index]);  // daughter in vertex
  //
  kink.SetMother(mother);
  kink.SetDaughter(daughter);
  kink.Update();
  kink.SetTPCRow0(GetLayerNumber(kink.GetR(),kink.GetPosition()[2]));
  kink.SetTPCncls(param0[index].GetNumberOfClusters(),0);
  kink.SetTPCncls(param1[index].GetNumberOfClusters(),1);
  kink.SetLabel(CookLabel(&mother,0.4, 0,kink.GetTPCRow0()),0);
  kink.SetLabel(CookLabel(&daughter,0.4, kink.GetTPCRow0(),nlayers),1);
  mother.SetLabel(kink.GetLabel(0));
  daughter.SetLabel(kink.GetLabel(1));

  return 1;
}
*/
/*
void IlcDCHtracker::UpdateKinkQuilctyM(IlcDCHseed * seed){
  //
  // update Kink quilcty information for mother after back propagation
  //
  int nlayers=fSector.GetNLayers();
  if (seed->GetKinkIndex(0)>=0) return; 
  for (Int_t ikink=0;ikink<3;ikink++){
    Int_t index = seed->GetKinkIndex(ikink);
    if (index>=0) break;
    index = TMath::Abs(index)-1;
    IlcESDkink * kink = fEvent->GetKink(index);
    //kink->fDCHdensity2[0][0]=-1;
    //kink->fDCHdensity2[0][1]=-1;
    kink->SetTPCDensity2(-1,0,0);
    kink->SetTPCDensity2(1,0,1);
    //
    Int_t layer0 = kink->GetTPCRow0() - 2 - Int_t( 0.5/ (0.05+kink->GetAngle(2)));
    if (layer0<fRecPar->GetMinNOfClusetrs()) layer0=fRecPar->GetMinNOfClusetrs();
    //
    Int_t layer1 = kink->GetTPCRow0() + 2 +  Int_t( 0.5/ (0.05+kink->GetAngle(2)));
    if (layer1>nlayers-fRecPar->GetMinNOfClusetrs()) layer1=nlayers-fRecPar->GetMinNOfClusetrs();
    //
    Int_t found,foundable,shared;
    seed->GetClusterStatistic(0,layer0, found, foundable,shared,kFALSE);
    if (foundable>5)   kink->SetTPCDensity2(Float_t(found)/Float_t(foundable),0,0);
    seed->GetClusterStatistic(layer1,nlayers-1, found, foundable,shared,kFALSE);
    if (foundable>5)   kink->SetTPCDensity2(Float_t(found)/Float_t(foundable),0,1);
  }
    
}

void IlcDCHtracker::UpdateKinkQuilctyD(IlcDCHseed * seed){
  //
  // update Kink quilcty information for daughter after refit
  //
  int nlayers=fSector.GetNLayers();
  if (seed->GetKinkIndex(0)<=0) return; 
  for (Int_t ikink=0;ikink<3;ikink++){
    Int_t index = seed->GetKinkIndex(ikink);
    if (index<=0) break;
    index = TMath::Abs(index)-1;
    IlcESDkink * kink = fEvent->GetKink(index);
    kink->SetTPCDensity2(-1,1,0);
    kink->SetTPCDensity2(-1,1,1);
    //
    Int_t layer0 = kink->GetTPCRow0() -2 - Int_t( 0.5/ (0.05+kink->GetAngle(2)));
    if (layer0<fRecPar->GetMinNOfClusetrs()) layer0=fRecPar->GetMinNOfClusetrs();
    //
    Int_t layer1 = kink->GetTPCRow0() +2 +  Int_t( 0.5/ (0.05+kink->GetAngle(2)));
    if (layer1>nlayers-fRecPar->GetMinNOfClusetrs()) layer1=nlayers-fRecPar->GetMinNOfClusetrs();
    //
    Int_t found,foundable,shared;
    seed->GetClusterStatistic(0,layer0, found, foundable,shared,kFALSE);
    if (foundable>5)   kink->SetTPCDensity2(Float_t(found)/Float_t(foundable),1,0);
    seed->GetClusterStatistic(layer1,nlayers-1, found, foundable,shared,kFALSE);
    if (foundable>5)   kink->SetTPCDensity2(Float_t(found)/Float_t(foundable),1,1);
  }
    
}


Int_t  IlcDCHtracker::CheckKinkPoint(IlcDCHseed*seed,IlcDCHseed &mother, IlcDCHseed &daughter, IlcESDkink &knk)
{
  //
  // check kink point for given track
  // if return value=0 kink point not found
  // otherwise seed0 correspond to mother particle
  //           seed1 correspond to daughter particle
  //           kink  parameter of kink point
  int nlayers=fSector.GetNLayers();
  IlcKink &kink=(IlcKink &)knk;

  Int_t middlelayer = (seed->GetFirstPoint()+seed->GetLastPoint())/2;
  Int_t first = seed->GetFirstPoint(); 
  Int_t last  = seed->GetLastPoint();
  const int nchecks=20;
  if (last-first< nchecks) return 0;          // shortest length - 2*30 = 60 pad-layers

  
  IlcDCHseed *seed1 = ReSeed(seed,middlelayer+ nchecks, kTRUE);  //middle of chamber
  if (!seed1) return 0;
  FollowProlongation(*seed1,seed->GetLastPoint()-nchecks);
  seed1->Reset(kTRUE);
  FollowProlongation(*seed1,nlayers-1);
  seed1->Reset(kTRUE);  
  last = seed1->GetLastPoint();
  //
  IlcDCHseed *seed0 = new IlcDCHseed(*seed);
  seed0->Reset(kFALSE);
  seed0->Reset();
  //
  IlcDCHseed  param0[nchecks];  // parameters along the track
  IlcDCHseed  param1[nchecks];  // parameters along the track
  IlcKink            kinks[nchecks];   // corresponding kink parameters

  Int_t layers[nchecks];
  for (Int_t ilayer=0; ilayer<nchecks;ilayer++){
    layers[ilayer] = first +((last-first)*ilayer)/(nchecks-1);
  }
  // store parameters along the track
  //
  for (Int_t ilayer=0;ilayer<nchecks;ilayer++){
    FollowBackProlongation(*seed0, layers[ilayer]);
    FollowProlongation(*seed1,layers[nchecks-1-ilayer]);       
    new(&param0[ilayer])     IlcDCHseed(*seed0);
    new(&param1[nchecks-1-ilayer])   IlcDCHseed(*seed1);
  }
  //
  // define kinks 
  for (Int_t ilayer=0; ilayer<nchecks-1;ilayer++){
    kinks[ilayer].SetMother(param0[ilayer]);
    kinks[ilayer].SetDaughter(param1[ilayer]);
    kinks[ilayer].Update();
  }
  //
  // choose kink with biggest change of angle
  Int_t index =-1;
  Double_t maxchange= 0;
  for (Int_t ilayer=1;ilayer<nchecks-1;ilayer++){
    if (TMath::Abs(kinks[ilayer].GetR())>fParam->GetOuterRadius()-fRecPar->GetDistFromEdge()) continue;
    if (TMath::Abs(kinks[ilayer].GetR())<fParam->GetInnerRadius()+fRecPar->GetDistFromEdge()) continue;
    Float_t quilcty = TMath::Abs(kinks[ilayer].GetAngle(2))/(fRecPar->GetScaleForQuilcty()+TMath::Abs(kinks[ilayer].GetR()-param0[ilayer].GetR()));//dddddddddddddd
    if ( quilcty > maxchange){
      maxchange = quilcty;
      index = ilayer;
      //
    }
  }
  delete seed0;
  delete seed1;
  if (index<0) return 0;
  //
  Int_t layer0    = GetLayerNumber(kinks[index].GetR(),kinks[index].GetPosition()[2]);   //layer 0 estimate
  seed0 = new IlcDCHseed(param0[index]);
  seed1 = new IlcDCHseed(param1[index]);
  seed0->Reset(kFALSE);
  seed1->Reset(kFALSE);
  seed0->ResetCovariance(10.);
  seed1->ResetCovariance(10.);
  FollowProlongation(*seed0,0);
  FollowBackProlongation(*seed1,nlayers-1);
  new (&mother) IlcDCHseed(*seed0);  // backup mother at position 0
  seed0->Reset(kFALSE);  
  seed1->Reset(kFALSE);
  seed0->ResetCovariance(10.);
  seed1->ResetCovariance(10.);
  //
  first = TMath::Max(layer0-nchecks,0);
  last  = TMath::Min(layer0+nchecks,nlayers-1);
  //
  for (Int_t ilayer=0; ilayer<nchecks;ilayer++){
    layers[ilayer] = first +((last-first)*ilayer)/(nchecks-1);
  }
  // store parameters along the track
  //
  for (Int_t ilayer=0;ilayer<nchecks;ilayer++){
    FollowBackProlongation(*seed0, layers[ilayer]);
    FollowProlongation(*seed1,layers[nchecks-1-ilayer]);       
    new(&param0[ilayer])     IlcDCHseed(*seed0);
    new(&param1[nchecks-1-ilayer])   IlcDCHseed(*seed1);
  }
  //
  // define kinks 
  for (Int_t ilayer=0; ilayer<nchecks-1;ilayer++){
    kinks[ilayer].SetMother(param0[ilayer]);
    kinks[ilayer].SetDaughter(param1[ilayer]);
    //    param0[ilayer].Dump();
    //param1[ilayer].Dump();
    kinks[ilayer].Update();
  }
  //
  // choose kink with biggest change of angle
  index =-1;
  maxchange= 0;
  for (Int_t ilayer=0;ilayer<nchecks;ilayer++){
    if (TMath::Abs(kinks[ilayer].GetR())>fParam->GetOuterRadius()-fRecPar->GetDistFromEdge()) continue;
    if (TMath::Abs(kinks[ilayer].GetR())<fParam->GetInnerRadius()+fRecPar->GetDistFromEdge()) continue;
    Float_t quilcty = TMath::Abs(kinks[ilayer].GetAngle(2))/(fRecPar->GetScaleForQuilcty()+TMath::Abs(kinks[ilayer].GetR()-param0[ilayer].GetR()));//ddddddddddddd
    if ( quilcty > maxchange){
      maxchange = quilcty;
      index = ilayer;
      //
    }
  }
  //
  //
  if (index==-1 || param0[index].GetNumberOfClusters()+param1[index].GetNumberOfClusters()<100){
    delete seed0;
    delete seed1;
    return 0;
  }
  //  Float_t anglesigma = TMath::Sqrt(param0[index].fC22+param0[index].fC33+param1[index].fC22+param1[index].fC33);
  
  kink.SetMother(param0[index]);
  kink.SetDaughter(param1[index]);
  kink.Update();
  layer0    = GetLayerNumber(kink.GetR(),kink.GetPosition()[2]);   
  kink.SetTPCRow0(layer0);
  kink.SetLabel(CookLabel(seed0,0.5,0,layer0),0);
  kink.SetLabel(CookLabel(seed1,0.5,layer0,nlayers-1),1);
  kink.SetIndex(-10,0);
  kink.SetIndex(int(param0[index].GetNumberOfClusters()+param1[index].GetNumberOfClusters()),1);
  kink.SetTPCncls(param0[index].GetNumberOfClusters(),0);
  kink.SetTPCncls(param1[index].GetNumberOfClusters(),1);
  //
  //
  //  new (&mother) IlcDCHseed(param0[index]);
  new (&daughter) IlcDCHseed(param1[index]);
  daughter.SetLabel(kink.GetLabel(1));  
  param0[index].Reset(kTRUE);
  FollowProlongation(param0[index],0);  
  new (&mother) IlcDCHseed(param0[index]);
  mother.SetLabel(kink.GetLabel(0));
  delete seed0;
  delete seed1;
  //
  return 1;
}

*/


IlcDCHseed*  IlcDCHtracker::ReSeed(IlcDCHseed *t)
{
  //
  // reseed - refit -  track
  //
  Int_t first = t->GetFirstPoint();
  //
  //
  IlcDCHseed * seed = MakeSeed(t,0.1,0.5,0.9);
  FollowBackProlongation(*t,fSector.GetNLayers()-1);
  t->Reset(kFALSE);
  FollowProlongation(*t,first);
  return seed;
}







//_____________________________________________________________________________
Int_t IlcDCHtracker::ReadSeeds(const TFile *inp) {
  //-----------------------------------------------------------------
  // This function reades track seeds.
  //-----------------------------------------------------------------
  TDirectory *savedir=gDirectory; 

  TFile *in=(TFile*)inp;
  if (!in->IsOpen()) {
     cerr<<"IlcDCHtracker::ReadSeeds(): input file is not open !\n";
     return 1;
  }

  in->cd();
  TTree *seedTree=(TTree*)in->Get("Seeds");
  if (!seedTree) {
     cerr<<"IlcDCHtracker::ReadSeeds(): ";
     cerr<<"can't get a tree with track seeds !\n";
     return 2;
  }
  IlcDCHtrack *seed=new IlcDCHtrack; 
  seedTree->SetBranchAddress("tracks",&seed);
  
  if (fSeeds==0) fSeeds=new TObjArray(15000);

  Int_t n=(Int_t)seedTree->GetEntries();
  for (Int_t i=0; i<n; i++) {
     seedTree->GetEvent(i);
     fSeeds->AddLast(new IlcDCHseed(*seed));
  }
  
  delete seed;
  delete seedTree; 
  savedir->cd();
  return 0;
}
/*
Int_t IlcDCHtracker::Clusters2Tracks (IlcESD *esd)
{
  //
  if (fSeeds) DeleteSeeds();
  fEvent = esd;
  Clusters2Tracks();
  if (!fSeeds) return 1;
  FillESD(fSeeds);
  return 0;
  //
}
*/

//_____________________________________________________________________________
Int_t IlcDCHtracker::Clusters2Tracks() {
  //-----------------------------------------------------------------
  // This is a track finder.
  //-----------------------------------------------------------------
  TDirectory *savedir=gDirectory; 
  TStopwatch timer;

  fIteration = 0;
  fSeeds = Tracking();
  fIteration = 1;
  for (Int_t i=0; i<fSeeds->GetEntriesFast(); i++) {
    IlcDCHseed *pt=(IlcDCHseed*)fSeeds->UncheckedAt(i);    
    if (!pt) continue;   
    pt->Reset(kFALSE);
    pt->SetRemoval(0);
  }
  PropagateBack(fSeeds); 
  fIteration = 2;
  for (Int_t i=0; i<fSeeds->GetEntriesFast(); i++) {
    IlcDCHseed *pt=(IlcDCHseed*)fSeeds->UncheckedAt(i);    
    if (!pt) continue;   
    if(pt->GetTgl()<0) pt->ChangeDirection();

    pt->Reset(kFALSE);
    pt->SetRemoval(0);
  }
  PropagateForward2(fSeeds);
  fIteration = 0;

  if (fDebug>=0){
    Info("Clusters2Tracks","Time for tracking: \t");timer.Print();timer.Start();
  }
  //activate again some tracks
  for (Int_t i=0; i<fSeeds->GetEntriesFast(); i++) {
    IlcDCHseed *pt=(IlcDCHseed*)fSeeds->UncheckedAt(i), &t=*pt;    
    if (!pt) continue;    
    Int_t nc=t.GetNumberOfClusters();
    if (nc<fRecPar->GetMinNOfClusetrs()) {
//      cout<<" rimossa per nc < 20 0\n";
      delete fSeeds->RemoveAt(i);
      continue;
    } 
    if((t.GetSigmaY2()+t.GetSigmaZ2())>1e10){
      delete fSeeds->RemoveAt(i);
      continue;
    }
    CookLabel(pt,0.1); 
    if (pt->GetRemoval()==10) {
      if (pt->GetDensityFirst(fRecPar->GetMinNOfClusetrs())>fRecPar->GetMinDensity() || pt->GetDensityFirst(int(fRecPar->GetMinNOfClusetrs()*1.5))>fRecPar->GetMinDensity() || pt->GetDensityFirst(fRecPar->GetMinNOfClusetrs()*2)>fRecPar->GetMinDensity())
	pt->Desactivate(10);  // make track again active
      else{
	pt->Desactivate(20); 	
	delete fSeeds->RemoveAt(i);
      }
    } 
  }
  //
  RemoveUsed2(fSeeds,0.85,0.85,0);
  //if (IlcDCHReconstructor::GetDoKinks()) FindKinks(fSeeds,fEvent);
  RemoveUsed2(fSeeds,0.3,0.2,fRecPar->GetMinNOfClusetrs());
 //  //
//   // refit short tracks
//   //
  Int_t nseed=fSeeds->GetEntriesFast();

  //
  Int_t found = 0;
  for (Int_t i=0; i<nseed; i++) {
    IlcDCHseed *pt=(IlcDCHseed*)fSeeds->UncheckedAt(i), &t=*pt;    
    if (!pt) continue;    
    Int_t nc=t.GetNumberOfClusters();
    if (fDebug>0) cout<<"nc= "<<t.GetNumberOfClusters()<<" "<<pt->IsActive()<<" "<<pt->GetRemoval()<<" ";
    if (nc<fRecPar->GetMinNOfClusetrs()) {
//    cout<<" rimossa per nc < 15\n";
      delete fSeeds->RemoveAt(i);
      continue;
    }
    CookLabel(pt,0.1); //For comparison only
    //if ((pt->IsActive() || (pt->fRemoval==10) )&& nc>50 &&pt->GetNumberOfClusters()>0.4*pt->fNFoundable){
    if ((pt->IsActive() || (pt->GetRemoval()==10) )){
      found++;      
      if (fDebug>0) cerr<<found<<'\r';      
      pt->SetLab2(i);
    }
    else
      delete fSeeds->RemoveAt(i);
    
    if (fDebug>0) cout<<endl;      
  }

  if (fDebug>0) 
  //RemoveOverlap(fSeeds,0.99,7,kTRUE);  
  SignShared(fSeeds);  
  //RemoveUsed(fSeeds,0.9,0.9,6);
  // 
  nseed=fSeeds->GetEntriesFast();
  found = 0;
  TTreeSRedirector &cstream = *fDebugStreamer;
  for (Int_t i=0; i<nseed; i++) {
    IlcDCHseed *pt=(IlcDCHseed*)fSeeds->UncheckedAt(i), &t=*pt;    
    if (!pt) continue;    
    Int_t nc=t.GetNumberOfClusters();
    if (nc<fRecPar->GetMinNOfClusetrs()) {
      delete fSeeds->RemoveAt(i);
      continue;
    }
    t.SetUniqueID(i);
    t.CookdEdx(0.02,0.6);
    //    CheckKinkPoint(&t,0.05);
    //if ((pt->IsActive() || (pt->fRemoval==10) )&& nc>50 &&pt->GetNumberOfClusters()>0.4*pt->fNFoundable){
    if ((pt->IsActive() || (pt->GetRemoval()==10) )){
      found++;
      pt->SetLab2(i);
      if (fDebug>0){
	cerr<<found<<'\r';      
        pt->Print("params");
        cout<<"chi2 "<<pt->GetChi2()<<" ncls "<<pt->GetNumberOfClusters()<<" P = "<<pt->GetP()<<" tgl "<<pt->GetTgl()<<" z "<<pt->GetZ()<<" dir "<<pt->GetDir()<<endl;
      }
      if (IlcDCHReconstructor::StreamLevel()>0){
	cstream<<"Found"
	       <<"seed.="<<pt<<"\n"; 
      }
    }
    else
      delete fSeeds->RemoveAt(i);
    
  }

  SortTracks(fSeeds, 1);
  
 
  //  fNTracks = found;
  if (fDebug>=0){
    Info("Clusters2Tracks","Time for overlap removal, track writing and dedx cooking: \t"); timer.Print();timer.Start();
  }
  //
  //  cerr<<"Number of found tracks : "<<"\t"<<found<<endl;  
  Info("Clusters2Tracks","Number of found tracks %d",found);  
  savedir->cd();
  //  UnloadClusters();
  //  
  return 0;
}

void IlcDCHtracker::Tracking(TObjArray * arr)
{
  //
  // tracking of the seeds
  //
  int nlayers=fSector.GetNLayers();
  ParallelTracking(arr,nlayers-1,int(0.5*nlayers));
  ParallelTracking(arr,int(0.5*nlayers),0);
}

TObjArray * IlcDCHtracker::Tracking(Int_t seedtype, Int_t i1, Int_t i2, Float_t cuts[4])
{
  //
  //
  //tracking routine
  TObjArray * arr = new TObjArray;
  // 
  TStopwatch timer;
  timer.Start();

  if (seedtype==3) MakeSeeds3(arr,i1,i2,cuts);
  if (seedtype==4) MakeSeeds5(arr,i1,i2,cuts);    
  if (seedtype==2) MakeSeeds2(arr,i1,i2,cuts);
  
  if (fDebug>0){
    Info("Tracking"," Seeding - seedtype-%d\t from layer %d\t to %d\t nseeds=%d\n",
	 seedtype,i1,i2,arr->GetEntriesFast());
    timer.Print();
    timer.Start();
  }
  Tracking(arr);  
  if (fDebug>0){
    timer.Print();
  }

  return arr;
}

TObjArray * IlcDCHtracker::Tracking()
{
  //
  //
  TStopwatch timer;
  timer.Start();
  Int_t nup=fSector.GetNLayers();

  TObjArray * seeds = new TObjArray;
  TObjArray * arr=0;
  
  Int_t gap =fRecPar->GetMinNOfClusetrs();
  Float_t cuts[4];
  cuts[0] = 0.05;
  cuts[1] = 1.5;
  cuts[2] = 100;
  cuts[3] = 5;
  Float_t fnumber  = 3.0;
  Float_t fdensity = 3.0;
  
  //  
  //find primaries  
  //
  //seed 3 
  //  cuts[0]   - fP4 cut
  //  cuts[2]   - zvertex cut
  bool oposstereo=false;

  for (Int_t delta = 0; delta<int(nup*0.3); delta+=3){
    //
    cuts[0]=0.060;
    arr = Tracking(3,nup-1-delta,nup-1-delta-gap,cuts);
    SumTracks(seeds,arr);   
    //SignClusters(seeds,fnumber,fdensity); 
    if(oposstereo){
      arr = Tracking(3,nup-1-delta,nup-1-delta-gap-1,cuts);
      SumTracks(seeds,arr);
    }   
    //  SignClusters(seeds,fnumber,fdensity); 
    //
    for (Int_t i=2;i<3;i+=2){
      // seed high pt tracks
      cuts[0]=0.050;
      arr = Tracking(3,nup-i-delta,nup-i-delta-gap,cuts);
      SumTracks(seeds,arr);   
      //SignClusters(seeds,fnumber,fdensity);        

      if(oposstereo){
	arr = Tracking(3,nup-i-delta,nup-i-delta-gap-1,cuts);
	SumTracks(seeds,arr);
      }   
      //gnClusters(seeds,fnumber,fdensity);        
    }
    SignClusters(seeds,fnumber,fdensity);        
  }
  fnumber  = 4;
  fdensity = 4.;

  //find primaries  
  double rmin=fParam->GetInnerRadius();

  for (Int_t delta = int(nup*0.3); delta<=nup-gap-15; delta+=10){
    //
    // seed high pt tracks
    cuts[0]=1./(0.5*rmin);

    arr = Tracking(3,nup-delta,nup-delta-gap,cuts);
    SumTracks(seeds,arr);   
    // SignClusters(seeds,fnumber,fdensity);            

    if(oposstereo){
      arr = Tracking(3,nup-delta,nup-delta-gap-1,cuts);
      SumTracks(seeds,arr);
    }   
    arr = Tracking(3,nup-delta-1,nup-delta-gap-1,cuts);
    SumTracks(seeds,arr);   
    // SignClusters(seeds,fnumber,fdensity);            
    
    if(oposstereo){
      arr = Tracking(3,nup-delta-1,nup-delta-gap-2,cuts);
      SumTracks(seeds,arr);
    }   
    SignClusters(seeds,fnumber,fdensity);            

    cuts[0]=1./(0.5*rmin);
    
    if(oposstereo){
      arr = Tracking(3,nup-delta-5+1,nup-delta-5-gap,cuts);
      SumTracks(seeds,arr);
    }   
    // SignClusters(seeds,fnumber,fdensity);            
    
    arr = Tracking(3,nup-delta-5,nup-delta-5-gap,cuts);
    SumTracks(seeds,arr);   
    arr = Tracking(3,nup-delta-5+1,nup-delta-5-gap+1,cuts);
    SumTracks(seeds,arr);   
    // SignClusters(seeds,fnumber,fdensity);            
    if(oposstereo){
      arr = Tracking(3,nup-delta-5,nup-delta-5-gap-1,cuts);
      SumTracks(seeds,arr);
    }   
    SignClusters(seeds,fnumber,fdensity);            
  }

  for (Int_t delta = 10 ; delta>=0; delta--){
    cuts[0]=1./(0.5*rmin);

    if(oposstereo){
      arr = Tracking(3,delta+gap+1,delta,cuts);
      SumTracks(seeds,arr);   
      SignClusters(seeds,fnumber,fdensity);            
    }

    arr = Tracking(3,delta+gap,delta,cuts);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity);            
    
    arr = Tracking(3,delta+gap-2,delta,cuts);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity);            
  }

  if (fDebug>0){
    Info("Tracking()","Primary seeding\t%d ",seeds->GetEntriesFast());
    timer.Print();
    timer.Start();
  }
  
  // find secondaries

  cuts[0] = 0.05;
  cuts[1] = 1.5;
  cuts[2] = 100.;
  cuts[3] = 3.5;
  fnumber  = 2.;
  fdensity = 2.;

  arr = Tracking(4,nup-1,nup-1-gap,cuts);
  SumTracks(seeds,arr);   
  SignClusters(seeds,fnumber,fdensity);   
  //
  arr = Tracking(4,nup-2,nup-2-gap,cuts);
  SumTracks(seeds,arr);   
  SignClusters(seeds,fnumber,fdensity);   
  //
  arr = Tracking(4,nup-3,nup-3-gap,cuts);
  SumTracks(seeds,arr);   
  SignClusters(seeds,fnumber,fdensity);   
  //


  for (Int_t delta = 3; delta<int(nup*0.3); delta+=5){
    //
    cuts[0] = 0.05;
    cuts[1] = 1.5;
    cuts[2] = 100.;
    cuts[3] = 3.5;
    arr = Tracking(4,nup-1-delta,nup-1-delta-gap,cuts);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity);   
  } 
  fnumber  = 1;
  fdensity = 1;
  //
  // change cuts
  fnumber  = 2.;
  fdensity = 2.;
  cuts[0]=0.05;

  // find secondaries
  for (Int_t delta = int(nup*0.3); delta<=nup-1-gap; delta+=10){
    //
    cuts[0] = 0.05;
    cuts[1] = 3.5;
    cuts[2] = 100.;
    cuts[3] = 3.5;
    arr = Tracking(4,nup-1-delta,nup-1-delta-gap,cuts);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity);   
  }
 
  if (fDebug>0){
    Info("Tracking()","Secondary seeding\t%d ",seeds->GetEntriesFast());
    timer.Print();
    timer.Start();
  }

/*
  for (Int_t delta = 0; delta<int(nup*0.3)-6; delta+=6){
    //
    cuts[0]=0.060;
    cuts[1] = 1.5;
    arr = Tracking(3,int(nup*0.3)+gap-delta,int(nup*0.3)-delta,cuts);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity); 
    //
    for (Int_t i=2;i<6;i+=2){
      // seed high pt tracks
      cuts[0]=0.030;
      cuts[1]=0.3;
      arr = Tracking(3,int(nup*0.3)+gap-delta-i,int(nup*0.3)-delta-i,cuts);
      SumTracks(seeds,arr);   
      SignClusters(seeds,fnumber,fdensity);        
    }
  }
//   if (fDebug>0){
    Info("Tracking()","Third seeding\t%d ",seeds->GetEntriesFast());
    timer.Print();
    timer.Start();
//   }
*/
  return seeds;
  //
      
}


TObjArray * IlcDCHtracker::TrackingSpecial()
{
  //
  // seeding adjusted for laser and cosmic tests - short tracks with big inclination angle
  // no primary vertex seeding tried
  //
  TStopwatch timer;
  timer.Start();
  Int_t nup=fSector.GetNLayers();

  TObjArray * seeds = new TObjArray;
  TObjArray * arr=0;
  
  Int_t   gap  = 15;
  Float_t cuts[4];
  Float_t fnumber  = 3.0;
  Float_t fdensity = 3.0;
  
  // find secondaries
  cuts[0] = IlcDCHReconstructor::GetMaxC();   // max curvature
  cuts[1] = 3.5;    // max tan(phi) angle for seeding
  cuts[2] = 3.;     // not used (cut on z primary vertex)     
  cuts[3] = 3.5;    // max tan(theta) angle for seeding

  for (Int_t delta = 0; nup-delta-gap-1>0; delta+=3){
    //
    arr = Tracking(4,nup-1-delta,nup-1-delta-gap,cuts);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity);   
  } 
 
  if (fDebug>0){
    Info("Tracking()","\n\nSecondary seeding\t%d\n\n",seeds->GetEntriesFast());
    timer.Print();
    timer.Start();
  }

  return seeds;
  //
      
}


void IlcDCHtracker::SumTracks(TObjArray *arr1,TObjArray *arr2) const
{
  //
  //sum tracks to common container
  //remove suspicious tracks
  Int_t nseed = arr2->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    IlcDCHseed *pt=(IlcDCHseed*)arr2->UncheckedAt(i);    
    if (pt){
      //
      // remove tracks with too big curvature
      //
      if (TMath::Abs(pt->GetC())>IlcDCHReconstructor::GetMaxC()){
	delete arr2->RemoveAt(i);
	continue;
      }
       // REMOVE VERY SHORT  TRACKS
      if(pt->GetLastPoint()>2*fRecPar->GetMinNOfClusetrs()){
	if (pt->GetNumberOfClusters()<fRecPar->GetMinNOfClusetrs()){
	  delete arr2->RemoveAt(i);
	  continue;
	}
      }else{
	if (pt->GetNumberOfClusters()<fRecPar->GetMinNOfClusetrs()-1){
	  delete arr2->RemoveAt(i);
	  continue;
	}	
      }
      // NORMAL ACTIVE TRACK
      if (pt->IsActive()){
	arr1->AddLast(arr2->RemoveAt(i));
	continue;
      }
      //remove not usable tracks
      if (pt->GetRemoval()!=10){
	delete arr2->RemoveAt(i);
	continue;
      }
     
      // ENABLE ONLY ENOUGH GOOD STOPPED TRACKS
      if (pt->GetDensityFirst(fRecPar->GetMinNOfClusetrs())>fRecPar->GetMinDensity() || pt->GetDensityFirst(int(1.5*fRecPar->GetMinNOfClusetrs()))>fRecPar->GetMinDensity() || pt->GetDensityFirst(2*fRecPar->GetMinNOfClusetrs())>fRecPar->GetMinDensity())
	arr1->AddLast(arr2->RemoveAt(i));
      else{      
	delete arr2->RemoveAt(i);
      }
    }
  }
  delete arr2;  
}



void  IlcDCHtracker::ParallelTracking(TObjArray * arr, Int_t rfirst, Int_t rlast)
{
  //
  // try to track in parralel
  Int_t nlayers=fSector.GetNLayers();

  Int_t nseed=arr->GetEntriesFast();
  //prepare seeds for tracking
  for (Int_t i=0; i<nseed; i++) {
    IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i), &t=*pt; 
    if (!pt) continue;
    if (!t.IsActive()) continue;
    // follow prolongation to the first layer
    FollowProlongation(t, rfirst+1);
  }


  //
  for (Int_t nr=rfirst; nr>=rlast; nr--){ 
    if(nseed>0) 
      if(IlcDebugLevelClass()>=7)
	IlcDebug(5,Form("Doit layer %i",nr));
    // make indexes with the cluster tracks for given       

    // find nearest cluster
    for (Int_t i=0; i<nseed; i++) {
      IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i), &t=*pt;      
      if (!pt) continue;
      if(IlcDebugLevelClass()>=7)
	IlcDebug(5,Form("Nseed=%i active=%i, nclusters=%i",i,int(pt->IsActive()),pt->GetNumberOfClusters()));

      if (nr==int(0.5*nlayers)) pt->UpdateReference();
      if (!pt->IsActive()) continue;
      UpdateClusters(t,nr);
    }
    // prolonagate to the nearest cluster - if founded
    for (Int_t i=0; i<nseed; i++) {
      IlcDCHseed *pt=(IlcDCHseed*)arr->UncheckedAt(i); 
      if (!pt) continue;
      if(IlcDebugLevelClass()>=7)
	IlcDebug(5,Form("Nseed=%i active=%i, nclusters=%i",i,int(pt->IsActive()),pt->GetNumberOfClusters()));
      if (!pt->IsActive()) continue; 
      FollowToNextCluster(*pt,nr);
    }
  }    
}

void IlcDCHtracker::PrepareForBackProlongation(TObjArray * arr,Float_t fac) const
{
  //
  //
  // if we use DCH track itself we have to "update" covariance
  //
  Int_t nseed= arr->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    IlcDCHseed *pt = (IlcDCHseed*)arr->UncheckedAt(i);
    if (pt) {
      pt->Modify(fac);
	
    }
    
  }


}
void IlcDCHtracker::PrepareForProlongation(TObjArray * arr, Float_t fac) const
{
  //
  //
  // if we use DCH track itself we have to "update" covariance
  //
  Int_t nseed= arr->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    IlcDCHseed *pt = (IlcDCHseed*)arr->UncheckedAt(i);
    if (pt) {
      pt->Modify(fac);
      pt->SetFirstPoint(pt->GetLastPoint()); 
    }
    
  }


}

Int_t IlcDCHtracker::PropagateBack(TObjArray * arr)
{
  //
  // make back propagation
  //
  Int_t nseed= arr->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    if(fDebug>15) IlcDebug(16,Form("back seed=%i from =%i",i,nseed));
    IlcDCHseed *pt = (IlcDCHseed*)arr->UncheckedAt(i);
    if (pt&& pt->GetKinkIndex(0)<=0) { 
      FollowBackProlongation(*pt,fSector.GetNLayers()-1);     
    }
    if (pt&& pt->GetKinkIndex(0)>0) {
      //      IlcESDkink * kink = fEvent->GetKink(pt->GetKinkIndex(0)-1);
      //pt->SetFirstPoint(kink->GetTPCRow0());
      FollowBackProlongation(*pt,fSector.GetNLayers()-1);  
    }
    
  }
  return 0;
}


Int_t IlcDCHtracker::PropagateForward2(TObjArray * arr)
{
  //
  // make forward propagation
  //
  Int_t nseed= arr->GetEntriesFast();
  //
  for (Int_t i=0;i<nseed;i++){
    if(fDebug>15) IlcDebug(16,Form("seed=%i from =%i",i,nseed));
    IlcDCHseed *pt = (IlcDCHseed*)arr->UncheckedAt(i);
    if (pt) { 
      FollowProlongation(*pt,0);
    }
  }
  return 0;
}


Int_t IlcDCHtracker::PropagateForward()
{
  //
  // propagate track forward
  cout<<"Wrong not used"<<endl;
  return 1;
  ParallelTracking(fSeeds,fSector.GetNLayers()-1,0);
  return 1;
}


//__________________________________________________________________________
Int_t IlcDCHtracker::CookLabel(IlcKalmanTrack *tr, Float_t wrong,Int_t first, Int_t last) const {
  //--------------------------------------------------------------------
  //This function "cooks" a track label. If label<0, this track is fake.
  //--------------------------------------------------------------------
  IlcDCHseed *t=dynamic_cast<IlcDCHseed *>(tr);
  std::map<int,int> labels;
  int numberofpoints=t->GetNumberOfClusters(),noc=0;
  for(int i=0;i<kMaxLayer;i++) for(int j=0;j<kMaxInLayer;j++) {
    if (i<first||i>last) continue;

    Int_t cindex = t->GetClusterIndex(i,j);
    if (cindex<0||(cindex&ClusterIndex::NotAccept)) continue;
    IlcCluster *cl = t->GetClusterPointer(i,j);
    if(!cl) continue;
    labels[cl->GetLabel(0)]++;
    noc++;
  }
  int maxlabel=-1,ntimes=0;
  for(std::map<int,int>::iterator it=labels.begin();it!=labels.end();it++){
    if((*it).second>ntimes){
      maxlabel=(*it).first;
      ntimes=(*it).second;
    }
  }

  if((noc!=numberofpoints)&&first==0&&last==kMaxLayer-1){
    IlcWarning(Form("Something wrong: number of found clusters=%i not equal to track=%i",noc,numberofpoints));
    int nn=0;
    for(int i=kMaxLayer-1;i>=0;i--) for(int j=0;j<kMaxInLayer;j++) {
      Int_t cindex = t->GetClusterIndex(i,j);
      if (cindex<0||(cindex&ClusterIndex::NotAccept)) continue;
      cout<<nn<<"lay "<<i<<" "<<j<<" "<<cindex<<" "<<t->GetClusterPointer(i,j)<<endl;
      nn++;
    }
    cout<<"TrackId "<<maxlabel<<" "<<ntimes<<endl;
  }
  // In fast recpoints first label correspond to real trackid(but with bigger sigma)
  /*  Int_t lflag=0;
  
  for(Int_t i=0;i<numberofpoints;i++){
    Int_t cindex = t->GetClusterIndex(i);
    IlcCluster *cl = GetCluster(cindex);
    
    if(cl->GetLabel(0)==maxlabel || 
       cl->GetLabel(1)==maxlabel || 
       cl->GetLabel(2)==maxlabel) lflag++;
  }
  */


  double fakeratio=1.-double(ntimes)/(noc>0?noc:1);
  if(maxlabel<0) maxlabel=1000000-maxlabel;

  t->SetFakeRatio(fakeratio);
  t->SetLabel(t->GetFakeRatio()>wrong?-maxlabel:maxlabel);

  return t->GetLabel();
}

Int_t  IlcDCHtracker::IlcDCHSector::GetLayerNumber(Double_t x,double z) const 
{
  //return pad layer number for this x
  Double_t rup,rlow,rmid;
  int nup=fN-1,nlow=0,nmid;
  
  if(fN<1) return 0;
  rup=fLayer[fN-1].GetRAtZ(z);
  if (x > rup) return fN-1;
  rlow=fLayer[0].GetRAtZ(z);
  if (x < rlow) return 0;
  while(nup>nlow+1){
    nmid=int(0.5*(nup+nlow));
    rmid=fLayer[nmid].GetRAtZ(z);//0.5*(fLayer[nmid].GetRAtZ(z)+fLayer[nmid+1].GetRAtZ(z));
    if(rmid<x) {nlow=nmid;rlow=rmid;}
    else {nup=nmid;rup=rmid;}
  };
  return fabs(x-rlow)<fabs(x-rup)?nlow:nup;
}

//_________________________________________________________________________
void IlcDCHtracker::IlcDCHSector::Setup(const IlcDCHParam *par,IlcDCHwireposition *wirepos) {
  //-----------------------------------------------------------------------
  // Setup inner sector
  //-----------------------------------------------------------------------
  fN=par->GetNLayers();
  fLayer=new IlcDCHLayer[fN];
  int nsuper=par->GetSuperLayerNum();
  int nrings=par->GetRingNum()-1;
  std::cout<<"Dimension below in cm!"<<std::endl;
  for(int i=0;i<nsuper;i++){
    for(int j=0;j<nrings;j++){
      wirepos->SelectWire(i,j,0);
      Float_t right[3],left[3];
      wirepos->WirePosAtEndcap(right,left);
      int layer=i*nrings+j;
      IlcDCHLayer& klay= fLayer[layer];
      klay.SetR(TMath::Hypot(left[0]+right[0],left[1]+right[1])/2);
      klay.SetStereoAngle(TMath::Sign(wirepos->GetWireEpsilon(),
				      (right[1]*left[0]-right[0]*left[1])));
      klay.SetUseZ((fabs(par->GetDrop())<1e-4));
      klay.SetMaxZ(fabs(left[2]));
      klay.SetFirstPhi(TMath::ATan2(left[1]+right[1],left[0]+right[0]));
      klay.SetNWires(par->GetSWireNum()+par->GetSDeltaWireNum()*i);

      if(i>=10 && i<18) klay.SetNWires(par->GetSWireNum()+par->GetSDeltaWireNum()*(i-4));
      if(i>=18 && i<24) klay.SetNWires(par->GetSWireNum()+par->GetSDeltaWireNum()*(i-6));

      klay.SetDeltaPhi(2*TMath::Pi()/klay.GetNWires());
      //cout<<layer<<" slayer "<<i+1<<" ring "<<j<<" Nwire "<<klay.GetNWires()<<" r "<<klay.GetR()<<" sa "<<klay.GetStereoAngle()<<" z "<<klay.GetMaxZ()<<endl;

      double prevphi=klay.GetFirstPhi(),newphi;
      for(int w=1;w<klay.GetNWires();w++){
	wirepos->SelectWire(i,j,w);
	wirepos->WirePosAtEndcap(right,left);
	if(fabs(TMath::Hypot(left[0]+right[0],left[1]+right[1])/2-klay.GetR())>1e-4)
	  cout<<"Different radius at same layer r="<<
	    TMath::Hypot(left[0]+right[0],left[1]+right[1])<<endl;
	if(fabs(TMath::Tan(TMath::Sign(wirepos->GetWireEpsilon(),
				       (right[1]*left[0]-right[0]*left[1])))
		-klay.GetTStereoAngle())>1e-4)
	  cout<<"Different stereoangle at same layer r="<<
	    TMath::Sign(wirepos->GetWireEpsilon(),
			(right[0]*left[1]-right[1]*left[0]))<<endl;
	newphi=TMath::ATan2(left[1]+right[1],left[0]+right[0]);
	double dphi=newphi-prevphi;
	if(dphi<-TMath::Pi()) dphi+=TMath::Pi()*2;
	if(dphi>TMath::Pi())  dphi-=TMath::Pi()*2;
	if(dphi<0&&w==1)klay.SetDeltaPhi(-klay.GetDeltaPhi());
	if(fabs(dphi-klay.GetDeltaPhi())>1e-6)
	  cout<<"Different dphi on same layer="<<
	    dphi<<" must be "<<klay.GetDeltaPhi()<<endl;
	prevphi=newphi;
	  
      }
      double nn[3]={-TMath::Sin(klay.GetStereoAngle())*TMath::Sin(klay.GetFirstPhi()),
		    TMath::Sin(klay.GetStereoAngle())*TMath::Cos(klay.GetFirstPhi()),
		    TMath::Cos(klay.GetStereoAngle())};
      cout<<layer<<" r "<<klay.GetR()<<" phi0 "<<klay.GetFirstPhi()<<" dphi "<<klay.GetDeltaPhi()<<" sangle "<<klay.GetStereoAngle()<<" z "<<klay.GetMaxZ()<<" nw "<<klay.GetNWires()<<" n=("<<nn[0]<<","<<nn[1]<<","<<nn[2]<<")"<<endl;

      
    }
  }
}

IlcDCHtracker::IlcDCHLayer::IlcDCHLayer() {
  //
  // default constructor
  fN=0;
  fNarray=0;
  fClustersArray=0;
  fTStereoAngle=0;
  fIsStereo=false;
  fIsUseZ=true;
}

IlcDCHtracker::IlcDCHLayer::~IlcDCHLayer(){
  //
  ResetClusters();
}



//_________________________________________________________________________
void 
IlcDCHtracker::IlcDCHLayer::InsertCluster(const IlcDCHcluster* c, UInt_t index) {
  //-----------------------------------------------------------------------
  // Insert a cluster into this pad layer in accordence with its y-coordinate
  //-----------------------------------------------------------------------
  if(fN==fNarray){
    cerr<<"IlcDCHLayer::InsertCluster(): Insert more than was declared "<<
      fN<<" from "<<fNarray<<"!\n"; return;
  }
  if (fN==kMaxClusterPerLayer) {
    cerr<<"IlcDCHLayer::InsertCluster(): Too many clusters !\n"; return;
  }
  if (fN==0) {fIndex[0]=index; fClusters[fN++]=c; return;}
  Int_t i=Find(c->GetZ());
  memmove(fClusters+i+1 ,fClusters+i,(fN-i)*sizeof(IlcDCHcluster*));
  memmove(fIndex   +i+1 ,fIndex   +i,(fN-i)*sizeof(UInt_t));
  fIndex[i]=index; fClusters[i]=c; fN++;
}

void IlcDCHtracker::IlcDCHLayer::ResetClusters() {
   //
   // reset clusters
   fN  = 0; 
   fNarray=0;
   if(fClustersArray) delete[] fClustersArray;
   fClustersArray=0;
}


//___________________________________________________________________
Int_t IlcDCHtracker::IlcDCHLayer::Find(Double_t z) const {
  //-----------------------------------------------------------------------
  // Return the index of the nearest cluster 
  //-----------------------------------------------------------------------
  if (fN==0) return 0;
  if (z <= fClusters[0]->GetZ()) return 0;
  if (z > fClusters[fN-1]->GetZ()) return fN;
  Int_t b=0, e=fN-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (z > fClusters[m]->GetZ()) b=m+1;
    else e=m; 
  }
  return m;
}



//___________________________________________________________________
IlcDCHcluster * IlcDCHtracker::IlcDCHLayer::FindNearest(Double_t phi, Double_t z, 
							Double_t roadphiR, Double_t roadz) const {
  //-----------------------------------------------------------------------
  // Return the index of the nearest cluster in z y 
  //-----------------------------------------------------------------------
  Float_t maxdistance = 2;
  if(!fIsUseZ) roadz=1e6;

  IlcDCHcluster *cl =0;
  for (Int_t i=Find(z-roadz); i<fN; i++) {
      IlcDCHcluster *c=(IlcDCHcluster*)(fClusters[i]);
      if (c->GetZ() > z+roadz) break;		
      double wphi=fabs(Phi_mpi_pi(GetWirePhi(c->GetW(),z)-phi))*fR-fabs(c->GetImP());
      if ( fabs(wphi)*fR >  roadphiR ) continue;
      Float_t distance = (c->GetZ()-z)*(c->GetZ()-z)/(roadz*roadz)+wphi*wphi/(roadphiR*roadphiR);
      if (maxdistance>distance) {
	maxdistance = distance;
	cl=c;       
      }
  }
  return cl;      
}

IlcDCHcluster * IlcDCHtracker::IlcDCHLayer::FindNearest2(Double_t phi,Double_t snp, Double_t z, Double_t roadphiR, Double_t roadz,UInt_t & index, int first) const 
{
  //-----------------------------------------------------------------------
  // Return the index of the nearest cluster in z y 
  //-----------------------------------------------------------------------
  Float_t maxdistance[10] = {2,2,2,2,2,2,2,2,2,2};
  IlcDCHcluster *cl[10] = {0,0,0,0,0,0,0,0,0,0};
  UInt_t ind[10]={0,0,0,0,0,0,0,0,0,0};

  if(!fIsUseZ) roadz=1e6;
  //PH Check boundaries. 510 is the size of fFastCluster
  Int_t iz1 = Int_t(((z-roadz)*fCoeffConv+1.)*255-0.5);
  if (iz1<0) iz1=0;
  if (iz1>=NFAST) iz1=NFAST-1;
  iz1 = TMath::Max(GetFastClusterZ(iz1)-1,0);
  Int_t iz2 = Int_t(((z+roadz)*fCoeffConv+1.)*255+0.5);
  if (iz2<0) iz2=0;
  if (iz2>=NFAST) iz2=NFAST-1;
  iz2 = TMath::Min(GetFastClusterZ(iz2)+1,fN);
  
  for (Int_t i=iz1; i<iz2; i++) {
      IlcDCHcluster *c=(IlcDCHcluster*)(fClusters[i]);
      if (c->GetZ() > z+roadz) break;
      double wphir=fabs(Phi_mpi_pi(GetWirePhi(c->GetW(),z)-phi)*fR*sqrt(1-snp*snp))-fabs(c->GetImP());
      if ( fabs(wphir) >  roadphiR ) continue;
      Float_t distance = (c->GetZ()-z)*(c->GetZ()-z)/(roadz*roadz)+wphir*wphir/(roadphiR*roadphiR);
      if (maxdistance[first]>distance) {
	for(int j=0;j<first+1;j++){
	  if(maxdistance[j]>distance){
	    for(int j2=first;j2>j;j2--){
	      maxdistance[j]=maxdistance[j-1];
	      cl[j]=cl[j-1];
	      ind[j]=ind[j-1];
	    }
	    maxdistance[j] = distance;
	    cl[j]=c;
	    ind[j]=i;
	    break;
	  }
	}
      }
  }
  index=ind[first];
  return cl[first];      
}

void IlcDCHtracker::IlcDCHLayer::FindNearest2(Double_t phi,Double_t snp,Double_t z, 
					      Double_t roadphiR, Double_t roadz,int first,
					      IlcDCHcluster* clar[],UInt_t inar[]) const{
  //-----------------------------------------------------------------------
  // Return the index of the nearest cluster in z y 
  //-----------------------------------------------------------------------
  Float_t maxdistance[10] = {2,2,2,2,2,2,2,2,2,2};
  IlcDCHcluster *cl[10] = {0,0,0,0,0,0,0,0,0,0};
  UInt_t ind[10]={0,0,0,0,0,0,0,0,0,0};
  
  if(!fIsUseZ) roadz=1e6;
  //PH Check boundaries. 510 is the size of fFastCluster
  Int_t iz1 = Int_t(((z-roadz)*fCoeffConv+1.)*255-0.5);
  if (iz1<0) iz1=0;
  if (iz1>=NFAST) iz1=NFAST-1;
  iz1 = TMath::Max(GetFastClusterZ(iz1)-1,0);
  Int_t iz2 = Int_t(((z+roadz)*fCoeffConv+1.)*255+0.5);
  if (iz2<0) iz2=0;
  if (iz2>=NFAST) iz2=NFAST-1;
  iz2 = TMath::Min(GetFastClusterZ(iz2)+1,fN);
  
  for (Int_t i=iz1; i<iz2; i++) {
      IlcDCHcluster *c=(IlcDCHcluster*)(fClusters[i]);
      if (c->GetZ() > z+roadz) break;
      double wphir=fabs(Phi_mpi_pi(GetWirePhi(c->GetW(),z)-phi)*fR*sqrt(1-snp*snp))-fabs(c->GetImP());
      Float_t distance = (c->GetZ()-z)*(c->GetZ()-z)/(roadz*roadz)+wphir*wphir/(roadphiR*roadphiR);
      // std::cout<<i<<" wphi "<<wphir<<" "<<distance<<" phi "<<Phi_0_2pi(phi)<<" "<<GetWirePhi(c->GetW(),z)<<" snp "<<snp
      // 	       <<" dist "<<fabs(Phi_mpi_pi(GetWirePhi(c->GetW(),z)-phi)*fR*sqrt(1-snp*snp))<<" "<<c->GetImP()<<endl;
      if ( fabs(wphir) >  roadphiR ) continue;
      if (maxdistance[first]>distance) {
	for(int j=0;j<first+1;j++){
	  if(maxdistance[j]>distance){
	    for(int j2=first;j2>j;j2--){
	      maxdistance[j2]=maxdistance[j2-1];
	      cl[j2]=cl[j2-1];
	      ind[j2]=ind[j2-1];
	    }
	    maxdistance[j] = distance;
	    cl[j]=c;
	    ind[j]=i;
	    break;
	  }
	}
      }
  }
  for(int i=0;i<first+1;i++){
    clar[i]=cl[i];
    inar[i]=ind[i];
  }
}

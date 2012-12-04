#ifndef ILCDCHTRACK_H
#define ILCDCHTRACK_H
/* Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 * See cxx source for full Copyright notice                               */
/* $Id: IlcDCHtrack.h,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */

//-------------------------------------------------------
//                    DCH Track Class
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//
// The track parameterization is fixed in the following way:                   
//      param0:   local Y-coordinate of a track (cm)                
//      param1:   local Z-coordinate of a track (cm)                
//      param2:   local sine of the track momentum azimuth angle    
//      param3:   tangent of the track momentum dip angle           
//      param4:   1/pt (1/(GeV/c))                                  
//
//-------------------------------------------------------

#include "IlcKalmanTrack.h"
#include <TMath.h>

#include "IlcDCHreco.h"
#include "IlcExternalTrackParam.h"
//class IlcESDtrack;
//class IlcESDVertex;
class IlcDCHcluster;
#include <iostream>
using namespace std; 

//_____________________________________________________________________________
class IlcDCHtrack : public IlcKalmanTrack {
public:
  IlcDCHtrack();
  IlcDCHtrack(Double_t x, Double_t alpha, const Double_t p[5], 
              const Double_t cov[15], Int_t index); 
  // IlcDCHtrack(const IlcESDtrack& t);
  IlcDCHtrack(const IlcDCHtrack& t);
  virtual ~IlcDCHtrack() {}

  Int_t Compare(const TObject *o) const;

  void SetdEdx(Double_t dedx) {fdEdx=dedx;}
  Double_t GetdEdx()  const {return fdEdx;}
  Double_t GetPIDsignal()  const {return GetdEdx();}

  Int_t GetClusterIndex(Int_t i,int j) const {return fIndex[i][j];}
  void  SetClusterIndex(Int_t i,Int_t idx,int j) {
    fIndex[i][j]=idx;
    if(idx>=0&&((idx&ClusterIndex::Layer)>>ClusterIndex::LayerOffset)!=UInt_t(i))
      cout<<"Set wrong cluster index "<<((idx&ClusterIndex::Layer)>>ClusterIndex::LayerOffset)<<" lay "<<i<<endl;
  }
  Int_t GetMaxNumberOfLayers() const { return kMaxLayer;}

  Double_t GetC() const {return IlcExternalTrackParam::GetC(GetBz());}
  Double_t GetSigma2C() const {
    Double_t cnv=GetBz()*kB2C;
    return GetSigma1Pt2()*cnv*cnv;
  }

  Double_t GetPredictedChi2(const IlcCluster *cluster) const;
  Double_t GetPredictedChi2(const IlcCluster *cluster2,
			    double rwire,double anglewire,double tanstereo,bool withZ) const;
  Int_t GetProlongation(Double_t xr, Double_t &y, Double_t & z) const;
  //  Bool_t PropagateTo(Double_t xr, Double_t rho=0.635e-3, Double_t x0=24.2096);
  Bool_t PropagateTo(Double_t xr, Double_t rho=0.531e-3, Double_t x0=21.1628);
  Bool_t Update(const IlcCluster *c, Double_t chi2, Int_t i);
  Bool_t Update(const IlcCluster* c2, Double_t chi2, Int_t i,
		double rwire,double anglewire,double tanstereo,bool withZ);

  Bool_t Rotate(Double_t alpha) {
    cout<<"Don't use it Rotate "<<GetAlpha()<<" "<<alpha<<endl;
    return IlcExternalTrackParam::Rotate(GetAlpha()+alpha,GetBz());
  }
  
  // Bool_t PropagateToVertex(const IlcESDVertex *v, 
  //                          Double_t rho=1.2e-3, Double_t x0=36.66);
  Double_t GetYat(Double_t x) const;


  void ResetClusters() {SetNumberOfClusters(0); SetChi2(0.);}
  void UpdatePoints();//update points 
  Float_t* GetPoints() {return fPoints;}

  Float_t Density(Int_t layer0, Int_t layer1); //calculate cluster density
  Float_t Density2(Int_t layer0, Int_t layer1); //calculate cluster density
  Double_t GetD(Double_t x=0, Double_t y=0) const {
     return IlcExternalTrackParam::GetD(x,y,GetBz());
  }
  IlcExternalTrackParam & GetReference(){ return fReference;}
  void UpdateReference(){ new (&fReference) IlcExternalTrackParam(*this);}
  Int_t   GetKinkIndex(Int_t i) const{ return fKinkIndexes[i];}
  Int_t*  GetKinkIndexes() { return fKinkIndexes;}
  Int_t   GetV0Index(Int_t i) const{ return fV0Indexes[i];}
  Int_t*  GetV0Indexes() { return fV0Indexes;}

  void SetFirstPoint(Int_t f) {fFirstPoint=f;}
  void SetLastPoint(Int_t f) {fLastPoint=f;}
  void SetRemoval(Int_t f) {fRemoval=f;}
  void SetLab2(Int_t f) {fLab2=f;}
  void SetKinkIndex(Int_t i, Int_t idx) {fKinkIndexes[i]=idx;}
  void SetBConstrain(Bool_t b) {fBConstrain=b;}
  void SetNFoundable(Int_t n) {fNFoundable=n;}
  void SetNShared(Int_t s) {fNShared=s;}

  Int_t GetFirstPoint() const {return fFirstPoint;}
  Int_t GetLastPoint() const {return fLastPoint;}
  Int_t GetRemoval() const {return fRemoval;}
  Int_t GetLab2() const {return fLab2;}
  Bool_t GetBConstrain() const {return fBConstrain;}
  Int_t GetNShared() const {return fNShared;}
  Int_t GetNFoundable() const {return fNFoundable;}

protected: 
  Double_t GetBz() const;
  Double_t fdEdx;           // dE/dx

  Int_t fIndex[kMaxLayer][kMaxInLayer];       // indices of associated clusters 
  Int_t fNInLayer[kMaxLayer];    //n clusters in layer
  Float_t fPoints[4];            //first, max dens layer  end points of the track and max density
  // MI addition
  Float_t fSdEdx;           // sigma of dedx 
  //
  Int_t   fNFoundable;      //number of foundable clusters - dead zone taken to the account
  Bool_t  fBConstrain;   // indicate seeding with vertex constrain
  Int_t   fLastPoint;     // last  cluster position     
  Int_t   fFirstPoint;    // first cluster position
  Int_t fRemoval;         // removal factor
  Int_t fTrackType;       // track type - 0 - normal - 1 - kink -  2 -V0  3- double found
  Int_t fLab2;            // index of corresponding track (kink, V0, double)
  Int_t   fNShared;       // number of shared points 
  IlcExternalTrackParam   fReference; // track parameters at the middle of the chamber
  Float_t  fKinkPoint[12];      //radius, of kink,  dfi and dtheta
  Int_t    fKinkIndexes[3];     // kink indexes - minus = mother + daughter
  Int_t    fV0Indexes[3];     // kink indexes - minus = mother + daughter

  ClassDef(IlcDCHtrack,3)   // Time Projection Chamber reconstructed tracks
};

#endif


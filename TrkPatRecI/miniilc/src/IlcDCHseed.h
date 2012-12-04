#ifndef ILCDCHSEED_H
#define ILCDCHSEED_H
/* Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 * See cxx source for full Copyright notice                               */


/* $Id: IlcDCHseed.h,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */

//-------------------------------------------------------
//   DCH seed
//   Class  needed for DCH parallel tracking 
//
//   Origin: 
//-------------------------------------------------------

#include <TError.h>

#include "IlcDCHtrack.h"
#include "IlcComplexCluster.h"
#include "IlcPID.h"

class TFile;
class IlcDCHParam;
class IlcDCHseed;
class IlcDCHcluster;
class IlcDCHTrackerPoint;
//class IlcESD;   
class TClonesArray;

class IlcDCHseed : public IlcDCHtrack {
  public:  
     IlcDCHseed();
     virtual ~IlcDCHseed();
     IlcDCHseed(const IlcDCHtrack &t);
     IlcDCHseed(const IlcDCHseed &s, Bool_t clusterOwner = kFALSE);
     IlcDCHseed(Double_t xr, Double_t alpha, const Double_t xx[5], 
                const Double_t cc[15], Int_t i);     
  // DCH specific propagation
     Double_t GetPredictedChi2(const IlcCluster *cluster2) const;
     Double_t GetPredictedChi2(const IlcCluster *cluster2,
     			       double rwire,double anglewire,double tanstereo,bool withZ) const;
     Bool_t Update(const IlcCluster* c2, Double_t chi2, Int_t i);
     Bool_t Update(const IlcCluster* c2, Double_t chi2, Int_t i,
		   double rwire,double anglewire,double tanstereo,bool withZ,int signImP=0,bool calcChi2=false);

  // seeding operation
     Int_t Compare(const TObject *o) const;
     void Reset(Bool_t all = kTRUE);
     IlcDCHTrackerPoint * GetTrackPoint(Int_t i,Int_t j);
     IlcDCHcluster * GetClusterFast(Int_t ilayer,int j){ return fClusterPointer[ilayer][j];}
     int GetClusterPointerI(Int_t ilayer, IlcDCHcluster* cl){
       for(int i=0;i<kMaxInLayer;i++)
	 if(fClusterPointer[ilayer][i]==cl||fClusterPointer[ilayer][i]==0) return i;
       return kMaxInLayer-1;
     };
     bool IsAlready(Int_t ilayer, IlcDCHcluster* cl){
       for(int i=0;i<kMaxInLayer;i++)
	 if(fClusterPointer[ilayer][i]==cl) return true;
       return false;
     }
     void SetClusterPointer(Int_t ilayer, IlcDCHcluster* cl,int j) {
	 fClusterPointer[ilayer][j]=cl;
     }
     void RebuildSeed(); // rebuild seed to be ready for storing
     Double_t GetDensityFirst(Int_t n);
     void GetClusterStatistic(Int_t first, Int_t last, Int_t &found, Int_t &foundable, Int_t &shared, Bool_t plus2);
     
     void Modify(Double_t factor);

     void SetErrorY2(Float_t sy2){fErrorY2=sy2;}
     void SetErrorZ2(Float_t sz2){fErrorZ2=sz2;}
     Float_t  CookdEdx(Double_t low=0.05, Double_t up=0.70, Int_t i1=0, Int_t i2=kMaxLayer, Bool_t onlyused = kFALSE);
     void CookPID();
     Double_t Bethe(Double_t bg);     // return bethe-bloch

     Bool_t IsActive() const { return !(fRemoval);}
     void Desactivate(Int_t reason){ fRemoval = reason;} 
     IlcDCHcluster* GetClusterPointer(Int_t i,int j) const {return fClusterPointer[i][j];}
     Float_t GetClusterChi2(Int_t i,int j) const {return fClusterChi2[i][j];}
     void SetClusterChi2(Int_t i,float chi2,int j) {fClusterChi2[i][j]=chi2;}
     Char_t GetCircular() const {return fCircular;}

     void SetCircular(Char_t c) {fCircular=c;}
     void SetIsSeeding(Bool_t s) {fIsSeeding=s;}
     void SetSeedType(Int_t s) {fSeedType=s;}
     void SetSeed1(Int_t s) {fSeed1=s;}
     void SetSeed2(Int_t s) {fSeed2=s;}
  //     void SetESD(IlcESDtrack* esd) {fEsd=esd;}
     void SetBSigned(Bool_t s) {fBSigned=s;}
     void SetSort(Int_t s) {fSort=s;}
     void SetOverlapLabel(Int_t i, Int_t l) {fOverlapLabels[i]=l;}
     void SetCurrentCluster(IlcDCHcluster* cl,Int_t clindex) {fCurrentCluster=cl;fCurrentClusterIndex=clindex;}
     void SetLayer(Int_t n) {fLayer=n;}
     void SetInDead(Bool_t s) {fInDead=s;}
     void SetPoints(TClonesArray* p) {fPoints=p;}
     void SetEPoints(TClonesArray* p) {fEPoints=p;}

     Double_t DCHrPID(Int_t i) const {return fDCHr[i];}
     Double_t* DCHrPIDs() {return fDCHr;}
     Bool_t GetIsSeeding() const {return fIsSeeding;}
     Int_t GetSeedType() const {return fSeedType;}
     Int_t GetSeed1() const {return fSeed1;}
     Int_t GetSeed2() const {return fSeed2;}
  //     IlcESDtrack* GetESD() {return fEsd;}
     Float_t GetSDEDX(Int_t i) const {return fSDEDX[i];}
     Int_t GetNCDEDX(Int_t i) const {return fNCDEDX[i];}
     Bool_t GetBSigned() const {return fBSigned;}
     Int_t GetSort() const {return fSort;}
     Int_t GetOverlapLabel(Int_t i) const {return fOverlapLabels[i];}
     IlcDCHcluster* GetCurrentCluster() const {return fCurrentCluster;}
     Int_t GetLayer() const {return fLayer;}
     Int_t GetCurrentClusterIndex() const {return fCurrentClusterIndex;}
     Bool_t GetInDead() const {return fInDead;}
     Float_t GetErrorY2() const {return fErrorY2;}
     Float_t GetErrorZ2() const {return fErrorZ2;}
     //
     //
 private:
     //     IlcDCHseed & operator = (const IlcDCHseed &)
     //  {::Fatal("= operator","Not Implemented\n");return *this;}
  //     IlcESDtrack * fEsd; //!
     IlcDCHcluster*   fClusterPointer[kMaxLayer][kMaxInLayer];  // array of cluster pointers  - 
     Float_t  fClusterChi2[kMaxLayer][kMaxInLayer];  // array of cluster chi2  - 
     Bool_t             fClusterOwner;         // indicates the track is owner of cluster
     TClonesArray * fPoints;              //!array with points along the track
     TClonesArray * fEPoints;             //! array with exact points - calculated in special macro not used in tracking
     //---CURRENT VALUES
     Int_t fLayer;                 //!current layer number  
     Float_t fErrorY2;           //!sigma of current cluster 
     Float_t fErrorZ2;           //!sigma of current cluster    
     IlcDCHcluster * fCurrentCluster; //!pointer to the current cluster for prolongation
     Int_t   fCurrentClusterIndex; //! index of the current cluster
     Bool_t  fInDead;            //! indicate if the track is in dead zone
     Bool_t  fIsSeeding;         //!indicates if it is proces of seeading
     Int_t   fSort;              //!indicate criteria for sorting
     Bool_t  fBSigned;        //indicates that clusters of this trackes are signed to be used
     //
     //
     Float_t fDEDX[4];         // dedx according padlayers
     Float_t fSDEDX[4];        // sdedx according padlayers
     Int_t   fNCDEDX[4];       // number of clusters for dedx measurment
     Double_t fDCHr[IlcPID::kSPECIES];   // rough PID according DCH   
     //
     Int_t   fSeedType;         //seeding type
     Int_t   fSeed1;            //first layer for seeding
     Int_t   fSeed2;            //last layer for seeding
     Int_t   fOverlapLabels[12];  //track labels and the length of the  overlap     
     Float_t fMAngular;           // mean angular factor
     Char_t   fCircular;           // indicates curlin track
     IlcDCHTrackerPoint  fTrackPoints[kMaxLayer][kMaxInLayer];  //track points - array track points
     ClassDef(IlcDCHseed,1)  
};





#endif



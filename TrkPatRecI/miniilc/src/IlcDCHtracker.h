#ifndef ILCDCHTRACKER_H
#define ILCDCHTRACKER_H
/* Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 * See cxx source for full Copyright notice                               */


/* $Id: IlcDCHtracker.h,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */

//-------------------------------------------------------
//                       DCH tracker
//   Parallel tracker 
//
//   Origin: 
//-------------------------------------------------------

#include <TError.h>
#include "IlcTracker.h"
#include "IlcDCHreco.h"
#include "IlcDCHcluster.h"
#include "IlcExternalTrackParam.h"
#include "IlcLog.h"

class TFile;
class IlcDCHrecoParam;
class IlcDCHseed;
class IlcDCHTrackerPoint;
//class IlcESD;   
class TTree;
//class IlcESDkink;
class TTreeSRedirector;
class IlcTrackPoint;
class IlcDCHwireposition;
class IlcDCHParam;
class IlcDCHcluster;

class IlcDCHtracker : public IlcTracker {
public:
  IlcDCHtracker():IlcTracker(){
    fSeeds=0; 
  }
  IlcDCHtracker(const IlcDCHParam *par,IlcDCHwireposition *wpos=0); 
  virtual ~IlcDCHtracker();
  void SetUseZCoordinate(){kUseZCoordinate=true;}
  //
  void SetIteration(Int_t iteration){fIteration = iteration;}
  // virtual Int_t Clusters2Tracks (IlcESD *esd);
  // virtual Int_t RefitInward (IlcESD *esd);
  virtual Int_t LoadClusters (TObjArray* array);
  virtual Int_t LoadClusters (TTree *clusterTree);
  Int_t ReadClusters(TObjArray *array, TTree *clusterTree) const;
  void   UnloadClusters();

  void   Transform(IlcCluster * cluster);
  //
  // void FillESD(TObjArray* arr);
  void DeleteSeeds();
  void SetDebug(Int_t debug){ fDebug = debug; IlcLog::SetGlobalDebugLevel(debug); if (debug==0) IlcLog::SetGlobalLogLevel(IlcLog::kFatal); }
  // void FindKinks(TObjArray * array, IlcESD * esd);
  // void FindV0s(TObjArray * array, IlcESD * esd);
  // void UpdateKinkQuilctyM(IlcDCHseed * seed);
  // void UpdateKinkQuilctyD(IlcDCHseed * seed);
  // Int_t CheckKinkPoint(IlcDCHseed*seed, IlcDCHseed &mother, IlcDCHseed &daughter, IlcESDkink &kink);
  // Int_t RefitKink(IlcDCHseed &mother, IlcDCHseed &daughter, IlcESDkink &kink);
  Int_t ReadSeeds(const TFile *in);
  TObjArray * GetSeeds(){return fSeeds;}
   //   
  IlcDCHcluster *GetCluster(Int_t index) const;
  Int_t Clusters2Tracks();
  virtual void  CookLabel(IlcKalmanTrack *t,Float_t wrong) const { CookLabel(t,wrong,0,kMaxLayer-1);}
  virtual Int_t CookLabel(IlcKalmanTrack *t,Float_t wrong, Int_t first, Int_t last) const; 
   
  Int_t FollowProlongation(IlcDCHseed& t, Int_t rf=0, Int_t step=1,int contin=0);
   Bool_t GetTrackPoint(Int_t index, IlcTrackPoint &p ) const; 

   Int_t FollowBackProlongation(IlcDCHseed& t, Int_t rf);
   Int_t FollowToNext(IlcDCHseed& t, Int_t nr);
  Int_t FindAndUpdate(IlcDCHseed& t, Int_t nr,double roady,double roadz);
   Int_t UpdateClusters(IlcDCHseed& t,  Int_t nr);
   Int_t FollowToNextCluster( IlcDCHseed& t, Int_t nr);

   Int_t PropagateBack(TObjArray * arr);
   // Int_t PropagateBack(IlcESD * event);
   Int_t PropagateForward();
   Int_t PropagateForward2(TObjArray * arr);

   void SortTracks(TObjArray * arr, Int_t mode) const;
  
   Double_t F1(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3); 
   Double_t F1old(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3); 
   Double_t F2(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3); 
   Double_t F2old(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3); 

   Double_t F3(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t z1,Double_t z2); 
   Double_t F3n(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t z1,Double_t z2, 
                Double_t c); 
   Bool_t GetProlongation(Double_t x1, Double_t x2, Double_t x[5], Double_t &y, Double_t &z);

 public:
// **************** Internal tracker class ********************** 
   class IlcDCHSector;
   class IlcDCHLayer {
   public:
     IlcDCHLayer();
     ~IlcDCHLayer();
     void InsertCluster(const IlcDCHcluster *c, UInt_t index);
     void ResetClusters();
     operator int() const {return fN;}
     void SetN(int n) {fNarray=n;}
     Int_t GetN() const {return fNarray;}
     const IlcDCHcluster* operator[](Int_t i) const {return fClusters[i];}
     UInt_t GetIndex(Int_t i) const {return fIndex[i];}

     inline Int_t Find(Double_t z) const; 
     IlcDCHcluster *  FindNearest(Double_t phi, Double_t z, Double_t roadphiR, Double_t roadz) const;
     IlcDCHcluster *  FindNearest2(Double_t phi,Double_t snp, Double_t z, Double_t roadphiR, Double_t roadz, UInt_t & index, int first=0) const;
     
     void FindNearest2(Double_t phi,Double_t snp, Double_t z, Double_t roadphiR, Double_t roadz,int first,IlcDCHcluster* clar[],UInt_t inar[]) const;
    
     void SetR(Double_t r) {fR=r;}
     void SetMaxZ(Double_t z) {fMaxZ=z;}
     void SetStereoAngle(Double_t stereoangle) {
       fTStereoAngle=TMath::Tan(stereoangle);
       if(fabs(fTStereoAngle)>1e-6) fIsStereo=true;
     }
     void SetUseZ(bool use) {fIsUseZ=use;}
     void SetFirstPhi(double phi) {fFirstPhi=phi;}
     void SetNWires  (double nw ) {fNWires  =nw ;}  
     void SetDeltaPhi(double dph) {fDeltaPhi=dph;}
     bool IsStereo()           const {return fIsStereo;}
     Double_t GetR()           const {return fR;}
     Double_t GetMaxZ()        const {return fMaxZ;}
     Double_t GetTStereoAngle() const {return fTStereoAngle;}
     Double_t GetStereoAngle()  const {return TMath::ATan(fTStereoAngle);}
     Double_t GetFirstPhi()     const {return fFirstPhi;}
     Double_t GetNWires  ()     const {return fNWires  ;}  
     Double_t GetDeltaPhi()     const {return fDeltaPhi;}
     
     Double_t GetWirePhi(int nw,double z=0)const {
       return (fFirstPhi+nw*fDeltaPhi)+TMath::ATan(z/fR*fTStereoAngle);
     }
     Double_t GetRAtZ(double z=0)const {
       return TMath::Hypot(z*fTStereoAngle,fR);
     }
     
     bool IsCrossedCell(double phi,double z,double snp,double tgl) const{
       return true;
       //below forr hexognal cells chamber
       double dphiw=TMath::ATan(z/fR*fTStereoAngle);
       double nw=Phi_0_2pi(phi-fFirstPhi-dphiw)/fDeltaPhi*3;
       if((int(nw)%3)==1){
	 double lphicell=fDeltaPhi/3*TMath::Sin(TMath::Pi()/3);
	 double nw2=(snp*fDeltaPhi/3-TMath::ATan((z+tgl*lphicell*fR)/fR*fTStereoAngle)+dphiw)/fDeltaPhi*3;
	 if(int((nw+nw2)/3+1)==int(nw/3+1)&&(int(nw+nw2+0.5)%3)!=0){
	   nw2=(-snp*fDeltaPhi/3-TMath::ATan((z-tgl*lphicell*fR)/fR*fTStereoAngle)+dphiw)/fDeltaPhi*3;
	   if(int((nw+nw2)/3+1)==int(nw/3+1)&&(int(nw+nw2+0.5)%3)!=0)
	     return false;
	 }
       }
       return true;
     }
     

     IlcDCHcluster* GetClusters() const {return fClustersArray;}
     void SetClusters(IlcDCHcluster* cl) {fClustersArray=cl;}
     void SetCluster(Int_t i, IlcDCHcluster cl) {fClustersArray[i]=cl;}
     IlcDCHcluster* GetCluster(Int_t i) const {return &fClustersArray[i];}
     Short_t GetFastClusterZ(Int_t i) const 
     {int ii=i;if(i<0) ii=0;if(i>=NFAST) ii=NFAST-1;return fFastClusterZ[ii];}
     void SetFastClusterZ(Int_t i, Short_t cl) 
     {int ii=i;if(i<0) ii=0;if(i>=NFAST) ii=NFAST-1;fFastClusterZ[ii]=cl;}
     Double_t GetCoeffConv() const {return fCoeffConv;}
     void SetCoeffConv(double coeff) {fCoeffConv=coeff;}

     const static int NFAST=510;
   private:  
     IlcDCHLayer & operator=(const IlcDCHLayer & );
     IlcDCHLayer(const IlcDCHLayer& /*r*/);           //dummy copy constructor
     Short_t fFastClusterZ[NFAST];   //index of the nearest cluster at given position
     Short_t fFastClusterPhi[NFAST];   //index of the nearest cluster at given position
     Double_t fCoeffConv;   //converstion coeff from Z to int(Z*fCoeffConv+255)
     Int_t fN;                                          //number of inserted clusters 
     const IlcDCHcluster *fClusters[kMaxClusterPerLayer]; //pointers to clusters
     Int_t fNarray;                                          //number of setted clusters 
     IlcDCHcluster *fClustersArray;                     // 
     UInt_t fIndex[kMaxClusterPerLayer];                  //indeces of clusters
     Double_t fR;                                 //radius of this layer
     Bool_t   fIsStereo;                          // is it a stereo layer
     Double_t fTStereoAngle;                       //stereo angle of this layer
     Double_t fMaxZ;                              //maximum Z of this layer
     Double_t fFirstPhi;                          //phi of wire id=0
     Double_t fNWires;                            // number of wires at layer
     Double_t fDeltaPhi;                          //dphi of wires 
     bool     fIsUseZ;                            //is it use or not use Z-coordinate
   };

// **************** Internal tracker class ********************** 
   class IlcDCHSector {
   public:
     IlcDCHSector() { fN=0; fLayer = 0; }
     ~IlcDCHSector() { delete[] fLayer; }
     IlcDCHLayer& operator[](Int_t i) const { return *(fLayer+i); }
     Int_t GetNLayers() const { return fN; }
     void Setup(const IlcDCHParam *par,IlcDCHwireposition *wirepos);
     Double_t GetR(Int_t l) const {return fLayer[l].GetR();}
     Double_t GetMaxZ(Int_t l) const { return fLayer[l].GetMaxZ();} 
     Int_t GetLayerNumber(Double_t  x,double z) const;
   private:
     IlcDCHSector & operator=(const IlcDCHSector & );
     IlcDCHSector(const IlcDCHSector &/*s*/);           //dummy copy contructor 
     Int_t fN;                        //number of pad layers 
     IlcDCHLayer *fLayer;                    //array of pad layers
     
   };

   Float_t OverlapFactor(IlcDCHseed * s1, IlcDCHseed * s2, Int_t &sum1, Int_t &sum2);
   void  SignShared(IlcDCHseed * s1, IlcDCHseed * s2);
   void  SignShared(TObjArray * arr);

   void  RemoveUsed(TObjArray * arr, Float_t factor1, Float_t factor2,  Int_t removilcndex);
   void  RemoveUsed2(TObjArray * arr, Float_t factor1, Float_t factor2, Int_t minimal);
   void  RemoveDouble(TObjArray * arr, Float_t factor1, Float_t factor2,  Int_t removilcndex);

   void  StopNotActive(TObjArray * arr, Int_t layer0, Float_t th0, Float_t th1, Float_t th2) const;
   void  StopNotActive(IlcDCHseed * seed, Int_t layer0, Float_t th0, Float_t th1, Float_t th2) const;
   Int_t AcceptCluster(IlcDCHseed * seed, IlcDCHcluster * cluster, Float_t factor, Float_t cory=1., Float_t corz=1.);

private:
  friend class IlcDCHtrack;

  IlcDCHtracker(const IlcDCHtracker& r);           //dummy copy constructor
  IlcDCHtracker &operator=(const IlcDCHtracker& r);//dummy assignment operator
  inline IlcDCHLayer &GetLayer(Int_t layer){return fSector[layer];};
  inline Bool_t     IsActive(Int_t /*layer*/){return kTRUE;/*fSector[layer].GetN()>0;*/};
  inline Double_t  GetRlayer(Int_t layer,double z=0) const{return fSector[layer].GetRAtZ(z);};
  inline Int_t GetLayerNumber(Double_t x,Double_t z) const{return fSector.GetLayerNumber(x,z);};
  Int_t GetLayerFromId(ULong_t clid) const;

   void MakeSeeds3(TObjArray * arr, Int_t i1, Int_t i2, Float_t cuts[4]); 
   void MakeSeeds3Stereo(TObjArray * arr, Int_t i1, Int_t i2, Float_t cuts[4]); 
   void MakeSeeds3StereoLine(TObjArray * arr, Int_t i1, Int_t i2, Float_t cuts[4]); 
   void MakeSeeds5(TObjArray * arr, Int_t i1, Int_t i2, Float_t cuts[4]);
   void MakeSeeds2(TObjArray * arr, Int_t i1, Int_t i2, Float_t cuts[4]);
  

   // void ReadSeeds(IlcESD *event, Int_t direction);  //read seeds from the event
   IlcDCHseed *MakeSeed(IlcDCHseed *t, Float_t r0, Float_t r1, Float_t r2); //reseed
   IlcDCHseed *ReSeed(IlcDCHseed *t, Float_t r0, Float_t r1, Float_t r2); //reseed
   IlcDCHseed *ReSeed(IlcDCHseed *t, Int_t r0, Bool_t forward); //reseed


  
   IlcDCHseed * ReSeed(IlcDCHseed *t);
   void UnsignClusters();
   void SignClusters(TObjArray * arr, Float_t fnumber=3., Float_t fdensity=2.);  

   void ParallelTracking(TObjArray * arr, Int_t rfirst, Int_t rlast);
   void Tracking(TObjArray * arr);
   TObjArray * Tracking(Int_t seedtype, Int_t i1, Int_t i2, Float_t cuts[4]);
   TObjArray * Tracking();
   TObjArray * TrackingSpecial();
   void SumTracks(TObjArray *arr1,TObjArray *arr2) const;
   void PrepareForBackProlongation(TObjArray * arr, Float_t fac) const;
   void PrepareForProlongation(TObjArray * arr, Float_t fac) const;

   void SetSampledEdx(IlcDCHseed */*t*/, Float_t /*q*/, Int_t /*i*/) {;}
   Int_t UpdateTrack(IlcDCHseed *t, Int_t accept); //update trackinfo


   IlcDCHSector fSector;  //pointer to loaded sectors;
   //
   // IlcESD * fEvent;      // output with esd tracks
   Int_t    fDebug;      // debug option        
   Int_t fNtracks;                     //current number of tracks
   TObjArray *fSeeds;                  //array of track seeds
   Int_t fIteration;                   // indicate iteration - 0 - flayerard -1 back - 2forward - back->forward
   
   const IlcDCHParam *fParam;          //pointer to the parameters
   const IlcDCHrecoParam *fRecPar;      //pointer to tracker parameters
   TTreeSRedirector *fDebugStreamer;     //!debug streamer
   bool kUseZCoordinate;               // to use or not use Z-ccordinate mesuarement of wire

   ClassDef(IlcDCHtracker,1) 
};

#endif


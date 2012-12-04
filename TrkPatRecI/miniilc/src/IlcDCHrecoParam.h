#ifndef ILCDCHRECOPARAMFI_H
#define ILCDCHRECOPARAMFI_H
#include "TObject.h"

class IlcDCHrecoParam:public TObject{
public:
  IlcDCHrecoParam();
  void SetRoadY(double roady)      {fRoadY=roady;}
  void SetRoadZ(double roadz)      {fRoadZ=roadz;}
  void SetSafetyZ(double safetyz)  {fSafetyZ=safetyz;}
  double GetRoadY()            const{return fRoadY;}
  double GetRoadZ()            const{return fRoadZ;}
  double GetSafetyZ()          const{return fSafetyZ;}
  void SetNearestByC       (double fnearestbyC  )  {fNearestByC       =fnearestbyC  ;}
  void SetNearestByTgl     (double fnearestbyTgl)  {fNearestByTgl     =fnearestbyTgl;}
  void SetNearestByDistance2(double fnearestbyDst2)  {fNearestByDistance2=fnearestbyDst2;}
  double GetNearestByC()       const{ return fNearestByC;  }          //
  double GetNearestByTgl()     const{return fNearestByTgl;}        //
  double GetNearestByDistance2()const{return fNearestByDistance2;}   //
  void SetMinNOfClusetrs(int minNC){fMinNOfClusetrs=minNC;}    
  int GetMinNOfClusetrs()const{return  fMinNOfClusetrs;}  
  void SetMinDensity(double den){fMinDensity=den;}
  double GetMinDensity() const{return fMinDensity;} 
  double SetSeedTrackSigma2(double seedsigma2){return  fSeedTrackSigma2=seedsigma2;}
  double GetSeedTrackSigma2()const{return  fSeedTrackSigma2;}
  double SetSeedBeamSigma2(double seedbeamsigma2){return  fSeedBeamSigma2=seedbeamsigma2;}
  double GetSeedBeamSigma2()const{return  fSeedBeamSigma2;}
  double GetSeedBeamSigmaZ2()const{return  fSeedBeamSigmaZ2;}
  void   SetCutTrackSigma2(double cutsigma2) {fCutTrackSigma2=cutsigma2;};
  double GetCutTrackSigma2() const{return fCutTrackSigma2;};
  void SetMinDCAr(double cut) {fMinDCAr=cut;};
  void SetMinDCAz(double cut) {fMinDCAz=cut;};
  double GetMinDCAr() const {return fMinDCAr;};
  double GetMinDCAz() const {return fMinDCAz;};
  void SetScaleForQuilcty(double dr) {fScaleForQuilcty=dr;};
  double GetScaleForQuilcty() const {return fScaleForQuilcty;};
  void SetDeltaZ(double dz) {fDeltaZ=dz;};
  double GetDeltaZ() const {return fDeltaZ;};
  void SetDistFromEdge(double dist)  {fDistFromEdge=dist;};
  double GetDistFromEdge() const {return fDistFromEdge;};
  double GetMinRofV0         () const {return fMinRofV0         ;}; 
  double GetMaxRofV0         () const {return fMaxRofV0         ;}; 
  double GetCutOnDistanceInV0() const {return fCutOnDistanceInV0;}; 
  void SetMinRofV0         (double rmin) {fMinRofV0=rmin     ;}; 
  void SetMaxRofV0         (double rmax) {fMaxRofV0=rmax     ;}; 
  void SetCutOnDistanceInV0(double dr)   {fCutOnDistanceInV0=dr;}; 
protected:
  double fRoadY; // road forr searching clusters
  double fRoadZ;
  double fSafetyZ; // Aditional Z forr good region
  // if (TMath::Abs(z)<(IlcDCHReconstructor::GetCtgRange()*t.GetX()+fSafetyZ)) t.fNFoundable++;

  //  Distance to check overlapping o tracks 
  double fNearestByC;          // Delta Curvature 
  double fNearestByTgl;        // Delta Tgl 
  double fNearestByDistance2;  // dy^2+dz^2
  // number of clusters
  int fMinNOfClusetrs;    // Minimum Number of Clusters forr acceptable tracks 
  double fMinDensity;      //minimum acceptable dencity of clusters in track
  // 
  double fSeedTrackSigma2; // default resolution of track position after make seed with 3 clusters
  double fSeedBeamSigma2; // default resolution addition to beam resolution during make seed with 2 clusters and beam point
                          // 25000*Curv*Curv+fSeedBeamSigma2
  double fSeedBeamSigmaZ2; // default resolution addition to beam resolution during make seed with 2 clusters and beam point
  double fCutTrackSigma2; // cut on resolution of track position after make seed

  //kinks
  double fMinDCAr; // cut on r near beam point when looking forr circular track
  double fMinDCAz; // cut on z near beam point when looking forr circular track
  double fDeltaZ; // cut on z in middle of DCH when merge tracks
  double fDistFromEdge; // distance from inner and outer radius of DCH forr radius of kink
  double fScaleForQuilcty; // dAngle/(fScaleForQuilcty+abs(Rkink-Rtrack))

  //v0
  double fMinRofV0          ; // Minimum radius of V0
  double fMaxRofV0          ; // Maximum radius of V0
  double fCutOnDistanceInV0 ; //Cut On Distance In V0 of tracks
  ClassDef(IlcDCHrecoParam,1);
};

#endif

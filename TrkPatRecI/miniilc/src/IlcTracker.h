#ifndef ILCTRACKER_H
#define ILCTRACKER_H
/* Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 * See cxx source for full Copyright notice                               */

/* $Id: IlcTracker.h,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */

//-------------------------------------------------------------------------
//                          class IlcTracker
//   that is the base for IlcTPCtracker, IlcVXDtrackerV2 and IlcTRDtracker
//       Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------
#include <TObject.h>

// class IlcMagF;
class IlcCluster;
class TTree;
class IlcKalmanTrack;
// class IlcESD;
class IlcTrackPoint;

class IlcTracker : public TObject {
public:
  IlcTracker();
  virtual ~IlcTracker(){}
  // virtual Int_t Clusters2Tracks(IlcESD *event)=0;
  // virtual Int_t PropagateBack(IlcESD *event)=0;
  // virtual Int_t RefitInward(IlcESD *event)=0;
  void SetVertex(const Double_t *xyz, const Double_t *ers=0) { 
     fX=xyz[0]; fY=xyz[1]; fZ=xyz[2];
     if (ers) { fSigmaX=ers[0]; fSigmaY=ers[1]; fSigmaZ=ers[2]; } 
  }

//protected:
  virtual Int_t LoadClusters(TTree *)=0;
  virtual void UnloadClusters()=0;
  virtual IlcCluster *GetCluster(Int_t index) const=0;
  virtual Bool_t GetTrackPoint(Int_t /* index */ , IlcTrackPoint& /* p */) const { return kFALSE;}
  virtual void  UseClusters(const IlcKalmanTrack *t, Int_t from=0) const;
  virtual void  CookLabel(IlcKalmanTrack *t,Float_t wrong) const; 
  Double_t GetX() const {return fX;}
  Double_t GetY() const {return fY;}
  Double_t GetZ() const {return fZ;}
  Double_t GetSigmaX() const {return fSigmaX;}
  Double_t GetSigmaY() const {return fSigmaY;}
  Double_t GetSigmaZ() const {return fSigmaZ;}

  // static void SetFieldMap(const IlcMagF* map, Bool_t uni);
  // static const IlcMagF *GetFieldMap() {return fgkFieldMap;}
  static Double_t GetBz(Float_t *r); 
  static Double_t GetBz(Double_t *r) {
    Float_t rr[]={r[0],r[1],r[2]};
    return GetBz(rr);
  }
  static Double_t GetBz() {return fgBz;}
  static void SetBz(double field) {fgBz=field;}
  static Bool_t UniformField() {return fgUniformField;}

protected:
  IlcTracker(const IlcTracker &atr);
private:
  IlcTracker & operator=(const IlcTracker & atr);

  static Bool_t fgUniformField;       // uniform field flag
  // static const IlcMagF *fgkFieldMap;  // field map
  static Double_t fgBz;               // Nominal Bz (kG)

  Double_t fX;  //X-coordinate of the primary vertex
  Double_t fY;  //Y-coordinate of the primary vertex
  Double_t fZ;  //Z-coordinate of the primary vertex
 
  Double_t fSigmaX; // error of the primary vertex position in X
  Double_t fSigmaY; // error of the primary vertex position in Y
  Double_t fSigmaZ; // error of the primary vertex position in Z

  ClassDef(IlcTracker,3) //abstract tracker
};

#endif



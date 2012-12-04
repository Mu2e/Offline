#ifndef ILCDCHCLUSTER_H
#define ILCDCHCLUSTER_H
/* Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 * See cxx source for full Copyright notice                               */

/* $Id: IlcDCHcluster.h,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */

//-------------------------------------------------------
//                    DCH Cluster Class
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------

#include "TMath.h"
#include "IlcCluster.h"

//_____________________________________________________________________________
class IlcDCHcluster : public IlcCluster {
public:
  IlcDCHcluster():IlcCluster(){
    fQ=0.0; fImP=-1.0; fSigmaYZ=0.0;fUsed=0;fZLocal=-1e30; fImPexact=-1.; fZLocalexact=0; fTrackTime=0.;
  }
  IlcDCHcluster(Int_t *lab, Float_t *hit) : IlcCluster(lab,hit) {
    fImP = hit[4];fUsed=0;fZLocal=-1e30; fImPexact=-1.; fZLocalexact=0; fTrackTime=0.;
  }

  void SetQ(Float_t q)			{fQ=q;}

  void SetImP(Float_t b)		{fImP=b;}
  void SetImPexact(Float_t b)		{fImPexact=b;}
  void SetSigmaImP2(Float_t simp2)      {fSigmaImP2=simp2;}
  void SetSigmaYZ(Float_t syz)		{fSigmaYZ=syz;}
  void SetLabelAll(Int_t *lab)		{fTracks[0]=lab[0]; fTracks[1]=lab[1]; fTracks[2]=lab[2];}
  virtual void SetX(Float_t x)		{fX=x;}
  void SetSigmaX2(Float_t sx2)		{fSigmaX2=sx2;}
  void SetZLocal(Float_t zl)		{fZLocal=zl;}
  void SetZLocalexact(Float_t zl)	{fZLocalexact=zl;}
  void SetTrackTime(Float_t tm)		{fTrackTime=tm;}
  void SetId(ULong_t lid)	  	{fId=lid;}
  void SetIndex(ULong_t lid)	  	{fIndex=lid;}
   
  Float_t GetQ()		const {return TMath::Abs(fQ);}
  Float_t GetImP()		const {return fImP;}
  Float_t GetImPexact()		const {return fImPexact;}
  Float_t GetSigmaImP2()	const {return fSigmaImP2;}
  Float_t GetSigmaYZ()		const {return fSigmaYZ;}
  virtual Float_t GetX()	const {return fX;}
  Float_t GetSigmaX2()		const {return fSigmaX2;}
  Float_t GetZLocal()		const {return (fZLocal<-1e-29?fZ:fZLocal);}
  Float_t GetZLocalexact()	const {return fZLocalexact;}
  Float_t GetTrackTime()	const {return fTrackTime;}
  
  Int_t  GetSL()  		const {return fId/100000; }
  Int_t  GetR()   		const {return (fId%100000)/1000; }
  Int_t  GetW()       		const {return fId%1000; }
  ULong_t GetId()		const {return fId;}
  ULong_t GetIndex()		const {return fIndex;}

  inline  void Use(Int_t inc=10){if (inc>0) fUsed+=inc; else fUsed=0;};
  Int_t IsUsed(Int_t th=10) const {return (fUsed>=th) ? 1 : 0;}
  virtual Bool_t IsSortable() const {return kTRUE;}; 
  static void SetSortType(int sort){fSortType=sort;}
  virtual Int_t Compare(const TObject* obj) const {
    if(fSortType==0)
      return (((IlcDCHcluster*)obj)->GetId()>fId)? -1:1; 
    return (((IlcDCHcluster*)obj)->GetTrackTime()>fTrackTime)? -1:1;
  }
 

private:
  Float_t   fQ ;	//Q of cluster (in ADC counts)
  Float_t   fImP ;	//track impact parameter (track distance from wire)
  Float_t   fImPexact ;	//track impact parameter (track distance from wire)
  Float_t   fSigmaImP2;	//Sigma Impact Parameter
  Float_t   fSigmaYZ;	//Sigma YZ mixing
  Float_t   fX ;	//X of cluster
  Float_t   fSigmaX2;	//Sigma X square of cluster
  Float_t   fZLocal ;	//Z along of wire
  Float_t   fZLocalexact ;// along of wire
  Float_t   fTrackTime ;// track time
  ULong_t   fId;
  ULong_t   fIndex;    //index in global array
  Char_t    fUsed;     //counter of usage  
  static int fSortType; //! sort type

  ClassDef(IlcDCHcluster,3)  // CLUCOU Drift Chamber clusters
};

#endif



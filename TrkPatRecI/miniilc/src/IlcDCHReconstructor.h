#ifndef ILCDCHRECONSTRUCTOR_H
#define ILCDCHRECONSTRUCTOR_H
/* Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 * See cxx source for full Copyright notice                               */

/* $Id: IlcDCHReconstructor.h,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */

//#include "IlcReconstructor.h"
#include "IlcDCHRecoParam.h"

class IlcDCHParam;


class IlcDCHReconstructor/*: public IlcReconstructor*/ {
public:
  IlcDCHReconstructor();
  virtual ~IlcDCHReconstructor() {if (fgkRecoParam) delete fgkRecoParam;};

  /*  virtual void         Reconstruct(IlcRunLoader* runLoader) const;
  virtual void         Reconstruct(IlcRunLoader* runLoader,
				   IlcRawReader* rawReader) const;
  virtual void         Reconstruct(TTree* digitsTree, TTree* clustersTree) const {
    IlcReconstructor::Reconstruct(digitsTree,clustersTree);
  }
  virtual void         Reconstruct(IlcRawReader* rawReader, TTree* clustersTree) const {
    IlcReconstructor::Reconstruct(rawReader,clustersTree);
  }

  virtual IlcTracker*  CreateTracker(IlcRunLoader* runLoader) const;
  virtual void         FillESD(IlcRunLoader* runLoader, IlcESD* esd) const;
  virtual void         FillESD(TTree* digitsTree, TTree* clustersTree, 
                               IlcESD* esd) const {
    IlcReconstructor::FillESD(digitsTree,clustersTree,esd);
  }
  virtual void         FillESD(IlcRawReader* rawReader, TTree* clustersTree, 
                               IlcESD* esd) const {
    IlcReconstructor::FillESD(rawReader,clustersTree,esd);
  }
  virtual void         FillESD(IlcRunLoader* runLoader, 
                               IlcRawReader* rawReader, IlcESD* esd) const {
    IlcReconstructor::FillESD(runLoader,rawReader,esd);
  }
  */
//   void SetRecoParam(IlcDCHRecoParam * param){ fgkRecoParam = param;}
//   static const IlcDCHRecoParam* GetRecoParam(){ return fgkRecoParam;}
  //
  static Double_t GetCtgRange()     { return fgkRecoParam->GetCtgRange();}
  static Double_t GetMaxSnpTracker(){ return fgkRecoParam->GetMaxSnpTracker();}
  static Double_t GetMaxSnpTrack()  { return fgkRecoParam->GetMaxSnpTrack();}
  static Double_t GetMaxC()         { return fgkRecoParam->GetMaxC();}
  static Bool_t   GetDoKinks()      { return fgkRecoParam->GetDoKinks();}

  static Int_t StreamLevel()               { return fgStreamLevel;}
  static void  SetStreamLevel(Int_t level) { fgStreamLevel = level;}

private:
  //IlcDCHParam*         GetDCHParam(/*IlcRunLoader* runLoader*/) const;
  static IlcDCHRecoParam *   fgkRecoParam; // reconstruction parameters
  static Int_t               fgStreamLevel; // flag for streaming      - for DCH reconstruction

  ClassDef(IlcDCHReconstructor, 0)   // class for the DCH reconstruction
};

#endif

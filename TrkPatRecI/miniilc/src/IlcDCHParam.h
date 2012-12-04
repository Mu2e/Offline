#ifndef ILCDCHPARAM_H
#define ILCDCHPARAM_H
/* Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 * See cxx source for full Copyright notice                               */

/* $Id: IlcDCHParam.h,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */

////////////////////////////////////////////////
//  Manager class for DCH parameters          //
////////////////////////////////////////////////

//#include "IlcDetectorParam.h"
#include <TNamed.h>
#include "TMath.h"

#include <TGeoMatrix.h>

//class IlcDCHParam : public IlcDetectorParam {
class IlcDCHParam : public TNamed  {
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //IlcDCHParam object to be possible change 
  //geometry and some other parameters of DCH   
  //used by IlcDCH and IlcDCHSector 
 
public:
  IlcDCHParam(); 
  virtual ~IlcDCHParam();
  

  virtual void SetDefault();          //set defaut TPCparam 

//--------------- DCH ----------------------------------

  void  SetInnerRadius(Float_t InnerRadius )  { fInnerRadius=InnerRadius;}
  void  SetOuterRadius(Float_t OuterRadius )  { fOuterRadius=OuterRadius;}
  void  SetInnerWallThickness(Float_t InnerWallThickness )  { fInnerWallThickness=InnerWallThickness;}
  void  SetOuterWallThickness(Float_t OuterWallThickness )  { fOuterWallThickness=OuterWallThickness;}
  void  SetEndCapWallThickness(Float_t EndCapWallThickness )  { fEndCapWallThickness=EndCapWallThickness;}
  void  SetFWireDiameter(Float_t FWireDiameter )  { fFWireDiameter=FWireDiameter;}
  void  SetSWireDiameter(Float_t SWireDiameter )  { fSWireDiameter=SWireDiameter;}
  void  SetEndCapType(Int_t EndCapType);
  void  SetExperiment(Int_t Experiment);
  void  SetExperimentSubVer(Int_t ExperimentSubVer);
  void  SetFWMaterialType(Int_t FWMaterialType);

  void  SetSWireNum(Int_t SWireNum) {fnum_wire_sense=SWireNum;}
  void  SetSDeltaWireNum(Int_t DeltaWireNum){ fdelta_num_wire_sense=DeltaWireNum;}
  void  SetSuperLayerNum(Int_t SuperLayerNum){fnsuperlayer=SuperLayerNum;}
  void  SetRingNum(Int_t RingNum) {fnring=RingNum;}
  void  SetDrop(Float_t Drop) {fdrop=Drop;}
  void  SetLength(Float_t Length) {flength=Length;}
  
  void  SetExtraEndCapDist(Float_t EndCapDist) {fextra_EndCap_dist=EndCapDist;}
  
  void  SetExtraDim(Float_t ExtraDim) {fextra_dim=ExtraDim;}
  void  SetEndCapWallThetaIn(Float_t EndCapWallThetaIn) {fEndCap_Wall_theta_inner=EndCapWallThetaIn;}
  void  SetEndCapWallThetaOut(Float_t EndCapWallThetaOut){fEndCap_Wall_theta_outer=EndCapWallThetaOut;}
  void  SetMaxEndCapDim(Float_t MaxEndCapDim){fMaxEndCapDim=MaxEndCapDim;}

//--------------- DCH ----------------------------------



//--------------- DCH ----------------------------------

  Float_t  GetInnerRadius() const {return fInnerRadius;}
  Float_t  GetOuterRadius() const {return fOuterRadius;} 
  Float_t  GetInnerWallThickness() const {return fInnerWallThickness;}
  Float_t  GetOuterWallThickness() const {return fOuterWallThickness;}
  Float_t  GetEndCapWallThickness() const {return fEndCapWallThickness;}
  Float_t  GetFWireDiameter() const {return fFWireDiameter;}
  Float_t  GetSWireDiameter() const {return fSWireDiameter;}
  Int_t    GetEndCapType() const {return fEndCapType;}
  Int_t    GetExperiment() const {return fExperiment;}
  Int_t    GetExperimentSubVer() const {return fExperimentSubVer;}
  Int_t    GetFWMaterialType() const {return fFWMaterialType;}
  
  Int_t    GetSWireNum() const {return fnum_wire_sense;}
  Int_t    GetNumberOfWires(int slay,int /*ring*/) const {
    Int_t nwires=GetSWireNum()+slay*GetSDeltaWireNum();
    
    if(slay>=10 && slay<18) nwires=GetSWireNum()+(slay-4)*GetSDeltaWireNum();
    if(slay>=18 && slay<24) nwires=GetSWireNum()+(slay-6)*GetSDeltaWireNum();
    
    return nwires;
  }
  Int_t    GetSDeltaWireNum()const { return fdelta_num_wire_sense;}
  Int_t    GetSuperLayerNum()const {return fnsuperlayer;}
  Int_t    GetRingNum() const {return fnring;}
  Int_t    GetNLayers() const {return (fnring-1)*fnsuperlayer;}
  Float_t  GetDrop() const {return fdrop;}
  Float_t  GetLength() const {return flength;}
  Float_t  GetExtraEndCapDist() const {return fextra_EndCap_dist;}
  Float_t  GetExtraDim() const { return fextra_dim;}
  Float_t  GetEndCapWallThetaIn() const {return fEndCap_Wall_theta_inner;}
  Float_t  GetEndCapWallThetaOut() const {return fEndCap_Wall_theta_outer;}
  Float_t  GetMaxEndCapDim() const {return fMaxEndCapDim;}


  Float_t  GetSigmaFrsTrk()    const {return 0.0055;}  
  Float_t  GetSigmaSecTrk()    const {return 0.01;}
  Float_t  GetTimeBetweenTrk() const {return 0.002;}  //in micro s
//  Float_t  GetStereoAngle()    const {return 0.120;}
  Float_t  GetMaxZResolution() const {return 0.045; /*0.014;*/}



//--------------- DCH ----------------------------------



protected :

//--------------- DCH ----------------------------------

  Float_t fInnerRadius;    // lower radius
  Float_t fOuterRadius;    // outer radius
  Float_t fInnerWallThickness;
  Float_t fOuterWallThickness;
  Float_t fEndCapWallThickness;
  Float_t fFWireDiameter;
  Float_t fSWireDiameter;
  Int_t   fEndCapType;
  Int_t   fExperiment;
  Int_t   fExperimentSubVer;
  Int_t   fFWMaterialType; // 0=Al 1=Kapton
  Int_t   fnum_wire_sense;
  Int_t   fdelta_num_wire_sense;
  Int_t   fnsuperlayer;
  Int_t   fnring;
  Float_t fdrop;
  Float_t flength;
  Float_t fextra_EndCap_dist;
  Float_t fextra_dim;
  Float_t fEndCap_Wall_theta_inner;
  Float_t fEndCap_Wall_theta_outer;
  Float_t fMaxEndCapDim;

//--------------- DCH ----------------------------------


  
private:
  IlcDCHParam(const IlcDCHParam &);
  IlcDCHParam & operator=(const IlcDCHParam &);
  ClassDef(IlcDCHParam,5)  //parameter  object for set:TPC
};

 



#endif  

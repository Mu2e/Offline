#ifndef ILCDCHGEOMETRY_H
#define ILCDCHGEOMETRY_H
/* Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 * See cxx source for full Copyright notice                               */ 
 
/* $Id: IlcDCHgeometry.h,v 1.2 2013/03/15 16:20:00 kutschke Exp $ */ 
 
/////////////////////////////////////////////////////////////////////////////// 
//                                                                           // 
//  DCH geometry class                                                       // 
//                                                                           // 
/////////////////////////////////////////////////////////////////////////////// 
 
// #include "IlcGeometry.h" 
#include "IlcDCHParam.h" 
#include "TObjArray.h"
#include "TVector3.h" 
#include "TMatrix.h" 
//class IlcRunLoader; 
class TGeoHMatrix; 
 
class IlcDCHgeometry : public TNamed /*: public IlcGeometry*/ { 
 
 public: 
 
  enum { kNplan = 4, kNcham = 3, kNsect =12 , kNdet = 144, kNdets = 5, kNlayers=10}; 
 
  IlcDCHgeometry(); 
  virtual ~IlcDCHgeometry(); 
 
  virtual void     CreateGeometry(Int_t *idtmed);
  virtual void     CreateGeometryCluCou(Int_t *idtmed);  
  virtual Int_t    IsVersion() { return 1; }; 
 virtual void     Init(); 
  //  virtual Bool_t   Impact(const TParticle* ) const { return kTRUE; }; 
  // virtual void     ImportGDMLGeometry(const char *filename="/afs/le.infn.it/user/g/giuter/ILC/Ilcv4-04-11/aaa.gdml"); 
 
//   virtual Bool_t   Local2Global(Int_t d, Double_t *local, Double_t *global) const; 
//   virtual Bool_t   Local2Global(Int_t p, Int_t c, Int_t s, Double_t *local, Double_t *global) const; 
//  virtual Bool_t   Global2Local(Int_t mode, Double_t *local, Double_t *global, Int_t* index) const; 
//  virtual Bool_t   Global2Detector(Double_t global[3], Int_t index[3]); 
//  virtual Bool_t   Rotate(Int_t d, Double_t *pos, Double_t *rot) const; 
//  virtual Bool_t   RotateBack(Int_t d, Double_t *rot, Double_t *pos) const; 
 
// 	  void     GroupChamber(Int_t iplan, Int_t icham, Int_t *idtmed); 
// 	  void     CreateFrame(Int_t *idtmed); 
// 	  void     CreateServices(Int_t *idtmed); 
 
//          Bool_t   ReadGeoMatrices();   
  TGeoHMatrix *    GetGeoMatrix(Int_t det)       { return (TGeoHMatrix *) fMatrixGeo->At(det);             } 
  TGeoHMatrix *    GetMatrix(Int_t det)          { return (TGeoHMatrix *) fMatrixArray->At(det);           } 
  TGeoHMatrix *    GetCorrectionMatrix(Int_t det){ return (TGeoHMatrix *) fMatrixCorrectionArray->At(det); } 
 
  void SetParam(IlcDCHParam *param){fDCHParam=param;}
  Double_t GetBaseExag(Int_t superlayer){return fBaseExag[superlayer];}
 

/*  
  static  Int_t    Nsect()   { return fgkNsect; }; 
  static  Int_t    Nplan()   { return fgkNplan; }; 
  static  Int_t    Ncham()   { return fgkNcham; }; 
  static  Int_t    Ndet()    { return fgkNdet;  }; 
  static  Int_t    Ntub(Int_t isp=0)    { return fgkNtub[isp];  };  //  
 
  static  Float_t  Rmin()    { return fgkRmin;  }; 
  static  Float_t  Rmax()    { return fgkRmax;  }; 
  static  Float_t  Zmax1()   { return fgkZmax1; }; 
  static  Float_t  Zmax2()   { return fgkZmax2; }; 
  static  Float_t  Step()   { return fgkstep; }; 
 
  static  Float_t  Cwidcha() { return (fgkSwidth2 - fgkSwidth1)  
                             / fgkSheight[0] * (fgkCH[0]+fgkVspace[0]/2); };
 
  static  Float_t  Cheight(Int_t isp=0) { return fgkCH[isp]+fgkVspace[isp]/2;      };  
  static  Float_t  CTheight(Int_t isp=0) { return fgkCH[isp];};
  static  Float_t  Cspace(Int_t isp=0)  { return fgkVspace[isp];  }; 
  static  Float_t  CraHght() { return fgkCamH;    }; 
  static  Float_t  CdrHght() { return 1;    }; 
  static  Float_t  CamHght() { return fgkCamH;    }; 
  static  Float_t  CroHght() { return 1;    }; 
  static  Float_t  CroWid()  { return 1;    }; 
  static  Float_t  MyThick() { return 1; }; 
  static  Float_t  DrThick() { return 1; }; 
  static  Float_t  AmThick() { return fgkAmThick; }; 
  static  Float_t  DrZpos()  { return 0;  }; 
  static  Float_t  RpadW()   { return 1;   }; 
  static  Float_t  CpadW()   { return 1;   }; 
 
  static Float_t Blen() {return fgkSlenTR2;}; 
  static Float_t TRlen() {return fgkSlenTR2;}; 
 
  static  Float_t  Sheight(Int_t isp=0)  { return fgkSheight[isp];  }; 
 */ 
  void     SetSMstatus(Int_t sm, Char_t status)     {  
            sm += 5; if (sm > kNsect-1) sm -= kNsect; fSMstatus[sm] = status; }; 
 
  virtual Bool_t   IsHole(Int_t /*iplan*/, Int_t /*icham*/, Int_t /*isect*/) const { return kFALSE; }; 
//   static  Int_t    GetDetectorSec(Int_t p, Int_t c); 
//   static  Int_t    GetDetector(Int_t p, Int_t c, Int_t s); 
//   virtual Int_t    GetPlane(Int_t d)   const; 
//   virtual Int_t    GetChamber(Int_t d) const; 
//   virtual Int_t    GetSector(Int_t d)  const; 
 
  Char_t   GetSMstatus(Int_t sm) const  { sm += 5; if (sm > 17) sm -= 18;  return fSMstatus[sm];  }; 
          Float_t  GetChamberWidth(Int_t p,Int_t s=0) const           { return fCWidth[p][s];}

          Float_t  GetChamberLength(Int_t p, Int_t c) const { return fClength[p][c]; };  
 
  //  virtual void     GetGlobal(const IlcRecPoint* , TVector3& , TMatrixF& ) const { };  
  // virtual void     GetGlobal(const IlcRecPoint* , TVector3& ) const { }; 
  
//   static  Double_t GetAlpha()  { return 2.0 * 3.14159265358979323846 / fgkNsect; };  
//  
//  static  Double_t GetAlpha1()  { return GetAlpha()-GetAlpha2(); }; 
//  static  Double_t GetAlpha2()  { return TMath::ATan(116./400.); }; 
 
 
 
//  IlcDCHgeometry* GetGeometry(IlcRunLoader* runLoader = nullptr); 
   
//   static Float_t   GetTime0(Int_t p)                        { return fgkTime0[p];    }; 
 
 
// virtual Float_t  GetRowPadSize(Int_t p, Int_t c, Int_t s) const { return fRowPadSize[p][c][s]; }; 
//  virtual Float_t  GetColPadSize(Int_t p)                   const { return fColPadSize[p];       }; 
//  virtual Float_t  GetTimeBinSize()                         const { return fTimeBinSize;         }; 
 
Float_t  GetDriftVelocity()                       const { return fDriftVelocity; }; 
 
//  double GetXLayer(int layer,int plane,int sec);
//  double GetXSector(int sec);
 
 
 private: 
 
 
 protected: 
  
   Double_t fBaseExag[24];
  IlcDCHParam         *fDCHParam;           //  pointer to DCH parameters 

/*    
  static const Int_t   fgkNsect;                            // Number of sectors in the full detector (6) 
  static const Int_t   fgkNplan;                            // Number of planes (5) 
  static const Int_t   fgkNcham;                            // Number of chambers for sector (4) 
  static const Int_t   fgkNdet;                             // Total number of detectors (3) 
  static const Int_t fgkNtub[2];                           // Tubes number in y direction (trap,box) 
  static const Float_t fgkRmin;                             // Minimal radius of the DCH 
  static const Float_t fgkRmax;                             // Maximal radius of the DCH 
  static const Float_t fgkCrack; 
 
  static const Float_t fgkZmax1;                            // Half-length of the DCH at outer radius 
  static const Float_t fgkZmax2;                            // Half-length of the DCH at inner radius 
 
 
 
 
  //BTR1 
  static const Float_t fgkSheight[2];                          // Height of the DCH-volume  (BTR1) 
  static const Float_t fgkBSheight; 
 
  static const Float_t fgkSwidth1;                          // Lower width of the DCH-volume (BTR1) 
  static const Float_t fgkSlenTR1;                          // Length of the DCH-volume (BTR1) 
   static const Float_t fgkSwidth2;                          // Upper width of the DCH-volume in spaceframe (BTR1-3) 
   
  static const Float_t fgkSlenTR2;                          // Length of mother wheel 
  static const Float_t fgkstep; 
  static const Float_t fgkTradius;                          // radius of tube 
  //BTR2 (TR sector) 
  static const Float_t fgkTRwidth1; 
 static const Float_t fgkTRwidth2; 
 static const Float_t fgkTRheight; 
 //Box1 (B sector) 
 static const Float_t fgkBwidth1; 
 static const Float_t fgkBwidth2; 
 static const Float_t fgkBheight; 
  //  
 // The sense wire(W) 
 
static  const Float_t fgkCwdrRmin; 
static  const Float_t fgkCwdrRmax; 
static  const Float_t fgkCwdrH;   
 
 static const Float_t fgkCamRmin; 
 static const Float_t fgkCamRmax; 
//CBox%d chamber 
 
static const Float_t fgkCBoxwidth; 
static const Float_t fgkCBoxheight; 
static const Float_t fgklenCBox; 
 
 static const Float_t fgkCalT; 
 static const Float_t fgkAmThick; 
 static const Float_t fgkAmZpos; 
  static const Float_t fgkCamH;                             // Height of the amplification    
static const Float_t fgkCH[2];                               // Total height of the chambers 
 
  static const Float_t fgkVspace[2];                           // Vertical spacing of the chambers 
  static const Float_t fgkHspace;                           // Horizontal spacing of the chambers 
 
  static const Float_t fgkCroW;                             // Additional width of the readout chamber frames 
 
  static const Float_t fgkCpadW;                            // Difference of outer chamber width and pad plane width 
  static const Float_t fgkRpadW;                            // Difference of outer chamber width and pad plane width 
 */ 
  Char_t               fSMstatus[kNsect];                   // Super module status byte 
 
  Float_t              fCWidth[kNplan][kNsect];                     // Outer widths of the chambers

  Float_t              fClength[kNplan][kNcham];            // Outer lengths of the chambers 
 
  Float_t              fRotA11[kNsect];                     // Matrix elements for the rotation 
  Float_t              fRotA12[kNsect];                     // Matrix elements for the rotation 
  Float_t              fRotA21[kNsect];                     // Matrix elements for the rotation 
  Float_t              fRotA22[kNsect];                     // Matrix elements for the rotation 
 
  Float_t              fRotB11[kNsect];                     // Matrix elements for the backward rotation 
  Float_t              fRotB12[kNsect];                     // Matrix elements for the backward rotation 
  Float_t              fRotB21[kNsect];                     // Matrix elements for the backward rotation 
  Float_t              fRotB22[kNsect];                     // Matrix elements for the backward rotation 
 
//   static const Double_t fgkTime0Base;                       // Base value for calculation of Time-position of pad 0 
//   static const Float_t  fgkTime0[kNplan];                   // Time-position of pad 0 
   
  Float_t              fChamberUAorig[3*kNdets][3];         // Volumes origin in 
  Float_t              fChamberUDorig[3*kNdets][3];         // the chamber 
  Float_t              fChamberUForig[3*kNdets][3];         // [3] = x, y, z 
  Float_t              fChamberUUorig[3*kNdets][3];         // 
 
  Float_t              fChamberUAboxd[3*kNdets][3];         // Volumes box 
  Float_t              fChamberUDboxd[3*kNdets][3];         // dimensions (half) 
  Float_t              fChamberUFboxd[3*kNdets][3];         // [3] = x, y, z 
  Float_t              fChamberUUboxd[3*kNdets][3];         //  
 
  TObjArray *          fMatrixArray;                        //! array of matrix - Transformation Global to Local 
  TObjArray *          fMatrixCorrectionArray;              //! array of Matrix - Transformation Cluster to  Tracking systerm 
  TObjArray *          fMatrixGeo;                          //! geo matrices 
 
 
//Float_t              fRowPadSize[kNplan][kNcham][kNsect];   // Pad size in z-direction 
//  Float_t              fColPadSize[kNplan];                 // Pad size in rphi-direction 
//  Float_t              fTimeBinSize;                        // Size of the time buckets 
 
Float_t              fDriftVelocity;                        //  Drift velocity (cm / mus) 
 
 
  ClassDef(IlcDCHgeometry,2)                                // DCH geometry class 
 
}; 
 
#endif 
 
 
 
 
 
 
 
 
 
 

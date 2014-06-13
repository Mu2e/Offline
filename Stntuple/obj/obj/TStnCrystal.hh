///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TStnCrystal_hh
#define TStnCrystal_hh

#include "Gtypes.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"
#include "TVector3.h"

#include "Stntuple/base/THexIndex.hh"
#include "Stntuple/base/THexagon.hh"
#include "Stntuple/obj/TCalHitData.hh"

// #ifndef __CINT__
// #include "CalorimeterGeom/inc/Crystal.hh"
// #else
// namespace mu2e {
//   class Crystal;
// };
// #endif

class TDisk;

class TStnCrystal: public TObject {
public:
  
protected:
  int                   fIndex;
  THexIndex             fHexIndex;

  TVector3              fCenter;
  double                fHexSize;

  TObjArray*            fListOfHits;    // keeps [not-owned] pointers to TCalDataBlock
					// display in XY view
  THexagon              fHexagon;

  int                   fFillStyle;
  int                   fFillColor;
  int                   fLineColor;
  int                   fNHits;		// number of (resolved) crystal hits

  float                 fEnergy;  	// total deposited energy
  TDisk*                fDisk;          // backward pointer to the disk

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TStnCrystal() {}
  TStnCrystal(THexIndex* Index, double X0, double Y0, double Z0, double HexSize); 

  virtual ~TStnCrystal();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TObjArray*     ListOfHits   () { return fListOfHits; }
  TCalHitData*   CalHitData(int I) { return (TCalHitData*) fListOfHits->UncheckedAt(I); }
  TVector3*      Center       () { return &fCenter   ; }
  int            Index        () { return fIndex     ; }
  int            NHits        () { return fNHits; }
  float          Energy       () { return fEnergy;     }
  double         Radius       () { return fHexagon.Radius(); }
  //  const mu2e::Crystal* Crystal() { return fCrystal;    }
  TDisk*         Disk         () { return fDisk; }

  THexIndex*     HexIndex     () { return &fHexIndex; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void  SetFillStyle(int Style) { fHexagon.fFillStyle = Style; }
  void  SetFillColor(int Color) { fHexagon.fFillColor = Color; }
  void  SetLineColor(int Color) { fHexagon.fLineColor = Color; }
  void  SetDisk     (TDisk* Disk) { fDisk = Disk; }

  void   AddHit(TCalHitData* CalHit) { 
    fListOfHits->Add(CalHit); 
    fEnergy += CalHit->Energy();
    fNHits++;
  }

  //  virtual void  Draw    (Option_t* option = "");

  virtual void  Paint   (Option_t* option = "");
  virtual void  PaintXY (Option_t* option = "");
  virtual void  PaintCal(Option_t* option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void   Clear(const char* Opt = "") ;
  virtual void   Print(const char* Opt = "") const ; // **MENU**

  ClassDef(TStnCrystal,0)
};


#endif

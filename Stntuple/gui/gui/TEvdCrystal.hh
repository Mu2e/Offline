///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdCrystal_hh
#define TEvdCrystal_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"

#include "Stntuple/base/TVisNode.hh"
#include "Stntuple/base/THexagon.hh"

#ifndef __CINT__
#include "CalorimeterGeom/inc/Crystal.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#else
namespace mu2e {
  class CaloCluster;
  class CaloCrystalHit;
  class CaloHit;
  class Crystal;
};
#endif

class TDisk;

class TEvdCrystal: public TObject {
public:
  
protected:
  const mu2e::Crystal*  fCrystal;
  mu2e::CaloCrystalHit* fCrystalHit;

  TClonesArray*         fListOfHits;
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
  TEvdCrystal() {}
  TEvdCrystal(const mu2e::Crystal* Cr, TDisk* Disk); 

  virtual ~TEvdCrystal();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TClonesArray*   ListOfHits   () { return fListOfHits; }
  int             NHits        () const { return fNHits;      }
  float           Energy       () const { return fEnergy;     }
  double          Radius       () const { return fHexagon.Radius(); }
  TDisk*          Disk         () const { return fDisk; }
  double          X0           () const { return fHexagon.X0(); }
  double          Y0           () const { return fHexagon.Y0(); }
  const THexagon* Hexagon      () const { return (const THexagon*) &fHexagon; }

  const mu2e::Crystal* Crystal() const { return fCrystal;    }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void  SetFillStyle(int Style) { fHexagon.fFillStyle = Style; }
  void  SetFillColor(int Color) { fHexagon.fFillColor = Color; }
  void  SetLineColor(int Color) { fHexagon.fLineColor = Color; }
  void  SetLineWidth(int Width) { fHexagon.fLineWidth = Width; }

  void   AddHit(const mu2e::CaloCrystalHit* CrystalHit);

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

  ClassDef(TEvdCrystal,0)
};


#endif

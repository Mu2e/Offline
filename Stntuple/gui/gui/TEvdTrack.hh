///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdTrack_hh
#define TEvdTrack_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"
#include "TEllipse.h"

class KalRep;
class TEvdTrkStrawHit;

class TEvdTrack: public TObject {
public:
  
protected:
  int                fNumber;
  const KalRep*      fKrep;

  TObjArray*         fListOfHits;

  static TEllipse*   fgEllipse;
public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdTrack();

  TEvdTrack(int Number, const KalRep* Krep);

  virtual ~TEvdTrack();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int NHits() { return fListOfHits->GetEntriesFast(); }

  TEvdTrkStrawHit* Hit(int I) { 
    return (TEvdTrkStrawHit*) fListOfHits->UncheckedAt(I); 
  }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void AddHit(TObject* Hit) { fListOfHits->Add(Hit); }
//-----------------------------------------------------------------------------
// drawing functions
//-----------------------------------------------------------------------------
  virtual void  PaintXY  (Option_t* option = "");
  virtual void  PaintRZ  (Option_t* option = "");
  virtual void  PaintCal (Option_t* option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive   (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveCal(Int_t px, Int_t py);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void  Paint(Option_t* Opt = "");
  virtual void  Clear(Option_t* Opt = "");
  virtual void  Print(Option_t* Opt = "") const ; // **MENU**

  ClassDef(TEvdTrack,0)
};


#endif

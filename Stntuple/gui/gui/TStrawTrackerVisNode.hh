// vis node displays one wedge

#ifndef TStrawTrackerVisNode_hh
#define TStrawTrackerVisNode_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"

#include "Stntuple/base/TVisNode.hh"
#include "Stntuple/gui/TEvdStraw.hh"

#ifndef __CINT__


#else

namespace mu2e {
};

#endif

class TEvdCalSection;
class TEvdCrystal;


class TStrawTrackerVisNode: public TVisNode {
public:
  enum {
    kPickClusters = 0,
    kPickTracks   = 1,
    kPickStrips   = 2
  };
protected:

  const mu2e::CalTimePeakCollection**          fCalTimePeakColl;  //

  TObjArray**        fListOfTracks;

  TEvdStrawTracker*  fStrawTracker; 	// 

  TObjArray*         fListOfStrawHits;

  const mu2e::CalTimePeak*  fTimePeak;
  
  Int_t              fDisplayHits;
  Int_t              fPickMode;

  int                fNCrystals;
  int                fFirst;

public:
					// ****** constructors and destructor

  TStrawTrackerVisNode() {}
  TStrawTrackerVisNode(const char* Name, const mu2e::Disk* Disk, int SectionID); 

  virtual ~TStrawTrackerVisNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  void  SetListOfTracks   (TObjArray** List) { fListOfTracks    = List; }

  void  SetPickMode   (Int_t Mode) { fPickMode    = Mode; }
  void  SetDisplayHits(Int_t Mode) { fDisplayHits = Mode; }

  void  Set(Int_t Side, Int_t Wedge) ; // **MENU**


  int   InitEvent();

  // ::Draw functions are called only for T**View objects

  //  virtual void  Draw    (Option_t* option = "");
  virtual void  Paint   (Option_t* option = "");
          void  PaintXY (Option_t* option = "");
          void  PaintRZ (Option_t* option = "");
          void  PaintCal(Option_t* option = "");
  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void   Clear(Option_t* Opt = "");
  virtual void   Print(Option_t* Opt = "") const ; // **MENU**

  ClassDef(TStrawTrackerVisNode,0)
};


#endif

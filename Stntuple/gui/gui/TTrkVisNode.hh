#ifndef TTrkVisNode_hh
#define TTrkVisNode_hh

#include "Gtypes.h"
#include "TClonesArray.h"

#include "Stntuple/base/TVisNode.hh"

#ifndef __CINT__
#include "KalmanTests/inc/KalRepCollection.hh"
#else
namespace mu2e {
  class KalRepCollection;
};
#endif

class TStnTrackBlock;

class TTrkVisNode: public TVisNode {
protected:

  TStnTrackBlock* fTrackBlock;
  TObjArray*      fListOfTracks;
  Color_t         fTrackColor;

  const mu2e::KalRepCollection* fKalRepCollection;

public:
					// ****** constructors and destructor

  TTrkVisNode(const char* Name = "", TStnTrackBlock* TrackBlock = 0);
  virtual ~TTrkVisNode();
					// ****** accessors

  TObjArray* GetListOfTracks() { return fListOfTracks; }
  Color_t    GetTrackColor  () { return fTrackColor;   }

					// ****** modifiers

  void  SetListOfTracks(TObjArray* List) { fListOfTracks = List; }
  void  SetTrackColor  (Color_t      c ) { fTrackColor   = c   ; }

//-----------------------------------------------------------------------------
// overloaded methods of TVisNode
//-----------------------------------------------------------------------------
  virtual int   InitEvent();
  virtual void  PaintXY(Option_t* option = "");
  virtual void  PaintRZ(Option_t* option = "");
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void  Paint  (Option_t* option = "");

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);

  ClassDef(TTrkVisNode,0)
};


#endif

#ifndef TCalView_hh
#define TCalView_hh


#include "TNamed.h"
#include "TPad.h"

class TCalView: public TNamed {
protected:
  Int_t               fPx1;
  Int_t               fPy1;
  Int_t               fPx2;
  Int_t               fPy2;
  Int_t               fSectionToDisplay;	// a disk or vane number
  TVirtualPad*        fPad;
public:

  TCalView(int Section = 0);
  virtual ~TCalView();

  TVirtualPad*  GetPad() { return fPad; }

  int    SectionToDisplay() { return fSectionToDisplay; }
  void   SetSectionToDisplay(int Section) { fSectionToDisplay = Section; }
  void   SetPad (TVirtualPad* Pad) { fPad = Pad; }
//-----------------------------------------------------------------------------
// menu
//-----------------------------------------------------------------------------
  void   SetMinClusterEnergy(float MinEnergy);  // *MENU*
  void   SetMinCrystalEnergy(float MinEnergy);  // *MENU*
  void   PrintClosestCrystal();                 // *MENU*
//-----------------------------------------------------------------------------
// overloaded virtual functions of TObject
//-----------------------------------------------------------------------------
  virtual void  Paint(Option_t* Option = "");
  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitive(Int_t px, Int_t py);

  ClassDef(TCalView,0)
};

#endif

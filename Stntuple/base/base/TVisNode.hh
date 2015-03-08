#ifndef TVisVode_hh
#define TVisVode_hh

#include "TObject.h"
#include "TString.h"

class TVisNode: public TObject {
protected:
  TString    fName;
  TObject*   fClosestObject;
  int        fSectionToDisplay;
  int        fDebugLevel;
public:
					// ****** constructors and destructor
  TVisNode(const char* name = "");
  virtual ~TVisNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TObject*  GetClosestObject() { return fClosestObject; }
  
  virtual const char* GetName() const { return fName.Data(); }

  int   SectionToDisplay() { return fSectionToDisplay; }

  int   DebugLevel() { return fDebugLevel; }

					// called by TEvdManager::DisplayEvent
  virtual int   InitEvent() = 0;

  virtual void  PaintXY(Option_t* option = "") = 0;
  virtual void  PaintRZ(Option_t* option = "") = 0;
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py) = 0;

  void SetSectionToDisplay(int Section) { fSectionToDisplay= Section; }
  void SetDebugLevel(int Level) { fDebugLevel = Level; }

  ClassDef(TVisNode,0)
};

#endif

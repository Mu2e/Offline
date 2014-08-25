///////////////////////////////////////////////////////////////////////////////
// vis node - one disk
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdCalSection_hh
#define TEvdCalSection_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"

#ifndef __CINT__

#include "CalorimeterGeom/inc/Disk.hh"

#else

namespace mu2e {
  class  Disk;
};

#endif

class TEvdCalSection: public TObject {
public:
  
protected:
  const mu2e::Disk*   fDisk;
  int                 fSectionID;

  TEllipse*           fEllipse[2];

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdCalSection() {}
  TEvdCalSection(const mu2e::Disk* Disk, int SectionID);

  virtual ~TEvdCalSection();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  const mu2e::Disk*   Disk     () { return fDisk; }
  int                 SectionID() { return fSectionID; }

  int NCrystals() { return fDisk->nCrystals(); }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------

  //  virtual void  Draw    (Option_t* option = "");

  virtual void  Paint   (Option_t* option = "");
  virtual void  PaintXY (Option_t* Option = "");
  virtual void  PaintRZ (Option_t* Option = "");
  virtual void  PaintCal(Option_t* Option = "");

  int   InitEvent();

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);

  //  virtual void   Print(const char* Opt = "") const ; // **MENU**

  ClassDef(TEvdCalSection,0)
};


#endif

///////////////////////////////////////////////////////////////////////////////
// Vis node - One section of CRV
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdCalSection_hh
#define TEvdCalSection_hh

#include "Rtypes.h" // (?)  - Gtypes is obsolete and contains nothing
#include "TClonesArray.h"
#include "TBox.h"

#ifndef __CINT__

#include "CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"

#else

namespace mu2e 
{
	class  CRSScintillatorShield;
};

#endif

class TEvdCrvSection : public TObject 
{
public:
	//-----------------------------------------------------------------------------
	// constructors and destructor
	//-----------------------------------------------------------------------------
	TEvdCrvSection() {}
	TEvdCrvSection(/*const mu2e::CRSScintillatorShield* Shield,*/ int SectionID);

	virtual ~TEvdCrvSection();
	//-----------------------------------------------------------------------------
	// accessors
	//-----------------------------------------------------------------------------
	//const mu2e::CRSScintillatorShield*	Shield()	{ return fShield; }
	int									SectionID()	{ return fSectionID; }

	//int NBars() { return fShield->(); } - No current support for getting the number of bars in a given section easily
	//-----------------------------------------------------------------------------
	// modifiers
	//-----------------------------------------------------------------------------

	//  virtual void  Draw    (Option_t* option = "");

	virtual void  Paint(Option_t* option = "");
	//virtual void  PaintXY(Option_t* Option = "");
	//virtual void  PaintRZ(Option_t* Option = "");
	//virtual void  PaintCal(Option_t* Option = "");
	virtual void  PaintCrv(Option_t* Option = "");

	int   InitEvent() { return 0; }

	//  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

	virtual Int_t DistancetoPrimitive(Int_t px, Int_t py);
	virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
	virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);

	//  virtual void   Print(const char* Opt = "") const ; // **MENU**

protected:
	//const mu2e::CRSScintillatorShield*  fShield;
	int									fSectionID;
	TBox*								fBox;


	ClassDef(TEvdCrvSection, 0)
};


#endif
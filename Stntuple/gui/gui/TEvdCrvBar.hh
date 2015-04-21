///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdCrvBar_hh
#define TEvdCrvBar_hh

#include <vector>
#include <algorithm>

#include "Rtypes.h" // (?)  - Gtypes is obsolete and contains nothing
#include "TClonesArray.h"
#include "TBox.h"

#include "Stntuple/base/TVisNode.hh"

#ifndef __CINT__

#include "CosmicRayShieldGeom/inc/CRSScintillatorBar.hh"
//#include "CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"
//#include "RecoDataProducts/inc/CrvRecoPulses.hh"

#include "RecoDataProducts/inc/CrvRecoPulsesCollection.hh" //


#else
namespace mu2e 
{
	class CRSScintillatorBar;
	class CrvRecoPulses;
        struct CrvRecoPulses::CrvSingleRecoPulse;
};
#endif

//class CRSScintillatorShield;  // (?)

class TEvdCrvBar : public TObject 
{
public:
	//-----------------------------------------------------------------------------
	// Constructors and destructor
	//-----------------------------------------------------------------------------
	TEvdCrvBar() {}
	TEvdCrvBar(const mu2e::CRSScintillatorBar& Bar, int SectionID/*CRSScintillatorShield* Shield*/);

	virtual ~TEvdCrvBar();
	//-----------------------------------------------------------------------------
	// Accessors
	//-----------------------------------------------------------------------------
	//int					NHits(int SiPM)	const	{ return fNHits[SiPM]; }
	//float					Height()		const	{ return fPulseHeight; }
	//std::vector<int>		PE(int SiPM)	const	{ return sipmPulses[SiPM]._PEs; }
	//int					maxPE(int SiPM);

	double		X0()	const	{ return fBox->GetBBoxCenter().GetX(); } //X Center of box?
	double		Y0()	const	{ return fBox->GetBBoxCenter().GetY(); } //Y Center of box?
	const TBox*	Box()	const	{ return fBox; }
	int			Bar()	const	{ return barIndex; }

	//const mu2e::CRSScintillatorBar Bar() const { return fBar; }
	//-----------------------------------------------------------------------------
	// Modifiers
	//-----------------------------------------------------------------------------
	
  void   AddPulse(const mu2e::CrvRecoPulses::CrvSingleRecoPulse &BarPulse, int SiPM);
	
	//void  SetFillStyle(int Style) { fBox->fFillStyle(Style); }
	//void  SetFillColor(int Color) { fBox->fFillColor(Color); }
	void	SetFillColor(int Color, int SiPM) { sipm[SiPM]->SetFillColor(Color); }
	void	SetFillStyle(int Style, int SiPM) { sipm[SiPM]->SetFillStyle(Style); }

	//  virtual void  Draw    (Option_t* option = "");
	//  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

	virtual void  Paint(Option_t* option = "");
	virtual void  PaintXY(Option_t* option = "");
	virtual void  PaintCrv(Option_t* option = "");

	virtual Int_t DistancetoPrimitive(Int_t px, Int_t py);
	virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
	virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);
	//-----------------------------------------------------------------------------
	// Overloaded methods of TObject
	//-----------------------------------------------------------------------------
	virtual void   Clear(const char* Opt = "");
	virtual void   Print(const char* Opt = "") const; // **MENU**

protected:
	//const mu2e::CRSScintillatorBar	fBar;
	
	int barIndex; //Index of the bar
	//int fNHits[4]; // Number of pulses for each SiPM
  std::vector<const mu2e::CrvRecoPulses::CrvSingleRecoPulse*> sipmPulses[4];	//Vector of pulses for each SiPM  // Read as a vector of pointers to constant objects
	
	
	// display in XY view
	TBox*					fBox;
	TBox*					sipm[4];

	//int						fFillStyle;
	//int						fFillColor;	

	//std::vector<float>	fPulseHeight[4];	// Height of pulse in volts ( ? ? Unnecessary ? ? )
	//std::vector<int>		sipmPEs[4];		// Height of pulse in PEs		( ? ? Unnecessary ? ? )


	ClassDef(TEvdCrvBar, 0)
};


#endif

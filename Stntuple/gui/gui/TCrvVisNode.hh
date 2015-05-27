#ifndef TCrvVisNode_hh
#define TCrvVisNode_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"	
#include "TPad.h"
#include "TBox.h"

#include "Stntuple/base/TVisNode.hh"
#include "Stntuple/gui/TEvdCrvBar.hh"


#ifndef __CINT__

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorBar.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorLayer.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorModule.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"
#include "RecoDataProducts/inc/CrvRecoPulsesCollection.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#else

namespace mu2e
{
	class CosmicRayShield;
	class CRSScintillatorBar;
	class CRSScintillatorLayer;
	class CRSScintillatorModule;
	class CRSScintillatorShield;
	class CrvRecoPulsesCollection;
	class CRSScintillatorBarIndex;
	class CrvRecoPulses;
	struct CrvRecoPulses::CrvSingleRecoPulse;
};

#endif

class TCrvVisNode : public TVisNode
{
public:
	//-----------------------------------------------------------------------------
	// Constructors and Destructor
	//-----------------------------------------------------------------------------
	TCrvVisNode() {}
	TCrvVisNode(const char* Name, /*const mu2e::CosmicRayShield* CRV,*/ int SectionID);

	virtual ~TCrvVisNode();

	//-----------------------------------------------------------------------------
	// Accessors
	//-----------------------------------------------------------------------------
	TObjArray* GetListOfTracks();

	int		SectionID() { return fSectionID; }
	TEvdCrvBar*  EvdCrvBar(int barIndex);

	//-----------------------------------------------------------------------------
	// Modifiers
	//-----------------------------------------------------------------------------
	void	SetRecoPulsesCollection(mu2e::CrvRecoPulsesCollection** List) { fCrvRecoPulsesCollection = List; }

	//void	SetMinPulseHeight(float PulseHeight) { fMinPulseHeight = PulseHeight; }
	void	SetMinPulsePEs(float minPulsePEs){ fMinPulsePEs = minPulsePEs; }
	void	SetTimeWindow(float timeLow, float timeHigh);

	//-----------------------------------------------------------------------------
	// Overloaded methods of TVisNode
	//-----------------------------------------------------------------------------
	virtual int		InitEvent();
	void			UpdateEvent();
	virtual void	PaintXY(Option_t* option = "");
	virtual void	PaintCrv(Option_t* option = "");
	virtual void	PaintRZ(Option_t* option = "");

	//-----------------------------------------------------------------------------
	// Overloaded methods of TObject
	//-----------------------------------------------------------------------------
	virtual void	Paint(Option_t* option = "");
	virtual void	Clear(Option_t* Opt = "");
	virtual void	Print(Option_t* Opt = "") const; // **MENU**

	virtual Int_t	DistancetoPrimitive(Int_t px, Int_t py);
	virtual Int_t	DistancetoPrimitiveXY(Int_t px, Int_t py);
	virtual Int_t	DistancetoPrimitiveRZ(Int_t px, Int_t py);

protected:
	mu2e::CrvRecoPulsesCollection**	fCrvRecoPulsesCollection;

	int				fSectionID;
	
	TObjArray*		fListOfEvdCrvBars;

	//double	fMinPulseHeight;
	float		fMinPulsePEs;
	float		ftimeLow;
	float		ftimeHigh;

	Int_t colorPalette[1000];

	int getCRVSection(int shieldNumber);


	ClassDef(TCrvVisNode, 0)
};


#endif

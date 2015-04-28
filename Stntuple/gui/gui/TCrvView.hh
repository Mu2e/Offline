#ifndef TCrvView_hh
#define TCrvView_hh


#include "TNamed.h"
#include "TPad.h"

class TCrvView : public TNamed {
protected:
	Int_t			fPx1;
	Int_t			fPy1;
	Int_t			fPx2;
	Int_t			fPy2;
	Int_t			fSectionToDisplay;
	TVirtualPad*	fPad;		 
public:
	TCrvView(){ fSectionToDisplay = 0; fPad = 0; };
	TCrvView(int Section);
	virtual ~TCrvView();

	TVirtualPad* GetPad() { return fPad; }

	int		SectionToDisplay()					{ return fSectionToDisplay; }
	void	SetSectionToDisplay(int Section)	{ fSectionToDisplay = Section; }
	void	SetPad(TVirtualPad* Pad)			{ fPad = Pad; }
	void	SetTimeWindow(float lowWindow, float highWindow);

//-----------------------------------------------------------------------------
// menu
//-----------------------------------------------------------------------------
	//void	SetMinPulseHeight(float MinHeight);	// *MENU*
	void	SetMinPulsePEs(float MinPEs);		// *MENU*
	void	PrintClosestBar();					// *MENU*
//-----------------------------------------------------------------------------
// overloaded virtual functions of TObject
//-----------------------------------------------------------------------------
	virtual void	Paint(Option_t* Option = "");
	virtual void	ExecuteEvent(Int_t event, Int_t px, Int_t py);
	virtual Int_t	DistancetoPrimitive(Int_t px, Int_t py);

	ClassDef(TCrvView, 0)
};

#endif
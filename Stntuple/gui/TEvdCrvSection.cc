#include "TVirtualX.h"
#include "TPad.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TLine.h"
#include "TArc.h"
#include "TArrow.h"
#include "TMath.h"
#include "TBox.h"
#include "TObjArray.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "Stntuple/gui/TEvdCluster.hh"
#include "Stntuple/gui/TEvdCrystal.hh"
#include "Stntuple/gui/TEvdCrvSection.hh"
#include "Stntuple/gui/TStnVisManager.hh"

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"

ClassImp(TEvdCrvSection)

//_____________________________________________________________________________
TEvdCrvSection::TEvdCrvSection(/*const mu2e::CRSScintillatorShield* Shield,*/ int SectionID) : TObject()
{
	//fShield = Shield;
	fSectionID = SectionID;

	//GeomHandle<CosmicRayShield> CRS;

	//switch (fSectionID)
	//{
	//case 0:
	//	//Right
	//	fBox = new TBox(); //Make a box to outline this CRV Section
	//	break;

	//case 1:
	//	//Left
	//	fBox = new TBox(); //Make a box to outline this CRV Section
	//	break;

	//case 2:
	//	//Top DS
	//	fBox = new TBox(); //Make a box to outline this CRV Section
	//	break;

	//case 3:
	//	//Downstream
	//	fBox = new TBox(); //Make a box to outline this CRV Section
	//	break;

	//case 4:
	//	//Upstream
	//	fBox = new TBox(); //Make a box to outline this CRV Section
	//	break;

	//case 5:
	//	//Cryo Hole Upstream
	//	break;
	//case 6:
	//	//Cryo Hole Downstream
	//	break;
	//case 7:
	//	//Cryo Hole Top
	//	break;
	//case 8:
	//	//Top TS
	//	fBox = new TBox(); //Make a box to outline this CRV Section
	//	break;
	//}

	fBox = new TBox(0, 0, 100, 100);

	fBox->SetLineColor(kBlack);
	fBox->SetFillColor(kBlue);
	fBox->SetFillStyle(3001);
}

//_____________________________________________________________________________
TEvdCrvSection::~TEvdCrvSection()
{
}

//_____________________________________________________________________________
void TEvdCrvSection::Paint(Option_t* option)
{
	// Paints one CRV section

	int iv;

	const char* view = TVisManager::Instance()->GetCurrentView();

	if (strstr(view, "crv") != 0)
	{
		sscanf(view, "crv,%i", &iv);
		if (iv == fSectionID)
			PaintCrv(option);
	}
}

void TEvdCrvSection::PaintCrv(Option_t* Option) {
	// Paints one CRV Section

	//fBox->Paint(Option);
	fBox->Draw("L");
}


//_____________________________________________________________________________
Int_t TEvdCrvSection::DistancetoPrimitive(Int_t px, Int_t py) {
	return 9999;
}

//_____________________________________________________________________________
Int_t TEvdCrvSection::DistancetoPrimitiveXY(Int_t px, Int_t py) {

	Int_t dist = 9999;

	static TVector3 global;
	//   static TVector3 local;

	//  Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

	global.SetXYZ(gPad->AbsPixeltoX(px), gPad->AbsPixeltoY(py), 0);

	return dist;
}

//_____________________________________________________________________________
Int_t TEvdCrvSection::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
	return 9999;
}

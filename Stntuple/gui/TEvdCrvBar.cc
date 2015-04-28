#include "TVector3.h"
#include "TBox.h"

#include "art/Framework/Principal/Handle.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "Stntuple/gui/TEvdCrvBar.hh"
#include "Stntuple/gui/TStnVisManager.hh"
#include "Stntuple/base/TObjHandle.hh"

#include "CosmicRayShieldGeom/inc/CRSScintillatorBar.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"
#include "RecoDataProducts/inc/CrvRecoPulses.hh"

ClassImp(TEvdCrvBar)

//_____________________________________________________________________________
TEvdCrvBar::TEvdCrvBar(const mu2e::CRSScintillatorBar& Bar, int SectionID) : TObject()
{
	//fBar = Bar;
	const CLHEP::Hep3Vector *pos;
	pos = &Bar.getPosition();
	barIndex = Bar.index().asInt();
	double halfWidth = Bar.getHalfWidth();
	double halfThickness = Bar.getHalfThickness();
	
	switch (SectionID)
	{
	case 0:
		//Right
		fBox = new TBox(pos->z() - halfWidth, pos->x() - halfThickness, pos->z() + halfWidth, pos->x() + halfThickness);
		sipm[0] = new TBox(pos->z() - halfWidth, pos->x() - halfThickness, pos->z(), pos->x());
		sipm[1] = new TBox(pos->z() - halfWidth, pos->x() + halfThickness, pos->z(), pos->x());
		sipm[2] = new TBox(pos->z() + halfWidth, pos->x() - halfThickness, pos->z(), pos->x());
		sipm[3] = new TBox(pos->z() + halfWidth, pos->x() + halfThickness, pos->z(), pos->x());

		break;
	case 1:
		//Left
		fBox = new TBox(pos->z() - halfWidth, pos->x() - halfThickness, pos->z() + halfWidth, pos->x() + halfThickness);
		sipm[0] = new TBox(pos->z() - halfWidth, pos->x() - halfThickness, pos->z(), pos->x());
		sipm[1] = new TBox(pos->z() - halfWidth, pos->x() + halfThickness, pos->z(), pos->x());
		sipm[2] = new TBox(pos->z() + halfWidth, pos->x() - halfThickness, pos->z(), pos->x());
		sipm[3] = new TBox(pos->z() + halfWidth, pos->x() + halfThickness, pos->z(), pos->x());
		break;
	case 2:
		//Top DS
		fBox = new TBox(pos->z() - halfWidth, pos->y() - halfThickness, pos->z() + halfWidth, pos->y() + halfThickness);
		sipm[0] = new TBox(pos->z() - halfWidth, pos->y() - halfThickness, pos->z(), pos->y());
		sipm[1] = new TBox(pos->z() - halfWidth, pos->y() + halfThickness, pos->z(), pos->y());
		sipm[2] = new TBox(pos->z() + halfWidth, pos->y() - halfThickness, pos->z(), pos->y());
		sipm[3] = new TBox(pos->z() + halfWidth, pos->y() + halfThickness, pos->z(), pos->y());
		break;
	case 3:
		//Downstream
		fBox = new TBox(pos->y() - halfWidth, pos->z() - halfThickness, pos->y() + halfWidth, pos->z() + halfThickness);
		sipm[0] = new TBox(pos->y() - halfWidth, pos->z() - halfThickness, pos->y(), pos->z());
		sipm[1] = new TBox(pos->y() - halfWidth, pos->z() + halfThickness, pos->y() , pos->z());
		sipm[2] = new TBox(pos->y() + halfWidth, pos->z() - halfThickness, pos->y(), pos->z());
		sipm[3] = new TBox(pos->y() + halfWidth, pos->z() + halfThickness, pos->y(), pos->z());
		break;
	case 4:
		//Upstream
		fBox = new TBox(pos->y() - halfWidth, pos->z() - halfThickness, pos->y() + halfWidth, pos->z() + halfThickness);
		sipm[0] = new TBox(pos->y() - halfWidth, pos->z() - halfThickness, pos->y(), pos->z());
		sipm[1] = new TBox(pos->y() - halfWidth, pos->z() + halfThickness, pos->y(), pos->z());
		sipm[2] = new TBox(pos->y() + halfWidth, pos->z() - halfThickness, pos->y(), pos->z());
		sipm[3] = new TBox(pos->y() + halfWidth, pos->z() + halfThickness, pos->y(), pos->z());
		break;
	case 5://Cryo
	case 6://Cryo
	case 7://Cryo
		break;
	case 8:
		//TopTS (May not be needed?)
		fBox = new TBox(pos->z() - halfWidth, pos->y() - halfThickness, pos->z() + halfWidth, pos->y() + halfThickness);
		sipm[0] = new TBox(pos->z() - halfWidth, pos->y() - halfThickness, pos->z(), pos->y() + halfThickness);
		sipm[1] = new TBox(0,0,0,0); //Not drawn, single end readout
		sipm[2] = new TBox(pos->z() + halfWidth, pos->y() - halfThickness, pos->z(), pos->y() + halfThickness);
		sipm[3] = new TBox(0,0,0,0); //Not drawn, single end readout
		break;
	}

	fBox->SetLineWidth(1);
	fBox->SetLineColor(kBlack);
	fBox->SetFillColor(kWhite);
	fBox->SetFillStyle(0); //Solid fill
}

//-----------------------------------------------------------------------------
TEvdCrvBar::~TEvdCrvBar() {
}

//-----------------------------------------------------------------------------
void TEvdCrvBar::Paint(Option_t* Option) {

	const char* view = TVisManager::Instance()->GetCurrentView();

	if (strstr(view, "crv") != 0)
		PaintCrv(Option);
	
	gPad->Modified();
}

//-----------------------------------------------------------------------------
void TEvdCrvBar::PaintXY(Option_t* Option) {
}

//-----------------------------------------------------------------------------
void TEvdCrvBar::PaintCrv(Option_t* Option) {
	
	sipm[0]->Paint(); //Draw SiPMs
	sipm[1]->Paint();
	sipm[2]->Paint();
	sipm[3]->Paint();

	fBox->Paint("L"); //Draw CrvBar with outline
}


//_____________________________________________________________________________
Int_t TEvdCrvBar::DistancetoPrimitive(Int_t px, Int_t py) {
	Int_t dist = 9999;

	static TVector3 global;
	//   static TVector3 local;

	//   Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

	global.SetXYZ(gPad->AbsPixeltoX(px), gPad->AbsPixeltoY(py), 0);

	return dist;
}

//_____________________________________________________________________________
Int_t TEvdCrvBar::DistancetoPrimitiveXY(Int_t px, Int_t py) {

	Int_t dist = 9999;

	static TVector3 global;
	//   static TVector3 local;

	//   Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

	global.SetXYZ(gPad->AbsPixeltoX(px), gPad->AbsPixeltoY(py), 0);

	return dist;
}

//_____________________________________________________________________________
Int_t TEvdCrvBar::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
	return 9999;
}

//_____________________________________________________________________________
void TEvdCrvBar::AddPulse(const mu2e::CrvRecoPulses::CrvSingleRecoPulse &BarPulse, int SiPM) {
        const mu2e::CrvRecoPulses::CrvSingleRecoPulse *barPulsePtr;
	barPulsePtr = &BarPulse;
	
	sipmPulses[SiPM].push_back(barPulsePtr);
}

//_____________________________________________________________________________
void TEvdCrvBar::Clear(Option_t* Opt) {
	
	for (int i = 0; i < 4; i++)
	{
		sipm[i]->SetFillColor(kWhite); //Return to white
		sipmPulses[i].clear(); //Remove associated pulses
	}
	//fNHits = 0;
	//fEnergy = 0.;
}


//_____________________________________________________________________________
void TEvdCrvBar::Print(Option_t* Opt) const {

	//TObjHandle*                  h;
	//const mu2e::CaloCrystalHit*  hit;

	//printf(" X0 = %10.3f Y0 = %10.3f  E = %10.3f  njits = %5i\n",
	//	X0(), Y0(), fEnergy, fNHits);

	//printf("----------------------------------------------------------------\n");
	//printf("CrystalID      Time   Energy    EnergyTot  NRoids               \n");
	//printf("----------------------------------------------------------------\n");
	//for (int i = 0; i<fNHits; i++) {

	//	h = (TObjHandle*) fListOfHits->At(i);
	//	hit = (const mu2e::CaloCrystalHit*) h->Object();

	//	printf("%7i  %10.3f %10.3f %10.3f %5i\n",
	//		hit->id(),
	//		hit->time(),
	//		hit->energyDep(),
	//		hit->energyDepTotal(),
	//		hit->numberOfROIdsUsed());
	//}
}

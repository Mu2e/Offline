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
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "RecoDataProducts/inc/CrvRecoPulses.hh"

ClassImp(TEvdCrvBar)

//_____________________________________________________________________________
TEvdCrvBar::TEvdCrvBar(const mu2e::CRSScintillatorBar& Bar, int SectionID) : TObject()
{
	//fBar = Bar;
	const CLHEP::Hep3Vector *pos;
	pos = &Bar.getPosition();
	barIndex = Bar.index();
	double halfWidth = Bar.getHalfWidth();
	double halfThickness = Bar.getHalfThickness();
	fSectionID = SectionID;

	//set some defaults which will be updated later
	fthreshold = 10;
	ftimeLow = 400;
	ftimeHigh = 1695;
	
	switch (fSectionID)
	{
	case 0:
		//Right
		fBox = new TBox(pos->z() - halfWidth, pos->x() - halfThickness, pos->z() + halfWidth, pos->x() + halfThickness);
		sipm[0] = new TBox(pos->z() - halfWidth, pos->x() - halfThickness, pos->z(), pos->x());
		sipm[1] = new TBox(pos->z() - halfWidth, pos->x() + halfThickness, pos->z(), pos->x());
		sipm[2] = new TBox(pos->z() + halfWidth, pos->x() - halfThickness, pos->z(), pos->x());
		sipm[3] = new TBox(pos->z() + halfWidth, pos->x() + halfThickness, pos->z(), pos->x());
		xDispLoc = pos->z();
		yDispLoc = pos->x();
		break;
	case 1:
		//Left
		fBox = new TBox(pos->z() - halfWidth, pos->x() - halfThickness, pos->z() + halfWidth, pos->x() + halfThickness);
		sipm[0] = new TBox(pos->z() - halfWidth, pos->x() - halfThickness, pos->z(), pos->x());
		sipm[1] = new TBox(pos->z() - halfWidth, pos->x() + halfThickness, pos->z(), pos->x());
		sipm[2] = new TBox(pos->z() + halfWidth, pos->x() - halfThickness, pos->z(), pos->x());
		sipm[3] = new TBox(pos->z() + halfWidth, pos->x() + halfThickness, pos->z(), pos->x());
		xDispLoc = pos->z();
		yDispLoc = pos->x();
		break;
	case 2:
		//Top DS
		fBox = new TBox(pos->z() - halfWidth, pos->y() - halfThickness, pos->z() + halfWidth, pos->y() + halfThickness);
		sipm[0] = new TBox(pos->z() - halfWidth, pos->y() - halfThickness, pos->z(), pos->y());
		sipm[1] = new TBox(pos->z() - halfWidth, pos->y() + halfThickness, pos->z(), pos->y());
		sipm[2] = new TBox(pos->z() + halfWidth, pos->y() - halfThickness, pos->z(), pos->y());
		sipm[3] = new TBox(pos->z() + halfWidth, pos->y() + halfThickness, pos->z(), pos->y());
		xDispLoc = pos->z();
		yDispLoc = pos->y();
		break;
	case 3:
		//Downstream
		if (Bar.id().getShieldNumber() == 18)
		{//	Apply Offset in Z for overlapping sections
			fBox = new TBox(pos->y() - halfWidth, (pos->z() + 60) - halfThickness, pos->y() + halfWidth, (pos->z() + 60) + halfThickness);
			sipm[0] = new TBox(pos->y() - halfWidth, (pos->z() + 60) - halfThickness, pos->y(), (pos->z() + 60));
			sipm[1] = new TBox(pos->y() - halfWidth, (pos->z() + 60) + halfThickness, pos->y(), (pos->z() + 60));
			sipm[2] = new TBox(pos->y() + halfWidth, (pos->z() + 60) - halfThickness, pos->y(), (pos->z() + 60));
			sipm[3] = new TBox(pos->y() + halfWidth, (pos->z() + 60) + halfThickness, pos->y(), (pos->z() + 60));
			xDispLoc = pos->y();
			yDispLoc = pos->z()+60;
		}
		else if (Bar.id().getShieldNumber() == 19)
		{// Apply Offset in Z for overlapping sections
			fBox = new TBox(pos->y() - halfWidth, (pos->z() - 60) - halfThickness, pos->y() + halfWidth, (pos->z() - 60) + halfThickness);
			sipm[0] = new TBox(pos->y() - halfWidth, (pos->z() - 60) - halfThickness, pos->y(), (pos->z() - 60));
			sipm[1] = new TBox(pos->y() - halfWidth, (pos->z() - 60) + halfThickness, pos->y(), (pos->z() - 60));
			sipm[2] = new TBox(pos->y() + halfWidth, (pos->z() - 60) - halfThickness, pos->y(), (pos->z() - 60));
			sipm[3] = new TBox(pos->y() + halfWidth, (pos->z() - 60) + halfThickness, pos->y(), (pos->z() - 60));
			xDispLoc = pos->y();
			yDispLoc = pos->z()-60;
		}
		else
		{
			fBox = new TBox(pos->y() - halfWidth, pos->z() - halfThickness, pos->y() + halfWidth, pos->z() + halfThickness);
			sipm[0] = new TBox(pos->y() - halfWidth, pos->z() - halfThickness, pos->y(), pos->z());
			sipm[1] = new TBox(pos->y() - halfWidth, pos->z() + halfThickness, pos->y(), pos->z());
			sipm[2] = new TBox(pos->y() + halfWidth, pos->z() - halfThickness, pos->y(), pos->z());
			sipm[3] = new TBox(pos->y() + halfWidth, pos->z() + halfThickness, pos->y(), pos->z());
			xDispLoc = pos->y();
			yDispLoc = pos->z();
		}
		break;
	case 4:
		//Upstream
		fBox = new TBox(pos->y() - halfWidth, pos->z() - halfThickness, pos->y() + halfWidth, pos->z() + halfThickness);
		sipm[0] = new TBox(pos->y() - halfWidth, pos->z() - halfThickness, pos->y(), pos->z() + halfThickness);
		sipm[1] = new TBox(0, 0, 0, 0); //Not drawn, single end readout
		sipm[2] = new TBox(pos->y() + halfWidth, pos->z() - halfThickness, pos->y(), pos->z() + halfThickness);
		sipm[3] = new TBox(0, 0, 0, 0); //Not drawn, single end readout
		xDispLoc = pos->y();
		yDispLoc = pos->z();
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
		xDispLoc = pos->z();
		yDispLoc = pos->y();
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
	//mu2e::GeomHandle< mu2e::CosmicRayShield > CRS;
	//const mu2e::CRSScintillatorBar bar = CRS->getBar(barIndex);
	//const CLHEP::Hep3Vector *pos;
	//pos = &bar.getPosition();

	//double x1, x2, y1, y2;
	//
	//gPad->GetRange(x1, y1, x2, y2);

	////Draw bars only within the visible range of the window
	//if ((pos->x() > x1) && (pos->x() < x2) && (pos->y() > y1) && (pos->y() > y2))
	//{
		sipm[0]->Paint(); //Draw SiPMs
		sipm[1]->Paint();
		sipm[2]->Paint();
		sipm[3]->Paint();

		fBox->Paint("L"); //Draw CrvBar with outline
	//}
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

const mu2e::CrvRecoPulses::CrvSingleRecoPulse* TEvdCrvBar::lastPulseInWindow(int SiPM, float threshold, float timeLow, float timeHigh)
{
	fthreshold = threshold;
	ftimeLow = timeLow;
	ftimeHigh = timeHigh;
	
	const mu2e::CrvRecoPulses::CrvSingleRecoPulse* tempPulse = 0; //Create a pointer to a 'null' pulse

	for (unsigned int i = 0; i < sipmPulses[SiPM].size(); i++)
	{
		// !!!!! THIS ASSUMES THE PULSES ARE ORGANIZED BY TIME !!!!!
		//					And they are...

		//Skip until we get to a pulse in the time window
		if (sipmPulses[SiPM][i]->_leadingEdge < timeLow || sipmPulses[SiPM][i]->_PEs < threshold)
			continue;
		//if the pulse occurs beyond the time window then break and return the previous pulse, if not then set the buffer (lastpulse)
		if (sipmPulses[SiPM][i]->_leadingEdge > timeHigh)
			break;

		tempPulse = sipmPulses[SiPM][i];		
	}

	return tempPulse;
}


//_____________________________________________________________________________
void TEvdCrvBar::Print(Option_t* Opt) {
	const mu2e::CrvRecoPulses::CrvSingleRecoPulse* tempPulse = 0;

	printf("----------------------------------------------------------------\n");
	printf("Bar	SiPM	Section	Time	PEs	\n");
	printf("----------------------------------------------------------------\n");
	
	for (int SiPM = 0; SiPM < 4; SiPM++)
	{
		if (fSectionID == 4 && (SiPM == 1 || SiPM == 3)) // Or if it is Top TS, but I will need to figure out how to determine this
			continue;

		tempPulse = lastPulseInWindow(SiPM, fthreshold, ftimeLow, ftimeHigh);
		if (tempPulse)
		{
			printf("thresh %f , tL %f , tH %f , %p \n", fthreshold, ftimeLow, ftimeHigh, tempPulse);
			printf("%i	%i	%i	%f	%i \n", barIndex.asInt(), SiPM, fSectionID, tempPulse->_leadingEdge, tempPulse->_PEs);
		}
	}

}

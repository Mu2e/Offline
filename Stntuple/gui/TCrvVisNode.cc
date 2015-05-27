#include <cmath>

#include "TVirtualX.h"
#include "TPad.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TBox.h"
#include "TObjArray.h"
#include "TColor.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "Stntuple/gui/TCrvVisNode.hh"
#include "Stntuple/gui/TStnVisManager.hh"
#include "Stntuple/gui/TEvdCrvBar.hh"

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorModule.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorLayer.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorBar.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "RecoDataProducts/inc/CrvRecoPulsesCollection.hh"


ClassImp(TCrvVisNode)

TCrvVisNode::TCrvVisNode(const char* Name, /*const mu2e::CosmicRayShield* CRV,*/ int SectionID) : TVisNode(Name)
{
	// This gradient palette will always be the same
	// but the time-scale is governed by the vis manager slider.
	// That way, a time window can be selected, and the colors will
	// go from blue-red (low-high) in the time window
	Double_t stop [] = { 0.00, 0.20, 0.40, 0.60, 0.70, 1.00 };
	Double_t r [] = { 0.00, 0.00, 0.00, 0.97, 0.97, 0.10 };
	Double_t g [] = { 0.97, 0.30, 0.40, 0.97, 0.00, 0.00 };
	Double_t b [] = { 0.97, 0.97, 0.00, 0.00, 0.00, 0.00 };
	Int_t FI = TColor::CreateGradientColorTable(6, stop, r, g, b, 1000);
	for (int i = 0; i < 1000; i++) colorPalette[i] = FI + i;

	fMinPulsePEs = 10; // Default number for minimum PEs

	mu2e::GeomHandle< mu2e::CosmicRayShield > CRS;
	int nmodules, nlayers, nbars;

	TEvdCrvBar *evd_bar;

	fSectionID = SectionID;
	
	fListOfEvdCrvBars = new TObjArray();

	switch (fSectionID)
	{
	case 0:
		//Right
		for (int shield = 0; shield <= 5; shield++) // Loop over all the shields in right, but skip short section
		{
			nmodules = CRS->getCRSScintillatorShield(shield).nModules();
			for (int module = 0; module < nmodules; module++) // Loop over all the modules in the shield
			{
				nlayers = CRS->getCRSScintillatorShield(shield).getModule(module).nLayers();
				for (int layer = 0; layer < nlayers; layer++) // Loop over all the layers in the module
				{
					nbars = CRS->getCRSScintillatorShield(shield).getModule(module).getLayer(layer).nBars();
					for (int bar = 0; bar < nbars; bar++) // Loop over all the bars in the layer
					{
						evd_bar = new TEvdCrvBar(CRS->getCRSScintillatorShield(shield).getModule(module).getLayer(layer).getBar(bar), fSectionID);
						fListOfEvdCrvBars->Add(evd_bar);
					}
				}
			}
			if (!shield) // Skip short shield (1)
				shield++;
		}
		break;

	case 1:
		//Left
		for (int shield = 6; shield <= 8; shield++) // Loop over all the shields in left
		{
			nmodules = CRS->getCRSScintillatorShield(shield).nModules();
			for (int module = 0; module < nmodules; module++) // Loop over all the modules in the shield
			{
				nlayers = CRS->getCRSScintillatorShield(shield).getModule(module).nLayers();
				for (int layer = 0; layer < nlayers; layer++) // Loop over all the layers in the module
				{
					nbars = CRS->getCRSScintillatorShield(shield).getModule(module).getLayer(layer).nBars();
					for (int bar = 0; bar < nbars; bar++) // Loop over all the bars in the layer
					{
						evd_bar = new TEvdCrvBar(CRS->getCRSScintillatorShield(shield).getModule(module).getLayer(layer).getBar(bar), fSectionID);
						fListOfEvdCrvBars->Add(evd_bar);
					}
				}
			}
		}
		break;

	case 2:
		//Top DS
		for (int shield = 10; shield <= 12; shield++) // Loop over all the shields in Top DS
		{
			nmodules = CRS->getCRSScintillatorShield(shield).nModules();
			for (int module = 0; module < nmodules; module++) // Loop over all the modules in the shield
			{
				nlayers = CRS->getCRSScintillatorShield(shield).getModule(module).nLayers();
				for (int layer = 0; layer < nlayers; layer++) // Loop over all the layers in the module
				{
					nbars = CRS->getCRSScintillatorShield(shield).getModule(module).getLayer(layer).nBars();
					for (int bar = 0; bar < nbars; bar++) // Loop over all the bars in the layer
					{
						evd_bar = new TEvdCrvBar(CRS->getCRSScintillatorShield(shield).getModule(module).getLayer(layer).getBar(bar), fSectionID);
						fListOfEvdCrvBars->Add(evd_bar);
					}
				}
			}
		}
		break;

	case 3:
		//Downstream
		nmodules = CRS->getCRSScintillatorShield(13).nModules(); // Get downstream shield
		for (int module = 0; module < nmodules; module++) // Loop over all the modules in the shield
		{
			nlayers = CRS->getCRSScintillatorShield(13).getModule(module).nLayers();
			for (int layer = 0; layer < nlayers; layer++) // Loop over all the layers in the module
			{
				nbars = CRS->getCRSScintillatorShield(13).getModule(module).getLayer(layer).nBars();
				for (int bar = 0; bar < nbars; bar++) // Loop over all the bars in the layer
				{
					evd_bar = new TEvdCrvBar(CRS->getCRSScintillatorShield(13).getModule(module).getLayer(layer).getBar(bar), fSectionID);
					fListOfEvdCrvBars->Add(evd_bar);
				}
			}
		}

		//Downstream lower sections - added 05/04/2015
		for (int shield = 18; shield <= 20; shield++)
		{
			nmodules = CRS->getCRSScintillatorShield(shield).nModules();
			for (int module = 0; module < nmodules; module++)
			{
				nlayers = CRS->getCRSScintillatorShield(shield).getModule(module).nLayers();
				for (int layer = 0; layer < nlayers; layer++)
				{
					nbars = CRS->getCRSScintillatorShield(shield).getModule(module).getLayer(layer).nBars();
					for (int bar = 0; bar < nbars; bar++)
					{
						evd_bar = new TEvdCrvBar(CRS->getCRSScintillatorShield(shield).getModule(module).getLayer(layer).getBar(bar), fSectionID);
						fListOfEvdCrvBars->Add(evd_bar);
					}
				}
			}
		}
		break;

	case 4:
		//Upstream
		nmodules = CRS->getCRSScintillatorShield(14).nModules(); // Get usptream shield
		for (int module = 0; module < nmodules; module++) // Loop over all the modules in the shield
		{
			nlayers = CRS->getCRSScintillatorShield(14).getModule(module).nLayers();
			for (int layer = 0; layer < nlayers; layer++) // Loop over all the layers in the module
			{
				nbars = CRS->getCRSScintillatorShield(14).getModule(module).getLayer(layer).nBars();
				for (int bar = 0; bar < nbars; bar++) // Loop over all the bars in the layer
				{
					evd_bar = new TEvdCrvBar(CRS->getCRSScintillatorShield(14).getModule(module).getLayer(layer).getBar(bar), fSectionID);
					fListOfEvdCrvBars->Add(evd_bar);
				}
			}
		}
		break;

	case 5: // Cryo Upstream
	case 6: // Cryo Downstream
	case 7: // Cryo Top
		// Cryo hole shields are not displayed for now
		break;

	case 8:
		//Top TS
		nmodules = CRS->getCRSScintillatorShield(9).nModules(); // Get Top TS shield
		for (int module = 0; module < nmodules; module++) // Loop over all the modules in the shield
		{
			nlayers = CRS->getCRSScintillatorShield(9).getModule(module).nLayers();
			for (int layer = 0; layer < nlayers; layer++) // Loop over all the layers in the module
			{
				nbars = CRS->getCRSScintillatorShield(9).getModule(module).getLayer(layer).nBars();
				for (int bar = 0; bar < nbars; bar++) // Loop over all the bars in the layer
				{
					evd_bar = new TEvdCrvBar(CRS->getCRSScintillatorShield(9).getModule(module).getLayer(layer).getBar(bar), fSectionID);
					fListOfEvdCrvBars->Add(evd_bar);
				}
			}
		}
		break;

	default:
		Warning("TCrvVisNode", Form("Unknown SectionID %i", SectionID));
		break;

	}
}

TCrvVisNode::~TCrvVisNode()
{
}

int TCrvVisNode::InitEvent()
{	
	ftimeLow = 400;
	ftimeHigh = 1695;

	TEvdCrvBar*              evd_bar;
	mu2e::GeomHandle<mu2e::CosmicRayShield> CRS;

	Clear();

	if (fCrvRecoPulsesCollection) //If we have a valid pointer to a collection (non-empty), then proceed to fill
	{
		//Loop over the RecoPulses in the collection
		for (mu2e::CrvRecoPulsesCollection::const_iterator icrpc = (*fCrvRecoPulsesCollection)->begin(), ecrpc = (*fCrvRecoPulsesCollection)->end(); icrpc != ecrpc; ++icrpc)
		{
			const mu2e::CRSScintillatorBar &CRVCounterBar = CRS->getBar(icrpc->first);

			if (getCRVSection(CRVCounterBar.id().getShieldNumber()) != fSectionID) //If the bar that we are looking at it is not a bar in this CRV Section, skip it
				continue;

			evd_bar = EvdCrvBar(icrpc->first.asInt());			//Set the bar pointer to the bar with map bar index

			const mu2e::CrvRecoPulses &crvRecoPulses = icrpc->second; //The set of pulses for this bar

			//Loop over each SiPM for this bar
			for (unsigned int SiPM = 0; SiPM < 4; SiPM++)
			{
			  const std::vector<mu2e::CrvRecoPulses::CrvSingleRecoPulse> &pulseVector = crvRecoPulses.GetRecoPulses(SiPM);

				//Loop over single pulses
				for (unsigned int i = 0; i < pulseVector.size(); i++)
				{
				  const mu2e::CrvRecoPulses::CrvSingleRecoPulse &pulse = pulseVector[i];
				  evd_bar->AddPulse(pulse, SiPM); // Add the pulse to the bar

				  if ((pulse._PEs > fMinPulsePEs) && (pulse._leadingEdge > ftimeLow) && (pulse._leadingEdge < ftimeHigh)) //If the pulse is above the threshold and within the initial time window, color the sipm
					{						
						evd_bar->SetFillStyle(1001, SiPM);
						evd_bar->SetFillColor(colorPalette[(int) ((pulse._leadingEdge - ftimeLow) / (ftimeHigh - ftimeLow) * 999)], SiPM);
					}
				} // Loop over single pulses
			} // Loop over SiPMs
		} // Loop over RecoPulses
	}
	printf("Finished TCrvVisNode::InitEvent() for section %i \n", fSectionID);

	return 0;
}

void TCrvVisNode::UpdateEvent()
{
	printf("Updating event... TCrvVisNode::UpdateEvent() for section %i \n", fSectionID);
	int          nbar;
	TEvdCrvBar  *evd_bar;

	nbar = fListOfEvdCrvBars->GetEntries();
	for (int ibar = 0; ibar<nbar; ibar++)
	{
		evd_bar = (TEvdCrvBar*) fListOfEvdCrvBars->At(ibar);
		for (unsigned int SiPM = 0; SiPM < 4; SiPM++)
		{
			const mu2e::CrvRecoPulses::CrvSingleRecoPulse* barPulse = evd_bar->lastPulseInWindow((int) SiPM, fMinPulsePEs, ftimeLow, ftimeHigh);
			if (barPulse) //If we have a valid pulse for the bar
			{
				evd_bar->SetFillColor(colorPalette[(int) ((barPulse->_leadingEdge - ftimeLow) / (ftimeHigh - ftimeLow) * 999)], SiPM);
			}
			else // Make the SiPM white since no pulses fall within the window
			{
				evd_bar->SetFillColor(kWhite, SiPM);
			}
		}
	}
	
}

void TCrvVisNode::Paint(Option_t* Option) {

	int iv;

	const char* view = TVisManager::Instance()->GetCurrentView();

	if (strstr(view, "crv") != 0)
	{
		sscanf(view, "crv,%i", &iv);
		if (iv == fSectionID)
			PaintCrv(Option);
	}

	gPad->Modified();
}


//-----------------------------------------------------------------------------
// There will probably be a number of VisNodes, one for each section in the CRV
//-----------------------------------------------------------------------------
void TCrvVisNode::PaintXY(Option_t* Option) {
	// draw crv
}

//_____________________________________________________________________________
void TCrvVisNode::PaintCrv(Option_t* Option) {
	// draw crv
	TEvdCrvBar* evd_bar;

	for (int i = 0; i < fListOfEvdCrvBars->GetEntries(); i++)
	{
		evd_bar = (TEvdCrvBar*) fListOfEvdCrvBars->At(i);

		evd_bar->Paint(Option);
	}
}


//_____________________________________________________________________________
void TCrvVisNode::PaintRZ(Option_t* option) {
	// draw crv
}

//_____________________________________________________________________________
Int_t TCrvVisNode::DistancetoPrimitive(Int_t px, Int_t py) {
	int              min_dist(9999), iv;
	static TVector3  global;

	TVisManager* vm = TVisManager::Instance();
	const char* view = vm->GetCurrentView();

	sscanf(view, "crv,%i", &iv);
	if (iv != fSectionID)                return min_dist;

	global.SetXYZ(gPad->AbsPixeltoX(px), gPad->AbsPixeltoY(py), 0);

	//  printf("px,py,X,Y = %5i %5i %10.3f %10.3f\n",px,py,global.X(),global.Y());

	min_dist = DistancetoPrimitiveXY(px, py);

	return min_dist;
}

//_____________________________________________________________________________
Int_t TCrvVisNode::DistancetoPrimitiveXY(Int_t px, Int_t py) {

	Int_t min_dist = 9999;

	static TVector3 global;
	static TVector3 locrv;
	double          dist;

	fClosestObject = NULL;
	
	global.SetXYZ(gPad->AbsPixeltoX(px), gPad->AbsPixeltoY(py), 0);

	double gx = global.X();
	double gy = global.Y();

	int          nbar;
	TEvdCrvBar  *evd_bar;

	nbar = fListOfEvdCrvBars->GetEntries();
	double gx0, gy0;

	printf("gx: %f , gy %f \n", gx, gy);
	for (int ibar = 0; ibar<nbar; ibar++) {
		evd_bar = (TEvdCrvBar*) fListOfEvdCrvBars->At(ibar);

		//px0 = gPad->XtoAbsPixel(evd_bar->X0());
		//py0 = gPad->YtoAbsPixel(evd_bar->Y0());
		gx0 = evd_bar->X0();
		gy0 = evd_bar->Y0();
		
		dist = sqrt(pow((gx0 - gx),2) + pow((gy0 - gy),2)); //actual distance
		if (dist < min_dist) {
			min_dist = dist;
			fClosestObject = evd_bar;
			printf("Found bar: %i (%f, %f) to be closest at (%f, %f) dist= %f \n", evd_bar->Bar().asInt(), gx0, gy0, gx, gy, dist);
		}
	}
	return min_dist;
}

//_____________________________________________________________________________
Int_t TCrvVisNode::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
	return 9999;
}


//_____________________________________________________________________________
void TCrvVisNode::Clear(Option_t* Opt) {
	int          nbar;
	TEvdCrvBar  *evd_bar;

	nbar = fListOfEvdCrvBars->GetEntries();
	for (int ibar = 0; ibar<nbar; ibar++) 
	{
		evd_bar = (TEvdCrvBar*) fListOfEvdCrvBars->At(ibar);
		evd_bar->Clear();
	}
}

//_____________________________________________________________________________
void TCrvVisNode::Print(Option_t* Opt) const {
	int          nbar;
	TEvdCrvBar  *evd_bar;

	nbar = fListOfEvdCrvBars->GetEntries();
	for (int ibar = 0; ibar<nbar; ibar++)
	{
		evd_bar = (TEvdCrvBar*) fListOfEvdCrvBars->At(ibar);
		evd_bar->Print(Opt);
	}
}

//_____________________________________________________________________________
TEvdCrvBar*  TCrvVisNode::EvdCrvBar(int barIndex)
{
	TEvdCrvBar* tempBarPtr = 0;
	for (int i = 0; i<fListOfEvdCrvBars->GetEntries(); i++)
	{
		tempBarPtr = (TEvdCrvBar*) fListOfEvdCrvBars->At(i);
		if (tempBarPtr->Bar().asInt() == barIndex)
				return tempBarPtr;
	}

	return tempBarPtr; //If it makes it through the loop without returning
}

void TCrvVisNode::SetTimeWindow(float timeLow, float timeHigh)
{
	ftimeLow = timeLow;
	ftimeHigh = timeHigh; 
}

//_____________________________________________________________________________
int TCrvVisNode::getCRVSection(int shieldNumber)
{
	int CRVSection = -1;
	switch (shieldNumber)
	{
	case 0:
	case 1:
	case 2:
	case 3:
	case 4:
	case 5:  CRVSection = 0; break;  //R
	case 6:
	case 7:
	case 8:  CRVSection = 1; break;  //L
	case 9:  CRVSection = 8; break;  //TS T
	case 10:
	case 11:
	case 12: CRVSection = 2; break;  //T
	case 13: CRVSection = 3; break;  //D
	case 14: CRVSection = 4; break;  //U
	case 15: CRVSection = 5; break;  //CU
	case 16: CRVSection = 6; break;  //CD
	case 17: CRVSection = 7; break;  //CT
	case 18:
	case 19:
	case 20: CRVSection = 3; break; //D lower shields - added 05/04/2015
	}
	return CRVSection;
}

//


#include "Stntuple/gui/TEvdMainFrame.hh"
#include "Stntuple/gui/TStnVisManager.hh"
#include "TGButton.h"

ClassImp(TEvdMainFrame)


//_____________________________________________________________________________
TEvdMainFrame::TEvdMainFrame() :
TGMainFrame(gClient->GetRoot(), 100, 200, kMainFrame)
{}

//_____________________________________________________________________________
TEvdMainFrame::TEvdMainFrame(const TGWindow* p, UInt_t w, UInt_t h, Int_t options) :
TGMainFrame(gClient->GetRoot(), w, h, options)
{}




//_____________________________________________________________________________
Bool_t TEvdMainFrame::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2) {
	// Handle menu items.

	TStnVisManager* vm = TStnVisManager::Instance();

	int message = GET_MSG(msg);

	switch (message) {
	case kC_COMMAND:
		int submessage = GET_SUBMSG(msg);
		switch (submessage) {
		case kCM_MENU:
			switch (parm1) {
				//-----------------------------------------------------------------------------
				//  SUBDETECTOR menu
				//-----------------------------------------------------------------------------
			case TStnVisManager::M_TRACKER_XY:
				vm->OpenTrkXYView();
				break;

			case TStnVisManager::M_TRACKER_RZ:
				vm->OpenTrkRZView();
				break;

			case TStnVisManager::M_CALORIMETER_XY:
				vm->OpenCalView();
				break;

			case TStnVisManager::M_CRV_XY:
				vm->OpenCrvView();
				break;

			case TStnVisManager::M_EXIT:
				vm->CloseWindow();
				break;
				//-----------------------------------------------------------------------------
				//  default
				//-----------------------------------------------------------------------------
			default:
				if (vm->DebugLevel() > 0) {
					printf(" *** TStnFrame::ProcessMessage msg: %li parm1: %li parm2: %li\n",
						msg, parm1, parm2);
				}
				break;
			}
		default:
			if (vm->DebugLevel() > 0) {
				printf(" *** TStnFrame::ProcessMessage msg: %li parm1: %li parm2: %li\n",
					msg, parm1, parm2);
			}
			break;
		}
	}
	return true;
}

//_____________________________________________________________________________
void TEvdMainFrame::HandleButtons()
{
	// Handle different buttons.

	TGButton *btn = (TGButton *) gTQSender;
	int id = btn->WidgetId();

	printf("Detected button %i clicked, handling...\n", id);

	switch (id)
	{
	
	case 1:
		TStnVisManager::Instance()->OpenTrkXYView();
		break;
	case 2:
		TStnVisManager::Instance()->OpenCalView();
		break;
	case 3:
		TStnVisManager::Instance()->OpenCrvView();
		break;
	default:
		printf("Unknown button clicked\n");
		break;
	}
}
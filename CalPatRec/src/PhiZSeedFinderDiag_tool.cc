#include "TH1.h"
#include "TH2.h"

#include <string.h>

#include "Offline/CalPatRec/inc/DeltaFinder_types.hh"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Offline/CalPatRec/inc/HlPrint.hh"
#include "Offline/CalPatRec/inc/PhiZSeedFinderAlg.hh"
#include "Offline/CalPatRec/inc/McPart_t.hh"

#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"

#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"

using namespace std;

namespace mu2e {

  using namespace PhiZSeedFinderTypes;

  class SimParticle;
  class StrawDigiMCCollection;

  class PhiZSeedFinderDiag: public ModuleHistToolBase {

    enum {
      kNEventHistSets  =  10,
    };

    struct EventHist_t {
      TH1F*  fEventNumber;
      TH1F*  fRunNumber;
    };

                                        // hits referred to here are the combo hits
    struct Hist_t {
      EventHist_t*  fEvent [kNEventHistSets ];
    };

  protected:

    bool           _mcDiag;

    int            _eventNumber;

    TObjArray      _listOfMcParticles; // list of particles with at least one ComboHit in the tracker

    Data_t*        _data;                 // shared data, passed from the caller
    Hist_t         _hist;

    std::unique_ptr<McUtilsToolBase>      _mcUtils;

    float          _ppi;                // proton pulse intensity

  public:

    PhiZSeedFinderDiag(const fhicl::Table<mu2e::PhiZSeedFinderTypes::Config>& config);
    ~PhiZSeedFinderDiag();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
    McPart_t*   mcPart   (int I) { return (McPart_t*) _listOfMcParticles.At(I); }
//-----------------------------------------------------------------------------
// other functions
//-----------------------------------------------------------------------------
    void        bookEventHistograms (EventHist_t*  Hist, art::TFileDirectory* Dir);

    void        fillEventHistograms (EventHist_t*  Hist);
//-----------------------------------------------------------------------------
// overriden virtual functions of the base class
//-----------------------------------------------------------------------------
  public:

    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override ;
    virtual int fillHistograms(void* Data, int Mode = -1 ) override ;
    virtual int debug         (void* Data, int Mode = -1 ) override ;
  };

//-----------------------------------------------------------------------------
// this routine is called once per job
//-----------------------------------------------------------------------------
  PhiZSeedFinderDiag::PhiZSeedFinderDiag(const fhicl::Table<mu2e::PhiZSeedFinderTypes::Config>& config):
    _mcDiag                (config().mcDiag()                ),
    _printOTracker         (config().printOTracker()         ),
    _printComboHits        (config().printComboHits()        ),
    _printGoodComboHits    (config().printGoodComboHits()    ),
    _printShcol            (config().printShcol()            ),
  {
    printf(" PhiZSeedFinderDiag::PhiZSeedFinderDiag : HOORAY Config! \n");

    if (_mcDiag != 0) _mcUtils = art::make_tool  <McUtilsToolBase>(config().mcUtils,"mcUtils");
    else              _mcUtils = std::make_unique<McUtilsToolBase>();
  }

//-----------------------------------------------------------------------------
  PhiZSeedFinderDiag::~PhiZSeedFinderDiag() {
  }

//-----------------------------------------------------------------------------
// this routine is called once per job (likely, from beginJob)
// TH1::AddDirectory makes sure one can have histograms with the same name
// in different subdirectories
//-----------------------------------------------------------------------------
  int PhiZSeedFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {

    TH1::AddDirectory(0);
    char folder_name[20];
//-----------------------------------------------------------------------------
// book event-level histograms - fill them once per event
//-----------------------------------------------------------------------------
    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[ 0] = 1;                // all events

    for (int i=0; i<kNEventHistSets; i++) {
      if (book_event_histset[i] != 0) {
        sprintf(folder_name,"evt_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist.fEvent[i] = new EventHist_t;
        bookEventHistograms(_hist.fEvent[i],&dir);
      }
    }
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
    return 0;
  }

//-----------------------------------------------------------------------------
  void  PhiZSeedFinderDiag::fillEventHistograms(EventHist_t* Hist) {

    int event_number = _data->event->event();

    Hist->fEventNumber->Fill(event_number);
  }

//-----------------------------------------------------------------------------
// main fill histograms function called once per event
// 'Mode' not used
//-----------------------------------------------------------------------------
  int PhiZSeedFinderDiag::fillHistograms(void* Data, int Mode) {
//-----------------------------------------------------------------------------
// event histograms - just one set
//-----------------------------------------------------------------------------
    fillEventHistograms(_hist.fEvent[0]);
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
    return 0;
  }

//-----------------------------------------------------------------------------
// debugLevel > 0: print seeds
// Mode = 1:
// Mode = 2: event
//-----------------------------------------------------------------------------
  int PhiZSeedFinderDiag::debug(void* Data, int Mode) {
  }

}

DEFINE_ART_CLASS_TOOL(mu2e::PhiZSeedFinderDiag)

///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TH1.h"
#include "TH2.h"

#include <string.h>

#include "CalPatRec/inc/DeltaFinder2_types.hh"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Mu2eUtilities/inc/McUtilsToolBase.hh"
#include "RecoDataProducts/inc/StrawHitPosition.hh"
#include "TrackerGeom/inc/Tracker.hh"

using namespace std;

namespace mu2e {
  using namespace DeltaFinder2Types;

  class SimParticle;
  class StrawDigiMCCollection;
  
  class DeltaFinderDiagCluster: public ModuleHistToolBase {

    enum {
      kNEventHistSets      = 10,
      kNHitClusterHistSets = 10
    };

    struct HitClusterHist_t {
      TH1F*  fSize;
      TH1F*  fUnique;
      TH1F*  fP;                          // momentum of the first cluster particle
      TH1F*  fPdgID;                      // of the first particle
    };

    struct EventHist_t {
      TH1F*  fEventNumber;
      TH1F*  fRunNumber;
      TH1F*  fNHits;
      TH1F*  fNOHits;
      TH1F*  fNHitClusters;
    };


    struct HitCluster_t {
      vector<HitData_t*>       _hits;
      int                      _unique;
      const mu2e::SimParticle* _simp;    // SimParticle corresponding to the first hit
      
      HitCluster_t() { 
	init(); 
      }

      void init() {
	_hits.clear();
	_unique = -1;
	_simp   = NULL;
      }
    };

    
    struct Hist_t {
      EventHist_t*       fEvent     [kNEventHistSets     ];
      HitClusterHist_t*  fHitCluster[kNHitClusterHistSets];
    };
  protected:

    bool                                  _mcDiag;
    int                                   _printShcol;
    int                                   _printOTracker;
    int                                   _printHitClusters;
    art::InputTag                         _stepPointMcCollTag;
    double                                _maxChi2W;

    std::unique_ptr<McUtilsToolBase>      _mcUtils;

    int                                   _eventNumber;
    //    const StrawDigiMCCollection*          _listOfMcStrawHits;
    
    std::vector<McPart_t*>                _list_of_mc_particles; // list_of_particles with hits in the tracker
    std::vector<McPart_t*>                _list_of_mc_part_hit ; // for each StrawHit, pointer to its McPart 

    Data_t*                               _data;                 // diag data, passed from the caller, cached

    Hist_t                                _hist;

    vector<HitCluster_t>                  _list_of_hclusters[kNStations][kNFaces][kNPanelsPerFace];

    int                                   _nsh;  // number of straw hits
    int                                   _nosh;  // number of straw hits in the oTracker structure
    int                                   _nhcl; // number of hit clusters
    
  public:
    
    DeltaFinderDiagCluster(const fhicl::ParameterSet& PSet);
    ~DeltaFinderDiagCluster();

  private:

    void        bookEventHistograms     (EventHist_t*      Hist, art::TFileDirectory* Dir);
    void        bookHitClusterHistograms(HitClusterHist_t* Hist, art::TFileDirectory* Dir);

    void        fillEventHistograms     (EventHist_t* Hist);
    void        fillHitClusterHistograms(HitClusterHist_t* Hist, HitCluster_t*        Mc );

    int         associateMcTruth();
    
    void        printHitData(const HitData_t* Hd, int Index);
    void        printHitClusters();
    void        printOTracker();

    // void        printStrawHit(const StrawHit* Sh, int Index);
    // void        printStrawHitCollection();
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
  DeltaFinderDiagCluster::DeltaFinderDiagCluster(const fhicl::ParameterSet& PSet):
    _mcDiag                (PSet.get<bool>         ("mcDiag"                       )),
    _printOTracker         (PSet.get<int>          ("printOTracker"                )),
    _printHitClusters      (PSet.get<int>          ("printHitClusters"             )),
    _stepPointMcCollTag    (PSet.get<string>       ("stepPointMcCollTag"           )),
    _maxChi2W              (PSet.get<double>       ("maxChi2W"                     ))
  {
    printf(" DeltaFinderDiagCluster::DeltaFinderDiagCluster : HOORAY! \n");

    if (_mcDiag != 0) _mcUtils = art::make_tool<McUtilsToolBase>(PSet.get<fhicl::ParameterSet>("mcUtils"));
    else              _mcUtils = std::make_unique<McUtilsToolBase>();
  }
  
//-----------------------------------------------------------------------------
  DeltaFinderDiagCluster::~DeltaFinderDiagCluster() {
  }
  
//-----------------------------------------------------------------------------
  void DeltaFinderDiagCluster::bookEventHistograms(EventHist_t* Hist, art::TFileDirectory* Dir) {

    Hist->fEventNumber     = Dir->make<TH1F>("event" , "Event Number" , 100, 0., 100000.);
    Hist->fRunNumber       = Dir->make<TH1F>("run"   , "Run Number"   , 100, 0., 100000.);
    Hist->fNHits           = Dir->make<TH1F>("nhits" , "nhits"        , 500, 0., 10000.);
    Hist->fNOHits          = Dir->make<TH1F>("nohits", "nohits"       , 500, 0., 10000.);
    Hist->fNHitClusters    = Dir->make<TH1F>("nhcl"  , "nhit clusters", 500, 0., 10000.);

  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiagCluster::bookHitClusterHistograms(HitClusterHist_t* Hist, art::TFileDirectory* Dir) {

    Hist->fSize    = Dir->make<TH1F>("size"  , "Cluster size",  25, 0.,  25.);
    Hist->fUnique  = Dir->make<TH1F>("unique", "unique"      ,  10, 0.,  10.);
    Hist->fP       = Dir->make<TH1F>("p"     , "momentum"    , 250, 0., 250.);
    Hist->fPdgID   = Dir->make<TH1F>("pdg_id", "PDG ID"      , 1000, -5000., 5000.);
  }

//-----------------------------------------------------------------------------
// this routine is called once per job (likely, from beginJob)
// TH1::AddDirectory makes sure one can have histograms with the same name
// in different subdirectories
//-----------------------------------------------------------------------------
  int DeltaFinderDiagCluster::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {

    TH1::AddDirectory(0);
    char folder_name[20];
//-----------------------------------------------------------------------------
// book event-level histograms - fill them once per event
//-----------------------------------------------------------------------------
    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[ 0] = 1;		// all events

    for (int i=0; i<kNEventHistSets; i++) {
      if (book_event_histset[i] != 0) {
	sprintf(folder_name,"evt_%i",i);
	art::TFileDirectory dir = Tfs->mkdir(folder_name);
	
	_hist.fEvent[i] = new EventHist_t;
	bookEventHistograms(_hist.fEvent[i],&dir);
      }
    }
//-----------------------------------------------------------------------------
// book MC histograms
//-----------------------------------------------------------------------------
    int book_hit_cluster_histset[kNHitClusterHistSets];
    for (int i=0; i<kNHitClusterHistSets; i++) book_hit_cluster_histset[i] = 0;

    book_hit_cluster_histset[  0] = 1;		// all hit clusters
    book_hit_cluster_histset[  1] = 1;		// unique hit clusters
    book_hit_cluster_histset[  2] = 1;		// merged hit clusters

    for (int i=0; i<kNHitClusterHistSets; i++) {
      if (book_hit_cluster_histset[i] != 0) {
	sprintf(folder_name,"hcl_%i",i);
	art::TFileDirectory dir = Tfs->mkdir(folder_name);

	_hist.fHitCluster[i] = new HitClusterHist_t;
	bookHitClusterHistograms(_hist.fHitCluster[i],&dir);
      }
    }
    return 0;
  }


//-----------------------------------------------------------------------------
  void  DeltaFinderDiagCluster::fillEventHistograms(EventHist_t* Hist) {

    int event_number = _data->event->event();
    
    Hist->fEventNumber->Fill(event_number);
    Hist->fNHits->Fill(_nsh);
    Hist->fNOHits->Fill(_nosh);
    Hist->fNHitClusters->Fill(_nhcl);

  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiagCluster::fillHitClusterHistograms(HitClusterHist_t* Hist, HitCluster_t* HCl) {

    Hist->fSize->Fill(HCl->_hits.size());
    Hist->fUnique->Fill(HCl->_unique);

    const SimParticle* simp = HCl->_simp;
    int pdg_id = _mcUtils->getPdgID(simp);
    double mom = _mcUtils->getStartMom(simp);

    Hist->fPdgID->Fill(pdg_id);
    Hist->fP    ->Fill(mom);
  }

//-----------------------------------------------------------------------------
// main fill histograms function called once per event
// 'Mode' not used
//-----------------------------------------------------------------------------
  int DeltaFinderDiagCluster::fillHistograms(void* Data, int Mode) {
    _data = (Data_t*) Data;
//-----------------------------------------------------------------------------
// start from precalculating MC-specific info
//-----------------------------------------------------------------------------
    int en = _data->event->event();
    if (_mcDiag) {
      if (_eventNumber != en) {
	_eventNumber       = en;
	//	_listOfMcStrawHits = _mcUtils->getListOfMcStrawHits(_data->event, _stepPointMcCollTag);
	associateMcTruth();
      }
    }

    _nsh = _data->shcol->size();

    _nosh = 0;
    _nhcl = 0;

    for (int station=0; station<kNStations; station++) {
      for (int face=0; face<kNFaces; face++) {
	for (int ip=0; ip<3; ip++) {
	  _nhcl += _list_of_hclusters[station][face][ip].size();

	  PanelZ_t* pz = &_data->oTracker[station][face][ip];
	  int nosh      = pz->fHitData.size();

	  _nosh += nosh;
	}
      }
    }
//-----------------------------------------------------------------------------
// event histograms - just one set
//-----------------------------------------------------------------------------
    fillEventHistograms(_hist.fEvent[0]);
//-----------------------------------------------------------------------------
// fill hit cluster histograms
//-----------------------------------------------------------------------------
    for (int station=0; station<kNStations; station++) {
      for (int face=0; face<kNFaces; face++) {
	for (int ip=0; ip<3; ip++) {
	  //	  int ncl = _list_of_hclusters[station][face][ip].size();
	  for (auto hcl : _list_of_hclusters[station][face][ip]) {
	    fillHitClusterHistograms(_hist.fHitCluster[0],&hcl);

	    if (hcl._unique == 1) fillHitClusterHistograms(_hist.fHitCluster[1],&hcl);
	    else                  fillHitClusterHistograms(_hist.fHitCluster[2],&hcl);
	  }
	}
      }
    }

    return 0;
  }

//-----------------------------------------------------------------------------
// for each DeltaSeed, create a list of SimParticle*'s parallel to its list of straw hits
//-----------------------------------------------------------------------------
  int DeltaFinderDiagCluster::associateMcTruth() {

    const StrawHit* sh0 = &_data->shcol->at(0);

    for (int station=0; station<kNStations; ++station) {
      for (int face=0; face<kNFaces; ++face) {
	for (int ip=0; ip<3; ++ip) {
	  vector<HitCluster_t>* list = & _list_of_hclusters[station][face][ip];

	  list->clear();
//-----------------------------------------------------------------------------
// count clusters in a given panel
//-----------------------------------------------------------------------------
	  PanelZ_t* pz = &_data->oTracker[station][face][ip];
	  int nh = pz->fHitData.size();

	  for (int ih=0; ih<nh; ih++) {
	    HitData_t* hd  = &pz->fHitData[ih];
	    uint16_t is    = hd->fStraw->id().straw();
//------------------------------------------------------------------------------
// first , see if this hit belongs to already found hit cluster
//-----------------------------------------------------------------------------
	    int added = 0;
	    for (auto hit_cluster = list->begin(); hit_cluster != list->end(); hit_cluster++) {
	      HitData_t* last_hit = hit_cluster->_hits.back();
	      uint16_t lstraw = last_hit->fStraw->id().straw();
	      if (is-lstraw < 3) {
//-----------------------------------------------------------------------------
// hit is close, check time and z
//-----------------------------------------------------------------------------
		float dt = hd->fHit->time()-last_hit->fHit->time();
		float dw = hd->fPos->wireDist()-last_hit->fPos->wireDist();

		double sigw1 =  hd->fPos->posRes(StrawHitPosition::wire);
		double sigw2 =  last_hit->fPos->posRes(StrawHitPosition::wire);

		double chi2w = dw*dw/(sigw1*sigw1+sigw2*sigw2);

		if ((fabs(dt) < 50) && (chi2w < _maxChi2W)) {
		  hit_cluster->_hits.push_back(hd);
		  added = 1;
		  break;
		}
	      }
	    }

	    if (added == 1)                                   continue;
//-----------------------------------------------------------------------------
// new hit cluster
//-----------------------------------------------------------------------------
	    HitCluster_t hcl;

	    hcl._hits.push_back(hd);

	    const StrawHit* sh            = hd->fHit;
	    int i0                        = sh-sh0;
	    const mu2e::SimParticle* simp = _mcUtils->getSimParticle(_data->event,i0);
	    hcl._simp                     = simp;

	    _list_of_hclusters[station][face][ip].push_back(hcl);
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// at this point, lists of hit clusters are formed, need to check whether 
// all hits in a cluster belong to the same particle
//-----------------------------------------------------------------------------
    for (int station=0; station<kNStations; station++) {
      for (int face=0; face<kNFaces; face++) {
	for (int ip=0; ip<3; ++ip) {
	  int ncl = _list_of_hclusters[station][face][ip].size();
	  for (int icl=0; icl<ncl; icl++) {
	    HitCluster_t* hcl = &_list_of_hclusters[station][face][ip][icl];
	    int nhits = hcl->_hits.size();

	    const StrawHit* sh           = hcl->_hits.front()->fHit;
	    int i0                       = sh-sh0;
	    const mu2e::SimParticle* mc0 = _mcUtils->getSimParticle(_data->event,i0);

	    hcl->_unique = 1;

	    for (int ih=1; ih<nhits; ih++) {
	      const StrawHit* sh           = hcl->_hits[ih]->fHit;
	      int i1                       = sh-sh0;
	      const mu2e::SimParticle* mc1 = _mcUtils->getSimParticle(_data->event,i1);

	      if (mc1 != mc0) {
		hcl->_unique = 0;
		break;
	      }
	    }
	  }
	}
      }
    }

    return 0;
  }
  
//-----------------------------------------------------------------------------
// debugLevel > 0: print seeds
//-----------------------------------------------------------------------------
  int DeltaFinderDiagCluster::debug(void* Data, int Mode) {

    _data = (Data_t*) Data;

    if (Mode != 2) return -1; // beginRun not handled yet

    int en = _data->event->event();
    if (_mcDiag) {
      if (_eventNumber != en) {
	_eventNumber       = en;
	//	_listOfMcStrawHits = _mcUtils->getListOfMcStrawHits(_data->event, _stepPointMcCollTag);
	//	InitMcDiag();
	associateMcTruth();
      }
    }

    if (_printHitClusters > 0) printHitClusters();
    if (_printOTracker    > 0) printOTracker   ();

    // if (_printShcol) printStrawHitCollection();

    return 0;
  }


//-----------------------------------------------------------------------------
  void DeltaFinderDiagCluster::printHitData(const HitData_t* Hd, int Index) {

    if (Index < 0) {
      printf("#-----------------------------------------------------------------------------");
      printf("------------------------------------------------------------------------------\n");
      printf("#      SHID  St:Pl P L Str     Time     dt        eDep       wdist     wres   ");
      printf("     PDG           ID       p      X        Y         Z   DeltaID radOK edepOK\n");
      printf("#-----------------------------------------------------------------------------");
      printf("------------------------------------------------------------------------------\n");
      return;
    }

    const StrawHit* sh0 = &_data->shcol->at(0);
    const StrawHit* sh  = Hd->fHit;
    int loc             = sh-sh0;

    const StrawHitPosition* shp = Hd->fPos;
    const StrawHitFlag*     shf = &_data->shfcol->at(loc);
   
    int radselOK        = (! shf->hasAnyProperty(StrawHitFlag::radsel));
    int edepOK          = (! shf->hasAnyProperty(StrawHitFlag::energysel));

    const SimParticle* sim(0);
    int                pdg_id(-9999), sim_id(-9999);
    float              mc_mom(-9999.);
	
    if (_mcDiag) {
      sim    = _mcUtils->getSimParticle(_data->event,loc);
      pdg_id = _mcUtils->getPdgID(sim);
      sim_id = _mcUtils->getID(sim);
      mc_mom = _mcUtils->getStartMom(sim);
    }

    const mu2e::Straw* straw = &_data->tracker->getStraw(sh->strawId());
    
    printf("%5i ",loc);
    printf("%5i" ,sh->strawId().asUint16());
	
    printf("  %2i:%2i %1i %1i %2i   %8.3f %7.3f  %9.6f   %8.3f %8.3f %10i   %10i %8.3f %8.3f %8.3f %9.3f %5i %5i %5i\n",
	   straw->id().getStation(),
	   straw->id().getPlane(),
	   straw->id().getPanel(),
	   straw->id().getLayer(),
	   straw->id().getStraw(),
	   sh->time(),
	   sh->dt(),
	   sh->energyDep(),
	   shp->wireDist(),
	   shp->posRes(StrawHitPosition::wire),
	   pdg_id,
	   sim_id,
	   mc_mom,
	   shp->pos().x(),
	   shp->pos().y(),
	   shp->pos().z(),
	   Hd->fDeltaIndex,
	   radselOK,
	   edepOK);
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void DeltaFinderDiagCluster::printOTracker() {

    printf(" >>> DeltaFinderDiagCluster::printOTracker\n");
    int nhitso = 0;
    for (int is=0; is<kNStations; is++) {
      for (int face=0; face<kNFaces; face++) {
	for (int ip=0; ip<kNPanelsPerFace; ip++) {
	  PanelZ_t* pz = &_data->oTracker[is][face][ip];
	  printf("#        --------------- station: %2i face: %2i panel: %2i nhits[0]:%3li\n",
		 is,face,ip, pz->fHitData.size());
	  printHitData(NULL,-1);
	  int nh = pz->fHitData.size();
	  for (int ih=0; ih<nh; ih++) {
	    printHitData(&pz->fHitData[ih],ih);
	  }
	  nhitso += nh;
	}
      }
    }

    printf(" nhits, nhitso : %6i %6i \n", (int) _data->shcol->size(),nhitso);
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void DeltaFinderDiagCluster::printHitClusters() {

    printf(" >>> DeltaFinderDiagCluster::printHitClusters\n");
    int nhitso = 0;
    for (int is=0; is<kNStations; is++) {
      for (int face=0; face<kNFaces; face++) {
	for (int ip=0; ip<kNPanelsPerFace; ip++) {
	  vector<HitCluster_t>* list = &_list_of_hclusters[is][face][ip];
	  int nhcl = list->size();
	  printf("#        --------------- station: %2i face: %2i panel: %2i N(clusters): %3i\n",
		 is,face,ip, nhcl);

	  for (int ihcl=0; ihcl<nhcl; ihcl++) {
	    HitCluster_t* hcl = &(*list)[ihcl];
	    int nh = hcl->_hits.size();
	    printf(" -------------------------- hit cluster # %2i nhits = %i unique=%i\n",ihcl,nh,hcl->_unique);
	    printHitData(NULL,-1);
	    for (int ih=0; ih<nh; ih++) {
	      printHitData(hcl->_hits[ih],ih);
	    }
	    nhitso += nh;
	  }
	}
      }
    }

    printf(" nhits, nhitso : %6i %6i \n", (int) _data->shcol->size(),nhitso);
  }
}

DEFINE_ART_CLASS_TOOL(mu2e::DeltaFinderDiagCluster)

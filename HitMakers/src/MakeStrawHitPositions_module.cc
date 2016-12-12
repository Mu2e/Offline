//
// A module to create strawhitpositions out of the strawhits
//
// $Id: MakeStrawHitPositions_module.cc,v 1.3 2014/04/28 13:34:43 brownd Exp $
// $Author: brownd $
// $Date: 2014/04/28 13:34:43 $
// 
//  Original Author: David Brown, LBNL
//  changes from G. Pezzullo, P. Murat

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

// art includes.
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// root
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TMarker.h"
#include "TList.h"
#include "TLegend.h"

// C++ includes.
#include <iostream>
#include <float.h>

//#include "CalPatRec/inc/THackData.hh"

using namespace std;

namespace mu2e {


  void printShP(StrawHitPosition shp){
    double x,y,z;
    x = shp.pos().x();
    y = shp.pos().y();
    z = shp.pos().z();
    
    printf("%5.3f  %5.3f  %5.3f\n", x,y,z );
  }

  class MakeStrawHitPositions : public art::EDProducer {

  public:
    explicit MakeStrawHitPositions(fhicl::ParameterSet const& pset);
    // Accept compiler written d'tor.

    void produce( art::Event& e);
    void beginJob();
    void beginRun(art::Run& Run);
    void printHits(const StrawHit& Hit, StrawHitPosition& Pos, int &banner);
  private:

    // Diagnostics level.
    int _diagLevel;
    int _printHits;
    // Name of the StrawHit collection
    string _shLabel;
    string _shfLabel;
    // diagnostics
    TH1F* _nhits;
    TH1F* _deltat;
    TH1F* _deltaE;
    TH1F* _dperp;
    TH1F* _sep;
    TH1F* _dTD;
    TH1F* _chisq;
    std::vector<TH2F*> _stations;

    const Tracker* _tracker;
    void drawStations(); 

    //    THackData* fHackData;
  };

  MakeStrawHitPositions::MakeStrawHitPositions(fhicl::ParameterSet const& pset) :

    // Parameters
    _diagLevel(pset.get<int>("diagLevel",0)),
    _printHits(pset.get<int>("printHits",0)),
    _shLabel(pset.get<string>("StrawHitCollectionLabel","makeSH")),
    _shfLabel(pset.get<string>("StrawHitFlagCollectionLabel","FlagBkgHits")),
    _nhits(0),_deltat(0),_deltaE(0),_dperp(0),_sep(0),_dTD(0),_chisq(0)
 {
    // Tell the framework what we make.
    produces<StrawHitPositionCollection>();

    //    fHackData = new THackData("HackData","Hack Data");
    //    gROOT->GetRootFolder()->Add(fHackData);
 
  }

   void MakeStrawHitPositions::beginJob(){
    // create diagnostics if requested
    if(_diagLevel > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _nhits = tfs->make<TH1F>("nhits","NHits",500,0,5000);
      _deltat = tfs->make<TH1F>("deltat","#Delta t;ns",100,-200.0,200.0);
      _deltaE = tfs->make<TH1F>("deltaE","#Delta E;MeV",100,0.0,0.05);
      _dperp = tfs->make<TH1F>("dperp","#Delta d;mm",100,0.0,500.0);
      _sep = tfs->make<TH1F>("sep","Face separation",6,-0.5,5.5);
      _dTD = tfs->make<TH1F>("dTD","#Delta Time Difference;mm",100,-100.0,100.0);
      _chisq = tfs->make<TH1F>("chisq","Chisquared",100,-1.0,50.0);
    }

  }

  void MakeStrawHitPositions::beginRun(art::Run& Run) {
    _tracker = &getTrackerOrThrow();
  }

  void MakeStrawHitPositions::printHits(const StrawHit& Hit, StrawHitPosition& Pos, int &banner) {
    double x,y,z;

    Straw const& straw = _tracker->getStraw(Hit.strawIndex());
    x = Pos.pos().x();
    y = Pos.pos().y();
    z = Pos.pos().z();

    if(banner==0){
      printf("-----------------------------------------------------------------------------------");
      printf("------------------------------\n");
      printf("   x       y     z   SHID    Station Panel Layer Straw     Flags      Time          dt       eDep \n");
      printf("-----------------------------------------------------------------------------------");
      printf("------------------------------\n");
      banner = false;
    }

    printf("%5.3f %5.3f %5.3f %5i  %5i  %5i   %5i   %5i   0x%s %8.3f   %8.3f   %9.6f\n",
	   x,y,z,
	   Hit.strawIndex().asInt(),
	   straw.id().getPlane(),
	   straw.id().getPanel(),
	   straw.id().getLayer(),
	   straw.id().getStraw(),
	   Pos.flag().hex().data(),
	   Hit.time(),
	   Hit.dt(),
	   Hit.energyDep() );

  }

  void MakeStrawHitPositions::produce(art::Event& event) {

    if ( _diagLevel > 0 ) cout << "MakeStrawHitPositions: produce() begin; event " << event.id().event() << endl;


    static bool first(true);
    if(_diagLevel >1 && first){
      first = false;
      drawStations();
   }

   // Handle to the conditions service
    ConditionsHandle<TrackerCalibrations> tcal("ignored");

    const StrawHitFlagCollection *_shfcol;
    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if(event.getByLabel(_shfLabel,shflagH))
      _shfcol = shflagH.product();
    else 
      _shfcol = 0;

    art::Handle<mu2e::StrawHitCollection> strawhitsH; 
    const StrawHitCollection* strawhits(0);
    if(event.getByLabel(_shLabel,strawhitsH))
      strawhits = strawhitsH.product();
    if(strawhits == 0){
      throw cet::exception("RECO")<<"mu2e::MakeStrawHitPositions: No StrawHit collection found for label " <<  _shLabel << endl;
    }
    // create a collection of StrawHitPosition, and intialize them using the time division
    size_t nsh = strawhits->size();
    if(_diagLevel > 0)_nhits->Fill(nsh);

    unique_ptr<StrawHitPositionCollection> shpos(new StrawHitPositionCollection);
    shpos->reserve(2*nsh);


 //01 - 13 - 2014 gianipez added some printout
    int banner(0);
    SHInfo shinfo;

    for(size_t ish=0;ish<nsh;++ish){
      StrawHit const& hit = strawhits->at(ish);
      Straw const& straw = _tracker->getStraw(hit.strawIndex());
      tcal->StrawHitInfo(straw,hit,shinfo);
      StrawHitFlag shflag = _shfcol->at(ish);
      StrawHitPosition shp(hit,straw,shinfo,shflag);

      if (_printHits>0) {
	printHits(hit,shp, banner);
	banner=1;
      }

      shpos->push_back(shp);
    }

    event.put(std::move(shpos));
  } // end MakeStrawHitPositions::produce.

  void MakeStrawHitPositions::drawStations() {
    // Get a reference to one of the L or T trackers.
    // Throw exception if not successful.
    const TTracker& tt = dynamic_cast<const TTracker&>(*_tracker);
    art::ServiceHandle<art::TFileService> tfs;
    int nsta = tt.nPlanes()/2;
    for(int ista=0;ista<nsta;++ista){
      char name[100];
      snprintf(name,100,"station%i",ista);
      _stations.push_back( tfs->make<TH2F>(name,name,100,-700,700,100,-700,700));
      _stations[ista]->SetStats(false);
      TList* flist = _stations[ista]->GetListOfFunctions();
      TLegend* sleg = new TLegend(0.1,0.6,0.3,0.9);
      flist->Add(sleg);
      for(int iplane=0;iplane<2;++iplane){
	const Plane& pln = tt.getPlane(2*ista+iplane);
	const std::vector<Panel>& panels = pln.getPanels();
	for(size_t ipnl=0;ipnl<panels.size();++ipnl){
	  int iface = ipnl%2;
	  const Panel& pnl = panels[ipnl];
	  CLHEP::Hep3Vector spos = pnl.straw0MidPoint();
	  CLHEP::Hep3Vector sdir = pnl.straw0Direction();
	  CLHEP::Hep3Vector end0 = spos - 100.0*sdir;
	  CLHEP::Hep3Vector end1 = spos + 100.0*sdir;
	  TLine* sline = new TLine(end0.x(),end0.y(),end1.x(),end1.y());
	  sline->SetLineColor(ipnl+1);
	  sline->SetLineStyle(2*iplane+iface+1);
	  flist->Add(sline);
	  TMarker* smark = new TMarker(end0.x(),end0.y(),8);
	  smark->SetMarkerColor(ipnl+1);
	  smark->SetMarkerSize(2);
	  flist->Add(smark);
	  char label[80];
	  double phi = spos.phi();
	  if(phi<0)phi+= 2*M_PI;
	  snprintf(label,80,"pln %i pnl %i phi= %f",iplane,(int)ipnl,180*phi/M_PI);
	  sleg->AddEntry(sline,label,"l");
	}
      }
    }
  }
} // end namespace mu2e

using mu2e::MakeStrawHitPositions;
DEFINE_ART_MODULE(MakeStrawHitPositions)


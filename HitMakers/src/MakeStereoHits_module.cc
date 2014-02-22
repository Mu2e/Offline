//
// A module to create simple stereo hits out of StrawHits.  This can work
// with either tracker.  StrawHit selection is done by flagging in an upstream module
//
// $Id: MakeStereoHits_module.cc,v 1.12 2014/02/22 13:37:36 brownd Exp $
// $Author: brownd $
// $Date: 2014/02/22 13:37:36 $
// 
//  Original Author: David Brown, LBNL
//  

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

// art includes.
#include "art/Persistency/Common/Ptr.h"
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
#
// C++ includes.
#include <iostream>
#include <float.h>

using namespace std;

namespace mu2e {

  class MakeStereoHits : public art::EDProducer {

  public:
    explicit MakeStereoHits(fhicl::ParameterSet const& pset);
    // Accept compiler written d'tor.

    void produce( art::Event& e);
    void beginJob();

  private:

    // Diagnostics level.
    int _diagLevel;
    // Name of the StrawHit collection
    string _shLabel;
  // Parameters
    double _maxDt; // maximum time separation between hits
    double _maxDE; // maximum deposited energy deference: this excludes inconsistent hits
    double _maxDPerp; // maximum transverse separation
    double _maxTDNsig; // maximum # of TimeDivision sigmas past the active edge of a wire to allow making stereo hits
    double _maxChisq; // maximum # of TimeDivision consistency chisquared to allow making stereo hits
    bool _writepairs; // write out the stereo pairs
    // diagnostics
    TH1F* _nhits;
    TH1F* _deltat;
    TH1F* _deltaE;
    TH1F* _dperp;
    TH1F* _sep;
    TH1F* _dTD;
    TH1F* _chisq;
    std::vector<TH2F*> _stations;
  };

  MakeStereoHits::MakeStereoHits(fhicl::ParameterSet const& pset) :

    // Parameters
    _diagLevel(pset.get<int>("diagLevel",0)),
    _shLabel(pset.get<string>("StrawHitCollectionLabel","makeSH")),
    _maxDt(pset.get<double>("maxDt",60.0)), // nsec
    _maxDE(pset.get<double>("maxDE",0.99)), // dimensionless, from 0 to 1
    _maxDPerp(pset.get<double>("maxDPerp",100)), // mm, maximum perpendicular distance between time-division points
    _maxTDNsig(pset.get<double>("maxTDNsig",5.0)),
    _maxChisq(pset.get<double>("maxChisquared",12.0)),
    _writepairs(pset.get<bool>("WriteStereoPairs",false)),
    _nhits(0),_deltat(0),_deltaE(0),_dperp(0),_sep(0),_dTD(0),_chisq(0)
 {
    // Tell the framework what we make.
    if(_writepairs)produces<StereoHitCollection>();
    produces<StrawHitPositionCollection>();
  }

   void MakeStereoHits::beginJob(){
    // create diagnostics if requested
    if(_diagLevel > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _nhits = tfs->make<TH1F>("nhits","NHits",500,0,5000);
      _deltat = tfs->make<TH1F>("deltat","#Delta t;ns",100,-200.0,200.0);
      _deltaE = tfs->make<TH1F>("deltaE","#Delta E;MeV",100,-0.2,0.2);
      _dperp = tfs->make<TH1F>("dperp","#Delta d;mm",100,0.0,500.0);
      _sep = tfs->make<TH1F>("sep","Face separation",6,-0.5,5.5);
      _dTD = tfs->make<TH1F>("dTD","#Delta Time Difference;mm",100,-100.0,100.0);
      _chisq = tfs->make<TH1F>("chisq","Chisquared",100,-1.0,50.0);
    }
  }

  void
  MakeStereoHits::produce(art::Event& event) {

    if ( _diagLevel > 0 ) cout << "MakeStereoHits: produce() begin; event " << event.id().event() << endl;

       // Get a reference to one of the L or T trackers.
       // Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();
    static bool first(true);
    if(_diagLevel >1 && first){
      first = false;
      const TTracker& tt = dynamic_cast<const TTracker&>(tracker);
      art::ServiceHandle<art::TFileService> tfs;
      unsigned nsta = tt.nDevices()/2;
      for(int ista=0;ista<nsta;++ista){
	char name[100];
	snprintf(name,100,"station%i",ista);
	_stations.push_back( tfs->make<TH2F>(name,name,100,-700,700,100,-700,700));
	_stations[ista]->SetStats(false);
	TList* flist = _stations[ista]->GetListOfFunctions();
	TLegend* sleg = new TLegend(0.1,0.6,0.3,0.9);
	flist->Add(sleg);
	for(int idev=0;idev<2;++idev){
	  const Device& dev = tt.getDevice(2*ista+idev);
	  const std::vector<Sector>& sectors = dev.getSectors();
	  for(size_t isec=0;isec<sectors.size();++isec){
	    int iface = isec%2;
	    const Sector& sec = sectors[isec];
	    CLHEP::Hep3Vector spos = sec.straw0MidPoint();
	    CLHEP::Hep3Vector sdir = sec.straw0Direction();
	    CLHEP::Hep3Vector end0 = spos - 100.0*sdir;
	    CLHEP::Hep3Vector end1 = spos + 100.0*sdir;
	    TLine* sline = new TLine(end0.x(),end0.y(),end1.x(),end1.y());
	    sline->SetLineColor(isec+1);
	    sline->SetLineStyle(2*idev+iface+1);
	    flist->Add(sline);
	    TMarker* smark = new TMarker(end0.x(),end0.y(),8);
	    smark->SetMarkerColor(isec+1);
	    smark->SetMarkerSize(2);
	    flist->Add(smark);
	    char label[80];
	    snprintf(label,80,"dev %i sec %i",idev,(int)isec);
	    sleg->AddEntry(sline,label,"l");
	  }
	}
      }
    }

   // Handle to the conditions service
    ConditionsHandle<TrackerCalibrations> tcal("ignored");

    art::Handle<mu2e::StrawHitCollection> strawhitsH; 
    const StrawHitCollection* strawhits(0);
    if(event.getByLabel(_shLabel,strawhitsH))
      strawhits = strawhitsH.product();
    if(strawhits == 0){
      throw cet::exception("RECO")<<"mu2e::MakeStereoHits: No StrawHit collection found for label " <<  _shLabel << endl;
    }
    // create a collection of StrawHitPosition, and intialize them using the time division
    size_t nsh = strawhits->size();
    if(_diagLevel > 0)_nhits->Fill(nsh);

    unique_ptr<StrawHitPositionCollection> shpos(new StrawHitPositionCollection);
    shpos->reserve(2*nsh);
    for(size_t ish=0;ish<nsh;++ish){
      StrawHit const& hit = strawhits->at(ish);
      Straw const& straw = tracker.getStraw(hit.strawIndex());
      SHInfo shinfo;
      tcal->StrawHitInfo(straw,hit,shinfo);
      shpos->push_back(StrawHitPosition(hit,straw,shinfo));
    }
    // create the stereo hits
    StereoHitCollection stereohits;
    stereohits.reserve(3*nsh);
    vector<double> minchisq(nsh,FLT_MAX);
    vector<int> minsep(nsh,SectorId::apart);
    vector<size_t> ibest(nsh);
    // double loop over selected straw hits
    for(size_t ish=0;ish<nsh;++ish){
      StrawHit const& sh1 = strawhits->at(ish);
      Straw const& straw1 = tracker.getStraw(sh1.strawIndex());
      StrawHitPosition const& shp1 = (*shpos)[ish];
      for(size_t jsh=ish+1;jsh<nsh;++jsh){
	StrawHit const& sh2 = strawhits->at(jsh);
	Straw const& straw2 = tracker.getStraw(sh2.strawIndex());
	StrawHitPosition const& shp2 = (*shpos)[jsh];
	SectorId::isep sep = straw1.id().getSectorId().separation(straw2.id().getSectorId());
	if(_diagLevel > 1 &&  sep != SectorId::same && sep < SectorId::apart) {
	  _deltat->Fill(sh1.time()-sh2.time());
	  _deltaE->Fill(sh1.energyDep() - sh2.energyDep());
	  _dperp->Fill((shp1.pos()-shp2.pos()).perp());
	}
	if( sep != SectorId::same && sep < SectorId::apart // hits are in the same station but not the same sector
	    && (sep <= minsep[ish] || sep <= minsep[jsh]) // this separation is at least as good as the current best for one of the hits
	    && fabs(sh1.time()-sh2.time()) < _maxDt // hits are close in time
	    && fabs(sh1.energyDep() - sh2.energyDep())/(sh1.energyDep()+sh2.energyDep()) < _maxDE   // deposited energy is roughly consistent (should compare dE/dx but have no path yet!)
	    && (shp1.pos()-shp2.pos()).perp() < _maxDPerp) { // transverse separation isn't too big
  // tentative stereo hit: this solves for the POCA
	  StereoHit sth(*strawhits,tracker,ish,jsh);
	  if(_diagLevel > 1 ) {
	    _dTD->Fill( straw1.getDetail().activeHalfLength()-fabs(sth.wdist1()) + straw1.getDetail().outerRadius()); 
	    _dTD->Fill( straw2.getDetail().activeHalfLength()-fabs(sth.wdist2()) + straw2.getDetail().outerRadius());
	  }
	  if( straw1.getDetail().activeHalfLength()-fabs(sth.wdist1()) > -straw1.getDetail().outerRadius() // stereo point is inside the active length
	      && straw2.getDetail().activeHalfLength()-fabs(sth.wdist2()) > -straw2.getDetail().outerRadius()) { // stereo point is inside the active length
	    // compute difference between stereo points and TD prediction
	    double chi1 = (shp1.wireDist()-sth.wdist1())/shp1.posRes(StrawHitPosition::phi);	      
	    double chi2 = (shp2.wireDist()-sth.wdist2())/shp2.posRes(StrawHitPosition::phi);
	    if(fabs(chi1) <_maxTDNsig && fabs(chi2) < _maxTDNsig)
	    {
	      // compute chisquared
	      double chisq = chi1*chi1+chi2*chi2; 
	      if(chisq < _maxChisq){
		sth.setChisquared(chisq);
		stereohits.push_back(sth);
// choose the best pair as:
// 1) take the pair with the minimum plane separation
// 2) otherwise, take the pair with the minimum chisquared
		if(sep < minsep[ish] || (sep == minsep[ish] && chisq < minchisq[ish])){
		  minsep[ish] = sep;
		  minchisq[ish] = chisq;
		  ibest[ish] = stereohits.size()-1;
		}
		if(sep < minsep[jsh] || (sep == minsep[jsh] && chisq < minchisq[jsh])){
		  minsep[jsh] = sep;
		  minchisq[jsh] = chisq;
		  ibest[jsh] = stereohits.size()-1;
		}
	      }
	    }
	  }
	}
      }
    }
// now, overwrite the positions for those hits which have stereosresolve the stereo hits to find the best position for each hit that particpates.  The algorithm is:
    for(size_t ish=0; ish<nsh;++ish){
      if(minsep[ish] < SectorId::apart)shpos->at(ish) = StrawHitPosition(stereohits,ibest[ish],ish);
      if(_diagLevel > 0){
	_sep->Fill(minsep[ish]);
	_chisq->Fill(minchisq[ish]);
      }
    }
    if(_writepairs){
      unique_ptr<StereoHitCollection> sthits(new StereoHitCollection(stereohits));
      event.put(std::move(sthits));
    }
    event.put(std::move(shpos));
  } // end MakeStereoHits::produce.

} // end namespace mu2e

using mu2e::MakeStereoHits;
DEFINE_ART_MODULE(MakeStereoHits)


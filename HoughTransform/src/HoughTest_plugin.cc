//
// An EDAnalyzer Module that runs the HoughTransform L-tracker code
//
// $Id: HoughTest_plugin.cc,v 1.2 2010/01/11 21:37:59 rhbob Exp $
// $Author: rhbob $ 
// $Date: 2010/01/11 21:37:59 $
//
// Original author R. Bernstein
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TSpectrum3.h"
#include "TMultiLayerPerceptron.h"
#include "TMLPAnalyzer.h"


// Mu2e includes.
#include "LTrackerGeom/inc/LTracker.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "HoughTransform/inc/HoughTransform.hh"
#include "HitCluster/inc/HitCluster.hh"
#include "GeneralUtilities/inc/RootNameTitleHelper.hh"
#include "GeneralUtilities/inc/pow.hh"

//CLHEP includes
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/Randomize.h"

using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::RandPoisson;
using namespace mu2e;
using namespace mu2e::houghtransform;
namespace mu2e {

  //--------------------------------------------------------------------
  //
  // 

    //the fitting function; these must be outside the class def'n or ROOT can't know about them (or I could make a static)
Double_t houghFitToRadius(Double_t *x, Double_t *par);

Double_t houghFitToRadius(Double_t *x, Double_t *par)
{
  Double_t background = par[0] + par[1]*x[0];
  //  Double_t peak = par[2] * (1./(TMath::Sqrt(2.*TMath::Pi())*par[4]) )* TMath::Exp(- 0.5*pow<2>((x[0] - par[3])/par[4]) );
  Double_t peak = par[2] * TMath::Exp(- 0.5*pow<2>((x[0] - par[3])/par[4]) );
  //  cout << "x, background, peak " << x[0] << " " << background <<" " << peak << endl;
  return background + peak; 
}

  class HoughTest : public edm::EDAnalyzer {
  public:
    explicit HoughTest(edm::ParameterSet const& pset) : 
      _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",10)),
      _nAnalyzed(0),
      _hRadius(0),
      _hTime(0),
      _hMultiplicity(0),
      _hDriftDist(0),
      _messageCategory("ToyHitInfo"){
    }
    virtual ~HoughTest() { }

    virtual void beginJob(edm::EventSetup const&);
    virtual void endJob();

    virtual void beginRun(edm::Run const &r, 
			  edm::EventSetup const& eSetup );

    virtual void beginLuminosityBlock(edm::LuminosityBlock const& lblock, 
				      edm::EventSetup const&);
 
    // This is called for each event.
    void analyze(const edm::Event& e, edm::EventSetup const&);

    RanecuEngine* inefficiencyEngine;
    RanecuEngine* accidentalEngine;
    RanecuEngine* regularEngine;


  private:

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Number of events analyzed.
    int _nAnalyzed;

    // Pointers to histograms to be filled.
    TH1F* _hRadius;
    TH1F* _hTime;
    TH1F* _hMultiplicity;
    TH1F* _hDriftDist;
    TH1F* _hxHit;
    TH1F* _hyHit;
    TH1F* _hzHit;
    TH1F* _hHitNeighbours;
    TH1F* _hCheckPointRadius;
    TH1F* _hRadiusDistribution2D;
    TH1F* _hDCADistribution2D;
    TH1F* _hChisqDistribution;
    TNtuple* _ntup;
    TF1* houghRadiusFit;
    TH2F* _hRadiusPeak;
    TH1F* _hNumberOfSpectrumPeaks;
    TH1F* _hNumberOfSpectrumPeaks2;

    // A category for the error logger.
    const std::string _messageCategory;

    // A helper function.
    int countHitNeighbours( Straw const& straw, 
			    StepPointMCCollection const* hits );

  };


  void HoughTest::beginJob(edm::EventSetup const& ){

    // Get access to the TFile service.
    edm::Service<edm::TFileService> tfs;

    // Create some 1D histograms.
    _hRadius       = tfs->make<TH1F>( "hRadius", "Radius of Hits;(mm)",          100,  0., 1000. );
    _hTime         = tfs->make<TH1F>( "hTime", "Pulse Height;(ns)",              100,  0.,  100. );
    _hMultiplicity = tfs->make<TH1F>( "hMultiplicity", "Hits per Event",         100,  0.,  100. );
    _hDriftDist    = tfs->make<TH1F>( "hDriftDist", "Crude Drift Distance;(mm)", 100,  0.,   3.  );
    _hxHit         = tfs->make<TH1F>( "hxHit",  "X of Hit;(mm)",                 
				      100,  -1000.,  1000. );
    _hyHit         = tfs->make<TH1F>( "hyHit",  "Y of Hit;(mm)",                 
				      100,  -1000.,  1000. );
    _hzHit         = tfs->make<TH1F>( "hzHit",  "Z of Hit;(mm)",                 
				      100,  -1400.,  1400. );

    _hHitNeighbours    = tfs->make<TH1F>( "hHitNeighbours",  "Number of hit neighbours",
					  10, 0., 10. );

    _hCheckPointRadius = tfs->make<TH1F>( "hCheckPointRadius",  "Radius of Reference point; (mm)",
					  100, 2.25, 2.75 );

    // Create an ntuple.
    _ntup           = tfs->make<TNtuple>( "ntup", "Hit ntuple", 
					  "evt:trk:sid:hx:hy:hz:wx:wy:wz:dca:time:dev:sec");
    _hRadiusDistribution2D = tfs->make<TH1F>("radiusDistribution2D","Radius Distribution in 2D search",100,0.,1000.);
    _hDCADistribution2D = tfs->make<TH1F>("DCADistribution2D","DCA Distribution in 2D Search",120,0.,1200.);

   //set up fit function in Hough Space
    houghRadiusFit = tfs->make<TF1>("houghRadiusFit",houghFitToRadius,0.,100.,5);

    //location of peak for each event (up to 100)
    _hRadiusPeak = tfs->make<TH2F>("radiusPeak","Radius Peak", 100,0.,100.,100,200.,400.);

    //and the number of peaks
    _hNumberOfSpectrumPeaks =  tfs->make<TH1F>("numberOfSpectrumPeaks" ,"Number of Radius Spectrum Peaks", 10,0.,10.);
    _hNumberOfSpectrumPeaks2 = tfs->make<TH1F>("numberOfSpectrumPeaks2","Number of Radius/DCA Spectrum Peaks", 10,0.,10.);

    //make a special random number generator for noise events
    accidentalEngine = new RanecuEngine();

    //make a special random number generator for inefficiency
    inefficiencyEngine = new RanecuEngine();

    //here's the regular engine
    regularEngine = new RanecuEngine();
    HepRandom::setTheEngine(regularEngine);

  }

  void HoughTest::endJob(){
  }



  void HoughTest::beginRun(edm::Run const& run,
                                 edm::EventSetup const& eSetup ){
  }

  void HoughTest::beginLuminosityBlock(edm::LuminosityBlock const& lblock,
					     edm::EventSetup const&){
  }


  void
  HoughTest::analyze(const edm::Event& evt, edm::EventSetup const&) {

    static int ncalls(0);
    ++ncalls;

    // Master geometry for the LTracker.
    GeomHandle<LTracker> ltracker;

    // Maintain a counter for number of events seen.
    ++_nAnalyzed;

    // Instance name of the module that created the hits of interest;
    static const string creatorName("g4run");


    // Ask the event to give us a "handle" to the requested hits.
    //    edm::Handle<StepPointMCCollection> hits;
    //evt.getByLabel(creatorName,hits);
    
    edm::Handle<StepPointMCCollection> hitsHandle;
    evt.getByLabel(creatorName,hitsHandle);
    StepPointMCCollection const* hits = hitsHandle.product();


    int nstraws = ltracker->getAllStraws().size();
    static double occupancy = 0.05;
    static int noiseMean = static_cast<int>(nstraws*occupancy);
    int numberOfNoiseHits = RandPoisson::shoot(accidentalEngine,noiseMean);
          numberOfNoiseHits = 0;

    //make a copy of the hits
    StepPointMCCollection plusNoise(*hitsHandle);
    plusNoise.reserve(plusNoise.size()+numberOfNoiseHits);
    //    cout << "size of plusNoise before " << plusNoise.size() << "and number of noiseHits =" << numberOfNoiseHits << endl;

    //add noise hits
    int trackIDnoise = 2;
    double eDepNoise = 0.;
    double timeNoise = 0.;
    Hep3Vector momentumNoise;
    HepRandom::setTheEngine(accidentalEngine);
    for (int i=0; i < numberOfNoiseHits; ++i)
      {
	//int istraw = static_cast<int>(nstraws*CLHEP::RandFlat::shoot());
	int istraw = static_cast<int>(nstraws*CLHEP::RandFlat::shoot(accidentalEngine));
	Straw const& straw = ltracker->getStraw( StrawIndex(istraw) );
	Hep3Vector mid = straw.getMidPoint(); //for the HT I only use x/y so the midpoint is perfect
	Hep3Vector w   = straw.getDirection();

		plusNoise.push_back(StepPointMC(trackIDnoise,istraw,eDepNoise,timeNoise,mid,momentumNoise));
      }
    //go back to normal engine
    HepRandom::setTheEngine(regularEngine);
    cout << "size of plusNoise after " << plusNoise.size() << endl;

    StepPointMCCollection const* plusNoiseHits = &plusNoise;
    //replace hits with plusNoise
    hits = &plusNoise;


    //remove hits according to inefficiency
    double inefficiency = 0.0;
    HepRandom::setTheEngine(inefficiencyEngine);

    //loop over all hits and delete the element

    //note to self, 1/5/10: I don't think I can just loop over hits and delete them because 
    //I'm going to foul up the iterator.  I think I need to loop over the original and then copy new elements
    //to it if they meet the criterion.  Rob says this is the right way for now.  Code it up post-Rice

    //behind the scenes somebody built me a copy constructor and assignment operators for vectors.  thanks!
    StepPointMCCollection inefficientHits(*hitsHandle);
    inefficientHits.clear();
    cout << "size of inefficientHits = " << inefficientHits.size() << endl;
    for( StepPointMCCollection::const_iterator i = hits->begin(); i!= hits->end(); ++i ) 
      {
	if (CLHEP::RandFlat::shoot(inefficiencyEngine) >= inefficiency)
	  {
	    // > inefficiency, add the element in inefficientHits
	    const StepPointMC& hit = *i;
	    inefficientHits.push_back(StepPointMC(hit));
	  }
      }
    cout << "size of inefficientHits after = " << inefficientHits.size() << endl;
    //and now make inefficientHits what gets passed to everyone else
    hits = &inefficientHits;

    // Fill histogram with number of hits per event.
    _hMultiplicity->Fill(hits->size());

 
    //    assert(2==1);

    // don't bother unless there are some minimum number of hits; start with 3. I might want this in the constructor as an argument...

    //so I have access to histos from results
    edm::Service<edm::TFileService> tfs;



	edm::LogVerbatim log(_messageCategory);
    log << "HoughTransforming event #: " 
	<< evt.id().event();

    cout << endl << "RHB HoughTransforming event # " << evt.id().event() << endl;

    int minHitsforHough = 3;
    cout << "minHitsforHough must be >3, size of hits = " << hits->size() << "and plusHits size = " << plusNoiseHits->size() << endl; 
    if (hits->size() > minHitsforHough)
      {
	mu2e::houghtransform::HoughTransform trialHough;
	mu2e::houghtransform::HoughTransform::houghCandidates houghTracks;

    // look for circles in hough space
	trialHough.foundHoughTracks(ltracker,plusNoiseHits,houghTracks);

	//loop over all found houghTracks, each of which is a candidate circle
	//and make a pair of histos that looks at accumulator space in radius and center

	
	RootNameTitleHelper radiusHisto("_histoRadius","Radius Accumulator Space Event",evt.id().event(),5);
	RootNameTitleHelper radiusErrorHisto("_histoErrorRadius","Radius Error Event",evt.id().event(),5);
	RootNameTitleHelper dcaHisto("_histoDCA","Distance of Closest Approch",evt.id().event(),5);
	RootNameTitleHelper centerHisto("_histoCenter","Center Accumulator Space Event",evt.id().event(),5);
	RootNameTitleHelper nStrawsHisto("_histoNStraws","Number of Straws in Circles, Event",evt.id().event(),5);
	RootNameTitleHelper radiusDCAHisto("_histoRadiusDCA","Radius vs. DCA, Event",evt.id().event(),5);
	RootNameTitleHelper radiusDCACenterHisto("_histoRadiusDCACenter","Radius vs. DCA vs Center, Event",evt.id().event(),5);

	TH1F* _histoRadius = tfs->make<TH1F>(radiusHisto.name(),radiusHisto.title(),
					     80,100.,500.);
	TH1F* _histoDCA = tfs->make<TH1F>(dcaHisto.name(),dcaHisto.title(),
					     110,0.,110.);
	TH2F* _histoRadiusDCA = tfs->make<TH2F>(radiusDCAHisto.name(),radiusDCAHisto.title(),
						80,100.,500.,120,-10.,110.);
	TH3F* _histoRadiusDCACenter = tfs->make<TH3F>(radiusDCACenterHisto.name(),radiusDCACenterHisto.title(),
						80,100.,500.,120,-10.,110.,100,0.,1000.);
	TH2F* _histoCenter = tfs->make<TH2F>(centerHisto.name(),centerHisto.title(), 
					     100,-1000.,1000.,100,-1000.,1000.);
	TH1F* _histoRadiusError = tfs->make<TH1F>(radiusErrorHisto.name(),radiusErrorHisto.title(), 
						  100,0.,100.);
	TH1F* _histoNStraws = tfs->make<TH1F>(nStrawsHisto.name(),nStrawsHisto.title(), 
						  100,0.,100.);

	cout << "number of Hough Tracks = " << houghTracks.size() << endl;

	for (int ithHough = 0; ithHough < houghTracks.size(); ++ithHough)
	  {
	    // only look at circles that have reasonable number of straws; three per passage-> 1 cluster, times three passages, 
	    //is enough to give three points (clusters) which is what you need to get a helix =9  
	    if (houghTracks[ithHough].numberOfStraws >=6 )
	      {
		Double_t distance = TMath::Sqrt(pow<2>(houghTracks[ithHough].x0) + pow<2>(houghTracks[ithHough].y0));
		//		cout << "distance = " << distance << endl;
		_histoRadius->Fill(static_cast<Double_t>(houghTracks[ithHough].radius),
				   static_cast<Double_t>(houghTracks[ithHough].numberOfStraws));
		_histoRadius->Fill(static_cast<Double_t>(houghTracks[ithHough].radius));
		_histoRadiusDCA->Fill(static_cast<Float_t>(houghTracks[ithHough].radius),
				      static_cast<Float_t>(houghTracks[ithHough].dca),
				      static_cast<Float_t>(houghTracks[ithHough].numberOfStraws));
		//		_histoRadiusDCA->Fill(250.,50.);
		_histoRadiusDCACenter->Fill(static_cast<Double_t>(houghTracks[ithHough].radius),
				      static_cast<Double_t>(houghTracks[ithHough].dca),
				      distance,
				      static_cast<Double_t>(houghTracks[ithHough].numberOfStraws));
		_histoDCA->Fill(houghTracks[ithHough].dca,houghTracks[ithHough].numberOfStraws);
		_histoCenter->Fill(houghTracks[ithHough].x0,houghTracks[ithHough].y0,houghTracks[ithHough].numberOfStraws);
	      }
	    _histoNStraws->Fill(houghTracks[ithHough].numberOfStraws);
	  }
		

	
	TSpectrum s(50);// calling this with (50) instead of () is a magic trick from Rene Brun. 
	s.Search(_histoRadius,1.0," ",0.50);
	Int_t nPeaks = s.GetNPeaks();
	_hNumberOfSpectrumPeaks->Fill(static_cast<double>(nPeaks));
	//TSpectrum2 doesn't work correctly, and examples in $ROOT_DIR don't either without turning 
	//off background and Markov smoothing; hence nobackgroundnomarkov. I don't think this matters
	//for us.  The examples in ROOT_DIR seem to work better but find a lot more ghost hits; with those
	//options it just finds random junk
	TSpectrum2 s2(50,1.0);  
	s2.Search(_histoRadiusDCA,2,"nobackgroundnomarkov",0.50);
	Int_t nPeaks2 = s2.GetNPeaks();
	cout << "nPeaks = " << nPeaks << endl << " and n2Peaks = " << nPeaks2 << endl;
	//make sure there's something to fit to!
	if (nPeaks > 0)
	  {
	    Float_t* firstPeak = s.GetPositionX();
	    for (int ithPeak = 0; ithPeak < nPeaks; ++ithPeak)
	      {
		cout <<"peak number " << ithPeak << " " << firstPeak[ithPeak] << endl;	
		_hRadiusPeak->Fill(static_cast<float>(ncalls) - 0.5,firstPeak[ithPeak]);
	      }
	  }
	_hNumberOfSpectrumPeaks2->Fill(static_cast<double>(nPeaks2));

	if (nPeaks2 > 0)
	  {
	    Float_t* firstPeakX = s2.GetPositionX();
	    Float_t* firstPeakY = s2.GetPositionY();
	    for (int ithPeak = 0; ithPeak < nPeaks2; ++ithPeak)
	      {
		cout <<"peak number " << ithPeak << " " << firstPeakX[ithPeak] << " " << firstPeakY[ithPeak] << endl;
		_hRadiusDistribution2D->Fill(firstPeakX[ithPeak]);
		_hDCADistribution2D->Fill(firstPeakY[ithPeak]);
	      }
	  }
	


      }
	
    
    
    // Loop over all plusNoiseHits.
    int n(0);
    StepPointMCCollection::const_iterator i = plusNoiseHits->begin();
    StepPointMCCollection::const_iterator e = plusNoiseHits->end();
    //scatter plot showing event
	RootNameTitleHelper eventPlot("_eventPlot","EventPlot",evt.id().event(),5);
	TH2F* _eventPlot= tfs->make<TH2F>(eventPlot.name(),eventPlot.title(), 
					     100,-1000.,1000.,100,-1000.,1000.);

	_eventPlot->SetMarkerStyle(20);
	_eventPlot->SetMarkerSize(0.5);
    for ( ; i!=e; ++i){
      
      // Aliases, used for readability.
      const StepPointMC& hit = *i;
      const Hep3Vector& pos = hit.position();
      const Hep3Vector& mom = hit.momentum();
      _eventPlot->Fill(pos[0],pos[1]);
      // Get the straw information.
      Straw const& straw = ltracker->getStraw( StrawIndex(hit.volumeId()) );
      //      cout << "straw, position of hit " << hit.volumeId() << " " << pos[0] << " " << pos[1] << endl;
      Hep3Vector mid = straw.getMidPoint();
      Hep3Vector w   = straw.getDirection();

      // Count how many nearest neighbours are also hit.
      int nNeighbours = countHitNeighbours( straw, plusNoiseHits );

      // Compute an estimate of the drift distance.
      TwoLinePCA pca( mid, w, pos, mom);

      // Check that the radius of the reference point in the local
      // coordinates of the straw.  Should be 2.5 mm.
      double s = w.dot(pos-mid);
      Hep3Vector point = pos - (mid + s*w);

      // I don't understand the distribution of the time variable.
      // I want it to be the time from the start of the spill.
      // It appears to be the time since start of tracking.

      // Fill some histograms
      _hRadius->Fill(pos.perp());
      _hTime->Fill(hit.time());
      _hHitNeighbours->Fill(nNeighbours);
      _hCheckPointRadius->Fill(point.mag());

      _hxHit->Fill(pos.x());
      _hyHit->Fill(pos.y());
      _hzHit->Fill(pos.z());

      _hDriftDist->Fill(pca.dca());

      /*
      // Fill the ntuple.
      nt[0]  = evt.id().event();
      nt[1]  = hit.trackId();
      nt[2]  = hit.volumeId();
      nt[3]  = pos.x();
      nt[4]  = pos.y();
      nt[5]  = pos.z();
      nt[6]  = mid.x();
      nt[7]  = mid.y();
      nt[8]  = mid.z();
      nt[9]  = pca.dca();
      nt[10] = hit.time();
      nt[11] = straw.Id().getDevice();
      nt[12] = straw.Id().getSector();

      _ntup->Fill(nt);
      */
      // Print out limited to the first few events.
      if ( ncalls < _maxFullPrint ){
	cout << "Readback hit: "
	     << evt.id().event()
	     << n++ <<  " "
	     << hit.trackId()    << "   "
	     << hit.volumeId() << " | "
	     << pca.dca()   << " "
	     << pos  << " "
	     << mom  << " "
	     << point.mag() 
	     << endl;
      }

    } // end loop over hits.
    
  } // end of ::analyze.

    //minuit fitting code



  // Count how many of this straw's nearest neighbours are hit.
  // If we have enough hits per event, it will make sense to make
  // a class to let us direct index into a list of which straws have hits.
  int HoughTest::countHitNeighbours( Straw const& straw, 
				    StepPointMCCollection const* hits ){
    
    int count(0);
    vector<StrawIndex> const& nearest = straw.nearestNeighboursByIndex();
    for ( vector<int>::size_type ihit =0;
	  ihit<nearest.size(); ++ihit ){

      StrawIndex idx = nearest[ihit];

      for( StepPointMCCollection::const_iterator 
	     i = hits->begin(),
	     e = hits->end(); i!=e ; ++i ) {
	const StepPointMC& hit = *i;
	if ( hit.volumeId() == idx.asInt() ){
	  ++count;
	  break;
	}
      }

    }
    return count;
  }
  
}


using mu2e::HoughTest;
DEFINE_FWK_MODULE(HoughTest);

//
// An EDProducer Module that runs the HoughTransform L-tracker code
//
// $Id: HoughTest_module.cc,v 1.11 2011/10/28 18:47:06 greenc Exp $
// $Author: greenc $
// $Date: 2011/10/28 18:47:06 $
//
// Original author R. Bernstein
//

#include "GeneralUtilities/inc/RootNameTitleHelper.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "HitCluster/inc/HitCluster.hh"
#include "HoughTransform/inc/HoughTransform.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMLPAnalyzer.h"
#include "TMath.h"
#include "TMultiLayerPerceptron.h"
#include "TNtuple.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TSpectrum3.h"
#include "RecoDataProducts/inc/HoughCircleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "cetlib/pow.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <iostream>
#include <string>

using CLHEP::Hep3Vector;
using cet::square;
using cet::sum_of_squares;

using namespace mu2e::houghtransform;
using namespace mu2e;
using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //

    //the fitting function; these must be outside the class def'n or ROOT can't know about them (or I could make a static)
Double_t houghFitToRadius(Double_t *x, Double_t *par);

Double_t houghFitToRadius(Double_t *x, Double_t *par)
{
  Double_t background = par[0] + par[1]*x[0];
  //  Double_t peak = par[2] * (1./(TMath::Sqrt(2.*TMath::Pi())*par[4]) )* TMath::Exp(- 0.5*square((x[0] - par[3])/par[4]) );
  Double_t peak = par[2] * TMath::Exp(- 0.5*square((x[0] - par[3])/par[4]) );
  return background + peak;
}
  static const bool singleHoughHisto = false;
  class HoughTest : public art::EDProducer {
  public:
    explicit HoughTest(fhicl::ParameterSet const& pset) :
      _maxFullPrint(pset.get<int>("maxFullPrint",10)),
      _nPeakSearch(pset.get<unsigned>("NPeakSearch")),
      _nAnalyzed(0),
      _hRadius(0),
      _hTime(0),
      _hMultiplicity(0),
      _hDriftDist(0),
      _messageCategory("HitInfo"),
      _hitCreatorName(pset.get<string>("hitCreatorName")),
      _useStepPointMC(pset.get<bool>("UseMCHits"))
    {
        produces<HoughCircleCollection>();
    }
    virtual ~HoughTest() { }

    virtual void beginJob();
    virtual void endJob();

    virtual void beginRun(art::Run const &r);

    virtual void beginSubRun(art::SubRun const& lblock);

    // This is called for each event.
    void produce(art::Event& e);


  private:

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;
    // number of peaks to look for
    unsigned _nPeakSearch;

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
                      // radius of center from center hough space -
                      // radius of center from radius vc. dca hough space
    TH1F* _hRadiusPeakEqual; // center and radius/dca peaks are equal in
                               // rank
    TH1F* _hRadiusPeakNNPlus; // center space peak is one rank higher than
                              // radius/dca space peak
    TH1F* _hRadiusPeakNNMinus; // center space peak is one rank lower than
                              // radius/dca space peak
    TH1F* _hRadiusPeakNotNN;// non of the above

    TH1F* _hRHitFromHTCenterMC; // radius of hits from peak hough center (MC hits)
    TH1F* _hRHitFromCircleMC; // residual of hits from peak hough track
    TH1F* _hRNoiseHitFromHTCenterMC; // radius of noise hits from peak hough center
    TH1F* _hRNoiseHitFromCircleMC; // residual of noise hits from peak hough track
    TH1F* _hRPhysHitFromHTCenterMC; // radius of physical hits from peak hough center
    TH1F* _hRPhysHitFromCircleMC; // residual of physical hits from peak hough track
    TH2F* _hRadius2v1; // radius from 2D vs. radius from 1D
    TH2F* _hYvsX; // Y vs. X of 1st peak
    TH1F* _hChisqDistribution;
    TNtuple* _ntup;
    TF1* houghRadiusFit;
    TH2F* _hRadiusPeak;
    TH1F* _hNumberOfSpectrumPeaks;
    TH1F* _hNumberOfSpectrumPeaks2;


    TH1F* _histoRadius;
    TH1F* _histoDCA;
    TH2F* _histoRadiusDCA;
    TH3F* _histoRadiusDCACenter;
    TH2F* _histoCenter;
    TH2F* _histoRXCenter;
    TH2F* _histoRYCenter;
    TH1F* _histoRadiusError;
    TH1F* _histoNStraws;
    TH2F* _eventPlot;

    // A category for the error logger.
    const std::string _messageCategory;

    //name of the module that created the hits to be used
    const std::string _hitCreatorName;

    bool _useStepPointMC; // use MC Hits

    // A helper function.
    int countHitNeighbours( Straw const& straw,
                            StepPointMCCollection const* hits );

    void bookEventHistos(art::EventNumber_t);
    void fillEventHistos(
         mu2e::houghtransform::HoughTransform::houghCircleStruct);
  };


  void HoughTest::beginJob(){

    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

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

    _hRHitFromHTCenterMC = tfs->make<TH1F>( "hRHitFromHTCenterMC",
                         "Radius of (MC) Hits from Peak Hough Center",100,0,400);

    _hRHitFromCircleMC = tfs->make<TH1F>( "hRHitFromCircleMC",
                         "Radius (MC) hits - radius Peak",100,-100,100);

    _hRNoiseHitFromHTCenterMC = tfs->make<TH1F>( "hRNoiseHitFromHTCenterMC",
                     "Radius of Noise (MC) Hits from Peak Hough Center",100,0,400);

    _hRNoiseHitFromCircleMC = tfs->make<TH1F>( "hRNoiseHitFromCircleMC",
                     "Radius Noise (MC) hits - radius Peak",100,-100,100);

    _hRPhysHitFromHTCenterMC = tfs->make<TH1F>( "hRPhysHitFromHTCenterMC",
                   "Radius of Physical (MC) Hits from Peak Hough Center",100,0,400);

    _hRPhysHitFromCircleMC = tfs->make<TH1F>( "hRPhysHitFromCircleMC",
                   "Radius Physical (MC) hits - radius Peak",100,-100,100);

    // Create an ntuple.
    _ntup           = tfs->make<TNtuple>( "ntup", "Hit ntuple",
                                          "evt:trk:sid:hx:hy:hz:wx:wy:wz:dca:time:dev:sec");
    _hRadiusDistribution2D = tfs->make<TH1F>("radiusDistribution2D","Radius Distribution in 2D search",100,0.,1000.);
    _hDCADistribution2D = tfs->make<TH1F>("DCADistribution2D","DCA Distribution in 2D Search",120,0.,1200.);
    _hRadius2v1 = tfs->make<TH2F>("hRadius2v1","First Radius from 2D vs. 1D ",
        100,0.,400.,100,0,400.);
    _hYvsX = tfs->make<TH2F>("hYvsX","First Radius from 2D vs. 1D ",
        100,-500,500,100,-500,500);

    double deltaRmin=-140;
    double deltaRmax=60;
    _hRadiusPeakEqual = tfs->make<TH1F>("hRadiusPeakEqual",
           "Delta R for center, center and rad/dca peaks equal rank",
           100,deltaRmin,deltaRmax);
    _hRadiusPeakNNPlus = tfs->make<TH1F>("hRadiusPeakNNPlus",
           "Delta R for center, center peak one rank higher than rad/dca",
           100,deltaRmin,deltaRmax);
    _hRadiusPeakNNMinus = tfs->make<TH1F>("hRadiusPeakNNMinus",
           "Delta R for center, center peak one rank lower than rad/dca",
           100,deltaRmin,deltaRmax);
    _hRadiusPeakNotNN = tfs->make<TH1F>("hRadiusPeakNotNN",
           "Delta R for center, center peak and rad/dca peak not neighbors in rank",
           100,deltaRmin,deltaRmax);

;

   //set up fit function in Hough Space
    houghRadiusFit = tfs->make<TF1>("houghRadiusFit",houghFitToRadius,0.,100.,5);

    //location of peak for each event (up to 100)
    _hRadiusPeak = tfs->make<TH2F>("radiusPeak","Radius Peak", 100,0.,100.,100,200.,400.);

    //and the number of peaks
    _hNumberOfSpectrumPeaks =  tfs->make<TH1F>("numberOfSpectrumPeaks" ,"Number of Radius Spectrum Peaks", 10,0.,10.);
    _hNumberOfSpectrumPeaks2 = tfs->make<TH1F>("numberOfSpectrumPeaks2","Number of Radius/DCA Spectrum Peaks", 10,0.,10.);

  }

  void HoughTest::endJob(){
  }



  void HoughTest::beginRun(art::Run const& run){
  }

  void HoughTest::beginSubRun(art::SubRun const& lblock){
  }


  void HoughTest::produce(art::Event& evt) {

    static int ncalls(0);
    ++ncalls;

    // Maintain a counter for number of events seen.
    ++_nAnalyzed;

    // the product
    auto_ptr<HoughCircleCollection> HoughResults(new HoughCircleCollection);


    // Ask the event to give us a handle to the requested hits.
    //    art::Handle<StepPointMCCollection> hits;
    //evt.getByLabel(creatorName,hits);
    art::Handle<StepPointMCCollection> hitsHandle;
    static const string collectionName("tracker");

    // someday hitClusters will be produced elsewhere, but for now,
    // we get them from StepPointMCs

    vector<mu2e::hitcluster::HitCluster> hitClusters;

    if (_useStepPointMC) {
      evt.getByLabel(_hitCreatorName,collectionName,hitsHandle);
       StepPointMCCollection const* hits = hitsHandle.product();
       if (!hits->size()) return; // return with empty product
       HoughTransform::MakeClusters(hits,hitClusters);
       if (!hitClusters.size()) return;
    }// else get them from somewhere else when the means to do so exists {}

    // Master geometry for the LTracker.
    GeomHandle<LTracker> ltracker;

    // a helper class for lots of the heavy lifting
    HoughTransform houghHelper(hitClusters);


    // double InitialRadius=-1;
/*
    // if were using a radius from a first pass, get it here
    if (_InitialRadiusCreater!="none") {
       art::Handle<HoughCircleCollection> hcHandle;
       evt.getByLabel(_InitialRadiusCreator,hcHandle);
       if (hcHandle->size()) {
          const HoughCircle& hc=hcHandle->at(0);
          InitialRadius=hc.Radius();
       } else {
          std::cout<<"Looking for initial radius, but no hough circles found."
            << std::endl;
       }
    }
*/

    // Book histogram on the first call regardless, and book new ones if
    //  unless disabled by singHoughHisto
    if ( ncalls == 1 || !singleHoughHisto) bookEventHistos(evt.id().event());

/*
    do 1st loop

    fill histograms ? or part of 1st loop?

      (set radius
      do 2nd loop
      fill histograms)

    make data products
*/


/*
    // if were using a radius from a first pass, get it here
    if (_InitialRadiusCreater!="none") {
       art::Handle<HoughCircleCollection> hcHandle;
       evt.getByLabel(_InitialRadiusCreator,hcHandle);
       if (hcHandle->size()) {
          const HoughCircle& hc=hcHandle->at(0);
          InitialRadius=hc.Radius();
       } else {
          std::cout<<"Looking for initial radius, but no hough circles found."
            << std::endl;
       }
    }
*/

    // Fill histogram with number of hits per event.
    _hMultiplicity->Fill(hitsHandle->size());


    // don't bother unless there are some minimum number of hits; start with 3. I might want this in the constructor as an argument...

    //so I have access to histos from results
    // used here? art::ServiceHandle<art::TFileService> tfs;



        mf::LogVerbatim log(_messageCategory);
    log << "HoughTransforming event #: "
        << evt.id().event();

    //cout << endl << "RHB HoughTransforming event # " << evt.id().event() << endl;



    size_t minHitsforHough = 3;
    if (hitClusters.size() > minHitsforHough)
      {
        mu2e::houghtransform::HoughTransform::houghCandidates houghTracks;

    // look for circles in hough space: findHoughTracks fills houghTracks
    // with all houghCircleStucts it can find.  The latter are all combinations
    // of 3 hits with their exact circle solution.
        houghHelper.foundHoughTracks(ltracker,houghTracks);

        cout << "ncalls = " << ncalls << endl;

        //loop over all found houghTracks, each of which is a candidate circle
        // and make a pair of histos that looks at accumulator space in radius
        // and center

        for (size_t ithHough = 0; ithHough < houghTracks.size(); ++ithHough)
          {
            // only look at circles that have reasonable number of straws; three per passage-> 1 cluster, times three passages,
            //is enough to give three points (clusters) which is what you need to get a helix =9
                 // fill the histos with the found circles
                  fillEventHistos(houghTracks[ithHough]);
           }



       // Look for peaks in the hough space
        TSpectrum s(_nPeakSearch);// calling this with (_nPeakSearch) instead of () is a magic trick from Rene Brun.
        s.Search(_histoRadius,1.0," ",0.95);
        Int_t nPeaks = s.GetNPeaks();
        _hNumberOfSpectrumPeaks->Fill(static_cast<double>(nPeaks));
        //TSpectrum2 doesn't work correctly, and examples in $ROOT_DIR don't either without turning
        //off background and Markov smoothing; hence nobackgroundnomarkov. I don't think this matters
        //for us.  The examples in ROOT_DIR seem to work better but find a lot more ghost hits; with those

        //options it just finds random junk
        TSpectrum2 s2(_nPeakSearch,1.0);
        s2.Search(_histoRadiusDCA,2,"nobackgroundnomarkov",0.95);
        Int_t nPeaks2 = s2.GetNPeaks();

        TSpectrum sCenterX(_nPeakSearch); // gets a histogram filled after we
        TSpectrum sCenterY(_nPeakSearch); // gets a histogram filled after we
                       // know R

        Float_t* rad1D=0; // radius from 1D plot
        Float_t* rad2D=0; // radius from 2D plot
        Float_t* dca2D=0; // Distance of closest approach from 2D plot
        Float_t* centX=0;  // center in X and Y
        Float_t* centY=0;

        //make sure there's something to fit to!
        if (nPeaks > 0)
          {
            rad1D = s.GetPositionX();
            for (int ithPeak = 0; ithPeak < nPeaks; ++ithPeak)
              {
                //                cout <<"peak number " << ithPeak << " " << firstPeak[ithPeak] << endl;
                _hRadiusPeak->Fill(static_cast<float>(ncalls) - 0.5,rad1D[ithPeak]);
              }
          }
        _hNumberOfSpectrumPeaks2->Fill(static_cast<double>(nPeaks2));

        if (nPeaks2 > 0)
          {
            rad2D = s2.GetPositionX();
            dca2D = s2.GetPositionY();
            for (int ithPeak = 0; ithPeak < nPeaks2; ++ithPeak)
              {
                //                cout <<"peak number " << ithPeak << " " << firstPeakX[ithPeak] << " " << firstPeakY[ithPeak] << endl;
                _hRadiusDistribution2D->Fill(rad2D[ithPeak]);
                _hDCADistribution2D->Fill(dca2D[ithPeak]);
              }

              if (nPeaks>0) _hRadius2v1->Fill(rad1D[0],rad2D[0]);

              // look for center in X in a band of best peak radius
              // Radius already found.  Same for Y
              int bin=_histoRXCenter->GetYaxis()->FindBin(rad2D[0]);
              int binl=max(bin-1,1);
              int binu=min(bin+1,_histoRXCenter->GetNbinsY());

              TH1D* xProj=_histoRXCenter->ProjectionX(Form("_px"),binl,binu);
              sCenterX.Search(xProj,1," ",0.95);

              TH1D* yProj=_histoRYCenter->ProjectionX(Form("_py"),binl,binu);
              sCenterY.Search(yProj,1," ",0.95);

          }//nPeaks2>0

          Int_t nPeaksCenterX = sCenterX.GetNPeaks();
          Int_t nPeaksCenterY = sCenterY.GetNPeaks();

          if (nPeaksCenterX > 0  && nPeaksCenterY >0 )
          {
               centX = sCenterX.GetPositionX();
               centY = sCenterY.GetPositionX();
              _hYvsX->Fill(centX[0],centY[0]);

          }


          // loop over all combinations of the peaks
          int nPeaksCenter=min(nPeaksCenterX,nPeaksCenterY);
          for (int ipc=0; ipc<nPeaksCenter; ipc++)  {

            // radius of the hough center (from detector axis) calculated
            // from center space peak
            double radial_center=sqrt(TMath::Power(centX[ipc],2)
                                     +TMath::Power(centY[ipc],2));

            for (int ipr=0; ipr<nPeaks2; ipr++) {

              // radius of the hough center (from detector axis) calculated
              // from radius-DCA space peak
           // this is useless - it has an ambiguity of whether the circle
           // encompases the origin
              double radial_radiusDCA=dca2D[ipr]+rad2D[ipr];

              // plots for matching peaks in the two different spaces
              TH1F* hrad=0;
              if (ipr==ipc) {
                 hrad=_hRadiusPeakEqual;
              } else if (ipr==ipc+1) {
                 hrad=_hRadiusPeakNNPlus;
              } else if (ipr==ipc-1) {
                 hrad=_hRadiusPeakNNMinus;
              } else {
                 hrad=_hRadiusPeakNotNN;
              }

              hrad->Fill(radial_center-radial_radiusDCA);

              // for now, "the" hough track is the one defined by
              // the center of the 1st peak in center space and the radius
              // of the 1st peak in radius-DCA space.
              if (ipc==0 && ipr==0) {


                 // create the data product
                        // currently "number of straws"=1.  This is a dummy,
                        // soon to be replaced with useful info
                 HoughCircle hc(centX[ipc],centY[ipc],rad2D[ipr],1);
                 HoughResults->push_back(hc);

                 if (_useStepPointMC) {

                    StepPointMCCollection const* hits = hitsHandle.product();
                    StepPointMCCollection::const_iterator ihit=hits->begin();
                    StepPointMCCollection::const_iterator endhit=hits->end();


                    for ( ; ihit!=endhit; ++ihit){

                      // Aliases, used for readability.
                       const StepPointMC& hit = *ihit;
                       const CLHEP::Hep3Vector& pos = hit.position();
                       double radFromCenter=sqrt(TMath::Power(centX[ipc]-pos.x(),2)
                                          +TMath::Power(centY[ipc]-pos.y(),2));
                       _hRHitFromHTCenterMC->Fill(radFromCenter);
                       _hRHitFromCircleMC->Fill(radFromCenter-rad2D[ipr]);
                       if (hit.trackId().asInt() == 2) {
                          _hRNoiseHitFromHTCenterMC->Fill(radFromCenter);
                          _hRNoiseHitFromCircleMC->Fill(radFromCenter-rad2D[ipr]);
                       } else {
                          _hRPhysHitFromHTCenterMC->Fill(radFromCenter);
                          _hRPhysHitFromCircleMC->Fill(radFromCenter-rad2D[ipr]);
                       }


                    } //MC hits

                 } //using MC hits


              } // ipc=ipr=0

            } // ipr
          } // ipc

      } // minHitsforHough



    // Loop over all hits
    int n(0);

    if (!singleHoughHisto)
      {
        art::ServiceHandle<art::TFileService> tfs;
        //scatter plot showing event
        RootNameTitleHelper eventPlot("_eventPlot","EventPlot",evt.id().event(),5);
        _eventPlot= tfs->make<TH2F>(eventPlot.name(),eventPlot.title(),
                                          100,-1000.,1000.,100,-1000.,1000.);

        _eventPlot->SetMarkerStyle(20);
        _eventPlot->SetMarkerSize(0.5);
      }
    // clear out multiple use histos
    if (singleHoughHisto)
      {
        _histoRadius->Reset();
        _histoDCA->Reset();
        _histoRadiusDCA->Reset();
        _histoRadiusDCACenter->Reset();
        _histoCenter->Reset();
      }


    if (_useStepPointMC) {
     // most of this needs to be switched over to Clusters

      StepPointMCCollection const* hits = hitsHandle.product();
      StepPointMCCollection::const_iterator ihit = hits->begin();
      const StepPointMCCollection::const_iterator endhit = hits->end();
      for ( ; ihit!=endhit; ++ihit){

        // Aliases, used for readability.
        const StepPointMC& hit = *ihit;
        const CLHEP::Hep3Vector& pos = hit.position();
        const CLHEP::Hep3Vector& mom = hit.momentum();
        if (!singleHoughHisto)
          {
            _eventPlot->Fill(pos[0],pos[1]);
          }
        // Get the straw information.
        Straw const& straw = ltracker->getStraw( StrawIndex(hit.volumeId()) );
        //      cout << "straw, position of hit " << hit.volumeId() << " " << pos[0] << " " << pos[1] << endl;
        CLHEP::Hep3Vector mid = straw.getMidPoint();
        CLHEP::Hep3Vector w   = straw.getDirection();

        // Count how many nearest neighbours are also hit.
        int nNeighbours = countHitNeighbours( straw, hits );

        // Compute an estimate of the drift distance.
        TwoLinePCA pca( mid, w, pos, mom);

        // Check that the radius of the reference point in the local
        // coordinates of the straw.  Should be 2.5 mm.
        double s = w.dot(pos-mid);
        CLHEP::Hep3Vector point = pos - (mid + s*w);

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
        nt[1]  = hit.trackId().asInt();
        nt[2]  = hit.volumeId();
        nt[3]  = pos.x();
        nt[4]  = pos.y();
        nt[5]  = pos.z();
        nt[6]  = mid.x();
        nt[7]  = mid.y();
        nt[8]  = mid.z();
        nt[9]  = pca.dca();
        nt[10] = hit.time();
        nt[11] = straw.id().getDevice();
        nt[12] = straw.id().getSector();

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
    } // if using MC hits

    evt.put(HoughResults);

  } // end of ::produce.

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
        if ( hit.volumeId() == idx.asUint() ){
          ++count;
          break;
        }
      }

    }
    return count;
  } //countHitNeighbors

  void HoughTest::bookEventHistos(art::EventNumber_t evtno)
  {
        RootNameTitleHelper radiusHisto("_histoRadius","Radius Accumulator Space Event",evtno,5);
        RootNameTitleHelper radiusErrorHisto("_histoErrorRadius","Radius Error Event",evtno,5);
        RootNameTitleHelper dcaHisto("_histoDCA","Distance of Closest Approch",evtno,5);
        RootNameTitleHelper centerHisto("_histoCenter","Center Accumulator Space Event",evtno,5);
        RootNameTitleHelper rXCenterHisto("_histoRXCenter","Radius vs. X Center Accumulator Space Event",evtno,5);
        RootNameTitleHelper rYCenterHisto("_histoRYCenter","Radius vs. Y Center Accumulator Space Event",evtno,5);
        RootNameTitleHelper nStrawsHisto("_histoNStraws","Number of Straws in Circles, Event",evtno,5);
        RootNameTitleHelper radiusDCAHisto("_histoRadiusDCA","Radius vs. DCA, Event",evtno,5);
        RootNameTitleHelper radiusDCACenterHisto("_histoRadiusDCACenter","Radius vs. DCA vs Center, Event",evtno,5);

        art::ServiceHandle<art::TFileService> tfs;

        _histoRadius = tfs->make<TH1F>(radiusHisto.name(),radiusHisto.title(),
                                       80,100.,500.);
        _histoDCA = tfs->make<TH1F>(dcaHisto.name(),dcaHisto.title(),
                                    110,0.,110.);
        _histoRadiusDCA = tfs->make<TH2F>(radiusDCAHisto.name(),radiusDCAHisto.title(),
                                          80,100.,500.,120,-10.,110.);
        _histoRadiusDCACenter = tfs->make<TH3F>(radiusDCACenterHisto.name(),radiusDCACenterHisto.title(),
                                                80,100.,500.,120,-10.,110.,100,0.,1000.);
        _histoCenter = tfs->make<TH2F>(centerHisto.name(),centerHisto.title(),
                                       100,-1000.,1000.,100,-1000.,1000.);
        _histoRXCenter = tfs->make<TH2F>(rXCenterHisto.name(),
                                         rXCenterHisto.title(),
                                       100,-1000.,1000.,80,100.,500.);
        _histoRYCenter = tfs->make<TH2F>(rYCenterHisto.name(),
                                         rYCenterHisto.title(),
                                       100,-1000.,1000.,80,100.,500.);
        _histoRadiusError = tfs->make<TH1F>(radiusErrorHisto.name(),radiusErrorHisto.title(),
                                            100,0.,100.);
        _histoNStraws = tfs->make<TH1F>(nStrawsHisto.name(),nStrawsHisto.title(),
                                        100,0.,100.);
  } //bookEventHistos

 void HoughTest::fillEventHistos(
         mu2e::houghtransform::HoughTransform::houghCircleStruct houghCircle)
 {
    _histoNStraws->Fill(houghCircle.numberOfStraws);
    // if there are enough straws, fill the other histograms
    if (houghCircle.numberOfStraws >=6 )
    {
       Double_t distance = TMath::Sqrt(sum_of_squares(houghCircle.x0, houghCircle.y0));
       _histoRadius->Fill(static_cast<Double_t>(houghCircle.radius));
       _histoRadiusDCA->Fill(static_cast<Float_t>(houghCircle.radius),
                  static_cast<Float_t>(houghCircle.dca),
                  static_cast<Float_t>(houghCircle.numberOfStraws));
       _histoRadiusDCACenter->Fill(static_cast<Double_t>(houghCircle.radius),
                  static_cast<Double_t>(houghCircle.dca), distance,
                  static_cast<Double_t>(houghCircle.numberOfStraws));
       _histoDCA->Fill(houghCircle.dca,houghCircle.numberOfStraws);
       _histoCenter->Fill(houghCircle.x0,houghCircle.y0,houghCircle.numberOfStraws);
       _histoRXCenter->Fill(houghCircle.x0,houghCircle.radius,houghCircle.numberOfStraws);
       _histoRYCenter->Fill(houghCircle.y0,houghCircle.radius,houghCircle.numberOfStraws);
    } // if enough straws

 } // fillEventHistos()

} // namespace HoughTest


using mu2e::HoughTest;
DEFINE_ART_MODULE(HoughTest);

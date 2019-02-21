// Straight Track finder- uses linear least squares 
// The purpose of this module is to fit to field-off cosmic tracks. These shall be utilised for the internal alignment of the trackers.
// $Id: StraightTrackFinder_module.cc
// $Author: S.middleton
// $Date: Nov. 2018 $


// ART:
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "GeneralUtilities/inc/Angles.hh"
#include "art/Utilities/make_tool.h"
//MU2E:
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/StraightTrackSeed.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "TrkReco/inc/StraightTrackFit.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TrkReco/inc/TrkTimeCalculator.hh"

//utils:
#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
//Track Fitting Headers:
#include "TrkPatRec/inc/StraightTrackFinder_types.hh"
#include "TrkReco/inc/StraightTrackFinderData.hh"
#include "TrkReco/inc/TrkFaceData.hh"

//CLHEP:

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"


//ROOT: 
#include "TH1F.h"
#include "TH2F.h"
#include "Math/VectorUtil.h"
#include "TVector2.h"

//C++:

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <utility>
#include <functional>
#include <float.h>
#include <vector>
#include <map>

using namespace std;
using namespace ROOT::Math::VectorUtil;

namespace{
    //create a compare struct to allow height ordering of the hits in an event. The highest hit will then seed the "cosmic" track....
    struct ycomp : public std::binary_function<mu2e::ComboHit,mu2e::ComboHit,bool> {
    bool operator()(mu2e::ComboHit const& p1, mu2e::ComboHit const& p2) { return p1._pos.y() > p2._pos.y(); }
  };//end ycomp

    struct zcomp : public std::binary_function<mu2e::ComboHit,mu2e::ComboHit,bool> {
    bool operator()(mu2e::ComboHit const& p1, mu2e::ComboHit const& p2) { return p1._pos.z() < p2._pos.z(); }
  }; //end zcomp

}//end namespace

namespace mu2e{

  class StraightTrackFinder : public art::EDProducer {
  public:
    explicit StraightTrackFinder(fhicl::ParameterSet const&);
    virtual ~StraightTrackFinder();
    virtual void beginJob();
    virtual void beginRun(art::Run& run);
    virtual void produce(art::Event& event );
  private:
  //config parameters:
    int				        n_TC_fails;
    int 				_diag,_debug;
    int                                 _printfreq;
    int 				_minnsh; // minimum # of strawHits to work wit
    TrkFitFlag				_saveflag;//write tracks that satisfy these flags
    unsigned				_maxniter;  // maximum # of iterations over outlier filtering + fitting TODO:USE This
    int 				_minNHitsTimeCluster; //min number of hits in a viable time cluster
    
    //Collection lists:
    art::ProductToken<ComboHitCollection> const _chToken;
    art::ProductToken<TimeClusterCollection> const _tcToken;
    

    //Define straw hit falgs
    StrawHitFlag  _hsel, _hbkg;
    TH1F* _niter, _chi2dXY;;
    //Define Track Fits
    StraightTrackFit     _tfit;
    TrkTimeCalculator _ttcalc;
    float             _t0shift; //TODO: use  
    StrawHitFlag      _outlier;

    std::unique_ptr<ModuleHistToolBase>   _hmanager;
    StraightTrackFinderTypes::Data_t      _data;
    StraightTrackFinderData               _stResult;
    
    
    void     OrderHitsY(StraightTrackFinderData& TrackData); //Order in height
    void     fillGoodHits(StraightTrackFinderData& TrackData);//apply "good" cut
    int     goodHitsTimeCluster(const TimeCluster TCluster, ComboHitCollection chcol);//select "good" time clusters
    
 
};//end private

//Constructor:
 StraightTrackFinder::StraightTrackFinder(fhicl::ParameterSet const& pset) :
    n_TC_fails  (pset.get<int>("n_TC_fails", 0)),
    _diag        (pset.get<int>("diagLevel",0)),
    _debug       (pset.get<int>("debugLevel",0)),
    _printfreq   (pset.get<int>   ("printFrequency", 101)),
    _minnsh      (pset.get<int>("minNStrawHits",4)),
    _saveflag    (pset.get<vector<string> >("SaveTrackFlag",vector<string>{"StraightTrackOK"})),
    _maxniter    (pset.get<unsigned>("MaxIterations",10)), // iterations over outlier removal
   _minNHitsTimeCluster(pset.get<int>("minNHitsTimeCluster", 1 )), 
   
    
    _chToken{consumes<ComboHitCollection>(pset.get<art::InputTag>("ComboHitCollection"))},
    _tcToken{consumes<TimeClusterCollection>(pset.get<art::InputTag>("TimeClusterCollection"))},
_hsel        (pset.get<std::vector<std::string> >("HitSelectionBits",std::vector<string>{"TimeDivision"})),
    _hbkg        (pset.get<std::vector<std::string> >("HitBackgroundBits",std::vector<std::string>{"Background"})),
    _tfit        (pset.get<fhicl::ParameterSet>("StraightTrackFit",fhicl::ParameterSet())), 
    _ttcalc      (pset.get<fhicl::ParameterSet>("T0Calculator",fhicl::ParameterSet())),
    _t0shift     (pset.get<float>("T0Shift",4.0)),//TODO:use
    _outlier     (StrawHitFlag::outlier)//, //TODO:check
   {
    produces<StraightTrackSeedCollection>();
    if (_debug != 0) _printfreq = 1;
    else    _hmanager = std::make_unique<ModuleHistToolBase>();
 }

//Destructor:
 StraightTrackFinder::~StraightTrackFinder(){}
 


/* ------------------------Begin JOb--------------------------//
//         Sets Up Historgram Book For Diag Plots            //
//----------------------------------------------------------*/

  void StraightTrackFinder::beginJob() {
    if (_debug != 0) std::cout<<"Beginning STFinder Job..."<<std::endl;
    art::ServiceHandle<art::TFileService> tfs;

    if (_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _niter = tfs->make<TH1F>( "niter" , "Number of Fit in XY Iteraions",201,-0.5,200.5);
      
      _hmanager->bookHistograms(tfs);
    }
   
  }
/* ------------------------Begin Run--------------------------//
//                   sets up the tracker                     //
//----------------------------------------------------------*/
  void StraightTrackFinder::beginRun(art::Run& ) {
    
// Sets up tracker information....adds this to the straighttrackfinder...
    if (_debug > 0)
    {
      std::cout << "BEGINNING STRAIGHT TRACK FINDING... " << std::endl;
    }
   mu2e::GeomHandle<mu2e::TTracker> th;
   const TTracker* tracker = th.get();
   _tfit.setTracker  (tracker);
   }


/* ------------------------PRODUCE--------------------------//
//                  Should produce ST seed Collection      //
//----------------------------------------------------------*/
  void StraightTrackFinder::produce(art::Event& event ) {
     if (_debug != 0) std::cout<<"Producing ST in  Finder..."<<std::endl;
     //Create output collection - this is collection of straght track seeds:
     unique_ptr<StraightTrackSeedCollection> seed_col(new StraightTrackSeedCollection());
     unsigned    STCounter(0);
     

     int _iev=event.id().event();
      if (_debug > 0){
          std::cout<<"ST Finder Event #"<<_iev<<std::endl;
      } //end "if"
      
     //fill in collection lists
     
     auto const& chH = event.getValidHandle(_chToken);
     const ComboHitCollection& chcol(*chH);

     auto const& tcH = event.getValidHandle(_tcToken);
     const TimeClusterCollection& tccol(*tcH);
     
     
    _data.event       = &event;
    _data.nseeds[0] = 0;
    
    
    //ensure all data types present. Fill data...
   
    _data.nTimePeaks  = tccol.size();
    _stResult._chcol  = &chcol; 
    _stResult._tccol  = &tccol;
   
    std::cout<<"Number of Time Peaks...: "<<tccol.size()<<std::endl;
    // create initial tracks from time clusters:

    for (size_t index=0;index< tccol.size();++index) {
      if(tccol.size()==0)		continue;
      int   nGoodTClusterHits(0);
      //get the t cluster
      const auto& tclust = tccol[index];
      // ensure only good clusters will be used for fit:
      nGoodTClusterHits     = goodHitsTimeCluster(tclust,chcol );
      if ( nGoodTClusterHits < _minNHitsTimeCluster)         continue;
      
      // create a place for the seed:
       StraightTrackSeed tseed ;

      //clear the variables in _stResult
      _stResult.clearTempVariables();

      //set variables used for searching the track candidate. pass these to the straighttrackfinderdata result (stResult)
      _stResult._tseed              = tseed;
      _stResult._timeCluster        = &tclust;
      _stResult._tseed._thits.setParent(chcol.parent());
      _stResult._tseed._t0          = tclust._t0;
      _stResult._tseed._timeCluster = art::Ptr<TimeCluster>(tcH,index);
 /*------------------------- FITTING STAGE 1 :---------------------------//
 Order hits in terms of height 
//------------------------- FITTING STAGE 1 :---------------------------*/
      OrderHitsY(_stResult); 
      
/*------------------------- FITTING STAGE 2 a :---------------------------//
Ensure good clusters used 
//------------------------- FITTING STAGE 2 a:---------------------------*/
     
     // if ( nGoodTClusterHits < _minNHitsTimeCluster)         continue;
/*------------------------- FITTING STAGE 2 b :---------------------------//
Ensure enough straw hits used
//------------------------- FITTING STAGE 2 a:---------------------------*/
      if (_debug != 0){
	 std::cout<<"#filtered SHits"<<_stResult._nFiltStrawHits<<"  Minimum allowed: "<<_minnsh<<std::endl;
      }
      
      
     
      if (_stResult._nFiltStrawHits < _minnsh)                  continue;
      
      _stResult._tseed._status.merge(TrkFitFlag::hitsOK);
      if(event.id().event() ==14 || event.id().event() ==78 || event.id().event() ==121){
	std::cout<< "hits ok"<<std::endl;
      }
      if (_diag) _stResult._diag.StraightTrackFitCounter = 0;
 /*------------------------- FITTING STAGE 3 :---------------------------//
 initial ST Pat Rec.. Call the "Fit" function and apply to the StraightTrackData 
//------------------------- FITTING STAGE 3 :---------------------------*/
      
      _tfit.BeginFit(_stResult);
      
      //Fill in diagnostics:
       if (_diag) {
	_stResult._diag.nShFitXY = _stResult._nXYSh;
	_stResult._diag.nChFitXY = _stResult._tseed._track.get_chisq_dof();
      } //end "diag"
    
      if (_stResult._tseed._status.hasAnyProperty(TrkFitFlag::StraightTrackOK)) { 
	
	      //tentatively store as a temporary result
	      StraightTrackFinderData tmpResult(_stResult);
	      

 /*------------------------- FITTING STAGE 4 :---------------------------//
Fill Seed Info into StraightTrackFinder_seed
//------------------------- FITTING STAGE 4 :---------------------------*/
	  
	      _stResult._tseed._status.merge(TrkFitFlag::StraightTrackOK);
              if (tmpResult._tseed.status().hasAnyProperty(_saveflag)){
		      //fille the hits in the seed collection
		      fillGoodHits(tmpResult);
		      StraightTrackSeedCollection* col = seed_col.get();
		      col->push_back(tmpResult._tseed);
                      
              }//end "saveflag"
      }//end "track ok" criteria
/*------------------------- FITTING STAGE 6 ----------------------------//
Fill Diagnostic Info.
//---------------------------------------------------------------------*/
      if (_diag > 0) {//TODO:Add these into stresult
            _hmanager->fillHistograms(&_data);
            
       } //end "diag"
       ++STCounter;
/*------------------------- FITTING STAGE 7 ----------------------------//
put into event data 
//---------------------------------------------------------------------*/
    }//end loop
  event.put(std::move(seed_col));    //here - more to add to event?
  std::cout<<"TC fails" << n_TC_fails<<std::endl;
  }//end produce


/*---------------------------Get GOOD HITS-----------------------//
//  Asks if an outlier  - removes if so  - is this necessary for cosmics  //
//---------------------------------------------------------------*/ 
  void StraightTrackFinder::fillGoodHits(StraightTrackFinderData& trackData){
    if (_debug != 0) {
	std::cout<<"Filling good hits..."<<std::endl;
    }
    ComboHit*     hit(0);
    for (unsigned f=0; f<trackData._chHitsToProcess.size(); ++f){
      hit = &trackData._chHitsToProcess[f];
      if (hit->_flag.hasAnyProperty(_outlier))     continue;
      
      ComboHit                thit(*hit);					
      trackData._tseed._thits.push_back(thit);
    }//end for loop

  }//end fill good function



  /*--------------Order Hits------------------------//
  // Need to organise tracks so that the hits in one XY plane are seeded from above. This will invlove ordering in terms of wire "Y" position. 
QUESTION: Are we using wire position or absolute position...wire+/-dres...?
 //------------------------------------------------*/
  void StraightTrackFinder::OrderHitsY(StraightTrackFinderData& TrackData){
    if (_debug != 0){
	 std::cout<<"Ordering Hits..."<<std::endl;
    }
    const vector<StrawHitIndex>& shIndices = TrackData._timeCluster->hits();
    mu2e::StraightTrackFinderData::ChannelID cx, co;

    int	    h;
    int     size  = shIndices.size();
    int     nFiltComboHits(0), nFiltStrawHits(0);
    
    ComboHitCollection ordChCol;
    ordChCol.reserve(size);
    //loop over hits
    
    for (int i=0; i<size; ++i) {
      h = shIndices[i];
      const ComboHit& ch  = (*_stResult._chcol)[h];
      ordChCol.push_back(ComboHit(ch));
      
    }//end loop
    //sort hits (here in y..need to reconsider for cosmics, should be "height")
    if (_debug != 0) std::cout<<"Number of ComboHits: "<<ordChCol.size()<<std::endl;
    std::sort(ordChCol.begin(), ordChCol.end(),ycomp());
    
    
    //Loop through order list now and update result:
    for (unsigned i=0; i<ordChCol.size(); ++i) {
      
      ComboHit& ch = ordChCol[i];
      //Add ordered hits back to the StraightTrackFinder Result:
      ComboHit hit(ch);
      _stResult._chHitsToProcess.push_back(hit);
      //set IDs:
      cx.Station                 = ch.strawId().station();
      cx.Plane                   = ch.strawId().plane() % 2;
      cx.Face                    = ch.strawId().face();
      cx.Panel                   = ch.strawId().panel();
      // get Height-ordered location
      TrackData.orderID(&cx, &co);
      //int os       = co.Station; 
      int of       = co.Face;
      //int op       = co.Panel;
      //sort wire position - reorder st hit results to reflect this Y based ordering...
      _stResult._chHitsWPos.push_back(XYWVec(hit.pos(),  of, hit.nStrawHits()));
      
      ++nFiltComboHits;
      nFiltStrawHits += ch.nStrawHits();
      
    }//end for loop
    TrackData._nFiltComboHits = nFiltComboHits;  //ComboHit counter
    TrackData._nFiltStrawHits = nFiltStrawHits;  //StrawHit counter
  }//end order function

  
  



/*--------------- GOOD HITS TIME CLUSTER ---------------------//
//             Gets the time clusters deemed "good" for fit   //
//------------------------------------------------------------*/
int  StraightTrackFinder::goodHitsTimeCluster(const TimeCluster TCluster, ComboHitCollection chcol){
    int   nhits         = TCluster.nhits();
    int   ngoodhits(0);
    
    double     minT(500.), maxT(2000.);// TODO: Check Where do these come from?
    for (int i=0; i<nhits; ++i){
      int          index   = TCluster.hits().at(i);
      
      ComboHit     sh      = chcol .at(index);
      
      
      if ( (sh.time() < minT) || (sh.time() > maxT) )  continue;

      // ++ngoodhits;
      ngoodhits += sh.nStrawHits();
    }

    return ngoodhits;
  } 
///////////////////////////////////////////////////
}//end mu2e namespace
using mu2e::StraightTrackFinder;
DEFINE_ART_MODULE(StraightTrackFinder);

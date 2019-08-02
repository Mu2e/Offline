// Cosmic Track finder- module calls seed fitting routine to begin cosmic track analysis
// $Id: CosmicTrackFinder_module.cc
// $Author: S.middleton
// $Date: Feb 2019 $
#include "TrkReco/inc/CosmicTrackFit.hh"
#include "TrkPatRec/inc/CosmicTrackFinder_types.hh"
#include "TrkReco/inc/CosmicTrackFinderData.hh"

//Mu2e General:
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/TimeClusterCollection.hh"
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
//MC:
// art
#include "canvas/Persistency/Common/Ptr.h"
// MC data
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
//MU2E:
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "TrkReco/inc/CosmicTrackFit.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TrkReco/inc/TrkTimeCalculator.hh"
//utils:
#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
//Track Fitting Headers:
#include "TrkPatRec/inc/CosmicTrackFinder_types.hh"
#include "TrkReco/inc/CosmicTrackFinderData.hh"
#include "TrkReco/inc/TrkFaceData.hh"
//For Drift:
#include "TrkReco/inc/PanelAmbigResolver.hh"
#include "TrkReco/inc/PanelStateIterator.hh"
// Mu2e BaBar
#include "BTrkData/inc/TrkStrawHit.hh"

//CLHEP:
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
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
using CLHEP::Hep3Vector;
using CLHEP::HepVector;

namespace{
    //create a compare struct to allow height ordering of the hits in an event. 
    struct ycomp : public std::binary_function<mu2e::ComboHit,mu2e::ComboHit,bool> {
    bool operator()(mu2e::ComboHit const& p1, mu2e::ComboHit const& p2) { return p1._pos.y() > p2._pos.y(); }
  };
  
    struct ycomp_digi : public std::binary_function<mu2e::StrawDigiMC,mu2e::StrawDigiMC,bool> {
    bool operator()(mu2e::StrawDigiMC const& p1, mu2e::StrawDigiMC const& p2) {
    art::Ptr<mu2e::StepPointMC> const& spmcp1 = p1.stepPointMC(mu2e::StrawEnd::cal); 
    art::Ptr<mu2e::StepPointMC> const& spmcp2 = p2.stepPointMC(mu2e::StrawEnd::cal);
    std::cout<<spmcp1->position().y() <<" "<< spmcp2->position().y()<<std::endl;
    return spmcp1->position().y() > spmcp2->position().y(); }
  };
  //For MC StrawDigis:
  struct ycomp_MC : public std::binary_function<mu2e::StrawDigiMC,mu2e::StrawDigiMC,bool> {
    bool operator()(art::Ptr<mu2e::StepPointMC> const& spmcp1, art::Ptr<mu2e::StepPointMC> const& spmcp2) {
    return spmcp1->position().y() > spmcp2->position().y(); }
  };

    struct zcomp : public std::binary_function<mu2e::ComboHit,mu2e::ComboHit,bool> {
    bool operator()(mu2e::ComboHit const& p1, mu2e::ComboHit const& p2) { return p1._pos.z() < p2._pos.z(); }
  }; 
  
  

}//end namespace

namespace mu2e{

  class CosmicTrackFinder : public art::EDProducer {
  public:
    explicit CosmicTrackFinder(fhicl::ParameterSet const&);
    virtual ~CosmicTrackFinder();
    virtual void beginJob();
    virtual void beginRun(art::Run& run);
    virtual void produce(art::Event& event );
  private:
  //config parameters:
    int 				_diag,_mcdiag, _debug;
    int                                 _printfreq;
    int 				_minnsh; // minimum # of strawHits in CH
    int 				_minnch; // minimum # of ComboHits for viable fit
    TrkFitFlag				_saveflag;//write tracks that satisfy these flags
    //unsigned				_maxniter;  // maximum # of iterations over outlier filtering + fitting 
    int 				_minNHitsTimeCluster; //min number of hits in a viable time cluster

    //Collection lists:
    art::ProductToken<ComboHitCollection> const _chToken;
    art::ProductToken<TimeClusterCollection> const _tcToken;
    art::ProductToken<StrawDigiMCCollection> const _mcToken;
      
    //Define Track Fits
    CosmicTrackFit     _tfit;
    //TrkTimeCalculator _ttcalc;
    //float             _t0shift; 
    StrawHitFlag      _outlier;
   
    std::unique_ptr<ModuleHistToolBase>   _hmanager;
    CosmicTrackFinderTypes::Data_t        _data;
    CosmicTrackFinderData                 _stResult;
    ProditionsHandle<StrawResponse> _strawResponse_h; 
    void     OrderHitsY(CosmicTrackFinderData& TrackData); //Order in height
    void     fillGoodHits(CosmicTrackFinderData& TrackData);//apply "good" cut
    void     fillPluginDiag(CosmicTrackFinderData& TrackData);
    int      goodHitsTimeCluster(const TimeCluster TCluster, ComboHitCollection chcol);//select "good" time clusters
};//end private

//Constructor:
 CosmicTrackFinder::CosmicTrackFinder(fhicl::ParameterSet const& pset) :
   
    _diag        (pset.get<int>("diagLevel",1)),
    _mcdiag      (pset.get<int>("mcdiagLevel",1)),
    _debug       (pset.get<int>("debugLevel",0)),
    _printfreq   (pset.get<int>   ("printFrequency", 101)),
    _minnsh      (pset.get<int>("minNStrawHits",2)),
    _minnch      (pset.get<int>("minNComboHits",4)),
    _saveflag    (pset.get<vector<string> >("SaveTrackFlag",vector<string>{"StraightTrackOK"})),
    _minNHitsTimeCluster(pset.get<int>("minNHitsTimeCluster", 1 )), 
    _chToken{consumes<ComboHitCollection>(pset.get<art::InputTag>("ComboHitCollection"))},
    _tcToken{consumes<TimeClusterCollection>(pset.get<art::InputTag>("TimeClusterCollection"))},
    _mcToken{consumes<StrawDigiMCCollection>(pset.get<art::InputTag>("StrawDigiMCCollection"))},
    _tfit        (pset.get<fhicl::ParameterSet>("CosmicTrackFit",fhicl::ParameterSet())), 
    //_ttcalc      (pset.get<fhicl::ParameterSet>("T0Calculator",fhicl::ParameterSet())),
    //_t0shift     (pset.get<float>("T0Shift",4.0)),
    _outlier     (StrawHitFlag::outlier)
  {
    produces<CosmicTrackSeedCollection>();
    if (_diag != 0) _hmanager = art::make_tool<ModuleHistToolBase>(pset.get<fhicl::ParameterSet>("diagPlugin"));
    else    _hmanager = std::make_unique<ModuleHistToolBase>();
 }

//Destructor:
 CosmicTrackFinder::~CosmicTrackFinder(){}

/* ------------------------Begin JOb--------------------------//
//         Sets Up Historgram Book For Diag Plots            //
//----------------------------------------------------------*/

  void CosmicTrackFinder::beginJob() {
    if (_debug != 0) std::cout<<"Beginning STFinder Job..."<<std::endl;
    art::ServiceHandle<art::TFileService> tfs;

    if (_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;  
      _hmanager->bookHistograms(tfs);
    } 
  }
/* ------------------------Begin Run--------------------------//
//                   sets up the tracker                     //
//----------------------------------------------------------*/
  void CosmicTrackFinder::beginRun(art::Run& run) {
    
// Sets up tracker information....adds this to the Cosmictrackfinder...
    if (_debug > 0)
    {
      std::cout << "BEGINNING Cosmic TRACK FINDING... " << std::endl;
    }
   mu2e::GeomHandle<mu2e::Tracker> th;
   const Tracker* tracker = th.get();
   auto const& srep = _strawResponse_h.get(run.id());
   _tfit.setTracker  (tracker);
   _tfit.setStrawResponse (srep);
   }

/* ------------------------PRODUCE--------------------------//
//                  Should produce ST seed Collection      //
//----------------------------------------------------------*/
  void CosmicTrackFinder::produce(art::Event& event ) {
     
     if (_debug != 0) std::cout<<"Producing Cosmic Track in  Finder..."<<std::endl;
     //Create output collection - this is collection of straght track seeds:
     unique_ptr<CosmicTrackSeedCollection> seed_col(new CosmicTrackSeedCollection());
     unsigned    STCounter(0);
     
     _stResult.clearMCVariables();
     int _iev=event.id().event();
      if (_debug > 0){
          std::cout<<"ST Finder Event #"<<_iev<<std::endl;
      } //end "if"
     
     //fill in collection lists
     auto const& chH = event.getValidHandle(_chToken);
     const ComboHitCollection& chcol(*chH);
     auto const& tcH = event.getValidHandle(_tcToken);
     const TimeClusterCollection& tccol(*tcH);

     if(_mcdiag>0){ 	 
           auto const& mcdH = event.getValidHandle(_mcToken);
           const StrawDigiMCCollection& mccol(*mcdH);
           _stResult._mccol =  &mccol;
           for(size_t imc=0; imc<mccol.size(); imc++){
                StrawDigiMC digi = mccol[imc];
           	_stResult._mcDigisToProcess.push_back(digi);
           }  
        }
    _data.event       = &event;
    _data.nTimePeaks  = tccol.size();
    _stResult._chcol  = &chcol; 
    _stResult._tccol  = &tccol;
    // create initial tracks from time clusters:
    for (size_t index=0;index< tccol.size();++index) {
      int   nGoodTClusterHits(0);
      //get the t cluster
      const auto& tclust = tccol[index];
      // ensure only good clusters will be used for fit:
      nGoodTClusterHits     = goodHitsTimeCluster(tclust,chcol );
      
/*------------------------- FITTING STAGE :---------------------------//
Ensure good clusters used 
//------------------------- FITTING STAGE---------------------------*/
      if ( nGoodTClusterHits < _minNHitsTimeCluster)         continue;
      // create a place for the seed:
       CosmicTrackSeed tseed ;
      //clear the variables in _stResult
      _stResult.clearTempVariables();
      //set variables used for searching the track candidate. pass these to the Cosmictrackfinderdata result (stResult)
      _stResult._tseed              = tseed;
      _stResult._timeCluster        = &tclust;
      _stResult._tseed._thits.setParent(chcol.parent());
      _stResult._tseed._t0          = tclust._t0;
      _stResult._tseed._timeCluster = art::Ptr<TimeCluster>(tcH,index);
 /*------------------------- FITTING STAGE :---------------------------//
 Order hits in terms of height 
//------------------------- FITTING STAGE  ---------------------------*/
      OrderHitsY(_stResult); 
/*------------------------- FITTING STAGE  ---------------------------//
Ensure enough straw hits used
//------------------------- FITTING STAGE ---------------------------*/
      if (_debug != 0){
	 std::cout<<"#filtered SHits"<<_stResult._nFiltStrawHits<<"  Minimum allowed: "<<_minnsh<<std::endl;
      }
      if (_stResult._nFiltComboHits < _minnch ) 	continue;
      if (_stResult._nFiltStrawHits < _minnsh)          continue;
     
      _stResult._nStrawHits = _stResult._nFiltStrawHits;
      _stResult._nComboHits = _stResult._nFiltComboHits;
      _stResult._tseed._status.merge(TrkFitFlag::hitsOK);
      
      if (_diag) _stResult._diag.CosmicTrackFitCounter = 0;
 /*------------------------- FITTING STAGE  :---------------------------//
 initial ST Pat Rec.. Call the "Fit" function and apply to the CosmicTrackData 
//------------------------- FITTING STAGE  :---------------------------*/
      if(_mcdiag > 0 && _stResult._chHitsToProcess.size() > 0){
           XYZVec first, last;
             //std::cout<<"In mc diag loop "<<std::endl;
             //StrawDigiMCCollection OrderMCDigi;
             
	     for(size_t ich= 0; ich<_stResult._chHitsToProcess.size(); ich++) {  
	        
           	std::vector<StrawDigiIndex> shids;           	
           	ComboHit const& chit = _stResult._chHitsToProcess.at(ich);
           	
           	if (chit.nCombo() >  chit.nStrawHits() ) continue;
                _stResult._chHitsToProcess.fillStrawDigiIndices(event,ich,shids);               
              
                StrawDigiMC const& mcd1 = _stResult._mccol->at(shids[0]);              
                //art::Ptr<StepPointMC> const& spmcp = mcd1.stepPointMC(StrawEnd::cal);
                //_stResult._mcDigisToProcess.push_back(mcd1);               
                
		if(ich ==0){
			first = _tfit.MCInitHit(mcd1) ; 	
		}
		if(ich == _stResult._chHitsToProcess.size()-1){
			last = _tfit.MCFinalHit(mcd1); 	
		} 
	     }
	     _tfit.MCDirection(first, last, _stResult); 
	     
	     
	     }
     ostringstream title;
     title << "Run: " << event.id().run()
     << "  Subrun: " << event.id().subRun()
     << "  Event: " << event.id().event()<<".root";
     
      _tfit.BeginFit(title.str().c_str(), _stResult, _data);

      //Fill in diagnostics:
       if (_diag) {
	_stResult._diag.nShFit = _stResult._nXYSh;
	_stResult._diag.nChFit = _stResult._nXYCh;
      } //end "diag"
     
      if (_stResult._tseed._status.hasAnyProperty(TrkFitFlag::StraightTrackOK) && _stResult._tseed._status.hasAnyProperty(TrkFitFlag::StraightTrackConverged)  ) { 
	       std::vector<CosmicTrackSeed>          track_seed_vec;
	      //tentatively store as a temporary result
	      CosmicTrackFinderData tmpResult(_stResult);
	     
 /*------------------------- FITTING STAGE  :---------------------------//
Fill Seed Info into CosmicTrackFinder_seed
//------------------------- FITTING STAGE  :---------------------------*/
	      _data.nseeds += 1;
	      _stResult._tseed._status.merge(TrkFitFlag::StraightTrackOK);
              if (tmpResult._tseed.status().hasAnyProperty(_saveflag)){
              
		      //fille the hits in the seed collection
		      fillGoodHits(tmpResult);
		      track_seed_vec.push_back(tmpResult._tseed);
		      if (_diag > 0) {
	      		fillPluginDiag(tmpResult);
	    	      }
		      CosmicTrackSeedCollection* col = seed_col.get();
		      //TODO - refit here to get outliers removed and optimization
		      if (track_seed_vec.size() == 0)     continue;
		      col->push_back(tmpResult._tseed);                 
              }//end "saveflag"
      }//end "track ok" criteria
/*------------------------- FITTING STAGE  ----------------------------//
Fill Diagnostic Info.
//---------------------------------------------------------------------*/
      if (_diag > 0) _hmanager->fillHistograms(&_data); 
       ++STCounter;
/*------------------------- FITTING STAGE  ----------------------------//
put into event data 
//---------------------------------------------------------------------*/
    }//end loop
  event.put(std::move(seed_col));    //here - more to add to event?
  
  }//end produce
 
      

/*---------------------------Get GOOD HITS-----------------------//
//  Asks if an outlier  - removes if so  - is this necessary for cosmics  //
//---------------------------------------------------------------*/ 
  void CosmicTrackFinder::fillGoodHits(CosmicTrackFinderData& trackData){
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
  void CosmicTrackFinder::OrderHitsY(CosmicTrackFinderData& TrackData){
    if (_debug != 0){
	 std::cout<<"Ordering Hits..."<<std::endl;
    }
    const vector<StrawHitIndex>& shIndices = TrackData._timeCluster->hits();
    mu2e::CosmicTrackFinderData::ChannelID cx, co;

    int	    h;
    int     size  = shIndices.size();
    int     nFiltComboHits(0), nFiltStrawHits(0);
    
    ComboHitCollection ordChCol;
    ordChCol.reserve(size);
 
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
    
      //Add ordered hits back to the CosmicTrackFinder Result:
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
  
//Currently not using Diag module but may "upgrade" so leave this here for now...
  void CosmicTrackFinder::fillPluginDiag(CosmicTrackFinderData& trackData) {
    
    int loc = _data.nseeds;
  //  printf(" N(seeds) > %i, IGNORE SEED\n",_data.maxSeeds());
    
      _data.nChPPanel[loc] = trackData._diag.nChPPanel;
      _data.nChHits[loc] = trackData._diag.nChHits;
      int nhits          = trackData._tseed._thits.size();
      _data.ntclhits[loc] = trackData._timeCluster->hits().size();
      _data.nhits[loc] = nhits;
      _data.nShFit[loc] = trackData._diag.nShFit;
      _data.nChFit[loc] = trackData._diag.nChFit;
      _data.niters[loc] = trackData._diag.niters;
      for (int i=0; i<trackData._diag.nChFit; ++i) {
        
	_data.hit_residualX[loc][i] = trackData._diag.hit_residualX[i];
	_data.hit_residualY[loc][i] = trackData._diag.hit_residualY[i];
    }
  
}


/*--------------- GOOD HITS TIME CLUSTER ---------------------//
//             Gets the time clusters deemed "good" for fit   //
//------------------------------------------------------------*/
int  CosmicTrackFinder::goodHitsTimeCluster(const TimeCluster TCluster, ComboHitCollection chcol){
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
using mu2e::CosmicTrackFinder;
DEFINE_ART_MODULE(CosmicTrackFinder);

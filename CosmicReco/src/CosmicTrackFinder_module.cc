// Author: S Middleton
// Purpose: Cosmic Track finder- module calls seed fitting routine to begin cosmic track analysis. The module can call the seed fit and drift fit. Producing a "CosmicTrackSeed" list.

#include "CosmicReco/inc/CosmicTrackFit.hh"
#include "CosmicReco/inc/CosmicTrackFinderData.hh"

#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
//Mu2e General:
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "TrackerGeom/inc/Tracker.hh"
// ART:
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
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
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "TrkReco/inc/TrkTimeCalculator.hh"
//utils:
#include "Mu2eUtilities/inc/ParametricFit.hh"
//For Drift:
#include "TrkReco/inc/PanelAmbigResolver.hh"
#include "TrkReco/inc/PanelStateIterator.hh"
#include "TrkReco/inc/TrkFaceData.hh"
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
    struct Config{
	      using Name=fhicl::Name;
	      using Comment=fhicl::Comment;
	      fhicl::Atom<int> mcdiag{Name("mcdiag"), Comment("set on for MC info"),2};
	      fhicl::Atom<int> debug{Name("debugLevel"), Comment("set to 1 for debug prints"),1};
              fhicl::Atom<int> printfreq{Name("printFrequency"), Comment("print frquency"), 101};
    	      fhicl::Atom<int> minnsh {Name("minNStrawHits"), Comment("minimum number of straw hits "),2};
    	      fhicl::Atom<int> minnch {Name("minNComboHits"), Comment("number of combohits allowed"),8};
    	      fhicl::Atom<TrkFitFlag> saveflag {Name("SaveTrackFlag"),Comment("if set to OK then save the track"), TrkFitFlag::helixOK};
              fhicl::Atom<int> minNHitsTimeCluster{Name("minNHitsTimeCluster"),Comment("minium allowed time cluster"), 1 };
    	      fhicl::Atom<art::InputTag> chToken{Name("ComboHitCollection"),Comment("tag for combo hit collection")};
	      fhicl::Atom<art::InputTag> tcToken{Name("TimeClusterCollection"),Comment("tag for time cluster collection")};
	      fhicl::Atom<art::InputTag> mcToken{Name("StrawDigiMCCollection"),Comment("StrawDigi collection tag")};
	      fhicl::Table<CosmicTrackFit::Config> tfit{Name("CosmicTrackFit"), Comment("fit")};
    };
    typedef art::EDProducer::Table<Config> Parameters;
    explicit CosmicTrackFinder(const Parameters& conf);
    virtual ~CosmicTrackFinder();
    virtual void beginJob();
    virtual void beginRun(art::Run& run);
    virtual void produce(art::Event& event ) override;
    
  private:
    
    Config _conf;

  //config parameters:
    int 				_mcdiag, _debug;
    int                                 _printfreq;
    int 				_minnsh; // minimum # of strawHits in CH
    int 				_minnch; // minimum # of ComboHits for viable fit
    TrkFitFlag				_saveflag;//write tracks that satisfy these flags
    
    int 				_minNHitsTimeCluster; //min number of hits in a viable time cluster
    float				_max_seed_chi2; ///maximum chi2 allowed for seed
    
    art::InputTag  _chToken;
    art::InputTag  _tcToken;
    art::InputTag  _mcToken;
    
    CosmicTrackFit     _tfit;
     
    StrawHitFlag      _outlier;
   
    CosmicTrackFinderData                 _stResult;
    ProditionsHandle<StrawResponse> _strawResponse_h; 
    void     OrderHitsY(CosmicTrackFinderData& TrackData); //Order in height
    void     OrderHitsYMC(CosmicTrackFinderData& TrackData, art::Event& event); //Same but for MCDigis
    void     fillGoodHits(CosmicTrackFinderData& TrackData);//apply "good" cut
    int      goodHitsTimeCluster(const TimeCluster TCluster, ComboHitCollection chcol);
   
};


 CosmicTrackFinder::CosmicTrackFinder(const Parameters& conf) :
   art::EDProducer(conf),
   	_mcdiag (conf().mcdiag()),
	_debug  (conf().debug()),
    	_printfreq  (conf().printfreq()),
    	_minnsh   (conf().minnsh()),
    	_minnch  (conf().minnch()),
    	_saveflag  (conf().saveflag()),
    	_minNHitsTimeCluster(conf().minNHitsTimeCluster()),
    	_chToken (conf().chToken()),
	_tcToken (conf().tcToken()),
	_mcToken (conf().mcToken()),
	_tfit (conf().tfit())
{
	    consumes<ComboHitCollection>(_chToken);
	    consumes<TimeClusterCollection>(_tcToken);
	    consumes<StrawDigiMCCollection>(_mcToken);
	    produces<CosmicTrackSeedCollection>();
	    
 }

 CosmicTrackFinder::~CosmicTrackFinder(){}

/* ------------------------Begin JOb--------------------------//
//         Sets Up Historgram Book For Diag Plots            //
//----------------------------------------------------------*/

  void CosmicTrackFinder::beginJob() {
   
    art::ServiceHandle<art::TFileService> tfs;
    
  }
/* ------------------------Begin Run--------------------------//
//                   sets up the tracker                     //
//----------------------------------------------------------*/
  void CosmicTrackFinder::beginRun(art::Run& run) {
    
   mu2e::GeomHandle<mu2e::Tracker> th;
   const Tracker* tracker = th.get();
  
   _stResult.run = &run;
   auto const& srep = _strawResponse_h.get(run.id());
   _tfit.setTracker  (tracker);
   _tfit.setStrawResponse (srep);
   }

  void CosmicTrackFinder::produce(art::Event& event ) {
     
     if (_debug != 0) std::cout<<"Producing Cosmic Track in  Finder..."<<std::endl;
     unique_ptr<CosmicTrackSeedCollection> seed_col(new CosmicTrackSeedCollection());
     
     _stResult.clearMCVariables();
     int _iev=event.id().event();
      if (_debug > 0){
          std::cout<<"ST Finder Event #"<<_iev<<std::endl;
      } 
     
     auto const& chH = event.getValidHandle<ComboHitCollection>(_chToken);
     const ComboHitCollection& chcol(*chH);
     auto  const& tcH = event.getValidHandle<TimeClusterCollection>(_tcToken);
     const TimeClusterCollection& tccol(*tcH);
     
     if(_mcdiag>0){ 	 
           auto const& mcdH = event.getValidHandle<StrawDigiMCCollection>(_mcToken);
           const StrawDigiMCCollection& mccol(*mcdH);
           _stResult._mccol =  &mccol;
       }
    
    
    _stResult.event   = &event;
    _stResult._chcol  = &chcol; 
    _stResult._tccol  = &tccol;
    
    
    for (size_t index=0;index< tccol.size();++index) {
      int   nGoodTClusterHits(0);
      const auto& tclust = tccol[index];
      nGoodTClusterHits     = goodHitsTimeCluster(tclust,chcol );
    
      if ( nGoodTClusterHits < _minNHitsTimeCluster)         continue;
      if (_debug > 0){
          std::cout<<"time clusters "<<_iev<<std::endl;
      }
       CosmicTrackSeed tseed ;
      _stResult.clearTempVariables();
      _stResult._tseed              = tseed;
      _stResult._timeCluster        = &tclust;
      _stResult._chHitsToProcess.setParent(chcol.parent());
      _stResult._tseed._panel_hits.setParent(chcol.parent());
      _stResult._tseed._t0          = tclust._t0;
      _stResult._tseed._timeCluster = art::Ptr<TimeCluster>(tcH,index);
      
      OrderHitsY(_stResult); 

      if (_debug != 0){
	 std::cout<<"#filtered SHits"<<_stResult._nFiltStrawHits<<" #filter CHits "<<_stResult._nFiltComboHits<<std::endl;
      }
      if (_stResult._nFiltComboHits < _minnch ) 	continue;
      if (_stResult._nFiltStrawHits < _minnsh)          continue;
      _stResult._tseed._status.merge(TrkFitFlag::Straight);
      _stResult._nStrawHits = _stResult._nFiltStrawHits;
      _stResult._nComboHits = _stResult._nFiltComboHits;
      _stResult._tseed._status.merge(TrkFitFlag::hitsOK);
      _stResult._tseed._panel_hits = _stResult._chHitsToProcess;
      
      if(_mcdiag > 0 && _stResult._chHitsToProcess.size() > 0){
		OrderHitsYMC(_stResult, event);
		for(auto const & mcd : _stResult._mcDigisToProcess){
			art::Ptr<StepPointMC> const& spmcp = mcd.stepPointMC(StrawEnd::cal);
            		XYZVec posN(spmcp->position().x(), spmcp->position().y(), spmcp->position().z());
                	
		}
			
             
	     }

     ostringstream title;
     title << "Run: " << event.id().run()
     << "  Subrun: " << event.id().subRun()
     << "  Event: " << event.id().event()<<".root";
     _tfit.BeginFit(title.str().c_str(), _stResult);

      if (_stResult._tseed._status.hasAnyProperty(TrkFitFlag::helixOK) && _stResult._tseed._status.hasAnyProperty(TrkFitFlag::helixConverged) && _stResult._tseed._track.converged == true ) { 
	       std::vector<CosmicTrackSeed>          track_seed_vec;
	       
	      fillGoodHits(_stResult);
	      
	      CosmicTrackFinderData tmpResult(_stResult);
	      _stResult._tseed._status.merge(TrkFitFlag::helixOK);
              if (tmpResult._tseed.status().hasAnyProperty(_saveflag)){
              	
		      std::vector<uint16_t> chindices;
		      if(tmpResult._tseed._track.converged == false) continue;
		      for(size_t ich= 0; ich<_stResult._chHitsToProcess.size(); ich++) { 
		   	 chindices.push_back(ich);
		   	
	              }
	              //get list of indices to straw level combohits
	              std::vector<ComboHitCollection::const_iterator> chids;  
		      tmpResult._chHitsToProcess.fillComboHits(event, chindices, chids); 
		      std::vector<ComboHitCollection::const_iterator> StrawLevelCHitIndices = chids;
		      for (auto const& it : chids){
		      	 //it = "_normal_iterator<const mu2e::ComboHit*, std::vector<mu2e::ComboHit> >"
		      	  const mu2e::ComboHit chit = it[0];
		      	  tmpResult._tseed._straw_chits.push_back(chit);
		      	  Straw const& straw = _tfit._tracker->getStraw(chit.strawId()); 
           		  tmpResult._tseed._straws.push_back(straw);
           		  
	      	      }
	
	      	      for(size_t ich= 0; ich<tmpResult._tseed._straw_chits.size(); ich++) {  
           	 
           		std::vector<StrawHitIndex> shitids;          	          		
           	        tmpResult._tseed._straw_chits.fillStrawHitIndices(event, ich,shitids);  
                        tmpResult._tseed._strawHitIdxs.push_back(ich);
               
			for(auto const& ids : shitids){ 
				size_t    istraw   = (ids);
			     	TrkStrawHitSeed tshs;
			     	tshs._index  = istraw;
			     	tshs._t0 = tclust._t0;
			     	tmpResult._tseed._trkstrawhits.push_back(tshs); 
	     		}  
	     		}
	              
                      if( _tfit.goodTrack(tmpResult._tseed._track) == false) continue;
		      _tfit.DriftFit(tmpResult);
		      
		      //Add tmp to seed list:
		      track_seed_vec.push_back(tmpResult._tseed);
		     
		      CosmicTrackSeedCollection* col = seed_col.get();
		      
		      if (track_seed_vec.size() == 0)     continue;
		      col->push_back(tmpResult._tseed);  
			          
              }
        }
    }
 
  event.put(std::move(seed_col));    
  }

  void CosmicTrackFinder::fillGoodHits(CosmicTrackFinderData& trackData){
    if (_debug != 0) {
	std::cout<<"Filling good hits..."<<std::endl;
    }
    ComboHit*     hit(0);
    for (unsigned f=0; f<trackData._chHitsToProcess.size(); ++f){
      hit = &trackData._chHitsToProcess[f];
      //if (hit->_flag.hasAnyProperty(_outlier))     continue;
      
      ComboHit                thit(*hit);					
      trackData._tseed._panel_hits.push_back(thit);
    }
  }
void CosmicTrackFinder::OrderHitsYMC(CosmicTrackFinderData& TrackData, art::Event& event){
    int     size  = _stResult._chHitsToProcess.size();
    StrawDigiMCCollection ordDigiCol;
    ordDigiCol.reserve(size);
 
    for (int i=0; i<size; ++i) {
     std::vector<StrawDigiIndex> shids;  
     _stResult._chHitsToProcess.fillStrawDigiIndices(event,i,shids);    
     StrawDigiMC const& mcd1 = _stResult._mccol->at(shids[0]);  
      ordDigiCol.push_back(mcd1); 
    }
   
    if (_debug != 0) std::cout<<"Number of Digis: "<<ordDigiCol.size()<<std::endl;
    std::sort(ordDigiCol.begin(), ordDigiCol.end(),ycomp_digi());

    for (unsigned i=0; i<ordDigiCol.size(); ++i) { 
      StrawDigiMC const& mcd1 = ordDigiCol[i]; 
     _stResult._mcDigisToProcess.push_back(mcd1);

    }

}

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
    }
   
    if (_debug != 0) std::cout<<"Number of ComboHits: "<<ordChCol.size()<<std::endl;
    std::sort(ordChCol.begin(), ordChCol.end(),ycomp());
    for (unsigned i=0; i<ordChCol.size(); ++i) { 
      ComboHit& ch = ordChCol[i];
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
    }
    TrackData._nFiltComboHits = nFiltComboHits;  //ComboHit counter
    TrackData._nFiltStrawHits = nFiltStrawHits;  //StrawHit counter
  }
  
int  CosmicTrackFinder::goodHitsTimeCluster(const TimeCluster TCluster, ComboHitCollection chcol){
    int   nhits         = TCluster.nhits();
    int   ngoodhits(0);
    double     minT(500.), maxT(2000.);
    for (int i=0; i<nhits; ++i){
      int          index   = TCluster.hits().at(i);
      ComboHit     sh      = chcol.at(index); 
      if ( (sh.time() < minT) || (sh.time() > maxT) )  continue;
  
      ngoodhits += sh.nStrawHits();
    }

    return ngoodhits;
  } 

///////////////////////////////////////////////////
}//end mu2e namespace
using mu2e::CosmicTrackFinder;
DEFINE_ART_MODULE(CosmicTrackFinder);

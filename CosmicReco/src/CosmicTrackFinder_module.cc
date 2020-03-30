// Author : S Middleton
// Data : March 2019
// Purpose: Cosmic Track finder- module calls seed fitting routine to begin cosmic track analysis. The module can call the seed fit and drift fit. Producing a "CosmicTrackSeed" list.

#include "CosmicReco/inc/CosmicTrackFit.hh"
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
#include "canvas/Persistency/Common/Ptr.h"

//MU2E:
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "TrkReco/inc/TrkTimeCalculator.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "RecoDataProducts/inc/LineSeed.hh"

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

    struct ycomp_iter : public std::binary_function<std::vector<ComboHit>::const_iterator, std::vector<ComboHit>::const_iterator, bool> {
      bool operator()(std::vector<ComboHit>::const_iterator p1, std::vector<ComboHit>::const_iterator p2) { return p1->_pos.y() > p2->_pos.y(); }
    };

}

namespace mu2e{

    class CosmicTrackFinder : public art::EDProducer {
    public:
	struct Config{
	      using Name=fhicl::Name;
	      using Comment=fhicl::Comment;
	      fhicl::Atom<int> debug{Name("debugLevel"), Comment("set to 1 for debug prints"),1};
	      fhicl::Atom<int> printfreq{Name("printFrequency"), Comment("print frquency"), 101};
	      fhicl::Atom<int> minnsh {Name("minNStrawHits"), Comment("minimum number of straw hits "),2};
	      fhicl::Atom<int> minnch {Name("minNComboHits"), Comment("number of combohits allowed"),8};
	      fhicl::Atom<TrkFitFlag> saveflag {Name("SaveTrackFlag"),Comment("if set to OK then save the track"), TrkFitFlag::helixOK};
	      fhicl::Atom<int> minNHitsTimeCluster{Name("minNHitsTimeCluster"),Comment("minium allowed time cluster"), 1 };
	      fhicl::Atom<art::InputTag> chToken{Name("ComboHitCollection"),Comment("tag for combo hit collection")};
	      fhicl::Atom<art::InputTag> shToken{Name("StrawHitCollection"),Comment("tag for straw hit collection")};
	      fhicl::Atom<art::InputTag> tcToken{Name("TimeClusterCollection"),Comment("tag for time cluster collection")};
	      fhicl::Atom<bool> DoDrift{Name("DoDrift"),Comment("turn on for drift fit")};
	      fhicl::Table<CosmicTrackFit::Config> tfit{Name("CosmicTrackFit"), Comment("fit")};
	};
	typedef art::EDProducer::Table<Config> Parameters;
	explicit CosmicTrackFinder(const Parameters& conf);
	virtual ~CosmicTrackFinder();
	virtual void beginJob() override;
	virtual void beginRun(art::Run& run) override;
	virtual void produce(art::Event& event ) override;
    
    private:
    
	Config _conf;

	int 				_debug;
	int                                 _printfreq;
	int 				_minnsh; // minimum # of strawHits in CH
	int 				_minnch; // minimum # of ComboHits for viable fit
	TrkFitFlag				_saveflag;//write tracks that satisfy these flags
	int 				_minNHitsTimeCluster; //min number of hits in a time cluster
	float				_max_seed_chi2; ///maximum chi2 allowed for seed

	art::InputTag  _chToken;
	art::InputTag  _shToken;
	art::InputTag  _tcToken;

	bool 	   _DoDrift;

	CosmicTrackFit     _tfit;

	StrawHitFlag      _dontuseflag;

	ProditionsHandle<StrawResponse> _strawResponse_h; 

	ProditionsHandle<Tracker> _alignedTracker_h;
        void     OrderHitsY(ComboHitCollection const&chcol, std::vector<StrawHitIndex> const&inputIdx, std::vector<StrawHitIndex> &outputIdxs);
        void     fillGoodHits(ComboHitCollection const& chcol, std::vector<StrawHitIndex> const& inputIdxs, std::vector<StrawHitIndex> &outputIdxs);
	int      goodHitsTimeCluster(const TimeCluster TCluster, ComboHitCollection chcol);
   
};


    CosmicTrackFinder::CosmicTrackFinder(const Parameters& conf) :
	art::EDProducer(conf),
	_debug  (conf().debug()),
	_printfreq  (conf().printfreq()),
	_minnsh   (conf().minnsh()),
	_minnch  (conf().minnch()),
	_saveflag  (conf().saveflag()),
	_minNHitsTimeCluster(conf().minNHitsTimeCluster()),
	_chToken (conf().chToken()),
        _shToken (conf().shToken()),
	_tcToken (conf().tcToken()),
	_DoDrift (conf().DoDrift()),
	_tfit (conf().tfit())
	{
		consumes<ComboHitCollection>(_chToken);
		consumes<ComboHitCollection>(_shToken);
		consumes<TimeClusterCollection>(_tcToken);
		produces<CosmicTrackSeedCollection>();
                produces<LineSeedCollection>();
	    
 	}

    CosmicTrackFinder::~CosmicTrackFinder(){}

    void CosmicTrackFinder::beginJob() {
   
	art::ServiceHandle<art::TFileService> tfs;
    }

    void CosmicTrackFinder::beginRun(art::Run& run) {
    }

    void CosmicTrackFinder::produce(art::Event& event ) {
      Tracker const& tracker = _alignedTracker_h.get(event.id());
      _tfit.setTracker(&tracker);
      StrawResponse const& srep = _strawResponse_h.get(event.id());

      if (_debug != 0) std::cout<<"Producing Cosmic Track in  Finder..."<<std::endl;
      unique_ptr<CosmicTrackSeedCollection> seed_col(new CosmicTrackSeedCollection());

      int _iev=event.id().event();
      if (_debug > 0){
        std::cout<<"ST Finder Event #"<<_iev<<std::endl;
      } 

      auto const& chH = event.getValidHandle<ComboHitCollection>(_chToken);
      const ComboHitCollection& chcol(*chH);
      auto const& shH = event.getValidHandle<ComboHitCollection>(_shToken);
      const ComboHitCollection& shcol(*shH);
      auto  const& tcH = event.getValidHandle<TimeClusterCollection>(_tcToken);
      const TimeClusterCollection& tccol(*tcH);

      unique_ptr<LineSeedCollection> lseed_col(new LineSeedCollection());

      std::vector<StrawHitIndex> panelHitIdxs;


      for (size_t index=0;index< tccol.size();++index) {
        int   nGoodTClusterHits(0);
        const auto& tclust = tccol[index];
        nGoodTClusterHits     = goodHitsTimeCluster(tclust,chcol);

        if ( nGoodTClusterHits < _minNHitsTimeCluster)         continue;
        if (_debug > 0){
          std::cout<<"time clusters "<<_iev<<std::endl;
        }
        CosmicTrackSeed initseed ;
        ComboHitCollection sortedcombohits;
        sortedcombohits.setParent(chcol.parent());
        initseed._panel_hits.setParent(chcol.parent());
        initseed._t0          = tclust._t0;
        initseed._timeCluster = art::Ptr<TimeCluster>(tcH,index);

        std::vector<StrawHitIndex> sortedShIndices;
        OrderHitsY(chcol,tclust.hits(),sortedShIndices); 

        int nFiltComboHits = 0;
        int nFiltStrawHits = 0;
        for (size_t i=0;i<sortedShIndices.size();i++){
          auto ch = chcol[sortedShIndices[i]];
          // FIXME for backwards compatibility
          sortedcombohits.push_back(ch);
          nFiltComboHits++;
          nFiltStrawHits += ch.nStrawHits();
        }


        if (_debug != 0){
          std::cout<<"#filtered SHits"<<nFiltStrawHits<<" #filter CHits "<<nFiltComboHits<<std::endl;
        }
        if (nFiltComboHits < _minnch ) 	continue;
        if (nFiltStrawHits < _minnsh)          continue;
        initseed._status.merge(TrkFitFlag::Straight);
        initseed._status.merge(TrkFitFlag::hitsOK);
        initseed._panel_hits = sortedcombohits;

        ostringstream title;
        title << "Run: " << event.id().run()
          << "  Subrun: " << event.id().subRun()
          << "  Event: " << event.id().event()<<".root";
        CosmicTrackSeed chifitseed = _tfit.BeginFit(title.str().c_str(), initseed, sortedcombohits);

        if (chifitseed._status.hasAnyProperty(TrkFitFlag::helixOK) && chifitseed._status.hasAnyProperty(TrkFitFlag::helixConverged) && chifitseed._track.converged == true ) { 
          std::vector<CosmicTrackSeed>          track_seed_vec;

          std::vector<StrawHitIndex> goodShIndices;
          fillGoodHits(chcol,sortedShIndices,goodShIndices);

          //FIXME for backwards compatibility
          for (size_t i=0;i<goodShIndices.size();i++){
            ComboHit thit = chcol[goodShIndices[i]];
            chifitseed._panel_hits.push_back(thit);
          }

          chifitseed._status.merge(TrkFitFlag::helixOK);
          if (chifitseed.status().hasAnyProperty(_saveflag)){

            std::vector<ComboHitCollection::const_iterator> chids;  
            chcol.fillComboHits(event, sortedShIndices, chids); 
            std::vector<ComboHitCollection::const_iterator> StrawLevelCHitIndices = chids;
            for (auto const& it : chids){
              chifitseed._straw_chits.push_back(it[0]);
            }

            for(size_t ich= 0; ich<chifitseed._straw_chits.size(); ich++) 					{  

              std::vector<StrawHitIndex> shitids;          	          		
              chifitseed._straw_chits.fillStrawHitIndices(event, ich, shitids);  
              chifitseed._strawHitIdxs.push_back(ich);

              for(auto const& ids : shitids){ 
                size_t    istraw   = (ids);
                TrkStrawHitSeed tshs;
                tshs._index  = istraw;
                tshs._t0 = tclust._t0;
                chifitseed._trkstrawhits.push_back(tshs); 
              }  
            }

            if( _tfit.goodTrack(chifitseed._track) == false){
              chifitseed._status.clear(TrkFitFlag::helixConverged);
              chifitseed._status.clear(TrkFitFlag::helixOK);
              continue;
            }
            ComboHitCollection tmpHits;
            if(_DoDrift){
              _tfit.DriftFit(chifitseed, srep);
              if( chifitseed._track.minuit_converged == false){
                chifitseed._status.clear( TrkFitFlag::helixConverged);
                chifitseed._status.clear(TrkFitFlag::helixOK);
                continue;
              }

              for(auto const &chit : chifitseed._straw_chits){

                if(!chit._flag.hasAnyProperty(StrawHitFlag::outlier)){

                  tmpHits.push_back(chit);
                }
              }
              chifitseed._straw_chits = tmpHits;
            }
            track_seed_vec.push_back(chifitseed);

            CosmicTrackSeedCollection* col = seed_col.get();

            if ((_DoDrift and tmpHits.size() == 0) or track_seed_vec.size() == 0)     continue;
            col->push_back(chifitseed);  

            LineSeed lseed;
            lseed._t0 = 0;
            lseed._converged = chifitseed._track.minuit_converged;
            lseed._seedSize = 0;
            lseed._timeCluster = chifitseed.timeCluster();
            std::vector<StrawHitIndex> strawHitIdxs;
            for (size_t i=0; i<chifitseed._straw_chits.size();i++){
              ComboHit const& chit0 = chifitseed._straw_chits[i];
              bool found = false;
              for (size_t j=0;j<shcol.size();j++){
                ComboHit const& chit1 = shcol[j];
                if (chit1.strawId() == chit0.strawId() && chit1.time() == chit0.time()){
                  strawHitIdxs.push_back(j);
                  found = true;
                  break;
                }
              }
              if (!found)
                std::cout << "COULD NOT FIND HIT ======================================================================================================" << std::endl;
            }
            lseed._strawHitIdxs = strawHitIdxs;
            lseed._seedInt = CLHEP::Hep3Vector(chifitseed._track.MinuitFitParams.A0,chifitseed._track.MinuitFitParams.B0,0);
            lseed._seedDir = CLHEP::Hep3Vector(chifitseed._track.MinuitFitParams.A1,chifitseed._track.MinuitFitParams.B1,1).unit();
            if (lseed._seedDir.y() > 0)
              lseed._seedDir *= -1;

            LineSeedCollection* lscol = lseed_col.get();
            lscol->push_back(lseed);

          }
        }
      }

      event.put(std::move(seed_col));    
      event.put(std::move(lseed_col));
    }

    void CosmicTrackFinder::fillGoodHits(ComboHitCollection const& chcol, std::vector<StrawHitIndex> const& inputIdxs, std::vector<StrawHitIndex> &outputIdxs){
	if (_debug != 0) {
		std::cout<<"Filling good hits..."<<std::endl;
	}
        for (size_t i=0;i<inputIdxs.size();i++){
          auto ch = &chcol[inputIdxs[i]];
          if (ch->_flag.hasAnyProperty(_dontuseflag)) continue;
          outputIdxs.push_back(inputIdxs[i]);
        }
    }

  void CosmicTrackFinder::OrderHitsY(ComboHitCollection const& chcol, std::vector<StrawHitIndex> const& inputIdxs, std::vector<StrawHitIndex> &outputIdxs){
    if (_debug != 0){
      std::cout<<"Ordering Hits..."<<std::endl;
    }

    std::vector<std::vector<ComboHit>::const_iterator> ordChColIters;
    for (size_t i=0;i<inputIdxs.size();i++){
      std::vector<ComboHit>::const_iterator tempiter = chcol.begin() + inputIdxs[i];
      ordChColIters.push_back(tempiter);
    }
    std::sort(ordChColIters.begin(), ordChColIters.end(),ycomp_iter());

    for (size_t i=0;i<ordChColIters.size();i++){
      auto thisiter = ordChColIters[i];
      outputIdxs.push_back(std::distance(chcol.begin(),thisiter));
    }
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

}
using mu2e::CosmicTrackFinder;
DEFINE_ART_MODULE(CosmicTrackFinder);

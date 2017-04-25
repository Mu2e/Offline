// $Id: FlagBkgHits_module.cc, without diagnostics $
// $Author: brownd & mpettee $ 
// $Date: 2016/11/30 $
//
// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
// Mu2e
#include "TrkReco/inc/ClusterStrawHits.hh"
#include "Mu2eUtilities/inc/MVATools.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// root 
#include "TMath.h"
#include "TH1F.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
 
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <float.h>
#include <vector>
#include <set>
#include <map>
#include <numeric>
using namespace std; 
using namespace boost::accumulators;

namespace mu2e 
{
  // structs for MVAs
  // select hits close to each other and/or to the cluster
  struct DeltaHitMVA {
    std::vector<Double_t> _pars;
    Double_t& _dphi; // phi diff to cluster
    Double_t& _drho; // rho diff to cluster
    Double_t& _dt;  // time diff to cluster
    DeltaHitMVA() : _pars(3,0.0),_dphi(_pars[0]),_drho(_pars[1]),_dt(_pars[2]){}
  };

  struct DeltaClusterMVA {
    std::vector<Double_t> _pars;
    Double_t& _prho;
    Double_t& _srho;
    Double_t& _zmin;
    Double_t& _zmax;
    Double_t& _zgap;
    Double_t& _ns;
    Double_t& _nsmiss;
    Double_t& _sphi;
    Double_t& _ngdhits;
    DeltaClusterMVA() : _pars(9,0.0),_prho(_pars[0]),_srho(_pars[1]),
    _zmin(_pars[2]),_zmax(_pars[3]),_zgap(_pars[4]),
    _ns(_pars[5]),_nsmiss(_pars[6]),
    _sphi(_pars[7]),_ngdhits(_pars[8]) {}
  };

  struct DeltaHitInfo {
    DeltaHitInfo(size_t hindex,double htime, double hphi, double hrho, double hz,double hgd,int hflag,bool stereo) :
      _hindex(hindex), _htime(htime), _hphi(hphi), _hrho(hrho), _hz(hz), _hgd(hgd), _hflag(hflag), _stereo(stereo) {}
    size_t _hindex;
    Float_t _htime, _hphi, _hrho, _hz, _hgd, _cdist;
    Int_t _hflag;
    Bool_t _stereo;
  };

  struct DeltaInfo {
    DeltaInfo() : _tcluster(0.0), _pcluster(0.0), _rcluster(0.0),
      _tmed(0.0), _pmed(0.0),_rmed(0.0),
      _tmean(0.0), _pmean(0.0), _rmean(0.0),
      _trms(0.0), _prms(0.0),_rrms(0.0),
      _nchits(0),_ngdhits(0), _ngdstereo(0),
      _zmin(0.0),_zmax(0.0),_zgap(0.0),
      _ismin(-1),_ismax(-1),_ns(-1),_nsmiss(-1),
      _idmin(-1),_idmax(-1),_nd(-1),_ndmiss(-1),
      _isdelta(false),_pmvaout(-100.0){}
    // cluster information
    double _tcluster, _pcluster, _rcluster;
    // accumulators for cluster center
    double _tmed, _pmed, _rmed;
    double _tmean, _pmean, _rmean;
    double _trms, _prms, _rrms;
    // summary information
    unsigned _nchits, _ngdhits, _ngdstereo;
    double _zmin, _zmax, _zgap;
    unsigned _ismin, _ismax, _ns, _nsmiss;
    unsigned _idmin, _idmax, _nd, _ndmiss;
    bool _isdelta;
    double _pmvaout;
    // information about the hits in the delta
    std::vector<DeltaHitInfo> _dhinfo;
  };

  class FlagBkgHits : public art::EDProducer
  {
    public:
      explicit FlagBkgHits(fhicl::ParameterSet const&);
      virtual ~FlagBkgHits();
      virtual void beginJob();
      virtual void beginRun(art::Run&);
      virtual void produce(art::Event& event ); 
      void endJob();
    private:
      // configuration parameters
      int _debug;
      int _printfreq;
      // event object labels
      art::InputTag _shlabel, _shplabel, _stlabel, _shflabel;
      // straw hit selection masks
      StrawHitFlag _stmask, _deltamask, _ismask;
      // delta-ray removal parameters
      art::InputTag _dhittype, _stereohitweights, _nonstereohitweights;
      art::InputTag _dclustertype, _nonstereoclusterweights, _stereoclusterweights;
      double _gdstereo, _gdnonstereo;
      bool _flagall;
      unsigned _mindh, _minns;
      double _stereoclusterhitfrac;
      double _stereoclustermvacut, _nonstereoclustermvacut;
      // input collections
      const StrawHitCollection* _shcol;
      const StrawHitPositionCollection* _shpcol;
      const StereoHitCollection* _stcol;
      const StrawHitFlagCollection* _shfcol;
      // output collections
      StrawHitFlagCollection* _bkgfcol;
      // internal helper functions
      void filterDeltas();
      void fillDeltaInfo(std::list<StrawHitCluster> const& clusters,std::vector<DeltaInfo>& dinfo);
      void fillDeltaSummary(DeltaInfo& delta);
      double deltaPhi(double phi1,double phi2) const;
      bool findData(const art::Event& evt);
      // clusterer
      ClusterStrawHits _clusterer;
      // delta removal diagnostics
      MVATools _stereohitMVA, _nonstereohitMVA; // assign hits to delta clusters
      DeltaHitMVA _dhmva; // input variables to TMVA for delta hit selection
      DeltaClusterMVA _dpmva;
      MVATools _stereoclusterMVA, _nonstereoclusterMVA; // classify delta clusters
      Int_t _iev;
      Bool_t _isdelta;
      Float_t _pphi, _pt, _prho;
      Float_t _zmin, _zmax, _zgap;
      Int_t _ns, _nd, _smin, _smax, _nsmiss, _ndmiss;
      Int_t _nchits,_ngdhits, _ngdstereo;
      Float_t _tmed, _pmed, _rmed;
      Float_t _tmean, _pmean, _rmean;
      Float_t _stime, _sphi, _srho;
      Float_t _pmvaout;
      // histograms
      TH1F *_nhits, *_niter, *_nchanged, *_nclusters;
  };

  FlagBkgHits::FlagBkgHits(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>("debugLevel",0)),
    _printfreq(pset.get<int>("printFrequency",101)),
    _shlabel(pset.get<art::InputTag>("StrawHitCollectionLabel","makeSH")),
    _shplabel(pset.get<art::InputTag>("StrawHitPositionCollectionLabel","MakeStereoHits")),
    _stlabel(pset.get<art::InputTag>("StereoHitCollectionLabel","MakeStereoHits")),
    _shflabel(pset.get<art::InputTag>("StrawHitFlagCollectionLabel","FlagStrawHits")),
    _stmask(StrawHitFlag::stereo),
    _deltamask(StrawHitFlag::delta),
    _ismask(StrawHitFlag::isolated),
    _gdstereo(pset.get<double>("StereoHitMVACut",0.55)),
    _gdnonstereo(pset.get<double>("NonStereoHitMVACut",0.55)),
    _flagall(pset.get<bool>("FlagAllHits",false)), // flag all hits in the cluster, regardless of MVA value
    _mindh(pset.get<unsigned>("MinDeltaHits",5)),
    _minns(pset.get<unsigned>("MinNStations",2)),
    _stereoclusterhitfrac(pset.get<double>("StereoClusterHitFraction",0.5)),
    _stereoclustermvacut(pset.get<double>("StereoClusterMVACut",0.8)),
    _nonstereoclustermvacut(pset.get<double>("NonStereoClusterMVACut",0.8)),
    _clusterer(pset.get<fhicl::ParameterSet>("ClusterStrawHits",fhicl::ParameterSet())),
    _stereohitMVA(pset.get<fhicl::ParameterSet>("StereoHitMVA",fhicl::ParameterSet())),
    _nonstereohitMVA(pset.get<fhicl::ParameterSet>("NonStereoHitMVA",fhicl::ParameterSet())),
    _stereoclusterMVA(pset.get<fhicl::ParameterSet>("StereoClusterMVA",fhicl::ParameterSet())),
    _nonstereoclusterMVA(pset.get<fhicl::ParameterSet>("NonStereoClusterMVA",fhicl::ParameterSet()))
  {
    produces<StrawHitFlagCollection>();
    // eventually, should also produce a collection of delta-rays and their properties, FIXME!!
  }

  FlagBkgHits::~FlagBkgHits(){}

  void FlagBkgHits::beginJob(){
    if(_debug > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _nhits=tfs->make<TH1F>("nhits","N Hits Used in Clustering",100,0,5000);
      _niter=tfs->make<TH1F>("niter","N Cluster Iterations",50,-0.5,49.5);
      _nchanged=tfs->make<TH1F>("nchanged","N Hits Changed in Last Iteration",100,-0.5,99.5);
      _nclusters=tfs->make<TH1F>("nclusters","N Clusters",200,0,1000);
    }
    // initialize the MVAs
    _stereohitMVA.initMVA();
    _nonstereohitMVA.initMVA();
    _stereoclusterMVA.initMVA();
    _nonstereoclusterMVA.initMVA();
     if(_debug > 0){
      cout << "Stereo Hit MVA : " << endl;
      _stereohitMVA.showMVA();
      cout << "Non-Stereo Hit MVA : " << endl;
      _nonstereohitMVA.showMVA();
      cout << "Stereo Cluster MVA : " << endl;
      _stereoclusterMVA.showMVA();
      cout << "Non-Stereo Cluster MVA : " << endl;
      _nonstereoclusterMVA.showMVA();
     }
  }

  void FlagBkgHits::beginRun(art::Run& ){}

  void FlagBkgHits::produce(art::Event& event ) {
    // event printout
    _iev=event.id().event();
    if(_debug > 0 && (_iev%_printfreq)==0)cout<<"FlagBkgHits: event="<<_iev<<endl;
    // find the data
    if(!findData(event)){
      throw cet::exception("RECO")<< "Missing input collection" << endl;
    }
    // create output
    unique_ptr<StrawHitFlagCollection> bkgfcol(new StrawHitFlagCollection);
    _bkgfcol = bkgfcol.get();
    // find clusters in time/phi/rho space
    std::vector<DeltaInfo> dinfo;
    // test of new straw hit clustering
    StrawHitClusterList clusters;
    _clusterer.findClusters(*_shcol,*_shpcol,*bkgfcol,clusters);
    if(_debug > 0){
      _nhits->Fill(clusters._nhits);
      _niter->Fill(clusters._niter);
      _nchanged->Fill(clusters._nchanged);
      _nclusters->Fill(clusters._clist.size());
    }
    fillDeltaInfo(clusters._clist,dinfo);
    // loop over clusters and classify them
    for(auto& delta : dinfo){
    // fill summary information
      fillDeltaSummary(delta);
    // if the cluster is a delta, flag its hits
    // (or optionally flag all the hits in the cluster)
      if(delta._isdelta){
	for(auto dhinfo : delta._dhinfo){
	  if(_flagall || dhinfo._hflag > 0 ){
	    bkgfcol->at(dhinfo._hindex).merge(_deltamask);
	  }
	}
      }
    // flag single-hit clusters as 'isolated'
      if(delta._dhinfo.size()==1)
	bkgfcol->at(delta._dhinfo[0]._hindex).merge(_ismask);
    }
    // put the background flag into the event
     event.put(std::move(bkgfcol));
  }

  void FlagBkgHits::endJob(){
    // does this cause the file to close?
    art::ServiceHandle<art::TFileService> tfs;
  }

  // find the input data objects 
  bool FlagBkgHits::findData(const art::Event& evt){
    bool retval(false);
    _shcol = 0; _shpcol = 0; _stcol = 0; _shfcol = 0;
    art::Handle<mu2e::StrawHitCollection> shH;
    if(evt.getByLabel(_shlabel,shH))
      _shcol = shH.product();
    art::Handle<mu2e::StrawHitPositionCollection> shposH;
    if(evt.getByLabel(_shplabel,shposH))
      _shpcol = shposH.product();
    art::Handle<mu2e::StereoHitCollection> stH;
    if(evt.getByLabel(_stlabel,stH))
      _stcol = stH.product();
    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if(evt.getByLabel(_shflabel,shflagH))
      _shfcol = shflagH.product();
    // don't require stereo hits, they are used only for diagnostics
    retval = _shcol != 0 && _shpcol != 0 && _shfcol != 0;

    return retval;
  }

  void FlagBkgHits::fillDeltaInfo(std::list<StrawHitCluster> const& clusters,std::vector<DeltaInfo>& dinfo){
    for (std::list<StrawHitCluster>::const_iterator icl=clusters.begin();icl!=clusters.end();++icl){
      StrawHitCluster const& cluster = *icl;
      DeltaInfo dp;
      dp._tcluster = dp._tmed = cluster.time();
      dp._pcluster = dp._pmed = cluster.pos().phi();
      dp._rcluster = dp._rmed = cluster.pos().perp();
      // find precise cluster position
      for(auto ich : cluster.hits()){
	    size_t ish = ich._index;
        StrawHit const& sh = _shcol->at(ish);
	    StrawHitPosition const& shp = _shpcol->at(ish);
	    double ct = sh.time();
	    double phi = dp._pcluster+deltaPhi(shp.pos().phi(),dp._pcluster);
	    double rho = shp.pos().perp();
	    double dphi = deltaPhi(phi,dp._pcluster);
	    _dhmva._dphi = fabs(dphi);
	    _dhmva._drho = fabs(rho - dp._rmed);
	    _dhmva._dt = fabs(ct - dp._tmed);
	    double gd(-1.0);
	    int iflag(0);
	    bool stereo = shp.flag().hasAllProperties(_stmask);
	    if(stereo){
	      gd = _stereohitMVA.evalMVA(_dhmva._pars);
	      if(gd > _gdstereo)iflag = 1;
    	} else {
    	  gd = _nonstereohitMVA.evalMVA(_dhmva._pars);
        if(gd > _gdnonstereo)iflag = 1;
    	}
    	DeltaHitInfo dhinfo(ish,ct,phi,rho,shp.pos().z(),gd,iflag,stereo);
        dhinfo._cdist = ich._dist;
    	dp._dhinfo.push_back(dhinfo);
      }
        fillDeltaSummary(dp);
        dinfo.push_back(dp);
    }
  }

  void FlagBkgHits::fillDeltaSummary(DeltaInfo& delta) {
    // tracker, to get StrawID later
    const TTracker& tracker = dynamic_cast<const TTracker&>(getTrackerOrThrow());
    unsigned nplanes = tracker.nPlanes();
    unsigned nstations = nplanes/2;
    // fill the BDT information about each hit, compared to the centroid
    // also accumulate information about selected hits
    accumulator_set<double, stats<tag::mean,tag::variance(lazy)> > tacc,pacc,racc;
    std::vector<bool> planes(nplanes,false);
    std::vector<bool> stations(nstations,false);
    std::vector<double> hz;
    hz.reserve(delta._dhinfo.size());
    delta._ngdhits = delta._ngdstereo = 0;
    delta._nchits = delta._dhinfo.size();
    for(auto ih : delta._dhinfo){
      DeltaHitInfo& dhinfo = ih;
      if(dhinfo._hflag > 0){
        tacc(dhinfo._htime);
        pacc(dhinfo._hphi);
        racc(dhinfo._hrho);
	const StrawHit& sh = _shcol->at(dhinfo._hindex);
	unsigned iplane = (unsigned)(tracker.getStraw(sh.strawIndex()).id().getPlaneId());
	unsigned istation = iplane/2;
	planes[iplane] = true;
	stations[istation] = true;
	hz.push_back(dhinfo._hz);
	++delta._ngdhits;
	if(dhinfo._stereo)
	  ++delta._ngdstereo;
      }
    }
    if(extract_result<tag::count>(tacc) > 0){
      delta._tmean = extract_result<tag::mean>(tacc);
      delta._pmean = extract_result<tag::mean>(pacc);
      delta._rmean = extract_result<tag::mean>(racc);
      delta._trms = extract_result<tag::variance>(tacc);
      delta._prms = extract_result<tag::variance>(pacc);
      delta._rrms = extract_result<tag::variance>(racc);
      delta._trms = (delta._trms>0.0) ? sqrt(delta._trms) : 0.0;
      delta._prms = (delta._prms>0.0) ? sqrt(delta._prms) : 0.0;
      delta._rrms = (delta._rrms>0.0) ? sqrt(delta._rrms) : 0.0;
      std::sort(hz.begin(),hz.end());
      delta._zmin = hz.front();
      delta._zmax = hz.back();
      // look for gaps in z
      delta._zgap = 0.0;
      for(unsigned iz=1;iz<hz.size();++iz){
	if(hz[iz]-hz[iz-1] > delta._zgap)delta._zgap = hz[iz]-hz[iz-1]; 
      }
      // count 'missing' planes between first and last
      delta._ismin = 0;
      delta._ismax = nplanes/2-1;
      delta._idmin = 0;
      delta._idmax = nplanes-1;
      delta._nsmiss = 0;
      delta._ndmiss = 0;
      while(!stations[delta._ismin])++delta._ismin;
      while(!stations[delta._ismax])--delta._ismax;
      delta._ns = delta._ismax-delta._ismin+1;
      for(unsigned is =delta._ismin;is<delta._ismax;++is){
	if(!stations[is])++delta._nsmiss;
      }
      while(!planes[delta._idmin])++delta._idmin;
      while(!planes[delta._idmax])--delta._idmax;
      delta._nd = delta._idmax-delta._idmin+1;
      for(unsigned id =delta._idmin;id<delta._idmax;++id){
	if(!planes[id])++delta._ndmiss;
      }
      // compute final delta selection based on the delta properties
      delta._isdelta = false;
      _dpmva._prho = delta._rmed;
      _dpmva._srho = delta._rrms;
      _dpmva._zmin = delta._zmin;
      _dpmva._zmax = delta._zmax;
      _dpmva._zgap = delta._zgap;
      _dpmva._ns = delta._ns;
      _dpmva._nsmiss = delta._nsmiss;
      _dpmva._sphi = delta._prms;
      _dpmva._ngdhits = delta._ngdhits;
      double sfrac = delta._ngdstereo/delta._ngdhits;
      bool clusterdelta;
      if(sfrac > _stereoclusterhitfrac){
	delta._pmvaout = _stereoclusterMVA.evalMVA(_dpmva._pars);
	clusterdelta = delta._pmvaout > _stereoclustermvacut;
      } else {
	delta._pmvaout = _nonstereoclusterMVA.evalMVA(_dpmva._pars);
	clusterdelta = delta._pmvaout > _nonstereoclustermvacut;
      }

      delta._isdelta = delta._ngdhits >= _mindh && delta._ns >= _minns && clusterdelta;
    }
  }

 double FlagBkgHits::deltaPhi(double phi1,double phi2) const {
    double dphi = phi2-phi1;
    if(dphi > M_PI){
      dphi -= 2*M_PI;
    } else if(dphi < -M_PI){
      dphi += 2*M_PI;
    }
    return dphi;
  }
    
}
using mu2e::FlagBkgHits;


DEFINE_ART_MODULE(FlagBkgHits);

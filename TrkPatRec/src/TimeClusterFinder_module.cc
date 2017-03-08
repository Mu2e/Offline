//
// TTracker time peak finder
//
// $Id: TimeClusterFinder_module.cc,v 1.3 2014/08/25 12:08:29 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/25 12:08:29 $
//
// Original author D. Brown and G. Tassielli
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
// Mu2e
#include "GeneralUtilities/inc/Angles.hh"
#include "Mu2eUtilities/inc/MVATools.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
// tracking
#include "TrkReco/inc/TrkUtilities.hh"
#include "TrkReco/inc/TrkTimeCalculator.hh"
// root
#include "TH1F.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp>
// C++
#include <memory>
#include <algorithm>
#include <utility>
using namespace std; 
using namespace boost::accumulators;
using CLHEP::Hep3Vector;

namespace mu2e {
// struct to run MVA
  struct TimePeakMVA {
    vector<Double_t> _pars;
    Double_t& _dt;
    Double_t& _dphi;
    Double_t& _rho;
    TimePeakMVA() : _pars(3,0.0), _dt(_pars[0]), _dphi(_pars[1]), _rho(_pars[2]) {}
  };

  // structs to facilitate peak finding
  typedef pair<Float_t,int> BinContent;
  struct BCComp : public binary_function<BinContent, BinContent, bool> {
     bool operator()(BinContent const& x, BinContent const& y) {
       return x.first > y.first; // sort descending
    }
  };

  class TimeClusterFinder : public art::EDProducer {
  public:
    enum ClusterAlgorithm {peak=0,scan};
    explicit TimeClusterFinder(fhicl::ParameterSet const& pset);
    virtual ~TimeClusterFinder();

    virtual void beginJob();

    // This is called for each event.
    void produce(art::Event & e);

  private:
// configuration parameters
    int		  _iev;
    int           _debug;
    int           _printfreq;
    ClusterAlgorithm _algo;
    // event object labels
    art::InputTag			_shTag;
    art::InputTag			_shpTag;
    art::InputTag			_shfTag;
    art::InputTag			_ccTag;
    
    StrawHitFlag  _hsel, _hbkg;
    double        _maxdt;

    // time spectrum parameters
    unsigned      _minnhits;
    double        _minpeakmva, _maxpeakdt;
    double        _tmin;
    double        _tmax;
    double        _tbin;
    TH1F	  _timespec;
    double        _ymin;
    bool          _refine;
    bool	  _usecc;
    double	  _ccmine, _ccwt;
    double	  _maxover;

    // cache of event objects
    const StrawHitCollection*             _shcol;
    const StrawHitFlagCollection*         _shfcol;
    const StrawHitPositionCollection*     _shpcol;
    const CaloClusterCollection*	  _cccol;
    // need the calo handle to so I can make a ptr later
    MVATools                              _peakMVA; // MVA for peak cleaning
    TimePeakMVA                           _pmva; // input variables to TMVA for peak cleaning

    // T0 calcuator
    TrkTimeCalculator _ttcalc;

    art::Handle<mu2e::StrawHitCollection> _strawhitsH;
// helper functions
    void initCluster	  (TimeCluster& tp); // fill peak information from the list of hits
    void refineCluster	  (TimeCluster& tp); // refine the peak information and hit list
    bool findData         (const art::Event& evt);
    void findClusters    (TimeClusterCollection*, art::Event const& evt);
    void findPeaks    (std::vector<double>& tctimes);
    void scanPeaks    (std::vector<double>& tctimes);
    void fillTimeSpectrum();
    void addCaloCluster(TimeCluster& tc,art::Event const& e);
    bool goodHit(StrawHitFlag const& flag) const;
    bool goodCaloCluster(CaloCluster const& cc) const;
    double clusterWeight(TimeCluster const& tc) const;

  };

  TimeClusterFinder::~TimeClusterFinder() {
  }
  
  TimeClusterFinder::TimeClusterFinder(fhicl::ParameterSet const& pset) :
    _debug             (pset.get<int>("debugLevel",0)),
    _printfreq         (pset.get<int>("printFrequency",101)),
    _algo	       (static_cast<ClusterAlgorithm>(pset.get<int>("ClusterAlgorithm",peak))),
    _shTag	       (pset.get<art::InputTag>("StrawHitCollection","makeSH")),
    _shpTag	       (pset.get<art::InputTag>("StrawHitPositionCollection","MakeStereoHits")),
    _shfTag	       (pset.get<art::InputTag>("StrawHitFlagCollection","FlagBkgHits")),
    _ccTag             (pset.get<art::InputTag>("caloClusterModuleLabel","MakeCaloCluster")),
    _hsel		(pset.get<std::vector<std::string> >("HitSelectionBits",vector<string>{"EnergySelection","TimeSelection","RadiusSelection"})),
    _hbkg              (pset.get<vector<string> >("HitBackgroundBits",vector<string>{"DeltaRay","Isolated"})),
    _maxdt             (pset.get<double>("DtMax",30.0)),
    _minnhits          (pset.get<unsigned>("MinNHits",10)),
    _minpeakmva        (pset.get<double>("MinTimePeakMVA",0.2)),
    _maxpeakdt         (pset.get<double>("MaxTimePeakDeltat",25.0)),
    _tmin              (pset.get<double>("tmin",500.0)),
    _tmax              (pset.get<double>("tmax",1700.0)),
    _tbin              (pset.get<double>("tbin",15.0)),
    _ymin              (pset.get<double>("ymin",5.0)),
    _refine	       (pset.get<bool>("RefineClusters",true)), // refine content of clusters
    _usecc	       (pset.get<bool>("UseCaloCluster",false)),
    _ccmine            (pset.get<double>("CaloClusterMinE",50.0)), // minimum energy to call a cluster 'good' (MeV)
    _ccwt              (pset.get<double>("CaloClusterWeight",10.0)), // Cluster weight in units of hits
    _maxover           (pset.get<double>("MaxOverlap",0.3)), // Maximum hit overlap to consider clusters as different
    _peakMVA           (pset.get<fhicl::ParameterSet>("PeakCleanMVA",fhicl::ParameterSet())),
    _ttcalc            (pset.get<fhicl::ParameterSet>("T0Calculator",fhicl::ParameterSet()))
  {
    // set # bins for time spectrum plot
    unsigned nbins = (unsigned)rint((_tmax-_tmin)/_tbin);
    _timespec = TH1F("timespec","time spectrum",nbins,_tmin,_tmax);
    // Tell the framework what we make
    produces<TimeClusterCollection>();
    produces<StrawHitFlagCollection>();
  }

  void TimeClusterFinder::beginJob(){
    _peakMVA.initMVA();
    if(_debug > 0){
      cout << "TimeClusterFinder MVA : " << endl; 
      _peakMVA.showMVA();
    }
  }

  void TimeClusterFinder::produce(art::Event & event ) {
    // debug printout
    _iev=event.id().event();
    if(_debug > 0 && (_iev%_printfreq)==0)cout<<"TimeClusterFinder: event="<<_iev<<endl;
    // find the data
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::TimeClusterFinder: data missing or incomplete"<< endl;
      return;
    }
    // create time peak collection
    std::unique_ptr<TimeClusterCollection>    tccol  (new TimeClusterCollection);
    // copy in the existing flags
    std::unique_ptr<StrawHitFlagCollection> flags(new StrawHitFlagCollection(*_shfcol));
    // find the clusters
    findClusters(tccol.get(),event);
    // set the flag for all hits associated to a time peak
    if(_debug > 0) cout << "Found " << tccol->size() << " Time Clusters " << endl;
    for (auto tpc : *tccol.get()){
      if(_debug > 1){
	cout << "Time Cluster time = " << tpc.t0().t0() << " +- " << tpc.t0().t0Err()
	<< " position = " << tpc._pos << endl;
      }
      for (auto shi : tpc._strawHitIdxs ) {
	flags->at(shi).merge(StrawHitFlag::tclust);
	if(_debug > 3){
	  cout << "Time Cluster hit at index " << shi << endl;
	}
      }
    }
    // put collections into the event
    event.put(std::move( tccol ));
    event.put(std::move( flags));
  }

   // find the input data objects
  bool TimeClusterFinder::findData(const art::Event& evt){
    _shcol  = 0;
    _shfcol = 0;
    _shpcol = 0;
    _cccol  = 0;

    auto shH = evt.getValidHandle<StrawHitCollection>(_shTag);
    _shcol = shH.product();
    auto shpH = evt.getValidHandle<StrawHitPositionCollection>(_shpTag);
    _shpcol = shpH.product();
    auto shfH = evt.getValidHandle<StrawHitFlagCollection>(_shfTag);
    _shfcol = shfH.product();

    if(_usecc){
      auto ccH = evt.getValidHandle<CaloClusterCollection>(_ccTag);
      _cccol = ccH.product();
    }
    return _shcol != 0 && _shfcol != 0 && _shpcol != 0 && (!_usecc || _cccol != 0);
  }

  void TimeClusterFinder::findClusters(TimeClusterCollection* tclusts,art::Event const& evt) {
    // first, fill a histogram with all the hit times
    fillTimeSpectrum();
  // find cluster seeds instead this spectrum
    vector<double> tctimes;
    switch (_algo ) {
      case peak : default:
	findPeaks(tctimes);
	break;
      case scan :
	scanPeaks(tctimes);
	break;
    }
// loop over the seeds and create time clusters from them
    for(auto tctime: tctimes) {
      if(_debug > 1)
	cout << "Peak Time " << tctime  << endl;
      TimeCluster tclust;
      // initial t0 is the peak position
      tclust._t0 = TrkT0(tctime,1.0);
     // associate all hits in the time window with this peak
      for(size_t istr=0; istr<_shcol->size(); ++istr){
	if(goodHit(_shfcol->at(istr))){
	  double time = _ttcalc.strawHitTime(_shcol->at(istr),_shpcol->at(istr));
	  if(fabs(time-tctime) < _maxdt){
	    tclust._strawHitIdxs.push_back(StrawHitIndex(istr));
	  }
	}
      }
      // if we're using the calo clusters, associate those too
      if(_usecc) addCaloCluster(tclust,evt);
      // initialize the t0 value using these hits
      initCluster(tclust);
      // refine the hit content and t0 value
      if(_refine)refineCluster(tclust);
      // final check on # of hits
      if(tclust._strawHitIdxs.size() >= _minnhits){
	bool overl(false);
	// compare this to the previous cluster: if they overlap too much, just take the best one
	for(auto itclust = tclusts->begin(); itclust < tclusts->end(); ++itclust) {
	  if(TrkUtilities::overlap(tclust,*itclust) > _maxover){
	    if(clusterWeight(tclust) > clusterWeight(*itclust)){ // new cluster is 'better', erase the old
	      tclusts->erase(itclust);
	      break;
	    } else {
	    // otherwise mark this cluster as overlapping and stop
	      overl = true;
	      break;
	    }
	  }
	}
	if(!overl)tclusts->push_back(tclust);
      }
    }
    // debug test of histogram
    if(_debug > 2){
      art::ServiceHandle<art::TFileService> tfs;
      TH1F* tspec = tfs->make<TH1F>(_timespec);
      char name[40];
      char title[100];
      snprintf(name,40,"tspec_%i",_iev);
      snprintf(title,100,"time spectrum event %i;nsec",_iev);
      tspec->SetNameTitle(name,title);
    }
  }

  void TimeClusterFinder::initCluster(TimeCluster& tclust) {
    // use medians to initialize robustly
    accumulator_set<double, stats<tag::min > > tmin;
    accumulator_set<double, stats<tag::max > > tmax;
    accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > tacc;
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > xacc;
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > yacc;
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > zacc;
    unsigned nstrs = tclust._strawHitIdxs.size();
    for(auto ish :tclust._strawHitIdxs) { 
      if(goodHit(_shfcol->at(ish))){
      // compute the corrected hit time
	Hep3Vector const& pos = _shpcol->at(ish).pos();
	double htime = _ttcalc.strawHitTime(_shcol->at(ish),_shpcol->at(ish));
	// weight inversely by the hit time resolution
	double wt = std::pow(1.0/_ttcalc.strawHitTimeErr(),2);
	tmin(htime);
	tmax(htime);
	tacc(htime,weight=wt);
	xacc(pos.x());
	yacc(pos.y());
	zacc(pos.z());
      }
    }
    // add cluster time
    if(tclust._caloCluster.isNonnull()){
      double ctime = _ttcalc.caloClusterTime(*tclust._caloCluster);
      double wt = std::pow(1.0/_ttcalc.caloClusterTimeErr(tclust._caloCluster->diskId()),2);
      tacc(ctime,weight=wt);
    }
    // set peak info
    static double invsqrt12(1.0/sqrt(12.0));
    tclust._t0._t0 = extract_result<tag::weighted_median>(tacc);
    tclust._t0._t0err = ( boost::accumulators::extract::max(tmax)-boost::accumulators::extract::min(tmin))*invsqrt12/sqrt(nstrs);
    tclust._pos = Hep3Vector(median(xacc),median(yacc),median(zacc));
  }

  void TimeClusterFinder::refineCluster(TimeCluster& tclust) {
    // iteratively filter the worst outlier
    double worstmva(100.0);
    double pphi(tclust._pos.phi());
    double ptime(tclust._t0.t0());
    do{
      // Use the MVA to find the least signal-like hit
      size_t iworst =0;
      worstmva = 100.0;
      for(size_t ips=0;ips<tclust._strawHitIdxs.size();++ips){
	unsigned ish = tclust._strawHitIdxs[ips];
	double dt = _ttcalc.strawHitTime(_shcol->at(ish),_shpcol->at(ish)) - ptime;
	double rho = _shpcol->at(ish).pos().perp();
	double phi = _shpcol->at(ish).pos().phi();
	double dphi = Angles::deltaPhi(phi,pphi);
	// compute MVA
	_pmva._dt = dt;
	_pmva._dphi = dphi;
	_pmva._rho = rho;
	double mvaout = _peakMVA.evalMVA(_pmva._pars);
	if(mvaout < worstmva){
	  worstmva = mvaout;
	  iworst = ips;
	}
      }
      // remove the worst hit
      if(worstmva < _minpeakmva){
	std::swap(tclust._strawHitIdxs[iworst],tclust._strawHitIdxs.back());
	tclust._strawHitIdxs.pop_back();
      }
      // re-compute the average phi and range
      accumulator_set<double, stats<tag::mean > > facc;
      accumulator_set<double, stats<tag::weighted_mean >, double > tacc;
      for(size_t ips=0;ips<tclust._strawHitIdxs.size();++ips){
	unsigned ish = tclust._strawHitIdxs[ips];
	double time = _ttcalc.strawHitTime(_shcol->at(ish),_shpcol->at(ish));
	double wt = std::pow(1.0/_ttcalc.strawHitTimeErr(),2);
	double phi = _shpcol->at(ish).pos().phi();
	Angles::deltaPhi(phi,pphi);
	tacc(time,weight=wt);
	facc(phi);
      }
      // add cluster time. 
      if(tclust._caloCluster.isNonnull()){
	double time = _ttcalc.caloClusterTime(*tclust._caloCluster);
	double wt = std::pow(1.0/_ttcalc.caloClusterTimeErr(tclust._caloCluster->diskId()),2);
	tacc(time,weight=wt);
      }
      pphi = extract_result<tag::mean>(facc);
      ptime = extract_result<tag::weighted_mean>(tacc);
    } while(tclust._strawHitIdxs.size() >= _minnhits && worstmva < _minpeakmva);
    // final pass: hard cut on dt 
    vector<size_t> toremove;
    accumulator_set<double, stats<tag::mean > > facc;
    accumulator_set<double, stats<tag::weighted_variance(lazy)>, double > terr;
    accumulator_set<double, stats<tag::mean > > racc;
    accumulator_set<double, stats<tag::mean > > zacc;
    for(size_t ips=0;ips<tclust._strawHitIdxs.size();++ips){
      unsigned ish = tclust._strawHitIdxs[ips];
      double dt = _ttcalc.strawHitTime(_shcol->at(ish),_shpcol->at(ish)) - ptime;
      double wt = std::pow(1.0/_ttcalc.strawHitTimeErr(),2);
      double phi = _shpcol->at(ish).pos().phi();
      double rho = _shpcol->at(ish).pos().perp();
      Angles::deltaPhi(phi,pphi);
      if(fabs(dt) < _maxpeakdt){
	terr(_ttcalc.strawHitTime(_shcol->at(ish),_shpcol->at(ish)),weight=wt);
	facc(phi);
	racc(rho);
	zacc(_shpcol->at(ish).pos().z());
      } else {
	toremove.push_back(ips);
      }
    }
    // add cluster time
    if(tclust._caloCluster.isNonnull()){
      double time = _ttcalc.caloClusterTime(*tclust._caloCluster);
      double wt = std::pow(1.0/_ttcalc.caloClusterTimeErr(tclust._caloCluster->diskId()),2);
      terr(time,weight=wt);
    }
    // actually remove the hits; must start from the back
    std::sort(toremove.begin(),toremove.end(),std::greater<size_t>());
    for(auto irm=toremove.begin();irm!=toremove.end();++irm){
      std::swap(tclust._strawHitIdxs[*irm],tclust._strawHitIdxs.back());
      tclust._strawHitIdxs.pop_back();
    }
    // update peak properties
    tclust._t0._t0 = extract_result<tag::weighted_mean>(terr);
    tclust._t0._t0err = sqrt(std::max(0.0,extract_result<tag::weighted_variance(lazy)>(terr))/extract_result<tag::count>(terr));
    pphi = extract_result<tag::mean>(facc);
    double prho = extract_result<tag::mean>(racc);
    double zpos = extract_result<tag::mean>(zacc);
    // update position
    tclust._pos = Hep3Vector(prho*cos(pphi),prho*sin(pphi),zpos);
  }

  void TimeClusterFinder::findPeaks(std::vector<double>& tctimes) {
    tctimes.clear();
    // repack the histogram into bin,content
    vector<BinContent> bcv;
    int nbins = _timespec.GetNbinsX()+1;
    vector<bool> used(nbins,false);
    for(int ibin=1;ibin < nbins; ++ibin)
      bcv.push_back(make_pair(_timespec.GetBinContent(ibin),ibin));
// sort by bin contents
    sort(bcv.begin(),bcv.end(),BCComp());
// loop over sorted contents; this starts with the highest bins and works down
    for(auto const& bc : bcv){
      if((!used[bc.second]) && bc.first >= _ymin){
	//local maximum.  Average around a few bins
	double tctime(0.0);
	double norm(0.0);
	for(int ibin = std::max(1,bc.second-2);ibin < std::min(nbins,bc.second+3); ++ibin){
	  if(fabs(_timespec.GetBinCenter(bc.second)-_timespec.GetBinCenter(ibin)) < _maxdt ){
	    norm += _timespec.GetBinContent(ibin);
	    tctime += _timespec.GetBinCenter(ibin)*_timespec.GetBinContent(ibin);
	    // mark these bins as used so they aren't found as separate peaks
	    used[ibin] = true;
	  }
	}
	tctime /= norm;
// check if sum is above threshold
	if(norm > _minnhits)
	  tctimes.push_back(tctime);
      }
    }
  }

  void TimeClusterFinder::scanPeaks(std::vector<double>& tctimes) {
    // implement this FIXME!!
  }

  bool TimeClusterFinder::goodHit(StrawHitFlag const& flag) const {
    return flag.hasAllProperties(_hsel) && !flag.hasAnyProperty(_hbkg);
  }

  bool TimeClusterFinder::goodCaloCluster(CaloCluster const& cc) const {
    return cc.energyDep() > _ccmine;
  }

  // give a weight to a time cluster for comparison purposes.
  double TimeClusterFinder::clusterWeight(TimeCluster const& tc) const {
    double wt = tc._strawHitIdxs.size();
    if(tc._caloCluster.isNonnull())wt += _ccwt;
    return wt;
  }
  
  void TimeClusterFinder::fillTimeSpectrum() {
    _timespec.Reset();
    // loop over straws hits and fill time spectrum plot for selected hits
    unsigned nstrs = _shcol->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      if(goodHit(_shfcol->at(istr))) {
	double time = _ttcalc.strawHitTime(_shcol->at(istr),_shpcol->at(istr));
	_timespec.Fill(time);
      }
    }
    // if we're using the calorimeter clusters, put those in too
    if(_usecc && _cccol != 0){
      for(auto icc = _cccol->begin();icc != _cccol->end(); ++icc){
	if(goodCaloCluster(*icc)){
	  double time = _ttcalc.caloClusterTime(*icc);
	  // weight the cluster WRT hits by an ad-hoc value.  This is more about signal/noise than resolution
	  _timespec.Fill(time, _ccwt);
	}
      }
    }
  }

  void TimeClusterFinder::addCaloCluster(TimeCluster& tc,art::Event const& evt) {
    // take the highest-energy cluster that's consistent with this time
    if(_cccol != 0){
      auto bestcc = _cccol->end();
      for(auto icc = _cccol->begin();icc != _cccol->end(); ++icc){
	if(goodCaloCluster(*icc)){
	  double time = _ttcalc.caloClusterTime(*icc);
	  if(fabs(tc._t0._t0-time) < _maxdt) {
	    if(bestcc == _cccol->end() || icc->energyDep() > bestcc->energyDep())
	      bestcc = icc;
	  }
	}
      }
      if(bestcc != _cccol->end()){
	size_t index = std::distance(_cccol->begin(),bestcc);
	auto ccH = evt.getValidHandle<CaloClusterCollection>(_ccTag);
	tc._caloCluster = art::Ptr<CaloCluster>(ccH,index);
      }
    }
  }

}  // end namespace mu2e

using mu2e::TimeClusterFinder;
DEFINE_ART_MODULE(TimeClusterFinder);

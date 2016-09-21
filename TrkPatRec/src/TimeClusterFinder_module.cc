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
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/TimeClusterCollection.hh"
// root
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGMsgBox.h"
#include "TTree.h"
#include "TMarker.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/error_of_mean.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
// C++
#include <memory>
#include <algorithm>
using namespace std; 
using namespace boost::accumulators;


namespace mu2e {
// struct to run MVA
   struct TimePeakMVA {
    std::vector<Double_t> _pars;
    Double_t& _dt;
    Double_t& _dphi;
    Double_t& _rho;
    TimePeakMVA() : _pars(3,0.0), _dt(_pars[0]), _dphi(_pars[1]), _rho(_pars[2]) {}
  };
// comparison functor for sorting peaks
  struct PeakSort : public binary_function<TimeCluster, TimeCluster, bool> {
    bool operator()(TimeCluster const& x, TimeCluster const& y) {
      return x._strawHitIdxs.size() > y._strawHitIdxs.size();
    }
  };
// comparison functor for sorting hitIndices by time
  struct hitTimeSort : public binary_function<hitIndex, hitIndex, bool> {
    hitTimeSort(const StrawHitCollection* shcol) : _shcol(shcol){}
    bool operator() (hitIndex const& x, hitIndex const& y){
      return _shcol->at(x).time() > _shcol->at(y).time();
    }
    const StrawHitCollection* _shcol;
  };

  class TimeClusterFinder : public art::EDProducer {
  public:

    explicit TimeClusterFinder(fhicl::ParameterSet const& pset);
    virtual ~TimeClusterFinder();

    virtual void beginJob();

    // This is called for each event.
    void produce(art::Event & e);

  private:

    // Start: run time parameters

    int           _debug;
    int           _printfreq;

    unsigned      _iev;

    // event object labels
    art::InputTag			_shTag;
    art::InputTag			_shpTag;
    art::InputTag			_shfTag;
    
    StrawHitFlag  _hsel, _hbkg;
    double        _maxdt;

    // time spectrum parameters
    int           _t0TypeCalculator;
    unsigned      _maxnpeak;
    unsigned      _minnhits;
    double        _minpeakmva, _maxpeakdt, _maxpeakdphi;
    double        _tmin;
    double        _tmax;
    double        _tbin;
    unsigned      _nbins;
    double        _ymin;
    double        _1dthresh,_tssigma;
    double        _nmnlWdNSigma;

    // cache of event objects
    const StrawHitCollection*             _shcol;
    const StrawHitFlagCollection*         _shfcol;
    const StrawHitPositionCollection*     _shpcol;

    MVATools                              _peakMVA; // MVA for peak cleaning
    TimePeakMVA                           _pmva; // input variables to TMVA for peak cleaning

    art::Handle<mu2e::StrawHitCollection> _strawhitsH;

    void initCluster	  (TimeCluster& tp); // fill peak information from the list of hits
    void refineCluster	  (TimeCluster& tp); // refine the peak information and hit list
    bool findData         (const art::Event& e);
    void findPeaks    (TimeClusterCollection*);
    void findPeaks    (TH1F const& tspect, std::vector<Float_t>& xpeak, std::vector<Float_t>& ypeak);
    bool goodHit(StrawHitFlag const& flag) const;

  };

  TimeClusterFinder::~TimeClusterFinder() {
  }
  
  TimeClusterFinder::TimeClusterFinder(fhicl::ParameterSet const& pset) :
    _debug             (pset.get<int>("debugLevel",0)),
    _printfreq         (pset.get<int>("printFrequency",101)),
    _shTag	 (pset.get<art::InputTag>("StrawHitCollection","makeSH")),
    _shpTag	 (pset.get<art::InputTag>("StrawHitPositionCollection","MakeStereoHits")),
    _shfTag	 (pset.get<art::InputTag>("StrawHitFlagCollection","FlagBkgHits")),
    _hsel              (pset.get<std::vector<std::string> >("HitSelectionBits")),
    _hbkg              (pset.get<vector<string> >("HitBackgroundBits",vector<string>{"DeltaRay","Isolated"})),
    _maxdt             (pset.get<double>("DtMax",30.0)),
    _t0TypeCalculator  (pset.get<int>("T0TypeCalculator",3)),
    _maxnpeak          (pset.get<unsigned>("MaxNPeaks",50)),
    _minnhits          (pset.get<unsigned>("MinNHits",10)),
    _minpeakmva        (pset.get<double>("MinTimePeakMVA",0.2)),
    _maxpeakdt         (pset.get<double>("MaxTimePeakDeltat",25.0)),
    _maxpeakdphi       (pset.get<double>("MaxTimePeakDeltaPhi",1.0)),
    _tmin              (pset.get<double>("tmin",500.0)),
    _tmax              (pset.get<double>("tmax",1700.0)),
    _tbin              (pset.get<double>("tbin",20.0)),
    _ymin              (pset.get<double>("ymin",8.0)),
    _1dthresh          (pset.get<double>("OneDPeakThreshold",5.0)),
    _nmnlWdNSigma      (pset.get<double>("NominalWidthNSigma",3.0)),
    _peakMVA           (pset.get<fhicl::ParameterSet>("PeakCleanMVA",fhicl::ParameterSet()))
  {
    // set # bins for time spectrum plot
    _nbins = (unsigned)rint((_tmax-_tmin)/_tbin);
    // Tell the framework what we make.
    produces<TimeClusterCollection>();
    produces<StrawHitFlagCollection>();

  }

  void TimeClusterFinder::beginJob(){
    _peakMVA.initMVA();
    if(_debug > 0){
      cout << "TimeClusterFinder MVA : " << endl; 
      _peakMVA.showMVA();
    }
    _iev = 0;
  }

  void TimeClusterFinder::produce(art::Event & event ) {

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
 
    // find the time peaks in the time spectrum of selected hits.  
    // This also estimates the t0 value and error, and associates hits
    findPeaks(tccol.get());

    // set the flag for all hits associated to a time peak
    if(_debug > 0) cout << "Found " << tccol->size() << " Time Clusters " << endl;
    for (auto tpc : *tccol.get()){
      if(_debug > 1){
	cout << "Time Cluster time = " << tpc.t0().t0() << " +- " << tpc.t0().t0Err()
	<< " position = " << tpc._pos << endl;
      }
      for (auto shi : tpc._strawHitIdxs ) {
	flags->at(shi).merge(StrawHitFlag::tclust);
	if(_debug > 2){
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

    auto shH = evt.getValidHandle<StrawHitCollection>(_shTag);
    _shcol = shH.product();
    auto shpH = evt.getValidHandle<StrawHitPositionCollection>(_shpTag);
    _shpcol = shpH.product();
    auto shfH = evt.getValidHandle<StrawHitFlagCollection>(_shfTag);
    _shfcol = shfH.product();

     return _shcol != 0 && _shfcol != 0 && _shpcol != 0;
   }

   void TimeClusterFinder::findPeaks(TimeClusterCollection* tclusts) {
     TH1F timespec("timespec","time spectrum",_nbins,_tmin,_tmax);
     // loop over straws hits and fill time spectrum plot for selected hits
     unsigned nstrs = _shcol->size();
     for(unsigned istr=0; istr<nstrs;++istr){
       if(goodHit(_shfcol->at(istr))) {
         double time = _shcol->at(istr).time();
         timespec.Fill(time);
       }
     }
     // find local maxima in the spectrum
     std::vector<Float_t> xpeaks,ypeaks;
     findPeaks(timespec,xpeaks,ypeaks);
     unsigned np = xpeaks.size();
     // Loop over peaks, looking only at those with a minimum peak value
     for (unsigned ip=0; ip<np; ++ip) {
       Float_t xp = xpeaks[ip];
       Float_t yp = ypeaks[ip];
       if(yp > _ymin){
       // time cluster
	 TimeCluster tclust;
	// initial t0 is the peak position
	 tclust._t0 = TrkT0(xp,1.0);
     // associate all hits in the time window with this peak
         for(size_t istr=0; istr<nstrs;++istr){
           if(goodHit(_shfcol->at(istr))){
             double time = _shcol->at(istr).time();
             if(fabs(time-xp) < _maxdt){
               tclust._strawHitIdxs.push_back(hitIndex(istr));
             }
           }
         }
	 // initialize the peak information
	 initCluster(tclust);
	 // refine the peak information
         refineCluster(tclust);
	 // final check on peak size
         if(tclust._strawHitIdxs.size() >= _minnhits)tclusts->push_back(tclust);
       }
     }
     // sort the peaks so that the largest comes first.  Not sure if this is really necessary
     static PeakSort psort;
     sort(tclusts->begin(),tclusts->end(),psort);
   }

   void TimeClusterFinder::initCluster(TimeCluster& tclust) {
     // use medians to initialize robustly
     accumulator_set<double, stats<tag::min > > tmin;
     accumulator_set<double, stats<tag::max > > tmax;
     accumulator_set<double, stats<tag::median(with_p_square_quantile) > > tacc;
     accumulator_set<double, stats<tag::median(with_p_square_quantile) > > xacc;
     accumulator_set<double, stats<tag::median(with_p_square_quantile) > > yacc;
     accumulator_set<double, stats<tag::median(with_p_square_quantile) > > zacc;
     unsigned nstrs = tclust._strawHitIdxs.size();
     for(unsigned istr=0; istr<nstrs;++istr){
       unsigned ish = tclust._strawHitIdxs[istr];
       if(goodHit(_shfcol->at(ish))){
	 tmin(_shcol->at(ish).time());
	 tmax(_shcol->at(ish).time());
	 tacc(_shcol->at(ish).time());
	 xacc(_shpcol->at(ish).pos().x());
	 yacc(_shpcol->at(ish).pos().y());
	 zacc(_shpcol->at(ish).pos().z());
       }
     }
     // set peak info
     static double invsqrt12(1.0/sqrt(12.0));
     tclust._t0._t0 = median(tacc);
     tclust._t0._t0err = ( boost::accumulators::extract::max(tmax)-boost::accumulators::extract::min(tmin))*invsqrt12/sqrt(nstrs);
     tclust._pos = CLHEP::Hep3Vector(median(xacc),median(yacc),median(zacc));
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
         double dt = _shcol->at(ish).time() - ptime;
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
       accumulator_set<double, stats<tag::mean > > tacc;
       for(size_t ips=0;ips<tclust._strawHitIdxs.size();++ips){
         unsigned ish = tclust._strawHitIdxs[ips];
         double time = _shcol->at(ish).time();
         double phi = _shpcol->at(ish).pos().phi();
         Angles::deltaPhi(phi,pphi);
         tacc(time);
         facc(phi);
       }
       pphi = extract_result<tag::mean>(facc);
       ptime = extract_result<tag::mean>(tacc);
     } while(tclust._strawHitIdxs.size() >= _minnhits && worstmva < _minpeakmva);
 // final pass: hard cut on dt and dphi
     vector<size_t> toremove;
     accumulator_set<double, stats<tag::mean > > facc;
     accumulator_set<double, stats<tag::error_of<tag::mean> > > tacc;
     accumulator_set<double, stats<tag::mean > > racc;
     accumulator_set<double, stats<tag::mean > > zacc;
     for(size_t ips=0;ips<tclust._strawHitIdxs.size();++ips){
       unsigned ish = tclust._strawHitIdxs[ips];
       double dt = _shcol->at(ish).time() - ptime;
       double phi = _shpcol->at(ish).pos().phi();
       double rho = _shpcol->at(ish).pos().perp();
       double dphi = Angles::deltaPhi(phi,pphi);
       if(fabs(dt) < _maxpeakdt && fabs(dphi) < _maxpeakdphi){
         tacc(_shcol->at(ish).time());
         facc(phi);
	 racc(rho);
	 zacc(_shpcol->at(ish).pos().z());
       } else {
	 toremove.push_back(ips);
       }
     }
     // actually remove the hits; must start from the back
     std::sort(toremove.begin(),toremove.end(),std::greater<size_t>());
     for(auto irm=toremove.begin();irm!=toremove.end();++irm){
       std::swap(tclust._strawHitIdxs[*irm],tclust._strawHitIdxs.back());
       tclust._strawHitIdxs.pop_back();
     }
     // update peak properties
     tclust._t0._t0 = extract_result<tag::mean>(tacc);
     tclust._t0._t0err = error_of<tag::mean>(tacc);
     pphi = extract_result<tag::mean>(facc);
     double prho = extract_result<tag::mean>(racc);
     double zpos = extract_result<tag::mean>(zacc);
// update position
     tclust._pos = CLHEP::Hep3Vector(prho*cos(pphi),prho*sin(pphi),zpos);
   }

   void TimeClusterFinder::findPeaks( TH1F const& tspect,std::vector<Float_t>& xpeak,std::vector<Float_t>& ypeak) {
     int ibin(1); // root starts bin numbers at 1
     int nbins = tspect.GetNbinsX();
     while(ibin < nbins+1){
       double y = tspect.GetBinContent(ibin);
       if(y >= _1dthresh){
         double xmax = tspect.GetBinCenter(ibin);
         double ymax = y;
         double yprev = y;
         // find contiguous bins above threshold
         int jbin = ibin+1;
         bool descending(false);
         while(jbin < nbins+1 && tspect.GetBinContent(jbin) >= _1dthresh){
           y =tspect.GetBinContent(jbin);
           descending |= yprev-y > sqrt(yprev);
           // don't follow next maximum
           if(descending && y-yprev > sqrt(y)){
             break;
           } else {
             if(y > ymax){
               ymax = y;
               xmax = tspect.GetBinCenter(jbin);
             }
             yprev = y;
             ibin = jbin;
             ++jbin;
           }
         }
         xpeak.push_back(xmax);
         ypeak.push_back(ymax);
       }
       ++ibin;
     }
   }

  bool TimeClusterFinder::goodHit(StrawHitFlag const& flag) const {
    return flag.hasAllProperties(_hsel) && !flag.hasAnyProperty(_hbkg);
  }

}  // end namespace mu2e

using mu2e::TimeClusterFinder;
DEFINE_ART_MODULE(TimeClusterFinder);

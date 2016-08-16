//
// TTracker time peak finder
//
// $Id: TimePeakFinder_module.cc,v 1.3 2014/08/25 12:08:29 tassiell Exp $
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
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
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
      return _shcol->at(x._index).time() > _shcol->at(y._index).time();
    }
    const StrawHitCollection* _shcol;
  };

  class TimePeakFinder : public art::EDProducer {
  public:

    explicit TimePeakFinder(fhicl::ParameterSet const& pset);
    virtual ~TimePeakFinder();

    virtual void beginJob();

    // This is called for each event.
    void produce(art::Event & e);

  private:

    // Start: run time parameters

    int           _debug;
    int           _printfreq;

    unsigned      _iev;

    // event object labels
    std::string   _shLabel;
    std::string   _shpLabel;
    std::string   _shfLabel;
    std::string   _mcdigislabel;
    
    StrawHitFlag  _hsel, _psel, _hbkg;
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

    void initPeak	  (TimeCluster& tp); // fill peak information from the list of hits
    void refinePeak	  (TimeCluster& tp); // refine the peak information and hit list
    bool findData         (const art::Event& e);
    void findPeaks    (TimeClusterCollection*);
    void findPeaks    (TH1F const& tspect, std::vector<Float_t>& xpeak, std::vector<Float_t>& ypeak);
    bool goodHit(StrawHitFlag const& flag) const;

  };

  TimePeakFinder::~TimePeakFinder() {
  }
  
  TimePeakFinder::TimePeakFinder(fhicl::ParameterSet const& pset) :
    _debug             (pset.get<int>("debugLevel",0)),
    _printfreq         (pset.get<int>("printFrequency",101)),
    _shLabel           (pset.get<std::string>("StrawHitCollectionLabel","makeSH")),
    _shpLabel          (pset.get<std::string>("StrawHitPositionCollectionLabel","MakeStereoHits")),
    _shfLabel          (pset.get<std::string>("StrawHitFlagCollectionLabel","FlagBkgHits")),
    _mcdigislabel      (pset.get<string>("StrawDigiMCLabel")),
    _hsel              (pset.get<std::vector<std::string> >("HitSelectionBits")),
    _psel              (pset.get<std::vector<std::string> >("PositionSelectionBits")),
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

  }

  void TimePeakFinder::beginJob(){
    _peakMVA.initMVA();
    if(_debug > 0){
      cout << "TimePeakFinder MVA : " << endl; 
      _peakMVA.showMVA();
    }
    // create diagnostics if requested 
    if(_diag > 0)createDiagnostics();
    _iev = 0;
  }

  void TimePeakFinder::produce(art::Event & event ) {

    _iev=event.id().event();
    if(_debug > 0 && (_iev%_printfreq)==0)cout<<"TimePeakFinder: event="<<_iev<<endl;
    // find the data
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::TimePeakFinder: data missing or incomplete"<< endl;
      return;
    }
    // create time peak collection
    std::unique_ptr<TimeClusterCollection>    thcc  (new TimeClusterCollection);
    // copy in the existing flags
    _flags = new StrawHitFlagCollection(*_shfcol);
    std::unique_ptr<StrawHitFlagCollection> flags(_flags );
 
    // find the time peaks in the time spectrum of selected hits.  
    // This also estimates the t0 value and error, and associates hits
    findPeaks(thcc.get());

    // set the flag for all hits associated to a time peak
    for (auto tpc : thcc.get()){
      for (auto shi : tpc._strawHitIdxs ) {
	flags->at(shi._index).merge(StrawHitFlag::trksel);
      }
    }
    // put collections into the event
    event.put(std::move( thcc ));
    event.put(std::move( flags));

  }

   // find the input data objects
  bool TimePeakFinder::findData(const art::Event& evt){
    _shcol  = 0;
    _shfcol = 0;
    _shpcol = 0;

    if(evt.getByLabel(_shLabel,_strawhitsH))
      _shcol = _strawhitsH.product();
    art::Handle<mu2e::StrawHitPositionCollection> shposH;
    if(evt.getByLabel(_shpLabel,shposH))
       _shpcol = shposH.product();
     art::Handle<mu2e::StrawHitFlagCollection> shflagH;
     if(evt.getByLabel(_shfLabel,shflagH))
       _shfcol = shflagH.product();

     return _shcol != 0 && _shfcol != 0 && _shpcol != 0;
   }

   void TimePeakFinder::findPeaks(TimeClusterCollection* tpeaks) {
     TH1F timespec("timespec","time spectrum",_nbins,_tmin,_tmax);
     // loop over straws hits and fill time spectrum plot for selected hits
     unsigned nstrs = _shcol->size();
     for(unsigned istr=0; istr<nstrs;++istr){
       if(goodHit(_shfcol->at(istr)) {
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
	 TimeCluster tpeak;
	// initial t0 is the peak position
	 tpeak._t0 = TrkT0(xp,1.0);
     // associate all hits in the time window with this peak
         for(size_t istr=0; istr<nstrs;++istr){
           if(goodHit(_shfcol->at(istr))){
             double time = _shcol->at(istr).time();
             if(fabs(time-xp) < _maxdt){
               tpeak._strawHitIdxs.push_back(hitIndex(istr));
             }
           }
         }
	 // initialize the peak information
	 initPeak(tpeak);
	 // refine the peak information
         refinePeak(tpeak);
	 // final check on peak size
         if(tpeak._strawHitIdxs.size() >= _minnhits)tpeaks.push_back(tpeak);
       }
     }
     // sort the peaks so that the largest comes first.  Not sure if this is really necessary
     static PeakSort psort;
     sort(tpeaks.begin(),tpeaks.end(),psort);
   }

   void TimePeakFinder::initPeak(TimeCluster const& tpeak) {
     // use medians to initialize robustly
     accumulator_set<double, stats<tag::median(with_p_square_quantile) >, stats<tag::max>, stats<tag::min> > tacc;
     accumulator_set<double, stats<tag::median(with_p_square_quantile) > > pacc;
     accumulator_set<double, stats<tag::median(with_p_square_quantile) > > racc;
     accumulator_set<double, stats<tag::median(with_p_square_quantile) > > zacc;
     unsigned nstrs = tpeak._strawHitIdxs.size();
     for(unsigned istr=0; istr<nstrs;++istr){
       unsigned ish = tpeak._strawHitIdxs[istr]._index;
       if(goodHit(_shfcol->at(ish))){
         double time = _shcol->at(ish).time();
         double rho = _shpcol->at(ish).pos().perp();
         double phi = _shpcol->at(ish).pos().phi();
         double dphi = Angles::deltaPhi(phi,pphi);
         tacc(time);
         pacc(phi);
         racc(rho);
	 zacc(_shpcol->at(ish).pos().z());
       }
     }
     // set peak info
     static double invsqrt12(1.0/sqrt(12.0));
     tpeak._t0._t0 = median(tacc);
     tpeak._t0._t0err = (max(tacc)-min(tacc))*invsqrt12/sqrt(count(tacc));
     double rho = median(racc);
     double phi = median(pacc);
    // use the rho and phi to define a position. 
     tpeak._pos = Hep3Vector(rho*cos(phi),rho*sin(phi),median(zacc));
   }

   void TimePeakFinder::refinePeak(TimeCluster& tpeak) {
     // iteratively filter the worst outlier
     double worstmva(100.0);
     double pphi(tpeak._pos.phi());
     double ptime(tpeak._t0);
     do{
 // Use the MVA to find the least signal-like hit
       size_t iworst =0;
       worstmva = 100.0;
       for(size_t ips=0;ips<tpeak._strawHitIdxs.size();++ips){
         unsigned ish = tpeak._strawHitIdxs[ips]._index;
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
	 // accumulate new averages
	 tacc(time);
         facc(phi);
       }
       // remove the worst hit
       if(worstmva < _minpeakmva){
         std::swap(tpeak._strawHitIdxs[iworst],tpeak._strawHitIdxs.back());
         tpeak._strawHitIdxs.pop_back();
       }
   // re-compute the average phi and range
       accumulator_set<double, stats<tag::mean > > facc;
       accumulator_set<double, stats<tag::mean > > tacc;
       for(size_t ips=0;ips<tpeak._strawHitIdxs.size();++ips){
         unsigned ish = tpeak._strawHitIdxs[ips]._index;
         double dt = _shcol->at(ish).time() - ptime;
         double phi = _shpcol->at(ish).pos().phi();
         double dphi = Angles::deltaPhi(phi,pphi);
         tacc(time);
         facc(phi);
       }
       pphi = extract_result<tag::mean>(facc);
       ptime = extract_result<tag::mean>(tacc);
     } while(tpeak._strawHitIdxs.size() >= _minnhits && worstmva < _minpeakmva);
 // final pass: hard cut on dt and dphi
     vector<size_t> toremove;
     accumulator_set<double, stats<tag::mean > > facc;
     accumulator_set<double, stats<tag::mean >, stats<tag::error_of<tag::mean> > > tacc;
     accumulator_set<double, stats<tag::mean > > racc;
     accumulator_set<double, stats<tag::mean > > zacc;
     for(size_t ips=0;ips<tpeak._strawHitIdxs.size();++ips){
       unsigned ish = tpeak._strawHitIdxs[ips]._index;
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
       std::swap(tpeak._strawHitIdxs[*irm],tpeak._strawHitIdxs.back());
       tpeak._strawHitIdxs.pop_back();
     }
     // update peak properties
     tpeak._t0._t0 = extract_result<tag::mean>(tacc);
     tpeak._t0._t0err = error_of<tag::mean>(tacc);
     pphi = extract_result<tag::mean>(facc);
     double prho = extract_result<tag::mean>(racc);
     double zpos = extract_result<tag::mean>(zacc);
// update position
     tpeak._pos = CLHEP::Hep3Vector(prho*cos(pphi),prho*sin(pphi),zpos);
   }

   void TimePeakFinder::findPeaks( TH1F const& tspect,std::vector<Float_t>& xpeak,std::vector<Float_t>& ypeak) {
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

  bool TimePeakFinder::goodHit(StrawHitFlag const& flag) const {
    return flag.hasAllProperties(_hsel) && flag.hasAnyProperty(_psel) && !flag.hasAnyProperty(_hbkg);
  }

}  // end namespace mu2e

using mu2e::TimePeakFinder;
DEFINE_ART_MODULE(TimePeakFinder);

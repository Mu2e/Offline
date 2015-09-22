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
// conditions
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
// Mu2e
#include "TrkPatRec/inc/TrkPatRec.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/TrackerHitTimeClusterCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
// root
#include "TFile.h"
#include "TH1F.h"
#include "TMVA/Reader.h"
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
// C++
#include <memory>
using namespace std; 
using namespace boost::accumulators;


namespace mu2e {

  class TimePeakFinder : public art::EDProducer {
  public:

    explicit TimePeakFinder(fhicl::ParameterSet const& pset);
    virtual ~TimePeakFinder();

    virtual void beginJob();
    void endJob();

    // This is called for each event.
    void produce(art::Event & e);

  private:

    // Start: run time parameters

    int _debug;
    int _diag;
    int _printfreq;

    unsigned _iev;

    // event object labels
    std::string _shLabel;
    std::string _shpLabel;
    std::string _shfLabel;
    std::string _mcdigislabel;
    StrawHitFlag _tsel;
    StrawHitFlag _tbkg;
    double _maxdt;
    // time spectrum parameters
    bool _findtpeak;
    unsigned _maxnpeak;
    unsigned _minnhits;
    bool _cleanpeaks;
    double _minpeakmva, _maxpeakdt, _maxpeakdphi;
    std::string _PMVAType; // type of MVA
    std::string _PMVAWeights; // file of MVA weights
    double _tmin;
    double _tmax;
    double _tbin;
    unsigned _nbins;
    double _ymin;
    double _1dthresh,_tssigma;
    double _nmnlWdNSigma;

    // cache of event objects
    const StrawHitCollection* _shcol;
    const StrawHitFlagCollection* _shfcol;
    StrawHitFlagCollection* _flags;
    const StrawHitPositionCollection* _shpcol;

    // cache of time peaks
    std::vector<TrkTimePeak> _tpeaks;
    //string _iname; // data instance name

    void initializeReaders();
    TMVA::Reader *_peakMVA; // MVA for peak cleaning
    TimePeakMVA _pmva; // input variables to TMVA for peak cleaning

    art::Handle<mu2e::StrawHitCollection> _strawhitsH;

    bool findData(const art::Event& e);
    void findTimePeaks();
    void findTimePeaks(TH1F const& tspect,std::vector<Float_t>& xpeak,std::vector<Float_t>& ypeak);
    void createTimePeak();
    void cleanTimePeak(TrkTimePeak& tp);
    void createDiagnostics();
    void fillTimeDiag();
    void fillPeakDiag(size_t ip, TrkTimePeak const& tp);
    const StrawDigiMCCollection *_mcdigis;
// time peak diag variables
    TTree* _tpdiag;
    Int_t _tpeventid, _peakid, _pmax, _nphits, _ncphits, _nchits;
    Float_t _ptime, _pdtimemax, _ctime, _cdtimemax;;
    Float_t _pphi, _cphi, _cphirange, _pdphimax, _cdphimax;
    vector<TimePeakHitInfo> _tphinfo;
  };

  TimePeakFinder::~TimePeakFinder() {
    if (_peakMVA!=NULL) { delete _peakMVA; }
  }
  
  TimePeakFinder::TimePeakFinder(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>("debugLevel",0)),
    _diag(pset.get<int>("diagLevel",0)),
    _printfreq(pset.get<int>("printFrequency",101)),
    _shLabel(pset.get<std::string>("StrawHitCollectionLabel","makeSH")),
    _shpLabel(pset.get<std::string>("StrawHitPositionCollectionLabel","MakeStereoHits")),
    _shfLabel(pset.get<std::string>("StrawHitFlagCollectionLabel","FlagBkgHits")),
    _mcdigislabel(pset.get<string>("StrawDigiMCLabel")),
    _tsel(pset.get<std::vector<std::string> >("TimeSelectionBits")),
    _tbkg(pset.get<vector<string> >("TimeBackgroundBits",vector<string>{"DeltaRay","Isolated"})),
    _maxdt(pset.get<double>("DtMax",30.0)),
    _findtpeak(pset.get<bool>("FindTimePeaks",true)),
    _maxnpeak(pset.get<unsigned>("MaxNPeaks",50)),
    _minnhits(pset.get<unsigned>("MinNHits",10)),
    _cleanpeaks(pset.get<bool>("CleanTimePeaks",true)),
    _minpeakmva(pset.get<double>("MinTimePeakMVA",0.5)),
    _maxpeakdt(pset.get<double>("MaxTimePeakDeltat",25.0)),
    _maxpeakdphi(pset.get<double>("MaxTimePeakDeltaPhi",1.0)),
    _PMVAType(pset.get<std::string>("TimePeakMVAType","MLP method")),
    _tmin(pset.get<double>("tmin",500.0)),
    _tmax(pset.get<double>("tmax",1700.0)),
    _tbin(pset.get<double>("tbin",20.0)),
    _ymin(pset.get<double>("ymin",8.0)),
    _1dthresh(pset.get<double>("OneDPeakThreshold",5.0)),
    _nmnlWdNSigma(pset.get<double>("NominalWidthNSigma",3.0)),
    _peakMVA(NULL)
  {
    // tag the data product instance by the direction and particle type found by this fitter
    produces<StrawHitFlagCollection>();
    // set # bins for time spectrum plot
    _nbins = (unsigned)rint((_tmax-_tmin)/_tbin);
    // location-independent files
    ConfigFileLookupPolicy configFile;
    std::string weights = pset.get<std::string>("PeakMVAWeights","TrkPatRec/test/TimePeak.weights.xml");
    _PMVAWeights = configFile(weights);

    // Tell the framework what we make.
    produces<TrackerHitTimeClusterCollection>();

  }

  void TimePeakFinder::beginJob(){

    initializeReaders();
    // create diagnostics if requested 
    if(_diag > 0)createDiagnostics();
    _iev = 0;
  }

  void TimePeakFinder::produce(art::Event & event ) {

    _iev=event.id().event();
    if((_debug > 0 && _iev%_printfreq)==0)cout<<"TimePeakFinder: event="<<_iev<<endl;
    // find the data
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::TimePeakFinder: data missing or incomplete"<< endl;
      return;
    }
    // copy in the existing flags
    _flags = new StrawHitFlagCollection(*_shfcol);
    std::unique_ptr<StrawHitFlagCollection> flags(_flags );

    // find the time peaks in the time spectrum of selected hits.  Otherwise, take all
    // selected hits as a peak
    _tpeaks.clear();
    if(_findtpeak){
      findTimePeaks();
    } else {
      createTimePeak();
    }
    if(_diag > 2 && _mcdigis !=0 && _mcdigis->size() > 0)fillTimeDiag();

    // create event objects and piut them into event
    std::unique_ptr<TrackerHitTimeClusterCollection> thcc(new TrackerHitTimeClusterCollection);
    for (std::vector<TrkTimePeak>::iterator itp=_tpeaks.begin(); itp!=_tpeaks.end(); ++itp) {
      thcc->push_back(TrackerHitTimeCluster());
      thcc->back()._meanTime=itp->_tpeak;
      thcc->back()._peakmax=itp->_peakmax;

      thcc->back()._minHitTime = thcc->back()._maxHitTime = _shcol->at( itp->_trkptrs.begin()->_index ).time();
      double sumOfSqr(0.0), sum(0.0);
      for (std::vector<hitIndex>::iterator hittpit=itp->_trkptrs.begin(); hittpit!=itp->_trkptrs.end(); ++hittpit) {
	double htime = _shcol->at(hittpit->_index).time();
	sum+=htime;
	sumOfSqr+=htime*htime;
	thcc->back()._selectedTrackerHits.push_back( art::Ptr<StrawHit> ( _strawhitsH, hittpit->_index ) );
	if ( htime < thcc->back()._minHitTime ) { thcc->back()._minHitTime=htime; }
	if ( htime > thcc->back()._maxHitTime ) { thcc->back()._maxHitTime=htime; }
	flags->at(hittpit->_index).merge(StrawHitFlag::timesel);
      }

      if (itp->_trkptrs.size()>1) {
	thcc->back()._sigma=sumOfSqr-sum*sum/((double)itp->_trkptrs.size());
	if (thcc->back()._sigma>0.0) {
	  thcc->back()._sigma=sqrt(thcc->back()._sigma/((double)(itp->_trkptrs.size()-1.0)));
	  thcc->back()._nominalWidth=_nmnlWdNSigma*thcc->back()._sigma;
	} else {
	  thcc->back()._sigma=0.0;
	  thcc->back()._nominalWidth = thcc->back()._maxHitTime-thcc->back()._minHitTime;
	}
      } else {
	thcc->back()._nominalWidth = thcc->back()._maxHitTime-thcc->back()._minHitTime;
      }

    }

    int iPeak(0);
    if (_debug>0) {
      std::cout<<"n Peaks "<<thcc->size()<<std::endl;
      for (TrackerHitTimeClusterCollection::iterator tpeak_it=thcc->begin(); tpeak_it!=thcc->end(); ++tpeak_it ) {
	std::cout<<"iPeak "<<std::endl;
	std::cout<<*tpeak_it<<std::endl;
	++iPeak;
      }
    }

    event.put(std::move(thcc));
    event.put(std::move(flags));

  } // end produce

  void TimePeakFinder::endJob(){
  }

  // find the input data objects
   bool TimePeakFinder::findData(const art::Event& evt){
     _shcol = 0;
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

    if(_diag > 0){
      art::Handle<StrawDigiMCCollection> mcdigisHandle;
      if(evt.getByLabel(_mcdigislabel,"StrawHitMC",mcdigisHandle))
	_mcdigis = mcdigisHandle.product();
    }

     return _shcol != 0 && _shfcol != 0 && _shpcol != 0;
   }

   void TimePeakFinder::findTimePeaks() {
     TH1F timespec("timespec","time spectrum",_nbins,_tmin,_tmax);
     // loop over straws hits and fill time spectrum plot for tight hits
     unsigned nstrs = _shcol->size();
     for(unsigned istr=0; istr<nstrs;++istr){
       if(_flags->at(istr).hasAllProperties(_tsel) && !_flags->at(istr).hasAnyProperty(_tbkg)){
         double time = _shcol->at(istr).time();
         timespec.Fill(time);
       }
     }
     std::vector<Float_t> xpeaks,ypeaks;
     findTimePeaks(timespec,xpeaks,ypeaks);
     unsigned np = xpeaks.size();
     // Loop over peaks, looking only at those with a minimum peak value
     for (unsigned ip=0; ip<np; ++ip) {
       Float_t xp = xpeaks[ip];
       Float_t yp = ypeaks[ip];
       TrkTimePeak tpeak(xp,yp);
       if(yp > _ymin){
     // associate all hits in the time window with this peak
         for(size_t istr=0; istr<nstrs;++istr){
           if(_flags->at(istr).hasAllProperties(_tsel) && !_flags->at(istr).hasAnyProperty(_tbkg)){
             double time = _shcol->at(istr).time();
             if(fabs(time-xp) < _maxdt){
               tpeak._trkptrs.push_back(istr);
             }
           }
         }
         // if requested, clean the peaks
         if(_cleanpeaks)cleanTimePeak(tpeak);
         if(tpeak._trkptrs.size() >= _minnhits)_tpeaks.push_back(tpeak);
       }
     }
     // sort the peaks so that the largest comes first
     sort(_tpeaks.begin(),_tpeaks.end(),greater<TrkTimePeak>());
     // if requested, fill diagnostics
    if(_diag>1 && _mcdigis != 0){
      for(size_t ip=0;ip<_tpeaks.size();++ip){
	TrkTimePeak const& tp = _tpeaks[ip];
	fillPeakDiag(ip,tp);
      }
    }
   }

   void TimePeakFinder::createTimePeak() {
     // find the median time
     accumulator_set<double, stats<tag::median(with_p_square_quantile) > > tacc;
     unsigned nstrs = _shcol->size();
     for(unsigned istr=0; istr<nstrs;++istr){
       if(_flags->at(istr).hasAllProperties(_tsel) && !_flags->at(istr).hasAnyProperty(_tbkg)){
         double time = _shcol->at(istr).time();
         tacc(time);
       }
     }
     unsigned np = boost::accumulators::extract::count(tacc);
     if(np >= _minnhits){
       double mtime  = median(tacc);
 // create a time peak from the full subset of selected hits
       TrkTimePeak tpeak(mtime,(double)nstrs);
       for(unsigned istr=0; istr<nstrs;++istr){
         if(_flags->at(istr).hasAllProperties(_tsel) && !_flags->at(istr).hasAnyProperty(_tbkg)){
           tpeak._trkptrs.push_back(istr);
         }
       }
       _tpeaks.push_back(tpeak);
     }
   }


   void TimePeakFinder::cleanTimePeak(TrkTimePeak& tpeak) {
     static const double twopi(2*M_PI);
     // iteratively filter outliers
     double worstmva(100.0);
     double pphi(1000.0);
     double ptime(-1000.0);
     do{
   // first, compute the average phi and range.  Take care for wrapping
       accumulator_set<double, stats<tag::mean > > facc;
       accumulator_set<double, stats<tag::mean > > tacc;
       for(size_t ips=0;ips<tpeak._trkptrs.size();++ips){
         unsigned ish = tpeak._trkptrs[ips]._index;
         double time = _shcol->at(ish).time();
         tacc(time);
         double phi = _shpcol->at(ish).pos().phi();
         if(extract_result<tag::count>(facc) > 0){
           double dphi = phi - extract_result<tag::mean>(facc);
           if(dphi > M_PI){
             phi -= twopi;
           } else if(dphi < -M_PI){
             phi += twopi;
           }
         }
         facc(phi);
       }
       pphi = extract_result<tag::mean>(facc);
       ptime = extract_result<tag::mean>(tacc);
 // find the least signal-like hit
       size_t iworst =0;
       worstmva = 100.0;
       for(size_t ips=0;ips<tpeak._trkptrs.size();++ips){
         unsigned ish = tpeak._trkptrs[ips]._index;
         double dt = _shcol->at(ish).time() - ptime;
         double rho = _shpcol->at(ish).pos().perp();
         double phi = _shpcol->at(ish).pos().phi();
         double dphi = phi - pphi;
         if(dphi > M_PI){
           dphi -= twopi;
         } else if(dphi < -M_PI){
           dphi += twopi;
         }
         // compute MVA
         _pmva._dt = dt;
         _pmva._dphi = dphi;
         _pmva._rho = rho;
         double mvaout = _peakMVA->EvaluateMVA(_PMVAType);
         if(mvaout < worstmva){
           worstmva = mvaout;
           iworst = ips;
         }
       }
       // remove the worst hit
       if(worstmva < _minpeakmva){
         std::swap(tpeak._trkptrs[iworst],tpeak._trkptrs.back());
         tpeak._trkptrs.pop_back();
       }
     } while(tpeak._trkptrs.size() >= _minnhits && worstmva < _minpeakmva);
 // final pass: hard cut on dt and dphi
     vector<size_t> toremove;
     accumulator_set<double, stats<tag::mean > > facc;
     accumulator_set<double, stats<tag::mean > > tacc;
     for(size_t ips=0;ips<tpeak._trkptrs.size();++ips){
       unsigned ish = tpeak._trkptrs[ips]._index;
       double dt = _shcol->at(ish).time() - ptime;
       double phi = _shpcol->at(ish).pos().phi();
       double dphi = phi - pphi;
       if(dphi > M_PI){
         dphi -= twopi;
         phi -= twopi;
       } else if(dphi < -M_PI){
         dphi += twopi;
         phi += twopi;
       }
       if(fabs(dt) < _maxpeakdt && fabs(dphi) < _maxpeakdphi){
         tacc(_shcol->at(ish).time());
         facc(phi);
       } else {
         toremove.push_back(ips);
       }
     }
     // actually remove the hits; must start from the back
     std::sort(toremove.begin(),toremove.end(),std::greater<size_t>());
     for(auto irm=toremove.begin();irm!=toremove.end();++irm){
       std::swap(tpeak._trkptrs[*irm],tpeak._trkptrs.back());
       tpeak._trkptrs.pop_back();
     }
     // update peak properties
     pphi = extract_result<tag::mean>(facc);
     ptime = extract_result<tag::mean>(tacc);
     tpeak._tpeak = ptime;
     tpeak._phi = pphi;
   }

   void TimePeakFinder::initializeReaders() {
     _peakMVA = new TMVA::Reader();
     _peakMVA->AddVariable("_dt",&_pmva._dt);
     _peakMVA->AddVariable("_dphi",&_pmva._dphi);
     _peakMVA->AddVariable("_rho",&_pmva._rho);
     _peakMVA->BookMVA(_PMVAType,_PMVAWeights);
   }

   void TimePeakFinder::findTimePeaks( TH1F const& tspect,std::vector<Float_t>& xpeak,std::vector<Float_t>& ypeak) {
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

  void TimePeakFinder::fillPeakDiag(size_t ip,TrkTimePeak const& tp) {
    _tpeventid = _iev;
    _peakid = ip;
    _pmax = tp._peakmax;
    _nphits = tp._trkptrs.size();
    _ncphits = 0;
    _pdtimemax = 0.0;
    _pdphimax = 0.0;
    _cdtimemax = 0.0;
    _cdphimax = 0.0;
    _tphinfo.clear();
    accumulator_set<double, stats<tag::mean > > facc;
    accumulator_set<double, stats<tag::mean > > tacc;
    for(size_t iph=0;iph<tp._trkptrs.size();++iph){
      unsigned ish = tp._trkptrs[iph]._index;
      StrawHit const& sh = _shcol->at(ish);
      StrawHitPosition const& shp = _shpcol->at(ish);
      double time = sh.time();
      tacc(time);
      CLHEP::Hep3Vector const& pos = shp.pos();
      double phi = pos.phi();
      if(extract_result<tag::count>(facc) > 0){
	double dphi = phi - extract_result<tag::mean>(facc);
	if(dphi > M_PI){
	  phi -= 2*M_PI;
	} else if(dphi < -M_PI){
	  phi += 2*M_PI;
	}
      }
      facc(phi);
    }
    _pphi =extract_result<tag::mean>(facc);
    _ptime =extract_result<tag::mean>(tacc);
    for(size_t iph=0;iph<tp._trkptrs.size();++iph){
      unsigned ish = tp._trkptrs[iph]._index;
      StrawHitPosition const& shp = _shpcol->at(ish);
      CLHEP::Hep3Vector const& pos = shp.pos();
      double phi = pos.phi();
      double dphi = phi - _pphi;
      if(dphi > M_PI)
	dphi -= 2*M_PI;
      else if(dphi < -M_PI)
	dphi += 2*M_PI;
      if(fabs(dphi) > _pdphimax)_pdphimax = fabs(dphi);

      double dt = _shcol->at(ish).time() - _ptime;
      if(fabs(dt) > _pdtimemax)_pdtimemax=fabs(dt);

      StrawDigiMC const& mcdigi = _mcdigis->at(ish);
      // use TDC channel 0 to define the MC match
      StrawDigi::TDCChannel itdc = StrawDigi::zero;
      if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
      art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
      art::Ptr<SimParticle> const& spp = spmcp->simParticle();
      int gid(-1);
      if(spp->genParticle().isNonnull())
	gid = spp->genParticle()->generatorId().id();

      bool conversion = (spp->pdgId() == 11 && gid == 2 && spmcp->momentum().mag()>90.0);

      if(conversion){
	_ncphits++;
	if(fabs(dt) > _cdtimemax)_cdtimemax = fabs(dt);
      }
      if(conversion && dphi > _cdphimax)_cdphimax = dphi;

      _pmva._dt = dt;
      _pmva._dphi = dphi;
      _pmva._rho = pos.perp();
      double mvaout = _peakMVA->EvaluateMVA(_PMVAType);

      TimePeakHitInfo tph;
      tph._dt = dt;
      tph._dphi = dphi;
      tph._rho = pos.perp();
      tph._mva = mvaout;
      tph._mcpdg = spp->pdgId();
      tph._mcgen = gid;
      tph._mcproc = spp->creationCode();
      _tphinfo.push_back(tph);
    }
    _tpdiag->Fill();
//    if(_ncphits > 0.5*_nchits)_icepeak = ip;
  }

  void TimePeakFinder::fillTimeDiag() {
    art::ServiceHandle<art::TFileService> tfs;
    TH1F *ctsp, *rtsp, *ttsp, *ltsp, *tdtsp, *ptsp;

    char rsname[100];
    char csname[100];
    char tsname[100];
    char lsname[100];
    char tdsname[100];
    char tpsname[100];
    snprintf(rsname,100,"rawtspectrum%i",_iev);
    snprintf(csname,100,"convtspectrum%i",_iev);
    snprintf(tsname,100,"tighttspectrum%i",_iev);
    snprintf(lsname,100,"loosetspectrum%i",_iev);
    snprintf(tdsname,100,"tightnodeltatspectrum%i",_iev);
    snprintf(tpsname,100,"protontspectrum%i",_iev);
    ttsp = tfs->make<TH1F>(tsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    ttsp->SetLineColor(kCyan);
    ltsp = tfs->make<TH1F>(lsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    ltsp->SetLineColor(kGreen);
    rtsp = tfs->make<TH1F>(rsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    rtsp->SetLineColor(kBlue);
    ctsp = tfs->make<TH1F>(csname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    ctsp->SetLineColor(kRed);
    ctsp->SetFillColor(kRed);
    ptsp = tfs->make<TH1F>(tpsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    ptsp->SetLineColor(kBlack);
    ptsp->SetFillColor(kBlack);
    tdtsp = tfs->make<TH1F>(tdsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    tdtsp->SetLineColor(kOrange);

    unsigned nstrs = _shcol->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      double time = _shcol->at(istr).time();
      bool conversion(false);
      bool proton(false);
      // summarize the MC truth for this strawhit
      if(_mcdigis != 0) {
	StrawDigiMC const& mcdigi = _mcdigis->at(istr);
	// use TDC channel 0 to define the MC match
	StrawDigi::TDCChannel itdc = StrawDigi::zero;
	if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
	art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
	art::Ptr<SimParticle> const& spp = spmcp->simParticle();
	int gid(-1);
	if(spp->genParticle().isNonnull())
	  gid = spp->genParticle()->generatorId().id();

	conversion = (spp->pdgId() == 11 && gid == 2 && spmcp->momentum().mag()>90.0);
	proton = spp->pdgId()==2212;
      }
      // fill plots
      rtsp->Fill(time);
      if(_flags->at(istr).hasAllProperties(_tsel)){
	ttsp->Fill(time);
      }
      if(_flags->at(istr).hasAllProperties(_tsel) && !_flags->at(istr).hasAnyProperty(_tbkg)){
	tdtsp->Fill(time);
      }
      if(_flags->at(istr).hasAllProperties(_tsel) && !_flags->at(istr).hasAnyProperty(_tbkg)){
	ltsp->Fill(time);
	if(proton)
	  ptsp->Fill(time);
      }
      if(conversion)
	ctsp->Fill(time);

    }
    // plot time peaks
    TList* flist = tdtsp->GetListOfFunctions();
    for(auto ipeak=_tpeaks.begin();ipeak!=_tpeaks.end();++ipeak){
      TMarker* smark = new TMarker(ipeak->_tpeak,ipeak->_peakmax,23);
      smark->SetMarkerColor(kRed);
      smark->SetMarkerSize(1.5);
      flist->Add(smark);
    }
  }

  void TimePeakFinder::createDiagnostics() {
    art::ServiceHandle<art::TFileService> tfs;
   // time peak diagnostics
    _tpdiag=tfs->make<TTree>("tpdiag","time peak diagnostics");
    _tpdiag->Branch("eventid",&_tpeventid,"eventid/I");
    _tpdiag->Branch("peakid",&_peakid,"peakid/I");
    _tpdiag->Branch("pmax",&_pmax,"pmax/I");
    _tpdiag->Branch("nphits",&_nphits,"nphits/I");
    _tpdiag->Branch("ncphits",&_ncphits,"ncphits/I");
    _tpdiag->Branch("nchits",&_nchits,"nchits/I");
    _tpdiag->Branch("ptime",&_ptime,"ptime/F");
    _tpdiag->Branch("ctime",&_ctime,"ctime/F");
    _tpdiag->Branch("pdtimemax",&_pdtimemax,"pdtimemax/F");
    _tpdiag->Branch("cdtimemax",&_cdtimemax,"cdtimemax/F");
    _tpdiag->Branch("pphi",&_pphi,"pphi/F");
    _tpdiag->Branch("cphi",&_cphi,"cphi/F");
    _tpdiag->Branch("cphirange",&_cphirange,"cphirange/F");
    _tpdiag->Branch("pdphimax",&_pdphimax,"pdphimax/F");
    _tpdiag->Branch("cdphimax",&_cdphimax,"cdphimax/F");
    _tpdiag->Branch("tphinfo",&_tphinfo);
 }



}  // end namespace mu2e

using mu2e::TimePeakFinder;
DEFINE_ART_MODULE(TimePeakFinder);

//
// Tracker time cluster finder
//
//
// Original author D. Brown and G. Tassielli
//
// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
// Mu2e
#include "GeneralUtilities/inc/Angles.hh"
#include "Mu2eUtilities/inc/MVATools.hh"
#include "Mu2eUtilities/inc/polyAtan2.hh"
// data
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
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

namespace {

  struct TimeCluMVA
  {
    vector<Float_t> _pars;
    Float_t& _dt;
    Float_t& _dphi;
    Float_t& _rho;
    Float_t& _nsh;
    Float_t& _plane;
    Float_t& _werr;
    Float_t& _wdist;
    
    TimeCluMVA() : _pars(7,0.0), _dt(_pars[0]), _dphi(_pars[1]), _rho(_pars[2]), _nsh(_pars[3]),
     _plane(_pars[4]), _werr(_pars[5]), _wdist(_pars[6]){}
//    TimeCluMVA() : _pars(5,0.0), _dt(_pars[0]), _nsh(_pars[1]),
//      _plane(_pars[2]), _werr(_pars[3]), _wdist(_pars[4]){}
  };
}

namespace mu2e {
   
  class TimeClusterFinder : public art::EDProducer
  {  
    public:
       

        struct Config
        {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;
            fhicl::Atom<art::InputTag>              comboHitCollection     {Name("ComboHitCollection"),     Comment("ComboHit collection {Name") };
            fhicl::Atom<art::InputTag>              strawHitFlagCollection {Name("StrawHitFlagCollection"), Comment("StrawHitFlag collection {Name") };
            fhicl::Atom<art::InputTag>              caloClusterCollection  {Name("CaloClusterCollection"),  Comment("Calo cluster collection {Name") };
            fhicl::Table<MVATools::Config>          tcMVA                  {Name("ClusterMVA"),             Comment("MVA for time cluster cleaning") }; 
            fhicl::Table<MVATools::Config>          tcCaloMVA              {Name("ClusterCaloMVA"),         Comment("MVA for time clsuter cleaning with calo") }; 
            fhicl::Sequence<std::string>            hsel                   {Name("HitSelectionBits"),       Comment("HitSelectionBits") }; 
            fhicl::Sequence<std::string>            hbkg                   {Name("HitBackgroundBits"),      Comment("HitBackgroundBits") }; 
            fhicl::Atom<bool>                       usecc                  {Name("UseCaloCluster"),         Comment("Use calorimeter cluster") }; 
            fhicl::Atom<bool>                       useccpos               {Name("UseCaloClusterPosition"), Comment("Use calorimeter cluster position") }; 
            fhicl::Atom<float>                      ccmine                 {Name("CaloClusterMinE"),        Comment("Minimum energy for calorimeter cluster") }; 
            fhicl::Atom<float>                      ccweight               {Name("CaloClusterWeight"),      Comment("Weight of cluster in tracker hits") }; 
            fhicl::Table<TrkTimeCalculator::Config> ttcalc                 {Name("T0Calculator"),           Comment("TimeTracker calculator config") };  
            fhicl::Atom<bool>                       testflag               {Name("TestFlag"),               Comment("Test hit flags") }; 
            fhicl::Atom<float>                      maxdt                  {Name("DtMax"),                  Comment("Maximum delta time for hit in cluster") }; 
            fhicl::Atom<unsigned>                   minnhits               {Name("MinNHits"),               Comment("Minimum number of hits for cluster") }; 
            fhicl::Atom<float>                      minkeepmva             {Name("MinKeepHitMVA"),          Comment("Minimum MVA score to keep in cluster") }; 
            fhicl::Atom<float>                      minaddmva              {Name("MinAddHitMVA"),           Comment("Minimum MVA score to add in cluster") }; 
            fhicl::Atom<float>                      maxdPhi                {Name("MaxdPhi"),                Comment("Maximum delta Phi for hit to be in cluster") }; 
            fhicl::Atom<float>                      tmin                   {Name("Tmin"),                   Comment("Time histogram start") }; 
            fhicl::Atom<float>                      tmax                   {Name("Tmax"),                   Comment("Time histogram end") }; 
            fhicl::Atom<float>                      tbin                   {Name("Tbin"),                   Comment("Time histogram bin width") }; 
            fhicl::Atom<float>                      pitch                  {Name("AveragePitch"),           Comment("Average helix pitch (= dz/dflight, =sin(lambda)") }; 
            fhicl::Atom<float>                      ymin                   {Name("Ymin"),                   Comment("Minimum hit in time histo bin for peak") }; 
            fhicl::Atom<bool>                       recover                {Name("RefineClusters"),         Comment("Apply hit refining algorithm") }; 
            fhicl::Atom<bool>                       refine                 {Name("PrefilterCluster"),       Comment("Apply hit pre-filtering algorithm") }; 
            fhicl::Atom<bool>                       preFilter              {Name("RecoverHits"),            Comment("Apply hit recovery algorithm") }; 
            fhicl::Atom<int>                        npeak                  {Name("PeakWidth"),              Comment("Time Peak Width") }; 
            fhicl::Atom<int>                        printfreq              {Name("printFrequency"),         Comment("Print frequency"), 100 }; 
            fhicl::Atom<int>                        debugLevel             {Name("debugLevel"),             Comment("Debut Level"), 0 }; 
        };

        explicit TimeClusterFinder(const art::EDProducer::Table<Config>& config);

        void beginJob() override;
        void produce(art::Event& e) override;

    
    private:
       typedef std::pair<Float_t,int> BinContent;
       typedef std::vector<StrawHitIndex>::iterator ISH;
       
       int                                             _iev; 
       const art::ProductToken<ComboHitCollection>     _chToken;
       const art::ProductToken<StrawHitFlagCollection> _shfToken;
       const art::ProductToken<CaloClusterCollection>  _ccToken;      
       const StrawHitFlagCollection* _shfcol;
       const ComboHitCollection*     _chcol;
       const CaloClusterCollection*  _cccol;
       StrawHitFlag                  _hsel;
       StrawHitFlag                  _hbkg;
       MVATools                      _tcMVA;     
       MVATools                      _tcCaloMVA; 
       bool                          _usecc, _useccpos;
       float                         _ccmine, _ccwt;
       TrkTimeCalculator             _ttcalc;
       bool                          _testflag;
       float                         _maxdt;
       unsigned                      _minnhits;
       float                         _minkeepmva, _minaddmva; 
       float                         _maxdPhi;
       float                         _tmin, _tmax, _tbin;
       float		             _pitch; 
       float                         _ymin;
       bool                          _refine;
       bool                          _preFilter;
       bool                          _recover;
       int                           _npeak;
       int                           _printfreq;
       int                           _debug;    
       TH1F                          _timespec;
       TimeCluMVA                    _pmva; // input variables to TMVA for cluster cleaning


      void findClusters(TimeClusterCollection& tccol);
      void findCaloSeeds(TimeClusterCollection& tccol, art::Handle<CaloClusterCollection> const& ccH);
      void fillTimeSpectrum();
      void initCluster(TimeCluster& tc);
      void prefilterCluster(TimeCluster& tc);
      void recoverHits(TimeCluster& tc);
      ISH  removeHit(TimeCluster& tc, ISH);
      void addHit(TimeCluster& tc,size_t iadd);
      void clusterMean(TimeCluster& tc);
      void refineCluster(TimeCluster& tc);
      void findPeaks(TimeClusterCollection& seeds);
      void assignHits(TimeClusterCollection& tccol );
      bool goodHit(const StrawHitFlag& flag) const;
  };

  
  TimeClusterFinder::TimeClusterFinder(const art::EDProducer::Table<Config>& config) :
     art::EDProducer{config},
     _chToken      { consumes<ComboHitCollection>(      config().comboHitCollection()) },
     _shfToken     { mayConsume<StrawHitFlagCollection>(config().strawHitFlagCollection()) },
     _ccToken      { mayConsume<CaloClusterCollection>( config().caloClusterCollection()) },
     _hsel         ( config().hsel()),
     _hbkg         ( config().hbkg()),
     _tcMVA        ( config().tcMVA()),
     _tcCaloMVA    ( config().tcCaloMVA()),
     _usecc        ( config().usecc()),
     _useccpos     ( config().useccpos()),
     _ccmine       ( config().ccmine()),
     _ccwt         ( config().ccweight()),
     _ttcalc       ( config().ttcalc()),
     _testflag     ( config().testflag()),
     _maxdt        ( config().maxdt()),
     _minnhits     ( config().minnhits()),
     _minkeepmva   ( config().minkeepmva()),
     _minaddmva    ( config().minaddmva()),
     _maxdPhi      ( config().maxdPhi()),
     _tmin         ( config().tmin()),
     _tmax         ( config().tmax()),
     _tbin         ( config().tbin()),
     _pitch        ( config().pitch()), 
     _ymin         ( config().ymin()),
     _refine       ( config().refine()),          
     _preFilter    ( config().preFilter()),        
     _recover      ( config().recover()),      
     _npeak        ( config().npeak()), 
     _printfreq    ( config().printfreq()),
     _debug        ( config().debugLevel())
    {
        unsigned nbins = (unsigned)rint((_tmax-_tmin)/_tbin);
        _timespec = TH1F("timespec","time spectrum",nbins,_tmin,_tmax);
        produces<TimeClusterCollection>();
    }

  void TimeClusterFinder::beginJob() {
    _tcMVA.initMVA();
    _tcCaloMVA.initMVA();
    if (_debug > 0)
    {
      std::cout << "TimeClusterFinder MVA : " << std::endl;
      _tcMVA.showMVA();
      std::cout << "TimeClusterFinder Calo MVA : " << std::endl;
      _tcCaloMVA.showMVA();
    }
  }


  //--------------------------------------------------------------------------------------------------------------
  void TimeClusterFinder::produce(art::Event & event ){
    _iev = event.id().event();

    if (_debug > 0 && (_iev%_printfreq)==0) std::cout<<"TimeClusterFinder: event="<<_iev<<std::endl;

    auto const& chH = event.getValidHandle(_chToken);
    _chcol = chH.product();

    art::Handle<CaloClusterCollection> ccH{}; // need to cache for later Ptr creation 
    if(_usecc){
      ccH = event.getHandle<CaloClusterCollection>(_ccToken);
      _cccol = ccH.product();
    }

    if(_testflag){
      auto shfH = event.getValidHandle(_shfToken);
      _shfcol = shfH.product();
      if(_shfcol->size() != _chcol->size())
	throw cet::exception("RECO")<<"TimeClusterFinder: inconsistent flag collection length " << endl;
    }

    std::unique_ptr<TimeClusterCollection> tccol(new TimeClusterCollection);
    // If requested, use calo clusters to for time cluster seeds
    if (_usecc) findCaloSeeds(*tccol,ccH);
    // find all the hit clusters
    findClusters(*tccol);

    if (_debug > 0) std::cout << "Found " << tccol->size() << " Time Clusters " << std::endl;

    if (_debug > 1){
      for(auto const& tc : *tccol) {
	std::cout << "Time Cluster time = " << tc.t0().t0() << " +- " << tc.t0().t0Err()
	  << " position = " << tc._pos << std::endl;
	if(_debug > 3){
	  for (auto shi : tc._strawHitIdxs ) {
	    std::cout << "Time Cluster hit at index " << shi << std::endl;
	  }
	}
      }
    }
    event.put(std::move(tccol));
  }


  //--------------------------------------------------------------------------------------------------------------
  void TimeClusterFinder::findClusters(TimeClusterCollection& tccol) {
    // find seed from hits
    fillTimeSpectrum();
    findPeaks(tccol);
    // associate hits to seeds
    assignHits(tccol);
    // loop over seeds and fill/refine information
    auto itc = tccol.begin();
    while(itc != tccol.end()){
      TimeCluster& tc = *itc;
      initCluster(tc);
      if (_preFilter) prefilterCluster(tc);
      if( tc.nStrawHits() >= _minnhits) {
	clusterMean(tc);
	if (_refine) refineCluster(tc);
	if (_recover) recoverHits(tc);
      }
      if (tc.nStrawHits() < _minnhits) {
	itc = tccol.erase(itc);
      } else
	++itc;
      //std::cout<<"Collection size final"<<tc._strawHitIdxs.size()<<std::endl;
    }
    // debug test of histogram
    if (_debug > 2) {
      art::ServiceHandle<art::TFileService> tfs;
      TH1F* tspec = tfs->make<TH1F>(_timespec);
      char name[40];
      char title[100];
      snprintf(name,40,"tspec_%i",_iev);
      snprintf(title,100,"time spectrum event %i;nsec",_iev);
      tspec->SetNameTitle(name,title);
    }
  }

  void TimeClusterFinder::findCaloSeeds(TimeClusterCollection& tccol, art::Handle<CaloClusterCollection>const& ccH) {
    for(size_t icalo=0; icalo < _cccol->size(); ++icalo){
      auto const& calo = (*_cccol)[icalo];
      if (calo.energyDep() > _ccmine){
	TimeCluster tc;
	tc._t0 = TrkT0(_ttcalc.caloClusterTime(calo,_pitch), _ttcalc.caloClusterTimeErr());
	tc._caloCluster = art::Ptr<CaloCluster>(ccH,icalo);
	tccol.push_back(tc);
      }
    }
  }

  //--------------------------------------------------------------------------------------------------------------
  void TimeClusterFinder::fillTimeSpectrum() {
    _timespec.Reset();
    for (unsigned istr=0; istr<_chcol->size();++istr) {
      if (_testflag && !goodHit((*_shfcol)[istr])) continue;
      ComboHit const& ch = (*_chcol)[istr];
      float time = _ttcalc.comboHitTime((*_chcol)[istr],_pitch);
      _timespec.Fill(time,ch.nStrawHits());
    }
  }

  void TimeClusterFinder::assignHits(TimeClusterCollection& tccol ) {
  // assign hits to the closest time peak
    for(size_t istr=0; istr<_chcol->size(); ++istr) {
      if ((!_testflag) || goodHit((*_shfcol)[istr])) {
	ComboHit const& ch =(*_chcol)[istr];
	float time = _ttcalc.comboHitTime(ch,_pitch);
	float mindt(1e5);
	auto besttc = tccol.end();
	// find the closest seed (if any)
	for (auto itc = tccol.begin(); itc != tccol.end(); ++itc) {
	  float dt = fabs(time - itc->_t0._t0);
	  // make an absolute cut, including error on the cluster t0
	  if (dt < _maxdt+itc->_t0._t0err && dt < mindt){
	    mindt = dt;
	    besttc = itc;
	  }
	}
	if(besttc != tccol.end())
	  besttc->_strawHitIdxs.push_back(istr);
      }
    }
  }

  //--------------------------------------------------------------------------------------------------------------
  void TimeClusterFinder::findPeaks(TimeClusterCollection& tccol) {
    int nbins = _timespec.GetNbinsX()+1;
    std::vector<bool> alreadyUsed(nbins,false);
    // blank out bins around input times (from calo clusters)
    for(auto const& tc : tccol ){ 
      int ibin = _timespec.FindBin(tc._t0._t0);
      for(int jbin = std::max(1,ibin-_npeak);jbin < std::min(nbins,ibin+_npeak+1); ++jbin)
	alreadyUsed[jbin] = true;
    }
    // loop over spectrum to find peaks 
    std::vector<BinContent> bcv;
    for (int ibin=1;ibin < nbins; ++ibin)
      if (_timespec.GetBinContent(ibin) >= _ymin) bcv.push_back(make_pair(_timespec.GetBinContent(ibin),ibin));
    std::sort(bcv.begin(),bcv.end(),[](const BinContent& x, const BinContent& y){return x.first > y.first;});

    for (const auto& bc : bcv) {
      if (alreadyUsed[bc.second]) continue;
      float nsh(0.0);
      float t0(0.0);
      for (int ibin = std::max(1,bc.second-_npeak);ibin < std::min(nbins,bc.second+_npeak+1); ++ibin) {
	nsh += _timespec.GetBinContent(ibin);
	t0 += _timespec.GetBinCenter(ibin)*_timespec.GetBinContent(ibin);
	alreadyUsed[ibin] = true;
      }
      t0 /= nsh;
      // if the count is enough, create a cluster
      if (nsh > _minnhits){
	TimeCluster tc;
	tc._t0 = TrkT0(t0,_tbin*0.5); // bin width
	tc._nsh = nsh;
	tccol.push_back(tc);
      }    
    }
  }

  //--------------------------------------------------------------------------------------------------------------
  void TimeClusterFinder::initCluster(TimeCluster& tc) {
    // use medians to initialize robustly
    accumulator_set<float, stats<tag::min > > tmin;
    accumulator_set<float, stats<tag::max > > tmax;
    accumulator_set<float, stats<tag::weighted_median(with_p_square_quantile) >, float > tacc, xacc, yacc, zacc;

    unsigned nstrs = tc._strawHitIdxs.size();
    tc._nsh = 0;
    for(auto ish :tc._strawHitIdxs) {
      if (_testflag && !goodHit((*_shfcol)[ish])) continue;
      ComboHit const& ch = (*_chcol)[ish];
      unsigned nsh = ch.nStrawHits();
      tc._nsh += nsh;
      const XYZVec& pos = ch.pos();
      float htime = _ttcalc.comboHitTime(ch,_pitch);
      float hwt = ch.nStrawHits();
      tmin(htime);
      tmax(htime);
      tacc(htime,weight=hwt);
      xacc(pos.x(),weight=hwt);
      yacc(pos.y(),weight=hwt);
      zacc(pos.z(),weight=hwt);
    }

    if (tc.hasCaloCluster()) {
    // don't update t0 if there's an assigned calo cluster
      if(_useccpos){
	xacc(tc._caloCluster->cog3Vector().x(),weight=_ccwt);
	yacc(tc._caloCluster->cog3Vector().y(),weight=_ccwt);
      }
    } else {
      static float invsqrt12(1.0/sqrt(12.0));
      tc._t0._t0 = extract_result<tag::weighted_median>(tacc);
      tc._t0._t0err = ( boost::accumulators::extract::max(tmax)-boost::accumulators::extract::min(tmin))*invsqrt12/sqrt(nstrs);
    }
    //
    tc._pos = XYZVec(extract_result<tag::weighted_median>(xacc),
        extract_result<tag::weighted_median>(yacc),
        extract_result<tag::weighted_median>(zacc));

    if (_debug > 0) std::cout<<"Init time peak "<<tc._t0._t0<<std::endl;
  }

  // prefilter based on a rough hemisphere cut and the initial robust position
  void TimeClusterFinder::prefilterCluster(TimeCluster& tc){
    bool changed(true);
    while (changed) {
      changed = false;
      float pphi = polyAtan2( tc._pos.y(), tc._pos.x());
      auto iworst = tc._strawHitIdxs.end();
      float maxadPhi(_maxdPhi);
      for( auto ips = tc._strawHitIdxs.begin(); ips != tc._strawHitIdxs.end(); ++ips){
	ComboHit const& ch = (*_chcol)[*ips];
	float phi   = polyAtan2(ch.pos().y(), ch.pos().x()); 
	float dphi  = Angles::deltaPhi(phi,pphi);
	float adphi = std::abs(dphi);
	if(adphi > maxadPhi ){
	  iworst = ips;
	  maxadPhi = adphi;
	}
      }
      if( iworst != tc._strawHitIdxs.end()){
	changed = true;
	removeHit(tc,iworst);
      }
    }
  }

  void TimeClusterFinder::recoverHits(TimeCluster& tc){
    bool changed(true);
    while (changed) {
      changed = false;
      float pphi = polyAtan2(tc._pos.y(), tc._pos.x());
      for(size_t ich=0;ich < _chcol->size(); ++ich){
	if ((!_testflag) || goodHit((*_shfcol)[ich])) {
	  if(std::find(tc._strawHitIdxs.begin(),tc._strawHitIdxs.end(),ich) == tc._strawHitIdxs.end()){
	    ComboHit const& ch = (*_chcol)[ich];
	    float cht = _ttcalc.comboHitTime(ch,_pitch);
	    _pmva._dt = fabs(cht - tc._t0._t0);
	    if(_pmva._dt < _maxdt+tc._t0._t0err){
	      float phi = polyAtan2(ch.pos().y(), ch.pos().x());//ch.phi();
	      float dphi = fabs(Angles::deltaPhi(phi,pphi));
	      if(dphi < _maxdPhi){ 
		_pmva._dphi = dphi;
		_pmva._rho = ch.pos().Perp2();
		_pmva._nsh = ch.nStrawHits();
		_pmva._plane = ch.strawId().plane();
		_pmva._werr = ch.wireRes();
		_pmva._wdist = fabs(ch.wireDist());

		float mvaout(-1.0);
		if (tc.hasCaloCluster())
		  mvaout = _tcCaloMVA.evalMVA(_pmva._pars);
		else
		  mvaout = _tcMVA.evalMVA(_pmva._pars);
		if (mvaout > _minaddmva) {
		  addHit(tc,ich);
		  changed = true;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  std::vector<StrawHitIndex>::iterator TimeClusterFinder::removeHit(TimeCluster& tc, ISH iworst) {
    ComboHit const& ch = (*_chcol)[*iworst];
    unsigned nsh = ch.nStrawHits();
    float denom = float(tc._nsh - nsh);
    // update time cluster properties 
    if(!tc.hasCaloCluster()){
      float cht = _ttcalc.comboHitTime(ch,_pitch);
      float newt0  = (tc._t0._t0*tc._nsh - cht*nsh)/denom;
      tc._t0._t0err = sqrt((tc._t0._t0err*tc._t0._t0err*tc._nsh - (cht-newt0)*(cht-tc._t0._t0)*nsh )/denom);
      tc._t0._t0 = newt0;
    }
    tc._pos.SetX((tc._pos.x()*tc._nsh - ch.pos().x()*nsh)/denom);
    tc._pos.SetY((tc._pos.y()*tc._nsh - ch.pos().x()*nsh)/denom);
    tc._pos.SetZ((tc._pos.z()*tc._nsh - ch.pos().x()*nsh)/denom);
    tc._nsh -= nsh;
    return tc._strawHitIdxs.erase(iworst);
  }

  void TimeClusterFinder::addHit(TimeCluster& tc,size_t iadd) {
    ComboHit const& ch = (*_chcol)[iadd];
    unsigned nsh = ch.nStrawHits();
    float denom = float(tc._nsh + nsh);
    // update time cluster properties 
    if(!tc.hasCaloCluster()){
      float cht = _ttcalc.comboHitTime(ch,_pitch);
      float newt0  = (tc._t0._t0*tc._nsh + cht*nsh)/denom;
      tc._t0._t0err = sqrt((tc._t0._t0err*tc._t0._t0err*tc._nsh + (cht-newt0)*(cht-tc._t0._t0)*nsh )/denom);
      tc._t0._t0 = newt0;
    }
    tc._pos.SetX((tc._pos.x()*tc._nsh + ch.pos().x()*nsh)/denom);
    tc._pos.SetY((tc._pos.y()*tc._nsh + ch.pos().x()*nsh)/denom);
    tc._pos.SetZ((tc._pos.z()*tc._nsh + ch.pos().x()*nsh)/denom);
    tc._nsh += nsh;
    tc._strawHitIdxs.push_back(iadd);
  }

  void TimeClusterFinder::clusterMean(TimeCluster& tc) {
    // compute properties using weighted mean
    accumulator_set<float, stats<tag::weighted_variance(lazy)>, float > terr;
    accumulator_set<float, stats<tag::weighted_mean >,float > xacc, yacc, zacc;
    for(StrawHitIndex ish : tc._strawHitIdxs) {
      ComboHit const& ch = (*_chcol)[ish];
      float hwt = ch.nStrawHits();
      float cht = _ttcalc.comboHitTime(ch,_pitch);
      terr(cht,weight=hwt);
      xacc(ch.pos().x(),weight=hwt);
      yacc(ch.pos().y(),weight=hwt);
      zacc(ch.pos().z(),weight=hwt);
    }
    if (tc.hasCaloCluster()) {
      if(_useccpos){
	xacc(tc._caloCluster->cog3Vector().x(),weight=_ccwt);
	yacc(tc._caloCluster->cog3Vector().y(),weight=_ccwt);
      }
    } else {
      tc._t0._t0 = extract_result<tag::weighted_mean>(terr);
      tc._t0._t0err = sqrtf(std::max(double(1.0),2.0*extract_result<tag::weighted_variance(lazy)>(terr))/extract_result<tag::count>(terr));
    }
    
    tc._pos = XYZVec(extract_result<tag::weighted_mean>(xacc),
	extract_result<tag::weighted_mean>(yacc),
	extract_result<tag::weighted_mean>(zacc));
  } 

  void TimeClusterFinder::refineCluster(TimeCluster& tc) {
    // mva filtering; remove worst hit iteratively
    bool changed = true;
    while (changed) {
      changed = false;
      auto iworst = tc._strawHitIdxs.end();
      float worstmva(100.0);
      float pphi = polyAtan2(tc._pos.y(), tc._pos.x());
      for (auto ips=tc._strawHitIdxs.begin();ips != tc._strawHitIdxs.end();++ips) {
        ComboHit const& ch = (*_chcol)[*ips];
        float cht = _ttcalc.comboHitTime(ch,_pitch);

        _pmva._dt = fabs(cht - tc._t0._t0);
        float phi = polyAtan2(ch.pos().y(), ch.pos().x());//ch.phi();
        float dphi = Angles::deltaPhi(phi,pphi);
        _pmva._dphi = fabs(dphi);
	_pmva._rho = ch.pos().Perp2();
	_pmva._nsh = ch.nStrawHits();
	_pmva._plane = ch.strawId().plane();
	_pmva._werr = ch.wireRes();
	_pmva._wdist = fabs(ch.wireDist());

	float mvaout(-1.0);
	if (tc.hasCaloCluster())
	   mvaout = _tcCaloMVA.evalMVA(_pmva._pars);
	else
	  mvaout = _tcMVA.evalMVA(_pmva._pars);
	if (mvaout < worstmva) {
	  worstmva = mvaout;
	  iworst = ips;
        }
      }

      if (worstmva < _minkeepmva) {
        changed = true;
	removeHit(tc,iworst);
      }
    }
  }

  bool TimeClusterFinder::goodHit(const StrawHitFlag& flag) const
  {
    return flag.hasAllProperties(_hsel) && !flag.hasAnyProperty(_hbkg);
  }

}

using mu2e::TimeClusterFinder;
DEFINE_ART_MODULE(TimeClusterFinder);

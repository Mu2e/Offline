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
#include "Mu2eUtilities/inc/polyAtan2.hh"
// data
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
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

namespace {

  struct TimePeakMVA
  {
    vector<Float_t> _pars;
    Float_t& _dt;
    Float_t& _dphi;
    Float_t& _rho;
    TimePeakMVA() : _pars(3,0.0), _dt(_pars[0]), _dphi(_pars[1]), _rho(_pars[2]) {}
  };
}

namespace mu2e {
  class TimeClusterFinder : public art::EDProducer {
    public:
      typedef std::pair<Float_t,int> BinContent;
      enum Mode{flag=0,filter};
      explicit TimeClusterFinder(fhicl::ParameterSet const& pset);

      void beginJob() override;
      void produce(art::Event& e) override;

    private:
      int               _iev;
      int               _debug;
      int               _printfreq;
      bool		_testflag;
      art::ProductToken<ComboHitCollection> const _chToken;
      art::ProductToken<StrawHitFlagCollection> const _shfToken;
      art::ProductToken<CaloClusterCollection> const _ccToken;
      const StrawHitFlagCollection *_shfcol;
      const ComboHitCollection *_chcol;
      const CaloClusterCollection *_cccol;
      StrawHitFlag      _hsel, _hbkg;
      float             _maxdt;
      unsigned          _minnhits;
      float             _minpeakmva, _maxpeakdt;
      float             _maxdPhi;
      float             _tmin, _tmax, _tbin;
      TH1F              _timespec;
      float             _ymin;
      bool              _refine;
      bool              _preFilter;
      bool              _finalDt;
      bool              _usecc, _useccpos;
      float             _ccmine, _ccwt;
      float             _maxover;
      MVATools          _tcMVA; // MVA for peak cleaning
      MVATools          _tcCaloMVA; // MVA for peak cleaning, with calo cluster
      TimePeakMVA       _pmva; // input variables to TMVA for peak cleaning
      TrkTimeCalculator _ttcalc;
      int               _npeak;


      void findClusters(TimeClusterCollection& tccol, art::Handle<CaloClusterCollection> const& ccH);
      void fillTimeSpectrum();
      void initCluster(TimeCluster& tc);
      void prefilterCluster(TimeCluster& tc);
      void removeHit(TimeCluster& tc,size_t iremove);
      void clusterMean(TimeCluster& tc);
      void refineCluster(TimeCluster& tc);
      void clusterDt(TimeCluster& tc);
      void findPeaks(TimeClusterCollection& seeds);
      void assignHits(TimeClusterCollection& tccol );
      bool goodHit(const StrawHitFlag& flag) const;
  };

  TimeClusterFinder::TimeClusterFinder(fhicl::ParameterSet const& pset) :
    _debug             (pset.get<int>("debugLevel",0)),
    _printfreq         (pset.get<int>("printFrequency",101)),
    _testflag(pset.get<bool>("TestFlag")),
    _chToken{consumes<ComboHitCollection>(pset.get<art::InputTag>("ComboHitCollection"))},
    _shfToken{mayConsume<StrawHitFlagCollection>(pset.get<art::InputTag>("StrawHitFlagCollection"))},
    _ccToken{mayConsume<CaloClusterCollection>(pset.get<art::InputTag>("CaloClusterCollection"))},
    _hsel              (pset.get<std::vector<std::string> >("HitSelectionBits",vector<string>{"EnergySelection","TimeSelection","RadiusSelection"})),
    _hbkg              (pset.get<vector<string> >("HitBackgroundBits",vector<string>{"Background"})),
    _maxdt             (pset.get<float>(  "DtMax",30.0)),
    _minnhits          (pset.get<unsigned>("MinNHits",10)),
    _minpeakmva        (pset.get<float>(  "MinTimePeakMVA",0.2)),
    _maxpeakdt         (pset.get<float>(  "MaxTimePeakDeltat",25.0)),
    _maxdPhi           (pset.get<float>(  "MaxdPhi",1.5)),
    _tmin              (pset.get<float>(  "tmin",450.0)),
    _tmax              (pset.get<float>(  "tmax",1700.0)),
    _tbin              (pset.get<float>(  "tbin",15.0)),
    _ymin              (pset.get<float>(  "ymin",5.0)),
    _refine            (pset.get<bool>(  "RefineClusters",true)),
    _preFilter         (pset.get<bool>(    "PrefilterCluster",true)),
    _finalDt	       (pset.get<bool>(    "FinalDtCut",true)),
    _usecc             (pset.get<bool>(    "UseCaloCluster",false)),
    _useccpos          (pset.get<bool>(    "UseCaloClusterPosition",false)),
    _ccmine            (pset.get<float>(  "CaloClusterMinE",50.0)),
    _ccwt              (pset.get<float>(  "CaloClusterWeight",5.0)),
    _tcMVA           (pset.get<fhicl::ParameterSet>("ClusterMVA",fhicl::ParameterSet())),
    _tcCaloMVA           (pset.get<fhicl::ParameterSet>("ClusterCaloMVA",fhicl::ParameterSet())),
    _ttcalc            (pset.get<fhicl::ParameterSet>("T0Calculator",fhicl::ParameterSet())),
    _npeak       (pset.get<int>("PeakWidth",1)) // # of bins
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
      event.getByToken(_ccToken, ccH);
      _cccol = ccH.product();
    }

    if(_testflag){
      auto shfH = event.getValidHandle(_shfToken);
      _shfcol = shfH.product();
      if(_shfcol->size() != _chcol->size())
	throw cet::exception("RECO")<<"TimeClusterFinder: inconsistent flag collection length " << endl;
    }

    std::unique_ptr<TimeClusterCollection> tccol(new TimeClusterCollection);
    findClusters(*tccol,ccH);

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
  void TimeClusterFinder::findClusters(TimeClusterCollection& tccol, art::Handle<CaloClusterCollection>const& ccH) {
    // find seeds from calo clusters
    if (_usecc) {
      for(size_t icalo=0; icalo < _cccol->size(); ++icalo){
	auto const& calo = (*_cccol)[icalo];
	if (calo.energyDep() > _ccmine){
	  TimeCluster tc;
	  tc._t0 = TrkT0(_ttcalc.caloClusterTime(calo), _ttcalc.caloClusterTimeErr(calo.diskId()));
	  tc._caloCluster = art::Ptr<CaloCluster>(ccH,icalo);
	  tccol.push_back(tc);
	}
      }
    }
    // find seed from hits
    fillTimeSpectrum();
    findPeaks(tccol);
    // associate hits to seeds
    assignHits(tccol);
    // loop over seeds and fill/refine information
    for (auto itc = tccol.begin(); itc < tccol.end(); ++itc){
      TimeCluster& tc = *itc;
      initCluster(tc);
      if (_refine && tc.nStrawHits() >= _minnhits) refineCluster(tc);
      if (_finalDt && tc.nStrawHits() >= _minnhits) clusterDt(tc);
// cleanup
      if (tc.nStrawHits() < _minnhits) {
	std::swap(tc,tccol.back());	
	tccol.pop_back();
      }
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

  //--------------------------------------------------------------------------------------------------------------
  void TimeClusterFinder::fillTimeSpectrum() {
    _timespec.Reset();
    for (unsigned istr=0; istr<_chcol->size();++istr) {
      if (_testflag && !goodHit((*_shfcol)[istr])) continue;
      ComboHit const& ch = (*_chcol)[istr];
      // if ( (ch.time() < _minT) || (ch.time() > _maxT) )  continue;
      float time = _ttcalc.comboHitTime((*_chcol)[istr]);
      _timespec.Fill(time,ch.nStrawHits());
    }
  }

  void TimeClusterFinder::assignHits(TimeClusterCollection& tccol ) {
  // assign hits to the closest time peak
    for(size_t istr=0; istr<_chcol->size(); ++istr) {
      if ((!_testflag) || goodHit((*_shfcol)[istr])) {
	ComboHit const& ch =(*_chcol)[istr];
	float time = _ttcalc.comboHitTime(ch);
	float mindt(_maxdt);
	auto besttc = tccol.end();
	// find the closest seed (if any)
	for (auto itc = tccol.begin(); itc != tccol.end(); ++itc) {
	  float dt = fabs(time - itc->_t0._t0);
	  if (dt < mindt){
	    mindt = dt;
	    besttc = itc;
	  }
	}
	if(besttc != tccol.end()) besttc->_strawHitIdxs.push_back(istr);
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
	static const float isqrt12 = 1.0/sqrt(12.0);
	tc._t0 = TrkT0(t0,_tbin*isqrt12);
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
      float htime = _ttcalc.comboHitTime(ch);
      float pwt = ch.nStrawHits();
      tmin(htime);
      tmax(htime);
      tacc(htime,weight=pwt);
      xacc(pos.x(),weight=pwt);
      yacc(pos.y(),weight=pwt);
      zacc(pos.z(),weight=pwt);
    }

    if (tc._caloCluster.isNonnull()) {
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
    // hemisphere filter
    if (_preFilter) prefilterCluster(tc);
  }

  // prefilter based on a rough hemisphere cut
  void TimeClusterFinder::prefilterCluster(TimeCluster& tc){
    clusterMean(tc);
    bool enoughhits(true);
    bool changed(true);
    //prefilter one after another
    while (enoughhits && changed) {
      float pphi = polyAtan2(tc._pos.y(), tc._pos.x());
      changed = false;
      size_t iworst(0);
      float maxadPhi(0);
      for (size_t ips = 0; ips <  tc._strawHitIdxs.size(); ++ips){
	StrawHitIndex ish = tc._strawHitIdxs[ips];
	ComboHit const& ch = (*_chcol)[ish];
	float phi   = polyAtan2(ch.pos().y(), ch.pos().x()); 
	float dphi  = Angles::deltaPhi(phi,pphi);
	float adphi = std::abs(dphi);
	if (adphi > maxadPhi) {iworst=ips; maxadPhi=adphi;}
      }
      if (maxadPhi>_maxdPhi){
	changed = true;
	removeHit(tc,iworst);
      }
      enoughhits = tc._nsh >= _minnhits;
    }
  }

  void TimeClusterFinder::removeHit(TimeCluster& tc,size_t iworst) {
    std::swap(tc._strawHitIdxs[iworst],tc._strawHitIdxs.back());
    ComboHit const& ch = (*_chcol)[tc._strawHitIdxs.back()];
    unsigned nsh = ch.nStrawHits();
    float denom = float(tc._nsh - nsh);
    float cht = _ttcalc.comboHitTime(ch);
    // update time cluster properties 
    if(!tc._caloCluster.isNonnull())tc._t0._t0 = (tc._t0._t0*tc._nsh - cht*nsh)/denom;
    // update t0 error too FIXME!
    tc._pos.SetX((tc._pos.x()*tc._nsh - ch.pos().x()*nsh)/denom);
    tc._pos.SetY((tc._pos.y()*tc._nsh - ch.pos().x()*nsh)/denom);
    tc._pos.SetZ((tc._pos.z()*tc._nsh - ch.pos().x()*nsh)/denom);
    tc._nsh -= nsh;
    tc._strawHitIdxs.pop_back();
  }

  void TimeClusterFinder::clusterMean(TimeCluster& tc) {
    // compute properties using weighted mean
    accumulator_set<float, stats<tag::weighted_variance(lazy)>, unsigned > terr;
    accumulator_set<float, stats<tag::weighted_mean >,unsigned > xacc, yacc, zacc;
    for(StrawHitIndex ish : tc._strawHitIdxs) {
      ComboHit const& ch = (*_chcol)[ish];
      unsigned nsh = ch.nStrawHits();
      float cht = _ttcalc.comboHitTime(ch);
      terr(cht,weight=nsh);
      xacc(ch.pos().x(),weight=nsh);
      yacc(ch.pos().y(),weight=nsh);
      zacc(ch.pos().z(),weight=nsh);
    }
    if (tc._caloCluster.isNonnull()) {
      if(_useccpos){
	xacc(tc._caloCluster->cog3Vector().x(),weight=_ccwt);
	yacc(tc._caloCluster->cog3Vector().y(),weight=_ccwt);
      }
    } else {
      tc._t0._t0 = extract_result<tag::weighted_mean>(terr);
      tc._t0._t0err = sqrtf(std::max(float(0.0),extract_result<tag::weighted_variance(lazy)>(terr))/extract_result<tag::count>(terr));
    }
    
    tc._pos = XYZVec(extract_result<tag::weighted_mean>(xacc),
	extract_result<tag::weighted_mean>(yacc),
	extract_result<tag::weighted_mean>(zacc));
  } 

  void TimeClusterFinder::refineCluster(TimeCluster& tc) {
    // convert to means
    clusterMean(tc);
    // mva filtering; remove worst hit iteratively
    bool changed = true;
    bool enoughhits = true;
    while (enoughhits && changed) {
      changed = false;
      size_t iworst(0);
      float worstmva(100.0);
      float pphi = polyAtan2(tc._pos.y(), tc._pos.x());
      for (size_t ips=0;ips<tc._strawHitIdxs.size();++ips) {
        unsigned ish = tc._strawHitIdxs[ips];
        ComboHit const& ch = (*_chcol)[ish];
        float cht = _ttcalc.comboHitTime(ch);

        float rho = sqrtf(ch.pos().Perp2());
        float phi = polyAtan2(ch.pos().y(), ch.pos().x());//ch.phi();
        float dphi = Angles::deltaPhi(phi,pphi);

        _pmva._dt = cht - tc._t0._t0;
        _pmva._dphi = dphi;
        _pmva._rho = rho; // change this to use radius^2 FIXME!

	float mvaout(-1.0);
	if (tc._caloCluster.isNonnull())
	   mvaout = _tcMVA.evalMVA(_pmva._pars);
	else
	   mvaout = _tcCaloMVA.evalMVA(_pmva._pars);
        if (mvaout < worstmva) {
          worstmva = mvaout;
          iworst = ips;
        }
      }

      if (worstmva < _minpeakmva) {
        changed = true;
	removeHit(tc,iworst);
      }
      enoughhits = tc._nsh >= _minnhits;
    }
  }

  void TimeClusterFinder::clusterDt(TimeCluster& tc) {
    // hard cut on dt
    std::vector<size_t> toremove;
    for (size_t ips = 0; ips <  tc._strawHitIdxs.size(); ++ips){
      StrawHitIndex ish = tc._strawHitIdxs[ips];
      ComboHit const& ch = (*_chcol)[ish];
      float cht = _ttcalc.comboHitTime(ch);
      float  dt = cht - tc._t0._t0;
      if (fabs(dt) > _maxpeakdt) {
	toremove.push_back(ips);
      }
    }
    clusterMean(tc);
  }

  bool TimeClusterFinder::goodHit(const StrawHitFlag& flag) const
  {
    return flag.hasAllProperties(_hsel) && !flag.hasAnyProperty(_hbkg);
  }


}

using mu2e::TimeClusterFinder;
DEFINE_ART_MODULE(TimeClusterFinder);

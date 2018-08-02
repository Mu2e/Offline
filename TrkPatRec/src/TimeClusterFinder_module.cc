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

  // should switch to float FIXME!
  struct TimePeakMVA
  {
    vector<Float_t> _pars;
    Float_t& _dt;
    Float_t& _dphi;
    Float_t& _rho;
    // should add number of hits in the cluster FIXME!!!
    TimePeakMVA() : _pars(3,0.0), _dt(_pars[0]), _dphi(_pars[1]), _rho(_pars[2]) {}
  };


  struct meanAccumulator
  {
    meanAccumulator(): sum(0),weight(0) {};
    void add(float m, float w) {sum +=m*w; weight +=w;}
    void remove(float m, float w) {sum -=m*w; weight -=w;}
    float mean() {return weight>0 ? sum/weight : 0;}
    float sum;
    float weight;
  };
}

namespace mu2e {
  class TimeClusterFinder : public art::EDProducer {
    public:
      enum ClusterAlgorithm {peak=0,scan};
      enum Mode{flag=0,filter};
      explicit TimeClusterFinder(fhicl::ParameterSet const& pset);

      void beginJob() override;
      void produce(art::Event& e) override;

    private:
      int               _iev;
      int               _debug;
      int               _printfreq;
      ClusterAlgorithm  _algo;
      bool		_testflag;
      art::ProductToken<ComboHitCollection> const _chToken;
      art::ProductToken<StrawHitFlagCollection> const _shfToken;
      art::InputTag const _ccTag;
      const StrawHitFlagCollection *_shfcol;
      const ComboHitCollection *_chcol;
      const CaloClusterCollection *_cccol;
      // FIXME: Caching an art::Handle is indicative of a poor design.
      // It should be rethought.
      art::Handle<CaloClusterCollection> _ccH{}; // invalid, used for creating Ptrs
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
      bool              _usecc;
      float             _ccmine, _ccwt;
      float             _maxover;
      MVATools          _peakMVA; // MVA for peak cleaning
      TimePeakMVA       _pmva; // input variables to TMVA for peak cleaning
      TrkTimeCalculator _ttcalc;
      int               _deltaNbins;

      typedef std::pair<Float_t,int> BinContent;

      void findClusters(TimeClusterCollection& tccol);
      void fillTimeSpectrum();
      void initCluster(TimeCluster& tc);
      void refineCluster(TimeCluster& tc);
      void findPeaks(std::vector<float>& tctimes);
      void scanPeaks(std::vector<float>& tctimes);
      bool goodHit(const StrawHitFlag& flag) const;
      void addCaloCluster(TimeCluster& tc);
  };

  TimeClusterFinder::TimeClusterFinder(fhicl::ParameterSet const& pset) :
    _debug             (pset.get<int>("debugLevel",0)),
    _printfreq         (pset.get<int>("printFrequency",101)),
    _algo              (static_cast<ClusterAlgorithm>(pset.get<int>("ClusterAlgorithm",peak))),
    _testflag(pset.get<bool>("TestFlag")),
    _chToken{consumes<ComboHitCollection>(pset.get<art::InputTag>("ComboHitCollection"))},
    _shfToken{mayConsume<StrawHitFlagCollection>(pset.get<art::InputTag>("StrawHitFlagCollection"))},
    _ccTag{pset.get<art::InputTag>("caloClusterModuleLabel","CaloClusterFast")},
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
    _usecc             (pset.get<bool>(    "UseCaloCluster",false)),
    _ccmine            (pset.get<float>(  "CaloClusterMinE",50.0)),
    _ccwt              (pset.get<float>(  "CaloClusterWeight",10.0)),
    _maxover           (pset.get<float>(  "MaxOverlap",0.3)),         // Maximum hit overlap to consider clusters as different
    _peakMVA           (pset.get<fhicl::ParameterSet>("PeakCleanMVA",fhicl::ParameterSet())),
    _ttcalc            (pset.get<fhicl::ParameterSet>("T0Calculator",fhicl::ParameterSet()))
    {
      mayConsume<CaloClusterCollection>(_ccTag);

      unsigned nbins = (unsigned)rint((_tmax-_tmin)/_tbin);
      _timespec = TH1F("timespec","time spectrum",nbins,_tmin,_tmax);
      _deltaNbins = int(_maxdt/_tbin)-1;

      produces<TimeClusterCollection>();

    }

  void TimeClusterFinder::beginJob()
  {
    _peakMVA.initMVA();
    if (_debug > 0)
    {
      std::cout << "TimeClusterFinder MVA : " << std::endl;
      _peakMVA.showMVA();
    }
  }


  //--------------------------------------------------------------------------------------------------------------
  void TimeClusterFinder::produce(art::Event & event ){
    _iev = event.id().event();

    if (_debug > 0 && (_iev%_printfreq)==0) std::cout<<"TimeClusterFinder: event="<<_iev<<std::endl;

    auto const& chH = event.getValidHandle(_chToken);
    _chcol = chH.product();

    if(_usecc){
      if(!event.getByLabel(_ccTag, _ccH))
        throw cet::exception("RECO")<<"TimeClusterFinder: No CaloCluster collection found for tag" <<  _ccTag << endl;
      _cccol = _ccH.product();
    }

    if(_testflag){
      auto shfH = event.getValidHandle(_shfToken);
      _shfcol = shfH.product();
      if(_shfcol->size() != _chcol->size())
        throw cet::exception("RECO")<<"TimeClusterFinder: inconsistent flag collection length " << endl;
    }

    std::unique_ptr<TimeClusterCollection> tccol(new TimeClusterCollection);
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
  void TimeClusterFinder::findClusters(TimeClusterCollection& tccol)
  {
    fillTimeSpectrum();
    vector<float> tctimes;
    switch (_algo ) {
      case peak : default:
        findPeaks(tctimes);
        break;
      case scan :
        scanPeaks(tctimes);
        break;
    }
    for (auto tctime: tctimes) {
      if(_debug > 1) std::cout << "Peak Time " << tctime  << std::endl;
      TimeCluster tc;
      tc._t0 = TrkT0(tctime,1.0);
      for(size_t istr=0; istr<_chcol->size(); ++istr) {
        if ((!_testflag) || goodHit((*_shfcol)[istr])) {
          ComboHit const& ch =(*_chcol)[istr];
          float time = _ttcalc.comboHitTime(ch);
          if (fabs(time-tctime) < _maxdt){
            tc._strawHitIdxs.push_back(StrawHitIndex(istr));
            tc._nsh += ch.nStrawHits();
          }
        }
      }

      if (_usecc) addCaloCluster(tc);

      initCluster(tc);
      if (_refine && tc.nStrawHits() >= _minnhits) refineCluster(tc);

      if (tc.nStrawHits() >= _minnhits) {
        bool overl(false);
        for (auto itc = tccol.begin(); itc < tccol.end(); ++itc)
        {
          if (TrkUtilities::overlap(tc,*itc) < _maxover) continue;
          if (tc.nStrawHits() > itc->nStrawHits())
          {
            tccol.erase(itc);
            break;
          } else {
            overl = true;
            break;
          }
        }
        if (!overl) tccol.push_back(tc);
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
      // if (ch.energyDep() > _maxElectronHitEnergy)         continue;
      // if ( (ch.time() < _minT) || (ch.time() > _maxT) )  continue;
      float time = _ttcalc.comboHitTime((*_chcol)[istr]);
      _timespec.Fill(time,ch.nStrawHits());
    }

    // weight the cluster WRT hits by an ad-hoc value.  This is more about signal/noise than resolution
    if (_usecc) {
      for(auto const& calo : *_cccol) {
        if (calo.energyDep() < _ccmine) continue;
        float time = _ttcalc.caloClusterTime(calo);
        _timespec.Fill(time, _ccwt);
      }
    }
  }

  //--------------------------------------------------------------------------------------------------------------
  void TimeClusterFinder::findPeaks(std::vector<float>& tctimes) {
    tctimes.clear();

    std::vector<BinContent> bcv;
    int nbins = _timespec.GetNbinsX()+1;
    std::vector<bool> alreadyUsed(nbins,false);

    for (int ibin=1;ibin < nbins; ++ibin)
      if (_timespec.GetBinContent(ibin) >= _ymin) bcv.push_back(make_pair(_timespec.GetBinContent(ibin),ibin));
    std::sort(bcv.begin(),bcv.end(),[](const BinContent& x, const BinContent& y){return x.first > y.first;});

    for (const auto& bc : bcv) {
      if (alreadyUsed[bc.second]) continue;

      float tctime(0.0);
      float norm(0.0);
      for (int ibin = std::max(1,bc.second-_deltaNbins);ibin < std::min(nbins,bc.second+_deltaNbins+1); ++ibin) {
        norm += _timespec.GetBinContent(ibin);
        tctime += _timespec.GetBinCenter(ibin)*_timespec.GetBinContent(ibin);
        alreadyUsed[ibin] = true;
      }
      tctime /= norm;
      if (norm > _minnhits) tctimes.push_back(tctime);
    }
  }

  //--------------------------------------------------------------------------------------------------------------
  void TimeClusterFinder::scanPeaks(std::vector<float>& tctimes) {
    //implement this FIXME!!
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
      // if (ch.energyDep() > _maxElectronHitEnergy)         continue;
      // if ( (ch.time() < _minT) || (ch.time() > _maxT) )  continue;
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
      float ctime = _ttcalc.caloClusterTime(*tc._caloCluster);
      tacc(ctime,weight=_ccwt);
      // add calo cluster position FIXME!!!
    }

    static float invsqrt12(1.0/sqrt(12.0));
    tc._t0._t0 = extract_result<tag::weighted_median>(tacc);
    tc._t0._t0err = ( boost::accumulators::extract::max(tmax)-boost::accumulators::extract::min(tmin))*invsqrt12/sqrt(nstrs);
    tc._pos = XYZVec(extract_result<tag::weighted_median>(xacc),
        extract_result<tag::weighted_median>(yacc),
        extract_result<tag::weighted_median>(zacc));

    if (_debug > 0) std::cout<<"Init time peak "<<tc._t0._t0<<std::endl;
  }

  void TimeClusterFinder::refineCluster(TimeCluster& tc)
  {
    tc._nsh = 0;

    float pphi(tc._pos.phi());
    float ptime(tc._t0.t0());

    bool enoughhits(true);
    bool changed(true);
    if (_preFilter) {
      //prefilter one after another
      while (enoughhits && changed) {
        changed = false;
        int iworst(-1);
        float maxadPhi(0);
        float worstphi(0),sumphi(0), sumt(0);
        unsigned nhits(0);
        for (size_t ips=0; ips<tc._strawHitIdxs.size(); ++ips) {
          unsigned ish = tc._strawHitIdxs[ips];
          ComboHit const& ch = (*_chcol)[ish];
          unsigned nsh = ch.nStrawHits();
          float phi   = ch.phi();  // should be float FIXME!
          float dphi  = Angles::deltaPhi(phi,pphi);
          float adphi = std::abs(dphi);
          sumphi += nsh*phi;
          sumt += nsh *_ttcalc.comboHitTime(ch);
          nhits += nsh;
          if (adphi > maxadPhi) {iworst=ips; maxadPhi=adphi; worstphi=phi;}
        }
        if (maxadPhi>_maxdPhi){
          changed = true;
// remove the worst hit
          std::swap(tc._strawHitIdxs[iworst],tc._strawHitIdxs.back());
          ComboHit const& ch = (*_chcol)[tc._strawHitIdxs.back()];
          unsigned nsh = ch.nStrawHits();
          nhits -= nsh;
          sumphi -= nsh * worstphi;
          sumt -= nsh *_ttcalc.comboHitTime(ch);
          tc._strawHitIdxs.pop_back();
        }
        enoughhits = nhits >= _minnhits;
        pphi  = sumphi/float(nhits);
        ptime = sumt/float(nhits);
      }
    }
    // mva filtering
    changed = true;
    while (enoughhits && changed) {
      changed = false;
      size_t iworst(0);
      float worstmva(100.0);
      float worstphi(0.0);
      float sumphi(0), sumt(0);
      unsigned nhits(0);
      for (size_t ips=0;ips<tc._strawHitIdxs.size();++ips) {
        unsigned ish = tc._strawHitIdxs[ips];
        ComboHit const& ch = (*_chcol)[ish];
        unsigned nsh = ch.nStrawHits();
        float time = _ttcalc.comboHitTime(ch);

        float rho = sqrtf(ch.pos().Perp2());
        float phi = ch.phi();
        float dphi = Angles::deltaPhi(phi,pphi);
        sumphi += nsh*phi;
        sumt += nsh*time;

        _pmva._dt = time - ptime;
        _pmva._dphi = dphi;
        _pmva._rho = rho; // change this to use radius^2 FIXME!

        float mvaout = _peakMVA.evalMVA(_pmva._pars);
        if (mvaout < worstmva) {
          worstmva = mvaout;
          worstphi = phi;
          iworst = ips;
        }
        nhits += nsh;
      }

      if (worstmva < _minpeakmva) {
        changed = true;
        std::swap(tc._strawHitIdxs[iworst],tc._strawHitIdxs.back());
        ComboHit const& ch = (*_chcol)[tc._strawHitIdxs.back()];
        unsigned nsh = ch.nStrawHits();
        nhits -= nsh;
        sumphi -= nsh * worstphi;
        sumt -= nsh *_ttcalc.comboHitTime(ch);
        tc._strawHitIdxs.pop_back();
      }
      enoughhits = nhits >= _minnhits;
      pphi  = sumphi/float(nhits);
      ptime = sumt/float(nhits);
    }

    if(enoughhits){
      // final pass: hard cut on dt
      std::vector<size_t> toremove;
      accumulator_set<float, stats<tag::weighted_mean >,unsigned > facc;
      accumulator_set<float, stats<tag::weighted_variance(lazy)>, unsigned > terr;
      accumulator_set<float, stats<tag::weighted_mean >,unsigned > racc;
      accumulator_set<float, stats<tag::weighted_mean >,unsigned > zacc;
      for(size_t ips=0;ips<tc._strawHitIdxs.size();++ips) {
        unsigned ish = tc._strawHitIdxs[ips];
        ComboHit const& ch = (*_chcol)[ish];
        unsigned nsh = ch.nStrawHits();
        float  time = _ttcalc.comboHitTime(ch);
        float  dt = time - ptime;
        float  phi = ch.phi();
        float rho = sqrtf(ch.pos().Perp2());
        Angles::deltaPhi(phi,pphi);
        if (fabs(dt) < _maxpeakdt) {
          terr(time,weight=nsh);
          facc(phi,weight=nsh);
          racc(rho,weight=nsh);
          zacc(ch.pos().z(),weight=nsh);
          tc._nsh += nsh;
        } else {
          toremove.push_back(ips);
        }
      }

      for (auto irm=toremove.rbegin();irm!=toremove.rend();++irm)
      {
        std::swap(tc._strawHitIdxs[*irm],tc._strawHitIdxs.back());
        tc._strawHitIdxs.pop_back();
      }

      tc._t0._t0 = extract_result<tag::weighted_mean>(terr);
      tc._t0._t0err = sqrtf(std::max(float(0.0),extract_result<tag::weighted_variance(lazy)>(terr))/extract_result<tag::count>(terr));
      pphi = extract_result<tag::weighted_mean>(facc);
      float prho = extract_result<tag::weighted_mean>(racc);
      float zpos = extract_result<tag::weighted_mean>(zacc);
      tc._pos = XYZVec(prho*cos(pphi),prho*sin(pphi),zpos);
    }
    //if (_debug > 0) std::cout<<"final time "<<tc._t0._t0<<std::endl;
  }

  void TimeClusterFinder::addCaloCluster(TimeCluster& tc)
  {
    auto bestcc = _cccol->end();
    for (auto icc = _cccol->begin();icc != _cccol->end(); ++icc)
    {
      if (icc->energyDep() > _ccmine) continue;
      float time = _ttcalc.caloClusterTime(*icc);
      if (fabs(tc._t0._t0-time) > _maxdt) continue;
      if (bestcc == _cccol->end() || icc->energyDep() > bestcc->energyDep()) bestcc = icc;
    }

    if (bestcc != _cccol->end())
    {
      size_t index = std::distance(_cccol->begin(),bestcc);
      tc._caloCluster = art::Ptr<CaloCluster>(_ccH,index);
    }
  }

  bool TimeClusterFinder::goodHit(const StrawHitFlag& flag) const
  {
    return flag.hasAllProperties(_hsel) && !flag.hasAnyProperty(_hbkg);
  }


}

using mu2e::TimeClusterFinder;
DEFINE_ART_MODULE(TimeClusterFinder);

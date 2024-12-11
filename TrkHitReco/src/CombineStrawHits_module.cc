//
// A module to combine straw hits in the same panel.  This improves resolution and
// reduces combinatoric for downstream modules.
//  Original Author: David Brown, LBNL
//
// Modified by B. Echenard (Caltech), assumes that the hits are ordered by panels
// Dave Brown confirmed this is the case

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "TMath.h"

#include <iostream>

namespace mu2e {

  class CombineStrawHits : public art::EDProducer {

    public:
      struct Config
      {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<int>              debug   { Name("debugLevel"),            Comment("Debug level"),0 };
        fhicl::Atom<art::InputTag>    CHC     { Name("ComboHitCollection"),    Comment("Input single-straw ComboHit collection") };
        fhicl::Atom<art::InputTag>    EWM     { Name("EventWindowMarker"),     Comment("EventWindowMarker")};
        fhicl::Atom<bool>             testflag{ Name("TestFlag"),              Comment("test input collection flag or not") };
        fhicl::Sequence<std::string>  shsel   { Name("StrawHitSelectionBits"), Comment("Input flag selection") };
        fhicl::Sequence<std::string>  shmask  { Name("StrawHitMask"),          Comment("Input flag anti-selection") };
        fhicl::Atom<bool>             useTOT  { Name("UseTOT"),                Comment("use tot to estimate drift time") };
        fhicl::Atom<float>            uerr    { Name("UError"),                Comment("intrinsic error along the straw (mm)") };
        fhicl::Atom<bool>             unsorted{ Name("Unsorted"),              Comment("Digi data unsorted by StrawId")};
        fhicl::Atom<int>              maxds   { Name("MaxDS"),                 Comment("maximum straw number difference") };
        fhicl::Atom<float>            maxdt   { Name("MaxDt"),                 Comment("maximum time separation between hits in ns") };
        fhicl::Atom<float>            maxwdchi{ Name("MaxWireDistDiffPull"),   Comment("maximum wire distance separation chi") };
        fhicl::Atom<bool>             filter  { Name("FilterHits"),            Comment("Filter hits (alternative is to just flag)") };
        fhicl::Atom<int>              minN    { Name("MinimumNHits"),          Comment("Minimum number of hits")};
        fhicl::Atom<int>              maxN    { Name("MaximumNHits"),          Comment("Maximum number of hits")};
        fhicl::Atom<float>            minE    { Name("MinimumEnergy"),         Comment("Minimum straw energy deposit (MeV)")};
        fhicl::Atom<float>            maxE    { Name("MaximumEnergy"),         Comment("Maximum straw energy deposit (MeV)")};
        fhicl::Atom<float>            minR    { Name("MinimumRadius"),         Comment("Minimum transverse radius (mm)")};
        fhicl::Atom<float>            maxR    { Name("MaximumRadius"),         Comment("Maximum transverse radius (mm)")};
        fhicl::Atom<float>            minT { Name("MinimumTime"),           Comment("Earliest StrawDigi time to process (nsec)")};
        fhicl::Atom<float>            maxT { Name("MaximumTime"),           Comment("Latest StrawDigi time to process (nsec)")};
        fhicl::OptionalAtom<float>    minTOff  { Name("MinimumTimeOffSpill"),           Comment("Earliest StrawDigi time to process (nsec)")};
        fhicl::OptionalAtom<float>    maxTOff  { Name("MaximumTimeOffSpill"),           Comment("Latest StrawDigi time to process (nsec)")};
        fhicl::Atom<bool>             checkWres{ Name("CheckWres"),            Comment("Check Wres for consistency") };
      };

      explicit CombineStrawHits(const art::EDProducer::Table<Config>& config);
      void produce( art::Event& e);

    private:
      void combine(EventWindowMarker const& ewm, ComboHitCollection const& chcOrig, ComboHitCollection& chcol);
      void combineHits(const ComboHitCollection& chcOrig, ComboHit& combohit);

      int           _debug;
      art::ProductToken<ComboHitCollection> const _chctoken;
      art::ProductToken<EventWindowMarker> const _ewmtoken;
      bool          _testflag;
      StrawHitFlag  _shsel;
      StrawHitFlag  _shmask;
      float         _maxdt;
      bool          _useTOT;
      float         _maxwdchi;
      float         _uvar;
      float         _minR2,_maxR2;
      float         _minT, _maxT, _minTOff, _maxTOff;
      bool _overrideminTOff, _overridemaxTOff;
      float         _minE, _maxE;
      int           _minN, _maxN;
      int           _maxds;
      bool          _unsorted;
      bool          _checkWres;
      bool          _filter;
      StrawIdMask   _mask;
  };

  CombineStrawHits::CombineStrawHits(const art::EDProducer::Table<Config>& config) :
    EDProducer{config},
    _debug(     config().debug()),
    _chctoken{consumes<ComboHitCollection>(config().CHC())},
    _ewmtoken{consumes<EventWindowMarker>(config().EWM())},
    _testflag(  config().testflag()),
    _shsel(     config().shsel()),
    _shmask(    config().shmask()),
    _maxdt(     config().maxdt()),
    _useTOT(    config().useTOT()),
    _maxwdchi(  config().maxwdchi()),
    _uvar(     config().uerr()*config().uerr()),
    _minR2(     config().minR()*config().minR()),
    _maxR2(     config().maxR()*config().maxR()),
    _minT(      config().minT()),
    _maxT(      config().maxT()),
    _overrideminTOff(      config().minTOff(_minTOff)),
    _overridemaxTOff(      config().maxTOff(_maxTOff)),
    _minE(      config().minE()),
    _maxE(      config().maxE()),
    _minN(      config().minN()),
    _maxN(      config().maxN()),
    _maxds(     config().maxds()),
    _unsorted( config().unsorted()),
    _checkWres( config().checkWres()),
    _filter(    config().filter()),
    _mask("uniquepanel")     // define the mask: ComboHits are made from straws in the same unique panel
    {
      produces<ComboHitCollection>();
    }

  void CombineStrawHits::produce(art::Event& event)
  {
    auto chcH = event.getValidHandle(_chctoken);
    const ComboHitCollection& chcOrig(*chcH);
    auto ewmH = event.getValidHandle(_ewmtoken);
    const EventWindowMarker& ewm(*ewmH);

    auto chcolNew = std::make_unique<ComboHitCollection>();
    chcolNew->reserve(chcOrig.size());
    chcolNew->setParent(chcH);

    if (_unsorted){
      // currently VST data is not sorted by panel number so we must sort manually
      // sort hits by panel
      ComboHitCollection chcolsort;
      chcolsort.reserve(chcOrig.size());
      chcolsort.setParent(chcOrig.parent());
      std::array<std::vector<uint16_t>,StrawId::_nupanels> panels;
      size_t nsh = chcOrig.size();
      for(uint16_t ish=0;ish<nsh;++ish){
        ComboHit const& ch = chcOrig[ish];
        // select hits based on flag
        panels[ch.strawId().uniquePanel()].push_back(ish);
      }
      for (uint16_t ipanel=0;ipanel<StrawId::_nupanels;ipanel++){
        for (uint16_t ish=0;ish<panels[ipanel].size();ish++){
          chcolsort.push_back(chcOrig.at(panels[ipanel][ish]));
        }
      }
      combine(ewm, chcolsort, *chcolNew);
    }else{
      combine(ewm, chcOrig, *chcolNew);
    }
    event.put(std::move(chcolNew));
  }


  void CombineStrawHits::combine(EventWindowMarker const& ewm, ComboHitCollection const& chcOrig, ComboHitCollection& chcol)
  {

    float minT = _minT;
    float maxT = _maxT;
    if (ewm.spillType() == EventWindowMarker::offspill && _overrideminTOff)
      minT = _minTOff;
    if (ewm.spillType() == EventWindowMarker::offspill && _overridemaxTOff)
      maxT = _maxTOff;

    std::vector<bool> isUsed(chcOrig.size(),false);
    for (size_t ich=0;ich<chcOrig.size();++ich) {
      if (isUsed[ich]) continue;
      isUsed[ich] = true;

      const ComboHit& hit1 = chcOrig[ich];
      if ( _testflag && hit1.flag().hasAnyProperty(StrawHitFlag::dead)) continue;
      if ( _testflag && (!hit1.flag().hasAllProperties(_shsel) || hit1.flag().hasAnyProperty(_shmask))) continue;
      ComboHit combohit;
      combohit.init(hit1,ich);
      int panel1 = hit1.strawId().uniquePanel();

      for (size_t jch=ich+1;jch<chcOrig.size();++jch) {
        if (isUsed[jch]) continue;
        const ComboHit& hit2 = chcOrig[jch];

        if (hit2.strawId().uniquePanel() != panel1) break;
        if (abs(hit2.strawId().straw()-hit1.strawId().straw())> _maxds ) continue; // hits are not sorted by straw number
        if ( _testflag && hit2.flag().hasAnyProperty(StrawHitFlag::dead)) continue;
        if ( _testflag && (!hit2.flag().hasAllProperties(_shsel) || hit2.flag().hasAnyProperty(_shmask)) ) continue;

        float dt = _useTOT ? fabs(hit1.correctedTime() - hit2.correctedTime()) : fabs(hit1.time() - hit2.time());
        if (dt > _maxdt) continue;

        float wderr = sqrtf(hit1.wireVar() + hit2.wireVar());
        float wdchi = fabs(hit1.wireDist() - hit2.wireDist())/wderr;
        if (wdchi > _maxwdchi) continue;

        bool ok = combohit.addIndex(jch);
        if (!ok){
          std::cout << "CombineStrawHits past limit" << std::endl;
        } else {
          isUsed[jch]= true;
        }
      }
      // clear the flag bits; they are reset later
      const static StrawHitFlag initialFlag("TimeDivision");
      combohit._flag = initialFlag;
      int nch = combohit.nCombo();
      if(nch  < _minN || nch > _maxN){
        if(_filter)break;
      } else
        combohit._flag.merge(StrawHitFlag::nhitsel);
      // actually combine the hits if necessar, and make the cuts
      if (nch > 1) combineHits(chcOrig, combohit);

      auto time = _useTOT ? combohit.correctedTime() : combohit.time();
      if (time < minT || time > maxT ){
        if(_filter)break;
      } else
        combohit._flag.merge(StrawHitFlag::timesel);

      auto energy = combohit.energyDep();
      if( energy > _maxE || energy < _minE ) {
        if(_filter)break;
      } else
        combohit._flag.merge(StrawHitFlag::energysel);

      auto r2 = combohit.pos().Perp2();
      if( r2 < _minR2 || r2 > _maxR2 ) {
        if(_filter)break;
      } else
        combohit._flag.merge(StrawHitFlag::radsel);
      combohit._mask = _mask;
      combohit._flag.merge(StrawHitFlag::panelcombo);
      chcol.push_back(std::move(combohit));
    }
  }


  void CombineStrawHits::combineHits(const ComboHitCollection& chcOrig, ComboHit& combohit)
  {
    // simple sums to speed up the trigger
    double eacc(0),ctacc(0),dtacc(0),twtsum(0),ptacc(0),wacc(0),wacc2(0),wwtsum(0);
    double etacc[2] = {0,0};
    XYZVectorF midpos;
    combohit._nsh = 0;
    if (_debug > 2) std::cout << "Combining " << combohit.nCombo() << " hits: ";

    for (size_t ich = 0; ich < combohit.nCombo(); ++ich)
    {
      size_t index = combohit.index(ich);
      if (_debug > 3)std::cout << index << ", ";
      if (index > chcOrig.size())
        throw cet::exception("RECO")<<"mu2e::CombineStrawHits: inconsistent index "<<std::endl;

      const ComboHit& ch = chcOrig[index];
      // simple average of basic quantities.  End times cannot take advantage of corrected time resolution
      combohit._nsh += ch.nStrawHits();
      eacc += ch.energyDep();
      etacc[StrawEnd::cal] += ch.endTime(StrawEnd::cal);
      etacc[StrawEnd::hv] += ch.endTime(StrawEnd::hv);
      dtacc += ch.driftTime();
      ptacc += ch.propTime();
      midpos += ch.centerPos();
      // weight time values by time error
      double twt = 1.0/(ch.timeVar());
      twtsum += twt;
      ctacc += ch.correctedTime()*twt;
      // weight position values by U (wire direction) error
      double wwt = 1.0/(ch.wireVar());
      wwtsum += wwt;
      wacc  += ch.wireDist()*wwt;
      wacc2 += ch.wireDist()*ch.wireDist()*wwt;
    }

    if(combohit.nStrawHits() != combohit.nCombo())
      throw cet::exception("RECO")<<"mu2e::CombineStrawHits: inconsistent count "<< std::endl;

    if(_debug > 2) std::cout << std::endl;

    double nsh = static_cast<double>(combohit.nStrawHits());
    combohit._edep       = eacc/nsh; // simple averages
    combohit._etime[StrawEnd::cal] = etacc[StrawEnd::cal]/nsh;
    combohit._etime[StrawEnd::hv] = etacc[StrawEnd::hv]/nsh;
    // these next should be removed if possible
    midpos /= nsh;
    combohit._dtime      = dtacc/nsh;
    combohit._ptime      = ptacc/nsh;
    combohit._time       = ctacc/twtsum;
    combohit._timevar    = 1.0/twtsum;
    combohit._wdist      = wacc/wwtsum;
    combohit._pos        = midpos + combohit._wdist*combohit.uDir();
    combohit._uvar       = 1.0/wwtsum + _uvar;
    combohit._vvar       = 1.0/twtsum;
    // compute U chisquared
    unsigned ndof = nsh - 1; // u direction only
    double chisq = wacc2 - wacc/wwtsum;
    combohit._qual       = TMath::Prob(chisq,ndof);
//-----------------------------------------------------------------------------
// for combohits made out of 2 and more straw hits:
// if combohit._ures is less than the sigma calculated on the individual hit positions,
// use the PDG prescription : scale the resolution to the sigma
//-----------------------------------------------------------------------------
    float wdist2 = wacc2/wwtsum;
    float nc2 = pow(ndof,2);
    float wvar = std::max(_uvar,(wdist2-combohit._wdist*combohit._wdist)/nc2);
    if (_checkWres)combohit._uvar = std::max(combohit._uvar,wvar);
  }
}

DEFINE_ART_MODULE(mu2e::CombineStrawHits)

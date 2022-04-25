//
// A module to combine straw hits in the same panel.  This improves resolution and
// reduces combinatoric for downstream modules.
//  Original Author: David Brown, LBNL
//
// Modified by B. Echenard (Caltech), assumes that the hits are ordered by panels
// Dave Brown confirmed this is the case

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>

namespace mu2e {

  class CombineStrawHits : public art::EDProducer {

    public:
      struct Config
      {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<int>              debug   { Name("debugLevel"),            Comment("Diagnosis level"),0 };
        fhicl::Atom<art::InputTag>    chTag   { Name("ComboHitCollection"),    Comment(" ") };
        fhicl::Atom<bool>             testflag{ Name("TestFlag"),              Comment("test flag or not") };
        fhicl::Atom<bool>             testrad { Name("TestRadius"),            Comment("test position radius") };
        fhicl::Sequence<std::string>  shsel   { Name("StrawHitSelectionBits"), Comment("flag selection") };
        fhicl::Sequence<std::string>  shmask  { Name("StrawHitMask"),          Comment("flag anti-selection") };
        fhicl::Atom<float>            maxdt   { Name("MaxDt"),                 Comment("maximum time separation between hits in ns") };
        fhicl::Atom<bool>             useTOT  { Name("UseTOT"),                Comment("use tot to estimate drift time") };
        fhicl::Atom<float>            maxwdchi{ Name("MaxWireDistDiffPull"),   Comment("maximum wire distance separation chi") };
        fhicl::Atom<float>            werr    { Name("WireError"),             Comment("intrinsic error on wire distance squared") };
        fhicl::Atom<float>            terr    { Name("TransError"),            Comment("intrinsic error transverse to wire (per straw)") };
        fhicl::Atom<float>            minR    { Name("MinimumRadius"),         Comment("Minimum transverse radius squared") };
        fhicl::Atom<float>            maxR    { Name("MaximumRadius"),         Comment("Maximum transverseradius squared") };
        fhicl::Atom<int>              maxds   { Name("MaxDS"),                 Comment("maximum straw number difference") };
        fhicl::Atom<bool>             isVSTdata{ Name("IsVSTData"),            Comment("Data from VST and needs sorting"), false };
      };

      explicit CombineStrawHits(const art::EDProducer::Table<Config>& config);
      void produce( art::Event& e);

    private:
      void combine(const ComboHitCollection* chcolOrig, ComboHitCollection& chcol);
      void combineHits(const ComboHitCollection* chcolOrig, ComboHit& combohit);

      int           _debug;
      art::InputTag _chTag;
      bool          _testflag;
      bool          _testrad;
      StrawHitFlag  _shsel;
      StrawHitFlag  _shmask;
      float         _maxdt;
      bool          _useTOT;
      float         _maxwdchi;
      float         _werr2;
      float         _terr;
      float         _minR2;
      float         _maxR2;
      int           _maxds;
      bool          _isVSTdata;
      StrawIdMask   _mask;
  };

  CombineStrawHits::CombineStrawHits(const art::EDProducer::Table<Config>& config) :
    EDProducer{config},
    _debug(     config().debug()),
    _chTag(     config().chTag()),
    _testflag(  config().testflag()),
    _testrad(   config().testrad()),
    _shsel(     config().shsel()),
    _shmask(    config().shmask()),
    _maxdt(     config().maxdt()),
    _useTOT(    config().useTOT()),
    _maxwdchi(  config().maxwdchi()),
    _werr2(     config().werr()*config().werr()),
    _terr(      config().terr()),
    _minR2(     config().minR()*config().minR()),
    _maxR2(     config().maxR()*config().maxR()),
    _maxds(     config().maxds()),
    _isVSTdata( config().isVSTdata()),
    _mask(      "uniquepanel")     // define the mask: ComboHits are made from straws in the same unique panel
    {
      consumes<ComboHitCollection>(_chTag);
      produces<ComboHitCollection>();
    }

  void CombineStrawHits::produce(art::Event& event)
  {
    // I have to get a Handle, not a ValidHandle, as a literal handle is needed to find the productID
    art::Handle<ComboHitCollection> chH;
    if(!event.getByLabel(_chTag, chH))
      throw cet::exception("RECO")<<"mu2e::CombineStrawHits: No ComboHit collection found for tag" <<  _chTag << std::endl;
    const ComboHitCollection* chcolOrig = chH.product();

    auto chcolNew = std::make_unique<ComboHitCollection>();
    chcolNew->reserve(chcolOrig->size());
    chcolNew->setParent(chH);

    if (_isVSTdata){
      // currently VST data is not sorted by panel number so we must sort manually
      // sort hits by panel
      ComboHitCollection* chcolsort = new ComboHitCollection();
      chcolsort->reserve(chcolOrig->size());
      chcolsort->setParent(chH);
      std::array<std::vector<uint16_t>,StrawId::_nupanels> panels;
      size_t nsh = chcolOrig->size();
      for(uint16_t ish=0;ish<nsh;++ish){
        ComboHit const& ch = (*chcolOrig)[ish];
        // select hits based on flag
        panels[ch.strawId().uniquePanel()].push_back(ish);
      }
      for (uint16_t ipanel=0;ipanel<StrawId::_nupanels;ipanel++){
        for (uint16_t ish=0;ish<panels[ipanel].size();ish++){
          chcolsort->push_back(chcolOrig->at(panels[ipanel][ish]));
        }
      }
      combine(chcolsort, *chcolNew);
    }else{
      combine(chcolOrig, *chcolNew);
    }
    event.put(std::move(chcolNew));
  }


  void CombineStrawHits::combine(const ComboHitCollection* chcolOrig, ComboHitCollection& chcol)
  {
    std::vector<bool> isUsed(chcolOrig->size(),false);
    for (size_t ich=0;ich<chcolOrig->size();++ich)
    {
      if (isUsed[ich]) continue;
      isUsed[ich] = true;

      const ComboHit& hit1 = (*chcolOrig)[ich];
      if ( _testflag && (!hit1.flag().hasAllProperties(_shsel) || hit1.flag().hasAnyProperty(_shmask)) ) continue;
      ComboHit combohit;
      combohit.init(hit1,ich);
      int panel1 = hit1.strawId().uniquePanel();

      for (size_t jch=ich+1;jch<chcolOrig->size();++jch)
      {
        if (isUsed[jch]) continue;
        const ComboHit& hit2 = (*chcolOrig)[jch];

        if (hit2.strawId().uniquePanel() != panel1) break;
        if (abs(hit2.strawId().straw()-hit1.strawId().straw())> _maxds ) continue; // hits are not sorted by straw number
        if ( _testflag && (!hit2.flag().hasAllProperties(_shsel) || hit2.flag().hasAnyProperty(_shmask)) ) continue;

        float dt = _useTOT ? fabs(hit1.correctedTime() - hit2.correctedTime()) : fabs(hit1.time() - hit2.time());
        if (dt > _maxdt) continue;

        float wderr = sqrtf(hit1.wireErr2() + hit2.wireErr2());
        float wdchi = fabs(hit1.wireDist() - hit2.wireDist())/wderr;
        if (wdchi > _maxwdchi) continue;

        bool ok = combohit.addIndex(jch);
        if (!ok) std::cout << "CombineStrawHits past limit" << std::endl;
        isUsed[jch]= true;
      }

      if (combohit.nCombo() > 1) combineHits(chcolOrig, combohit);

      float r2     = combohit.pos().Perp2();
      bool goodrad = r2 < _maxR2 && r2 > _minR2;
      if (goodrad) combohit._flag.merge(StrawHitFlag::radsel);
      if (!_testrad || goodrad) chcol.push_back(std::move(combohit));
    }
  }


  void CombineStrawHits::combineHits(const ComboHitCollection* chcolOrig, ComboHit& combohit)
  {
    combohit._mask = _mask;
    combohit._flag.merge(StrawHitFlag::panelcombo);

    float eacc(0),tacc(0),dtacc(0),ptacc(0),placc(0),werracc(0),wacc(0),wacc2(0),weights(0);
    XYZVectorF midpos;
    combohit._nsh = 0;
    if (_debug > 2) std::cout << "Combining " << combohit.nCombo() << " hits: ";

    for (size_t ich = 0; ich < combohit.nCombo(); ++ich)
    {
      size_t index = combohit.index(ich);
      if (_debug > 3)std::cout << index << ", ";
      if (index > chcolOrig->size())
        throw cet::exception("RECO")<<"mu2e::CombineStrawHits: inconsistent index "<<std::endl;

      const ComboHit& ch = (*chcolOrig)[index];
      combohit._flag.merge(ch.flag());
      float wt = 1.0/(ch.wireErr2());
      eacc += ch.energyDep();
      tacc += ch.time();// time is an unweighted average
      dtacc += ch.driftTime();
      ptacc += ch.propTime();
      placc += ch.pathLength();
      werracc += ch.wireRes();
      weights += wt;
      wacc  += ch.wireDist()*wt;
      wacc2 += ch.wireDist()*ch.wireDist()*wt;
      midpos += ch.centerPos(); // simple average for position
      combohit._nsh += ch.nStrawHits();
    }

    if(combohit.nStrawHits() < combohit.nCombo())
      throw cet::exception("RECO")<<"mu2e::CombineStrawHits: inconsistent count "<< std::endl;

    if(_debug > 2) std::cout << std::endl;

    midpos /= combohit._nsh;
    combohit._edep       = eacc/float(combohit.nCombo());
    combohit._time       = tacc/float(combohit.nCombo());
    combohit._dtime      = dtacc/float(combohit.nCombo());
    combohit._ptime      = ptacc/float(combohit.nCombo());
    combohit._pathlength = placc/float(combohit.nCombo());
    combohit._wdist      = wacc/weights;
    combohit._pos        = midpos + combohit._wdist*combohit._wdir;
    combohit._wres       = sqrt(1.0/weights + _werr2);
    combohit._tres       = _terr/sqrt(combohit._nsh); // error proportional to # of straws (roughly)
    //      float wvar           = sqrtf((wacc2/weights-wacc/weights*wacc/weights));//define quality as variance/average ratio
    combohit._qual       = 1.0; // quality isn't used and should be removed FIXME
  }
}

DEFINE_ART_MODULE(mu2e::CombineStrawHits)

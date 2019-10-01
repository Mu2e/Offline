//
// A module to combine straw hits in the same panel.  This improves resolution and
// reduces combinatorix for downstream modules.
//  Original Author: David Brown, LBNL
//

// Mu2e includes.
#include "RecoDataProducts/inc/ComboHit.hh"
// art includes.
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp>
using namespace boost::accumulators;
// C++ includes.
#include <iostream>
#include <float.h>

using namespace std;

namespace mu2e {

  class CombineStrawHits : public art::EDProducer {

  public:
    explicit CombineStrawHits(fhicl::ParameterSet const& pset);

    void produce( art::Event& e);

  private:
    // utility functions
    void combineHits(ComboHit& combohit);
    // configuration
    int _debug;
    // event object Tags
    art::InputTag _chTag;
    // input event collection cache
    const ComboHitCollection* _chcol;
    // Parameters
    bool _testflag; //test flag or not
    bool _testrad; // test position radius
    StrawHitFlag _shsel; // flag selection
    StrawHitFlag _shmask; // flag anti-selection
    float _maxdt; // maximum time separation between hits
    bool _useTOT; // use tot to estimate drift time
    float _maxwdchi; // maximum wire distance separation chi
    float _werr2; // intrinsic error on wire distance squared
    float _terr; // intrinsic error transverse to wire (per straw)
    float _minR2, _maxR2; // transverse radius (squared)
    int _maxds; // maximum straw number difference
    StrawIdMask _mask;
  };

  CombineStrawHits::CombineStrawHits(fhicl::ParameterSet const& pset) :
    art::EDProducer{pset},
    // Parameters
    _debug(pset.get<int>("debugLevel",0)),
    _chTag(pset.get<art::InputTag>("ComboHitCollection")),
    _testflag(pset.get<bool>("TestFlag")),
    _testrad(pset.get<bool>("TestRadius")),
    _shsel(pset.get<vector<string> >("StrawHitSelectionBits",vector<string>{"EnergySelection","TimeSelection"} )),
    _shmask(pset.get<vector<string> >("StrawHitMaskBits",vector<string>{} )),
    _maxdt(pset.get<float>("MaxDt",40.0)), // nsec 
    _useTOT(pset.get<bool>("UseTOT",false)), // use TOT corrected time
    _maxwdchi(pset.get<float>("MaxWireDistDiffPull",4.0)), //units of resolution sigma
    _terr(pset.get<float>("TransError",8.0)), //mm
    _maxds(pset.get<int>("MaxDS",3)) // how far away 2 straws can be, in 0-95 numbering (including layers!!)
  {
    consumes<ComboHitCollection>(_chTag);
    float werr = pset.get<float>("WireError",10.0); // mm
    _werr2 = werr*werr;
    produces<ComboHitCollection>();
    // define the mask: straws are in the same unique panel
    std::vector<StrawIdMask::field> fields{StrawIdMask::plane,StrawIdMask::panel};
    _mask = StrawIdMask(fields);
    float minR = pset.get<float>("minimumRadius",395); // mm
    _minR2 = minR*minR;
    float maxR = pset.get<float>("maximumRadius",650); // mm
    _maxR2 = maxR*maxR;
  }

  void CombineStrawHits::produce(art::Event& event)
  {
    // find event data.  Note I have to get a Handle, not a ValidHandle,
    // as a literal handle is needed to find the productID
    art::Handle<ComboHitCollection> chH;
    if(!event.getByLabel(_chTag, chH))
      throw cet::exception("RECO")<<"mu2e::CombineStrawHits: No ComboHit collection found for tag" <<  _chTag << endl;
    _chcol = chH.product();

    // create output
    auto chcol = std::make_unique<ComboHitCollection>();
    chcol->reserve(_chcol->size());
    // reference the parent in the new collection
    chcol->setParent(chH);

    // sort hits by panel
    std::array<std::vector<uint16_t>,StrawId::_nupanels> panels;
    size_t nsh = _chcol->size();
    for(uint16_t ish=0;ish<nsh;++ish){
      ComboHit const& ch = (*_chcol)[ish];
      // select hits based on flag
      if((!_testflag) || (ch.flag().hasAllProperties(_shsel) && (!ch.flag().hasAnyProperty(_shmask))) ){
        panels[ch.strawId().uniquePanel()].push_back(ish);
      }
    }
    // loop over panels
    for(auto const& phits : panels ) {
      // keep track of which hits are used as part of a combo hit
      std::vector<bool> used(phits.size(),false);
      // loop over hit pairs in this panel
      for(size_t ihit=0;ihit < phits.size(); ++ihit){
        if(!used[ihit]){
          used[ihit] = true;
          ComboHit const& hit1 = (*_chcol)[phits[ihit]];
          // create a combo hit for every hit; initialize it with this hit
          ComboHit combohit;
          combohit.init(hit1,phits[ihit]);
          // loop over other hits in this panel
          for(size_t jhit=ihit+1;jhit < phits.size(); ++jhit){
            if(!used[jhit]){
              ComboHit const& hit2 = (*_chcol)[phits[jhit]];
              // require straws be near each other
              int ds = abs( (int)hit1.strawId().straw()-(int)hit2.strawId().straw());
              if(ds > 0 && ds <= _maxds ){
                // require times be consistent
                float dt;
                if (_useTOT)
		  dt = fabs(hit1.correctedTime() - hit2.correctedTime());
                else
                  dt = fabs(hit1.time() - hit2.time());
                if(dt < _maxdt){
                  // compute the chi of the differnce in wire positions
                  float wderr = sqrtf(hit1.wireErr2() + hit2.wireErr2());
                  float wdchi = fabs(hit1.wireDist() - hit2.wireDist())/wderr;
                  if(wdchi < _maxwdchi){
                    // add a neural net selection here someday for Offline use  FIXME!
                    // these hits match: add the 2nd to the combo hit
                    bool ok = combohit.addIndex(phits[jhit]);
                    if(!ok)std::cout << "CombineStrawHits past limit" << std::endl;
                    used[jhit]= true;
                  } // consistent positions along wire
                }// consistent times
              }// straw proximity
            } // 2nd hit not used
          } // 2nd panel hit
          // compute floating point info for this combo hit and save it
          if(combohit.nCombo() > 1)combineHits(combohit);
          // radius test
          float r2 = combohit.pos().Perp2();
          bool goodrad = r2 < _maxR2 && r2 > _minR2;
          if(goodrad) combohit._flag.merge(StrawHitFlag::radsel);
          if(!_testrad || goodrad)
            chcol->push_back(std::move(combohit));
        } // 1st hit not used
      } // 1st panel hit
    } // panels
    // store data in the event
    event.put(std::move(chcol));
  }

  // compute the properties of this combined hit
  void CombineStrawHits::combineHits(ComboHit& combohit) {
    // if there's only 1 hit, take the info from the orginal collections
    // This is because the boost accumulators sometimes don't work for low stats
    // init from the 0th hit
    combohit._mask = _mask;
    // add the flag
    combohit._flag.merge(StrawHitFlag::panelcombo);
    accumulator_set<float, stats<tag::mean> > eacc;
    accumulator_set<float, stats<tag::mean> > tacc;
    accumulator_set<float, stats<tag::mean> > dtacc;
    accumulator_set<float, stats<tag::mean> > ptacc;
    accumulator_set<float, stats<tag::mean> > placc;
    accumulator_set<float, stats<tag::weighted_variance(lazy)>, float> wacc;
    accumulator_set<float, stats<tag::mean> > werracc;
    XYZVec midpos;
    combohit._nsh = 0;
    if(_debug > 2)std::cout << "Combining " << combohit.nCombo() << " hits: ";
    for(size_t ich = 0; ich < combohit.nCombo(); ++ich){
      // get back the original information
      size_t index = combohit.index(ich);
      if(_debug > 3)std::cout << index << ", ";
      if(index > _chcol->size())
        throw cet::exception("RECO")<<"mu2e::CombineStrawHits: inconsistent index "<< endl;
      ComboHit const& ch = (*_chcol)[index];
      combohit._flag.merge(ch.flag());
      eacc(ch.energyDep());
      tacc(ch.time());// time is an unweighted average
      dtacc(ch.driftTime());
      ptacc(ch.propTime());
      placc(ch.pathLength());
      float wt = 1.0/(ch.wireErr2());
      wacc(ch.wireDist(),weight=wt); // wire position is weighted
      werracc(ch.wireRes());
      midpos += ch.centerPos(); // simple average for position
      combohit._nsh += ch.nStrawHits();
    }
    if(combohit.nStrawHits() < combohit.nCombo())
      throw cet::exception("RECO")<<"mu2e::CombineStrawHits: inconsistent count "<< endl;
    if(_debug > 2)std::cout << std::endl;
    combohit._time = extract_result<tag::mean>(tacc);
    combohit._dtime = extract_result<tag::mean>(dtacc);
    combohit._ptime = extract_result<tag::mean>(ptacc);
    combohit._pathlength = extract_result<tag::mean>(placc);
    combohit._edep = extract_result<tag::mean>(eacc);
    combohit._wdist = extract_result<tag::weighted_mean>(wacc);
    midpos /= combohit._nsh;
    combohit._pos = midpos + combohit._wdist*combohit._wdir;
    combohit._wres = sqrt(1.0/extract_result<tag::sum_of_weights>(wacc) + _werr2);
    combohit._tres = _terr/sqrt(combohit._nsh); // error proportional to # of straws (roughly)
    // for now, define the quality as the ratio of the variance to the average
    float wvar = sqrtf(std::max(extract_result<tag::variance>(wacc),float(0.0)));
    combohit._qual = wvar/extract_result<tag::mean>(werracc);
  }
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::CombineStrawHits)

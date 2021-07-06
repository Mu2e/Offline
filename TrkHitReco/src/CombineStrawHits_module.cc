//
// A module to combine straw hits in the same panel.  This improves resolution and
// reduces combinatoric for downstream modules.
//  Original Author: David Brown, LBNL
//
// Modified by B. Echenard (Caltech), assumes that the hits are ordered by panels
// Dave Brown confirmed this is the case

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>

namespace mu2e {

  class CombineStrawHits : public art::EDProducer {

     public:
       explicit CombineStrawHits(fhicl::ParameterSet const& pset);
       void produce( art::Event& e);

     private:
       void combine(const ComboHitCollection* chcolOrig, ComboHitCollection& chcol);
       void combineHits(const ComboHitCollection* chcolOrig, ComboHit& combohit);

       int           _debug;
       art::InputTag _chTag;
       bool          _testflag;      // test flag or not
       bool          _testrad;       // test position radius
       StrawHitFlag  _shsel;         // flag selection
       StrawHitFlag  _shmask;        // flag anti-selection
       float         _maxdt;         // maximum time separation between hits
       bool          _useTOT;        // use tot to estimate drift time
       float         _maxwdchi;      // maximum wire distance separation chi
       float         _werr2;         // intrinsic error on wire distance squared
       float         _terr;          // intrinsic error transverse to wire (per straw)
       float         _minR2, _maxR2; // transverse radius (squared)
       int           _maxds;         // maximum straw number difference
       StrawIdMask   _mask;
  };



  CombineStrawHits::CombineStrawHits(fhicl::ParameterSet const& pset) :
    art::EDProducer{pset},
    _debug(     pset.get<int>("debugLevel",0)),
    _chTag(     pset.get<art::InputTag>("ComboHitCollection")),
    _testflag(  pset.get<bool>("TestFlag")),
    _testrad(   pset.get<bool>("TestRadius")),
    _shsel(     pset.get<std::vector<std::string>> ("StrawHitSelectionBits",std::vector<std::string>{"EnergySelection","TimeSelection"} )),
    _shmask(    pset.get<std::vector<std::string>>("StrawHitMask",std::vector<std::string>{} )),
    _maxdt(     pset.get<float>("MaxDt",40.0)), // nsec 
    _useTOT(    pset.get<bool>("UseTOT",false)), // use TOT corrected time
    _maxwdchi(  pset.get<float>("MaxWireDistDiffPull",4.0)), //units of resolution sigma
    _terr(      pset.get<float>("TransError",8.0)), //mm
    _maxds(     pset.get<int>("MaxDS",3)), // how far away 2 straws can be, in 0-95 numbering (including layers!!)
    _mask("uniquepanel")// define the mask: ComboHits are made from straws in the same unique panel
  {
      float werr = pset.get<float>("WireError",10.0); // mm
      _werr2 = werr*werr;
      float minR = pset.get<float>("minimumRadius",395); // mm
      _minR2 = minR*minR;
      float maxR = pset.get<float>("maximumRadius",650); // mm
      _maxR2 = maxR*maxR;
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

      combine(chcolOrig, *chcolNew);
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
              if (hit2.strawId().straw()-hit1.strawId().straw()> _maxds ) break;
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
      XYZVec midpos;
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
      float wvar           = sqrtf((wacc2/weights-wacc/weights*wacc/weights));//define quality as variance/average ratio
      combohit._qual       = wvar/(werracc/float(combohit.nCombo()));
  }
} 

DEFINE_ART_MODULE(mu2e::CombineStrawHits)

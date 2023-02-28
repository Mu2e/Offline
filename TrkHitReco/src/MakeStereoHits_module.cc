//
// A module to create simple stereo hits out of StrawHits. StrawHit selection is done by flagging in an upstream module
//
//
//  Original Author: David Brown, LBNL
//

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA_XYZ.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
using namespace boost::accumulators;

#include <iostream>
#include <float.h>
using namespace std;

namespace {

  struct StereoMVA
  {
    StereoMVA() : _pars(4,0.0),_dt(_pars[0]),_chisq(_pars[1]),_rho(_pars[2]),_ndof(_pars[3]){}

    std::vector <float> _pars;
    float& _dt;
    float& _chisq;
    float& _rho;
    float& _ndof;
  };
}

namespace mu2e {
  class MakeStereoHits : public art::EDProducer {
    public:
      struct Config {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<int>                 debug    { Name("debugLevel"), Comment("Debug level")};
        fhicl::Atom<art::InputTag>       chTag    { Name("ComboHitCollection"),   Comment("Input ComboHit collection") };
        fhicl::Sequence<std::string>     shsel    { Name("StrawHitSelectionBits"),    Comment("Mask for selecting hits") };
        fhicl::Sequence<std::string>     shrej    { Name("StrawHitRejectionBits"),   Comment("Mask for rejecting hits") };
        fhicl::Atom<float>               maxDt    { Name("MaxDt"),   Comment("Maximum time separation") };
        fhicl::Atom<bool>                useTOT   { Name("UseTOT"),   Comment("Use corrected time") };
        fhicl::Atom<float>               maxDPerp { Name("MaxDPerp"),   Comment("maximum perpendicular distance") };
        fhicl::Atom<float>               minwdot  { Name("MinWdot"),   Comment("minimum cos of angle between straws") };
        fhicl::Atom<float>               maxChisq { Name("MaxChisquared"),   Comment("position matching") };
        fhicl::Atom<float>               minMVA   { Name("MinMVA"),   Comment("MVA cut") };
        fhicl::Atom<float>               wfac     { Name("ZErrorFactor"),   Comment("error component due to z separation") };
        fhicl::Atom<float>               tfac     { Name("TErrorFactor"),   Comment("error component due to transverse separation") };
        fhicl::Atom<float>               minRho   { Name("MinRho"),   Comment("Minimum transverse radius of combination") };
        fhicl::Atom<float>               maxRho   { Name("MaxRho"),   Comment("Maximum transverse radius of combination") };
        fhicl::Atom<bool>                doMVA    { Name("DoMVA"),   Comment("Use MVA") };
        fhicl::Atom<unsigned>            maxfsep  { Name("MaxFaceSeparation"),   Comment("max separation between faces in a station") };
        fhicl::Atom<bool>                testflag { Name("TestFlag"),   Comment("") };
        fhicl::Atom<std::string>         smask    { Name("SelectionMask"),   Comment("define the mask to select hits") };
        fhicl::Table<MVATools::Config>   MVA      { Name("MVA"), Comment("MVA Configuration") };
      };

      explicit MakeStereoHits(const art::EDProducer::Table<Config>& config);
      void produce( art::Event& e);
      virtual void beginJob();
      virtual void beginRun(art::Run & run);
    private:
      typedef std::vector<uint16_t> ComboHits;

      int            _debug;
      art::InputTag _chTag;

      StrawHitFlag   _shsel;     // input flag selection
      StrawHitFlag   _shrej;     // input flag rejection
      float         _maxDt;      // maximum time separation between hits
      bool          _useTOT;     // use TOT to estimate drift time
      float         _maxDPerp;   // maximum transverse separation
      float         _minwdot;    // minimum dot product of straw directions
      float         _maxChisq;   // maximum chisquared to allow making stereo hits
      float         _minMVA;     // minimum MVA output
      float         _wfac;       // resolution factor along the wire
      float         _tfac;       // resolution transverse to the wire
      float  _minrho, _maxrho;   // transverse radius range
      bool          _doMVA;      // do MVA eval or simply use chi2 cut
      unsigned      _maxfsep;    // max face separation
      bool          _testflag;   // test the flag or not
      StrawIdMask   _smask;      // mask for combining hits

      MVATools _mvatool;
      StereoMVA _vmva;

      std::array<std::vector<StrawId>,StrawId::_nupanels > _panelOverlap;   // which panels overlap each other
      void genMap();
      void finalize(ComboHit& combohit,ComboHitCollection const& inchcol);
  };

  MakeStereoHits::MakeStereoHits(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    _debug(config().debug()),
    _chTag(config().chTag()),
    _shsel(config().shsel()),
    _shrej(config().shrej()),
    _maxDt(config().maxDt()),
    _useTOT(config().useTOT()),
    _maxDPerp(config().maxDPerp()),
    _minwdot(config().minwdot()),
    _maxChisq(config().maxChisq()),
    _minMVA(config().minMVA()),
    _wfac(config().wfac()),
    _tfac(config().tfac()),
    _minrho(config().minRho()),
    _maxrho(config().maxRho()),
    _doMVA(config().doMVA()),
    _maxfsep(config().maxfsep()),
    _testflag(config().testflag()),
    _smask(config().smask()),
    _mvatool(config().MVA())
    {}

  void MakeStereoHits::beginJob() {
    if(_doMVA){
      _mvatool.initMVA();
      if (_debug > 0){
        std::cout << "MakeStereoHits MVA parameters: " << std::endl;
       _mvatool.showMVA();
      }
    }
  }

  void MakeStereoHits::beginRun(art::Run & run) {
    genMap();
  }

  void MakeStereoHits::produce(art::Event& event) {
    // I have to get a Handle, not a ValidHandle, as a literal handle is needed to find the productID
    art::Handle<ComboHitCollection> chH;
    if(!event.getByLabel(_chTag, chH))
      throw cet::exception("RECO")<<"mu2e::MakeStereoHits: No ComboHit collection found for tag" <<  _chTag << std::endl;
    auto const& inchcol = *chH.product();
    auto chcol = std::make_unique<ComboHitCollection>();
    chcol->reserve(inchcol.size());
    chcol->setParent(chH);
    // sort hits by unique panel.  This should be built in by construction upstream FIXME!!
    std::array<std::vector<uint16_t>,StrawId::_nupanels> phits;
    size_t nch = inchcol.size();
    if(_debug > 2)cout << "MakeStereoHits found " << nch << " Input hits" << endl;
    std::vector<bool> used(nch,false);
    for(uint16_t ihit=0;ihit<nch;++ihit){
      ComboHit const& ch = inchcol[ihit];
      // select hits based on flag
      if( (!_testflag) ||( ch.flag().hasAllProperties(_shsel) && (!ch.flag().hasAnyProperty(_shrej))) ){
        phits[ch.strawId().uniquePanel()].push_back(ihit);
      }
    }
    if(_debug > 3){
      for (unsigned ipan=0; ipan < StrawId::_nupanels; ++ipan) {
        if(phits[ipan].size() > 0 ){
          cout << "Panel " << ipan << " has " << phits[ipan].size() << " hits "<< endl;
        }
      }
    }
    //  Loop over all hits.  Every one must appear somewhere in the output
    for (size_t ihit=0;ihit<nch;++ihit) {
      if(used[ihit])continue;
      used[ihit] = true;
      // create an output combo hit for every hit; initialize it with this hit
      ComboHit const& ch1 = inchcol[ihit];
      ComboHit combohit;
      combohit.init(ch1,ihit);
      // zero values that accumulate in pairs
      combohit._qual = 0.0;
      combohit._pos = XYZVectorF(0.0,0.0,0.0);
      // loop over the panels which overlap this hit's panel
      for (auto sid : _panelOverlap[ch1.strawId().uniquePanel()]) {
        // loop over hits in the overlapping panel
        for (auto jhit : phits[sid.uniquePanel()]) {
          const ComboHit& ch2 = inchcol[jhit];
          if(_debug > 4) cout << " comparing hits " << ch1.strawId().uniquePanel() << " and " << ch2.strawId().uniquePanel();
          if (!used[jhit] ){
            float dt;
            if (_useTOT)
              dt = fabs(ch1.correctedTime()-ch2.correctedTime());
            else
              dt = fabs(ch1.time()-ch2.time());
            if(_debug > 4) cout << " dt = " << dt;
            if (dt < _maxDt){
              float wdot = ch1.wdir().Dot(ch2.wdir());
              XYZVectorF dp = ch1.pos()-ch2.pos();
              float dperp = sqrt(dp.perp2());
              // negative crosings are in opposite quadrants and longitudinal separation isn't too big
              if(_debug > 4) cout << " wdot = " << wdot << " dperp = " << dperp;
              if (wdot > _minwdot && dperp < _maxDPerp ) {
                // solve for the POCA.
                TwoLinePCA_XYZ pca(ch1.pos(),ch1.wdir(),ch2.pos(),ch2.wdir());
                if(pca.closeToParallel()){
                  cet::exception("RECO")<<"mu2e::StereoHit: parallel wires" << std::endl;
                }
                // check the intersection is inside the tracker active volume
                float rho = pca.point1().Rho();
                if(_debug > 4) cout << " rho = " << rho;
                if(rho < _maxrho && rho > _minrho ){
                  // compute chisquared; include error for particle angle
                  // should be a cumulative linear regression FIXME!
                  float terr = _tfac*fabs(ch1.pos().z()-ch2.pos().z());
                  float terr2 = terr*terr;
                  float dw1 = pca.s1();
                  float dw2 = pca.s2();
                  float chisq = dw1*dw1/(ch1.wireVar()+terr2) + dw2*dw2/(ch2.wireVar()+terr2);
                  if(_debug > 3) cout << " chisq = " << chisq;
                  if (chisq < _maxChisq){
                    if(_debug > 3) cout << " added ";
                    // if we get to here, try to add the hit
                    // accumulate the chisquared
                    if(combohit.addIndex(jhit)) {
                      // average z
                      combohit._qual += chisq;
                      combohit._pos += XYZVectorF(pca.point1().x(),pca.point1().y(),0.5*(pca.point1().z()+pca.point2().z()));
                    } else
                      throw cet::exception("RECO")<< "MakeStereoHits can't add hit" << std::endl;
                    used[jhit] = true;
                  }
                }
              }
            }
          }
          if(_debug > 3) cout << endl;
        }
      }
      finalize(combohit,inchcol);
      chcol->push_back(std::move(combohit));
    }
    event.put(std::move(chcol));
  }

  void MakeStereoHits::finalize(ComboHit& combohit,ComboHitCollection const& inchcol) {
    combohit._mask = _smask;
    if(combohit.nCombo() > 1){
      combohit._flag.merge(StrawHitFlag::stereo);
      combohit._flag.merge(StrawHitFlag::radsel);
      accumulator_set<float, stats<tag::weighted_mean>, unsigned > eacc;
      accumulator_set<float, stats<tag::weighted_mean>, unsigned > tacc;
      accumulator_set<float, stats<tag::weighted_mean>, unsigned > dtacc;
      accumulator_set<float, stats<tag::weighted_mean>, unsigned > placc;
      accumulator_set<float, stats<tag::min > > zmin;
      accumulator_set<float, stats<tag::max > > zmax;
      combohit._nsh = 0;
      for(size_t ich = 0; ich < combohit.nCombo(); ++ich){
        size_t index = combohit.index(ich);
        ComboHit const& ch = inchcol[index];
        combohit._flag.merge(ch.flag());
        eacc(ch.dEdx(),weight=ch.pathLength()*ch.nStrawHits());
        tacc(ch.time(),weight=ch.nStrawHits());
        dtacc(ch.driftTime(),weight=ch.nStrawHits());
        placc(ch.pathLength(),weight=ch.nStrawHits());
        zmin(ch.pos().z());
        zmax(ch.pos().z());
        combohit._nsh += ch.nStrawHits();
      }
      float maxz = extract_result<tag::max>(zmax);
      float minz = extract_result<tag::min>(zmin);
      combohit._time = extract_result<tag::weighted_mean>(tacc);
      combohit._dtime = extract_result<tag::weighted_mean>(dtacc);
      combohit._pathlength = extract_result<tag::weighted_mean>(placc);
      combohit._dedx = extract_result<tag::weighted_mean>(eacc);
      float dz = (maxz-minz);
      combohit._wdist = dz;
      combohit._tres = dz*_tfac;
      combohit._wres = dz*_wfac;
      combohit._qual /= (combohit.nCombo()-1);// normalize by # of pairs
      combohit._pos /= (combohit.nCombo()-1);
      combohit._wdir = XYZVectorF(0.0,0.0,1.0);
    } else {
      size_t index = combohit.index(0);
      combohit._pos = inchcol[index].pos();// put back original position
    }
  }

  // generate the overlap map
  void MakeStereoHits::genMap() {
    static bool init(false);
    if(!init){
      init = true;
      // initialize
      const Tracker& tt(*GeomHandle<Tracker>());
      // establihit the extent of a panel using the longest straw (0)
      Straw const& straw = tt.getStraw(StrawId(0,0,0));
      float phi0 = (straw.getMidPoint()-straw.halfLength()*straw.getDirection()).phi();
      float phi1 = (straw.getMidPoint()+straw.halfLength()*straw.getDirection()).phi();
      float lophi = std::min(phi0,phi1);
      float hiphi = std::max(phi0,phi1);
      float phiwidth = hiphi-lophi;
      if (phiwidth>M_PI) phiwidth = 2*M_PI-phiwidth;
      if(_debug > 1)std::cout << "Panel Phi width = " << phiwidth << std::endl;
      // loop over all unique panels
      for(size_t ipla = 0;ipla < StrawId::_nplanes; ++ipla) {
        for(int ipan=0;ipan<StrawId::_npanels;++ipan){
          StrawId sid(ipla,ipan,0);
          uint16_t upan = sid.uniquePanel();
          Straw const& straw = tt.getStraw(StrawId(ipla,ipan,0));
          float phi = straw.getMidPoint().phi();
          if(_debug > 1)std::cout << "Plane " << ipla << " Panel " << ipan << " phi = " << phi << " z = " << straw.getMidPoint().z() << endl;
          // loop over nearby panels and check for an overlap
          size_t minpla = (size_t)std::max(0,(int)ipla-1);
          size_t maxpla = (size_t)std::min(StrawId::_nplanes-1,(int)ipla+1);
          for(size_t jpla = minpla; jpla <= maxpla;++jpla){
            for(int jpan=0;jpan<StrawId::_npanels;++jpan){
              StrawId osid(jpla,jpan,0);
              Straw const& ostraw = tt.getStraw(StrawId(jpla,jpan,0));
              if(_smask.equal(osid,sid) && osid.uniqueFace() != sid.uniqueFace()
                  && (unsigned)abs(osid.uniqueFace() - sid.uniqueFace()) <= _maxfsep ) {
                float dphi = fabs(phi - ostraw.getMidPoint().phi());
                if (dphi > M_PI) dphi = 2*M_PI-dphi;
                if (dphi < phiwidth) _panelOverlap[upan].push_back(osid);
              }
            }
          }
        }
      }
      if (_debug >1) {
        for(uint16_t ipan = 0; ipan < StrawId::_nupanels; ++ipan) {
          std::cout << "Unique Panel " << ipan << " Overlaps with the panels: ";
          for(auto sid : _panelOverlap[ipan])
            std::cout << sid.uniquePanel() << ", ";
          std::cout << std::endl;
        }
      }
    }
  }

}

using mu2e::MakeStereoHits;
DEFINE_ART_MODULE(MakeStereoHits)


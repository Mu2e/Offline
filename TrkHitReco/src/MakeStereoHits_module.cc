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

#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/TrkHitReco/inc/CombineStereoPoints.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
// boost
//#include <boost/accumulators/accumulators.hpp>
//#include <boost/accumulators/statistics/mean.hpp>
//#include <boost/accumulators/statistics/stats.hpp>
//#include <boost/accumulators/statistics/weighted_variance.hpp>
//#include <boost/accumulators/statistics/max.hpp>
//#include <boost/accumulators/statistics/min.hpp>
//using namespace boost::accumulators;
#include "TMath.h"

#include <iostream>
#include <limits>
#include <cfloat>
#include <list>

namespace mu2e {
  class MakeStereoHits : public art::EDProducer {
    public:
      struct Config {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<int>                 debug    { Name("DebugLevel"), Comment("Debug level")};
        fhicl::Atom<art::InputTag>       CHC      { Name("ComboHitCollection"),   Comment("Input ComboHit collection") };
        fhicl::Sequence<std::string>     shsel    { Name("StrawHitSelectionBits"),    Comment("Mask for selecting hits") };
        fhicl::Sequence<std::string>     shrej    { Name("StrawHitRejectionBits"),   Comment("Mask for rejecting hits") };
        fhicl::Atom<float>               maxDt    { Name("MaxDt"),   Comment("Maximum time separation (ns)") };
        fhicl::Atom<float>               maxDPerp { Name("MaxDPerp"),   Comment("maximum perpendicular distance (mm)") };
        fhicl::Atom<float>               maxwdot  { Name("MaxWdot"),   Comment("maximum cos of angle between straws") };
        fhicl::Atom<float>               maxChisq { Name("MaxChisquared"),   Comment("position matching") };
        fhicl::Atom<float>               uvres    { Name("UVRes"),   Comment("Resolution in U,V (X,Y) plane") };
        fhicl::Atom<float>               minRho   { Name("MinRho"),   Comment("Minimum transverse radius of combination (mm)") };
        fhicl::Atom<float>               maxRho   { Name("MaxRho"),   Comment("Maximum transverse radius of combination (mm)") };
        fhicl::Atom<float>               minE     { Name("MinimumEnergy"),         Comment("Minimum straw energy deposit (MeV)")};
        fhicl::Atom<float>               maxE     { Name("MaximumEnergy"),         Comment("Maximum straw energy deposit (MeV)")};
        fhicl::Atom<unsigned>            maxfsep  { Name("MaxFaceSeparation"),   Comment("max separation between faces") };
        fhicl::Atom<float>               maxDz    { Name("MaxDz"),   Comment("max Z separation between panels (mm)") };
        fhicl::Atom<bool>                testflag { Name("TestFlag"),   Comment("Test input hit flags") };
        fhicl::Atom<bool>                filter   { Name("FilterHits"),            Comment("Filter hits (alternative is to just flag)") };
        fhicl::Atom<bool>                sline    { Name("StereoLine"),            Comment("Fit stereohit ComboHit daughters to a line") };
        fhicl::Atom<unsigned>            slinendof{ Name("StereoLineNDOF"),            Comment("Minimum NDOF for stereo line fit to be used in output ComboHits") };
        fhicl::Atom<std::string>         smask    { Name("SelectionMask"),   Comment("define the mask to select hits") };
      };

      explicit MakeStereoHits(const art::EDProducer::Table<Config>& config);
      void produce( art::Event& e);
      virtual void beginJob();
      virtual void beginRun(art::Run & run);
    private:
      typedef std::vector<uint16_t> ComboHits;

      int            _debug;
      art::ProductToken<ComboHitCollection> const _chctoken;

      StrawHitFlag   _shsel;     // input flag selection
      StrawHitFlag   _shrej;     // input flag rejection
      float         _maxDt;      // maximum time separation between hits
      float         _maxDPerp;   // maximum transverse separation
      float         _maxwdot;    // minimum dot product of straw directions
      float         _maxChisq;   // maximum chisquared to allow making stereo hits
      float         _uvvar;      // intrinsic variance in UV plane
      float         _minrho, _maxrho;   // transverse radius range
      float         _minE, _maxE; // edep cuts
      unsigned      _maxfsep;    // max face separation
      float         _maxDz;      // max z sepration
      bool          _testflag;   // test the flag or not
      bool          _filter;
      bool          _sline;      // fit to a line
      unsigned      _slinendof;  // minimum NDOF to use the sline fit when producing output ComboHits
      StrawIdMask   _smask;      // mask for combining hits

      std::array<std::vector<StrawId>,StrawId::_nupanels > _panelOverlap;   // which panels overlap each other
      void genMap();
      void fillComboHit(ComboHit& ch, CombineStereoPoints const& cpts, ComboHitCollection const& inchcol) const;
  };

  MakeStereoHits::MakeStereoHits(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    _debug(config().debug()),
    _chctoken{consumes<ComboHitCollection>(config().CHC())},
    _shsel(config().shsel()),
    _shrej(config().shrej()),
    _maxDt(config().maxDt()),
    _maxDPerp(config().maxDPerp()),
    _maxwdot(config().maxwdot()),
    _maxChisq(config().maxChisq()),
    _uvvar(config().uvres()*config().uvres()),
    _minrho(config().minRho()),
    _maxrho(config().maxRho()),
    _minE(config().minE()),
    _maxE(config().maxE()),
    _maxfsep(config().maxfsep()),
    _maxDz(config().maxDz()),
    _testflag(config().testflag()),
    _filter(config().filter()),
    _sline(config().sline()),
    _slinendof(config().slinendof()),
    _smask(config().smask())
    {
      produces<ComboHitCollection>();
    }

  void MakeStereoHits::beginJob() {
  }

  void MakeStereoHits::beginRun(art::Run & run) {
    genMap();
  }

  void MakeStereoHits::produce(art::Event& event) {
    auto chcH = event.getValidHandle(_chctoken);
    const ComboHitCollection& inchcol(*chcH);
    auto chcol = std::make_unique<ComboHitCollection>();
    chcol->reserve(inchcol.size());
    chcol->setParent(chcH);
    // sort hits by unique panel.  This should be built in by construction upstream TODO!!
    std::array<std::vector<uint16_t>,StrawId::_nupanels> phits;
    size_t nch = inchcol.size();
    if(_debug > 2)std::cout << "MakeStereoHits found " << nch << " Input hits" << std::endl;
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
          std::cout << "Panel " << ipan << " has " << phits[ipan].size() << " hits "<< std::endl;
        }
      }
    }
    //  Loop over all hits.  Every one must appear somewhere in the output
    for (size_t ihit=0;ihit<nch;++ihit) {
      if(used[ihit])continue;
      used[ihit] = true;
      ComboHit const& ch1 = inchcol[ihit];
      // initialize new hit
      ComboHit combohit;
      combohit._sid = _smask.maskStrawId(ch1.strawId());
      // create a combiner seeded on this hit
      CombineStereoPoints cpts(_uvvar);
      cpts.addPoint(StereoPoint(ch1.pos(),ch1.uDir(),ch1.uVar(),ch1.vVar()),ihit);
      if( (!_testflag) ||( ch1.flag().hasAllProperties(_shsel) && (!ch1.flag().hasAnyProperty(_shrej))) ){
        // loop over the panels which overlap this hit's panel
        for (auto sid : _panelOverlap[ch1.strawId().uniquePanel()]) {
          // loop over hits in the overlapping panel
          for (auto jhit : phits[sid.uniquePanel()]) {
            const ComboHit& ch2 = inchcol[jhit];
            if (!used[jhit] && cpts.nPoints() < ComboHit::MaxNCombo  && ( (!_testflag) ||( ch2.flag().hasAllProperties(_shsel) && (!ch2.flag().hasAnyProperty(_shrej)))) ){
              if(_debug > 3) std::cout << " comparing hits in panels " << ch1.strawId().uniquePanel() << " and " << ch2.strawId().uniquePanel() << std::endl;
              float dt;
              dt = fabs(ch1.correctedTime()-ch2.correctedTime());
              if(_debug > 4) std::cout << " dt = " << dt;
              if (dt < _maxDt){
                XYZVectorF dp = ch1.pos()-ch2.pos();
                auto dperp = dp.Rho();
                // negative crosings are in opposite quadrants and longitudinal separation isn't too big
                if(_debug > 4) std::cout << " dperp = " << dperp;
                if (dperp < _maxDPerp ) {
                  // compute chisq WRT this point
                  StereoPoint twodpt(ch2.pos(),ch2.uDir(),ch2.uVar(),ch2.vVar());
                  double dchi2 = cpts.dChi2(twodpt);
                  if(_debug > 3) std::cout << " dchisq = " << dchi2;
                  if(dchi2 < _maxChisq){
                    // provisionally add this hit
                    cpts.addPoint(twodpt,jhit);
                    auto rho = cpts.point().point().pos2().R();
                    if(_debug > 4) std::cout << " rho = " << rho << std::endl;
                    if(rho < _maxrho && rho > _minrho ){
                      if(_debug > 3) std::cout << " added index";
                      used[jhit] = true;
                    } else {
                      // remove the point
                      cpts.removePoint(jhit);
                    }
                  }
                }
              }
            }
            if(_debug > 3) std::cout << std::endl;
          }
        }
      }
      if(cpts.nPoints() >1){
        // create the combined hit
        fillComboHit(combohit,cpts,inchcol);
      } else {
        // only 1 combo hit: either failed pre-selection or no matching hits.  Just reference the initial hit
        combohit.init(ch1,ihit);
      }
      // test the final hit
      auto rho = combohit.pos().Rho();
      if(rho < _maxrho && rho > _minrho ){
        combohit._flag.merge(StrawHitFlag::radsel);
      } else {
        combohit._flag.clear(StrawHitFlag::radsel);
      }
      auto energy = combohit.energyDep();
      if( energy < _maxE && energy > _minE ) {
        combohit._flag.merge(StrawHitFlag::energysel);
      } else {
        combohit._flag.clear(StrawHitFlag::energysel);
      }
      // update the level
      combohit._mask = _smask;
      // final test
      if( (!_filter) || ( combohit.flag().hasAllProperties(_shsel) &&
            (!combohit.flag().hasAnyProperty(_shrej))) ) chcol->push_back(combohit);
    }
    event.put(std::move(chcol));
  }

  void MakeStereoHits::fillComboHit(ComboHit& combohit, CombineStereoPoints const& cpts, ComboHitCollection const& inchcol) const {
    if(_debug > 1){
      std::cout << "Combining " << cpts.nPoints() << " hits" << std::endl;
    }
    bool sth = _sline;
    if(sth){
    // solve for the line
      StereoLine sline;
      sth = cpts.stereoLine(sline) && sline.ndof_ >= _slinendof;
      // if the fit succeeded and has enough NDOF use it to fill the ComboHit
      if(sth){
        combohit._pos = sline.pos(sline.z0());
        // create a 2-D point from the upper component of this
        TwoDPoint spt(sline.pars().Sub<TwoDPoint::SVEC>(0), sline.cov().Sub<TwoDPoint::SMAT>(0,0));
        // use this covariance to define the u direction
        combohit._udir = spt.udir();
        // project the covariances
        combohit._uvar = spt.uvar();
        combohit._vvar = spt.vvar();
        // convert slopes from the fit into a direction and direction errors
        combohit._hdir = sline.dir();
        TwoDPoint::SMAT dmat = sline.cov().Sub<TwoDPoint::SMAT>(StereoLine::dxdz,StereoLine::dxdz);
        // variance is on cos(theta), so a flat prior
        TwoDPoint::SVEC R = sline.pars().Sub<TwoDPoint::SVEC>(StereoLine::dxdz);
        auto rmag = R(0)*R(0) + R(1)*R(1);
        TwoDPoint::SVEC dCostdR =  R * std::pow(1.0 + rmag,-1.5);
        TwoDPoint::SVEC dPhidR = TwoDPoint::SVEC(-R(1),R(0)) * (1.0/rmag);
        combohit._hcostvar = ROOT::Math::Similarity(dCostdR,dmat);
        combohit._hphivar = ROOT::Math::Similarity(dPhidR,dmat);
        // fit quality
        combohit._qual = TMath::Prob(sline.chisq(),sline.ndof());
        // flag this hit as having stereo line information
        combohit._flag.merge(StrawHitFlag::sline);
      }
    }
    // otherwise, fall back to the 2D projection
    // fill position and variance from combined info
    if(!sth){
      auto const& pt = cpts.point();
      combohit._pos = pt.pos3();
      combohit._udir = pt.udir();
      combohit._uvar = pt.uvar();
      combohit._vvar = pt.vvar();
      combohit._qual = cpts.consistency();
    }

    // fill the remaining variables
    double twtsum(0), zmin(std::numeric_limits<float>::max()), zmax(std::numeric_limits<float>::lowest());
    for(auto iwt = cpts.weights().begin(); iwt != cpts.weights().end(); ++iwt){
      auto ihit = iwt->first;
      if(combohit.addIndex(ihit)) {
        auto const& ch = inchcol[ihit];
        unsigned nsh = ch.nStrawHits();
        combohit._nsh += nsh;
        combohit._flag.merge(ch.flag());
        double z = ch._pos.Z();
        zmin = std::min(zmin,z);
        zmax = std::max(zmax,z);
        combohit._flag.merge(StrawHitFlag::stereo);
        combohit._flag.clear(StrawHitFlag::panelcombo);
        combohit._edep += ch.energyDep()*nsh;
        // the following have unclear meaning for stereo hits, but we fill them anyways
        double twt = 1.0/ch.timeVar();
        for(size_t iend=0;iend<2;++iend){
          combohit._etime[iend] += ch._etime[iend]*twt;
          combohit._tot[iend] += ch._tot[iend]*twt;
        }
        twtsum += twt;
        combohit._time += ch.correctedTime()*twt;
        // to remove
        combohit._dtime += ch._dtime*twt;
        combohit._ptime += ch._ptime*twt;
      } else {
        std::cout << "MakeStereoHits past limit" << std::endl;
      }
    }
    // define w error from range
    static const double inv12 = 1.0/12.0;
    double dz = zmax-zmin;
    combohit._wvar = inv12*dz*dz;
    // simple average for edep
    combohit._edep /= combohit._nsh;
    // average time
    combohit._time /= twtsum;
    combohit._dtime /= twtsum;
    combohit._ptime  /= twtsum;
    combohit._timevar = 1.0/twtsum;
    for(size_t iend=0;iend<2;++iend){
      combohit._etime[iend] /=twtsum;
      combohit._tot[iend] /= twtsum;
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
          if(_debug > 1)std::cout << "Plane " << ipla << " Panel " << ipan << " phi = " << phi << " z = " << straw.getMidPoint().z() << std::endl;
          // loop over nearby panels and check for an overlap
          for(size_t jpla = 0;jpla < StrawId::_nplanes; ++jpla) {
            for(int jpan=0;jpan<StrawId::_npanels;++jpan){
              StrawId osid(jpla,jpan,0);
              if(osid != sid && _smask.equal(osid,sid) && (unsigned)abs(osid.uniqueFace() - sid.uniqueFace()) <= _maxfsep ) {
                Straw const& ostraw = tt.getStraw(StrawId(jpla,jpan,0));
                float dphi = fabs(fmod(phi - ostraw.getMidPoint().phi(),2*M_PI));
                if (dphi > M_PI) dphi = 2*M_PI-dphi;
                if(_debug > 1)std::cout << "Test Plane " << jpla << " Panel " << jpan << " dphi " << dphi << std::endl;
                if (dphi < phiwidth){
                  // insure the straws aren't parallel and are close enough in Z
                  double wdot = fabs(straw.direction().dot(ostraw.direction()));
                  double dz = fabs((straw.origin()-ostraw.origin()).z());
                  if(_debug > 1)std::cout << "Dz " << dz << " wdot " << wdot << std::endl;
                  if(wdot < _maxwdot && dz < _maxDz ){
                    if(_debug > 1)std::cout << "Added overlapping panel " << std::endl;
                    _panelOverlap[upan].push_back(osid);
                  }
                }
              }
            }
          }
        }
      }
      if (_debug >0) {
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


//
// A module to create simple stereo hits out of StrawHits. StrawHit selection is done by flagging in an upstream module
//
// $Id: MakeStereoHits_module.cc,v 1.23 2014/09/18 08:42:47 brownd Exp $
// $Author: brownd $
// $Date: 2014/09/18 08:42:47 $
// 
//  Original Author: David Brown, LBNL
//  

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeneralUtilities/inc/TwoLinePCA_XYZ.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "Mu2eUtilities/inc/MVATools.hh"
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

    std::vector <Double_t> _pars;
    Double_t& _dt; 
    Double_t& _chisq; 
    Double_t& _rho;  
    Double_t& _ndof; 
  };
}

namespace mu2e {
  class MakeStereoHits : public art::EDProducer {
    public:
      enum matchtype{panel=0,station};
      enum mode{filter=0,flag};
      explicit MakeStereoHits(fhicl::ParameterSet const& pset);
      void produce( art::Event& e);
      virtual void beginJob();
      virtual void beginRun(art::Run & run);
    private:
      typedef std::vector<uint16_t> ComboHits;

      int            _debug;
      art::InputTag  _shTag;
      art::InputTag  _chTag;
      art::InputTag  _shfTag;
      ComboHitCollection const* _chcol;

      StrawHitFlag   _shsel;      // flag selection
      StrawHitFlag   _shmask;     // flag anti-selection 
      double         _maxDt;      // maximum time separation between hits
      double         _maxDE;      // maximum deposited energy deference: this excludes inconsistent hits
      double         _maxDZ;      // maximum longitudinal separation
      double         _maxDPerp;   // maximum transverse separation
      double         _minDdot;    // minimum dot product of straw directions
      double         _rmax;      // maximum radius for stereo to be in active volume
      double         _maxChi;     // maximum # of TimeDivision sigmas past the active edge of a wire to allow making stereo hits
      double         _maxChisq;   // maximum chisquared to allow making stereo hits
      double         _minMVA;     // minimum MVA output
      double         _wfac;       // resolution factor along the wire
      double         _tfac;       // resolution transverse to the wire
      bool           _doMVA;      // do MVA eval or simply use chi2 cut
      matchtype	    _match;	  // match either to panel or station
      mode	    _mode; // what mode to operat
      StrawIdMask _mask; 

      MVATools _mvatool;
      StereoMVA _vmva; 

      std::array<std::vector<StrawId>,StrawId::_nupanels > _panelOverlap;   // which panels overlap each other
      void genMap();    
      void finalize(ComboHit& combohit);
  };

  MakeStereoHits::MakeStereoHits(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>(           "debugLevel",0)),
    _chTag(pset.get<art::InputTag>("ComboHitCollection","CombineStrawHits")),
    _shsel(pset.get<std::vector<std::string> >("StrawHitSelectionBits",std::vector<std::string>{"EnergySelection","TimeSelection"} )),
    _shmask(pset.get<std::vector<std::string> >("StrawHitMaskBits",std::vector<std::string>{} )),
    _maxDt(pset.get<double>(   "maxDt",40.0)), // nsec
    _maxDE(pset.get<double>(   "maxDE",0.99)), // dimensionless, from 0 to 1
    _maxDZ(pset.get<double>(   "maxDZ",1000.)), // mm, maximum longitudinal distance between straws
    _maxDPerp(pset.get<double>("maxDPerp",500.)), // mm, maximum perpendicular distance between time-division points
    _minDdot(pset.get<double>( "minDdot",0.6)), // minimum angle between straws
    _rmax(pset.get<double>(   "maxRadius",700.0)), // maximum radius (mm)
    _maxChisq(pset.get<double>("maxChisquared",100.0)), // position matching
    _minMVA(pset.get<double>(  "minMVA",0.6)), // MVA cut
    _wfac(pset.get<double>(    "ZErrorFactor",0.3)), // error component due to z separation
    _tfac(pset.get<double>(    "ZErrorFactor",1.0)), // error component due to z separation
    _doMVA(pset.get<bool>(  "doMVA",true)),
    _match(static_cast<matchtype>(pset.get<int>("MatchType",station))),
    _mode(static_cast<mode>(pset.get<int>("Mode",filter))),
    _mvatool(pset.get<fhicl::ParameterSet>("MVATool",fhicl::ParameterSet()))
    {
      _maxChi = sqrt(2.0*_maxChisq);// add some buffer
      // define the mask: straws are in the same unique panel
      std::vector<StrawIdMask::field> fields;
      if(_match==panel)
	fields.push_back(StrawIdMask::panel);
      else
	fields.push_back(StrawIdMask::station);
      _mask = StrawIdMask(fields);
      produces<ComboHitCollection>();
    }

  void MakeStereoHits::beginJob()
  {
    if(_doMVA){
      _mvatool.initMVA();    
      if (_debug > 0) std::cout << "MakeStereoHits MVA parameters: " << std::endl;
      if (_debug > 0) _mvatool.showMVA();
    }
  }

  void MakeStereoHits::beginRun(art::Run & run)
  {
    genMap();
  }

  void MakeStereoHits::produce(art::Event& event) 
  {
// find input
//    auto chH = event.getValidHandle<ComboHitCollection>(_chTag);
//    Can't use ValidHandles so do this manually: what a hack!!
    art::Handle<ComboHitCollection> chH;
    if(!event.getByLabel(_chTag, chH))
      throw cet::exception("RECO")<<"mu2e::MakeStereoHits: No ComboHit collection found for tag" <<  _chTag << endl;
    _chcol = chH.product();
    // setup output
    std::unique_ptr<ComboHitCollection> chcol(new ComboHitCollection());
    chcol->reserve(_chcol->size());
    // reference the parent in the new collection
    chcol->setParent(chH);
    // sort hits by unique panel.  This should be built in by construction upstream FIXME!!
    std::array<std::vector<uint16_t>,StrawId::_nupanels> phits;
    size_t nch = _chcol->size();
    std::vector<bool> used(nch,false);
    for(uint16_t ihit=0;ihit<nch;++ihit){
      ComboHit const& ch = (*_chcol)[ihit];
      // select hits based on flag
      if(_mode == filter ||( ch.flag().hasAllProperties(_shsel) && (!ch.flag().hasAnyProperty(_shmask))) ){
	phits[ch.sid().uniquePanel()].push_back(ihit);
      }
    }
    //  Loop over all hits.  Every one must appear somewhere in the output 
    for (size_t ihit=0;ihit<nch;++ihit) {
      // create an output combo hit for every hit; initialize it with this hit
      ComboHit const& ch1 = (*_chcol)[ihit];
      ComboHit combohit;
      combohit.init(ch1,ihit);
      // loop over the panels which overlap this hit's panel
      for (auto sid : _panelOverlap[ch1.sid().uniquePanel()]) {
      // loop over hits in the overlapping panel
	for (auto jhit : phits[sid.uniquePanel()]) {
	  const ComboHit& ch2 = (*_chcol)[jhit];
	  if (used[jhit] ) continue;
	  double dt = fabs(ch1.time()-ch2.time());
	  if (dt > _maxDt) continue;
	  float de = std::min(1.0f,std::abs((ch1.energyDep() - ch2.energyDep())/(ch1.energyDep()+ch2.energyDep())));
	  if (de > _maxDE ) continue;
	  int sep = abs(ch1.sid().uniquePanel() - ch2.sid().uniquePanel());
	    // hits are in the same station but not the same panel
	  if ( sep == 0 || sep >= 12) continue;
	  double ddot = ch1.wdir().Dot(ch2.wdir());
	  XYZVec dp = ch1.pos()-ch2.pos();
	  double dperp = sqrt(dp.perp2());
	  double dz = fabs(dp.z());
	    // negative crosings are in opposite quadrants and longitudinal separation isn't too big
	  if (ddot < _minDdot || dz > _maxDZ || dperp > _maxDPerp ) continue;
	  // solve for the POCA.
	  TwoLinePCA_XYZ pca(ch1.pos(),ch1.wdir(),ch2.pos(),ch2.wdir());
	  if(pca.closeToParallel()){  
	    cet::exception("RECO")<<"mu2e::StereoHit: parallel wires" << std::endl;
	  }
	  // check the points are inside the tracker active volume
	  double rho1 = sqrt(pca.point1().Perp2());
	  double rho2 = sqrt(pca.point2().Perp2());
	  if(rho1 > _rmax || rho2 > _rmax) continue;
	  // compute chisquared; include error for particle angle
	  // should be a cumulative linear regression FIXME!
	  float terr = _tfac*fabs(ch1.pos().z()-ch2.pos().z());
	  float terr2 = terr*terr;
	  float dw1 = ch1.wireDist()-pca.s1();
	  float dw2 = ch2.wireDist()-pca.s2();
	  float chisq = dw1*dw1/(ch1.wireErr2()+terr2) + dw2*dw2/(ch2.wireErr2()+terr2);
	  if (chisq > _maxChisq) continue;
	  // accumulate the chisquared
	  combohit._qual += chisq;
	  // if we get to here, add the hit
	  bool ok = combohit.addIndex(jhit);
	  if(!ok)std::cout << "MakeStereoHits can't add hit" << std::endl;
	  used[jhit] = true;
	}	
      }
      finalize(combohit);
      chcol->push_back(combohit);
      used[ihit] = true;
    }
    // Finish up

    event.put(std::move(chcol));
  } 

  void MakeStereoHits::finalize(ComboHit& combohit) {
    combohit._mask = _mask;
    combohit._flag.merge(StrawHitFlag::stereo);
    accumulator_set<float, stats<tag::mean> > eacc;
    accumulator_set<float, stats<tag::mean> > tacc;
    accumulator_set<float, stats<tag::min > > zmin;
    accumulator_set<float, stats<tag::max > > zmax;
    accumulator_set<float, stats<tag::weighted_variance(lazy)>, float> xacc,yacc,zacc;
    combohit._nsh = 0;
    for(size_t ich = 0; ich < combohit.nCombo(); ++ich){
      size_t index = combohit.index(ich);
      ComboHit const& ch = (*_chcol)[index];
      float wt = 1.0/(ch.wireErr2());
      xacc(ch.pos().x(),weight=wt); // wire position is weighted
      yacc(ch.pos().y(),weight=wt); // wire position is weighted
      zacc(ch.pos().z(),weight=wt); // wire position is weighted
      combohit._flag.merge(ch.flag());
      eacc(ch.energyDep());
      tacc(ch.time());// time is an unweighted average
      zmin(ch.pos().z());
      zmax(ch.pos().z());
      combohit._nsh += ch.nStrawHits();
    }
    float xp = extract_result<tag::weighted_mean>(xacc);
    float yp = extract_result<tag::weighted_mean>(yacc);
    float zp = extract_result<tag::weighted_mean>(zacc);
    float maxz = extract_result<tag::max>(zmax);
    float minz = extract_result<tag::min>(zmin);
    combohit._pos = XYZVec(xp,yp,zp);
    combohit._time = extract_result<tag::mean>(tacc);
    combohit._edep = extract_result<tag::mean>(eacc);
    combohit._wdist = combohit._pos.z() - 0.5*(maxz+minz); 
    float dz = (maxz-minz);
    float terr = dz*_tfac;
    // component of error
    combohit._tres = sqrt(1.0/extract_result<tag::sum_of_weights>(xacc) + terr*terr);
    combohit._wres = dz*_wfac;
    combohit._qual/= combohit.nCombo();// normalize quality
    combohit._wdir = XYZVec(0.0,0.0,1.0);
  }
// generate the overlap map
  void MakeStereoHits::genMap() {
    static bool init(false);
    if(!init){
      init = true;
      // initialize
      const TTracker& tt(*GeomHandle<TTracker>());
      // establihit the extent of a panel using the longest straw (0)
      Straw const& straw = tt.getStraw(StrawId(0,0,0));
      double phi0 = (straw.getMidPoint()-straw.getHalfLength()*straw.getDirection()).phi();
      double phi1 = (straw.getMidPoint()+straw.getHalfLength()*straw.getDirection()).phi();
      double lophi = std::min(phi0,phi1);
      double hiphi = std::max(phi0,phi1);
      double phiwidth = hiphi-lophi;
      if (phiwidth>M_PI) phiwidth = 2*M_PI-phiwidth;
      if(_debug > 0)std::cout << "Panel Phi width = " << phiwidth << std::endl;
      // loop over all unique panels
      for(size_t ipla = 0;ipla < StrawId::_nplanes; ++ipla) {
	for(int ipan=0;ipan<StrawId::_npanels;++ipan){
	  StrawId sid(ipla,ipan,0);
	  uint16_t upan = sid.uniquePanel();
	  Straw const& straw = tt.getStraw(StrawId(ipla,ipan,0));
	  float phi = straw.getMidPoint().phi();
	  // loop over nearby panels and check for an overlap
	  size_t minpla(ipla), maxpla(ipla);
	  if(_match==station){
	    minpla = (size_t)std::max(0,(int)ipla-1);
	    maxpla = (size_t)std::min(StrawId::_nplanes-1,(int)ipla+1);
	  }
	  for(size_t jpla = minpla; jpla <= maxpla;++jpla){
	    for(int jpan=0;jpan<StrawId::_npanels;++jpan){
	      StrawId osid(jpla,jpan,0);
	      if(osid.uniquePanel() != sid.uniquePanel() && osid.station() == sid.station()){
		Straw const& ostraw = tt.getStraw(StrawId(jpla,jpan,0));
		float dphi = fabs(phi - ostraw.getMidPoint().phi());
		if (dphi > M_PI) dphi = 2*M_PI-dphi;
		if (dphi < phiwidth) _panelOverlap[upan].push_back(osid);
	      }
	    }
	  }
	}
      }
      if (_debug >0)
      {
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


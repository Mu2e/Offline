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
#include "GeneralUtilities/inc/TwoLinePCA.hh"
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
using namespace boost::accumulators;

#include <iostream>
#include <float.h>
using namespace std;
using CLHEP::Hep3Vector;

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
      double         _radBuf;      // radial buffer for stereo to be in active volume
      double         _maxChi;     // maximum # of TimeDivision sigmas past the active edge of a wire to allow making stereo hits
      double         _maxChisq;   // maximum chisquared to allow making stereo hits
      double         _minMVA;     // minimum MVA output
      double         _wres;       // resolution to assign along the wire
      double         _minR;       // min radius flag
      double         _maxR;       // max radius flag
      double         _minE;       // min hit energy
      double         _maxE;       // max hit energy
      bool           _doMVA;      // do MVA eval or simply use chi2 cut
      matchtype	    _match;	  // match either to panel or station
      StrawIdMask _mask; 

      MVATools _mvatool;
      StereoMVA _vmva; 

      std::vector <std::vector<StrawId> > _panelOverlap;   // which panels overlap each other
      void genMap();    
      void combineHits(ComboHit& combohit);
  };

  MakeStereoHits::MakeStereoHits(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>(           "debugLevel",0)),
    _chTag(pset.get<art::InputTag>("ComboHitCollection","makeSH")),
    _shsel(pset.get<std::vector<std::string> >("StrawHitSelectionBits",std::vector<std::string>{"EnergySelection","TimeSelection"} )),
    _shmask(pset.get<std::vector<std::string> >("StrawHitMaskBits",std::vector<std::string>{} )),
    _maxDt(pset.get<double>(   "maxDt",40.0)), // nsec
    _maxDE(pset.get<double>(   "maxDE",0.99)), // dimensionless, from 0 to 1
    _maxDZ(pset.get<double>(   "maxDZ",1000.)), // mm, maximum longitudinal distance between straws
    _maxDPerp(pset.get<double>("maxDPerp",500.)), // mm, maximum perpendicular distance between time-division points
    _minDdot(pset.get<double>( "minDdot",0.6)), // minimum angle between straws
    _radBuf(pset.get<double>(   "radBuf",20.0)), // radial buffer (mm)
    _maxChisq(pset.get<double>("maxChisquared",20.0)), // position matching
    _minMVA(pset.get<double>(  "minMVA",0.6)), // MVA cut
    _wres(pset.get<double>(    "LongitudinalResolution",20.0)), // estimated resolution of stereo reco
    _minR(pset.get<double>("minimumRadius",395.0)), // mm
    _maxR(pset.get<double>("maximumRadius",650.0)), // mm
    _minE(pset.get<double>("minimumEnergy",0.0)), // Minimum deposited straw energy (MeV)
    _maxE(pset.get<double>("maximumEnergy",0.0035)), // MeV
    _doMVA(pset.get<bool>(  "doMVA",true)),
    _match(static_cast<matchtype>(pset.get<int>("MatchType",station))),
    _mvatool(pset.get<fhicl::ParameterSet>("MVATool",fhicl::ParameterSet())),
    _panelOverlap(StrawId::_nupanels)
    {
      _maxChi = sqrt(_maxChisq);
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
    _mvatool.initMVA();    
    if (_debug > 0) std::cout << "MakeStereoHits MVA parameters: " << std::endl;
    if (_debug > 0) _mvatool.showMVA();
  }

  void MakeStereoHits::beginRun(art::Run & run)
  {
    genMap();
  }

  void MakeStereoHits::produce(art::Event& event) 
  {
    Tracker const& tracker = getTrackerOrThrow();
    const TTracker& tt = dynamic_cast<const TTracker&>(tracker);
    float rmax = tt.rOut() + _radBuf;
// find input
//    auto chH = event.getValidHandle<ComboHitCollection>(_chTag);
//    Can't use ValidHandles so do this manually: what a hack!!
    art::Handle<ComboHitCollection> chH;
    if(!event.getByLabel(_chTag, chH))
      throw cet::exception("RECO")<<"mu2e::MakeStereoHits: No ComboHit collection found for tag" <<  _chTag << endl;
    _chcol = chH.product();
    // setup output
    std::unique_ptr<ComboHitCollection> chcol(new ComboHitCollection(*_chcol));
    chcol->reserve(_chcol->size());
    // reference the parent in the new collection
    chcol->setParent(chH);
    // sort hits by unique panel.  This should be built in by construction upstream FIXME!!
    std::array<std::vector<uint16_t>,StrawId::_nupanels> phits;
    size_t nsh = _chcol->size();
    std::vector<bool> used(nsh,false);
    for(uint16_t ihit=0;ihit<nsh;++ihit){
      ComboHit const& ch = (*_chcol)[ihit];
      // select hits based on flag
      if(ch.flag().hasAllProperties(_shsel) && (!ch.flag().hasAnyProperty(_shmask)) ){
	phits[ch.sid().uniquePanel()].push_back(ihit);
      }
    }
    //  Loop over all hits.  Every one must appear somewhere in the output 
    for (size_t ihit=0;ihit<nsh;++ihit)
    {
      // create an output combo hit for every hit; initialize it with this hit
      ComboHit const& ch1 = (*_chcol)[ihit];
      ComboHit combohit;
      combohit.init(ch1,ihit);
      if (used[ihit] || !ch1.flag().hasAllProperties(_shsel) || ch1.flag().hasAnyProperty(_shmask) ) continue;
      // loop over the panels which overlap this hit's panel
      for (auto sid2 : _panelOverlap[ch1.sid().uniquePanel()])
      {
      // loop over hits in the overlapping panel
	for (auto jhit : phits[sid2.uniquePanel()])
	{
	  const ComboHit& ch2 = (*_chcol)[jhit];
	  if (used[jhit] || !ch2.flag().hasAllProperties(_shsel) || ch2.flag().hasAnyProperty(_shmask) ) continue;
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
	  // solve for the POCA.  must translate FIXME!
	  Hep3Vector p1(ch1.pos().x(),ch1.pos().y(),ch1.pos().z());
	  Hep3Vector w1(ch1.wdir().x(),ch1.wdir().y(),ch1.wdir().z());
	  Hep3Vector p2(ch2.pos().x(),ch2.pos().y(),ch2.pos().z());
	  Hep3Vector w2(ch2.wdir().x(),ch2.wdir().y(),ch2.wdir().z());
	  TwoLinePCA pca(p1,w1,p2,w2);
	  if(pca.closeToParallel()){  
	    cet::exception("RECO")<<"mu2e::StereoHit: parallel wires" << std::endl;
	  }
	  // check the points are inside the tracker active volume
	  XYZVec pos(0.5*(pca.point1().x()+pca.point2().x()),
	      0.5*(pca.point1().y()+pca.point2().y()),
	      0.5*(pca.point1().z()+pca.point2().z()));
	  double rho1 = pca.point1().perp();
	  double rho2 = pca.point2().perp();
	  if(rho1 > rmax || rho2 > rmax) continue;
	  // compute chisquared
	  float chi1 = (ch1.wireDist()-pca.s1())/ch1.posRes(ComboHit::wire);
	  float chi2 = (ch2.wireDist()-pca.s2())/ch2.posRes(ComboHit::wire);
	  if (fabs(chi1) >_maxChi || fabs(chi2) > _maxChi) continue;
	  double chisq = chi1*chi1+chi2*chi2; 
	  if (chisq > _maxChisq) continue;
	  // compute MVA
	  double mvaout(-1.0);
	  if (_doMVA)
	  {
	    _vmva._dt = dt;
	    _vmva._chisq = chisq;
	    _vmva._rho = sqrt(pos.Perp2());
	    _vmva._ndof = 2;
	    mvaout = _mvatool.evalMVA(_vmva._pars);
	      if (mvaout < _minMVA) continue;
	  }
	  // if we get to here, add the hit
	  bool ok = combohit.addIndex(jhit);
	  if(!ok)std::cout << "MakeStereoHits past limit" << std::endl;
	  used[jhit] = true;
	}	
      }
      // finalize the position and quality computation and save the hit
      combineHits(combohit);
      chcol->push_back(combohit);
      used[ihit] = true;
    } 

    event.put(std::move(chcol));
  } 

  void MakeStereoHits::genMap()
  {
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
    // loop over all unique panels
    for(size_t ipla = 0;ipla < StrawId::_nplanes; ++ipla) {
      for(int ipan=0;ipan<StrawId::_npanels;++ipan){
	StrawId sid2(ipla,ipan,0);
	uint16_t upan = sid2.uniquePanel();
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
	    StrawId osid2(jpla,jpan,0);
	    if(osid2.uniquePanel() != sid2.uniquePanel() && osid2.station() == sid2.station()){
	      Straw const& ostraw = tt.getStraw(StrawId(ipla,ipan,0));
	      float dphi = fabs(phi - ostraw.getMidPoint().phi());
	      if (dphi > M_PI) dphi = 2*M_PI-dphi;
	      if (dphi < phiwidth) _panelOverlap[upan].push_back(osid2);
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
  // compute the properties of this combined hit
  void MakeStereoHits::combineHits(ComboHit& combohit) {
  // this should be a constrained linear regression fit FIXME!
    if(combohit._nsh > 1){
      combohit._mask = _mask;
      combohit._flag.merge(StrawHitFlag::stereo);
      accumulator_set<float, stats<tag::mean> > eacc;
      accumulator_set<float, stats<tag::mean> > tacc;
      accumulator_set<float, stats<tag::weighted_variance(lazy)>, float> xacc,yacc,zacc;
      accumulator_set<float, stats<tag::mean> > wtacc;
      accumulator_set<float, stats<tag::mean> > werracc;
      XYZVec midpos;
      combohit._nsh = 0;
      for(unsigned ich = 0; ich < combohit.nCombo(); ++ich){
	// get back the original information
	ComboHit const& ch = (*_chcol)[combohit.index(ich)];
	combohit._flag.merge(ch.flag());
	eacc(ch.energyDep());
	tacc(ch.time());// time is an unweighted average
	float wt = 1.0/(ch.wireErr2());
	xacc(ch.pos().x(),weight=wt); // wire position is weighted
	yacc(ch.pos().y(),weight=wt); // wire position is weighted
	zacc(ch.pos().z(),weight=wt); // wire position is weighted
	wtacc(wt);
	werracc(ch.wireRes());
	combohit._nsh += ch.nStrawHits();
      }
      combohit._time = extract_result<tag::mean>(tacc);
      combohit._edep = extract_result<tag::mean>(eacc);
      combohit._wdist = 0.0;
      combohit._wdir = XYZVec(0.0,0.0,1.0);
      float xp = extract_result<tag::weighted_mean>(xacc);
      float yp = extract_result<tag::weighted_mean>(yacc);
      float zp = extract_result<tag::weighted_mean>(zacc);
      combohit._pos = XYZVec(xp,yp,zp);
      combohit._tres = 1.0/sqrt(extract_result<tag::mean>(wtacc));
      combohit._wres = _wres; 
      combohit._qual = 1.0; //FIXME!!
    }
  }
}

using mu2e::MakeStereoHits;
DEFINE_ART_MODULE(MakeStereoHits)


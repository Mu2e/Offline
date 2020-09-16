//
// A module to create simple stereo hits out of StrawHits. StrawHit selection is done by flagging in an upstream module
//
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
#include "art_root_io/TFileService.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eUtilities/inc/TwoLinePCA_XYZ.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
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
      float         _maxDt;      // maximum time separation between hits
      bool          _useTOT;     // use TOT to estimate drift time
      float         _maxDPerp;   // maximum transverse separation
      float         _minDdot;    // minimum dot product of straw directions
      float _minR2, _maxR2; // transverse radius (squared) 
      float         _maxChisq;   // maximum chisquared to allow making stereo hits
      float         _minMVA;     // minimum MVA output
      float         _wfac;       // resolution factor along the wire
      float         _tfac;       // resolution transverse to the wire
      bool           _doMVA;      // do MVA eval or simply use chi2 cut
      unsigned      _maxfsep;	  // max face separation
      bool	    _testflag; // test the flag or not
      StrawIdMask _smask; // define matches inside a station

      MVATools _mvatool;
      StereoMVA _vmva; 

      std::array<std::vector<StrawId>,StrawId::_nupanels > _panelOverlap;   // which panels overlap each other
      void genMap();    
      void finalize(ComboHit& combohit);
  };

  MakeStereoHits::MakeStereoHits(fhicl::ParameterSet const& pset) :
    art::EDProducer{pset},
    _debug(pset.get<int>(           "debugLevel",0)),
    _chTag(pset.get<art::InputTag>("ComboHitCollection")),
    _shsel(pset.get<std::vector<std::string> >("StrawHitSelectionBits",std::vector<std::string>{"EnergySelection","TimeSelection"} )),
    _shmask(pset.get<std::vector<std::string> >("StrawHitMaskBits",std::vector<std::string>{} )),
    _maxDt(pset.get<float>(   "maxDt",40.0)), // nsec //FIXME tune with TOT
    _useTOT(pset.get<bool>("UseTOT",false)), // use TOT to estimate drift time
    _maxDPerp(pset.get<float>("maxDPerp",500.)), // mm, maximum perpendicular distance between time-division points
    _minDdot(pset.get<float>( "minDdot",0.6)), // minimum angle between straws
    _maxChisq(pset.get<float>("maxChisquared",5.0)), // position matching
    _minMVA(pset.get<float>(  "minMVA",0.6)), // MVA cut
    _wfac(pset.get<float>(    "ZErrorFactor",0.3)), // error component due to z separation
    _tfac(pset.get<float>(    "ZErrorFactor",1.0)), // error component due to z separation
    _doMVA(pset.get<bool>(  "doMVA",false)),
    _maxfsep(pset.get<unsigned>("MaxFaceSeparation",3)), // max separation between faces in a station
    _testflag(pset.get<bool>("TestFlag")),
    _mvatool(pset.get<fhicl::ParameterSet>("MVATool",fhicl::ParameterSet()))
    {
      float minR = pset.get<float>("minimumRadius",395); // mm
      _minR2 = minR*minR;
      float maxR = pset.get<float>("maximumRadius",650); // mm
      _maxR2 = maxR*maxR;
      // define the mask: straws are in the same unique panel
      std::vector<StrawIdMask::field> fields;
      fields.push_back(StrawIdMask::station);
      _smask = StrawIdMask(fields);
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

  void MakeStereoHits::produce(art::Event& event) {
// find input: I have to get a Handle, not ValidHandle, to get the productID
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
    if(_debug > 1)cout << "MakeStereoHits found " << nch << " Input hits" << endl;
    std::vector<bool> used(nch,false);
    for(uint16_t ihit=0;ihit<nch;++ihit){
      ComboHit const& ch = (*_chcol)[ihit];
      // select hits based on flag
      if( (!_testflag) ||( ch.flag().hasAllProperties(_shsel) && (!ch.flag().hasAnyProperty(_shmask))) ){
	phits[ch.strawId().uniquePanel()].push_back(ihit);
      }
    }
    if(_debug > 2){
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
      ComboHit const& ch1 = (*_chcol)[ihit];
      ComboHit combohit;
      combohit.init(ch1,ihit);
      // zero values that accumulate in pairs
      combohit._qual = 0.0;
      combohit._pos = XYZVec(0.0,0.0,0.0);
      // loop over the panels which overlap this hit's panel
      for (auto sid : _panelOverlap[ch1.strawId().uniquePanel()]) {
      // loop over hits in the overlapping panel
	for (auto jhit : phits[sid.uniquePanel()]) {
	  const ComboHit& ch2 = (*_chcol)[jhit];
	  if(_debug > 3) cout << " comparing hits " << ch1.strawId().uniquePanel() << " and " << ch2.strawId().uniquePanel();
	  if (!used[jhit] ){
            float dt;
            if (_useTOT)
	      dt = fabs(ch1.correctedTime()-ch2.correctedTime());
            else
	      dt = fabs(ch1.time()-ch2.time());
	    if(_debug > 3) cout << " dt = " << dt;
	    if (dt < _maxDt){
	      float ddot = ch1.wdir().Dot(ch2.wdir());
	      XYZVec dp = ch1.pos()-ch2.pos();
	      float dperp = sqrt(dp.perp2());
	      // negative crosings are in opposite quadrants and longitudinal separation isn't too big
	      if(_debug > 3) cout << " ddot = " << ddot << " dperp = " << dperp;
	      if (ddot > _minDdot && dperp < _maxDPerp ) {
		// solve for the POCA.
		TwoLinePCA_XYZ pca(ch1.pos(),ch1.wdir(),ch2.pos(),ch2.wdir());
		if(pca.closeToParallel()){  
		  cet::exception("RECO")<<"mu2e::StereoHit: parallel wires" << std::endl;
		}
		// check the points are inside the tracker active volume; these are all the same as the
		float rho2 = pca.point1().Perp2();
		if(_debug > 3) cout << " rho2 = " << rho2;
		if(rho2 < _maxR2 && rho2 > _minR2 ){
		  // compute chisquared; include error for particle angle
		  // should be a cumulative linear regression FIXME!
		  float terr = _tfac*fabs(ch1.pos().z()-ch2.pos().z());
		  float terr2 = terr*terr;
		  float dw1 = pca.s1();
		  float dw2 = pca.s2();
		  float chisq = dw1*dw1/(ch1.wireErr2()+terr2) + dw2*dw2/(ch2.wireErr2()+terr2);
		  if(_debug > 3) cout << " chisq = " << chisq;
		  if (chisq < _maxChisq){
		    if(_debug > 3) cout << " added ";
		    // if we get to here, try to add the hit
		    // accumulate the chisquared
		    if(combohit.addIndex(jhit)) {
		      // average z 
		      combohit._qual += chisq;
		      combohit._pos += XYZVec(pca.point1().x(),pca.point1().y(),0.5*(pca.point1().z()+pca.point2().z()));	    
		    } else
		      std::cout << "MakeStereoHits can't add hit" << std::endl;
		    used[jhit] = true;
		  }	
		}
	      }
	    }
	  }
	  if(_debug > 3) cout << endl;
	}
      }
      finalize(combohit);
      chcol->push_back(std::move(combohit));
    }
    event.put(std::move(chcol));
  } 

  void MakeStereoHits::finalize(ComboHit& combohit) {
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
	ComboHit const& ch = (*_chcol)[index];
	combohit._flag.merge(ch.flag());
	eacc(ch.energyDep(),weight=ch.nStrawHits());
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
      combohit._edep = extract_result<tag::weighted_mean>(eacc);
      float dz = (maxz-minz);
      combohit._wdist = dz; 
      combohit._tres = dz*_tfac;
      combohit._wres = dz*_wfac;
      combohit._qual /= (combohit.nCombo()-1);// normalize by # of pairs
      combohit._pos /= (combohit.nCombo()-1);
      combohit._wdir = XYZVec(0.0,0.0,1.0);
    } else {
      size_t index = combohit.index(0);
      combohit._pos = (*_chcol)[index].pos();// put back original position
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
      if(_debug > 0)std::cout << "Panel Phi width = " << phiwidth << std::endl;
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
	      if(_smask.equal(osid,sid) && osid.uniqueFace() != sid.uniqueFace() && (unsigned)abs(osid.uniqueFace() - sid.uniqueFace()) <= _maxfsep ) {
		float dphi = fabs(phi - ostraw.getMidPoint().phi());
		if (dphi > M_PI) dphi = 2*M_PI-dphi;
		if (dphi < phiwidth) _panelOverlap[upan].push_back(osid);
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


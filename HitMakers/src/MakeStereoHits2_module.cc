// Can we speed up StereoHit constructor?
// We create additional position and flag collections, do everything into one module?
//
// Need to fix MVA tool to use simple cuts
//
// A module to create simple stereo hits out of StrawHits. StrawHit selection is done by flagging in an upstream module
//
// $Id: MakeStereoHits2_module.cc,v 1.23 2014/09/18 08:42:47 brownd Exp $
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
#include "TrackerGeom/inc/Tracker.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitPosition.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "Mu2eUtilities/inc/MVATools.hh"

#include <iostream>
#include <float.h>


//Notes: Main contributors: MakeStereoHits and MVATools. 
//       Could we cache the panel separation?
//       StereoHit <-> Single Hit position 

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


  class MakeStereoHits2 : public art::EDProducer {

     public:

       explicit MakeStereoHits2(fhicl::ParameterSet const& pset);

       void produce( art::Event& e);
       virtual void beginJob();
       virtual void beginRun(art::Run & run);


     private:

       typedef std::vector<size_t> PanelHits;
       typedef std::vector<PanelHits> StationPanels;
       typedef std::vector< StationPanels > TrackerStations;

       int            _debug;
       art::InputTag  _shTag;
       art::InputTag  _shpTag;
       art::InputTag  _shfTag;

       StrawHitFlag   _shsel;      // flag selection
       StrawHitFlag   _shmask;     // flag anti-selection 
       double         _maxDt;      // maximum time separation between hits
       double         _maxDE;      // maximum deposited energy deference: this excludes inconsistent hits
       double         _maxDZ;      // maximum longitudinal separation
       double         _maxDPerp;   // maximum transverse separation
       double         _minDdot;    // minimum dot product of straw directions
       double         _minDL;      // minimum distance from end of active straw;
       double         _maxChi;     // maximum # of TimeDivision sigmas past the active edge of a wire to allow making stereo hits
       double         _maxChisq;   // maximum # of TimeDivision consistency chisquared to allow making stereo hits
       double         _minMVA;     // minimum MVA output
       double         _wres;       // resolution to assign along the wire
       double         _minR;       // min radius flag
       double         _maxR;       // max radius flag
       double         _minE;       // min hit energy
       double         _maxE;       // max hit energy
       bool           _doMVA;      // do MVA eval or simply use chi2 cut
       bool           _bestpair;   // require both hits agree as to the best pair
       bool           _writepairs; // write out the stereo pairs

       MVATools _mvatool;
       StereoMVA _vmva; 

       size_t _nStations; 
       size_t _nPanels; 
       std::vector <std::pair<int,int> > _panelOverlap;  

       void genMap();    
       double longRes(const StereoHit& sthit) const; 
       bool betterPair(const StereoHit& newpair, const StereoHit& oldpair) const;
       void reportInsertedHits(const TrackerStations& stax);
  };

  MakeStereoHits2::MakeStereoHits2(fhicl::ParameterSet const& pset) :
     _debug(pset.get<int>(           "debugLevel",0)),
     _shTag(pset.get<art::InputTag>( "StrawHitCollection","makeSH")),
     _shpTag(pset.get<art::InputTag>("StrawHitPositionCollection","makeSH")),
     _shfTag(pset.get<art::InputTag>("StrawHitFlagCollection","makeSH")),
     _shsel(pset.get<std::vector<std::string> >("StrawHitSelectionBits",std::vector<std::string>{"EnergySelection","TimeSelection"} )),
     _shmask(pset.get<std::vector<std::string> >("StrawHitMaskBits",std::vector<std::string>{} )),
     _maxDt(pset.get<double>(   "maxDt",40.0)), // nsec
     _maxDE(pset.get<double>(   "maxDE",0.99)), // dimensionless, from 0 to 1
     _maxDZ(pset.get<double>(   "maxDZ",1000.)), // mm, maximum longitudinal distance between straws
     _maxDPerp(pset.get<double>("maxDPerp",500.)), // mm, maximum perpendicular distance between time-division points
     _minDdot(pset.get<double>( "minDdot",0.6)), // minimum angle between straws
     _minDL(pset.get<double>(   "minDL",-20.0)), // extent along straw
     _maxChisq(pset.get<double>("maxChisquared",80.0)), // position matching
     _minMVA(pset.get<double>(  "minMVA",0.6)), // MVA cut
     _wres(pset.get<double>(    "LongitudinalResolution",20.0)), // estimated resolution of stereo reco
     _minR(pset.get<double>("minimumRadius",395.0)), // mm
     _maxR(pset.get<double>("maximumRadius",650.0)), // mm
     _minE(pset.get<double>("minimumEnergy",0.0)), // Minimum deposited straw energy (MeV)
     _maxE(pset.get<double>("maximumEnergy",0.0035)), // MeV
     _doMVA(pset.get<bool>(  "doMVA",true)),
     _bestpair(pset.get<bool>(  "BestStereoPair",false)),
     _writepairs(pset.get<bool>("WriteStereoPairs",true)),
     _mvatool(pset.get<fhicl::ParameterSet>("MVATool",fhicl::ParameterSet()))
  {
      _maxChi = sqrt(_maxChisq);

      if (_writepairs) produces<StereoHitCollection>();
      produces<StrawHitPositionCollection>();
      produces<StrawHitFlagCollection>();
  }

  void MakeStereoHits2::beginJob()
  {
     _mvatool.initMVA();    
     if (_debug > 0) std::cout << "MakeStereoHits2 MVA parameters: " << std::endl;
     if (_debug > 0) _mvatool.showMVA();
  }

  void MakeStereoHits2::beginRun(art::Run & run)
  {
      genMap();
  }


  void MakeStereoHits2::produce(art::Event& event) 
  {    
     const TTracker& tt(*GeomHandle<TTracker>());
     
     art::Handle<StrawHitCollection> strawHitsHandle;
     event.getByLabel(_shTag,strawHitsHandle);
     const StrawHitCollection& _shcol(*strawHitsHandle);
     size_t nsh = _shcol.size();

     art::Handle<StrawHitPositionCollection> strawHitPosHandle;
     event.getByLabel(_shpTag,strawHitPosHandle);
     const StrawHitPositionCollection& _shpcol(*strawHitPosHandle);

     art::Handle<StrawHitFlagCollection> strawHitFlagsHandle;
     event.getByLabel(_shfTag,strawHitFlagsHandle);
     const StrawHitFlagCollection& _shfcol(*strawHitFlagsHandle);


     std::unique_ptr<StrawHitFlagCollection> shfcol(new StrawHitFlagCollection);      
     std::unique_ptr<StrawHitPositionCollection> shpcol(new StrawHitPositionCollection(_shpcol));
     std::unique_ptr<StereoHitCollection> stereohits(new StereoHitCollection);
     shfcol->reserve(nsh);
     stereohits->reserve(3*nsh);
     std::vector<int> ibest(nsh,-1);

     size_t nres = std::max(size_t(32),nsh/10);
     PanelHits hpstax;
     hpstax.reserve(nres);
     StationPanels pstax(2*_nPanels,hpstax);
     TrackerStations stax(_nStations,pstax);

     for (size_t ish=0;ish<nsh;++ish)
     {
         StrawHitFlag shf = _shfcol.at(ish);
         shf.merge(shpcol->at(ish).flag());
         if (!shf.hasAllProperties(_shsel) || shf.hasAnyProperty(_shmask) ) continue;

	 const StrawHit& hit = _shcol.at(ish);
	 const Straw& straw  = tt.getStraw(hit.strawIndex());

         size_t iplane  = straw.id().getPlane();
         size_t ipnl    = straw.id().getPanel();
         size_t station = iplane/2;

         // define a 'global' panel for the station.  This changes sign with odd-even stations
	 size_t jpnl = (station%2==0) ? ipnl + (iplane%2)*_nPanels : ipnl + (1-iplane%2)*_nPanels;
	 if ( _debug > 2) std::cout << "Inserting hit " << ish << " into station " << station << " panel " << jpnl << std::endl;
         stax[station][jpnl].push_back(ish);      
     }


    // should be able to rewrite the loop to avoid the loop on station and overlaps
    
    // loop over stations
    for (const StationPanels& pstax : stax)
    {                                
       for (const auto& ij : _panelOverlap)
       {
          int ipnl = ij.first;
          int jpnl = ij.second;
          if (pstax[jpnl].empty()) continue;	      

          for (size_t ish : pstax[ipnl])
          {                                    		
             const StrawHit& sh1 = _shcol.at(ish);
	     const Straw& straw1 = tt.getStraw(sh1.strawIndex());
	     const StrawHitPosition& shp1 = _shpcol.at(ish);
             double t1 = sh1.time();

             for (size_t jsh : pstax[jpnl])
             {                              
 	        double dt = fabs(t1-_shcol.at(jsh).time());
                if (dt > _maxDt) continue;

                const StrawHit& sh2 = _shcol.at(jsh);
	        const Straw& straw2 = tt.getStraw(sh2.strawIndex());
	        const StrawHitPosition& shp2 = _shpcol.at(jsh);

                float de = std::min(1.0f,std::abs((sh1.energyDep() - sh2.energyDep())/(sh1.energyDep()+sh2.energyDep())));
                if (de > _maxDE ) continue;


                PanelId::isep sep = straw1.id().getPanelId().separation(straw2.id().getPanelId());
                // hits are in the same station but not the same panel
                if ( sep == PanelId::same || sep >= PanelId::apart) continue;


                double ddot = straw1.direction().dot(straw2.direction());
	        CLHEP::Hep3Vector dp = shp1.pos()-shp2.pos();
	        double dperp = dp.perp();
	        double dz = fabs(dp.z());

                // negative crosings are in opposite quadrants and longitudinal separation isn't too big
                if (ddot < _minDdot || dz > _maxDZ || dperp > _maxDPerp ) continue;


	        // tentative stereo hit: this solves for the POCA
	        StereoHit sth(_shcol,tt,ish,jsh);
	        double dl1 = straw1.getDetail().activeHalfLength()-fabs(sth.wdist1());
	        double dl2 = straw2.getDetail().activeHalfLength()-fabs(sth.wdist2());
                if (dl1 < _minDL || dl2 < _minDL) continue;

		double chi1(0.0), chi2(0.0);
		unsigned ndof(0);
		if (shp1.flag().hasAllProperties(StrawHitFlag::tdiv))
                {
		  chi1 = (shp1.wireDist()-sth.wdist1())/shp1.posRes(StrawHitPosition::wire);
		  ++ndof;
		}
		if (shp2.flag().hasAllProperties(StrawHitFlag::tdiv)){
		  chi2 = (shp2.wireDist()-sth.wdist2())/shp2.posRes(StrawHitPosition::wire);
		  ++ndof;
		}
		if (fabs(chi1) >_maxChi || fabs(chi2) > _maxChi) continue;

		double chisq = chi1*chi1+chi2*chi2; 
		if (chisq > _maxChisq) continue;
		sth.setChisquared(chisq);

                double mvaout(-1.0);
                if (_doMVA)
                {
                   _vmva._dt = dt;
		   _vmva._chisq = chisq;
		   _vmva._rho = sth.pos().perp();
		   _vmva._ndof = ndof;
		   mvaout = _mvatool.evalMVA(_vmva._pars);
                   if (mvaout < _minMVA) continue;
                }
		
		sth.setMVAOut(mvaout);
		stereohits->push_back(sth);
		size_t isth = stereohits->size()-1;

		if (ibest[ish] < 0 || betterPair(sth,stereohits->at(ibest[ish]))) ibest[ish] = isth;
		if (ibest[jsh] < 0 || betterPair(sth,stereohits->at(ibest[jsh]))) ibest[jsh] = isth;
	     } 
	  } 
       } 
    } 
    

    // overwrite the positions for stereoHits and update flags
    for (size_t ish=0; ish<nsh;++ish)
    {
       if (ibest[ish] > 0)
       {
	  const StereoHit& sthit = stereohits->at(ibest[ish]);
	  // optionally require this be the best pairing for both hits
	  if( (!_bestpair) || ibest[sthit.hitIndex1()] == ibest[sthit.hitIndex2()])
          {
	     StrawHitPosition& shpos = shpcol->at(ish);
	     shpos._flag.merge(StrawHitFlag::stereo);
	     shpos._stindex = ibest[ish];
	     shpos._pos = sthit.pos();
	     shpos._phi = sthit.pos().phi();
	     shpos._wres = longRes(sthit);
          }
       }

       //recreate full flag list and update the position / energy flag. The energy flag must be cleared first
       const StrawHitPosition& shp = shpcol->at(ish);
       double rad = shp.pos().perp();
       double energy = _shcol.at(ish).energyDep();

       StrawHitFlag flag(_shfcol.at(ish));  //merge old flag 
       flag.merge(_shpcol.at(ish).flag());  //merge old straw hitposition flag
       flag.merge(shpcol->at(ish)._flag);   //merge stereo flag     
       flag.clear(StrawHitFlag::energysel); //clear energysel flag
       if (rad > _minR && rad < _maxR) flag.merge(StrawHitFlag::radsel);
       if (energy > _minE && energy < _maxE) flag.merge(StrawHitFlag::energysel);
       shfcol->push_back(std::move(flag));
    }



    if (_writepairs) event.put(std::move(stereohits));
    event.put(std::move(shpcol));
    event.put(std::move(shfcol));
  } 

  
  // estimate the resolution on the stereo hit position projection along a wire direction
  double MakeStereoHits2::longRes(StereoHit const& sthit) const 
  {
     return _wres; // this should be a real calculation FIXME!!!
  }


  bool MakeStereoHits2::betterPair(StereoHit const& newpair, StereoHit const& oldpair) const 
  {
     // choose the best pair as:
     // 1) take the pair with the minimum plane separation
     // 2) otherwise, take the pair with the maximum MVA output
     if(newpair.panelSeparation() == oldpair.panelSeparation() ) {
       return newpair.mvaout() > oldpair.mvaout();
     } else {
       return newpair.panelSeparation() < oldpair.panelSeparation();
     }
  }


  void MakeStereoHits2::genMap()
  {
      const TTracker& tt(*GeomHandle<TTracker>());

      _nStations = tt.nPlanes()/2;
      _nPanels = tt.getPlane(0).nPanels(); 

      const Straw& straw = tt.getStraw(StrawId(0,0,0,0));
      double phi0 = (straw.getMidPoint()-straw.getHalfLength()*straw.getDirection()).phi();
      double phi1 = (straw.getMidPoint()+straw.getHalfLength()*straw.getDirection()).phi();
      double lophi = std::min(phi0,phi1);
      double hiphi = std::max(phi0,phi1);
      double phiwidth = hiphi-lophi;
      if (phiwidth>M_PI) phiwidth = 2*M_PI-phiwidth;

      std::vector<double> panphi(2*_nPanels);

      for(int ipla=0;ipla<2;++ipla)
      {
         Plane const& plane = tt.getPlane(ipla);
         for(int ipan=0;ipan<plane.nPanels();++ipan)
         {
	    Panel const& panel = plane.getPanel(ipan);
	    // expand to station-wide 'panel' number.  This changes sign with station
	    int jpan = ipan + (ipla%2)*_nPanels;
	    panphi[jpan] = panel.straw0MidPoint().phi();
         }
      }
      
      for (size_t iphi = 0;iphi<2*_nPanels;++iphi)
      {
	  double phi = panphi[iphi];
	  for (size_t jphi=iphi+1;jphi<2*_nPanels;++jphi)
          {
	     double dphi = fabs(phi - panphi[jphi]);
	     if (dphi > M_PI) dphi = 2*M_PI-dphi;
	     if (dphi < phiwidth) _panelOverlap.emplace_back(std::make_pair<int,int>(iphi,jphi));
	  }
      }

      if (_debug >0)
      {
         // for testing, assume every panel overlaps
         if (_debug > 9) 
         {          
            _panelOverlap.clear();
            for (size_t iphi = 0;iphi<2*_nPanels;++iphi)
	      for (size_t jphi=iphi+1;jphi<2*_nPanels;++jphi) 
                _panelOverlap.emplace_back(std::make_pair<int,int>(iphi,jphi));
         }
                  
         std::cout << "panel phi width = " << phiwidth << "  panel phi positions = ";
         for(auto pphi : panphi) std::cout << pphi << " ";
         std::cout << std::endl;
 
         for (auto& kv : _panelOverlap)         
	   std::cout << "Panel " << kv.first << " Overlaps with panel "<<kv.second<<std::endl;         
      }
  
  }


  void MakeStereoHits2::reportInsertedHits(const TrackerStations& stax)
  {
      if (stax.empty()) return;
      if ( !stax[0].empty()) std::cout << "Built PanelHits size = " << stax[0][0].size() << std::endl;
      std::cout << "Built StationPanels size = " << stax[0].size() << std::endl;

      std::cout << "Built TrackerStations size = " << stax.size() << std::endl;
      for( size_t station =0; station < stax.size(); ++station)     
	  for(size_t ipnl = 0; ipnl < stax[station].size(); ++ipnl) 
	    for(size_t ish : stax[station][ipnl]) 	   
 	       std::cout << "Inserted hit " << ish << " into station " << station << " panel " << ipnl << std::endl;	 	

  }


} 

using mu2e::MakeStereoHits2;
DEFINE_ART_MODULE(MakeStereoHits2)


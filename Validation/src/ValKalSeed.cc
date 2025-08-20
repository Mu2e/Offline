#include "Offline/Validation/inc/ValKalSeed.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/KalSeedMC.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include <cmath>

namespace mu2e {

  ValKalSeed::ValKalSeed(std::string name) : _name(name) {
    // build the VDId -> SId map by hand
    _vdmap[VirtualDetectorId(VirtualDetectorId::TT_FrontHollow)] = SurfaceId("TT_Front");
    _vdmap[VirtualDetectorId(VirtualDetectorId::TT_Mid)] = SurfaceId("TT_Mid");
    _vdmap[VirtualDetectorId(VirtualDetectorId::TT_MidInner)] = SurfaceId("TT_Mid");
    _vdmap[VirtualDetectorId(VirtualDetectorId::TT_Back)] = SurfaceId("TT_Back");
    _vdmap[VirtualDetectorId(VirtualDetectorId::TT_OutSurf)] = SurfaceId("TT_Outer");
    _vdmap[VirtualDetectorId(VirtualDetectorId::TT_InSurf)] = SurfaceId("TT_Inner");
  }

  int ValKalSeed::declare(const art::TFileDirectory& tfs) {
    _hVer = tfs.make<TH1D>("Ver", "Version Number", 101, -0.5, 100.0);
    _hN = tfs.make<TH1D>("NSeed", "N KalSeed", 11, -0.5, 10.0);
    _hNStraw = tfs.make<TH1D>("NHit", "N Straw Hits", 101, -0.5, 100.0);
    _hNSeg = tfs.make<TH1D>("NSeg", "N KalSegment", 100, 0.0, 500.0);
    _hNInter = tfs.make<TH1D>("NInter", "N Intersections", 51, -0.5, 50.0);
    _hTraj = tfs.make<TH1D>("Traj","Trajectory Type", 5,0.5,5.5);
    _hStatus = tfs.make<TH1D>("Status", "Status", 32, -0.5, 31.5);
    _hchi2 = tfs.make<TH1D>("Chi2N", "Chi2/DOF", 100, 0.0, 10.0);
    _hhasCal = tfs.make<TH1D>("hasCal", "CalCluster attached", 2, -0.5, 1.5);
    _hactiveCal = tfs.make<TH1D>("activeCal", "Calo Hit Active", 2, -0.5, 1.5);
    _hfitCon = tfs.make<TH1D>("FitConn", "Fit CL", 100, 0.0, 1.0);
    _hfitConC = tfs.make<TH1D>("FitConnC", "Fit CL CPR", 100, 0.0, 1.0);
    _hfitConT = tfs.make<TH1D>("FitConnT", "Fit CL TPR", 100, 0.0, 1.0);
    _hp = tfs.make<TH1D>("p", "p", 200, -110., 110.);
    _hp2 = tfs.make<TH1D>("p2", "p", 300, -500., 500.);
    _hpC = tfs.make<TH1D>("pC", "p CPR", 100, 0., 110.);
    _hpT = tfs.make<TH1D>("pT", "p TPR", 100, 0., 110.);
    _hpce = tfs.make<TH1D>("pce", "p CE", 100, 95.0, 110.);
    _hpcep = tfs.make<TH1D>("pcep", "p CE+", 100, 82.0, 97.);
    _hsignedp = tfs.make<TH1D>("signedp", "signedp", 200, -110., 110.);
    _hsignedp2 = tfs.make<TH1D>("signedp2", "signedp2", 300, -500., 500.);
    _hpe = tfs.make<TH1D>("pe", "p error", 100, 0.0, 1.0);
    _hRho = tfs.make<TH1D>("rho", "Transverse radius", 100, 0.0, 800.);
    _hPhi = tfs.make<TH1D>("phi", "phi", 100, -M_PI, M_PI);
    _hCost = tfs.make<TH1D>("cost", "cos(Theta)", 100, -1.0, 1.0);
    _ht0 = tfs.make<TH1D>("t0", "t0", 100, 0.0,1800);
    _ht0e = tfs.make<TH1D>("t0e", "t0 error", 100, 0.0,5.0);
    _ht02 = tfs.make<TH1D>("t02", "t0", 100, 0.0,1e5);
    _hCuts = tfs.make<TH1D>("Cuts", "Cut series", 8, 0.5, 8.5);
    _hCCdisk = tfs.make<TH1D>("CCdisk", "Calo Disk", 2, -0.5, 1.5);
    _hCCEoverP = tfs.make<TH1D>("CCEoverP", "Calo E Over Track P", 100, 0.0, 1.5);
    _hCCDt = tfs.make<TH1D>("CCDt", "Calo time residual", 100, -5.0, 5.0);
    _hCCDOCA = tfs.make<TH1D>("CCDOCA", "Calo DOCA to Track", 100, -100.0, 100.0);
    _hCCcdepth = tfs.make<TH1D>("CCcdepth", "Calo POCA Depth", 100, -100.0, 500.0);
    _hCCtz = tfs.make<TH1D>("CCtz", "Calo POCA Track Z", 100, 1000.0, 5000.0);
    _hHDrift = tfs.make<TH1D>("HDrift", "Hit Drift Radius;drift radius (mm)", 100, 0.0, 3.0);
    _hHDOCA = tfs.make<TH1D>("HDOCA", "Hit Wire DOCA;DOCA (mm)", 100, -5.0, 5.0);
    _hHEDep = tfs.make<TH1D>("HEDep", "Hit Energy Deposition;EDep (KeV)", 100, 0, 5.0);
    _hHPanel = tfs.make<TH1D>("HPanel", "Hit Unique Panel", 216, -0.5, 215.5);
    _hSRadLen = tfs.make<TH1D>("SRadLen", "Fractional Straw Radiation Length", 100, 0, 1.0e-3);
    _hSRadLenSum = tfs.make<TH1D>( "SRadLenSum", "Sum Fractional Straw Radiation Length", 100, 0, 0.04);
    int ibin = 1;
    _hCuts->GetXaxis()->SetBinLabel(ibin++, "MC Primary");  // bin 1, first visible
    _hCuts->GetXaxis()->SetBinLabel(ibin++, "MC Momentum");
    _hCuts->GetXaxis()->SetBinLabel(ibin++, "KalmanOK");  // 3
    _hCuts->GetXaxis()->SetBinLabel(ibin++, "Fit Quality");
    _hCuts->GetXaxis()->SetBinLabel(ibin++, "Livegate");  // 5
    _hCuts->GetXaxis()->SetBinLabel(ibin++, "Reco pitch");
    _hCuts->GetXaxis()->SetBinLabel(ibin++, "Cosmic Rejection");
    _hCuts->GetXaxis()->SetBinLabel(ibin++, "Momentum window");  // 8
    _hCuts->SetMinimum(0.0);
    _hTraj->GetXaxis()->SetBinLabel(1, "LoopHelix");
    _hTraj->GetXaxis()->SetBinLabel(2, "CentralHelix");
    _hTraj->GetXaxis()->SetBinLabel(3, "KinematicLine");
    _hPRes = tfs.make<TH1D>("PRes", "Relative Momentum resolution at tracker mid", 200, -3.0, 3.0);
    _hPResA = tfs.make<TH1D>("PResA", "Absolute Momentum resolution", 200, -5.0, 3.0);

    return 0;
  }

  int ValKalSeed::fill(const KalSeedCollection& coll, art::Event const& event) {
    auto const& ptable = GlobalConstantsHandle<ParticleDataList>();
    // increment this by 1 any time the defnitions of the histograms or the
    // histogram contents change, and will not match previous versions
    _hVer->Fill(12.0);

    _hN->Fill(coll.size());
    for (auto const& ks : coll) {
      _hNStraw->Fill(ks.hits().size());
      _hNSeg->Fill(ks.segments().size());
      _hNInter->Fill(ks.intersections().size());
      if(ks.loopHelixFit())_hTraj->Fill(1);
      if(ks.centralHelixFit())_hTraj->Fill(2);
      if(ks.kinematicLineFit())_hTraj->Fill(3);
      const TrkFitFlag& tff = ks.status();
      bool isCPR = tff.hasAllProperties(TrkFitFlag::CPRHelix);
      bool isTPR = tff.hasAllProperties(TrkFitFlag::TPRHelix);

      //    for(TrkFitFlagDetail::bit_type i=0; i<f.size(); i++)
      //  if(f.hasAnyProperty(i)) _hStatus->Fill(i);

      for (auto sn : tff.bitNames()) {
        if (tff.hasAnyProperty(TrkFitFlag(sn.first)))
          _hStatus->Fill(std::log2(sn.second));
      }

      _hchi2->Fill(ks.chisquared() /ks.nDOF());
      int q = ks.hasCaloCluster();
      _hhasCal->Fill(q);
      q = ks.caloHit()._flag.hasAllProperties(StrawHitFlag::active);
      _hactiveCal->Fill(q);
      _hfitCon->Fill(ks.fitConsistency());
      if (isCPR) _hfitConC->Fill(ks.fitConsistency());
      if (isTPR) _hfitConT->Fill(ks.fitConsistency());

      // sample fit at an appropriate intersection
      VirtualDetectorId vdid = VirtualDetectorId(VirtualDetectorId::TT_Mid);
      if(!ks.loopHelixFit()){
        if(ks.intersections(SurfaceId("TT_Outer")).size() > 0){
          vdid = VirtualDetectorId::TT_OutSurf;
        } else if(ks.intersections(SurfaceId("TT_Back")).size() > 0){
          vdid = VirtualDetectorId::TT_Back;
        } else
          vdid = VirtualDetectorId::TT_FrontHollow;
      }
      //
      double t0 = ks.t0Val();
      double t0var = ks.t0Var();
      _ht0->Fill(t0);
      _ht02->Fill(t0);
      _ht0e->Fill(sqrt(t0var));

      // stop at the first found intersection
      // p of MC associated particle at this intersection
      double p_pri(0.0);
      double p_mc = mcTrkP(event,vdid,p_pri);
      SurfaceId sid = _vdmap[vdid];
      auto kintercol = ks.intersections(sid);
      double ksCharge = ptable->particle(ks.particle()).charge();
      auto ikinter = ks.t0Segment(t0);
      if(ikinter != ks.segments().end()){
        auto mom3 = ikinter->momentum3();
        double p = mom3.R();
        double recoCharge(1.);
        auto   kseg = ks.segments().back();
        if(ks.centralHelixFit()) recoCharge = kseg.centralHelix().charge();
        _hp->Fill(p*recoCharge);
        _hp2->Fill(p*recoCharge);
        if (isCPR) _hpC->Fill(p);
        if (isTPR) _hpT->Fill(p);
        _hpce->Fill(p);
        _hpcep->Fill(p);
        _hsignedp->Fill(p*ksCharge);
        _hsignedp2->Fill(p*ksCharge);
        _hpe->Fill(ikinter->momerr());
        _hRho->Fill(ikinter->position3().Rho());
        _hPhi->Fill(mom3.Phi());
        double cost = cos(mom3.Theta());
        _hCost->Fill(cost);

        // fill the cut series; this needs updating TODO
        bool d0cut =true;
        // the first of the cut series, number of events
        _hCuts->Fill(1.0);
        // MC CE found
       double td = 1.0/tan(mom3.Theta());
        // note: these are crude and arbitrary cuts just for validation,
        // do not intepret these as a correct analysis selection!
        _hPResA->Fill(p - p_pri);
        if (p_mc > 90.0) {
          _hCuts->Fill(2.0);
          if (tff.hasAllProperties(TrkFitFlag("KalmanOK"))) {
            _hCuts->Fill(3.0);
            if (ks.fitConsistency() > 0.002) {
              _hCuts->Fill(4.0);
              if (t0 > 700.0 && t0 < 1695.0) {
                _hCuts->Fill(5.0);
                if (td > 0.577 && td < 1.0) {
                  _hCuts->Fill(6.0);
                  if (d0cut) {
                    _hCuts->Fill(7.0);
                    _hPRes->Fill(p - p_mc);
                  }
                }
              }
            }
          }
        }
      }
      // associated cluster info
      if (ks.hasCaloCluster()) {
        auto const& chs = ks.caloHit();
        _hCCdisk->Fill(chs.caloCluster()->diskID());
        // get momentum from the last segment
        auto const& ss = ks.segments().back();  // KalSegment
        double p = ss.mom();
        _hCCEoverP->Fill(ks.caloCluster()->energyDep() / p);
        _hCCDt->Fill(chs._udt);
        _hCCDOCA->Fill(chs._udoca);
        _hCCcdepth->Fill(chs._cdepth);
        _hCCtz->Fill(chs._cpos.Z());
      }
      // Assocated TrkStrawHit info
      for (auto const& tshs : ks.hits()) {
        if (tshs.strawHitState()>WireHitState::inactive) {
          _hHDrift->Fill(tshs.driftRadius());
          _hHDOCA->Fill(tshs.wireDOCA());
          _hHEDep->Fill(1000 * tshs.energyDep());
          _hHPanel->Fill(tshs.strawId().uniquePanel());
        }
      }
      // Assocated Material info
      float radlensum(0.0);
      for (auto const& ts : ks.straws()) {
        if (ts.active()) {
          _hSRadLen->Fill(ts._radlen);
          radlensum += ts._radlen;
        }
      }
      _hSRadLenSum->Fill(radlensum);
    }
    return 0;
  }

  double ValKalSeed::mcTrkP(art::Event const& event, VirtualDetectorId const& vdid,double& p_pri) {
    // first, find the SimParticle associated
    std::vector<art::Handle<PrimaryParticle>> ppv = event.getMany<PrimaryParticle>();
    if(ppv.size()>0){
      auto const& pp = *ppv[0];
      auto psimp = pp.primarySimParticles()[0];
      p_pri = psimp->startMomXYZT().Vect().R();
      // find first KalSeedMC associated with MC primary.
      // This isnt horribly efficient and can fail when there are multiple fits of the
      // same SimParticle, but Validation isn't designed to allow true comparisons
      std::vector<art::Handle<KalSeedMCCollection>> ksmchv = event.getMany<KalSeedMCCollection>();
      for(auto ksmch : ksmchv) {
        for(auto const& ksmc : *ksmch) {
          for(auto const& spp : ksmc.simParticles()){
            if(spp._spkey == psimp->id()){
              // find the correct VDStep
              for(auto const& vdstep : ksmc.vdSteps()){
                if(vdstep._vdid == vdid){
                  return vdstep._mom.R();
                }
              }
            }
          }
        }
      }
    }
    // failure
    return 0.0;
  }
}

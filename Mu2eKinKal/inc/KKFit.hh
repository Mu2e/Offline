#ifndef Mu2eKinKal_KKFit_hh
#define Mu2eKinKal_KKFit_hh
//
// helper class for constructing KinKal Fits using Mu2e data
//
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHitCluster.hh"
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawXing.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawMaterial.hh"
#include "Offline/Mu2eKinKal/inc/KKCaloHit.hh"
#include "Offline/Mu2eKinKal/inc/KKFitUtilities.hh"
#include "Offline/Mu2eKinKal/inc/KKFitSettings.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
// art includes
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Handle.h"
// mu2e Includes
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/DataProducts/inc/StrawIdMask.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/KalSeedAssns.hh"
#include "Offline/RecoDataProducts/inc/KalIntersection.hh"
// geometry
#include "Offline/KinKalGeom/inc/Tracker.hh"
// KinKal includes
#include "KinKal/Fit/Status.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Trajectory/SensorLine.hh"
// Other
#include "cetlib_except/exception.h"
#include <memory>
#include <cmath>
#include <limits>
#include <algorithm>
namespace mu2e {
  using KinKal::SensorLine;
  using KinKal::TimeRange;
  using KinKal::Status;
  using KinKal::BFieldMap;
  using StrawHitIndexCollection = std::vector<StrawHitIndex>;
  using Mu2eKinKal::KKFitConfig;
  using CCHandle = art::Handle<CaloClusterCollection>;
  template <class KTRAJ> class KKFit {
    public:
      // fit configuration
      using KKTRK = KKTrack<KTRAJ>;
      using KKTRKPTR = std::unique_ptr<KKTRK>;
      using PTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using PCA = KinKal::PiecewiseClosestApproach<KTRAJ,SensorLine>;
      using TCA = KinKal::ClosestApproach<KTRAJ,SensorLine>;
      using KKHIT = KinKal::Measurement<KTRAJ>;
      using KKMAT = KinKal::Material<KTRAJ>;
      using KKSTRAWHIT = KKStrawHit<KTRAJ>;
      using KKSTRAWHITPTR = std::shared_ptr<KKSTRAWHIT>;
      using KKSTRAWHITCOL = std::vector<KKSTRAWHITPTR>;
      using KKSTRAWHITCLUSTER = KKStrawHitCluster<KTRAJ>;
      using KKSTRAWHITCLUSTERPTR = std::shared_ptr<KKSTRAWHITCLUSTER>;
      using KKSTRAWHITCLUSTERCOL = std::vector<KKSTRAWHITCLUSTERPTR>;
      using KKSTRAWHITCLUSTERER = KKStrawHitClusterer<KTRAJ>;
      using KKSTRAWXING = KKStrawXing<KTRAJ>;
      using KKSTRAWXINGPTR = std::shared_ptr<KKSTRAWXING>;
      using KKSTRAWXINGCOL = std::vector<KKSTRAWXINGPTR>;
      using KKCALOHIT = KKCaloHit<KTRAJ>;
      using KKCALOHITPTR = std::shared_ptr<KKCALOHIT>;
      using KKCALOHITCOL = std::vector<KKCALOHITPTR>;
      //      using KKPANELHIT = KKPanelHit<KTRAJ>;
      using MEAS = KinKal::Hit<KTRAJ>;
      using MEASPTR = std::shared_ptr<MEAS>;
      using MEASCOL = std::vector<MEASPTR>;
      using EXING = KinKal::ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
      using EXINGCOL = std::vector<EXINGPTR>;
      enum SaveTraj {none=0, full, detector, t0seg};
      // construct from fit configuration objects
      explicit KKFit(KKFitConfig const& fitconfig);
      // helper functions used to create components of the fit
      bool makeStrawHits(Tracker const& tracker,StrawResponse const& strawresponse, BFieldMap const& kkbf, KKStrawMaterial const& smat,
          PTRAJ const& ptraj, ComboHitCollection const& chcol, StrawHitIndexCollection const& strawHitIdxs,
          KKSTRAWHITCOL& hits, KKSTRAWXINGCOL& exings) const;
      SensorLine caloAxis(CaloCluster const& cluster, Calorimeter const& calo) const; // should come from CaloCluster TODO
      bool makeCaloHit(CCPtr const& cluster, Calorimeter const& calo, PTRAJ const& pktraj, KKCALOHITCOL& hits) const;
      // extend a track with a new configuration, optionally searching for and adding hits and straw material
      void extendTrack(Config const& config, BFieldMap const& kkbf, Tracker const& tracker,
          StrawResponse const& strawresponse, KKStrawMaterial const& smat, ComboHitCollection const& chcol,
          Calorimeter const& calo, CCHandle const& cchandle,
          KKTRK& kktrk) const;
      // extend the fit to the surfaces specified in the config
      void extendFit(KKTRK& kktrk);
      // save the complete fit trajectory as a seed
      KalSeed createSeed(KKTRK const& kktrk, TrkFitFlag const& seedflag, Calorimeter const& calo) const;
      TimeRange range(KKSTRAWHITCOL const& strawhits, KKCALOHITCOL const& calohits, KKSTRAWXINGCOL const& strawxings) const; // time range from a set of hits and element Xings
      bool useCalo() const { return usecalo_; }
      bool correctMaterial() const { return matcorr_; }
      bool addMaterial() const { return addmat_; }
      bool addHits() const { return addhits_; }
      auto const& strawHitClusterer() const { return shclusterer_; }
    private:
      void fillTrackerInfo(Tracker const& tracker) const;
      void addStrawHits(Tracker const& tracker,StrawResponse const& strawresponse, BFieldMap const& kkbf, KKStrawMaterial const& smat,
          KKTRK const& kktrk, ComboHitCollection const& chcol, KKSTRAWHITCOL& hits) const;
      void addStraws(Tracker const& tracker, KKStrawMaterial const& smat, KKTRK const& kktrk, KKSTRAWHITCOL const& addhits, KKSTRAWXINGCOL& addexings) const;
      void addCaloHit(Calorimeter const& calo, KKTRK& kktrk, CCHandle cchandle, KKCALOHITCOL& hits) const;
      void sampleFit(KKTRK const& kktrk,KalIntersectionCollection& inters) const; // sample fit at the surfaces specified in the config
      void extendFit(KKTRK& kktrk) const;
      int printLevel_;
      unsigned minNStrawHits_;
      bool matcorr_, addhits_, addmat_, usecalo_; // flags
      KKSTRAWHITCLUSTERER shclusterer_; // functor to cluster KKStrawHits
      // CaloHit configuration
      double caloDt_; // calo time offset; should come from proditions TODO!
      double caloPosRes_; // calo cluster transverse position resolution; should come from proditions or CaloCluster TODO!
      double caloTimeRes_; // calo cluster time resolution; should come from proditions or CaloCluster TODO!
      double caloPropSpeed_; // effective light propagation speed in a crystal (including reflections).  Should come from prodtions TODO
      double minCaloEnergy_; // minimum CaloCluster energy
      double maxCaloDt_; // maximum track-calo time difference
      double maxCaloDoca_; // maximum track-calo DOCA
      // straw hit creation and processing
      double tprec_; // TPOCA calculation nominal precision
      StrawHitFlag addsel_, addrej_; // selection and rejection flags when adding hits
      // parameters controlling adding hits
      float maxStrawHitDoca_, maxStrawHitDt_, maxStrawDoca_, maxStrawDocaCon_, maxStrawUposBuff_;
      int maxDStraw_; // maximum distance from the track a strawhit can be to consider it for adding.
      // cached info computed from the tracker, used in hit adding; these must be lazy-evaluated as the tracker doesn't exist on construction
      mutable double strawradius_;
      mutable double ymin_, ymax_, umax_; // panel-level info
      mutable double rmin_, rmax_; // plane-level info
      mutable double spitch_;
      mutable bool needstrackerinfo_ = true;

      SaveTraj savetraj_; // trajectory saving option
  };

  template <class KTRAJ> KKFit<KTRAJ>::KKFit(KKFitConfig const& fitconfig) :
    printLevel_(fitconfig.printLevel()),
    minNStrawHits_(fitconfig.minNStrawHits()),
    matcorr_(fitconfig.matCorr()),
    addhits_(fitconfig.addHits()),
    addmat_(fitconfig.addMaterial()),
    usecalo_(fitconfig.useCaloCluster()),
    shclusterer_(StrawIdMask(fitconfig.strawHitClusterLevel()),fitconfig.strawHitClusterDeltaStraw(),fitconfig.strawHitClusterDeltaT()),
    caloDt_(fitconfig.caloDt()),
    caloPosRes_(fitconfig.caloPosRes()),
    caloTimeRes_(fitconfig.caloTimeRes()),
    caloPropSpeed_(fitconfig.caloPropSpeed()),
    minCaloEnergy_(fitconfig.minCaloEnergy()),
    maxCaloDt_(fitconfig.maxCaloDt()),
    maxCaloDoca_(fitconfig.maxCaloDoca()),
    tprec_(fitconfig.tpocaPrec()),
    addsel_(fitconfig.addHitSelect()),
    addrej_(fitconfig.addHitReject()),
    maxStrawHitDoca_(fitconfig.maxStrawHitDOCA()),
    maxStrawHitDt_(fitconfig.maxStrawHitDt()),
    maxStrawDoca_(fitconfig.maxStrawDOCA()),
    maxStrawDocaCon_(fitconfig.maxStrawDOCAConsistency()),
    maxStrawUposBuff_(fitconfig.maxStrawUposBuff()),
    maxDStraw_(fitconfig.maxDStraw())
  {
    if (fitconfig.saveTraj() == "T0") {
        savetraj_ = t0seg;
    } else if (fitconfig.saveTraj() == "Full") {
        savetraj_ = full;
    } else if (fitconfig.saveTraj() == "Detector") {
        savetraj_ = detector;
    } else if (fitconfig.saveTraj() == "None") {
        savetraj_ = none;
    } else {
      throw cet::exception("RECO")<<"mu2e::KKFit: unknown trajectory option "<< fitconfig.saveTraj() << endl;
    }
  }

  template <class KTRAJ> bool KKFit<KTRAJ>::makeStrawHits(Tracker const& tracker,StrawResponse const& strawresponse,BFieldMap const& kkbf, KKStrawMaterial const& smat,
      PTRAJ const& ptraj, ComboHitCollection const& chcol, StrawHitIndexCollection const& strawHitIdxs,
      KKSTRAWHITCOL& hits, KKSTRAWXINGCOL& exings) const {
    unsigned ngood(0);
    // loop over the individual straw combo hits
    for(auto strawidx : strawHitIdxs) {
      const ComboHit& combohit(chcol.at(strawidx));
      if(combohit.mask().level() != StrawIdMask::uniquestraw){
        throw cet::exception("RECO")<<"mu2e::KKFit: ComboHit error"<< endl;
      }
      const Straw& straw = tracker.getStraw(combohit.strawId());
      auto wline = Mu2eKinKal::hitLine(combohit,straw,strawresponse); // points from the signal to the straw center
      double htime = combohit.correctedTime();
      // use the hit z position to estimate the particle time.
      auto hpos = wline.position3(htime);
      auto ppos = ptraj.position3(htime);
      auto vel = ptraj.velocity(htime);
      double ptime = htime;
      if (vel.Z() != 0)
        ptime += (hpos.Z()-ppos.Z())/vel.Z();
      CAHint hint(ptime,htime);
      // compute PCA between the seed trajectory and this straw
      PCA pca(ptraj, wline, hint, tprec_ );
      // create the hit.  Note these may initially be unusable
      hits.push_back(std::make_shared<KKSTRAWHIT>(kkbf, pca, combohit, straw, strawidx, strawresponse));
      // create the material crossing, including this reference
      if(matcorr_) exings.push_back(std::make_shared<KKSTRAWXING>(hits.back(),smat));
      if(hits.back()->hitState().usable())ngood++;
    }
    return ngood >= minNStrawHits_;
  }

  template <class KTRAJ> SensorLine KKFit<KTRAJ>::caloAxis(CaloCluster const& cluster, Calorimeter const& calo) const {
    // move cluster COG into the tracker frame.  COG is at the front face of the disk
    CLHEP::Hep3Vector cog = calo.geomUtil().mu2eToTracker(calo.geomUtil().diskFFToMu2e( cluster.diskID(), cluster.cog3Vector()));
    // project this along the crystal axis to the SIPM, which is at the back.  This is the point the time measurement corresponds to
    VEC3 ffcog(cog);
    double lcrystal = calo.caloInfo().getDouble("crystalZLength"); // text-keyed lookup is very inefficient FIXME!
    VEC3 crystalF2B = VEC3(0.0,0.0,lcrystal); // this should come directly from the calogeometry, TODO
    VEC3 sipmcog = ffcog + crystalF2B;
    // create the SensorLine trajectory from this information: signal goes towards the sipm
    return SensorLine(sipmcog,ffcog,cluster.time()+caloDt_,caloPropSpeed_);
  }


  template <class KTRAJ> bool KKFit<KTRAJ>::makeCaloHit(CCPtr const& cluster, Calorimeter const& calo, PTRAJ const& pktraj, KKCALOHITCOL& hits) const {
    bool retval(false);
    auto caxis = caloAxis(*cluster,calo);
    // find the time the seed traj passes the middle of the crystal to form the hint
    auto pmid = caxis.position3(caxis.timeAtMidpoint());
    double zt = Mu2eKinKal::zTime(pktraj,pmid.Z(),pktraj.range().end());
    CAHint hint( zt, caxis.measurementTime());
    // compute a preliminary PCA between the seed trajectory and the cluster axis
    auto pca = PCA(pktraj, caxis, hint, tprec_ );
    if(pca.usable() && fabs(pca.doca()) < maxCaloDoca_ && fabs(pca.deltaT()) < maxCaloDt_){
      // check that CA is within the active volume of the calorimeter
      double dz = pca.sensorPoca().Z() - caxis.position3(caxis.measurementTime()).Z();
      if( dz > -caxis.length() -maxCaloDoca_ && dz < maxCaloDoca_) {
//        double tvar = cluster->timeErr()*cluster->timeErr(); the returned value seems unphysically smalL, ~70 ps.
        double tvar = caloTimeRes_*caloTimeRes_; // temporary kludge, this number comes from MDC2020 sim studies.  FIXME
        double wvar = caloPosRes_*caloPosRes_;
        hits.push_back(std::make_shared<KKCALOHIT>(cluster,pca,tvar,wvar));
        retval = true;
      }
    }
    return retval;
  }

  template <class KTRAJ> void KKFit<KTRAJ>::extendTrack(Config const& exconfig, BFieldMap const& kkbf, Tracker const& tracker,
      StrawResponse const& strawresponse, KKStrawMaterial const& smat, ComboHitCollection const& chcol,
      Calorimeter const& calo, CCHandle const& cchandle,
      KKTRK& kktrk) const {
    KKSTRAWHITCOL addstrawhits;
    KKCALOHITCOL addcalohits;
    KKSTRAWXINGCOL addstrawxings;
    if(addhits_)addStrawHits(tracker, strawresponse, kkbf, smat, kktrk, chcol, addstrawhits );
    if(matcorr_ && addmat_)addStraws(tracker, smat, kktrk, addstrawhits, addstrawxings);
    if(addhits_ && usecalo_ && kktrk.caloHits().size()==0)addCaloHit(calo, kktrk, cchandle, addcalohits);
    if(printLevel_ > 1){
      std::cout << "KKTrk extension adding "
        << addstrawhits.size() << " StrawHits and "
        << addcalohits.size() << " CaloHits and "
        << addstrawxings.size() << " Straw Xings" << std::endl;
    }
    kktrk.extendTrack(exconfig,addstrawhits,addstrawxings,addcalohits);
  }

  template <class KTRAJ> void KKFit<KTRAJ>::addStrawHits(Tracker const& tracker,StrawResponse const& strawresponse, BFieldMap const& kkbf, KKStrawMaterial const& smat,
      KKTRK const& kktrk, ComboHitCollection const& chcol, KKSTRAWHITCOL& addhits) const {
    auto const& ftraj = kktrk.fitTraj();
    // build the set of existing hits
    std::set<StrawHitIndex> oldhits;
    for(auto const& strawhit : kktrk.strawHits())oldhits.insert(strawhit->strawHitIndex());
    for( size_t ich=0; ich < chcol.size();++ich){
      if(oldhits.find(ich)==oldhits.end()){      // make sure this hit wasn't already found
        ComboHit const& strawhit = chcol[ich];
        if(strawhit.flag().hasAllProperties(addsel_) && (!strawhit.flag().hasAnyProperty(addrej_))){
          double zt = Mu2eKinKal::zTime(ftraj,strawhit.pos().Z(),strawhit.correctedTime());
          if(fabs(strawhit.correctedTime()-zt) < maxStrawHitDt_) {      // compare the measured time with the estimate from the fit
            const Straw& straw = tracker.getStraw(strawhit.strawId());
            auto wline = Mu2eKinKal::hitLine(strawhit,straw,strawresponse);
            double psign = wline.direction().Dot(straw.wireDirection());  // wire distance is WRT straw center, in the nominal wire direction
            double htime = wline.measurementTime() - (straw.halfLength()-psign*strawhit.wireDist())/wline.speed(wline.timeAtMidpoint());
            CAHint hint(zt,htime);
            // compute PCA between the trajectory and this straw
            PCA pca(ftraj, wline, hint, tprec_ );
            if(fabs(pca.doca()) < maxStrawHitDoca_){ // add test of chi TODO
              addhits.push_back(std::make_shared<KKSTRAWHIT>(kkbf, pca, strawhit, straw, ich, strawresponse));
            }
          }
        }
      }
    }
  }

  template <class KTRAJ> void KKFit<KTRAJ>::addStraws(Tracker const& tracker, KKStrawMaterial const& smat, KKTRK const& kktrk,
      KKSTRAWHITCOL const& addhits, KKSTRAWXINGCOL& addexings) const {
    // this algorithm assumes the track never hits the same straw twice.  That could be violated by reflecting tracks, and could be addressed
    // by including the time of the Xing as part of its identity.  That would slow things down so it remains to be proven it's a problem  TODO
    // build the set of existing straws
    auto const& ftraj = kktrk.fitTraj();
    // pre-compute some tracker info if needed
    if(needstrackerinfo_)fillTrackerInfo(tracker);
    // list the IDs of existing straws: this speeds the search
    std::set<StrawId> oldstraws;
    for(auto const& strawxing : kktrk.strawXings())oldstraws.insert(strawxing->strawId());
    // create materials for straw hits just added
    for(auto const& strawhit : addhits){
      if(oldstraws.find(strawhit->strawId()) == oldstraws.end()){
        addexings.push_back(std::make_shared<KKSTRAWXING>(strawhit,smat));
        oldstraws.insert(strawhit->strawId());
      }
    }
    // Go hierarchically through planes and panels to find new straws hit by this track.
    for(auto const& plane : tracker.planes()){
      if(tracker.planeExists(plane.id())) {
        double plz = plane.origin().z();
        // find the track position in at this plane's central z.
        double zt = Mu2eKinKal::zTime(ftraj,plz,ftraj.range().begin());
        auto plpos = ftraj.position3(zt);
        // rough check on the point radius
        double rho = plpos.Rho();
        if(rho > rmin_ && rho < rmax_){
          // loop over panels in this plane
          for(auto panel_p : plane.panels()){
            auto const& panel = *panel_p;
            // linearly correct position for the track direction due to difference in panel-plane Z position
            auto tdir = ftraj.direction(zt);
            double dz = panel.origin().z()-plz;
            auto papos = plpos + (dz/tdir.Z())*tdir;
            // convert this position into panel coordinates
            CLHEP::Hep3Vector cpos(papos.X(),papos.Y(),papos.Z()); // clumsy translation
            auto pposv = panel.dsToPanel()*cpos;
            // translate the y position into a rough straw number
            int istraw = static_cast<int>(rint( (pposv.y()-ymin_)*spitch_));
            // require this be within the (integral) straw buffer.  This just reduces the number of calls to PCA
            if(istraw >= -maxDStraw_ && istraw < static_cast<int>(panel.nStraws()) + maxDStraw_ ){
              unsigned istrmin = static_cast<unsigned>(std::max(istraw-maxDStraw_,0));
              // largest straw is the innermost; use that to test length
              if(fabs(pposv.x()) < panel.getStraw(istrmin).halfLength() + maxStrawUposBuff_ ) {
                unsigned istrmax = static_cast<unsigned>(std::min(istraw+maxDStraw_,static_cast<int>(panel.nStraws())-1));
                // loop over straws
                for(unsigned istr = istrmin; istr <= istrmax; ++istr){
                  auto const& straw = panel.getStraw(istr);
                  // add strawExists test TODO
                  // make sure we haven't already seen this straw
                  if(oldstraws.find(straw.id()) == oldstraws.end()){
                    auto sline = Mu2eKinKal::strawLine(straw,zt); // line down the straw axis center
                    CAHint hint(zt,zt);
                    // compute PCA between the trajectory and this straw
                    PCA pca(ftraj, sline, hint, tprec_ );
                    // require consistency with this track passing through this straw
                    double du = fabs((pca.sensorPoca().Vect()-VEC3(straw.wirePosition(0.0))).Dot(VEC3(straw.wireDirection(0.0))));
                    double doca = fabs(pca.doca());
                    double dsig = std::max(0.0,doca-strawradius_)/sqrt(pca.docaVar());
                    if(doca < maxStrawDoca_ && dsig < maxStrawDocaCon_ && du < straw.halfLength() + maxStrawUposBuff_){
                      addexings.push_back(std::make_shared<KKSTRAWXING>(pca.localClosestApproach(),smat,straw));
                      oldstraws.insert(straw.id());
                    }
                  } // not existing straw cut
                } // straws loop
              } // straw length cut
            } // straw index cut
          } // panels loop
        } // radius cut
      } // plane exists cut
    }  // planes loop
  } // end function

  template <class KTRAJ> void KKFit<KTRAJ>::addCaloHit(Calorimeter const& calo, KKTRK& kktrk, CCHandle cchandle, KKCALOHITCOL& hits) const {
    double crystalLength = calo.caloInfo().getDouble("crystalZLength");
    auto const& ftraj = kktrk.fitTraj();
    auto cccol = cchandle.product();
    double edep(-1.0);
    std::shared_ptr<KKCALOHIT> chitptr;
    // loop over disks to decide which are worth testing
    std::array<bool,2> test{false,false};
    for(unsigned idisk=0; idisk < 2; ++idisk){
      auto ffpos = calo.geomUtil().mu2eToTracker(calo.disk(idisk).geomInfo().frontFaceCenter());
      double rmin = calo.disk(idisk).geomInfo().innerEnvelopeR() - maxCaloDoca_;
      double rmax = calo.disk(idisk).geomInfo().outerEnvelopeR() + maxCaloDoca_;
      // test at both faces; if the track is in the right area, test the clusters on this disk
      // Replace this with an intersection with the calo face TODO
      for(int iface=0; iface<2; ++iface){
        double zt = Mu2eKinKal::zTime(ftraj,ffpos.z()+iface*crystalLength,ftraj.range().end());
        auto tpos = ftraj.position3(zt);
        double rho = tpos.Rho();
        test[idisk] |= rho > rmin && rho < rmax;
      }
    }
    // now loop over crystals and find the best match
    for(size_t icc=0;icc < cccol->size(); ++icc){
      auto const& cc = (*cccol)[icc];
      auto idisk = static_cast<size_t>(cc.diskID());
      if (test[idisk] && cc.energyDep() > minCaloEnergy_ && cc.energyDep() > edep){
        // create PCA from this cluster and the traj
        auto caxis = caloAxis(cc,calo);
        // find the time the seed traj passes the middle of the crystal to form the hint
        auto pmid = caxis.position3(caxis.timeAtMidpoint());
        double zt = Mu2eKinKal::zTime(ftraj,pmid.Z(),ftraj.range().end());
        CAHint hint( zt, caxis.timeAtMidpoint());
        // compute closest approach between the fit trajectory and the cluster axis
        auto pca = PCA(ftraj, caxis, hint, tprec_ );
        if(pca.usable() && fabs(pca.doca()) < maxCaloDoca_ && fabs(pca.deltaT()) < maxCaloDt_){
          // check that the position is within the active position of the crystal
          double dz = pca.sensorPoca().Z() - caxis.position3(caxis.measurementTime()).Z();
          if( dz > -caxis.length() -maxCaloDoca_ && dz < maxCaloDoca_) {
            art::Ptr<CaloCluster> ccPtr = art::Ptr<CaloCluster>(cchandle,icc);
//            double tvar = cc.timeErr()*cc.timeErr();
            double tvar = caloTimeRes_*caloTimeRes_;
            double wvar = caloPosRes_*caloPosRes_;
            chitptr = std::make_shared<KKCALOHIT>(ccPtr,pca,tvar,wvar);
            edep = cc.energyDep();
          }
        }
      }
    }
    // Add the best match (if there is one)
    if(!(chitptr == nullptr))hits.push_back(chitptr);
  }

  template <class KTRAJ> void KKFit<KTRAJ>::fillTrackerInfo(Tracker const& tracker) const {
    // pre-compute some info derived from the tracker
    strawradius_ = tracker.strawOuterRadius();
    auto const& frontplane = tracker.planes().front();
    auto const& firstpanel = frontplane.getPanel(0);
    auto const& innerstraw = firstpanel.getStraw(0);
    auto const& outerstraw = firstpanel.getStraw(StrawId::_nstraws-1);
    // compute limits: add some buffer for the finite size of the straw
    auto DStoP = firstpanel.dsToPanel();
    auto innerstraw_origin = DStoP*innerstraw.origin();
    auto outerstraw_origin = DStoP*outerstraw.origin();
    ymin_ = innerstraw_origin.y();
    ymax_ = outerstraw_origin.y();
    umax_ = innerstraw.halfLength() + strawradius_; // longest possible straw
    // plane-level variables: these add some buffer
    rmin_ = innerstraw_origin.y() - maxDStraw_*strawradius_;
    rmax_ = outerstraw.wireEnd(StrawEnd::cal).mag() + maxDStraw_*strawradius_;
    spitch_ = (StrawId::_nstraws-1)/(ymax_-ymin_);
    needstrackerinfo_= false;
  }


  template <class KTRAJ> TimeRange KKFit<KTRAJ>::range(KKSTRAWHITCOL const& strawhits, KKCALOHITCOL const& calohits, KKSTRAWXINGCOL const& strawxings) const{
    double tmin = std::numeric_limits<float>::max();
    double tmax = -tmin;
    for( auto const& strawhit : strawhits) {
      tmin = std::min(tmin,strawhit->time());
      tmax = std::max(tmax,strawhit->time());
    }
    for( auto const& calohit : calohits) {
      tmin = std::min(tmin,calohit->time());
      tmax = std::max(tmax,calohit->time());
    }
    for( auto const& strawxing : strawxings) {
      tmin = std::min(tmin,strawxing->time());
      tmax = std::max(tmax,strawxing->time());
    }
    return TimeRange(tmin,tmax);
  }

   template <class KTRAJ> void KKFit<KTRAJ>::extendFit(KKTRK& kktrk) {
     // extend the fit upstream and downstream (up and down for cosmics) to the specified surfaces;
 //TODO

   }

  template <class KTRAJ> KalSeed KKFit<KTRAJ>::createSeed(KKTRK const& kktrk, TrkFitFlag const& seedflag, Calorimeter const& calo) const {
    TrkFitFlag fflag(seedflag);  // initialize the flag with the seed fit flag
    if(kktrk.fitStatus().usable()){
      fflag.merge(TrkFitFlag::kalmanOK);
      fflag.merge(TrkFitFlag::seedOK);
    }
    if(kktrk.fitStatus().status_ == Status::converged) fflag.merge(TrkFitFlag::kalmanConverged);
    if(matcorr_)fflag.merge(TrkFitFlag::MatCorr);
    if(kktrk.config().bfcorr_ )fflag.merge(TrkFitFlag::BFCorr);
    // explicit T0 is needed for backwards-compatibility; sample from the appropriate trajectory piece
    auto const& fittraj = kktrk.fitTraj();
    double t0val = fittraj.t0();
    auto const& t0piece = fittraj.nearestPiece(t0val);
    double t0sig = sqrt(t0piece.params().covariance()(KTRAJ::t0_,KTRAJ::t0_));
    HitT0 t0(t0val,t0sig);
    // create the shell for the output
    KalSeed kseed(kktrk.fitParticle(),fflag);
    auto const& fstatus = kktrk.fitStatus();
    kseed._chisq = fstatus.chisq_.chisq();
    kseed._ndof = fstatus.chisq_.nDOF();
    kseed._fitcon = fstatus.chisq_.probability();
    size_t igap;
    double maxgap,avggap;
    fittraj.gaps(maxgap,igap,avggap);
    kseed._maxgap = maxgap;
    kseed._avggap = avggap;
    // loop over track components and store them
    kseed._hits.reserve(kktrk.strawHits().size());
    for(auto const& strawhit : kktrk.strawHits() ) {
      Residual utres = strawhit->refResidual(Mu2eKinKal::tresid);
      Residual udres = strawhit->refResidual(Mu2eKinKal::dresid);
      // compute unbiased residuals; this can fail if the track has marginal coverage
      if(kktrk.fitStatus().usable()) {
        try {
          udres = strawhit->residual(Mu2eKinKal::dresid);
          utres = strawhit->residual(Mu2eKinKal::tresid);
        } catch (std::exception const& error) {
         // std::cout << "Unbiased KKStrawHit residual calculation failure, nDOF = " << fstatus.chisq_.nDOF() << std::endl;
        }
      }
      kseed._hits.emplace_back(strawhit->strawHitIndex(),strawhit->hit(),
          strawhit->closestApproach().tpData(),
          strawhit->unbiasedClosestApproach().tpData(),
          utres, udres,
          strawhit->refResidual(Mu2eKinKal::tresid),
          strawhit->refResidual(Mu2eKinKal::dresid),
          strawhit->fillDriftInfo(strawhit->closestApproach()),
          strawhit->hitState(),
          strawhit->straw());
    }
    if(kktrk.caloHits().size() > 0){
      auto const& calohit = kktrk.caloHits().front(); // for now take the front: not sure if there will ever be >1 TODO
      auto const& ca = calohit->closestApproach();
      StrawHitFlag hflag;
      if(calohit->active()){
        hflag.merge(StrawHitFlag::active);
        hflag.merge(StrawHitFlag::doca);
      }
      // calculate the unbiased time residual
      Residual ctres = calohit->refResidual(0);
      if(kktrk.fitStatus().usable()){
        try {
          ctres = calohit->residual(0);
        } catch (std::exception const& error) {
          // std::cout << "Unbiased KKCaloHit residual calculation failure, nDOF = " << fstatus.chisq_.nDOF() << std::endl;
        }
      }
      // calculate the cluster depth = distance along the crystal axis from the POCA to the back face of this disk (where the SiPM sits)
      double backz = calo.geomUtil().mu2eToTracker(calo.disk(calohit->caloCluster()->diskID()).geomInfo().backFaceCenter()).z();
      // calculate the distance from POCA to the SiPM, along the crystal (Z) direction, and projected along the track
      float clen = backz-ca.sensorPoca().Z();
      float trklen = clen/ca.particleTraj().direction(ca.particleToca()).Z();
      kseed._chit = TrkCaloHitSeed(calohit->caloCluster(),hflag,
          clen,trklen,
          ca.tpData(),
          calohit->unbiasedClosestApproach().tpData(),
          ctres,ca.particleTraj().momentum3(ca.particleToca()));
    }
    kseed._straws.reserve(kktrk.strawXings().size());
    for(auto const& sxing : kktrk.strawXings()) {
      // compute energy loss
      double dm, paramomvar, perpmomvar;
      sxing->materialEffects(dm, paramomvar,perpmomvar);
      double radfrac = sxing->radiationFraction();
      // find the gas path
      double gaspath(0.0);
      for(auto const& matxing : sxing->matXings()){
        if(matxing.dmat_.name().compare(5,3,"gas",3)==0)gaspath = matxing.plen_;
      }
      auto const& tpocad = sxing->closestApproach();
      kseed._straws.emplace_back(sxing->strawId(),
          tpocad.doca(),
          tpocad.particleToca(),
          tpocad.sensorToca(),
          gaspath,
          radfrac,
          dm,
          sxing->active() );
    }
    // save the fit segments as requested
    if(kktrk.fitStatus().usable()){
      if (savetraj_ == full){
        kseed._segments.reserve(fittraj.pieces().size());
        for (auto const& traj : fittraj.pieces() ){
          // skip zero-range segments.  By convention, sample the state at the mid-time
          if(traj->range().range() > 0.0) kseed._segments.emplace_back(*traj,traj->range().mid());
        }
      } else if (savetraj_ == detector ) {
        // only save segments inside the tracker volume. First, find the time limits for that
        double tmin = std::numeric_limits<float>::max();
        double tmax = -tmin;
        static const SurfaceId tt_front("TT_Front");
        static const SurfaceId tt_mid("TT_Mid");
        static const SurfaceId tt_back("TT_Back");
        for(auto const& interpair : kktrk.intersections()) {
          auto const& sid = std::get<0>(interpair);
          if(sid == tt_front || sid == tt_mid || sid == tt_back){
            auto const& inter = std::get<1>(interpair);
            tmin = std::min(tmin,inter.time_);
            tmax = std::max(tmax,inter.time_);
          }
        }
        // extend as needed to the calohit. Eventually we should also extend to the calorimeter, but that needs to be an extrapolation
        if(kktrk.caloHits().size() > 0){
          auto const& calohit = kktrk.caloHits().front();
          if(calohit->active()){
            tmin = std::min(tmin,calohit->time());
            tmax = std::max(tmax,calohit->time());
          }
        }
        if(tmin > tmax){
          std::cout << "mu2e::KKFit: tracker intersections missing"<< std::endl;
          tmin = fittraj.range().begin();
          tmax = fittraj.range().end();
        }
        kseed._segments.reserve(fittraj.pieces().size());// this will be oversized
        for (auto const& traj : fittraj.pieces() ){
          // skip segments outside the tracker volume range
          if(traj->range().range() > 0.0 && (traj->range().inRange(tmin) || traj->range().inRange(tmax) || (traj->range().begin() > tmin && traj->range().end() < tmax)) ) kseed._segments.emplace_back(*traj,traj->range().mid());
        }
      } else if (savetraj_ == t0seg ) {
        kseed._segments.emplace_back(t0piece,t0val);
      }
    } else {
      kseed._segments.emplace_back(t0piece,t0val); // save the t0 piece even for failed fits
    }
    // sample the fit at the locations recorded in the fit
    sampleFit(kktrk,kseed._inters);
    // remove unused storage
    kseed._segments.shrink_to_fit();
    return kseed;
  }

  template <class KTRAJ> void KKFit<KTRAJ>::sampleFit(KKTRK const& kktrk,KalIntersectionCollection& inters) const {
    // add IPA and ST Xings. Sample just before the transit to include the effect of the material
    static const double epsilon(1.0e-6);
    auto const& ftraj = kktrk.fitTraj();
    for(auto const& ipaxing : kktrk.IPAXings()){
      double stime = ipaxing->time() - epsilon;
      auto const& ktraj = ftraj.nearestPiece(stime);
      double dmom,paramomvar,perpmomvar;
      ipaxing->materialEffects(dmom,paramomvar,perpmomvar);
      inters.emplace_back(ktraj.stateEstimate(ipaxing->time()),XYZVectorF(ktraj.bnom()),ipaxing->surfaceId(),ipaxing->intersection(),dmom);
    }
    for(auto const& stxing : kktrk.STXings()){
      double stime = stxing->time() - epsilon;
      auto const& ktraj = ftraj.nearestPiece(stime);
      double dmom,paramomvar,perpmomvar;
      stxing->materialEffects(dmom,paramomvar,perpmomvar);
      inters.emplace_back(ktraj.stateEstimate(stxing->time()),XYZVectorF(ktraj.bnom()),stxing->surfaceId(),stxing->intersection(),dmom);
    }
    for(auto const& crvxing : kktrk.CRVXings()){
      double stime = crvxing->time() - epsilon;
      auto const& ktraj = ftraj.nearestPiece(stime);
      double dmom,paramomvar,perpmomvar;
      crvxing->materialEffects(dmom,paramomvar,perpmomvar);
      inters.emplace_back(ktraj.stateEstimate(crvxing->time()),XYZVectorF(ktraj.bnom()),crvxing->surfaceId(),crvxing->intersection(),dmom);
    }
    // record other intersections saved in the track
    for(auto const& interpair : kktrk.intersections()) {
      auto const& sid = std::get<0>(interpair);
      auto const& inter = std::get<1>(interpair);
      auto const& ktraj = ftraj.nearestPiece(inter.time_);
      inters.emplace_back(ktraj.stateEstimate(inter.time_),XYZVectorF(ktraj.bnom()),sid,inter);
    }
    // sort by time TODO
  }
}
#endif

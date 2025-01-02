#ifndef Mu2eKinKal_KKTrack_hh
#define Mu2eKinKal_KKTrack_hh
//
// subclass of KinKal Track specialized for Mu2e
//
#include "KinKal/Fit/Track.hh"
#include "KinKal/Detector/ParameterHit.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHitCluster.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawXing.hh"
#include "Offline/Mu2eKinKal/inc/KKShellXing.hh"
#include "Offline/Mu2eKinKal/inc/KKCaloHit.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "KinKal/Geometry/Intersection.hh"
#include <tuple>
namespace mu2e {

  using KinKal::Config;
  using KinKal::BFieldMap;
  using KinKal::TimeDir;
  using KinKal::Intersection;

  template <class KTRAJ> class KKTrack : public KinKal::Track<KTRAJ> {
    public:
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
      using MEAS = KinKal::Hit<KTRAJ>;
      using MEASPTR = std::shared_ptr<MEAS>;
      using MEASCOL = std::vector<MEASPTR>;
      using EXING = KinKal::ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
      using EXINGCOL = std::vector<EXINGPTR>;
      using KKIPAXING = KKShellXing<KTRAJ,KinKal::Cylinder>;
      using KKIPAXINGPTR = std::shared_ptr<KKIPAXING>;
      using KKIPAXINGCOL = std::vector<KKIPAXINGPTR>;
      using KKSTXING = KKShellXing<KTRAJ,KinKal::Annulus>;
      using KKSTXINGPTR = std::shared_ptr<KKSTXING>;
      using KKSTXINGCOL = std::vector<KKSTXINGPTR>;
      using KKCRVXING = KKShellXing<KTRAJ,KinKal::Rectangle>;
      using KKCRVXINGPTR = std::shared_ptr<KKCRVXING>;
      using KKCRVXINGCOL = std::vector<KKCRVXINGPTR>;
      using KKINTER = std::tuple<SurfaceId,KinKal::Intersection>;
      using KKINTERCOL = std::vector<KKINTER>;
      using TRACK = KinKal::Track<KTRAJ>;
      // construct from configuration, fit environment, and hits and materials
      KKTrack(Config const& config, BFieldMap const& bfield, KTRAJ const& seedtraj, PDGCode::type tpart, KKSTRAWHITCLUSTERER const& shclusterer,
          KKSTRAWHITCOL const& strawhits, KKSTRAWXINGCOL const& strawxings, KKCALOHITCOL const& calohits, std::array<double, KinKal::NParams()> constraints = {0});
      // extend the track according to new configuration, hits, and/or exings
      void extendTrack(Config const& config,
          KKSTRAWHITCOL const& strawhits, KKSTRAWXINGCOL const& strawxings, KKCALOHITCOL const& calohits );
      // extend the track to cover a new set of material Xings.  This will reuse the existing config object
      void extendTrack(EXINGCOL const& xings);
      // add IPA Xing
      void addIPAXing(KKIPAXINGPTR const& ipaxing,TimeDir const& tdir);
      // add ST Xing
      void addSTXing(KKSTXINGPTR const& stxing,TimeDir const& tdir);
      // add TCRV Xing
      void addTCRVXing(KKCRVXINGPTR const& crvxing,TimeDir const& tdir);
      // add intersections
      void addIntersection(SurfaceId const& sid, Intersection const& inter) { inters_.emplace_back(sid,inter); }

      // accessors
      PDGCode::type fitParticle() const { return tpart_;}
      KKSTRAWHITCOL const& strawHits() const { return strawhits_; }
      KKSTRAWHITCLUSTERCOL const& strawHitClusters() const { return strawhitclusters_; }
      KKSTRAWXINGCOL const& strawXings() const { return strawxings_; }
      KKIPAXINGCOL const& IPAXings() const { return ipaxings_; }
      KKSTXINGCOL const& STXings() const { return stxings_; }
      KKCRVXINGCOL const& CRVXings() const { return crvxings_; }
      KKINTERCOL const& intersections() const { return inters_; }
      KKCALOHITCOL const& caloHits() const { return calohits_; }
      void printFit(std::ostream& ost=std::cout,int detail=0) const;
    private:
      // record the particle type
      PDGCode::type tpart_;
      KKSTRAWHITCLUSTERER shclusterer_; // clustering functor
      KKSTRAWHITCOL strawhits_;  // straw hits used in this fit
      KKSTRAWXINGCOL strawxings_;  // straw material crossings used in this fit
      KKIPAXINGCOL ipaxings_;  // ipa material crossings used in extrapolation
      KKSTXINGCOL stxings_;  // stopping target material crossings used in extrapolation
      KKCRVXINGCOL crvxings_; // crv crossings using in extrapolation
      KKINTERCOL inters_; // other recorded intersections
      KKSTRAWHITCLUSTERCOL strawhitclusters_;  // straw hit clusters used in this fit
      KKCALOHITCOL calohits_;  // calo hits used in this fit
      // utility function to convert to generic types
      void convertTypes( KKSTRAWHITCOL const& strawhits, KKSTRAWXINGCOL const& strawxings,KKCALOHITCOL const& calohits,
          MEASCOL& hits, EXINGCOL& exings);
      // add hits to clusters
      void addHitClusters(KKSTRAWHITCOL const& strawhits,KKSTRAWXINGCOL const& strawxings, MEASCOL& hits);
  };

  template <class KTRAJ> KKTrack<KTRAJ>::KKTrack(Config const& config, BFieldMap const& bfield, KTRAJ const& seedtraj, PDGCode::type tpart,
      KKSTRAWHITCLUSTERER const& shclusterer,
      KKSTRAWHITCOL const& strawhits,
      KKSTRAWXINGCOL const& strawxings,
      KKCALOHITCOL const& calohits,
      std::array<double, KinKal::NParams()> constraints) :
    KinKal::Track<KTRAJ>(config,bfield,seedtraj), tpart_(tpart), shclusterer_(shclusterer),
    strawhits_(strawhits),
    strawxings_(strawxings),
    calohits_(calohits) {
      MEASCOL hits; // polymorphic container of hits
      EXINGCOL exings; // polymorphic container of detector element crossings
      // add the hits to clusters, as required
      addHitClusters(strawhits_,strawxings_,hits);
      if(this->config().plevel_ > 0){
        std::cout << "created " << strawhitclusters_.size() << " StrawHitClusters " << std::endl;
        for (auto const& shcluster : strawhitclusters_) {
          shcluster->print(std::cout,this->config().plevel_);
        }
      }
      convertTypes(strawhits_, strawxings_, calohits_,  hits,exings);

      std::array<bool,KinKal::NParams()> mask = {false};
      bool constraining = false;
      for (size_t i=0;i<KinKal::NParams();i++){
        if (constraints[i] > 0){
          mask[i] = true;
          constraining = true;
        }
      }
      if (constraining){
        KinKal::Parameters cparams = seedtraj.params();
        for (size_t ipar=0;ipar<KinKal::NParams();ipar++){
          for (size_t jpar=0;jpar<KinKal::NParams();jpar++){
            cparams.covariance()[ipar][jpar] = 0.0;
          }
        }
        for(size_t ipar=0; ipar < KinKal::NParams(); ipar++){
          if (mask[ipar])
            cparams.covariance()[ipar][ipar] = constraints[ipar]*constraints[ipar];
          else
            cparams.covariance()[ipar][ipar] = 1.0; // otherwise inversion fails
        }
        hits.push_back(std::make_shared<KinKal::ParameterHit<KTRAJ>>(seedtraj.range().mid(),seedtraj,cparams,mask));
      }

      this->fit(hits,exings);
    }

  template <class KTRAJ> void KKTrack<KTRAJ>::addHitClusters(KKSTRAWHITCOL const& strawhits,KKSTRAWXINGCOL const& strawxings,MEASCOL& hits) {
    if(shclusterer_.clusterLevel() != StrawIdMask::none){
      for(auto const& strawhitptr : strawhits){
        bool added(false);
        for(auto& shclusterptr : strawhitclusters_) {
          if(shclusterptr->canAddHit(strawhitptr,shclusterer_)){
            shclusterptr->addHit(strawhitptr,shclusterer_);
            added = true;
            break;
          }
        }
        if(!added){
          strawhitclusters_.emplace_back(std::make_shared<KKSTRAWHITCLUSTER>(strawhitptr));
          hits.emplace_back(std::static_pointer_cast<MEAS>(strawhitclusters_.back()));
        }
      }
      // now add material between clusters
      for(auto& shclusterptr : strawhitclusters_) {
        if(shclusterptr->strawHits().size()>1){
          auto trange = shclusterptr->timeRange();
          for(auto const& sxing : strawxings ) {
            if(trange.inRange(sxing->time())){
              shclusterptr->addXing(sxing);
            }
          }
          // also check existing material
          for(auto const& sxing : strawxings_ ) {
            if(trange.inRange(sxing->time())){
              shclusterptr->addXing(sxing);
            }
          }
        }
      }
    }
  }

  template <class KTRAJ> void KKTrack<KTRAJ>::convertTypes(
      KKSTRAWHITCOL const& strawhits,
      KKSTRAWXINGCOL const& strawxings,
      KKCALOHITCOL const& calohits,
      MEASCOL& hits, EXINGCOL& exings) {
    hits.reserve(strawhits_.size() + calohits_.size());
    exings.reserve(strawxings_.size());
    for(auto const& strawhit : strawhits)hits.emplace_back(std::static_pointer_cast<MEAS>(strawhit));
    for(auto const& calohit : calohits)hits.emplace_back(std::static_pointer_cast<MEAS>(calohit));
    for(auto const& strawxing : strawxings)exings.emplace_back(std::static_pointer_cast<EXING>(strawxing));
  }

  template <class KTRAJ> void KKTrack<KTRAJ>::extendTrack(Config const& config,
      KKSTRAWHITCOL const& strawhits, KKSTRAWXINGCOL const& strawxings, KKCALOHITCOL const& calohits) {
    // convert the hits and Xings to generic types and extend the track
    MEASCOL hits; // polymorphic container of hits
    EXINGCOL exings; // polymorphic container of detector element crossings
    addHitClusters(strawhits,strawxings,hits);
    if(strawhits.size() > 0 && this->config().plevel_ > 0){
      unsigned nhit(0);
      std::cout << "extended " << strawhits.size() << " hits into " << strawhitclusters_.size() << " StrawHitClusters " << std::endl;
      for (auto const& shcluster : strawhitclusters_) {
        for(auto const& strawhit : strawhits ) {
          for (auto const& sh : shcluster->strawHits()) {
            if(strawhit->strawId() == sh->strawId()){
              shcluster->print(std::cout,this->config().plevel_);
            }
          }
        }
        nhit += shcluster->strawHits().size();
      }
      if(nhit != strawhits_.size()+strawhits.size()) std::cout << "cluster hit sum doesn't match " << nhit << " " << strawhits_.size() << std::endl;
    }
    convertTypes(strawhits,strawxings,calohits,hits,exings);
    this->extend(config,hits,exings);
    // store the new hits
    strawhits_.reserve(strawhits_.size()+strawhits.size());
    calohits_.reserve(calohits_.size()+calohits.size());
    strawxings_.reserve(strawxings_.size()+strawxings.size());
    for(auto const& strawhit : strawhits)strawhits_.emplace_back(strawhit);
    for(auto const& calohit : calohits)calohits_.emplace_back(calohit);
    for(auto const& strawxing : strawxings)strawxings_.emplace_back(strawxing);
  }

  template <class KTRAJ> void KKTrack<KTRAJ>::extendTrack( EXINGCOL const& exings) {
    MEASCOL nohits;
    this->extend(this->config(),nohits,exings);
  }

  template <class KTRAJ> void KKTrack<KTRAJ>::printFit(std::ostream& ost,int printlevel) const {
    if(printlevel > 1) std::cout << "Seed Helix " << this->seedTraj() << std::endl;
    TRACK::print(ost,0);
    ost << "Fit with " << strawhits_.size() << " StrawHits and " << calohits_.size() << " CaloHits and " << strawxings_.size() << " Straw Xings" << std::endl;
    if(printlevel > 2){
      for(auto const& strawhit : strawhits_) strawhit->print(std::cout,2);
      for(auto const& calohit : calohits_) calohit->print(std::cout,2);
      for(auto const& strawxing :strawxings_) strawxing->print(std::cout,2);
      for(auto const& strawhitcluster :strawhitclusters_) strawhitcluster->print(std::cout,2);
    }
  }

  template <class KTRAJ> void KKTrack<KTRAJ>::addIPAXing(KKIPAXINGPTR const& ipaxingptr,TimeDir const& tdir) {
    // convert to a generic Xing
    std::shared_ptr<KinKal::ElementXing<KTRAJ>> exptr = std::static_pointer_cast<KinKal::ElementXing<KTRAJ>>(ipaxingptr);
    // extrapolate the fit throug this xing
    if(!this->extrapolate(exptr,tdir))throw cet::exception("RECO")<<"mu2e::KKTrack: Shell extrapolation failure"<< std::endl;
    // store the xing
    ipaxings_.push_back(ipaxingptr);
  }

  template <class KTRAJ> void KKTrack<KTRAJ>::addSTXing(KKSTXINGPTR const& stxingptr,TimeDir const& tdir) {
    // convert to a generic Xing
    std::shared_ptr<KinKal::ElementXing<KTRAJ>> exptr = std::static_pointer_cast<KinKal::ElementXing<KTRAJ>>(stxingptr);
    // extrapolate the fit through this xing
    if(!this->extrapolate(exptr,tdir))throw cet::exception("RECO")<<"mu2e::KKTrack: Shell extrapolation failure"<< std::endl;
    // store the xing
    stxings_.push_back(stxingptr);
  }

  template <class KTRAJ> void KKTrack<KTRAJ>::addTCRVXing(KKCRVXINGPTR const& crvxingptr,TimeDir const& tdir) {
    // convert to a generic Xing
    std::shared_ptr<KinKal::ElementXing<KTRAJ>> exptr = std::static_pointer_cast<KinKal::ElementXing<KTRAJ>>(crvxingptr);
    // extrapolate the fit throug this xing
    if(!this->extrapolate(exptr,tdir))throw cet::exception("RECO")<<"mu2e::KKTrack: Shell extrapolation failure"<< std::endl;
    // store the xing
    crvxings_.push_back(crvxingptr);
  }



}
#endif

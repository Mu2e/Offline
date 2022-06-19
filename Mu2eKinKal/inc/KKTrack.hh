#ifndef Mu2eKinKal_KKTrack_hh
#define Mu2eKinKal_KKTrack_hh
//
// subclass of KinKal Track specialized for Mu2e
//
#include "KinKal/Fit/Track.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHitCluster.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawXing.hh"
#include "Offline/Mu2eKinKal/inc/KKCaloHit.hh"
namespace mu2e {

  using KinKal::Config;
  using KinKal::BFieldMap;

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
      using TRACK = KinKal::Track<KTRAJ>;
      // construct from configuration, fit environment, and hits and materials
      KKTrack(Config const& config, BFieldMap const& bfield, KTRAJ const& seedtraj, PDGCode::type tpart, KKSTRAWHITCLUSTERER const& shclusterer,
          KKSTRAWHITCOL const& strawhits, KKSTRAWXINGCOL const& strawxings, KKCALOHITCOL const& calohits );
      // extend the track according to new configuration, hits, and/or exings
      void extendTrack(Config const& config,
          KKSTRAWHITCOL const& strawhits, KKSTRAWXINGCOL const& strawxings, KKCALOHITCOL const& calohits );
      // accessors
      PDGCode::type fitParticle() const { return tpart_;}
      KKSTRAWHITCOL const& strawHits() const { return strawhits_; }
      KKSTRAWHITCLUSTERCOL const& strawHitClusters() const { return strawhitclusters_; }
      KKSTRAWXINGCOL const& strawXings() const { return strawxings_; }
      KKCALOHITCOL const& caloHits() const { return calohits_; }
      void printFit(std::ostream& ost=std::cout,int detail=0) const;
    private:
      // record the particle type
      PDGCode::type tpart_;
      KKSTRAWHITCLUSTERER shclusterer_; // clustering functor
      KKSTRAWHITCOL strawhits_;  // straw hits used in this fit
      KKSTRAWXINGCOL strawxings_;  // straw material crossings used in this fit
      KKSTRAWHITCLUSTERCOL strawhitclusters_;  // straw hit clusters used in this fit
      KKCALOHITCOL calohits_;  // calo hits used in this fit
      // utility function to convert to generic types
      void convertTypes( KKSTRAWHITCOL const& strawhits, KKSTRAWXINGCOL const& strawxings,KKCALOHITCOL const& calohits,
          MEASCOL& hits, EXINGCOL& exings);
      // add hits to clusters
      void addHitClusters(KKSTRAWHITCOL const& strawhits,MEASCOL& hits);
  };

  template <class KTRAJ> KKTrack<KTRAJ>::KKTrack(Config const& config, BFieldMap const& bfield, KTRAJ const& seedtraj, PDGCode::type tpart,
      KKSTRAWHITCLUSTERER const& shclusterer,
      KKSTRAWHITCOL const& strawhits,
      KKSTRAWXINGCOL const& strawxings,
      KKCALOHITCOL const& calohits) :
    KinKal::Track<KTRAJ>(config,bfield,seedtraj), tpart_(tpart), shclusterer_(shclusterer),
    strawhits_(strawhits),
    strawxings_(strawxings),
    calohits_(calohits) {
      MEASCOL hits; // polymorphic container of hits
      EXINGCOL exings; // polymorphic container of detector element crossings
      // add the hits to clusters, as required
      addHitClusters(strawhits_,hits);
      if(this->config().plevel_ > 0){
        std::cout << "created " << strawhitclusters_.size() << " StrawHitClusters " << std::endl;
        for (auto const& shcluster : strawhitclusters_) {
          shcluster->print(std::cout,this->config().plevel_);
        }
      }
      convertTypes(strawhits_, strawxings_, calohits_,  hits,exings);
      this->fit(hits,exings);
    }

  template <class KTRAJ> void KKTrack<KTRAJ>::addHitClusters(KKSTRAWHITCOL const& strawhits,MEASCOL& hits) {
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
    KKSTRAWHITCLUSTERCOL strawhitcluster;
    // convert the hits and Xings to generic types and extend the track
    MEASCOL hits; // polymorphic container of hits
    EXINGCOL exings; // polymorphic container of detector element crossings
    addHitClusters(strawhits,hits);
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

  template <class KTRAJ> void KKTrack<KTRAJ>::printFit(std::ostream& ost,int printlevel) const {
    if(printlevel > 1) std::cout << "Seed Helix " << this->seedTraj() << std::endl;
    TRACK::print(ost,0);
    ost << "Fit with " << strawhits_.size() << " StrawHits and " << calohits_.size() << " CaloHits and " << strawxings_.size() << " Straw Xings" << std::endl;
    if(printlevel > 2){
      for(auto const& strawhit : strawhits_) strawhit->print(std::cout,2);
      for(auto const& calohit : calohits_) calohit->print(std::cout,2);
      for(auto const& strawxing :strawxings_) strawxing->print(std::cout,2);
    }
  }

}
#endif

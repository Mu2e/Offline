#ifndef Mu2eKinKal_KKTrack_hh
#define Mu2eKinKal_KKTrack_hh
//
// subclass of KinKal Track specialized for Mu2e
//
#include "KinKal/Fit/Track.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHitGroup.hh"
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
      using KKSTRAWHITGROUP = KKStrawHitGroup<KTRAJ>;
      using KKSTRAWHITGROUPPTR = std::shared_ptr<KKSTRAWHITGROUP>;
      using KKSTRAWHITGROUPCOL = std::vector<KKSTRAWHITGROUPPTR>;
      using KKSTRAWHITGROUPER = KKStrawHitGrouper<KTRAJ>;
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
      KKTrack(Config const& config, BFieldMap const& bfield, KTRAJ const& seedtraj, PDGCode::type tpart, KKSTRAWHITGROUPER const& shgrouper,
          KKSTRAWHITCOL const& strawhits, KKSTRAWXINGCOL const& strawxings, KKCALOHITCOL const& calohits );
      // extend the track according to new configuration, hits, and/or exings
      void extendTrack(Config const& config,
          KKSTRAWHITCOL const& strawhits, KKSTRAWXINGCOL const& strawxings, KKCALOHITCOL const& calohits );
      // accessors
      PDGCode::type fitParticle() const { return tpart_;}
      KKSTRAWHITCOL const& strawHits() const { return strawhits_; }
      KKSTRAWHITGROUPCOL const& strawHitGroups() const { return strawhitgroups_; }
      KKSTRAWXINGCOL const& strawXings() const { return strawxings_; }
      KKCALOHITCOL const& caloHits() const { return calohits_; }
      void printFit(std::ostream& ost=std::cout,int detail=0) const;
    private:
      // record the particle type
      PDGCode::type tpart_;
      KKSTRAWHITGROUPER shgrouper_; // grouping functor
      KKSTRAWHITCOL strawhits_;  // straw hits used in this fit
      KKSTRAWXINGCOL strawxings_;  // straw material crossings used in this fit
      KKSTRAWHITGROUPCOL strawhitgroups_;  // straw hit groups used in this fit
      KKCALOHITCOL calohits_;  // calo hits used in this fit
      // utility function to convert to generic types
      void convertTypes( KKSTRAWHITCOL const& strawhits, KKSTRAWXINGCOL const& strawxings,KKCALOHITCOL const& calohits,
          MEASCOL& hits, EXINGCOL& exings);
      // add hits to groups
      void addHitGroups(KKSTRAWHITCOL const& strawhits,MEASCOL& hits);
  };

  template <class KTRAJ> KKTrack<KTRAJ>::KKTrack(Config const& config, BFieldMap const& bfield, KTRAJ const& seedtraj, PDGCode::type tpart,
      KKSTRAWHITGROUPER const& shgrouper,
      KKSTRAWHITCOL const& strawhits,
      KKSTRAWXINGCOL const& strawxings,
      KKCALOHITCOL const& calohits) :
    KinKal::Track<KTRAJ>(config,bfield,seedtraj), tpart_(tpart), shgrouper_(shgrouper),
    strawhits_(strawhits),
    strawxings_(strawxings),
    calohits_(calohits) {
      MEASCOL hits; // polymorphic container of hits
      EXINGCOL exings; // polymorphic container of detector element crossings
      // add the hits to groups, as required
      addHitGroups(strawhits_,hits);
      if(this->config().plevel_ > 0){
        std::cout << "created " << strawhitgroups_.size() << " StrawHitGroups " << std::endl;
        for (auto const& shgroup : strawhitgroups_) {
          shgroup->print(std::cout,this->config().plevel_);
        }
      }
      convertTypes(strawhits_, strawxings_, calohits_,  hits,exings);
      this->fit(hits,exings);
    }

  template <class KTRAJ> void KKTrack<KTRAJ>::addHitGroups(KKSTRAWHITCOL const& strawhits,MEASCOL& hits) {
    if(shgrouper_.groupLevel() != StrawIdMask::none){
      for(auto const& strawhitptr : strawhits){
        bool added(false);
        for(auto& shgroupptr : strawhitgroups_) {
          if(shgroupptr->canAddHit(strawhitptr,shgrouper_)){
            shgroupptr->addHit(strawhitptr,shgrouper_);
            added = true;
            break;
          }
        }
        if(!added){
          strawhitgroups_.emplace_back(std::make_shared<KKSTRAWHITGROUP>(strawhitptr));
          hits.emplace_back(static_pointer_cast<MEAS>(strawhitgroups_.back()));
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
    for(auto const& strawhit : strawhits)hits.emplace_back(static_pointer_cast<MEAS>(strawhit));
    for(auto const& calohit : calohits)hits.emplace_back(static_pointer_cast<MEAS>(calohit));
    for(auto const& strawxing : strawxings)exings.emplace_back(static_pointer_cast<EXING>(strawxing));
  }

  template <class KTRAJ> void KKTrack<KTRAJ>::extendTrack(Config const& config,
      KKSTRAWHITCOL const& strawhits, KKSTRAWXINGCOL const& strawxings, KKCALOHITCOL const& calohits) {
    KKSTRAWHITGROUPCOL strawhitgroup;
    // convert the hits and Xings to generic types and extend the track
    MEASCOL hits; // polymorphic container of hits
    EXINGCOL exings; // polymorphic container of detector element crossings
    addHitGroups(strawhits,hits);
    if(strawhits.size() > 0 && this->config().plevel_ > 0){
      unsigned nhit(0);
      std::cout << "extended " << strawhits.size() << " hits into " << strawhitgroups_.size() << " StrawHitGroups " << std::endl;
      for (auto const& shgroup : strawhitgroups_) {
        for(auto const& strawhit : strawhits ) {
          for (auto const& sh : shgroup->strawHits()) {
            if(strawhit->strawId() == sh->strawId()){
              shgroup->print(std::cout,this->config().plevel_);
            }
          }
        }
        nhit += shgroup->strawHits().size();
      }
      if(nhit != strawhits_.size()+strawhits.size()) std::cout << "group hit sum doesn't match " << nhit << " " << strawhits_.size() << std::endl;
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

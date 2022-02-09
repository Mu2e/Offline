#ifndef Mu2eKinKal_KKTrack_hh
#define Mu2eKinKal_KKTrack_hh
//
// subclass of KinKal Track specialized for Mu2e
//
#include "KinKal/Fit/Track.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
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
      // construct from configuration, fit environment, and hits and materials
      KKTrack(Config const& config, BFieldMap const& bfield, KTRAJ const& seedtraj, PDGCode::type tpart,
          KKSTRAWHITCOL const& strawhits, KKCALOHITCOL const& calohits, KKSTRAWXINGCOL const& strawxings);
      // extend the track according to new configuration, hits, and/or exings
      void extendTrack(Config const& config,
          KKSTRAWHITCOL const& strawhits, KKCALOHITCOL const& calohits, KKSTRAWXINGCOL const& strawxings);
      // accessors
      PDGCode::type fitParticle() const { return tpart_;}
      KKSTRAWHITCOL const& strawHits() const { return strawhits_; }
      KKSTRAWXINGCOL const& strawXings() const { return strawxings_; }
      KKCALOHITCOL const& caloHits() const { return calohits_; }
    private:
      // record the particle type
      PDGCode::type tpart_;
      KKSTRAWHITCOL strawhits_;  // straw hits used in this fit
      KKCALOHITCOL calohits_;  // calo hits used in this fit
      KKSTRAWXINGCOL strawxings_;  // straw material crossings used in this fit
      // utility function to convert to generic types
      void convertTypes( KKSTRAWHITCOL const& strawhits, KKCALOHITCOL const& calohits, KKSTRAWXINGCOL const& strawxings,
          MEASCOL& hits, EXINGCOL& exings);
  };
  template <class KTRAJ> KKTrack<KTRAJ>::KKTrack(Config const& config, BFieldMap const& bfield, KTRAJ const& seedtraj, PDGCode::type tpart,
      KKSTRAWHITCOL const& strawhits, KKCALOHITCOL const& calohits, KKSTRAWXINGCOL const& strawxings) :
    KinKal::Track<KTRAJ>(config,bfield,seedtraj), tpart_(tpart), strawhits_(strawhits), calohits_(calohits), strawxings_(strawxings)
  {
    MEASCOL hits; // polymorphic container of hits
    EXINGCOL exings; // polymorphic container of detector element crossings
    convertTypes(strawhits_, calohits_, strawxings_,hits,exings);
    this->fit(hits,exings);
  }

  template <class KTRAJ> void KKTrack<KTRAJ>::convertTypes( KKSTRAWHITCOL const& strawhits, KKCALOHITCOL const& calohits, KKSTRAWXINGCOL const& strawxings,
      MEASCOL& hits, EXINGCOL& exings) {
    hits.reserve(strawhits_.size() + calohits_.size());
    exings.reserve(strawxings_.size());
    for(auto const& strawhit : strawhits)hits.emplace_back(static_pointer_cast<MEAS>(strawhit));
    for(auto const& calohit : calohits)hits.emplace_back(static_pointer_cast<MEAS>(calohit));
    for(auto const& strawxing : strawxings)exings.emplace_back(static_pointer_cast<EXING>(strawxing));
  }

  template <class KTRAJ> void KKTrack<KTRAJ>::extendTrack(Config const& config,
      KKSTRAWHITCOL const& strawhits, KKCALOHITCOL const& calohits, KKSTRAWXINGCOL const& strawxings) {
    // convert the hits and Xings to generic types and extend the track
    MEASCOL hits; // polymorphic container of hits
    EXINGCOL exings; // polymorphic container of detector element crossings
    convertTypes(strawhits, calohits, strawxings,hits,exings);
    this->extend(config,hits,exings);
    // store the new hits
    strawhits_.reserve(strawhits_.size()+strawhits.size());
    calohits_.reserve(calohits_.size()+calohits.size());
    strawxings_.reserve(strawxings_.size()+strawxings.size());
    for(auto const& strawhit : strawhits)strawhits_.emplace_back(strawhit);
    for(auto const& calohit : calohits)calohits_.emplace_back(calohit);
    for(auto const& strawxing : strawxings)strawxings_.emplace_back(strawxing);
  }

}
#endif

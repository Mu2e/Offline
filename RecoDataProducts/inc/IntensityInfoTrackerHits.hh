//
// Class to collect the info needed for monitoring POT / stop muons
// bvitali May 2021
//


#ifndef RecoDataProducts_IntensityInfoTrackerHits_hh
#define RecoDataProducts_IntensityInfoTrackerHits_hh

#include <vector>

namespace mu2e {

  class IntensityInfoTrackerHits
  {
  public:
    IntensityInfoTrackerHits() {}
    IntensityInfoTrackerHits(unsigned short nTrackerHits):
      nTrackerHits_(nTrackerHits)
    {}


    void setNTrackerHits   (unsigned short tmp) {nTrackerHits_ = tmp;}

    unsigned short nTrackerHits () const { return nTrackerHits_; }

  private:
    unsigned short  nTrackerHits_ = 0;
  };
  typedef std::vector<mu2e::IntensityInfoTrackerHits> IntensityInfosTrackerHits;
}

#endif

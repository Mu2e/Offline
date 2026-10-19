//
// Class to collect the info needed for monitoring POT / stop muons
// bvitali May 2021
//


#ifndef RecoDataProducts_IntensityInfoTimeCluster_hh
#define RecoDataProducts_IntensityInfoTimeCluster_hh

#include <vector>

namespace mu2e {

  class IntensityInfoTimeCluster
  {
  public:
    IntensityInfoTimeCluster() {}
    IntensityInfoTimeCluster(unsigned short nTCs):
      nTCs_(nTCs)
    {}


    void setNTCs(unsigned short tmp) {nTCs_ = tmp;}

    unsigned short nTCs() const { return nTCs_; }

  private:
    unsigned short  nTCs_ = 0;
  };

  typedef std::vector<mu2e::IntensityInfoTimeCluster> IntensityInfosTimeCluster;
}

#endif

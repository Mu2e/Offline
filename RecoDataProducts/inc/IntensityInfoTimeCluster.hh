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
    IntensityInfoTimeCluster(unsigned short nProtonTCs):
      nProtonTCs_(nProtonTCs)
    {}


    void setNProtonTCs     (unsigned short tmp) {nProtonTCs_   = tmp;}

    unsigned short nProtonTCs   () const { return nProtonTCs_  ; }

  private:
    unsigned short  nProtonTCs_   = 0;
  };

  typedef std::vector<mu2e::IntensityInfoTimeCluster> IntensityInfosTimeCluster;
}

#endif

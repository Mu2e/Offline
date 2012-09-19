#include "RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"


#include <iostream>
//#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
#define AGDEBUG(stuff)

namespace mu2e {

  ExtMonFNALRecoClusterCollection::ExtMonFNALRecoClusterCollection(unsigned nplanes)
    : storage_()
    , planes_(nplanes)
  {}

  void ExtMonFNALRecoClusterCollection::insert(const ExtMonFNALRecoCluster& c) {
    storage_.push_back(c);
    planes_[c.plane()].push_back(storage_.size() - 1);
    AGDEBUG("storage_.size() = "<<storage_.size()<<", planes_.size()="<<planes_.size() 
            <<", planes_[c.plane()].indexes_.back()="<< planes_[c.plane()].indexes_.back()
            <<", planes_[c.plane()].parent_="<<planes_[c.plane()].parent_
            <<", this="<<this
            );
  }

} // namespace mu2e

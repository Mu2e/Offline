// Reco clusters grouped by plane.
//
// $Id: ExtMonFNALRecoClusterCollection.hh,v 1.3 2012/11/01 23:44:32 gandr Exp $
// $Author: gandr $
// $Date: 2012/11/01 23:44:32 $
//
// Original author Andrei Gaponenko
//

#ifndef RecoDataProducts_ExtMonFNALRecoClusterCollection_hh
#define RecoDataProducts_ExtMonFNALRecoClusterCollection_hh

#include <vector>

#include "canvas/Persistency/Common/Wrapper.h"

#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"

namespace mu2e {

  class ExtMonFNALRecoClusterCollection {
  public:

    //----------------------------------------------------------------
    // A transient class to provide per-plane view of the clusters
    class PlaneClusters {
      typedef ExtMonFNALRecoClusterCollection Parent;
      friend class ExtMonFNALRecoClusterCollection;
      const Parent *parent_;
      unsigned plane_;
    public:
      unsigned int size() const { return parent_->planes_[plane_].size(); }
      bool empty() const { return !size(); }
      const ExtMonFNALRecoCluster& operator[](unsigned int i) const {
        return parent_->storage_[ parent_->planes_[plane_][i] ];
      }

      explicit PlaneClusters(const Parent* p, unsigned plane) : parent_(p), plane_(plane) {}
    };

    //----------------------------------------------------------------
    // Main collection interface
    explicit ExtMonFNALRecoClusterCollection(unsigned nplanes);

    unsigned nplanes() const { return planes_.size(); }

    PlaneClusters clusters(unsigned plane) const { return PlaneClusters(this, plane); }

    void insert(const ExtMonFNALRecoCluster& c);


  private:

    typedef std::vector<ExtMonFNALRecoCluster> StorageImpl;
    StorageImpl storage_;
    std::vector<std::vector<unsigned int> > planes_;

    // For persistency
    template<class T> friend class art::Wrapper;
    ExtMonFNALRecoClusterCollection() {}

  public:

    //----------------------------------------------------------------
    // Extra public stuff needed by art::Ptr
    typedef StorageImpl::const_iterator const_iterator;
    typedef StorageImpl::value_type value_type;
    const_iterator begin() const { return storage_.begin(); }
    const_iterator end() const { return storage_.end(); }
    const_iterator cbegin() const { return storage_.cbegin(); }
    const_iterator cend() const { return storage_.cend(); }
    const ExtMonFNALRecoCluster& operator[](std::size_t globalIndex) const { return storage_[globalIndex]; }

    // clients need a correct element index to create art::Ptr
    std::size_t globalIndex(unsigned plane, unsigned clusterInPlane) const { return planes_[plane][clusterInPlane]; }

    std::size_t size() const { return storage_.size(); }
  };

} // namespace mu2e

//----------------
// This specialization is needed to support read-back of art::Ptr<ExtMonFNALRecoCluster>
namespace art {
  template<> struct has_setPtr<mu2e::ExtMonFNALRecoClusterCollection> { static const bool value = true; };
}

//----------------

#endif /* RecoDataProducts_ExtMonFNALRecoClusterCollection_hh */

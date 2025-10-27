// Ed Callaghan
// Container for CaloDigiWrappers, with functionality to resolve overlapping digitization windows
// September 2025

#ifndef CaloMC_CaloDigiWrapperCollection_hh
#define CaloMC_CaloDigiWrapperCollection_hh

// stl
#include <algorithm>
#include <map>
#include <memory>
#include <vector>

// mu2e
#include "Offline/CaloMC/inc/CaloDigiWrapper.hh"
#include "Offline/TrackerMC/inc/TimeBasedBucket.hh"

namespace mu2e{
  using CDWC_iterator = std::vector<CaloDigiWrapper>::iterator;
  using CDWC_const_iterator = std::vector<CaloDigiWrapper>::const_iterator;

  class CaloDigiWrapperCollection{
    public:
      using SiPMID_t = int;
      using sample_t = int;
      using pos_t = size_t;
      // forwarded calls to underlying container
      size_t size() const;
      CDWC_iterator begin();
      CDWC_const_iterator begin() const;
      CDWC_iterator end();
      CDWC_const_iterator end() const;
      CaloDigiWrapper& operator[](size_t);
      const CaloDigiWrapper& operator[](size_t) const;
      CaloDigiWrapper&       front();
      const CaloDigiWrapper& front() const;
      CaloDigiWrapper&       back();
      const CaloDigiWrapper& back() const;

      // insertion
      void Append(const CaloDigiCollection&);
      void Append(const CaloDigiWrapper&);

      // accessor
      std::unique_ptr<CaloDigiCollection> GetDigis() const;

      // identify sets of  digis with overlapping digitization windows, and
      // reduce each such set to a single digi, representing their "sum"
      void ResolveCollisions(CaloDigiWrapperCollection&);

    protected:
      std::vector<CaloDigiWrapper> _wrappers;
      void ResolveCollision(CaloDigiWrapperCollection&,
                            CaloDigiWrapperCollection&);

    private:
      /**/
  };
} // namespace mu2e

#endif

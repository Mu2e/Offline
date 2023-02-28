#ifndef RecoDataProducts_ComboHit_hh
#define RecoDataProducts_ComboHit_hh
//
// Class to describe a combination of one or more StrawHits.
// These can be aggregated together through the dedicated collection
//
// Original author David Brown Dec 2017
//
// Mu2e includes
#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/StrawIdMask.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include <stdint.h>
// art includes
#ifndef __ROOTCLING__
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#endif
// C++ includes
#include <array>
#include <vector>
namespace mu2e {

  struct ComboHit {
    enum edir{wire=0,trans};
    constexpr static size_t MaxNCombo = 8; // needs tuning FIXME!
    using PIArray = std::array<uint16_t,MaxNCombo>; // array of indices into parent collection
    // General accessors that apply to all kinds of combo hits
    XYZVectorF const& pos() const { return _pos; }
    XYZVectorF const& wdir() const { return _wdir; } // along wire
    XYZVectorF const& sdir() const { return _sdir; } // perp to wire and Z
    XYZVectorF centerPos() const { return _pos - _wdist*_wdir; }
//
    float posRes(edir dir) const;
    float dEdx() const { return _dedx; }
    float energyDep() const { return _dedx*_pathlength; }
    float pathLength() const { return _pathlength; }
    float phi() const { return _pos.phi();} // legacy function
    float timeRes() const { return _timeres; }
    float timeVar() const { return _timeres*_timeres; }
    float correctedTime() const { return _time; }
    float qual() const { return _qual; }
    StrawHitFlag const& flag() const { return _flag; }
    StrawId const& strawId() const { return _sid; }
    float wireRes() const { return _wres; }
    float transRes() const { return _tres; }
    float transVar() const { return _tres*_tres; }
    float wireVar() const { return _wres*_wres; }
    float wireDist() const { return _wdist; }
    // book-keeping accessors
    uint16_t nCombo() const { return _ncombo; }
    uint16_t nStrawHits() const { return _nsh; }
    StrawIdMask const& mask() const { return _mask;}
    void init(ComboHit const& other, uint16_t index);
    uint16_t index(uint16_t ish=0) const { return _pind.at(ish); }
    bool addIndex(uint16_t shi); // append an index to the
    PIArray const& indexArray() const { return _pind; }
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;
    // Accessors that only make sense for single-straw Combo hits
    float TOT(StrawEnd end=StrawEnd::cal)       const { return _tot[end];}
     // compatibility constructor (deprecated)
    float endTime(StrawEnd end=StrawEnd::cal)     const { return _ttdc[end];}
    auto const& TOTs() const { return _tot; }
    auto const& endTimes() const { return _ttdc; }
    StrawEnd const& earlyEnd() const { return _eend; } // End with earliest tdc time
    StrawEnd lateEnd() const { return _eend.otherEnd(); } // End with later tdc time
    // Accessors for hits used in helices
    float helixPhi() const { return _hphi;}
    // interface returning calibration info: these should be refactored to
    // use StrawResponse TODO
    float driftTime() const { return _dtime; }
    float propTime() const { return _ptime; }
    // legacy functions
    //  No new code should use these accessors, they should be removed soon TODO
//    ComboHit(const ComboHit&, StrawHitIndex, double);
    float time() const { return _ttdc[_eend]; }
    CLHEP::Hep3Vector posCLHEP() const { return GenVector::Hep3Vec(pos()); }
    // persistent payload
    // vector information.  These are stored explicitly even though they are reducible to a
    // smaller payload as the processing time cost of reconstituting them is higher than the memory access time cost
    XYZVectorF _pos; // best estimate of the position of this hit in space
    XYZVectorF _wdir; // 'direction' of this hit, used to define error elipsoid axis
    XYZVectorF _sdir; // straw radial direction, perp to Z and wire direction
    float _wres = -1.0;
    float _tres = -1.0; // resolution along and transverse to the 'wire' direction
    float _wdist = 0.0; // distance from wire center along the wire direction: this can be derived from _pos and _wdir
    float _time = 0.0; // best estimate of time the physical particle created this hit: aggregate and calibrated
    float _timeres = -1.0; // estimated resolution of time measurement
    float _dedx = 0.0; // average energy loss per unit length
    float _pathlength = 0.0; // straw gas volume path length estimate
    float _qual = 0.0;; // quality of hit or combination
    StrawHitFlag _flag; // condition of this hit
    StrawId _sid; // straw identifier; for composites, not all fields are complete, use in conjunction with mask
    StrawIdMask _mask; // mask of valid StrawId fields
    StrawEnd _eend; // early tdc end
    // low-level quantities needed for calibration.  These only make sense for single-straw hits
    TrkTypes::TOTTimes  _tot = {0.0, 0.0 };   // TOT times in ns from each end
    TrkTypes::TDCTimes _ttdc = {0.0, 0.0 }; // threshold crossing times in ns from each end
    // bookkeeping info
    uint16_t _ncombo = 0; // number of associated input objects
    uint16_t _nsh = 0; // number of underlying straw hits
    PIArray _pind = {0,0,0,0,0,0,0,0}; // Indices back to parent objects
    // information specific to hits associated with a helix
    float _hphi = 0.0; // azimuth relative to a helix center
    float _xyWeight = 0.0;       // weight used to perform the x-y circle fit
    float _zphiWeight = 0.0;     // weight used to perfom the z-phi linear fit
    // low-level derived data that should move to StrawResponse
    float _dtime = 0.0; // TOT based drift time estimate
    float _ptime = 0.0; // prop time estimate
  };
  // ComboHitCollection is a non-trivial subclass of vector which includes navigation of nested ComboHits
  class ComboHitCollection : public std::vector<mu2e::ComboHit> {
    public:
      ComboHitCollection(bool sorted=false) : _sorted(sorted) {}
      typedef std::vector<ComboHitCollection::const_iterator> CHCIter;
      // fill a vector of indices to the underlying digis used in a given ComboHit
      // This function is called recursively, so the the vector must be empty on the top-most call
#ifndef __ROOTCLING__
      void fillStrawDigiIndices(art::Event const& event, uint16_t chindex, std::vector<StrawHitIndex>& shids) const;
      // similarly fill to the StrawHit level
      void fillStrawHitIndices(art::Event const& event, uint16_t chindex, std::vector<StrawHitIndex>& shids) const;
      // do this for all the hits in the collection
      void fillStrawHitIndices(art::Event const& event, std::vector<std::vector<StrawHitIndex> >& shids) const;
      // translate a collection of ComboHits into the lowest-level (straw) combo hits.  This function is recursive
      void fillComboHits(art::Event const& event, std::vector<uint16_t> const& indices, CHCIter& iters) const;
      // fill a vector of iterators to the ComboHits 1 layer below a given ComboHit.  This is NOT RECURSIVE
      // return value says whether there's a layer below or not (if not, output is empty)
      bool fillComboHits(art::Event const& event, uint16_t chindex, CHCIter& iters) const;
      // recover the parent collection handle from the event
      void setParentHandle(art::Event const& event, art::Handle<ComboHitCollection>& phandle) const;
      // set the parent Id given a handle to the parent collection
      void setParent(art::Handle<ComboHitCollection> const& phandle);
      // or directly from the product ID
#endif
      void setParent(art::ProductID const& par){ _parent = par; }
      // accessors
      art::ProductID const& parent() const { return _parent; }
      bool sorted() const { return _sorted; }
      uint16_t nStrawHits() const;
    private:
      // reference back to the input ComboHit collection this one references
      // This can be used to chain back to the original StrawHit indices
      art::ProductID _parent;
      bool _sorted; // record if this collection was sorted
  };
  inline std::ostream& operator<<( std::ostream& ost,
                                   ComboHit const& hit){
    hit.print(ost,false);
    return ost;
  }
}
#endif

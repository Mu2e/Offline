#ifndef RecoDataProducts_ComboHit_hh
#define RecoDataProducts_ComboHit_hh
//
// Class to describe a combination of one or more StrawHits.
// These can be aggregated together through the dedicated collection
//
// Original author David Brown Dec 2017
//
// Mu2e includes
#include "DataProducts/inc/StrawEnd.hh"
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawIdMask.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include <stdint.h>
// root includes
#include "Rtypes.h"
// art includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
// C++ includes
#include <array>
#include <vector>
namespace mu2e {

  struct ComboHit {
    enum edir{wire=0,trans};
    constexpr static size_t MaxNCombo = 8; // needs tuning FIXME!
    typedef std::array<uint16_t,MaxNCombo> PIArray; // array of indices into parent collection
    ComboHit();
    // compatibility constructor (deprecated)
    ComboHit(const ComboHit&, StrawHitIndex, double);
    // accessors
    XYZVec centerPos() const { return _pos - _wdist*_wdir; }
    XYZVec const& pos() const { return _pos; }
    XYZVec const& wdir() const { return _wdir; }
// CLHEP-versions of these for backwards compatibilty
    CLHEP::Hep3Vector centerPosCLHEP() const { return Geom::Hep3Vec(_pos - _wdist*_wdir); }
    CLHEP::Hep3Vector posCLHEP() const { return Geom::Hep3Vec(_pos); }
    CLHEP::Hep3Vector wdirCLHEP() const { return Geom::Hep3Vec(_wdir); }
//
    Float_t posRes(edir dir) const;
    Float_t energyDep() const { return _edep; }
    Float_t phi() const { return _pos.phi();}
    Float_t helixPhi() const { return _hphi;}
    Float_t time() const { return _time; }
    Float_t driftTime() const { return _dtime; }
    Float_t propTime() const { return _ptime; }
    Float_t correctedTime() const { return _time - _ptime - _dtime; }
    Float_t specificIonization() const { return _edep/_pathlength; }
    Float_t pathLength() const { return _pathlength; }
    Float_t qual() const { return _qual; }
    StrawHitFlag const& flag() const { return _flag; }
    StrawEnd const& driftEnd() const { return _tend; } // which end was used for time
    StrawId const& strawId() const { return _sid; }
    Float_t wireRes() const { return _wres; }
    Float_t transRes() const { return _tres; }
    Float_t transErr2() const { return _tres*_tres; }
    Float_t wireErr2() const { return _wres*_wres; }
    Float_t wireDist() const { return _wdist; }
    uint16_t nCombo() const { return _ncombo; }
    uint16_t nStrawHits() const { return _nsh; }
    StrawIdMask const& mask() const { return _mask;}
    void init(ComboHit const& other, uint16_t index);
    uint16_t index(uint16_t ish=0) const;
    bool addIndex(uint16_t shi); // append an index to the
    PIArray const& indexArray() const { return _pind; }
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;
    //
    XYZVec _pos; // position of this hit
    XYZVec _wdir; // 'direction' of this hit, used to define error elipsoid axis
    XYZVec _sdir;           // straw radial direction, perp to Z and wire direction
    Float_t _wres, _tres; // resolution along and transverse to the 'wire' direction
    Float_t _wdist; // distance from wire center along this direction (agregate)
    Float_t _time, _edep, _qual; // derived StrawHit (agregate) info
    Float_t _dtime; // drift time estimate
    Float_t _ptime; // prop time estimate
    Float_t _pathlength; // path length estimate
    Float_t _hphi; // azimuth relative to a helix center
    Float_t _xyWeight;       // weight used to perform the x-y circle fit
    Float_t _zphiWeight;     // weight used to perfom the z-phi linear fit
    uint16_t _ncombo; // number of associated input objects
    uint16_t _nsh; // number of underlying straw hits
    PIArray _pind; // Indices back to parent objects
    StrawHitFlag _flag; // flag condition of this hit (agregate)
    StrawId _sid; // straw identifier; some fields may not be complete, use in conjunction with mask
    StrawIdMask _mask; // mask for valid StrawId fields
    StrawEnd _tend; // end used to define time measruement
  };
  // ComboHitCollection is a non-trivial subclass of vector which includes navigation of nested ComboHits
  class ComboHitCollection : public std::vector<mu2e::ComboHit> {
    public:
      ComboHitCollection(bool sorted=false) : _sorted(sorted) {}
      typedef std::vector<ComboHitCollection::const_iterator> CHCIter;
      // fill a vector of indices to the underlying digis used in a given ComboHit
      // This function is called recursively, so the the vector must be empty on the top-most call
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



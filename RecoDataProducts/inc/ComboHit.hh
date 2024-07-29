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
#include "Math/SMatrix.h"
// art includes
#ifndef __ROOTCLING__
#include "art/Framework/Principal/Handle.h"
#endif
#include "canvas/Persistency/Common/ProductPtr.h"
// C++ includes
#include <array>
#include <vector>
namespace mu2e {

  struct ComboHit {
    enum edir{wire=0,trans,z}; // should switch to  UVW TODO
    enum UVDir{UDir=0,VDir,WDir}; //  UVW
    constexpr static size_t MaxNCombo = 8;
    using PIArray = std::array<uint16_t,MaxNCombo>; // array of indices into parent collection
    // General accessors that apply to all kinds of combo hits
    // UVW coordinate system
    XYZVectorF uDir() const { return XYZVectorF(_udir.X(),_udir.Y(),0.0); }
    XYZVectorF vDir() const { return XYZVectorF(-_udir.Y(),_udir.X(),0.0); }
    XYZVectorF wDir() const { return XYZVectorF(0.0,0.0,1.0); } // by definition along Z
    // native 2D versionsof the above
    XYVectorF const& uDir2D() const { return _udir; } // along wire
    XYVectorF vDir2D() const { return XYVectorF(-_udir.Y(),_udir.X()); } // perp to wire
    // hit direction, =wDir except for stereo hits, where it can be the result of a fit
    XYZVectorF const& hDir() const { return _hdir; }
    // position
    XYZVectorF const& pos() const { return _pos; }
    float uPos() const { return _pos.Dot(uDir()); }
    float vPos() const { return _pos.Dot(vDir()); }
    XYZVectorF centerPos() const { return _pos - _wdist*uDir(); } // wire center position
    // resolution accessors
    float posRes(edir dir) const;
    // position and slope variances
    float uVar() const { return _uvar; }
    float vVar() const { return _vvar; }
    float wVar() const { return _wvar; }
    float hcostVar() const { return _hcostvar; }
    float hphiVar() const { return _hphivar; }
    // resolutions from the variances
    float uRes() const { return sqrt(uVar()); }
    float vRes() const { return sqrt(vVar()); }
    float wRes() const { return sqrt(wVar()); }
    float hcostRes() const { return sqrt(hcostVar()); }
    float hphiRes() const { return sqrt(hphiVar()); }
    // other info
    float energyDep() const { return _edep; }
    float qual() const { return _qual; }
    StrawHitFlag const& flag() const { return _flag; }
    StrawId const& strawId() const { return _sid; }
    // book-keeping accessors
    auto nCombo() const { return _ncombo; }
    auto nStrawHits() const { return _nsh; }
    auto const& mask() const { return _mask;}
    void init(ComboHit const& other, size_t index);
    auto index(size_t ish=0) const { return _pind.at(ish); }
    bool addIndex(size_t shi); // append an index to the
    auto const& indexArray() const { return _pind; }
    // general timing info, valid for all hits
    float timeVar() const { return _timevar; }
    float timeRes() const { return sqrt(_timevar); }
    float correctedTime() const { return _time; }
    // info for single-panel Combo hits
    float TOT(StrawEnd end=StrawEnd::cal)       const { return _tot[end];}
    float wireDist() const { return _wdist; }
    // compatibility constructor (deprecated)
    float endTime(StrawEnd end=StrawEnd::cal)     const { return _etime[end];}
    auto const& TOTs() const { return _tot; }
    auto const& endTimes() const { return _etime; }
    StrawEnd const& earlyEnd() const { return _eend; } // End with earliest tdc time
    StrawEnd lateEnd() const { return _eend.otherEnd(); } // End with later tdc time
    // Accessors for hits used in helices
    float helixPhi() const { return _hphi;}
    // other
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;
    // interface returning calibration info: these should be refactored to
    // use StrawResponse TODO
    float driftTime() const { return _dtime; }
    float propTime() const { return _ptime; }
    // legacy functions
    //  No new code should use these accessors, they should be removed soon TODO
    //    ComboHit(const ComboHit&, StrawHitIndex, double);
    float wireRes() const { return uRes(); }
    float wireVar() const { return uVar(); }
    float transRes() const { return vRes(); }
    float transVar() const { return vVar(); }
    float phi() const { return _pos.phi();}
    float time() const { return _etime[_eend]; }
    CLHEP::Hep3Vector posCLHEP() const { return GenVector::Hep3Vec(pos()); }
    // persistent payload
    // vector information.  These are stored explicitly even though they are reducible to a
    // smaller payload as the processing time cost of reconstituting them is higher than the memory access time cost
    XYZVectorF _pos; // best estimate of the position of this hit in space
    XYVectorF _udir; // always perp to Z, defined as the semi-major direction of the covariance matrix
    XYZVectorF _hdir; // direction of the hit
    float _uvar= 0.0, _vvar = 0.0, _wvar = 0.0; // diagonals of position covariance;
    float _hcostvar = 0.0, _hphivar = 0.0; // diagonals of hit direction covariance
    float _wdist = 0.0; // distance from wire center along the wire direction
    float _time = 0.0; // best estimate of time the physical particle created this hit: aggregate and calibrated
    float _timevar = 0.0; // estimated variance of time measurement
    float _edep = 0.0; // average energy deposition
    float _qual = 0.0; // quality of hit or combination
    StrawHitFlag _flag; // condition of this hit
    StrawId _sid; // straw identifier; for composites, not all fields are complete, use in conjunction with mask
    StrawIdMask _mask; // mask of valid StrawId fields
    StrawEnd _eend; // early tdc end
    // low-level quantities needed for calibration.  These only make sense for single-straw hits
    TrkTypes::TOTTimes  _tot = {0.0, 0.0 };   // TOT times in ns from each end
    TrkTypes::TDCTimes _etime = {0.0, 0.0 }; // threshold crossing times in ns from each end
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
      enum Sort {unsorted=0,zsort,sidsort,timesort}; // define sort state of the contents
      using CHCIter = std::vector<ComboHitCollection::const_iterator>;
      using CHCPTR = art::ProductPtr<ComboHitCollection>;
      using SHIV = std::vector<StrawHitIndex>;
      ComboHitCollection(Sort sort=unsorted) : _sort(sort) {}
      // fill a vector of indices to the underlying digis used in a given ComboHit
      // This function is called recursively, so the the vector must be empty on the top-most call
#ifndef __ROOTCLING__
      // find the parent at a given level
      // if parent and grandparent are the same level, will select the grandparent unless stopatfirst set
      CHCPTR parent(StrawIdMask::Level level, bool stopatfirst=false) const;
      void fillStrawDigiIndices( size_t chindex, SHIV& shids, bool stopatfirst=false) const;
      // Fill indices to the specified level.  Return value is the collection to whic
      // the indices apply.  first, given all my hits
      ComboHitCollection const* fillStrawHitIndices( SHIV& shiv, StrawIdMask::Level clevel=StrawIdMask::uniquestraw, bool stopatfirst=false) const;
      // given a specific hit (index) in myself
      ComboHitCollection const* fillStrawHitIndices( size_t chindex, SHIV& shiv, StrawIdMask::Level clevel=StrawIdMask::uniquestraw, bool stopatfirst=false) const;
      // given a vector of indices
      ComboHitCollection const* fillStrawHitIndices(SHIV const& inshiv, SHIV& outshiv, StrawIdMask::Level clevel=StrawIdMask::uniquestraw, bool stopatfirst=false) const;
      // the following are deprecated in favor of the more-efficient and self-checking functions above
      // translate a collection of ComboHits into the lowest-level (straw) combo hits.  This function is recursive
      void fillComboHits( std::vector<uint16_t> const& indices, CHCIter& iters) const;
      // fill a vector of iterators to the ComboHits 1 layer below a given ComboHit.  This is NOT RECURSIVE
      void fillComboHits( size_t chindex, CHCIter& iters) const;
      // set the parent Id given a handle to the parent collection
      void setParent(art::Handle<ComboHitCollection> const& phandle);
      void setParent(art::ValidHandle<ComboHitCollection> const& phandle);
      // or directly from the productPtr
      void setParent(CHCPTR const& parent);
      // or set to be the same as another collection
      void setSameParent(ComboHitCollection const& other);
      // if argument has its own parent, set this parent to match.
      // otherwise set this parent to be collection in argument
      void setAsSubset(CHCPTR const& other);
      void setAsSubset(art::Handle<ComboHitCollection> const& ohandle);
      void setAsSubset(art::ValidHandle<ComboHitCollection> const& ohandle);
      // optionally specify what level to make as parent
      // if parent and grandparent are the same level, will select the grandparent unless stopatfirst set
      void setAsSubset(CHCPTR const& optr, StrawIdMask::Level level, bool stopatfirst=false);
      void setAsSubset(art::Handle<ComboHitCollection> const& ohandle, StrawIdMask::Level level, bool stopatfirst=false);
      void setAsSubset(art::ValidHandle<ComboHitCollection> const& ohandle, StrawIdMask::Level level, bool stopatfirst=false);
#endif
      // accessors
      auto const& parent() const { return _parent; }
      StrawIdMask::Level level() const;
      auto sort() const { return _sort; }
      unsigned nStrawHits() const;
    private:
      // reference back to the input ComboHit collection this one references
      CHCPTR _parent; // pointer to the parent object
      Sort _sort; // record how this collection was sorted
  };
  inline std::ostream& operator<<( std::ostream& ost,
      ComboHit const& hit){
    hit.print(ost,false);
    return ost;
  }
}
#endif

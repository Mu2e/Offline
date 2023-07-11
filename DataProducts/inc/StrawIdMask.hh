#ifndef DataProducts_StrawIdMask_hh
#define DataProducts_StrawIdMask_hh
//
//  This is a helper class that, conjunction with StrawId, allows matching tracker elements at various levels
//  Original author: Dave Brown (LBNL)
//
#include "Offline/DataProducts/inc/StrawId.hh"
#include <vector>
#include <string>
namespace mu2e {
  class StrawIdMask {
    // levels of StrawId to mask.  Note that 'straw' and 'panel' refer to ALL straws and panels with that value (ie
    // straw 18 of every panel in the detector.  uniquestraw and uniquepanel refer to individual straws and panels
    public:
    enum Level{none=-1,tracker,station,plane,panel,straw,uniquestraw,uniquepanel};
    // compute the mask associated with a given level
    static uint16_t levelMask(Level fval);
    // specify which levels to compare on construction.
    StrawIdMask(Level lev=none) : _level(lev), _mask(levelMask(_level)) {}
    // identity is obvious
    bool operator ==(StrawIdMask const& other ) const { return other._mask == _mask; } // must match exactly to be equal
    bool operator < (StrawIdMask const& other ) const { return other._level < _level; } // needed for sorting, no physical significance
    explicit StrawIdMask(std::string const& asstring); // valid values: 'plane', 'panel' (match ALL panels with the same value)
    // 'uniquepanel' (match exactly 1 panel), 'straw' (will match ALL straws with the same value), 'uniquestraw' (will match exactly 1 straw)
    // accessors
    Level level() const { return _level; }
    std::string levelName() const { return levelName(_level); }
    uint16_t mask() const { return _mask; }
    // compare StrawIds given this objects mask.  The 2 match if their masked levels are the same.
    // This is the main function of this class
    bool equal (StrawId const& sid1, StrawId const& sid2) const {
      return (level() != none) && (sid1.asUint16() & _mask) == (sid2.asUint16() & _mask); }
    bool notequal (StrawId const& sid1, StrawId const& sid2) const {
      return (level() == none) || (! equal(sid1,sid2)); }
    // return a truncated StrawID.  This is returned by value
    StrawId maskStrawId(StrawId const& sid) const {
      return StrawId(sid.asUint16() & _mask); }
    static std::string levelName(Level level);
    private:
    Level _level;
    uint16_t _mask;
  };
}
#endif


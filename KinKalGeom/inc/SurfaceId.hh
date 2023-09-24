//
//  Define identifiers for the reconstruction surfaces
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKalGeom_SurfaceId_hh
#define KinKalGeom_SurfaceId_hh
#include <map>
#include <vector>
#include <string>
#include "Offline/GeneralUtilities/inc/EnumToStringSparse.hh"

namespace mu2e {
  class SurfaceIdDetail {
    public:
      enum enum_type {
        unknown =-1,
        TT_Front=0, TT_Mid, TT_Back, TT_Inner, TT_Outer, // tracker VD equivalents
        DS_Front=80, DS_Back, DS_Inner, DS_Outer,
        IPA, OPA, TSDA, // Absorbers in the DS.  These are the inner surfaces for the OPA and TSDA
        ST_Front=100,ST_Back, ST_Inner, ST_Outer, ST_Foils, // stopping target bounding surfaces and foils
        TCRV=200 // CRV test planes
      };

    static std::string const& typeName();
    static std::map<enum_type,std::string> const& names();
  };
  using SurfaceIdEnum = EnumToStringSparse<SurfaceIdDetail>;
  class SurfaceId {
    public:
      using enum_type = SurfaceIdDetail::enum_type;
      // copy the constructors
      SurfaceId() : index_(0) {}
      SurfaceId(std::string const& name, int index=0) : sid_(name), index_(index) {}
      SurfaceId(enum_type sid, int index=0) : sid_(sid), index_(index) {}

      // forward some accessors
      auto const& id() const { return sid_; }
      int index() const { return index_; }
      auto const& name() const { return sid_.name(); }

      bool indexMatch(SurfaceId const& other) const { return index_ == other.index_ || index_ < 0 || other.index_ < 0; }
      bool indexCompare(SurfaceId const& other) const { return index_<0 || other.index_ < 0 ? false : index_ < other.index_; }
      bool operator == (SurfaceId const& other ) const { return sid_ == other.sid_ && indexMatch(other) ; }
      bool operator != (SurfaceId const& other ) const { return sid_ != other.sid_ || !indexMatch(other) ; }
      bool operator < (SurfaceId const& other ) const { return sid_ == other.sid_ ? indexCompare(other) : sid_ < other.sid_; }
    private:
      SurfaceIdEnum sid_;
      int index_; // index.  Negative value is a wild card
  };
  using SurfaceIdCollection = std::vector<SurfaceId>;
}
#endif

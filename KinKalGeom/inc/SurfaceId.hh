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
        ST_Front=100,ST_Back, ST_Inner, ST_Outer, ST_Foil, // stopping target bounding surfaces and foils
        TCRV_Plane=200 // CRV test plane
      };

    static std::string const& typeName();
    static std::map<enum_type,std::string> const& names();
  };
  using SurfaceIdEnum = EnumToStringSparse<SurfaceIdDetail>;
  class SurfaceId {
    public:
      // copy the constructors
      SurfaceId() : index_(0) {}
      SurfaceId(std::string const& name, unsigned index=0) : sid_(name), index_(index) {}
      SurfaceId(SurfaceIdDetail::enum_type sid, unsigned index=0) : sid_(sid), index_(index) {}

      bool operator == (SurfaceId const& other ) const { return sid_ == other.sid_ && index_ == other.index_; }
      bool operator != (SurfaceId const& other ) const { return sid_ != other.sid_ || index_ != other.index_; }
      bool operator < (SurfaceId const& other ) const { return sid_ == other.sid_ ? index_ < other.index_ : sid_ < other.sid_; }
    private:
      SurfaceIdEnum sid_;
      unsigned index_; // index
  };
  using SurfaceIdCollection = std::vector<SurfaceId>;
}
#endif

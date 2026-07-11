//
//  Define identifiers for surfaces used in simulation and reconstruction. These can be virtual or physical
//  original author: David Brown (LBNL) 2023
//
#ifndef DataProducts_SurfaceId_hh
#define DataProducts_SurfaceId_hh
#include <cstddef>
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
        IPA_Legacy, //introduced for backwards compatibility with MDS1
        DS_CryoInner, DS_CryoOuter, DS_ShieldInner, DS_ShieldOuter, DS_Coil,
        IPA=90, IPA_Front, IPA_Back,
        OPA=95, TSDA, // Absorbers in the DS
        ST_Front=100,ST_Back, ST_Inner, ST_Outer, ST_Foils, ST_Wires, // stopping target bounding surfaces and components
        TCRV=200, // CRV test plane (deprecated)
        CRV_EX=201, CRV_T1, CRV_T2, // keep values 201-203 for extracted geometry compatibility
        // CRV_T1 and CRV_T2 are Run 2 CRV-TS sectors 9 and 10 in crv_counters_v10.txt.
        // Run 2 CRV sectors (from crv_counters_v10.txt, prefixed CRV_ by CosmicRayShieldMaker)
        CRV_T3=204, CRV_T4, CRV_T5,                           // CRV-Top
        CRV_R1=210, CRV_R2, CRV_R3, CRV_R4, CRV_R5, CRV_R6,   // CRV-Right
        CRV_L1=220, CRV_L2, CRV_L3,                           // CRV-Left
        CRV_E1=230, CRV_E2,                                   // CRV-TS Extension
        CRV_U =240,                                           // CRV-Upstream
        CRV_D1=250, CRV_D2, CRV_D3, CRV_D4,                   // CRV-Downstream
        CRV_C1=260, CRV_C2,                                   // CRV-Cryo-Outer
        CRV_C3=262,                                           // legacy MDC2020/crv_counters_v09; retire once those datasets are no longer used
        CRV_M1=270, CRV_M2, CRV_M3, CRV_M4, CRV_M5, CRV_M6, CRV_M7, CRV_M8, // CRV-Muon-Taggers (Mu2e/Offline PR #1864)
        CRV_StrongBack=280,                                   // CRV module Al strongback (tracker-side support plate)
        DS_HatchConcrete=300                                  // detector-area hatch concrete block approximation
      };

    // Update this counter whenever you add/remove surface IDs from the enum above.
    static constexpr std::size_t nSurfaceIds = 63;

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
      int index_; // index.  Negative value is a wild card for matching
  };
  using SurfaceIdCollection = std::vector<SurfaceId>;
  // printout
  std::ostream& operator<<(std::ostream& ost, const SurfaceId& s );
}
#endif

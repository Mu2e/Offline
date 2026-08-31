#include <cstddef>
#include <type_traits>
#include <utility>

#include "Offline/DataProducts/inc/SurfaceId.hh"

using namespace std;

namespace mu2e {

  std::string const& SurfaceIdDetail::typeName() {
    static const std::string type("SurfaceIdEnum");
    return type;
  }

  namespace {
    using SurfaceIdName = std::pair<SurfaceIdEnum::enum_type, char const*>;

    constexpr SurfaceIdName surfaceIdNames[] = {
    std::make_pair(SurfaceIdEnum::unknown,   "unknown"  ),
    std::make_pair(SurfaceIdEnum::TT_Front, "TT_Front"),
    std::make_pair(SurfaceIdEnum::TT_Mid, "TT_Mid"),
    std::make_pair(SurfaceIdEnum::TT_Back, "TT_Back"),
    std::make_pair(SurfaceIdEnum::TT_Inner, "TT_Inner"),
    std::make_pair(SurfaceIdEnum::TT_Outer, "TT_Outer"),
    std::make_pair(SurfaceIdEnum::DS_Front, "DS_Front"),
    std::make_pair(SurfaceIdEnum::DS_Back, "DS_Back"),
    std::make_pair(SurfaceIdEnum::DS_Inner, "DS_Inner"),
    std::make_pair(SurfaceIdEnum::DS_Outer, "DS_Outer"),
    std::make_pair(SurfaceIdEnum::IPA_Legacy, "IPA_Legacy"),
    std::make_pair(SurfaceIdEnum::DS_CryoInner, "DS_CryoInner"),
    std::make_pair(SurfaceIdEnum::DS_CryoOuter, "DS_CryoOuter"),
    std::make_pair(SurfaceIdEnum::DS_ShieldInner, "DS_ShieldInner"),
    std::make_pair(SurfaceIdEnum::DS_ShieldOuter, "DS_ShieldOuter"),
    std::make_pair(SurfaceIdEnum::DS_Coil, "DS_Coil"),
    std::make_pair(SurfaceIdEnum::IPA, "IPA"),
    std::make_pair(SurfaceIdEnum::IPA_Front, "IPA_Front"),
    std::make_pair(SurfaceIdEnum::IPA_Back, "IPA_Back"),
    std::make_pair(SurfaceIdEnum::OPA, "OPA"),
    std::make_pair(SurfaceIdEnum::TSDA, "TSDA"),
    std::make_pair(SurfaceIdEnum::ST_Front, "ST_Front"),
    std::make_pair(SurfaceIdEnum::ST_Back, "ST_Back"),
    std::make_pair(SurfaceIdEnum::ST_Inner, "ST_Inner"),
    std::make_pair(SurfaceIdEnum::ST_Outer, "ST_Outer"),
    std::make_pair(SurfaceIdEnum::ST_Foils, "ST_Foils"),
    std::make_pair(SurfaceIdEnum::ST_Wires, "ST_Wires"),
    std::make_pair(SurfaceIdEnum::TCRV, "TCRV"),
    std::make_pair(SurfaceIdEnum::CRV_EX, "CRV_EX"),
    std::make_pair(SurfaceIdEnum::CRV_T1, "CRV_T1"),
    std::make_pair(SurfaceIdEnum::CRV_T2, "CRV_T2"),
    std::make_pair(SurfaceIdEnum::CRV_T3, "CRV_T3"),
    std::make_pair(SurfaceIdEnum::CRV_T4, "CRV_T4"),
    std::make_pair(SurfaceIdEnum::CRV_T5, "CRV_T5"),
    std::make_pair(SurfaceIdEnum::CRV_R1, "CRV_R1"),
    std::make_pair(SurfaceIdEnum::CRV_R2, "CRV_R2"),
    std::make_pair(SurfaceIdEnum::CRV_R3, "CRV_R3"),
    std::make_pair(SurfaceIdEnum::CRV_R4, "CRV_R4"),
    std::make_pair(SurfaceIdEnum::CRV_R5, "CRV_R5"),
    std::make_pair(SurfaceIdEnum::CRV_R6, "CRV_R6"),
    std::make_pair(SurfaceIdEnum::CRV_L1, "CRV_L1"),
    std::make_pair(SurfaceIdEnum::CRV_L2, "CRV_L2"),
    std::make_pair(SurfaceIdEnum::CRV_L3, "CRV_L3"),
    std::make_pair(SurfaceIdEnum::CRV_E1, "CRV_E1"),
    std::make_pair(SurfaceIdEnum::CRV_E2, "CRV_E2"),
    std::make_pair(SurfaceIdEnum::CRV_U,  "CRV_U"),
    std::make_pair(SurfaceIdEnum::CRV_D1, "CRV_D1"),
    std::make_pair(SurfaceIdEnum::CRV_D2, "CRV_D2"),
    std::make_pair(SurfaceIdEnum::CRV_D3, "CRV_D3"),
    std::make_pair(SurfaceIdEnum::CRV_D4, "CRV_D4"),
    std::make_pair(SurfaceIdEnum::CRV_C1, "CRV_C1"),
    std::make_pair(SurfaceIdEnum::CRV_C2, "CRV_C2"),
    std::make_pair(SurfaceIdEnum::CRV_C3, "CRV_C3"),
    std::make_pair(SurfaceIdEnum::CRV_C4, "CRV_C4"),
    std::make_pair(SurfaceIdEnum::CRV_M1, "CRV_M1"),
    std::make_pair(SurfaceIdEnum::CRV_M2, "CRV_M2"),
    std::make_pair(SurfaceIdEnum::CRV_M3, "CRV_M3"),
    std::make_pair(SurfaceIdEnum::CRV_M4, "CRV_M4"),
    std::make_pair(SurfaceIdEnum::CRV_M5, "CRV_M5"),
    std::make_pair(SurfaceIdEnum::CRV_M6, "CRV_M6"),
    std::make_pair(SurfaceIdEnum::CRV_M7, "CRV_M7"),
    std::make_pair(SurfaceIdEnum::CRV_M8, "CRV_M8"),
    std::make_pair(SurfaceIdEnum::CRV_StrongBack, "CRV_StrongBack"),
    std::make_pair(SurfaceIdEnum::DS_HatchConcrete, "DS_HatchConcrete")
    };

    constexpr std::size_t nSurfaceIdNames = sizeof(surfaceIdNames)/sizeof(surfaceIdNames[0]);
    static_assert(nSurfaceIdNames == SurfaceIdDetail::nSurfaceIds,
        "SurfaceId enum and name map must stay in sync");
  }

  static const std::map<SurfaceIdEnum::enum_type,std::string> nam(
      surfaceIdNames, surfaceIdNames + nSurfaceIdNames);

  std::map<SurfaceIdEnum::enum_type,std::string> const& SurfaceIdDetail::names(){
    return nam;
  }

  std::ostream& operator<<(std::ostream& ost, const SurfaceId& s ) {
    ost << s.id() << ":" << s.index();
    return ost;
  }

}

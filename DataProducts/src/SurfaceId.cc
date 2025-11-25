#include <type_traits>
#include <utility>

#include "Offline/DataProducts/inc/SurfaceId.hh"

using namespace std;

namespace mu2e {

  std::string const& SurfaceIdDetail::typeName() {
    static const std::string type("SurfaceIdEnum");
    return type;
  }

  static const std::map<SurfaceIdEnum::enum_type,std::string> nam{
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
    std::make_pair(SurfaceIdEnum::TCRV, "TCRV"),

    std::make_pair(SurfaceIdEnum::EMC_Disk_0_SurfIn, "EMC_Disk_0_SurfIn"),
    std::make_pair(SurfaceIdEnum::EMC_Disk_0_SurfOut, "EMC_Disk_0_SurfOut"),
    std::make_pair(SurfaceIdEnum::EMC_Disk_1_SurfIn, "EMC_Disk_1_SurfIn"),
    std::make_pair(SurfaceIdEnum::EMC_Disk_1_SurfOut, "EMC_Disk_1_SurfOut"),
    std::make_pair(SurfaceIdEnum::EMC_Disk_0_EdgeIn, "EMC_Disk_0_EdgeIn"),
    std::make_pair(SurfaceIdEnum::EMC_Disk_0_EdgeOut "EMC_Disk_0_EdgeOut"),
    std::make_pair(SurfaceIdEnum::EMC_Disk_1_EdgeIn, "EMC_Disk_1_EdgeIn"),
    std::make_pair(SurfaceIdEnum::EMC_Disk_1_EdgeOut, "EMC_Disk_1_EdgeOut"),

    std::make_pair(SurfaceIdEnum::EMC_0_FrontIn, "EMC_0_FrontIn"),
    std::make_pair(SurfaceIdEnum::EMC_0_FrontOut, "EMC_0_FrontOut"),
    std::make_pair(SurfaceIdEnum::EMC_1_FrontIn, "EMC_1_FrontIn"),
    std::make_pair(SurfaceIdEnum::EMC_1_FrontOut, "EMC_1_FrontOut"),
    std::make_pair(SurfaceIdEnum::EMC_0_FrontIn, "EMC_2_FrontIn"),
    std::make_pair(SurfaceIdEnum::EMC_0_FrontOut, "EMC_2_FrontOut"),
    std::make_pair(SurfaceIdEnum::EMC_1_FrontIn, "EMC_3_FrontIn"),
    std::make_pair(SurfaceIdEnum::EMC_1_FrontOut, "EMC_3_FrontOut"),

  };
  std::map<SurfaceIdEnum::enum_type,std::string> const& SurfaceIdDetail::names(){
    return nam;
  }

  std::ostream& operator<<(std::ostream& ost, const SurfaceId& s ) {
    ost << s.id() << ":" << s.index();
    return ost;
  }

}

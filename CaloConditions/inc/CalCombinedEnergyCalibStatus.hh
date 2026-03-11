#ifndef CaloConditions_CalCombinedEnergyCalibStatus_hh
#define CaloConditions_CalCombinedEnergyCalibStatus_hh

//
// Status codes for combined calorimeter energy calibrations.
//
// 0   - updated: cosmic + source, consistent
// 1   - fallback: methods inconsistent, kept old value
// 2   - fallback: all methods statistically invalid
// 101 - updated using cosmic only
// 102 - updated using source only
//

namespace mu2e {

namespace CalCombinedEnergyCalibStatus {

constexpr int Updated = 0;
constexpr int FallbackInconsistentMethods = 1;
constexpr int FallbackAllMethodsInvalid = 2;
constexpr int UpdatedUsingCosmicOnly = 101;
constexpr int UpdatedUsingSourceOnly = 102;

inline const char* description(int statusCode) {
  switch (statusCode) {
    case Updated:
      return "updated: cosmic + source, consistent";
    case FallbackInconsistentMethods:
      return "fallback: methods inconsistent, kept old value";
    case FallbackAllMethodsInvalid:
      return "fallback: all methods statistically invalid";
    case UpdatedUsingCosmicOnly:
      return "updated using cosmic only";
    case UpdatedUsingSourceOnly:
      return "updated using source only";
    default:
      return "unknown status";
  }
}

}  // namespace CalCombinedEnergyCalibStatus

}  // namespace mu2e

#endif

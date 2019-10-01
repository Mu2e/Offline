#ifndef TrackerConditions_Mu2eMaterialConfig_hh
#define TrackerConditions_Mu2eMaterialConfig_hh
//
// Initialize Mu2eMaterial from fcl
//
#include <string>
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {

  struct Mu2eMaterialConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0,1,2")};
    fhicl::Atom<std::string> elements{
      Name("elements"), Comment("elements") };
    fhicl::Atom<std::string> isotopes{
      Name("isotopes"), Comment("isotopes") };
    fhicl::Atom<std::string> materials{
      Name("materials"), Comment("materials") };
    fhicl::Atom<std::string> strawGasMaterialName{
      Name("strawGasMaterialName"), Comment("strawGasMaterialName") };
    fhicl::Atom<std::string> strawWallMaterialName{
      Name("strawWallMaterialName"), Comment("strawWallMaterialName") };
    fhicl::Atom<std::string> strawWireMaterialName{
      Name("strawWireMaterialName"), Comment("strawWireMaterialName") };
    fhicl::Atom<double> dahlLynchScatteringFraction{
      Name("dahlLynchScatteringFraction"), Comment("dahlLynchScatteringFraction") };
    fhicl::Atom<double> intersectionTolerance{
      Name("intersectionTolerance"), Comment("intersectionTolerance") };
    fhicl::Atom<double> strawElementOffset{
      Name("strawElementOffset"), Comment("strawElementOffset") };
    fhicl::Atom<double> maximumIntersectionRadiusFraction{
      Name("maximumIntersectionRadiusFraction"), Comment("maximumIntersectionRadiusFraction") };

  };

}

#endif

#ifndef SimulationConfig_BookkeeperConfig_hh
#define SimulationConfig_BookkeeperConfig_hh

//
// Configuration for the MVACatalog ProditionsEntitiy
//

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {
  struct BookkeeperConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity: 0 or 1")}; 
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")}; 
    fhicl::Atom<double> muBeamPerPOT{
      Name("muBeamPerPOT"), Comment("muBeamPerPOT")};
    fhicl::Atom<double> flashPerMuBeam{
      Name("flashPerMuBeam"), Comment("flashPerMuBeam")};
  };
}

#endif

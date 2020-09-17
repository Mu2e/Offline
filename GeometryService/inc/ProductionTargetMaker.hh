#ifndef PRODUCTIONTARGETMAKER_HH
#define PRODUCTIONTARGETMAKER_HH

#include <memory>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace mu2e { class SimpleConfig; }
namespace mu2e { 
  class ProductionTarget; 
  class ProductionTargetMu2eII;
}

namespace mu2e {
  class ProductionTargetMaker {
  public:


    static std::unique_ptr<ProductionTarget>       make(const SimpleConfig& config, double solenoidOffset);
    static std::unique_ptr<ProductionTargetMu2eII> makeMu2eII(const SimpleConfig& config, double solenoidOffset);
    static const int tier1{1};
    // version 2 is the low density hayman
    static const int hayman_v_2_0{3};
    //enums for identifying Mu2eII target configurations, starting from hayman_v_2_0
    enum {mu2eii_conveyor = 4, mu2eii_rotating = 5};
 
    static std::unique_ptr<ProductionTarget> makeTier1(const SimpleConfig& config, double solenoidOffset);
    static std::unique_ptr<ProductionTarget> makeHayman_v_2_0(const SimpleConfig& config, double solenoidOffset);

    //different target object for Mu2e-II designs
    static std::unique_ptr<ProductionTargetMu2eII> makeMu2eIIConveyor(const SimpleConfig& config, double solenoidOffset);
    static std::unique_ptr<ProductionTargetMu2eII> makeMu2eIIRotating(const SimpleConfig& config, double solenoidOffset);



  };
}
#endif/*PRODUCTIONTARGETMAKER_HH*/

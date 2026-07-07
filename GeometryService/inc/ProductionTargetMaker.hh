#ifndef PRODUCTIONTARGETMAKER_HH
#define PRODUCTIONTARGETMAKER_HH

#include <memory>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <functional>

namespace mu2e { class SimpleConfig; }
namespace mu2e { class ProductionTarget; }

namespace mu2e {
  class ProductionTargetMaker {
  public:

    // Registry map for model makers — allows scaling without modifying make()
    using MakerFunction = std::function<std::unique_ptr<ProductionTarget>(const SimpleConfig&, double)>;
    static const std::map<std::string, MakerFunction>& getMakerRegistry();

    static std::unique_ptr<ProductionTarget> make(const SimpleConfig& config, double solenoidOffset);
    static const int tier1{1};
    // version 2 is the low density hayman
    static const int hayman_v_2_0{3};
    static const int stickman_v_1_0{4};


    static std::unique_ptr<ProductionTarget> makeTier1(const SimpleConfig& config, double solenoidOffset);
    static std::unique_ptr<ProductionTarget> makeHayman_v_2_0(const SimpleConfig& config, double solenoidOffset);
    static std::unique_ptr<ProductionTarget> makeStickman_v_1_0(const SimpleConfig& config, double solenoidOffset);



  };
}
#endif/*PRODUCTIONTARGETMAKER_HH*/

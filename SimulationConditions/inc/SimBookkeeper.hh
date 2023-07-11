#ifndef SimulationConditions_SimBookkeeper_hh
#define SimulationConditions_SimBookkeeper_hh

//
// SimBookkeeper is the ProditionsEntry for
// simulation proditions
//
// As of Oct 2020 this includes:
// 1) simulation stage efficiencies (DbTable: SimEfficiencies2)
// and nothing else
//

// Mu2e includes
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

  class SimBookkeeper : public ProditionsEntity {

  public:
    constexpr static const char* cxname = {"SimBookkeeper"};

    SimBookkeeper() : ProditionsEntity(cxname) {}
    // accessors
    double const getEff(const std::string& name) const {
      for (const auto& i_eff : _effs) {
        if (i_eff.first == name) {
          return i_eff.second;
        }
      }
      return -1;
    }

    // setters
    void addEff(std::string name, double new_val) {
      _effs[name] = new_val;
    }

    const std::string print() const;
    void print(std::ostream& os) const;

    // typedefs
    typedef std::shared_ptr<SimBookkeeper> ptr_t;
    typedef std::shared_ptr<const SimBookkeeper> cptr_t;

  private:
    // data
    std::map<std::string, double> _effs;
  };
}

#endif

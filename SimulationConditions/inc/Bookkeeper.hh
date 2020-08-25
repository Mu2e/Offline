#ifndef SimulationConditions_Bookkeeper_hh
#define SimulationConditions_Bookkeeper_hh

//
// Bookkeeper is the ProditionsEntry for
//

// Mu2e includes
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

  class Bookkeeper : virtual public ProditionsEntity {

  public:
    Bookkeeper() : _name("Bookkeeper") {}
    // accessors
    std::string const& name() const { return _name; }
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

    const std::string print() const {
      std::stringstream out;
      print(out);
      return out.str();
    }

    void print(std::ostream& os) const {
      os << "Efficiencies in " << _name << ":" << std::endl;
      for (const auto& i_eff : _effs) {
        os << i_eff.first << " = " << i_eff.second << std::endl;
      }
      os << std::endl;
    }

    // typedefs
    typedef std::shared_ptr<Bookkeeper> ptr_t;
    typedef std::shared_ptr<const Bookkeeper> cptr_t;

  private:
    // data
    std::string _name;
    std::map<std::string, double> _effs;
  };
}

#endif

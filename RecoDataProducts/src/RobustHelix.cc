#include "RecoDataProducts/inc/RobustHelix.hh"
namespace mu2e {
  std::ostream& operator<<(std::ostream& os, mu2e::RobustHelix const& helix) {
    os << " RobustHelix Center = " << helix.center() 
      << " Radius = " << helix.radius()
      << " Lambda = " << helix.lambda() << " PhiZ0 = " << helix.fz0(); return os; }
}

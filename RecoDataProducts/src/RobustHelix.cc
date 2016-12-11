#include "RecoDataProducts/inc/RobustHelix.hh"
namespace mu2e {
  std::ostream& operator<<(std::ostream& os, mu2e::RobustHelix const& helix) {
    os << " RobustHelix Center Radius = " << helix.rcent() << " Center Azimuth " << helix.fcent()
      << " Radius = " << helix.radius()
      << " Lambda = " << helix.lambda() << " PhiZ0 = " << helix.fz0(); return os; }
}

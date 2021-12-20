
#include "Offline/GlobalConstantsService/inc/ParticleData.hh"

namespace mu2e {

  std::ostream& operator<<(std::ostream& output, const ParticleData& pd) {
    output << std::setw(11) << pd.id();
    output << std::setw(22) << pd.name();
    output << std::setw(35) << pd.codeName();
    output << std::setw(10) << pd.charge();
    output << std::setw(12) << pd.mass();
    output << std::setw(12) << pd.lifetime();
    output << std::endl;
    return output;
  }

  void ParticleData::print(std::ostream& ostr) const { 
    ostr << *this; 
  }

} // namespace mu2e

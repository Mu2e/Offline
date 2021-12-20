
#include "Offline/GlobalConstantsService/inc/ParticleData.hh"

namespace mu2e {

std::ostream& operator<<( std::ostream& output, 
                          const ParticleData& pd ) {
  output << std::setw(11) << pd._id;
  output << std::setw(22) << pd._name;
  output << std::setw(35) << pd._codeName;
  output << std::setw(10) << pd._charge;
  output << std::setw(12) << pd._mass;
  output << std::setw(12) << pd._lifetime;
  output << std::endl;
  return output;            
}

}

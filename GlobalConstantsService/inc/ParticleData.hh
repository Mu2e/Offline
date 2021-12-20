#ifndef GlobalConstantsService_ParticleData_hh
#define GlobalConstantsService_ParticleData_hh
//
// Partcile data loaded from ParticleDataList
//
// the content is
//  - PDG ID
//  - displayname (includes special chars like +-*)
//  - progamming name (includes only alphanumeric and _)
//        this should be the same as  DataProducts/inc/PDGCode.hh
//  - charge (units of proton charge)
//  - mass (MeV)
//  - lifetime (ns)
//

#include <iostream>
#include <iomanip>

#include "Offline/GlobalConstantsService/inc/ParticleData.hh"

namespace mu2e {


  class ParticleData {

  public:

    ParticleData(int id, const std::string& name, const std::string& codeName,
                 float charge, float mass, float lifetime): 
      _id(id),_name(name),_codeName(codeName),_charge(charge),
      _mass(mass),_lifetime(lifetime) {}

    int id() const { return _id; }
    const std::string& name() const { return _name; }
    const std::string& displayName() const { return _name; }
    const std::string& codeName() const { return _codeName; }
    float charge() const { return _charge; }
    float mass() const { return _mass; }
    float lifetime() const { return _lifetime; }

    void print( std::ostream& ostr=std::cout) const {
      ostr << *this;
    }

    friend std::ostream& operator<<( std::ostream& output, 
                                     const ParticleData& pd );

  private:
    int _id;
    std::string _name;
    std::string _codeName;
    float _charge;
    float _mass;
    float _lifetime;

  };  // ParticleData

} //end namespace mu2e

#endif /* GlobalConstantsService_ParticleData_hh */

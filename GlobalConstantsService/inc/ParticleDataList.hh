#ifndef GlobalConstantsService_ParticleDataList_hh
#define GlobalConstantsService_ParticleDataList_hh
//
// A simple particle data list
//
// this list was created from 
// - HepPDT list
// - mass_width_2008.mc
// - g4nuclei.tbl
// - aspects of the geant particle list
// - aspects of DataProducts/inc/PDGCode.hh
//
// Units are standard units: proton_charge, MeV, and ns
//
// the text file columns and accessors are
//  1) PDG ID
//  2) displayname (includes special chars like +-*)
//  3) code name (includes only alphanumeric and _)
//         this should be the same as  DataProducts/inc/PDGCode.hh
//  4) alias - another common name that can be used to look up the particle
//  5) charge (units of proton charge)
//  6) mass (MeV)
//  7) lifetime (ns)
//


#include <string>
#include <iostream>
#include <map>

#include "Offline/Mu2eInterfaces/inc/ConditionsEntity.hh"
#include "Offline/GlobalConstantsService/inc/ParticleData.hh"

namespace mu2e {

  class SimpleConfig;

  class ParticleDataList : virtual public ConditionsEntity{

  public:

    ParticleDataList( SimpleConfig const& config );

    //lookup by PDG ID
    const ParticleData& particle( int id ) const;
    // lookup by name - the display name, code name or alias
    const ParticleData& particle( std::string const& name ) const;

    // the entire table as a map with pdgID as key
    std::map<int,ParticleData> const& list() const { return _list;}
    // the name lookup map with name as key and PDG id as mapped_type
    std::map<std::string,int> const& names() const { return _names;}

    void printTable( std::ostream& ostr=std::cout) const;

  private:

    // The actual particle data table.
    // we add an entry if the particle id corresponds to a nonexistent nuclei
    mutable std::map<int,ParticleData> _list;

    // map for lookup by name
    mutable std::map<std::string,int> _names;

  };  // ParticleDataList

} //end namespace mu2e

#endif /* GlobalConstantsService_ParticleDataList_hh */

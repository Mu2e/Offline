#ifndef Mu2e_ParticleDataTable_hh
#define Mu2e_ParticleDataTable_hh
//
// Mu2e wrapper around HepPDT::ParticleDataTable 
//
//   $Id: ParticleDataTable.hh,v 1.1 2010/03/19 01:18:41 kutschke Exp $
//   $Author: kutschke $
//   $Date: 2010/03/19 01:18:41 $
//
//  Original author Rob Kutschke
//
// Notes:
// 1) For each particle the following info is stored:
//      pdgid, name, mass, width, lifetime. 
//     plus errors on the mass and width. See HepPDT::ParticleData.
//
// 2) Units are non-standard for Mu2e:
//      mass in GeV
//      lifetimes in seconds
//
// 3) Why do we have a wrapper?  
//      a) For safety.  See bullet 4.
//      b) So that this class can be managed by the Conditions system, this class
//         needs to inherit from ConditionsEntity.
//
// 4) What is the safety concern?
//    Some of the HepPDT methods to access ParticleData objects
//    will, if the requested particle cannot be found, return
//    a null pointer.   This leaves it the user's job to check
//    for a null pointer.  The wrapper will do the check
//    and will throw if the pointer is null.  We can add an non-checking
//    version if it turns out to be important for speed.
//
// 5) Why not use the ProcessUnknownID method to do the throw?
//    Because that method is not called for the accessors that
//    take a particle name as an argument.  So using this would
//    leave those accessors unprotected.
//
// 6) Why not use inheritance and then override only the unsafe accessors?
//    HepPDT::ParticleDataTable does not have a virtual destructor.
//
// 7) I removed all of the accessors and iterators that returned non-const things.
//


// Mu2e includes.
#include "ConditionsService/inc/ConditionsEntity.hh"

// External includes.
#include "HepPDT/ParticleDataTable.hh"

namespace mu2e {

  // Forward declarations
  class SimpleConfig;

  class ParticleDataTable : public ConditionsEntity{
  
  public:

    typedef HepPDT::ParticleDataTable::PDTMap::const_iterator     const_iterator;
    typedef HepPDT::ParticleDataTable::PDTNameMap::const_iterator const_iteratorByName;

    ParticleDataTable( SimpleConfig const& config );
    ParticleDataTable( std::string const& name, std::string const& tableFilename );

    // Accept the compiler supplied destructor.

    /// Access particle information via ParticleID or particle name
    HepPDT::ParticleData const& particle( HepPDT::ParticleID ) const;
    HepPDT::ParticleData const& particle( std::string const& name ) const;

    // Duplicate accessors with [] syntax.
    HepPDT::ParticleData const& operator[] ( HepPDT::ParticleID id ) const{
      return particle(id);
    }
    HepPDT::ParticleData const& operator[] ( std::string const& name ) const{
      return particle(name);
    }

    /// Size of the particle data table and iterators over it.
    int             size()  const { return pdt.size(); }
    const_iterator  begin() const { return pdt.begin(); }
    const_iterator  end()   const { return pdt.end(); }

    /// Size of the map of particle names and iterators over it.
    int                   sizeNameMap()  const { return pdt.sizeNameMap(); }
    const_iteratorByName  beginNameMap() const { return pdt.beginNameMap(); }
    const_iteratorByName  endNameMap()   const { return pdt.endNameMap(); }

    /// Return the name of this particle data table
    //  Cannot return const& since pdt returns by value.
    std::string tableName() const { return pdt.tableName(); }

    // Access the table directly to get at other functions that are not forwarded.
    HepPDT::ParticleDataTable const& table() const { return pdt;}

  private:

    // The actual particle data table.
    HepPDT::ParticleDataTable pdt;

    // The name of the file from which the data was loaded.
    std::string _tableFilename;

    // ---  copying; forbidden:
    ParticleDataTable( const HepPDT::ParticleDataTable & orig );
    ParticleDataTable& operator=( const HepPDT::ParticleDataTable & );

    // A helper function to load the table from a file.
    void loadTableFromFile();

  };  // ParticleDataTable

} //end namespace mu2e

#endif

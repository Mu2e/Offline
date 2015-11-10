#ifndef GlobalConstantsService_ParticleDataTable_hh
#define GlobalConstantsService_ParticleDataTable_hh
//
// Mu2e wrapper around HepPDT::ParticleDataTable
//
//  Original author Rob Kutschke
//
// Notes:
// 1) For each particle the following info is stored:
//      pdgid, name, mass, width
//     plus errors on the mass and width. See HepPDT::ParticleData.
//     Lifetimes are computed from the stored width.
//
// 2) In the native file format, units are non-standard for Mu2e:
//      - masses and widths in GeV
//      - lifetimes in seconds
//    So convert the masses and widths to MeV.  This converts
//    The lifetimes to kilo-seconds!
//    We are waiting for a change from HepPDT to be able to get the
//    lifetimes in ns.
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
#include "Mu2eInterfaces/inc/ConditionsEntity.hh"
#include "cetlib/maybe_ref.h"

// External includes.
#include "HepPDT/ParticleDataTable.hh"

namespace mu2e {

  // Forward declarations
  class SimpleConfig;

  class ParticleDataTable : virtual public ConditionsEntity{

  public:

    // Return type of accessors
    typedef cet::maybe_ref<HepPDT::ParticleData const>            maybe_ref;

    typedef HepPDT::ParticleDataTable::PDTMap::const_iterator     const_iterator;
    typedef HepPDT::ParticleDataTable::PDTNameMap::const_iterator const_iteratorByName;

    ParticleDataTable( SimpleConfig const& config );
    ParticleDataTable( std::string const& name, std::string const& tableFilename );

    // Accept the compiler supplied destructor.  Copying forbidden - see below.

    /// Access particle information via ParticleID or particle name
    maybe_ref particle( HepPDT::ParticleID const & id ) const ;
    maybe_ref particle( std::string const& name ) const;

    // Duplicate accessors with [] syntax.
    maybe_ref operator[] ( HepPDT::ParticleID id ) const{
      return particle(id);
    }
    maybe_ref operator[] ( std::string const& name ) const{
      return particle(name);
    }

    /// Size of the particle data table and iterators over it.
    int             size()  const { return _pdt.size(); }
    const_iterator  begin() const { return _pdt.begin(); }
    const_iterator  end()   const { return _pdt.end(); }

    /// Size of the map of particle names and iterators over it.
    int                   sizeNameMap()  const { return _pdt.sizeNameMap(); }
    const_iteratorByName  beginNameMap() const { return _pdt.beginNameMap(); }
    const_iteratorByName  endNameMap()   const { return _pdt.endNameMap(); }

    /// Return the name of this particle data table
    //  Cannot return const& since pdt returns by value.
    std::string tableName() const { return _pdt.tableName(); }

    // Access the table directly to get at other functions that are not forwarded.
    HepPDT::ParticleDataTable const& table() const { return _pdt;}

    //    void printTable( std::ostream& ostr);

  private:

    // The actual particle data table.
    // we shall add an entry if the particle id corresponds to a nonexistent nuclei
    mutable HepPDT::ParticleDataTable _pdt;

    // The name of the file from which the data was loaded.
    std::string _tableFilename;

    // Tempoary hack.  The name of the auxillary particle data file,
    // which has better values for the masses and widths but which is
    // missing anti-particles.
    std::string _auxillaryFilename;

    // The name of the particle data file based on geant4 nuclei data
    std::string _geant4PDTFilename;

    // Keep track if the units were changed or not.
    bool _unitsChanged;

    // variable controlling the amount of printout
    int  _verbosityLevel;

    // ---  copying; forbidden:
    ParticleDataTable( const HepPDT::ParticleDataTable & orig );
    ParticleDataTable& operator=( const HepPDT::ParticleDataTable & );

    // A helper function to load the table from a file.
    void loadTableFromFile();

    // Change the units of the masses and widths from GeV to MeV.
    // This changes the units of lifetimes to kilo-seconds!
    void changeUnits();

    // Improve masses and lifetimes.
    void improveData();

    // Add geant4 nuclei data
    void addGeant4Data();

  };  // ParticleDataTable

} //end namespace mu2e

#endif /* GlobalConstantsService_ParticleDataTable_hh */

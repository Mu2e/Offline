//
// Mu2e wrapper around HepPDT::ParticleDataTable 
//
//   $Id: ParticleDataTable.cc,v 1.1 2010/03/19 01:18:41 kutschke Exp $
//   $Author: kutschke $
//   $Date: 2010/03/19 01:18:41 $
//


#include <iostream>
#include <fstream>

// Framework include files.
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e include files
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

// External include files.
#include "HepPDT/TableBuilder.hh"

using namespace std;

namespace mu2e {

  // Default construct an empty table and then fill it.
  ParticleDataTable::ParticleDataTable( std::string const& name,
                                        std::string const& tableFilename ):
    pdt(name),
    _tableFilename(tableFilename){

    loadTableFromFile();
  }

  ParticleDataTable::ParticleDataTable( SimpleConfig const& config ):
    pdt("Mu2eParticleData"){
    
    _tableFilename = config.getString("particleDataTable.filename",
                                      "ConditionsService/data/mass_width_2008.mc");    
    
    loadTableFromFile();
  }

  void ParticleDataTable::loadTableFromFile(){

    // Construct the table builder.
    HepPDT::TableBuilder  tb(pdt);

    // Build the table from the data file.
    ifstream in(_tableFilename.c_str());
    if ( !in ) {
      throw cms::Exception("FILE")
        << "Unable to open particle data file: " 
        << _tableFilename << "\n";
    }
    HepPDT::addPDGParticles( in, tb);

  }

  // Accessor by ID that checks for null pointer.
  HepPDT::ParticleData const& ParticleDataTable::particle( HepPDT::ParticleID id ) const{
    HepPDT::ParticleData const* p = pdt.particle(id);
    if ( p == 0 ){
      throw cms::Exception("RANGE")
        << "Could not find requested particle in the ParticleDataTable.  " 
        << "Requested paricle id code was: "
        << id.pid() << "\n";
    }
    return *p;
  }

  // Accessor by name that checks for null pointer.
  HepPDT::ParticleData const& ParticleDataTable::particle( std::string const& name ) const{
    HepPDT::ParticleData const* p = pdt.particle(name);
    if ( p == 0 ){
      throw cms::Exception("RANGE")
        << "Could not find requested particle in the ParticleDataTable.  " 
        << "Requested particle name was: "
        << name << "\n";
    }
    return *p;
  }

}  // end namespace mu2e

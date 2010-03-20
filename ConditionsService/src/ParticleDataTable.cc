//
// Mu2e wrapper around HepPDT::ParticleDataTable 
//
//   $Id: ParticleDataTable.cc,v 1.2 2010/03/20 00:59:39 kutschke Exp $
//   $Author: kutschke $
//   $Date: 2010/03/20 00:59:39 $
//
//
// 1) In changeUnits I want to loop over all elements in the table and change some
//    units, in situ. The minor complication is that there is no non-const iterator
//    over the table.  But we can stitch one together using two components.  First,
//    there is an accessor by ParticleID to give a non-const pointer to any element
//    in the table.  Also, there is a const iterator over the table.
//
// 2) It appears that data in the table is not guarantteed to be accessible until 
//    the builder goes out of scope.  However it appears that sometimes it is 
//    possible to access the table before that.
//
#include <fstream>

// Framework include files.
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e include files
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"

// External include files.
#include "HepPDT/TableBuilder.hh"
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

using HepPDT::ParticleData;
using HepPDT::ParticleID;
using HepPDT::Measurement;
using CLHEP::GeV;
using CLHEP::second;

namespace mu2e {

  // Default construct an empty table and then fill it.
  ParticleDataTable::ParticleDataTable( std::string const& name,
                                        std::string const& tableFilename ):
    _pdt(name),
    _tableFilename(tableFilename),
    _unitsChanged(false){

    loadTableFromFile();
  }

  ParticleDataTable::ParticleDataTable( SimpleConfig const& config ):
    _pdt("Mu2eParticleData"),
    _unitsChanged(false){
    
    _tableFilename = config.getString("particleDataTable.filename",
                                      "ConditionsService/data/mass_width_2008.mc");    
    
    loadTableFromFile();
  }

  void ParticleDataTable::loadTableFromFile(){

    // See note 2.
    {
      // Construct the table builder.
      HepPDT::TableBuilder  tb(_pdt);

      // Build the table from the data file.
      ifstream in(_tableFilename.c_str());
      if ( !in ) {
        throw cms::Exception("FILE")
          << "Unable to open particle data file: " 
          << _tableFilename << "\n";
      }
      HepPDT::addPDGParticles( in, tb);
    }

    // Make sure masses and widths are in MeV.
    changeUnits();
  }

  // Accessor by ID that checks for null pointer.
  ParticleData const& ParticleDataTable::particle( ParticleID id ) const{
    ParticleData const* p = _pdt.particle(id); 
    if ( p == 0 ){
      throw cms::Exception("RANGE")
        << "Could not find requested particle in the ParticleDataTable.  " 
        << "Requested paricle id code was: "
        << id.pid() << "\n";
    }
    return *p;
  }

  // Accessor by name that checks for null pointer.
  ParticleData const& ParticleDataTable::particle( std::string const& name ) const{
    ParticleData const* p = _pdt.particle(name);
    if ( p == 0 ){
      throw cms::Exception("RANGE")
        << "Could not find requested particle in the ParticleDataTable.  " 
        << "Requested particle name was: "
        << name << "\n";
    }
    return *p;
  }

  void ParticleDataTable::changeUnits(){

    // Electron mass, in MeV.
    double eMassMeV = 0.510999;

    // Electron mass from the table.
    double eMass = particle(PDGCode::e_minus).mass().value();

    // Ratio of two measures of the mass.
    double r = eMass/eMassMeV;

    // Tolerance on the ratio.
    double tolerance = 1.e-7;

    // Decide if the table is in MeV, GeV or something else?
    int units(0);
    if ( std::abs(r-1.0) < tolerance ){
      units = 1;
    } else if ( std::abs(r-0.001) < tolerance ){
      units = 2;
    }

    if ( units == 0 ){
      edm::LogWarning("CONDITIONS") 
        << "Did not recognize the units of masses in the particle data table.\n"
        << "The electron mass appears to be: "
        << eMass;
    }

    // Units are GeV, so change them to MeV.
    if ( units == 2 ){
      _unitsChanged = true;

      edm::LogWarning("CONDITIONS") 
        << "The HepPDT particle data table has masses in GeV. Changing to MeV.\n"
        << "  ( This leaves the lifetimes in a screwed up state: they are in kilo-seconds."
        << "    This does not affect Geant4 which has its own table of lifetimes. )";

      for ( HepPDT::ParticleDataTable::const_iterator i=_pdt.begin(), e=_pdt.end();
            i!=e; ++i ){

        // Get non-const reference to the particle data.  See Note 1.
        ParticleData const& tmp = i->second;
        ParticleData& particle = *_pdt.particle(tmp.ID());
          
        // Extract properties with dimensions of mass.
        Measurement mass   = particle.mass();
        Measurement width  = particle.totalWidth();
        double lowerCutoff = particle.lowerCutoff();
        double upperCutoff = particle.upperCutoff();

        // Compute same quantities in new units.
        Measurement newmass (  mass.value()*GeV,  mass.sigma()*GeV );
        Measurement newwidth( width.value()*GeV, width.sigma()*GeV );

        // Reset the properties in the table.  
        particle.setMass(newmass);
        particle.setTotalWidth(newwidth);
        particle.setLowerCutoff( lowerCutoff*GeV );
        particle.setUpperCutoff( upperCutoff*GeV );

      }
    }

  } // end changeUnits.

}  // end namespace mu2e

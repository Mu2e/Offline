#ifndef ToyDP_GenId_hh
#define ToyDP_GenId_hh
//
// An enum-matched-to-names class for generator Id's.
//
//
// $Id: GenId.hh,v 1.13 2011/05/09 16:33:05 onoratog Exp $
// $Author: onoratog $ 
// $Date: 2011/05/09 16:33:05 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) I put the enum_type in the class.  If I do not, and if I
//    repeat this model in another enum class, then there will
//    be ambiguities about lastEnum.
// 2) I do want both the name() operator and the << operator.
//    If only <<, then users will need to make temp ostringsteams
//    if they want to do their own formatting.
// 3) There are some notes on alternate implementations 
//    at the end of the .cc file.
// 4) Root stores enum types as 32 bit ints.

#include <iostream>
#include <vector>
#include <string>

namespace mu2e {

  class GenId {

  public:

    // Need to keep the enum and the _name member in sync.
    enum enum_type {
      unknown,       particleGun,       conversionGun,
      cosmicToy,     cosmicDYB,         cosmic,            dio1,
      dio2,          dio3,              pionCapture,
      muonCapture,   muonDecayInFlight, ejectedProtonGun, 
      piEplusNuGun,  primaryProtonGun,  fromG4BLFile,      ePlusfromStoppedPi,
      ejectedNeutronGun, 
      lastEnum
    };
  
    // Keep this in sync with the enum. Used in GenId.cc
#define GENID_NAMES                                                                \
      "unknown",      "particleGun",       "conversionGun",                        \
      "cosmicToy",    "cosmicDYB",         "cosmic",           "dio1",             \
      "dio2",         "dio3",              "pionCapture",                          \
      "muonCapture",  "muonDecayInFlight", "ejectedProtonGun",                     \
      "piEplusNuGun", "primaryProtonGun",  "fromG4BLFile"    , "ePlusfromStoppedPi", "ejectedNeutronGun"

  public:

    // The most important c'tor and accessor methods are first.
    explicit GenId( enum_type id):
      _id(id)
    {};
    
    enum_type Id() const { return _id;}
    
    const std::string name() const { 
      return std::string( _name[_id] );
    }
    
    // ROOT requires a default c'tor.
    GenId():
      _id(unknown){
    };
    
    // Need this to interface with HepMC::GenEvent which stores an int.
    explicit GenId( int id):
      _id(static_cast<enum_type>(id)){
      if ( !isValid() ){
        // throw or something
        exit(-1);
      }
    }
    
    virtual ~GenId(){}

    bool operator==(const GenId g) const{
      return ( _id == g._id );
    }

    bool operator==(const GenId::enum_type g) const{
      return ( _id == g );
    }

    // Accessor for the version.
    static int version() { return _version; }
 
    // Static version of the name method.
    static const std::string name( enum_type id ){ 
      return std::string( _name[id] );
    }
    
    // Check validity of an Id. Unknown is defined to be valid.
    static bool isValid( enum_type id){
      if ( id <  unknown  ) return false;
      if ( id >= lastEnum ) return false;
      return true;
    }

    static void printAll( std::ostream& ost);

    static void printAll(){
      printAll(std::cout);
    }

    // Member function version of static functions.
    bool isValid() const{
      return isValid(_id);
    }

    // List of names corresponding to the enum.
    const static char* _name[];

  private:

    // The one and only per-instance member datum.
    enum_type _id;

    // Can this make sense?  What happens if I read in two different
    // files that have different versions?  Should I use cvs version instead?
    // This is really an edm question not a question for the class itself.
    static const uint32_t _version = 1000;
    
  };
  
  // Shift left (printing) operator.
  inline std::ostream& operator<<(std::ostream& ost,
                                  const GenId& id ){
    ost << "( "
        << id.Id() << ": "
        << id.name()
        << " )";
    return ost;
  }

}

#endif

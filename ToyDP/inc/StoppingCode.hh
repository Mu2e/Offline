#ifndef ToyDP_StoppingCode_hh
#define ToyDP_StoppingCode_hh
//
// An enum-matched-to-names class for stopping codes from G4.
//
// $Id: StoppingCode.hh,v 1.1 2010/11/10 23:42:27 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/11/10 23:42:27 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) The enum_type is in the class.  There are other classes
//    that follow this model and lastEnum must be specific to each.
//
// 2) I do want both the name() operator and the << operator.
//    If only <<, then users will need to make temp ostringsteams
//    if they want to do their own formatting.
//
// 3) There are some notes on alternate implementations in
//    ToyDP/inc/GenId.hh
//
// 4) Root stores enum types as 32 bit ints.

#include <iostream>
#include <vector>
#include <string>

namespace mu2e {

  class StoppingCode {

  public:

    // Need to keep the enum and the _name member in sync.
    // Add new elements just before lastEnum; do not insert new elements
    // prior to this - it will break backwards compatibility.
    // 
    // No real values yet - just a placeholder.
    //
    enum enum_type {
      unknown,
      lastEnum
    };
  
    // Keep this in sync with the enum. Used in StoppingCode.cc
#define STOPPINGCODE_NAMES                                                     \
      "unknown"

  public:

    // The most important c'tor and accessor methods are first.
    explicit StoppingCode( enum_type id):
      _id(id)
    {};

    // Constructor from an int; should not be needed often.
    explicit StoppingCode( int id):
      _id(static_cast<enum_type>(id)){
      if ( !isValid() ){
        // throw or something
        exit(-1);
      }
    }

    // ROOT requires a default c'tor.
    StoppingCode():
      _id(unknown){
    };

    // Accept compiler supplied d'tor, copy c'tor and assignment operator.
    
    enum_type id() const { return _id;}
    
    const std::string name() const { 
      return std::string( _name[_id] );
    }
    
    bool operator==(const StoppingCode g) const{
      return ( _id == g._id );
    }

    bool operator==(const StoppingCode::enum_type g) const{
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

    static void printAll( std::ostream& ost = std::cout);

    // Member function version of static functions.
    bool isValid() const{
      return isValid(_id);
    }

  private:

    // The one and only per-instance member datum.
    enum_type _id;

    // List of names corresponding to the enum.
    const static char* _name[];

    // Can this make sense?  What happens if I read in two different
    // files that have different versions?  Should I use cvs version instead?
    // This is really an edm question not a question for the class itself.
    static const uint32_t _version = 1000;
    
  };
  
  // Shift left (printing) operator.
  inline std::ostream& operator<<(std::ostream& ost,
                                  const StoppingCode& id ){
    ost << "( "
        << id.id() << ": "
        << id.name()
        << " )";
    return ost;
  }

}

#endif

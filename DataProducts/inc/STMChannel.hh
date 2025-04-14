#ifndef DataProducts_STMChannel_hh
#define DataProducts_STMChannel_hh
//
// An enum-matched-to-names class for STM channels
// based on GenId.hh
//

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "canvas/Utilities/InputTag.h"

namespace mu2e {

  class STMChannel {

  public:

    // Need to keep the enum and the _name member in sync.
    // Note 1. Unknown has a non-standard numerica value since we start channel numbering at 0
    // Note 2. We restrict to 8-bits because this may be in every STMDigi
    enum enum_type : int8_t { HPGe, LaBr, unknown, lastEnum };

#ifndef SWIG
    // Keep this in sync with the enum. Used in STMChannel.cc
#define STMCHANNEL_NAMES \
    "HPGe", "LaBr", "unknown"
#endif

  public:

    // The most important c'tor and accessor methods are first.
    STMChannel( enum_type id):
      _id(id)
    {}

    enum_type id() const { return _id;}

    const std::string name() const {
      return std::string( _name[_id] );
    }

    // ROOT requires a default c'tor.
    STMChannel():
      _id(unknown){
    }

    virtual ~STMChannel() = default;

    bool operator==(const STMChannel g) const{
      return ( _id == g._id );
    }

    bool operator==(const STMChannel::enum_type g) const{
      return ( _id == g );
    }

    bool operator!=(const STMChannel::enum_type g) const{
      return ( _id != g );
    }

    bool operator!=(const STMChannel g) const{
      return ( _id != g._id );
    }

    bool operator<(const STMChannel g) const{
      return ( _id < g._id );
    }

    // Accessor for the version.
    static int version() { return _version; }

    // Static version of the name method.
    static const std::string name( enum_type id ){
      return std::string( _name[id] );
    }

    // Check validity of an Id. Unknown is defined to be valid.
    // Unknown has a non-standard numeric value
    static bool isValid( enum_type id){
      if ( id <  0  ) return false;
      if ( id >= unknown ) return false;
      return true;
    }

    // Return the STMChannel that corresponds to this name.
    static STMChannel findByName ( std::string const& name);

    static void printAll( std::ostream& ost);

    static void printAll(){
      printAll(std::cout);
    }

    // Member function version of static functions.
    bool isValid() const{
      return isValid(_id);
    }

#ifndef SWIG
    // List of names corresponding to the enum.
    const static char* _name[];
#endif

    // Number of valid codes, not including lastEnum, but including "unknown".
    static std::size_t size(){
      return lastEnum;
    }

  private:

    // The one and only per-instance member datum.
    enum_type _id;

    // Can this make sense?  What happens if I read in two different
    // files that have different versions?  Should I use cvs version instead?
    // This is really an edm question not a question for the class itself.
    static const unsigned _version = 1000;

  };

  // Shift left (printing) operator.
  inline std::ostream& operator<<(std::ostream& ost,
                                  const STMChannel& id ){
    ost << "( "
        << static_cast<int>(id.id()) << ": "
        << id.name()
        << " )";
    return ost;
  }

}

#endif /* DataProducts_STMChannel_hh */

#ifndef Compression_CompressionLevel_hh
#define Compression_CompressionLevel_hh
//
// An enum-matched-to-names class used to indicate the level of compression
// requested. The class contains enum entries for the different possible
// levels of compression.
//
// Original author Rob Kutschke
// Modified for Compression: Andy Edmonds, Dec 2020
//
// Notes:
//
// 1) The enum_type is in the class.  There are other classes
//    that follow this model and lastEnum must be specific to each.
//
// 2) I do want both the name() operator and the << operator.
//    If only <<, then users will need to make temp ostringsteams
//    if they want to do their own formatting.
//
// 3) There are some notes on alternate implementations in
//    Compression/inc/ProcessCode.hh and Compression/inc/GenId.hh
//
// 4) Root stores enum types as 32 bit ints.

#include <iostream>
#include <string>
#include <vector>

namespace mu2e {

  class CompressionLevel {

  public:

    // Need to keep the enum and the _name member in sync.
    // Add new elements just before lastEnum; do not insert new elements
    // prior to this - it will break backwards compatibility.
    enum enum_type { unknown, kNoCompression, kSimParticleCompression, kFullCompression,
      lastEnum
    };

    // Keep this list of names in sync with the enum. Used in CompressionLevel.cc
    // lastEnum does not appear in this list of names.
#define COMPRESSIONLEVEL_NAMES  "unknown", "noCompression", "simParticleCompression", "fullCompression"

  public:

    // The most important c'tor and accessor methods are first.
    CompressionLevel( enum_type id):
      _id(id)
    {}

    // ROOT requires a default c'tor.
    CompressionLevel():
      _id(unknown){
    }

    // Accept compiler supplied d'tor, copy c'tor and assignment operator.

    enum_type id() const { return _id;}

    // Return the name that corresponds to this enum.
    std::string name() const {
      return std::string( _name[_id] );
    }

    // Return the process code that corresponds to this name.
    // Return CompressionLevel(unknown) if there is no such string.
    static CompressionLevel findByName ( std::string const& name);

    // This operator implements:
    //   CompressionLevel a;
    //   enum_type b;
    // a = b;
    CompressionLevel& operator=(CompressionLevel::enum_type const& c){
      _id = c;
      return *this;
    }

    // This operator implements:
    //   CompressionLevel a;
    //   enum_type b = a;
    operator CompressionLevel::enum_type ()const{
      return _id;
    }

    // Tests for equality.
    bool operator==(const CompressionLevel g) const{
      return ( _id == g._id );
    }

    bool operator==(const CompressionLevel::enum_type g) const{
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

    // Number of valid codes, not including lastEnum, but including "unknown".
    static std::size_t size(){
      return lastEnum;
    }

    // Print all valid codes and their text strings.
    static void printAll( std::ostream& ost = std::cout);

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
    static const unsigned _version = 1000;

  };

  // Shift left (printing) operator.
  inline std::ostream& operator<<(std::ostream& ost,
                                  const CompressionLevel& id ){
    ost << "( "
        << id.id() << ": "
        << id.name()
        << " )";
    return ost;
  }

}

#endif /* Compression_CompressionLevel_hh */

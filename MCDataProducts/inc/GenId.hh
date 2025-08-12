#ifndef MCDataProducts_GenId_hh
#define MCDataProducts_GenId_hh
//
// An enum-matched-to-names class for generator Id's.
//
//
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

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace mu2e {

  class GenId {

  public:

    // Need to keep the enum and the _name member in sync.
    enum enum_type {
      unknown,       particleGun,       CeEndpoint,
      cosmicToy,     cosmicDYB,         cosmic,          obsolete1, //6
      dioTail,       obsolete2,         obsolete3,       obsolete4,           ExternalRPC, //11
      muonCapture,   muonDecayInFlight, ejectedProtonGun, //14
      piEplusNuGun,  primaryProtonGun,  fromG4BLFile,      ePlusfromStoppedPi, //18
      ejectedNeutronGun, ejectedPhotonGun, nuclearCaptureGun, InternalRPC, //22
      extMonFNALGun, fromStepPointMCs, stoppedMuonGun, PiCaptureCombined, //26
      MARS, StoppedParticleReactionGun, bremElectronGun, muonicXRayGun, //30
      fromSimParticleStartPoint, fromSimParticleCompact, StoppedParticleG4Gun, //33
      CaloCalib, InFlightParticleSampler, muplusDecayGun, StoppedMuonXRayGammaRayGun, //37
      cosmicCRY,  pbarFlat, fromAscii, ExternalRMC, InternalRMC, CeLeadingLog, cosmicCORSIKA, //44
      MuCapProtonGenTool, MuCapDeuteronGenTool, DIOGenTool, MuCapNeutronGenTool, // 48
      MuCapPhotonGenTool, MuCapGammaRayGenTool, CeLeadingLogGenTool, MuplusMichelGenTool,// 52
      gammaPairProduction, antiproton, Mu2eXGenTool,//55
      lastEnum //56
    };

#ifndef SWIG
    // Keep this in sync with the enum. Used in GenId.cc
#define GENID_NAMES                                                     \
    "unknown",      "particleGun",       "CeEndpoint",               \
      "cosmicToy",    "cosmicDYB",         "cosmic",           "obsolete1",  \
      "dioTail", "obsolete2",  "obsolete3", "obsolete4",           "ExternalRPC", \
      "muonCapture",  "muonDecayInFlight", "ejectedProtonGun",          \
      "piEplusNuGun", "primaryProtonGun",  "fromG4BLFile"    , "ePlusfromStoppedPi", \
      "ejectedNeutronGun", "ejectedPhotonGun", "nuclearCaptureGun", "InternalRPC", \
      "extMonFNALGun", "fromStepPointMCs", "stoppedMuonGun", "PiCaptureCombined", \
      "MARS", "StoppedParticleReactionGun","bremElectronGun", "muonicXRayGun", \
      "fromSimParticleStartPoint", "fromSimParticleCompact", "StoppedParticleG4Gun", \
      "CaloCalib", "InFlightParticleSampler","muplusDecayGun", "StoppedMuonXRayGammaRayGun", \
      "CosmicCRY", "pbarFlat","fromAscii","ExternalRMC","InternalRMC","CeLeadingLog", "CosmicCORSIKA", \
    "MuCapProtonGenTool", "MuCapDeuteronGenTool", "DIOGenTool", "MuCapNeutronGenTool", \
      "MuCapPhotonGenTool", "MuCapGammaRayGenTool","CeLeadingLogGenTool","MuplusMichelGenTool", \
      "gammaPairProduction", "antiproton", "Mu2eXGenTool"
#endif

  public:

    // The most important c'tor and accessor methods are first.
    GenId( enum_type id):
      _id(id)
    {}

    // c'tor from a string.
    // If the name is not recognzied or has the value "unknown", it throws an exception.
    GenId( std::string const& name){
      _id = findByName(name).id();
    }

    enum_type id() const { return _id;}

    const std::string name() const {
      return std::string( _name[_id] );
    }

    // ROOT requires a default c'tor.
    GenId():
      _id(unknown){
    }

    virtual ~GenId(){}

    bool operator==(const GenId g) const{
      return ( _id == g._id );
    }

    bool operator==(const GenId::enum_type g) const{
      return ( _id == g );
    }

    bool operator!=(const GenId::enum_type g) const{
      return ( _id != g );
    }

    bool operator!=(const GenId g) const{
      return ( _id != g._id );
    }

    bool isCosmic() const {
      return (_id == cosmicToy || _id == cosmicDYB || _id == cosmic || _id == cosmicCRY || _id == cosmicCORSIKA);
    }

    bool isConversion() const {
      return _id == GenId::CeEndpoint || _id == GenId::CeLeadingLog;
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

    // Return the GenId that corresponds to this name.
    static GenId findByName ( std::string const& name, bool throwIfUnknown = true, bool throwIfUndefined = true);

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
                                  const GenId& id ){
    ost << "( "
        << id.id() << ": "
        << id.name()
        << " )";
    return ost;
  }

}

#endif /* MCDataProducts_GenId_hh */

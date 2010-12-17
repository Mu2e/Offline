#ifndef ToyDP_StoppingCode_hh
#define ToyDP_StoppingCode_hh
//
// An enum-matched-to-names class for physics processes from G4,
// plus mu2e defined reasons for stopping particles.
//
// $Id: StoppingCode.hh,v 1.2 2010/12/17 22:21:43 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/12/17 22:21:43 $
//
// Original author Rob Kutschke
//
// Notes:
// 0) The names starting with mu2e represent actions taken by SteppingAction.
//    The other codes are the names of G4 processes.
//
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
    enum enum_type {
      unknown,                AlphaInelastic,          annihil,             AntiLambdaInelastic,
      AntiNeutronInelastic,   AntiOmegaMinusInelastic, AntiProtonInelastic, AntiSigmaMinusInelastic,
      AntiSigmaPlusInelastic, AntiXiMinusInelastic,    AntiXiZeroInelastic, CHIPSNuclearCaptureAtRest,
      compt,                  conv,                    Decay,               DeuteronInelastic,
      eBrem,                  eIoni,                   ElectroNuclear,      hBrems,
      hElastic,               hIoni,                   hPairProd,           ionIoni,
      KaonMinusInelastic,     KaonPlusInelastic,       KaonZeroLInelastic,  KaonZeroSInelastic,
      LambdaInelastic,        msc,                     muBrems,             muIoni,
      muMinusCaptureAtRest,   muMsc,                   muPairProd,          nCapture,
      NeutronInelastic,       nFission,                nKiller,             OmegaMinusInelastic,
      phot,                   PhotonInelastic,         PionMinusInelastic,  PionPlusInelastic,
      PositronNuclear,        ProtonInelastic,         SigmaMinusInelastic, SigmaPlusInelastic,
      StepLimiter,            Transportation,          TritonInelastic,     XiMinusInelastic,
      XiZeroInelastic,        mu2eLowEKine,            mu2eHallAir,         mu2eMaxSteps,
      lastEnum
    };
  
    // Keep this list of names in sync with the enum. Used in StoppingCode.cc
    // lastEnum does not appear in this list of names.
#define STOPPINGCODE_NAMES                                                                                   \
    "unknown",                "AlphaInelastic",          "annihil",             "AntiLambdaInelastic",       \
    "AntiNeutronInelastic",   "AntiOmegaMinusInelastic", "AntiProtonInelastic", "AntiSigmaMinusInelastic",   \
    "AntiSigmaPlusInelastic", "AntiXiMinusInelastic",    "AntiXiZeroInelastic", "CHIPSNuclearCaptureAtRest", \
    "compt",                  "conv",                    "Decay",               "DeuteronInelastic",         \
    "eBrem",                  "eIoni",                   "ElectroNuclear",      "hBrems",                    \
    "hElastic",               "hIoni",                   "hPairProd",           "ionIoni",                   \
    "KaonMinusInelastic",     "KaonPlusInelastic",       "KaonZeroLInelastic",  "KaonZeroSInelastic",        \
    "LambdaInelastic",        "msc",                     "muBrems",             "muIoni",                    \
    "muMinusCaptureAtRest",   "muMsc",                   "muPairProd",          "nCapture",                  \
    "NeutronInelastic",       "nFission",                "nKiller",             "OmegaMinusInelastic",       \
    "phot",                   "PhotonInelastic",         "PionMinusInelastic",  "PionPlusInelastic",         \
    "PositronNuclear",        "ProtonInelastic",         "SigmaMinusInelastic", "SigmaPlusInelastic",        \
    "StepLimiter",            "Transportation",          "TritonInelastic",     "XiMinusInelastic",          \
    "XiZeroInelastic",        "mu2eLowEKine",            "mu2eHallAir",         "mu2eMaxSteps"

  public:

    // The most important c'tor and accessor methods are first.
    StoppingCode( enum_type id):
      _id(id)
    {}

    // Constructor from an int; should not be needed often.  This checks for validity and throws.
    explicit StoppingCode( int id);

    // ROOT requires a default c'tor.
    StoppingCode():
      _id(unknown){
    }

    // Accept compiler supplied d'tor, copy c'tor and assignment operator.
    
    enum_type id() const { return _id;}

    // Return the name that corresponds to this enum.
    const std::string name() const { 
      return std::string( _name[_id] );
    }

    // Return the stopping code that corresponds to this name.
    // Return StoppingCode(unknown) if there is no such string.
    static StoppingCode findByName ( std::string const& name);

    // Stopping Code a;
    // enum_type b;
    // a = b;
    StoppingCode& operator=(StoppingCode::enum_type const& c){
      _id = c;
      return *this;
    }

    // Stopping Code a;
    // enum_type b = a;
    operator StoppingCode::enum_type ()const{
      return _id;
    }

    // Tests for equality.
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

    // Number of valid codes, not including lastEnum, but including "unknown".
    static size_t size(){
      return lastEnum;
    }

    // Print all valid codes and their text strings.
    static void printAll( std::ostream& ost = std::cout);

    // Return a list of codes that are mu2e specific.
    static std::vector<StoppingCode> mu2eCodes();

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

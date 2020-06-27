#ifndef MCDataProducts_ProcessCode_hh
#define MCDataProducts_ProcessCode_hh
//
// An enum-matched-to-names class used to indicate why a SimParticle was
// created and why it stopped. The class contains enum entries for all
// physics processes known in G4; it also contains an enum entry to indicate
// that the particle is a primary particle and other enum entries to
// indicate that a particle was killed in one of the user actions written by G4.
//
// $Id: ProcessCode.hh,v 1.19 2014/07/29 20:32:11 genser Exp $
// $Author: genser $
// $Date: 2014/07/29 20:32:11 $
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
//    MCDataProducts/inc/GenId.hh
//
// 4) Root stores enum types as 32 bit ints.

#include <iostream>
#include <string>
#include <vector>

namespace mu2e {

  class ProcessCode {

  public:

    // Need to keep the enum and the _name member in sync.
    // Add new elements just before lastEnum; do not insert new elements
    // prior to this - it will break backwards compatibility.
    enum enum_type {
      unknown,                AlphaInelastic,          annihil,             AntiLambdaInelastic, // 3
      AntiNeutronInelastic,   AntiOmegaMinusInelastic, AntiProtonInelastic, AntiSigmaMinusInelastic, // 7
      AntiSigmaPlusInelastic, AntiXiMinusInelastic,    AntiXiZeroInelastic, CHIPSNuclearCaptureAtRest, // 11
      compt,                  conv,                    Decay,               DeuteronInelastic, // 15
      eBrem,                  eIoni,                   ElectroNuclear,      hBrems, // 19
      hElastic,               hIoni,                   hPairProd,           ionIoni,  // 23
      KaonMinusInelastic,     KaonPlusInelastic,       KaonZeroLInelastic,  KaonZeroSInelastic, // 27
      LambdaInelastic,        msc,                     muBrems,             muIoni,  // 31
      muMinusCaptureAtRest,   muMsc,                   muPairProd,          nCapture,  // 35
      NeutronInelastic,       nFission,                nKiller,             OmegaMinusInelastic, // 39
      phot,                   PhotonInelastic,         PionMinusInelastic,  PionPlusInelastic, // 43
      PositronNuclear,        ProtonInelastic,         SigmaMinusInelastic, SigmaPlusInelastic, // 47
      StepLimiter,            Transportation,          TritonInelastic,     XiMinusInelastic, // 51
      XiZeroInelastic,        mu2eLowEKine,            mu2eKillerVolume,    mu2eMaxSteps,  // 55
      mu2ePrimary,            muMinusConversionAtRest, hadElastic,          CoulombScat, // 59
      nuclearStopping,        mu2eMaxGlobalTime,       TNuclearCapture,     muMinusAtomicCapture, // 63
      MuAtomDecay,            Rayl,                    ionInelastic,        He3Inelastic, // 67
      alphaInelastic,         AntiHe3InelasticProcess, AntiAlphaInelasticProcess, AntiDeuteronInelastic, // 71
      dInelastic,             tInelastic,              RadioactiveDecay,    CHIPS_Inelastic, // 75
      NotSpecified,           hFritiofCaptureAtRest,   hBertiniCaptureAtRest, AntiTritonInelasticProcess, // 79
      anti_He3Inelastic,      anti_alphaInelastic,     anti_deuteronInelastic, anti_lambdaInelastic, // 83
      anti_neutronInelastic,  anti_omega_MinusInelastic, anti_protonInelastic, anti_sigma_PlusInelastic,  // 87
      anti_sigma_MinusInelastic, anti_tritonInelastic, anti_xi_MinusInelastic, anti_xi0Inelastic,  // 91
      kaon_PlusInelastic,     kaon_MinusInelastic,     kaon0LInelastic,     kaon0SInelastic, // 95
      lambdaInelastic,        neutronInelastic,        omega_MinusInelastic, pi_PlusInelastic,  // 99
      pi_MinusInelastic,      protonInelastic,         sigma_PlusInelastic, sigma_MinusInelastic, // 103
      sigma0Inelastic,        xi_MinusInelastic,       xi0Inelastic,        positronNuclear, // 107
      electronNuclear,        photonNuclear,           antilambdaInelastic, DecayWithSpin, // 111
      ionElastic,             EMCascade,               DIO,                 NuclearCapture, // 115
      muonNuclear,            GammaToMuPair,           AnnihiToMuPair,      ee2hadr, // 119
      G4MinEkineCuts,         G4MaxTimeCuts,           OpAbsorption,        OpBoundary, // 123
      Scintillation,          inelastic,               G4ErrorEnergyLoss,   G4ErrorStepLengthLimit, // 127
      G4ErrorMagFieldLimit,   ePairProd,               FieldPropagator,     Mu2eRecorderProcess,  // 131
      mu2eProtonInelastic,    RadioactiveDecayBase,
      lastEnum,
      // An alias for backward compatibility
      mu2eHallAir = mu2eKillerVolume
    };

    // Keep this list of names in sync with the enum. Used in ProcessCode.cc
    // lastEnum does not appear in this list of names.
#define PROCESSCODE_NAMES                                                                                   \
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
    "XiZeroInelastic",        "mu2eLowEKine",            "mu2eKillerVolume",    "mu2eMaxSteps",              \
    "mu2ePrimary",            "muMinusConversionAtRest", "hadElastic",          "CoulombScat",               \
    "nuclearStopping",        "mu2eMaxGlobalTime",       "TNuclearCapture",     "muMinusAtomicCapture",      \
    "MuAtomDecay",            "Rayl",                    "ionInelastic",        "He3Inelastic",              \
    "alphaInelastic",         "AntiHe3InelasticProcess", "AntiAlphaInelasticProcess", "AntiDeuteronInelastic", \
    "dInelastic",             "tInelastic",              "RadioactiveDecay",    "CHIPS_Inelastic",           \
    "NotSpecified",           "hFritiofCaptureAtRest",   "hBertiniCaptureAtRest", "AntiTritonInelasticProcess", \
    "anti_He3Inelastic",      "anti_alphaInelastic",     "anti_deuteronInelastic", "anti_lambdaInelastic", \
    "anti_neutronInelastic",  "anti_omega-Inelastic",    "anti_protonInelastic",   "anti_sigma+Inelastic", \
    "anti_sigma-Inelastic",   "anti_tritonInelastic",    "anti_xi-Inelastic",      "anti_xi0Inelastic", \
    "kaon+Inelastic",         "kaon-Inelastic",          "kaon0LInelastic",        "kaon0SInelastic", \
    "lambdaInelastic",        "neutronInelastic",        "omega-Inelastic",        "pi+Inelastic", \
    "pi-Inelastic",           "protonInelastic",         "sigma+Inelastic",        "sigma-Inelastic", \
    "sigma0Inelastic",        "xi-Inelastic",            "xi0Inelastic",           "positronNuclear", \
    "electronNuclear",        "photonNuclear",           "anti-lambdaInelastic",   "DecayWithSpin", \
    "ionElastic",             "EMCascade",               "DIO",                    "NuclearCapture", \
    "muonNuclear",            "GammaToMuPair",           "AnnihiToMuPair",         "ee2hadr", \
    "G4MinEkineCuts",         "G4MaxTimeCuts",           "OpAbsorption",           "OpBoundary", \
    "Scintillation",          "inelastic",               "G4ErrorEnergyLoss",      "G4ErrorStepLengthLimit", \
    "G4ErrorMagFieldLimit",   "ePairProd",               "FieldPropagator",        "Mu2eRecorderProcess", \
    "mu2eProtonInelastic",    "RadioactiveDecayBase"

  public:

    // The most important c'tor and accessor methods are first.
    ProcessCode( enum_type id):
      _id(id)
    {}

    // ROOT requires a default c'tor.
    ProcessCode():
      _id(unknown){
    }

    // Accept compiler supplied d'tor, copy c'tor and assignment operator.

    enum_type id() const { return _id;}

    // Return the name that corresponds to this enum.
    std::string name() const {
      return std::string( _name[_id] );
    }

    // Return the process code that corresponds to this name.
    // Return ProcessCode(unknown) if there is no such string.
    static ProcessCode findByName ( std::string const& name);

    // This operator implements:
    //   ProcessCode a;
    //   enum_type b;
    // a = b;
    ProcessCode& operator=(ProcessCode::enum_type const& c){
      _id = c;
      return *this;
    }

    // This operator implements:
    //   ProcessCode a;
    //   enum_type b = a;
    operator ProcessCode::enum_type ()const{
      return _id;
    }

    // Tests for equality.
    bool operator==(const ProcessCode g) const{
      return ( _id == g._id );
    }

    bool operator==(const ProcessCode::enum_type g) const{
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

    // Return a list of codes that are mu2e specific.
    static std::vector<ProcessCode> mu2eCodes();

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
                                  const ProcessCode& id ){
    ost << "( "
        << id.id() << ": "
        << id.name()
        << " )";
    return ost;
  }

}

#endif /* MCDataProducts_ProcessCode_hh */

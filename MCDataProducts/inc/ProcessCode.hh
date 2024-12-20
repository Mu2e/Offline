#ifndef MCDataProducts_ProcessCode_hh
#define MCDataProducts_ProcessCode_hh
//
// An enum-matched-to-names class used to indicate why a SimParticle was
// created and why it stopped. The class contains enum entries for all
// physics processes known in G4; it also contains an enum entry to indicate
// that the particle is a primary particle and other enum entries to
// indicate that a particle was killed in one of the user actions written by G4.
//
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
      mu2ePrimary,            mu2eSpecialCutsProcess,  hadElastic,          CoulombScat, // 59
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
      G4ErrorMagFieldLimit,   ePairProd,               mu2eFieldPropagator, mu2eRecorderProcess,  // 131
      mu2eProtonInelastic,    RadioactiveDecayBase,    B_PlusInelastic,     B_MinusInelastic, // 135
      B0Inelastic,            Bc_PlusInelastic,        Bc_MinusInelastic,   Bs0Inelastic, // 139
      D_PlusInelastic,        D_MinusInelastic,        D0Inelastic,         Ds_PlusInelastic, // 143
      Ds_MinusInelastic,      anti_B0Inelastic,        anti_Bs0Inelastic,   anti_D0Inelastic, // 147
      anti_lambda_bInelastic, anti_lambda_c_PlusInelastic, anti_omega_b_MinusInelastic, anti_omega_c0Inelastic, // 151
      anti_xi_b_MinusInelastic, anti_xi_b0Inelastic,   anti_xi_c_PlusInelastic, anti_xi_c0Inelastic, // 155
      lambda_bInelastic,      lambda_c_PlusInelastic,  omega_b_MinusInelastic, omega_c0Inelastic, // 159
      xi_b_MinusInelastic,    xi_b0Inelastic,          xi_c_PlusInelastic,  xi_c0Inelastic, // 163
      // stopped-muon physics processes, specific to Mu2e
      truncated,       mu2eMuonCaptureAtRest,  mu2eMuonDecayAtRest,       mu2eCeMinusEndpoint, // 167
      mu2eCeMinusLeadingLog,   mu2eCePlusEndpoint,  mu2eDIOLeadingLog, mu2eInternalRMC,  // 171
      mu2eExternalRMC,         mu2eFlateMinus,      mu2eFlatePlus, mu2eFlatPhoton, // 175
      mu2eCePlusLeadingLog, mu2ePionCaptureAtRest, mu2eExternalRPC, mu2eInternalRPC, // 179
      mu2eCaloCalib, mu2ePienu, mu2eunused7, mu2eunused8, // 183
      uninitialized, NoProcess, GammaGeneralProc, // 186
      mu2eGammaConversion, Radioactivation, nCaptureHP, nFissionHP, mu2eAntiproton, // 191
      lastEnum,
      // An alias for backward compatibility
      mu2eHallAir = mu2eKillerVolume
    };

#ifndef SWIG
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
    "mu2ePrimary",            "mu2eSpecialCutsProcess",  "hadElastic",          "CoulombScat",               \
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
    "G4ErrorMagFieldLimit",   "ePairProd",               "mu2eFieldPropagator",    "mu2eRecorderProcess", \
    "mu2eProtonInelastic",    "RadioactiveDecayBase",      "B+Inelastic",         "B-Inelastic", \
    "B0Inelastic",            "Bc+Inelastic",            "Bc-Inelastic",        "Bs0Inelastic", \
    "D+Inelastic",            "D-Inelastic",             "D0Inelastic",         "Ds+Inelastic", \
    "Ds-Inelastic",           "anti_B0Inelastic",        "anti_Bs0Inelastic",   "anti_D0Inelastic", \
    "anti_lambda_bInelastic", "anti_lambda_c+Inelastic", "anti_omega_b-Inelastic", "anti_omega_c0Inelastic", \
    "anti_xi_b-Inelastic",    "anti_xi_b0Inelastic",     "anti_xi_c+Inelastic", "anti_xi_c0Inelastic", \
    "lambda_bInelastic",      "lambda_c+Inelastic",      "omega_b-Inelastic",   "omega_c0Inelastic", \
    "xi_b-Inelastic",         "xi_b0Inelastic",          "xi_c+Inelastic",      "xi_c0Inelastic", \
    "truncated", "mu2eMuonCaptureAtRest", "mu2eMuonDecayAtRest",  "mu2eCeMinusEndpoint", \
    "mu2eCeMinusLeadingLog", "mu2eCePlusEndpoint",  "mu2eDIOLeadingLog", "mu2eInternalRMC", \
    "mu2eExternalRMC",  "mu2eFlateMinus",      "mu2eFlatePlus", "mu2eFlatPhoton", \
  "mu2eCePlusLeadingLog", "mu2ePionCaptureAtRest", "mu2eExternalRPC", "mu2eInternalRPC", \
    "mu2eCaloCalib", "mu2ePienu", "mu2eunused7", "mu2eunused8", \
      "uninitialized", "NoProcess", "GammaGeneralProc", \
      "mu2eGammaConversion","Radioactivation", "nCaptureHP", "nFissionHP", "mu2eAntiproton"
#endif

  public:

    // The most important c'tor and accessor methods are first.
    ProcessCode( enum_type id):
      _id(id)
    {}

    // ROOT requires a default c'tor.
    ProcessCode():
      _id(uninitialized){
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

#ifndef SWIG
    // This operator implements:
    //   ProcessCode a;
    //   enum_type b;
    // a = b;
    ProcessCode& operator=(ProcessCode::enum_type const& c){
      _id = c;
      return *this;
    }
#endif

#ifndef SWIG
    // This operator implements:
    //   ProcessCode a;
    //   enum_type b = a;
    operator ProcessCode::enum_type ()const{
      return _id;
    }
#endif

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
      if ( id <  unknown  || id >= lastEnum ) return false;
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

#ifndef SWIG
    // List of names corresponding to the enum.
    const static char* _name[];
#endif

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

#ifndef DataProducts_CompressedPDGCode_hh
#define DataProducts_CompressedPDGCode_hh
//
// Enum-matched-to-String class defining an enum of compressed PDG codes.
// This is typically used to fill histograms with all PDG ids,
// but compressed to a few categories
// If the list is modified, update the min and max values
// used to create histograms
//

#include "Offline/GeneralUtilities/inc/EnumToStringSparse.hh"
#include <map>
#include <string>

namespace mu2e {

  class CompressedPDGCodeDetail
  {
  public:
    enum enum_type
      {
        anti_n0     = -7,
        anti_proton = -6,
        K_minus     = -5,
        pi_minus    = -4,
        anti_lepton = -3,
        mu_minus    = -2,
        e_minus     = -1,
        gamma       =  0,
        e_plus      =  1,
        mu_plus     =  2,
        lepton      =  3,
        pi_plus     =  4,
        K_plus      =  5,
        proton      =  6,
        n0          =  7,
        pi0         =  8,
        K_S0        =  9,
        K_L0        = 10,
        ud_meson    = 11,
        s_meson     = 12,
        cb_meson    = 13,
        ud_baryon   = 14,
        s_baryon    = 15,
        cb_baryon   = 16,
        other       = 17,
        nuclei      = 18,
        unknown     = 19
      };  // end enum CompressedPDGCode_type

    // needed to make histograms
    static constexpr int minBin = -7;
    static constexpr int maxBin = 19;

    static std::string const& typeName();
    static std::map<enum_type,std::string> const& names();

  }; // end class CompressedPDGCodeDetail

  typedef EnumToStringSparse<CompressedPDGCodeDetail> CompressedPDGCode;

} // end namespace mu2e.

#endif /* DataProducts_CompressedPDGCode_hh */

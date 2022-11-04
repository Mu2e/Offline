//
// Enum-matched-to-String class defining an enum of compressed PDG codes.
// This is typically used to fill histograms with all PDG ids,
// but compressed to a few categories
//

#include <type_traits>
#include <utility>

#include "Offline/DataProducts/inc/CompressedPDGCode.hh"

namespace mu2e {

  std::string const& CompressedPDGCodeDetail::typeName() {
    static const std::string type("CompressedPDGCode");
    return type;
  }

  static const std::map<CompressedPDGCodeDetail::enum_type,std::string> nam{
    {CompressedPDGCodeDetail::anti_n0       , "anti_n0" },      // -7
    {CompressedPDGCodeDetail::anti_proton   , "anti_proton" },  // -6
    {CompressedPDGCodeDetail::K_minus       , "K_minus" },      // -5
    {CompressedPDGCodeDetail::pi_minus      , "pi_minus" },     // -4
    {CompressedPDGCodeDetail::anti_lepton   , "anti_lepton" },  // -3
    {CompressedPDGCodeDetail::mu_minus      , "mu_minus" },     // -2
    {CompressedPDGCodeDetail::e_minus       , "e_minus" },      // -1
    {CompressedPDGCodeDetail::gamma         , "gamma" },        //  0
    {CompressedPDGCodeDetail::e_plus        , "e_plus" },       //  1
    {CompressedPDGCodeDetail::mu_plus       , "mu_plus" },      //  2
    {CompressedPDGCodeDetail::lepton        , "lepton" },       //  3
    {CompressedPDGCodeDetail::pi_plus       , "pi_plus" },      //  4
    {CompressedPDGCodeDetail::K_plus        , "K_plus" },       //  5
    {CompressedPDGCodeDetail::proton        , "proton" },       //  6
    {CompressedPDGCodeDetail::n0            , "n0" },           //  7
    {CompressedPDGCodeDetail::pi0           , "pi0" },          //  8
    {CompressedPDGCodeDetail::K_S0          , "K_S0" },         //  9
    {CompressedPDGCodeDetail::K_L0          , "K_L0" },         // 10
    {CompressedPDGCodeDetail::ud_meson      , "ud_meson" },     // 11
    {CompressedPDGCodeDetail::s_meson       , "s_meson" },      // 12
    {CompressedPDGCodeDetail::cb_meson      , "cb_meson" },     // 13
    {CompressedPDGCodeDetail::ud_baryon     , "ud_baryon" },    // 14
    {CompressedPDGCodeDetail::s_baryon      , "s_baryon" },     // 15
    {CompressedPDGCodeDetail::cb_baryon     , "cb_baryon" },    // 16
    {CompressedPDGCodeDetail::other         , "other" },        // 17
    {CompressedPDGCodeDetail::nuclei        , "nuclei" },       // 18
    {CompressedPDGCodeDetail::unknown       , "unknown" }       // 19
     }; // end initialization of static variable nam.
  std::map<CompressedPDGCodeDetail::enum_type,std::string> const& CompressedPDGCodeDetail::names(){
    return nam;
  }

} // end of namespace mu2e

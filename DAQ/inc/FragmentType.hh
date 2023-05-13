#ifndef mu2e_artdaq_core_Overlays_FragmentType_hh
#define mu2e_artdaq_core_Overlays_FragmentType_hh
#include "artdaq-core/Data/Fragment.hh"

namespace mu2e {

namespace detail {
enum FragmentType : artdaq::Fragment::type_t
{
	EMPTY = artdaq::Fragment::EmptyFragmentType,
	MISSED = artdaq::Fragment::FirstUserFragmentType,
	//DTC = artdaq::Fragment::FirstUserFragmentType + 1,  // DEPRECATED
	//MU2E = artdaq::Fragment::FirstUserFragmentType + 2, // DEPRECATED
	//MU2EEVENT = artdaq::Fragment::FirstUserFragmentType + 3, // DEPRECATED
	TRK = artdaq::Fragment::FirstUserFragmentType + 4,     // Tracker fragment
	CAL = artdaq::Fragment::FirstUserFragmentType + 5,     // Calorimeter fragment
	CRV = artdaq::Fragment::FirstUserFragmentType + 6,     // Cosmic Ray Veto fragment
	DBG = artdaq::Fragment::FirstUserFragmentType + 7,     // Debug Packet Fragment
	DTCEVT = artdaq::Fragment::FirstUserFragmentType + 8, // DTC Event Fragment
	INVALID  // Should always be last.
};

// Safety check.
static_assert(artdaq::Fragment::isUserFragmentType(FragmentType::INVALID - 1), "Too many user-defined fragments!");
}  // namespace detail

using detail::FragmentType;

std::unordered_map<FragmentType, std::string> const names{
	{FragmentType::MISSED, "MISSED"},
	//{FragmentType::DTC, "DTC"},   // DEPRECATED
	//{FragmentType::MU2E, "MU2E"}, // DEPRECATED
	//{FragmentType::MU2EEVENT, "MU2EEVENT"},
	{FragmentType::TRK, "TRK"},
	{FragmentType::CAL, "CAL"},
	{FragmentType::CRV, "CRV"},
	{FragmentType::DBG, "DBG"},
	{FragmentType::DTCEVT, "DTCEVT"}
};

FragmentType toFragmentType(std::string t_string);
std::string fragmentTypeToString(FragmentType val);

/**
 * \brief Create a list of all Fragment types defined by this package, in the format that RawInput expects
 * \return A list of all Fragment types defined by this package, in the format that RawInput expects
 */
std::map<artdaq::Fragment::type_t, std::string> makeFragmentTypeMap();
}  // namespace mu2e
#endif /* mu2e_artdaq_Overlays_FragmentType_hh */

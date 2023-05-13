#ifndef artdaq_core_Data_Mu2eEventFragment_hh
#define artdaq_core_Data_Mu2eEventFragment_hh

#include "Offline/DAQ/inc/DTC_Packets.h"
#include <memory>
#include "artdaq-core/Data/Fragment.hh"
#include "cetlib_except/exception.h"

// #include <ostream>
// #include <vector>

// Implementation of "DTCEventFragment", an artdaq::Fragment overlay class

namespace mu2e {
class DTCEventFragment;
}

/**
 * \brief The artdaq::DTCEventFragment class represents a Fragment which contains one or more DTC_Events
 */
class mu2e::DTCEventFragment
{
public:
	/// The current version of the DTCEventFragment
	static constexpr uint8_t CURRENT_VERSION = 1;

	/**
	 * \param f The Fragment object to use for data storage
	 *
	 * The constructor simply sets its const private member "artdaq_Fragment_"
	 * to refer to the artdaq::Fragment object
	 */
	explicit DTCEventFragment(artdaq::Fragment const& f)
		: artdaq_Fragment_(f) {}

	virtual ~DTCEventFragment()
	{
	}

	DTCLib::DTC_Event getData()
	{
		if (event_ptr_ == nullptr)
		{
			event_ptr_.reset(new DTCLib::DTC_Event(artdaq_Fragment_.dataBeginBytes()));
			event_ptr_->SetupEvent();
		}
		return *event_ptr_.get();
	}

	std::vector<DTCLib::DTC_DataBlock> getSubsystemData(DTCLib::DTC_Subsystem subsys)
	{
		auto data = getData();
		return data.GetSubsystemDataBlocks(subsys);
	}

protected:
private:
	DTCEventFragment(DTCEventFragment const&) = delete;             // DTCEventFragment should definitely not be copied
	DTCEventFragment(DTCEventFragment&&) = delete;                  // DTCEventFragment should not be moved, only the underlying Fragment
	DTCEventFragment& operator=(DTCEventFragment const&) = delete;  // DTCEventFragment should definitely not be copied
	DTCEventFragment& operator=(DTCEventFragment&&) = delete;       // DTCEventFragment should not be moved, only the underlying Fragment

	artdaq::Fragment const& artdaq_Fragment_;
	std::unique_ptr<DTCLib::DTC_Event> event_ptr_{nullptr};
};

#endif /* artdaq_core_Data_Mu2eEventFragment_hh */

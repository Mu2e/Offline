#include "Offline/DAQ/inc/CRVFragment.hh"

std::unique_ptr<mu2e::CRVFragment::CRVROCStatusPacket> mu2e::CRVFragment::GetCRVROCStatusPacket(size_t blockIndex) const
{
	auto dataPtr = dataAtBlockIndex(blockIndex);
	if (dataPtr == nullptr) return nullptr;

	std::unique_ptr<CRVROCStatusPacket> output(nullptr);
	output.reset(new CRVROCStatusPacket(*reinterpret_cast<CRVROCStatusPacket const*>(dataPtr->GetData())));
	return output;
}

std::vector<mu2e::CRVFragment::CRVHitReadoutPacket> mu2e::CRVFragment::GetCRVHitReadoutPackets(size_t blockIndex) const
{
	auto dataPtr = dataAtBlockIndex(blockIndex);
	if (dataPtr == nullptr) return std::vector<CRVHitReadoutPacket>();

	auto crvRocHdr = reinterpret_cast<CRVROCStatusPacket const*>(dataPtr->GetData());
        size_t nHits = 0;
        if(2*crvRocHdr->ControllerEventWordCount>sizeof(CRVROCStatusPacket))
	       nHits = (2*crvRocHdr->ControllerEventWordCount-sizeof(CRVROCStatusPacket)) / sizeof(CRVHitReadoutPacket);

	std::vector<CRVHitReadoutPacket> output(nHits);

	memcpy(&output[0], reinterpret_cast<CRVHitReadoutPacket const*>(crvRocHdr + 1), nHits * sizeof(CRVHitReadoutPacket));

	return output;
}

#include "Offline/DAQ/inc/CRVFragmentTmp.hh"

std::unique_ptr<mu2e::CRVFragmentTmp::CRVROCStatusPacket> mu2e::CRVFragmentTmp::GetCRVROCStatusPacket(size_t blockIndex) const
{
	auto dataPtr = dataAtBlockIndex(blockIndex);
	if (dataPtr == nullptr) return nullptr;

	std::unique_ptr<CRVROCStatusPacket> output(nullptr);
	output.reset(new CRVROCStatusPacket(*reinterpret_cast<CRVROCStatusPacket const*>(dataPtr->GetData())));
	return output;
}

std::vector<mu2e::CRVFragmentTmp::CRVHitReadoutPacket> mu2e::CRVFragmentTmp::GetCRVHitReadoutPackets(size_t blockIndex) const
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


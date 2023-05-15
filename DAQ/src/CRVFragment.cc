#include "Offline/DAQ/inc/CRVFragment.hh"

std::unique_ptr<mu2e::CRVFragment::CRVROCStatusPacket> mu2e::CRVFragment::GetCRVROCStatusPacket(size_t blockIndex) const
{
	auto dataPtr = dataAtBlockIndex(blockIndex);
	if (dataPtr == nullptr) return nullptr;

	std::unique_ptr<CRVROCStatusPacket> output(nullptr);
	output.reset(new CRVROCStatusPacket(*reinterpret_cast<CRVROCStatusPacket const*>(dataPtr->GetData())));
	return output;
}

std::vector<mu2e::CRVFragment::CRVHit> mu2e::CRVFragment::GetCRVHits(size_t blockIndex) const
{
	auto dataPtr = dataAtBlockIndex(blockIndex);
	if (dataPtr == nullptr) return std::vector<CRVHit>();

	auto crvRocHdr = reinterpret_cast<CRVROCStatusPacket const*>(dataPtr->GetData());
        size_t eventSize = 2*crvRocHdr->ControllerEventWordCount;
        size_t pos = sizeof(CRVROCStatusPacket);

        std::vector<mu2e::CRVFragment::CRVHit> output;
        while(pos<eventSize)
        {
          output.resize(output.size()+1);

	  memcpy(&output.back().first, reinterpret_cast<const uint8_t*>(dataPtr->GetData())+pos, sizeof(CRVHitInfo));
          pos += sizeof(CRVHitInfo);

          size_t nWaveformSamples = output.back().first.NumSamples;
          output.back().second.resize(nWaveformSamples);
	  memcpy(&output.back().second[0], reinterpret_cast<const uint8_t*>(dataPtr->GetData())+pos, nWaveformSamples*sizeof(CRVHitWaveformSample));
          pos += sizeof(CRVHitWaveformSample)*nWaveformSamples;
        }

	return output;
}

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
        size_t totalEventSize = 2*crvRocHdr->ControllerEventWordCount;
        size_t usedEventSize = sizeof(CRVROCStatusPacket);

        std::vector<mu2e::CRVFragment::CRVHit> output;
        auto hitInfoPtr = reinterpret_cast<CRVHitInfo const*>(crvRocHdr+1);
        while(usedEventSize<totalEventSize)
        {
          output.resize(output.size()+1);

	  memcpy(&output.back().first, hitInfoPtr, sizeof(CRVHitInfo));
          size_t nWaveformSamples = output.back().first.NumSamples;
          auto waveformPtr = reinterpret_cast<CRVHitWaveformSample const*>(hitInfoPtr+1);

          output.back().second.resize(nWaveformSamples);
	  memcpy(&output.back().second[0], waveformPtr, nWaveformSamples*sizeof(CRVHitWaveformSample));
          hitInfoPtr = reinterpret_cast<CRVHitInfo const*>(waveformPtr+nWaveformSamples);

          usedEventSize+=sizeof(CRVHitInfo) + nWaveformSamples*sizeof(CRVHitWaveformSample);
        }

	return output;
}

#include "CaloMC/inc/CaloWFExtractor.hh"
#include <algorithm>

namespace mu2e {

   void CaloWFExtractor::extract(const std::vector<int>& wf, std::vector<unsigned>& starts, std::vector<unsigned>& stops) const
   {
        unsigned timeSample(nBinsPeak_);
        while (timeSample+nBinsPeak_ < wf.size()) 
        {
	    // find the local maximum over threshold
	    if (wf[timeSample] < minPeakADC_) {++timeSample; continue;}

            auto it1 =  wf.begin()+timeSample-nBinsPeak_, it2=wf.begin()+timeSample+nBinsPeak_+1;
            if (std::max_element(it1,it2) != wf.begin()+timeSample) {++timeSample; continue;}

	    // find the starting / stopping point of the peak (stop = first value under threshold)
            auto     stopIter    = std::find_if(wf.begin()+timeSample, wf.end(),[this](const auto& a){return a < this->minPeakADC_;}); 
            unsigned sampleStop  = std::distance(wf.begin(),stopIter);
            unsigned sampleStart = (timeSample > bufferDigi_) ? timeSample - bufferDigi_ : 0;

	    starts.push_back(sampleStart);
	    stops.push_back(sampleStop);

	    //fast forward to end of waveform to search for next one 
	    timeSample = sampleStop+1;   
        }

        // Concatenate peaks
        unsigned iprev(0), icurrent(1);
        while (icurrent < starts.size())
        {
            if (stops[iprev] > starts[icurrent]) {stops[iprev]=stops[icurrent]; starts[icurrent]=stops[icurrent]=wf.size()+10;}
            else                                       {iprev = icurrent;}
            ++icurrent;
        }

        auto pred = [&wf](const unsigned a) {return a>wf.size();};
        starts.erase(std::remove_if(starts.begin(),starts.end(),pred),starts.end());
        stops.erase(std::remove_if(stops.begin(), stops.end(), pred),stops.end());
   }
}




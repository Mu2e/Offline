#include "Offline/CaloMC/inc/CaloWFExtractor.hh"
#include <algorithm>
#include <iostream>

namespace mu2e {

   void CaloWFExtractor::extract(const std::vector<int>& wf, std::vector<size_t>& starts, std::vector<size_t>& stops) const
   {
        size_t timeSample(nBinsPeak_+startOffset_);
        while (timeSample+nBinsPeak_ < wf.size())
        {
            // find starting point
            if (wf[timeSample] < minPeakADC_) {++timeSample; continue;}

            size_t imax(timeSample-nBinsPeak_);
            for (auto i = timeSample-nBinsPeak_; i<=timeSample+nBinsPeak_;++i) {if (wf[i]>wf[imax]) imax=i;}
            if (timeSample != imax)  {++timeSample; continue;}

            // find the starting / stopping point of the peak (stop = first value under threshold)
            size_t sampleStart = (timeSample > bufferDigi_) ? timeSample - bufferDigi_ : 0;
            size_t sampleStop(timeSample);
            ++sampleStop;
            while (sampleStop < wf.size() && wf[sampleStop] >= minPeakADC_) ++sampleStop;

            starts.push_back(sampleStart);
            stops.push_back(sampleStop);

            //fast forward to end of waveform to search for next one
            timeSample = sampleStop+1;
        }


        // Concatenate peaks and remove unused values (flag value to remove past wf.size() since the latter is a legitimate value)
        size_t iprev(0), icurrent(1);
        while (icurrent < starts.size())
        {
            if (stops[iprev] >= starts[icurrent]) {stops[iprev]=stops[icurrent]; starts[icurrent]=stops[icurrent]=wf.size()+1;}
            else                                  {iprev = icurrent;}
            ++icurrent;
        }

        auto pred = [&wf](const auto a) {return a>wf.size();};
        starts.erase(std::remove_if(starts.begin(),starts.end(),pred),starts.end());
        stops.erase(std::remove_if(stops.begin(), stops.end(), pred),stops.end());
   }
}


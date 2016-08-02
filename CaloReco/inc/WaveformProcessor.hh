#ifndef WaveformProcessor_HH
#define WaveformProcessor_HH

#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

  class WaveformProcessor {

    public:
      
      WaveformProcessor(fhicl::ParameterSet const& param) {};      
      virtual ~WaveformProcessor() {};

      virtual void   initialize() = 0;
      virtual void   reset() = 0;
      virtual void   extract(std::vector<double> &xInput, std::vector<double> &yInput) = 0;
      virtual void   plot(std::string pname) = 0;

      virtual int    nPeaks()                     const = 0;
      virtual double chi2()                       const = 0;
      virtual int    ndf()                        const = 0;
      virtual double amplitude(unsigned int i)    const = 0;
      virtual double amplitudeErr(unsigned int i) const = 0;
      virtual double time(unsigned int i)         const = 0;
      virtual double timeErr(unsigned int i)      const = 0;
      virtual bool   isPileUp(unsigned int i)     const = 0;
  
  };

}

#endif

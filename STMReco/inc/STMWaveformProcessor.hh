#ifndef STMWaveformProcessor_HH
#define STMWaveformProcessor_HH

#include <vector>
#include <string>

namespace mu2e {

  class STMWaveformProcessor {

     public:
        virtual ~STMWaveformProcessor() {};

        virtual void     initialize() = 0;
        virtual void     reset() = 0;
        virtual void     extract(const std::vector<double>& xInput, const std::vector<double>& yInput) = 0;
        virtual void     plot   (const std::string& pname) const = 0;

        virtual int      nPeaks()                     const = 0;
        virtual double   chi2()                       const = 0;
        virtual int      ndf()                        const = 0;
        virtual double   amplitude(unsigned int i)    const = 0;
        virtual double   amplitudeErr(unsigned int i) const = 0;
        virtual double   time(unsigned int i)         const = 0;
        virtual double   timeErr(unsigned int i)      const = 0;
        virtual bool     isPileUp(unsigned int i)     const = 0;
   };

}

#endif

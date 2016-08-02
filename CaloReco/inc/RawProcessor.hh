#ifndef RawProcessor_HH
#define RawProcessor_HH


#include "CaloReco/inc/WaveformProcessor.hh"
#include "fhiclcpp/ParameterSet.h"
#include <vector>


namespace mu2e {


  class RawProcessor : public WaveformProcessor {
     

     public:
     
                    RawProcessor(fhicl::ParameterSet const& param);
        virtual    ~RawProcessor() {};
	
        virtual void   initialize();
	virtual void   reset();
        virtual void   extract(std::vector<double> &xInput, std::vector<double> &yInput);
        virtual void   plot(std::string pname);

        virtual int    nPeaks()                     const {return nPeaks_;}
        virtual double chi2()                       const {return chi2_;}
        virtual int    ndf()                        const {return ndf_;}
        virtual double amplitude(unsigned int i)    const {return resAmp_.at(i);}
        virtual double amplitudeErr(unsigned int i) const {return resAmpErr_.at(i);}
        virtual double time(unsigned int i)         const {return resTime_.at(i);}
        virtual double timeErr(unsigned int i)      const {return resTimeErr_.at(i);}  
        virtual bool   isPileUp(unsigned int i)     const {return nPeaks_ > 1;}

	


    private:
       
       int                 windowPeak_;
       double              minPeakAmplitude_;
       double              shiftTime_;
       double              scaleFactor_;
       int                 diagLevel_;
       
       int                 nPeaks_;
       double              chi2_;
       int                 ndf_;
       std::vector<double> res_;
       std::vector<double> resAmp_;
       std::vector<double> resAmpErr_;
       std::vector<double> resTime_;
       std::vector<double> resTimeErr_;
       std::vector<double> xvec_;
       std::vector<double> yvec_;

  };

}
#endif

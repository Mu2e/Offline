#ifndef TrkHitReco_PeakFitFunction_hh
#define TrkHitReco_PeakFitFunction_hh

// Class holding functions which describe the ideal waveforms.  These are used to construct fits.

#include "TrkHitReco/inc/PeakFitParams.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
#include <functional>
#include "Rtypes.h"
#include "Math/ParamFunctor.h"

class TF1;
// Contains functions to be used for fitting by the PeakFit class

namespace mu2e {

  namespace TrkHitReco {

    // define the specific fit configuration among all options
    struct FitConfig 
    {
        enum fitOption : size_t {earlyPeak=0,latePeak,floatPedestal,floatWidth,truncateADC,findPeaks,nOptions};

        unsigned _options;
        int      _debug;  
        unsigned _maxnit; 

        bool hasOption(fitOption option) const { return (_options & (1<<option)) != 0; }
        void setOption(fitOption option) { _options |= (1<<option); }


        FitConfig() : _options(0), _debug(0), _maxnit(1) {}
        FitConfig(unsigned maxnit, int debug) : _options(0), _debug(debug), _maxnit(maxnit) {}
        FitConfig(const std::vector<fitOption>& options) : FitConfig() 
        {
	   for (auto iopt : options) setOption(iopt);	
        }
    };



    class PeakFitFunction : public ROOT::Math::ParamFunctor {

      public:

        // construct the model from the Straw conditions, the fit configuration
        PeakFitFunction(const StrawResponse& srep); 
        ~PeakFitFunction();

        //initialize the object once config is done
        void init(const FitConfig& fitConfig);
        // Returns the truncated response in ADC units
        void truncateResponse(Float_t& currentFunctionValue) const;
        
        // the actual fit function, taking the explicit parameters as input
        Double_t fitModel(Double_t time, PeakFitParams const& params) const;
        // the root version of same.  This calls down to the above
        Double_t fitModelRoot(Double_t* x, Double_t* p) const;
        // provide a TF1 using the fit function above
        TF1* fitModelTF1() const { return _tf1; }
        // Method for creating a TF1 using the fit function.  This OVERWRITES the state of the TF1
        void resetTF1(TF1& rootPeakFitFunction);
        // overwrite the TF1 based on a PeakFitParams object
        void resetTF1(PeakFitParamsLimits const& params);
        // override root base class function
        Double_t operator()(Double_t* x, Double_t* p);
        // Model for early peak approximated as exponential decay
        Float_t earlyPeak(const Double_t t,const Double_t charge) const;
        // Model for a unconvolved single peak
        // Shaping power set to 1
        Float_t unConvolvedSinglePeak(const Double_t t) const;
        // Convolution of single peak with uniform distribution
        Float_t convolvedSinglePeak(const Double_t t, const Double_t sigma) const;
     
      private:
        typedef std::function<Double_t(Double_t,PeakFitParams const&)> FitFunction;
        typedef std::function<Double_t(Double_t *,Double_t *)> FitFunctionRoot;

        const StrawResponse& _srep;      
        FitConfig _fitConfig;  
        TF1* _tf1;  

        void createTF1();

    }; 
  }
} 
#endif

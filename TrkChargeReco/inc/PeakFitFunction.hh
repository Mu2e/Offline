#ifndef TrkChargeReco_PeakFitFunction_hh
#define TrkChargeReco_PeakFitFunction_hh

// Class holding functions which describe the ideal waveforms.  These are used to construct fits.

#include "TrkChargeReco/inc/PeakFitParams.hh"
#include "TrackerConditions/inc/StrawElectronics.hh"
#include <functional>
#include "Rtypes.h"
#include "Math/ParamFunctor.h"

class TF1;
// Contains functions to be used for fitting by the PeakFit class

namespace mu2e {

  namespace TrkChargeReco {

// define the specific fit configuration among all options
    struct FitConfig {
      enum fitOption : size_t {earlyPeak=0,latePeak,floatPedestal,floatWidth,truncateADC,findPeaks,nOptions};
      // pack option fields into individual bits
      unsigned _options;
      int _debug; // debug level
      unsigned _maxnit; // maximum # of iterations
      // decode and encode functions
      bool hasOption(fitOption option) const { return (_options & (option<<1)) != 0; }
      void setOption(fitOption option) { _options |= (option<<1); }
      // default constructor zeros the options
      FitConfig() : _options(0), _debug(0), _maxnit(1) {}
      // construct from a vector of options
      FitConfig(std::vector<fitOption> const& options) : FitConfig() {
	for(auto iopt : options){
	  setOption(iopt);
	}
      }
    };

    class PeakFitFunction : public ROOT::Math::ParamFunctor {

      public:

// construct the model from the Straw conditions, the fit configuration
      PeakFitFunction(StrawElectronics const& strawele, FitConfig const& fitConfig); //
      ~PeakFitFunction();
// the actual fit function, taking the explicit parameters as input
      Double_t fitModel(Double_t time, PeakFitParams const& params) const;
// the root version of same.  This calls down to the above
      Double_t fitModelRoot(Double_t* x, Double_t* p) const;
// provide a TF1 using the fit function above
      TF1* fitModelTF1() const { return _tf1; }
// Method for creating a TF1 using the fit function.  This OVERWRITES the state of the TF1
      void resetTF1(TF1& rootPeakFitFunction) const;
// override root base class function
      Double_t operator()(Double_t* x, Double_t* p);
      private:
// provide the explicit function for this model
      typedef std::function<Double_t(Double_t,PeakFitParams const&)> FitFunction;
// provide the (c-style) function needed by root fits using TF1
      typedef std::function<Double_t(Double_t *,Double_t *)> FitFunctionRoot;
// general electronics parameters
      StrawElectronics const& _strawele;      
// configuration difinig the fit function
      FitConfig _fitConfig; // local cache of configuration
      TF1* _tf1; // root fitting object
// Helper functions, used in the above
      void createTF1();
      // Model for early peak approximated as exponential decay
      Float_t earlyPeak(const Double_t t,const Double_t charge) const;
      // Model for a unconvolved single peak
      // Shaping power set to 1
      Float_t unConvolvedSinglePeak(const Double_t t) const;
      // Convolution of single peak with uniform distribution
      Float_t convolvedSinglePeak(const Double_t t, const Double_t sigma) const;
      // Returns the truncated response in ADC units
      void truncateResponse(Float_t& currentFunctionValue) const;


    }; // PeakFitFunction class
  } //TrkChargeReco namespace
} // mu2e namespace
#endif

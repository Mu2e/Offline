#ifndef CaloNoiseARFitter_HH
#define CaloNoiseARFitter_HH

// Utility to perform auto-regressive fit of waveform

#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "SeedService/inc/SeedService.hh"
#include "CLHEP/Random/RandGaussQ.h"
#include <vector>


namespace mu2e {

  class CaloNoiseARFitter 
  {
     public:
        CaloNoiseARFitter(CLHEP::HepRandomEngine& engine, unsigned nParFit, int diagLevel);

        void                  setWaveform(const std::vector<double>& wf);
        void                  fitARCoeff();
        void                  generateWF(std::vector<double>& wf);

     private:
        void                  buildFullNoise(unsigned iRO);
        double                ARsigma(double sigma0);
        double                calcSigma(const std::vector<double>& values);
        void                  generateWF(std::vector<double>& wf, double sigma0);

        unsigned              nparFit_;
        std::vector<double>   param_;
        double                sigmaAR_;
        int                   status_;
        CLHEP::RandGaussQ     randGauss_;
        int                   diagLevel_;
   };
}
#endif

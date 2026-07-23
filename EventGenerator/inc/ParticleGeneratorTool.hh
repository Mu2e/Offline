#ifndef EventGenerator_ParticleGeneratorTool_hh
#define EventGenerator_ParticleGeneratorTool_hh

#include <memory>
#include <vector>
#include <string>

#include "CLHEP/Vector/LorentzVector.h"

#include "art/Utilities/ToolConfigTable.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SpectrumConfig.hh"
#include "Offline/Mu2eUtilities/inc/BinnedSpectrum.hh"

#include "Offline/GeneralUtilities/inc/RSNTIO.hh"

namespace mu2e {

  class ParticleGeneratorTool {
  public:

    struct Kinematic {
      PDGCode::type pdgId;
      ProcessCode creationCode;
      CLHEP::HepLorentzVector fourmom;
    };

    virtual void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string& materialName,
                                      const bool isPrimary) = 0;

    virtual std::vector<Kinematic> generate() = 0;

    // This interface should be removed when we retire ntuple-based muon resampling
    virtual void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) = 0;

    virtual std::unique_ptr<SpectrumConfig> spectrumConfig() {
      auto config = std::make_unique<SpectrumConfig>();
      return config;
    }

    virtual ~ParticleGeneratorTool() noexcept = default;

    double calculateBinnedSpectrumEnergyFraction(fhicl::ParameterSet pset,
                                                 const bool correct_full_integral = true, // correct for missing low tail to full integral
                                                 std::string var_low = "elow", // default to energy spectrum variables
                                                 std::string var_high = "ehi",
                                                 double full_var_low = 0., // if elow == ehi, it does the full spectrum
                                                 double full_var_high = 0.) const {

      // Initialize the spectra with and without the (possible) energy restriction
      BinnedSpectrum spectrum(pset);
      pset.erase(var_low);
      pset.erase(var_high);
      pset.put(var_low, full_var_low);
      pset.put(var_high, full_var_high);
      BinnedSpectrum fullSpectrum(pset);

      // Calculate the integrals
      double integral = 0.;
      double fullIntegral = 0.;
      for(size_t ibin=0;ibin < spectrum.getNbins();++ibin){
        integral += spectrum.getPDF(ibin);
      }
      for(size_t ibin=0;ibin < fullSpectrum.getNbins();++ibin){
        fullIntegral += fullSpectrum.getPDF(ibin);
      }
      // The "full" integral may be missing the lowest tail component in some cases (e.g. DIO)
      // --> apply a linear interpolation correction if requested
      const double xmin = fullSpectrum.getXMin();
      if(correct_full_integral && xmin > full_var_low) {
        const double pdfmin  = fullSpectrum.getPDF(0);
        const double binsize = fullSpectrum.getBinWidth();
        fullIntegral += 0.5*(pdfmin/binsize)*(xmin - full_var_low); // 0 to <pdfmin> over x = <spectrum low edge> to <xmin> linear interpolation
      }

      // Evaluate the fraction being sampled
      const double fraction = (fullIntegral > 0.) ? integral / fullIntegral : 0.;

      return fraction;
    }

    bool _isPrimary = true; // flag to indicate if this is for primary generation or not
  };
}

#endif

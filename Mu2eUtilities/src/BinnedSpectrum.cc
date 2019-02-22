#include "Mu2eUtilities/inc/BinnedSpectrum.hh"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "Mu2eUtilities/inc/CzarneckiSpectrum.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "Mu2eUtilities/inc/ConversionSpectrum.hh"
#include "Mu2eUtilities/inc/EjectedProtonSpectrum.hh"
#include "Mu2eUtilities/inc/MuonCaptureSpectrum.hh"
#include "Mu2eUtilities/inc/PionCaptureSpectrum.hh"
#include "DataProducts/inc/PDGCode.hh"

namespace mu2e {

  BinnedSpectrum::BinnedSpectrum(const fhicl::ParameterSet& psphys) :
    _fixMax(psphys.get<bool>("FixMax",false)),
    _finalBin(false)
  {
    const std::string spectrumShape(psphys.get<std::string>("spectrumShape"));
    if (spectrumShape == "Czarnecki") {
      // should be total energy
      PDGCode::type pdgId = PDGCode::type(psphys.get<int>("pdgId"));
      const double mass(GlobalConstantsHandle<ParticleDataTable>()->particle(pdgId).ref().mass().value());
      double elow = psphys.get<double>("elow", mass);
      double ehi  = GlobalConstantsHandle<PhysicsParams>()->getEndpointEnergy();
      double res = psphys.get<double>("spectrumResolution");
      this->initialize< CzarneckiSpectrum >(elow, ehi, res);
    }
    else if (spectrumShape == "flat") {
      double elow = psphys.get<double>("elow");
      double ehi  = psphys.get<double>("ehi");
      this->initialize<SimpleSpectrum>(elow, ehi, ehi-elow, SimpleSpectrum::Spectrum::Flat );
    }
    else if (spectrumShape == "monoenergetic") {
      this->initialize(psphys.get<double>("ehi"));
    }
    else if (spectrumShape == "CeEndpoint") {
      // think this is total energy
      double endpoint = psphys.get<double>("ehi", GlobalConstantsHandle<PhysicsParams>()->getEndpointEnergy());
      this->initialize(endpoint);
    }
    else if (spectrumShape == "ceLeadingLog") {
//-----------------------------------------------------------------------------
// ehi determines the conversion electron energy (tabulated in the .FCL file)
      // should be total energy
      double elow = psphys.get<double>("elow",0);
      double ehi  = psphys.get<double>("ehi" );
                                        // for radiatively corrected spectrum, elow and ehi are derivatives 
      double bin   = psphys.get<double>("spectrumResolution");
      // int    ratio = *ehi/bin;
      // *ehi         = (ratio+1.)*bin;

      _finalBin = true;
      this->initialize<ConversionSpectrum>(elow,ehi,bin,ehi,bin);
    }else if (spectrumShape == "ejectedProtons") {
      // should be kinetic energy
      double elow = 0.;
      double ehi = 105.; // cut off at muon mass
      double bin = (ehi - elow)/psphys.get<unsigned>("nbins");
      this->initialize<EjectedProtonSpectrum>(elow, ehi, bin);
    }else if (spectrumShape == "RMC") {
      double elow = psphys.get<double>("elow");
      double ehi = psphys.get<double>("ehi");
      double res = psphys.get<double>("spectrumResolution");
      bool kMaxUserSet = psphys.get<bool>  ("kMaxUserSet",false);
      double kMaxUser = psphys.get<double>("kMaxUser",0);
      const double bindingEnergyFit{0.464};
      const double recoilEnergyFit {0.220};
      const double deltaMassFit    {3.121};
      const double mmu = GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::mu_minus).ref().mass().value();
      const double kMaxMax =mmu - bindingEnergyFit - recoilEnergyFit - deltaMassFit;
      this->initialize<MuonCaptureSpectrum>(elow, ehi, res, kMaxUserSet, kMaxUser, kMaxMax);
    }else if (spectrumShape == "Bistirlich") {
      double elow = psphys.get<double>("elow");
      double ehi = psphys.get<double>("ehi");
      double res = psphys.get<double>("spectrumResolution");
      this->initialize<PionCaptureSpectrum>(elow, ehi, res);
    }else if (spectrumShape == "tabulated") {
//-----------------------------------------------------------------------------
// P.Murat: assume that tabulated are the bin centers
// set 'BinCenter' to false if it is the left edges
//-----------------------------------------------------------------------------
      // use elow and ehi if you want to create a truncated version of the table
      double elow = psphys.get<double>("elow",0);
      double ehi = psphys.get<double>("ehi",0);
      this->initialize(loadTable<2>( ConfigFileLookupPolicy()( psphys.get<std::string>("spectrumFileName"))),psphys.get<bool>("BinCenter", false),elow,ehi);
    
      if(_xmin < 0.0) throw cet::exception("BADCONFIG")
        << "BinnedSpectrum: negative energy endpoint "<< _xmin  <<"\n";
    }
    else {
      throw cet::exception("BADCONFIG")
        << "BinnedSpectrum: unknown spectrum shape "<<spectrumShape<<"\n";
    }
  }
}

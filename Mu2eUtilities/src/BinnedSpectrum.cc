#include <exception>
#include <string>

#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/exception.h"

#include "Offline/Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Offline/Mu2eUtilities/inc/ConversionSpectrum.hh"
#include "Offline/Mu2eUtilities/inc/CzarneckiSpectrum.hh"
#include "Offline/Mu2eUtilities/inc/EjectedProtonSpectrum.hh"
#include "Offline/Mu2eUtilities/inc/MuonCaptureSpectrum.hh"
#include "Offline/Mu2eUtilities/inc/PionCaptureSpectrum.hh"
#include "Offline/Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "Offline/Mu2eUtilities/inc/Mu2eXSpectrum.hh"
namespace mu2e {

  BinnedSpectrum::BinnedSpectrum(const fhicl::ParameterSet& psphys) :
    _fixMax(psphys.get<bool>("FixMax",false)),
    _finalBin(false)
  {
    const std::string spectrumShape(psphys.get<std::string>("spectrumShape"));
    if (spectrumShape == "Czarnecki") {
      // should be total energy
      PDGCode::type pdgId = PDGCode::type(psphys.get<int>("pdgId"));
      const double mass(GlobalConstantsHandle<ParticleDataList>()->particle(pdgId).mass());
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
    } else if (spectrumShape == "Mu2eX"){ // Sophie Middleton, 2021
      double elow = psphys.get<double>("elow",100);
      double ehi  = psphys.get<double>("ehi",105 );
      double bin = psphys.get<double>("spectrumResolution");
      this->initialize<Mu2eXSpectrum>(elow,ehi,bin,ehi,bin);
    } else if (spectrumShape == "ejectedProtons") {
      // should be kinetic energy
      double elow = 0.;
      // cut off at muon mass
      double ehi =  GlobalConstantsHandle<ParticleDataList>()->particle(PDGCode::mu_minus).mass();
      double bin = (ehi - elow)/psphys.get<unsigned>("nbins");
      this->initialize<EjectedProtonSpectrum>(elow, ehi, bin);
    }else if (spectrumShape == "RMC") {
      double elow = psphys.get<double>("elow");
      double ehi = psphys.get<double>("ehi");
      double res = psphys.get<double>("spectrumResolution");
      bool kMaxUserSet = psphys.get<bool>  ("kMaxUserSet",false);
      double kMaxUser = psphys.get<double>("kMaxUser",0);
      const double bindingEnergyFit = GlobalConstantsHandle<PhysicsParams>()->getRMCbindingEnergyFit("Al");
      const double recoilEnergyFit  = GlobalConstantsHandle<PhysicsParams>()->getRMCrecoilEnergyFit("Al");
      const double deltaMassFit     = GlobalConstantsHandle<PhysicsParams>()->getRMCdeltaMassFit("Al");
      const double mmu = GlobalConstantsHandle<ParticleDataList>()->particle(PDGCode::mu_minus).mass();
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

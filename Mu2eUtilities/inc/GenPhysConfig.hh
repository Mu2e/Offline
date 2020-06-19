#ifndef Mu2eUtilities_GenPhysConfig_hh_
#define Mu2eUtilities_GenPhysConfig_hh_

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

namespace mu2e { 

  struct GenPhysConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;

    fhicl::Atom<int> pdgId{Name("pdgId"), Comment("PDG ID number for this particle")};
    fhicl::Atom<std::string> spectrumVariable{Name("spectrumVariable"), Comment("The variable the defined spectrum is of")};
    fhicl::Atom<std::string> genId{Name("genId"), Comment("Generator ID string for this physics")};
    fhicl::Atom<std::string> spectrumShape{Name("spectrumShape"), Comment("Shape of the spectrum")};
    fhicl::OptionalAtom<std::string> spectrumFileName{Name("spectrumFileName"), Comment("File name to get spectrum from")};

    // For BinnedSpectrum class
    fhicl::Atom<bool> FixMax{Name("FixMax"), Comment("align bins so that top edge is exactly required xmax"), false};
    fhicl::OptionalAtom<double> elow{Name("elow"), Comment("Lowest energy of spectrum")};
    fhicl::OptionalAtom<double> ehi{Name("ehi"), Comment("Highest energy of spectrum")};
    fhicl::OptionalAtom<double> spectrumResolution{Name("spectrumResolution"), Comment("Resolution of BinnedSpectrum")};
    fhicl::OptionalAtom<unsigned> nbins{Name("nbins"), Comment("Number of bins in the spectrum")};
    fhicl::Atom<bool> kMaxUserSet{Name("kMaxUserSet"), Comment("True/false"), false};
    fhicl::OptionalAtom<double> kMaxUser{Name("kMaxUser"), Comment("Value of kMax for RMC spectrum")};
    fhicl::Atom<bool> BinCenter{Name("BinCenter"), Comment("True/false whether the tabulated data is defined at bin centers"), false};
  };

}

#endif

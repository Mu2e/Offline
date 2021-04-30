#ifndef Mu2eG4_constructProductionTargetMon_hh
#define Mu2eG4_constructProductionTargetMon_hh

//
// Free function. Approach borrowed from constructPS
// Constructs the downstream production target scanning monitor.
// Parent volume is the air in the target hall. Probably?
//

//Mu2e includes
#include "Mu2eG4Helper/inc/VolumeInfo.hh"
#include "G4LogicalVolume.hh"

namespace mu2e {

    class SimpleConfig;
    class PTMonPWC;

    void constructProductionTargetMon(VolumeInfo const& parent, SimpleConfig const& _config);

    // helper methods

    void constructTargetHallPWC(VolumeInfo const& motherVolume, 
                              const PTMonPWC* pwc, 
                              SimpleConfig const& _config, 
                              bool const doSurfaceCheck, 
                              int const verbosity);

    void insertOuterFrame(VolumeInfo const& container, 
                        const PTMonPWC* pwc, 
                        SimpleConfig const& _config,
                        bool const doSurfaceCheck,
                        int const verbosity);

    void insertWindows(VolumeInfo const& container, 
                     const PTMonPWC* pwc, 
                     SimpleConfig const& _config,
                     bool const doSurfaceCheck,
                     int const verbosity);

    void insertOuterGasBlocks(VolumeInfo const& container, 
                            const PTMonPWC* pwc, 
                            G4Material* gasMaterial, 
                            SimpleConfig const& _config,
                            bool const doSurfaceCheck,
                            int const verbosity);

    void insertVerticalProfileWires(VolumeInfo const& container, 
                                  const PTMonPWC* pwc, 
                                  G4Material* gasMaterial, 
                                  std::string const& wireNameSuffix, 
                                  SimpleConfig const& _config,
                                  bool const doSurfaceCheck,
                                  int const verbosity);

    void insertHorizontalProfileWires(VolumeInfo const& container, 
                                    const PTMonPWC* pwc, 
                                    G4Material* gasMaterial, 
                                    std::string const& wireNameSuffix, 
                                    SimpleConfig const& _config,
                                    bool const doSurfaceCheck,
                                    int const verbosity);

} // namespace mu2e


#endif /* Mu2eG4_constructProductionTargetMon_hh */

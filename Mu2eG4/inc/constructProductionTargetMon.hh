#ifndef Mu2eG4_constructProductionTargetMon_hh
#define Mu2eG4_constructProductionTargetMon_hh

//
// Free function. Approach borrowed from constructPS
// Constructs the downstream production target scanning monitor.
// Parent volume is the air in the target hall. Probably?
//

//Mu2e includes
#include "Mu2eG4Helper/inc/VolumeInfo.hh"

namespace mu2e {

    class VolumeInfo;
    class SimpleConfig;
    class G4LogicalVolume;
    class PTMonPWC;

    void constructProductionTargetMon(VolumeInfo const& parent, SimpleConfig const& _config);

    // helper methods
    
    void constructTargetHallPWC(VolumeInfo const& motherVolume, PTMonPWC* pwc, bool const doSurfaceCheck, int const verbosity);

    void insertOuterFrame(VolumeInfo const& container, 
                        PTMonPWC* pwc, 
                        bool const doSurfaceCheck,
                        int const verbosity);

    void insertWindows(G4LogicalVolume* windowLogical, 
                     VolumeInfo const& container, 
                     PTMonPWC* pwc, 
                     bool const doSurfaceCheck,
                     int const verbosity);

    void insertOuterGasBlocks(VolumeInfo const& container, 
                            PTMonPWC* pwc, 
                            G4Material* gasMaterial, 
                            bool const doSurfaceCheck,
                            int const verbosity);

    void insertVerticalProfileWires(VolumeInfo const& container, 
                                  PTMonPWC* pwc, 
                                  G4Material* gasMaterial, 
                                  std::string const& wireNameSuffix, 
                                  bool const doSurfaceCheck,
                                  int const verbosity);

    void insertHorizontalProfileWires(VolumeInfo const& container, 
                                    PTMonPWC* pwc, 
                                    G4Material* gasMaterial, 
                                    std::string const& wireNameSuffix, 
                                    bool const doSurfaceCheck,
                                    int const verbosity);

} // namespace mu2e


#endif /* Mu2eG4_constructProductionTargetMon_hh */

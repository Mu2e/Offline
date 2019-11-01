// Cosmic rays generator using CORSIKA

#ifndef Sources_inc_CosmicCORSIKA_hh
#define Sources_inc_CosmicCORSIKA_hh

#include <vector>


#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandFlat.h"

#include "fhiclcpp/types/ConfigurationTable.h"

#include "Mu2eUtilities/inc/VectorVolume.hh"


namespace art
{
  class Run;
}

struct Config
{
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;
  fhicl::Sequence<std::string> showerInputFiles{Name("fileNames"),Comment("List of CORSIKA binary output paths")};
  fhicl::Atom<int> firstEventNumber{Name("firstEventNumber"), Comment("Reference point on the Y axis"), 0};
  fhicl::Atom<int> firstSubRunNumber{Name("firstSubRunNumber"), Comment("Reference point on the Y axis"), 0};
  fhicl::Atom<int> maxEvents{Name("maxEvents"), Comment("Reference point on the Y axis"), 0};
  fhicl::Atom<unsigned int> runNumber{Name("runNumber"), Comment("Reference point on the Y axis"), 0};
  fhicl::Atom<std::string> module_label{Name("module_label"), Comment("Reference point on the Y axis"), ""};
  fhicl::Atom<std::string> module_type{Name("module_type"), Comment("Reference point on the Y axis"), ""};
  fhicl::Atom<bool> projectToTargetBox{Name("projectToTargetBox"), Comment("Store only events that cross the target box"), false};
  fhicl::Atom<float> showerAreaExtension{Name("showerAreaExtension"), Comment("Reference point on the Y axis"), 10000};
  fhicl::Atom<float> tOffset{Name("tOffset"), Comment("Time offset"), 0};
  fhicl::Atom<float> fluxConstant{Name("fluxConstant"), Comment("Flux constant"), 1.8e4};
  fhicl::Atom<float> targetBoxXmin{Name("targetBoxXmin"), Comment("Extension of the generation plane"), -5000};
  fhicl::Atom<float> targetBoxXmax{Name("targetBoxXmax"), Comment("Extension of the generation plane"), 5000};
  fhicl::Atom<float> targetBoxYmin{Name("targetBoxYmin"), Comment("Extension of the generation plane"), -5000};
  fhicl::Atom<float> targetBoxYmax{Name("targetBoxYmax"), Comment("Extension of the generation plane"), 5000};
  fhicl::Atom<float> targetBoxZmin{Name("targetBoxZmin"), Comment("Extension of the generation plane"), -5000};
  fhicl::Atom<float> targetBoxZmax{Name("targetBoxZmax"), Comment("Extension of the generation plane"), 5000};

};
typedef fhicl::WrappedTable<Config> Parameters;

class GenParticleCollection;

namespace mu2e {

  class CosmicCORSIKA{

    public:
      CosmicCORSIKA(const Config& conf);
      // CosmicCORSIKA(art::Run &run, CLHEP::HepRandomEngine &engine);
      ~CosmicCORSIKA();
      const float getLiveTime();
      const unsigned int getNumShowers();

      const std::map<unsigned int, int> corsikaToPdgId = {
        {1, 22}, // gamma
        {2, -11}, // e+
        {3, 11}, // e-
        {5, -13}, // mu+
        {6, 13}, // mu-
        {7, 111}, // pi0
        {8, 211}, // pi+
        {9, -211}, // pi-
        {10, 130}, // K0_L
        {11, 321}, // K+
        {12, -321}, // K-
        {13, 2112}, // n
        {14, 2212}, // p
        {15, -2212}, // pbar
        {16, 310}, // K0_S
        {17, 221}, // eta
        {18, 3122}, // Lambda
        {19, 3222}, // Sigma+
        {20, 3212}, // Sigma0
        {21, 3112}, // Sigma-
        {22, 3322}, // Cascade0
        {23, 3312}, // Cascade-
        {24, 3334}, // Omega-
        {25, -2112}, // nbar
        {26, -3122}, // Lambdabar
        {27, -3112}, // Sigma-bar
        {28, -3212}, // Sigma0bar
        {29, -3222}, // Sigma+bar
        {30, -3322}, // Cascade0bar
        {31, -3312}, // Cascade+bar
        {32, -3334}, // Omega+bar

        {50, 223}, // omega
        {51, 113}, // rho0
        {52, 213}, // rho+
        {53, -213}, // rho-
        {54, 2224}, // Delta++
        {55, 2214}, // Delta+
        {56, 2114}, // Delta0
        {57, 1114}, // Delta-
        {58, -2224}, // Delta--bar
        {59, -2214}, // Delta-bar
        {60, -2114}, // Delta0bar
        {61, -1114}, // Delta+bar
        {62, 10311}, // K*0
        {63, 10321}, // K*+
        {64, -10321}, // K*-
        {65, -10311}, // K*0bar
        {66, 12}, // nu_e
        {67, -12}, // nu_ebar
        {68, 14}, // nu_mu
        {69, -14}, // nu_mubar

        {116, 421}, // D0
        {117, 411}, // D+
        {118, -411}, // D-bar
        {119, -421}, // D0bar
        {120, 431}, // D+_s
        {121, -431}, // D-_sbar
        {122, 441}, // eta_c
        {123, 423}, // D*0
        {124, 413}, // D*+
        {125, -413}, // D*-bar
        {126, -423}, // D*0bar
        {127, 433}, // D*+_s
        {128, -433}, // D*-_s

        {130, 443}, // J/Psi
        {131, -15}, // tau+
        {132, 15}, // tau-
        {133, 16}, // nu_tau
        {134, -16}, // nu_taubar

        {137, 4122}, // Lambda+_c
        {138, 4232}, // Cascade+_c
        {139, 4132}, // Cascade0_c
        {140, 4222}, // Sigma++_c
        {141, 4212}, // Sigma+_c
        {142, 4112}, // Sigma0_c
        {143, 4322}, // Cascade'+_c
        {144, 4312}, // Cascade'0_c
        {145, 4332}, // Omega0_c
        {149, -4122}, // Lambda-_cbar
        {150, -4232}, // Cascade-_cbar
        {151, -4132}, // Cascade0_cbar
        {152, -4222}, // Sigma--_cbar
        {153, -4212}, // Sigma-_cbar
        {154, -4112}, // Sigma0_cbar
        {155, -4322}, // Cascade'-_cbar
        {156, -4312}, // Cascade'0_cbar
        {157, -4332}, // Omega0_cbar
        {161, 4224}, // Sigma*++_c
        {162, 1214}, // Sigma*+_c
        {163, 4114}, // Sigma*0_c

        {171, -4224}, // Sigma*--_cbar
        {172, -1214}, // Sigma*-_cbar
        {173, -4114}, // Sigma*0_cbar
        {176, 511}, // B0
        {177, 521}, // B+
        {178, -521}, // B-bar
        {179, -511}, // B0bar
        {180, 531}, // B0_s
        {181, -531}, // B0_sbar
        {182, 541}, // B+_c
        {183, -541}, // B-_cbar
        {184, 5122}, // Lambda0_b
        {185, 5112}, // Sigma-_b
        {186, 5222}, // Sigma+_b
        {187, 5232}, // Cascade0_b
        {188, 5132}, // Cascade-_b
        {189, 5332}, // Omega-_b
        {190, -5112}, // Lambda0_bbar
        {191, -5222}, // Sigma+_bbar
        {192, -5112}, // Sigma-_bbar
        {193, -5232}, // Cascade0_bbar
        {194, -5132}, // Cascade+_bbar
        {195, -5332} // Omega+_bbar
      };

      virtual bool generate(GenParticleCollection &);
      void openFile(FILE *f);

    private:

      bool genEvent(std::map<std::pair<int,int>, GenParticleCollection> &particles_map);
      float wrapvarBoxNo(const float var, const float low, const float high, int &boxno);
      static float wrapvar(const float var, const float low, const float high);
      std::vector<CLHEP::Hep3Vector> _targetBoxIntersections;
      std::vector<CLHEP::Hep3Vector> _worldIntersections;
      std::map<std::pair<int,int>, GenParticleCollection> _particles_map;

      GlobalConstantsHandle<ParticleDataTable> pdt;

      const float _GeV2MeV = CLHEP::GeV / CLHEP::MeV;
      const float _ns2s = CLHEP::ns / CLHEP::s;
      const float _cm2mm = CLHEP::cm / CLHEP::mm;

      CLHEP::Hep3Vector _cosmicReferencePointInMu2e;
      float _fluxConstant = 1.8e4;
      float _tOffset = 0; ///< Time offset of sample, defaults to zero (no offset) [s]
      bool _projectToTargetBox = false;
      float _showerAreaExtension = 0; ///< Extend distribution of corsika particles in x,z by this much (e.g. 1000 will extend 10 m in -x, +x, -z, and +z) [mm]

      std::string _refPointChoice;
      bool _geomInfoObtained = false;
      float _targetBoxXmin = 0;
      float _targetBoxXmax = 0;
      float _targetBoxYmin = 0;
      float _targetBoxYmax = 0;
      float _targetBoxZmin = 0;
      float _targetBoxZmax = 0;

      FILE *in = nullptr;
      int _loops = 0; // number of loops, necessary to skip the garbage data between blocks
      float _garbage;

      unsigned int _primaries = 0;
  };  // CosmicCORSIKA

}

#endif
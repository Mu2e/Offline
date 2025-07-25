//
// Instantiation of RegrowKalSeed for LoopHelix fits
//
// Original author: D. Brown (LBNL) 4/18/2025
#include "Offline/Mu2eKinKal/inc/RegrowKalSeed_module.hh"
namespace mu2e {
  using KinKal::VEC3;
  using KinKal::DMAT;
  using KinKal::DVEC;
  using KinKal::TimeDir;
  using MatEnv::DetMaterial;
  using KKConfig = Mu2eKinKal::KinKalConfig;
  using Mu2eKinKal::KKFinalConfig;
  using KKFitConfig = Mu2eKinKal::KKFitConfig;
  using KKModuleConfig = Mu2eKinKal::KKModuleConfig;
  struct RegrowLoopHelixConfig {
    fhicl::Table<KKFitConfig> kkfitSettings { Name("KKFitSettings") };
    fhicl::Table<KKConfig> fitSettings { Name("FitSettings") };
    fhicl::Table<KKConfig> extSettings { Name("ExtensionSettings") };
    fhicl::OptionalTable<KKExtrapConfig> Extrapolation { Name("Extrapolation") };

  };

  class RegrowLoopHelix : public RegrowKalSeed<KinKal::LoopHelix> {
    public:
      using Parameters = art::EDProducer::Table<RegrowLoopHelixConfig>;
      using KTRAJ = KinKal::LoopHelix;
      using PTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using KKTRK = KKTrack<KTRAJ>;
      using KKTRKCOL = OwningPointerCollection<KKTRK>;
      using KKSTRAWHIT = KKStrawHit<KTRAJ>;
      using KKSTRAWHITPTR = std::shared_ptr<KKSTRAWHIT>;
      using KKSTRAWHITCOL = std::vector<KKSTRAWHITPTR>;
      using KKSTRAWXING = KKStrawXing<KTRAJ>;
      using KKSTRAWXINGPTR = std::shared_ptr<KKSTRAWXING>;
      using KKSTRAWXINGCOL = std::vector<KKSTRAWXINGPTR>;
      using KKIPAXING = KKShellXing<KTRAJ,KinKal::Cylinder>;
      using KKIPAXINGPTR = std::shared_ptr<KKIPAXING>;
      using KKIPAXINGCOL = std::vector<KKIPAXINGPTR>;
      using KKSTXING = KKShellXing<KTRAJ,KinKal::Annulus>;
      using KKSTXINGPTR = std::shared_ptr<KKSTXING>;
      using KKSTXINGCOL = std::vector<KKSTXINGPTR>;
      using KKCALOHIT = KKCaloHit<KTRAJ>;
      using KKCALOHITPTR = std::shared_ptr<KKCALOHIT>;
      using KKCALOHITCOL = std::vector<KKCALOHITPTR>;
      using KKFIT = KKFit<KTRAJ>;

      using MEAS = KinKal::Hit<KTRAJ>;
      using MEASPTR = std::shared_ptr<MEAS>;
      using MEASCOL = std::vector<MEASPTR>;
      using EXING = KinKal::ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
      using EXINGCOL = std::vector<EXINGPTR>;

      using KKMaterialConfig = KKMaterial::Config;

      explicit RegrowLoopHelis(const Parameters& settings);
      void beginRun(art::Run& run) override;
      void produce(art::Event& event) override;
      void endJob() override;
    private:
      produces<KKTRKCOL>();
      produces<KalSeedCollection>();
  };
}
DEFINE_ART_MODULE(mu2e::RegrowLoopHelix)

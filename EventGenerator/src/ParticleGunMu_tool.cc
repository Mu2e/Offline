//-----------------------------------------------------------------------------
// Nov 2023 P.Murat: more or less universal particle gun plugin
//-----------------------------------------------------------------------------
#include "art/Utilities/ToolMacros.h"
#include <memory>
#include <vector>

#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandFlat.h"

#include "Offline/EventGenerator/inc/ParticleGeneratorTool.hh"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"

#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/Sequence.h"

namespace mu2e {
  class ParticleGunMu : public ParticleGeneratorTool {
  private:
    int               _nParticles;
    PDGCode::type     _pdgCode;
    double            _m;
    double            _mean;
    double            _pmin;
    double            _pmax;
    double            _czmin;
    double            _czmax;
    double            _phimin;
    double            _phimax;
    std::string       _shape;          // "box" or "sphere"
    CLHEP::Hep3Vector _center;
    CLHEP::Hep3Vector _dimensions;        // for box: dX/2, dY/2, dZ/2

    // BinnedSpectrum    _spectrum;

    std::unique_ptr<CLHEP::RandFlat>        _randFlat        ;
    std::unique_ptr<RandomUnitSphere>       _randomUnitSphere;
    std::unique_ptr<CLHEP::RandGeneral>     _randSpectrum    ;

  public:
    struct PhysConfig {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int>               nParticles  {Name("nParticles"  ), Comment("N particles"   ), 1};
      fhicl::Atom<int>               pdgCode     {Name("pdgCode"     ), Comment("PDG code"      )};
      fhicl::Atom<double>            pmin        {Name("pmin"        ), Comment("pmin"          )};
      fhicl::Atom<double>            pmax        {Name("pmax"        ), Comment("pmax"          )};
      fhicl::Atom<std::string>       mode        {Name("mode"        ), Comment("mode: flat, or smth else")};
      fhicl::Atom<double>            czmin       {Name("czmin"       ), Comment("cos(theta) min"),  -1};
      fhicl::Atom<double>            czmax       {Name("czmax"       ), Comment("cos(theta) max"),   1};
      fhicl::Atom<double>            phimin      {Name("phimin"      ), Comment("phi min"       ),   0};
      fhicl::Atom<double>            phimax      {Name("phimax"      ), Comment("phi max"       ), 360};
      fhicl::Atom<std::string>       shape       {Name("shape"       ), Comment("shape: box/sphere")};
      fhicl::Sequence<float>         center      {Name("center"      ), Comment("X,Y,Z"         )};
      fhicl::Sequence<float>         dimensions  {Name("dimensions"  ), Comment("dimensions: dX/2,dY/2,dZ/2 or rX,dY,rZ")};
    };
    typedef art::ToolConfigTable<PhysConfig> Parameters;

    explicit ParticleGunMu(Parameters const& conf);

    std::vector<ParticleGeneratorTool::Kinematic> generate() override;
    void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop);

    virtual GenId         genId      () override { return GenId::particleGun;               }
    virtual ProcessCode   processCode() override { return ProcessCode::mu2eMuonDecayAtRest; }

    void   getMom(CLHEP::HepLorentzVector* Mom);
    void   getXYZ(CLHEP::Hep3Vector*       Xyz);

    virtual void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string&) override {
      _randFlat         = std::make_unique<CLHEP::RandFlat>(eng);
      _randomUnitSphere = std::make_unique<RandomUnitSphere>(eng);
      //      _randSpectrum     = std::make_unique<CLHEP::RandGeneral>(eng, _spectrum.getPDF(), _spectrum.getNbins());
    }
  };


//-----------------------------------------------------------------------------
  ParticleGunMu::ParticleGunMu(Parameters const& conf) :
    _nParticles  (conf().nParticles  ()),
    _pdgCode     (static_cast<PDGCode::type>(conf().pdgCode())),
    _pmin        (conf().pmin        ()),
    _pmax        (conf().pmax        ()),
    _czmin       (conf().czmin       ()),
    _czmax       (conf().czmax       ()),
    _phimin      (conf().phimin      ()),
    _phimax      (conf().phimax      ()),
    _shape       (conf().shape       ())

    // _mass(GlobalConstantsHandle<ParticleDataList>()->particle(_pdgId).mass()),
    // _spectrum(BinnedSpectrum(conf().spectrum.get<fhicl::ParameterSet>()))
  {
    GlobalConstantsHandle<ParticleDataList> pdt;
    _m  = pdt->particle(_pdgCode).mass();

    std::vector v(conf().center());
    _center.set(v[0],v[1],v[2]);

    std::vector d(conf().dimensions());
    _dimensions.set(d[0],d[1],d[2]);

  }

//-----------------------------------------------------------------------------
  std::vector<ParticleGeneratorTool::Kinematic> ParticleGunMu::generate() {
    std::vector<ParticleGeneratorTool::Kinematic>  res;

    for (int i=0; i<_nParticles; i++) {
      double p = _pmin + (_pmax-_pmin) * _randFlat->fire(0,1);
      double e = sqrt(p*p+_m*_m);

      CLHEP::HepLorentzVector p4(_randomUnitSphere->fire(p), e);

      ParticleGeneratorTool::Kinematic k{_pdgCode, processCode(),p4};
      res.emplace_back(k);
    }

    return res;
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void ParticleGunMu::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
  }

//-----------------------------------------------------------------------------
  void ParticleGunMu::getXYZ(CLHEP::Hep3Vector* Xyz) {

    double x[3];
    for (int i=0; i<3; ++i) {
      x[i] = _center[i]+_randFlat->fire(-1.,1.)*_dimensions[i];
    }

    Xyz->set(x[0],x[1],x[2]);
  }

//-----------------------------------------------------------------------------
// generate 4-momentum of the next particle
//-----------------------------------------------------------------------------
  void ParticleGunMu::getMom(CLHEP::HepLorentzVector* Mom) {

    double p = _pmin+_randFlat->fire(0.,1.)*(_pmax-_pmin);
    double e = sqrt(p*p+_m*_m);

    CLHEP::Hep3Vector p3 = _randomUnitSphere->fire(p);

    Mom->set(p3[0],p3[1],p3[2],e);
  }
}

DEFINE_ART_CLASS_TOOL(mu2e::ParticleGunMu)

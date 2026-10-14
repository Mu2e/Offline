///////////////////////////////////////////////////////////////////////////////
// AprHelixTruthAna: per-event TTree comparing helix seeds (e.g. the trigger APR
// helices) with the MC truth of the main signal particle, to study helix finding
// efficiency and wrong-lambda / wrong-circle helices. Optionally records the
// decisions of trigger paths of the current process (triggerPaths), so trigger
// efficiencies can be split by the truth of the signal particle.
//
// The signal particle is the highest-momentum e+- or mu+- with >= minHits straw hits
// (via the StrawDigiMCs). Its truth at the earliest StrawGasStep gives p, pT, pz, the
// circle center and radius, and lambda = pz/(c B). Its combo hits in the
// helix finder's input collection give the z span and the largest z gap.
// For each helix: lambda, radius, center, N(straw hits), N(signal straw hits),
// z span and largest z gap of its hits.
///////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"

#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/Mu2eUtilities/inc/TriggerResultsNavigator.hh"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "fhiclcpp/types/Sequence.h"

#include "CLHEP/Units/PhysicalConstants.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

namespace mu2e {

  class AprHelixTruthAna : public art::EDAnalyzer {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> hsTag {Name("helixSeedTag"), Comment("helix seed collection")};
      fhicl::Atom<art::InputTag> chTag {Name("comboHitTag"), Comment("helix finder input combo hits")};
      fhicl::Atom<art::InputTag> sdmcTag {Name("strawDigiMCTag"), Comment("StrawDigiMC collection")};
      fhicl::Atom<int> minHits {Name("minHits"), Comment("min straw hits of the signal particle"), 5};
      fhicl::Sequence<std::string> paths {Name("triggerPaths"), Comment("trigger paths of this process to record (pass = 1)"), std::vector<std::string>{}};
    };
    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit AprHelixTruthAna(const Parameters& config);
    void beginJob() override;
    void beginRun(const art::Run&) override;
    void analyze(const art::Event& event) override;

  private:
    // z span, largest gap and N(distinct z planes) of a set of z values
    static void zInfo(std::vector<float> z, float& zmin, float& zmax, float& gap, int& nz);

    art::InputTag _hsTag, _chTag, _sdmcTag;
    int _minHits;
    std::vector<std::string> _paths;
    std::vector<int> _pass; // per entry of _paths
    float _bzConv = 0.; // mm per MeV/c
    TTree* _tree = nullptr;

    // event
    int _run, _subrun, _event;
    bool _hasHelices; // the helix collection exists (its producer ran)
    // truth of the signal particle
    int _pdg, _nsh; // pdg, N(straw hits)
    float _p, _pt, _pz, _rtrue, _lambdatrue, _xctrue, _yctrue, _d0true;
    int _nchsig, _nzsig; // signal combo hits in the input collection, distinct z
    float _zminsig, _zmaxsig, _gapsig;
    // helices
    std::vector<float> _hlambda, _hradius, _hxc, _hyc, _hd0, _hzmin, _hzmax, _hgap, _hchi2xy, _hchi2zphi;
    std::vector<int> _hnsh, _hnshsig, _hnz;
    // circle chi2/ndof of the helix's hits for the reco and the true circle (finder error model),
    // and N(signal combo hits in the input collection) within 3 sigma of each circle
    std::vector<float> _hchi2creco, _hchi2ctrue;
    std::vector<int> _hnsig3reco, _hnsig3true;
    int _nsig3true;
    // signal combo hits: position, N(straw hits), bit mask of the helices containing it
    std::vector<float> _sx, _sy, _sz;
    std::vector<int> _snsh, _shel;
    static float circleChi2(const ComboHit& hit, float xc, float yc, float r);
  };

  AprHelixTruthAna::AprHelixTruthAna(const Parameters& config) :
    art::EDAnalyzer{config}, _hsTag(config().hsTag()), _chTag(config().chTag()), _sdmcTag(config().sdmcTag()), _minHits(config().minHits()), _paths(config().paths()), _pass(_paths.size(), 0) {}

  void AprHelixTruthAna::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    _tree = tfs->make<TTree>("helices", "helix seeds vs signal truth");
    _tree->Branch("run", &_run);
    _tree->Branch("subrun", &_subrun);
    _tree->Branch("event", &_event);
    _tree->Branch("hasHelices", &_hasHelices);
    _tree->Branch("pdg", &_pdg);
    _tree->Branch("nsh", &_nsh);
    _tree->Branch("p", &_p);
    _tree->Branch("pt", &_pt);
    _tree->Branch("pz", &_pz);
    _tree->Branch("rtrue", &_rtrue);
    _tree->Branch("lambdatrue", &_lambdatrue);
    _tree->Branch("xctrue", &_xctrue);
    _tree->Branch("yctrue", &_yctrue);
    _tree->Branch("d0true", &_d0true);
    _tree->Branch("nchsig", &_nchsig);
    _tree->Branch("nzsig", &_nzsig);
    _tree->Branch("zminsig", &_zminsig);
    _tree->Branch("zmaxsig", &_zmaxsig);
    _tree->Branch("gapsig", &_gapsig);
    _tree->Branch("hlambda", &_hlambda);
    _tree->Branch("hradius", &_hradius);
    _tree->Branch("hxc", &_hxc);
    _tree->Branch("hyc", &_hyc);
    _tree->Branch("hd0", &_hd0);
    _tree->Branch("hzmin", &_hzmin);
    _tree->Branch("hzmax", &_hzmax);
    _tree->Branch("hgap", &_hgap);
    _tree->Branch("hchi2xy", &_hchi2xy);
    _tree->Branch("hchi2zphi", &_hchi2zphi);
    _tree->Branch("hnsh", &_hnsh);
    _tree->Branch("hnshsig", &_hnshsig);
    _tree->Branch("hnz", &_hnz);
    _tree->Branch("hchi2creco", &_hchi2creco);
    _tree->Branch("hchi2ctrue", &_hchi2ctrue);
    _tree->Branch("hnsig3reco", &_hnsig3reco);
    _tree->Branch("hnsig3true", &_hnsig3true);
    _tree->Branch("nsig3true", &_nsig3true);
    _tree->Branch("sx", &_sx);
    _tree->Branch("sy", &_sy);
    _tree->Branch("sz", &_sz);
    _tree->Branch("snsh", &_snsh);
    _tree->Branch("shel", &_shel);
    for (size_t i = 0; i < _paths.size(); ++i) _tree->Branch(("trig_" + _paths[i]).c_str(), &_pass[i]);
  }

  void AprHelixTruthAna::beginRun(const art::Run&) {
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    const double bz = bfmgr->getBField(det->toMu2e(CLHEP::Hep3Vector(0., 0., 0.))).z();
    _bzConv = bz * CLHEP::c_light / 1000.; // MeV/c per mm
  }

  // squared circle residual in units of the AgnosticHelixFinder circle error
  float AprHelixTruthAna::circleChi2(const ComboHit& hit, float xc, float yc, float r) {
    const float dx = hit.pos().x() - xc, dy = hit.pos().y() - yc;
    const float d2 = dx*dx + dy*dy;
    const float dxn = dx * hit.vDir().x() + dy * hit.vDir().y();
    const float costh2 = dxn * dxn / d2;
    const float err2 = hit.wireVar() * (1.f - costh2) + hit.transVar() * costh2;
    const float res = r - std::sqrt(d2);
    return res * res / err2;
  }

  void AprHelixTruthAna::zInfo(std::vector<float> z, float& zmin, float& zmax, float& gap, int& nz) {
    zmin = zmax = gap = 0.;
    nz = 0;
    if (z.empty()) return;
    std::sort(z.begin(), z.end());
    zmin = z.front();
    zmax = z.back();
    nz = 1;
    for (size_t i = 1; i < z.size(); ++i) {
      const float dz = z[i] - z[i-1];
      gap = std::max(gap, dz);
      if (dz > 10.) ++nz; // a new plane face
    }
  }

  void AprHelixTruthAna::analyze(const art::Event& event) {
    _run = event.run(); _subrun = event.subRun(); _event = event.event();
    if (!_paths.empty()) {
      art::Handle<art::TriggerResults> trH;
      event.getByLabel(art::InputTag("TriggerResults", "", moduleDescription().processName()), trH);
      for (size_t i = 0; i < _paths.size(); ++i) {
        _pass[i] = 0;
        if (trH.isValid()) {
          TriggerResultsNavigator nav(trH.product());
          if (nav.validPath(_paths[i])) _pass[i] = nav.accepted(_paths[i]);
        }
      }
    }
    const auto& sdmcs = *event.getValidHandle<StrawDigiMCCollection>(_sdmcTag);
    art::Handle<ComboHitCollection> chH;
    event.getByLabel(_chTag, chH);
    static const ComboHitCollection noHits;
    const auto& chcol = chH.isValid() ? *chH : noHits;
    // the helix finder may not have run (an earlier filter in its path failed)
    art::Handle<HelixSeedCollection> hsH;
    event.getByLabel(_hsTag, hsH);
    _hasHelices = hsH.isValid();
    static const HelixSeedCollection empty;
    const auto& hscol = _hasHelices ? *hsH : empty;

    // signal particle: the highest-momentum e+- or mu+- (at its earliest straw hit) with >= minHits straw hits
    // (the conversion electron in CE samples; pileup loopers have more hits but low momentum)
    std::map<art::Ptr<SimParticle>, std::pair<int, std::pair<float, float>>> info; // N(hits), (earliest time, p)
    for (const auto& mc : sdmcs) {
      const auto& sgs = mc.earlyStrawGasStep();
      if (sgs.isNull()) continue;
      const auto& simp = sgs->simParticle();
      if (std::abs(simp->pdgId()) != 11 && std::abs(simp->pdgId()) != 13) continue; // e+- and mu+-
      auto& in = info[simp];
      if (in.first == 0 || sgs->time() < in.second.first) in.second = {float(sgs->time()), float(sgs->momvec().mag())};
      ++in.first;
    }
    art::Ptr<SimParticle> sig;
    _nsh = 0;
    float pbest = 0.;
    for (const auto& [simp, in] : info) {
      if (in.first < _minHits || in.second.second <= pbest) continue;
      pbest = in.second.second;
      sig = simp;
      _nsh = in.first;
    }

    // truth at the earliest step of the signal particle
    _pdg = 0; _p = _pt = _pz = _rtrue = _lambdatrue = _xctrue = _yctrue = _d0true = 0.;
    float tmin = 1.e12;
    for (const auto& mc : sdmcs) {
      const auto& sgs = mc.earlyStrawGasStep();
      if (sgs.isNull() || sgs->simParticle() != sig || sgs->time() >= tmin) continue;
      tmin = sgs->time();
      const auto mom = sgs->momvec();
      const auto pos = sgs->position();
      _pdg = sig->pdgId();
      _p = mom.mag(); _pt = mom.perp(); _pz = mom.z();
      const float q = (sig->pdgId() == 11 || sig->pdgId() == 13 || sig->pdgId() == -211 || sig->pdgId() == -2212) ? -1. : 1.;
      _rtrue = _pt / std::fabs(_bzConv);
      _lambdatrue = _pz / std::fabs(_bzConv);
      // circle center: pos + R * q * sign(B) * (py, -px)/pT
      const float s = q * (_bzConv > 0 ? 1. : -1.);
      _xctrue = pos.x() + _rtrue * s * mom.y() / _pt;
      _yctrue = pos.y() - _rtrue * s * mom.x() / _pt;
      _d0true = std::sqrt(_xctrue*_xctrue + _yctrue*_yctrue) - _rtrue;
    }

    // signal combo hits in the helix finder input collection: hits with a majority of signal straws
    std::vector<float> zsig;
    std::vector<size_t> sigHits; // signal combo hits in the input collection
    std::vector<StrawHitIndex> shids;
    for (size_t ich = 0; ich < chcol.size(); ++ich) {
      shids.clear();
      chcol.fillStrawHitIndices(ich, shids);
      int nsig = 0;
      for (auto sh : shids) {
        const auto& sgs = sdmcs.at(sh).earlyStrawGasStep();
        if (sgs.isNonnull() && sgs->simParticle() == sig) ++nsig;
      }
      if (2 * nsig > (int)shids.size()) { zsig.push_back(chcol[ich].pos().z()); sigHits.push_back(ich); }
    }
    _nchsig = zsig.size();
    _sx.clear(); _sy.clear(); _sz.clear(); _snsh.clear(); _shel.clear();
    for (auto ich : sigHits) {
      const auto& hit = chcol[ich];
      _sx.push_back(hit.pos().x()); _sy.push_back(hit.pos().y()); _sz.push_back(hit.pos().z());
      _snsh.push_back(hit.nStrawHits());
      int mask = 0;
      for (size_t ih = 0; ih < hscol.size() && ih < 30; ++ih) {
        for (const auto& hh : hscol[ih].hits()) {
          if ((hh.pos() - hit.pos()).Mag2() < 1.e-6) { mask |= (1 << ih); break; }
        }
      }
      _shel.push_back(mask);
    }
    _nsig3true = 0;
    for (auto ich : sigHits) if (circleChi2(chcol[ich], _xctrue, _yctrue, _rtrue) < 9.f) ++_nsig3true;
    zInfo(zsig, _zminsig, _zmaxsig, _gapsig, _nzsig);

    // helices
    _hlambda.clear(); _hradius.clear(); _hxc.clear(); _hyc.clear(); _hd0.clear(); _hzmin.clear(); _hzmax.clear();
    _hgap.clear(); _hchi2xy.clear(); _hchi2zphi.clear(); _hnsh.clear(); _hnshsig.clear(); _hnz.clear();
    _hchi2creco.clear(); _hchi2ctrue.clear(); _hnsig3reco.clear(); _hnsig3true.clear();
    for (const auto& hs : hscol) {
      const auto& h = hs.helix();
      _hlambda.push_back(h.lambda());
      _hradius.push_back(h.radius());
      _hxc.push_back(h.centerx());
      _hyc.push_back(h.centery());
      _hd0.push_back(h.rcent() - h.radius());
      _hchi2xy.push_back(h.chi2dXY());
      _hchi2zphi.push_back(h.chi2dZPhi());
      const auto& hits = hs.hits();
      std::vector<float> z;
      int nsh = 0, nshsig = 0;
      for (size_t i = 0; i < hits.size(); ++i) {
        z.push_back(hits[i].pos().z());
        shids.clear();
        hits.fillStrawHitIndices(i, shids);
        for (auto sh : shids) {
          ++nsh;
          const auto& sgs = sdmcs.at(sh).earlyStrawGasStep();
          if (sgs.isNonnull() && sgs->simParticle() == sig) ++nshsig;
        }
      }
      float zmin, zmax, gap;
      int nz;
      zInfo(z, zmin, zmax, gap, nz);
      _hzmin.push_back(zmin); _hzmax.push_back(zmax); _hgap.push_back(gap); _hnz.push_back(nz);
      _hnsh.push_back(nsh); _hnshsig.push_back(nshsig);
      float c2r = 0., c2t = 0.;
      for (size_t i = 0; i < hits.size(); ++i) {
        c2r += circleChi2(hits[i], h.centerx(), h.centery(), h.radius());
        c2t += circleChi2(hits[i], _xctrue, _yctrue, _rtrue);
      }
      const float ndof = std::max<float>(1.f, hits.size() - 3.f);
      _hchi2creco.push_back(c2r / ndof);
      _hchi2ctrue.push_back(c2t / ndof);
      int n3r = 0, n3t = 0;
      for (auto ich : sigHits) {
        if (circleChi2(chcol[ich], h.centerx(), h.centery(), h.radius()) < 9.f) ++n3r;
        if (circleChi2(chcol[ich], _xctrue, _yctrue, _rtrue) < 9.f) ++n3t;
      }
      _hnsig3reco.push_back(n3r);
      _hnsig3true.push_back(n3t);
    }
    _tree->Fill();
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::AprHelixTruthAna)

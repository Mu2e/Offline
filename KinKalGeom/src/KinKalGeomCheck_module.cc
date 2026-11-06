//
// Check the KinKalGeom passive material: (1) no two solid elements (DS shells, passive planes, CRV sectors) share a
// volume and no two passive planes are coplanar, of the same material and thickness, with overlapping or touching
// footprints, and no line crosses a layer of two CRV sectors; (2) given a GDML file of the G4 geometry, points sampled
// inside each element lie in G4 volumes of a matching material, the CRV layer spans match the G4 scintillator, and the
// column density along straight lines from the tracker to the CRV sectors, outside the DS bore, agrees between
// KinKalGeom and G4. G4 material crossed by those lines outside every KinKalGeom element is
// tallied per G4 volume. The job fails (throws at beginRun) on (1), and on (2) if failOnG4Mismatch is set.
//
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Run.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include "Offline/KinKalGeom/inc/KKMaterial.hh"
#include "Offline/KinKalGeom/inc/KinKalGeomCheck.hh"
#include "Offline/KinKalGeom/inc/Tracker.hh"
#include "Offline/KinKalGeom/inc/DetectorSolenoid.hh"
#include "Offline/KinKalGeom/inc/CRV.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoMaterial.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

namespace mu2e {
  class KinKalGeomCheck : public art::EDAnalyzer {
    public:
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      struct Config {
        fhicl::Atom<int> debug{Name("debug"), Comment("1: list the elements, 2: also every bad ray"), 0};
        fhicl::Atom<unsigned> seed{Name("seed"), Comment("random seed"), 12345};
        fhicl::Atom<unsigned> nOverlapSample{Name("nOverlapSample"), Comment("points per element for the overlap test"), 20000};
        fhicl::Atom<double> minOverlapVolume{Name("minOverlapVolume"), Comment("overlap volume (mm^3) above which a pair is reported"), 1000.0};
        fhicl::Atom<double> coplanarTolerance{Name("coplanarTolerance"), Comment("distance (mm) under which planes are coplanar, and footprint gap under which they touch"), 0.1};
        fhicl::Atom<bool> allowCRVSectorOverlaps{Name("allowCRVSectorOverlaps"), Comment("skip CRV sector pairs in the volume test: their staggered layers interleave, and their layer fractions are checked instead"), true};
        fhicl::Atom<double> layerSumStep{Name("layerSumStep"), Comment("step (mm) across the shared band of overlapping CRV sectors for the layer-fraction test"), 1.0};
        fhicl::Atom<std::string> crvScintillator{Name("crvScintillator"), Comment("G4 material of the CRV scintillator bars"), "G4_POLYSTYRENE"};
        fhicl::Atom<double> minLayerAgreement{Name("minLayerAgreement"), Comment("fraction of points on CRV layer center planes where KinKal and G4 agree on scintillator"), 0.95};
        fhicl::Atom<bool> failOnOverlaps{Name("failOnOverlaps"), Comment("throw on overlaps or coplanar planes"), true};
        fhicl::Atom<std::string> gdmlFile{Name("gdmlFile"), Comment("GDML file of the G4 geometry; empty skips the G4 comparison"), ""};
        fhicl::Sequence<double,3> gdmlOffset{Name("gdmlOffset"), Comment("Mu2e origin in the GDML world frame (mm)"), std::array<double,3>{0.0,0.0,0.0}};
        fhicl::Sequence<std::string> materialAliases{Name("materialAliases"), Comment("'KinKalName=G4Name' pairs of matching materials"),
          std::vector<std::string>{"CRVStrongbackAl=G4_Al","DSShieldAl=G4_Al","DSCoilMix=DS1CoilMix","DSCoilMix=DS2CoilMix"}};
        fhicl::Atom<unsigned> nMaterialSample{Name("nMaterialSample"), Comment("points per element for the material test"), 2000};
        fhicl::Atom<double> minMaterialFraction{Name("minMaterialFraction"), Comment("fraction of an element's points that must be in matching G4 material"), 0.8};
        fhicl::Atom<unsigned> nRays{Name("nRays"), Comment("tracker-to-CRV lines"), 5000};
        fhicl::Atom<double> minDensity{Name("minDensity"), Comment("G4 materials below this density (g/cm^3) are ignored"), 0.01};
        fhicl::Atom<double> coverageStep{Name("coverageStep"), Comment("step (mm) of the coverage sampling along G4 steps"), 5.0};
        fhicl::Atom<double> rayAbsTolerance{Name("rayAbsTolerance"), Comment("column density difference (g/cm^2) always accepted"), 2.0};
        fhicl::Atom<double> rayRelTolerance{Name("rayRelTolerance"), Comment("relative column density difference accepted"), 0.1};
        fhicl::Atom<double> maxBadRayFraction{Name("maxBadRayFraction"), Comment("fraction of lines allowed outside the tolerance"), 0.05};
        fhicl::Atom<bool> failOnG4Mismatch{Name("failOnG4Mismatch"), Comment("throw when the G4 comparison fails"), false};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit KinKalGeomCheck(Parameters const& config);
      void beginRun(art::Run const& run) override;
      void analyze(art::Event const&) override {}
    private:
      using VEC3 = KinKal::VEC3;
      bool checkInternal(KinKalGeom const& kkg, KKGeomCheck::Elements const& elems) const;
      bool checkG4(KinKalGeom const& kkg, KKGeomCheck::Elements const& elems) const;
      std::string g4Material(VEC3 const& detpos, std::string* volname = nullptr) const;
      VEC3 toGDML(VEC3 const& detpos) const;  // detector coordinates (mm) -> GDML world (cm)
      Config config_;
      std::multimap<std::string,std::string> aliases_;
  };

  KinKalGeomCheck::KinKalGeomCheck(Parameters const& config) : art::EDAnalyzer(config), config_(config()) {
    for(auto const& alias : config_.materialAliases()){
      auto const eq = alias.find('=');
      if(eq == std::string::npos) throw cet::exception("CONFIG") << "KinKalGeomCheck: material alias '" << alias << "' is not 'KinKalName=G4Name'\n";
      aliases_.emplace(alias.substr(0,eq),alias.substr(eq+1));
    }
  }

  void KinKalGeomCheck::beginRun(art::Run const&) {
    GeomHandle<KinKalGeom> kkg_h;
    auto const elems = KKGeomCheck::elements(*kkg_h);
    GeomHandle<DetectorSystem> det;
    if(config_.debug() > 0){
      std::cout << "KinKalGeomCheck: " << elems.size() << " elements (Mu2e coordinates, mm)" << std::endl;
      for(auto const& el : elems){
        auto const c = det->toMu2e(CLHEP::Hep3Vector(el.center.X(),el.center.Y(),el.center.Z()));
        std::cout << "  " << std::setw(22) << std::left << el.name << std::right << " " << std::setw(16) << el.material
          << " center (" << c.x() << "," << c.y() << "," << c.z() << ") w " << el.wdir;
        if(el.shape == KKGeomCheck::Element::box) std::cout << " halves u " << el.uhalf << " v " << el.vhalf << " w " << el.whalf << " udir " << el.udir;
        else std::cout << " r " << el.rin << "-" << el.rout << " half length " << el.whalf;
        std::cout << std::endl;
      }
    }
    bool const internalOK = checkInternal(*kkg_h,elems);
    bool const g4OK = config_.gdmlFile().empty() ? true : checkG4(*kkg_h,elems);
    std::cout << "KinKalGeomCheck: internal " << (internalOK ? "PASS" : "FAIL")
      << ", G4 " << (config_.gdmlFile().empty() ? "not run" : (g4OK ? "PASS" : "FAIL")) << std::endl;
    if((!internalOK && config_.failOnOverlaps()) || (!g4OK && config_.failOnG4Mismatch()))
      throw cet::exception("KINKALGEOM") << "KinKalGeomCheck failed\n";
  }

  bool KinKalGeomCheck::checkInternal(KinKalGeom const& kkg, KKGeomCheck::Elements const& elems) const {
    auto const ovls = KKGeomCheck::overlaps(elems,config_.nOverlapSample(),config_.seed(),config_.minOverlapVolume(),config_.allowCRVSectorOverlaps());
    std::cout << "KinKalGeomCheck: " << ovls.size() << " overlapping element pairs" << std::endl;
    for(auto const& ovl : ovls)
      std::cout << "  OVERLAP " << elems[ovl.first].name << " / " << elems[ovl.second].name << " volume " << ovl.volume/1000.0 << " cm^3" << std::endl;
    auto const copl = KKGeomCheck::coplanarPlanes(kkg,config_.coplanarTolerance());
    auto const& planes = kkg.passiveMaterialPlanes();
    std::cout << "KinKalGeomCheck: " << copl.size() << " coplanar same-material plane pairs" << std::endl;
    for(auto const& cp : copl)
      std::cout << "  COPLANAR " << planes[cp.first].sid_ << " / " << planes[cp.second].sid_ << " " << planes[cp.first].material_
        << " in-plane overlaps " << cp.uoverlap << ", " << cp.voverlap << " mm" << std::endl;
    // a line through the band where two sectors' staggered layers interleave crosses each layer in only one of them
    auto const sums = KKGeomCheck::sectorLayerSums(kkg,config_.coplanarTolerance(),config_.layerSumStep());
    unsigned nbadsum = 0;
    std::cout << "KinKalGeomCheck: " << sums.size() << " CRV sector pairs with overlapping envelopes (summed layer fraction, min and max)" << std::endl;
    for(auto const& sp : sums){
      bool const bad = sp.maxsum > 1.0 + 1.0e-9;
      if(bad) ++nbadsum;
      std::cout << (bad ? "  BAD " : "  ") << kkg.CRV()->sectorName(sp.first) << " / " << kkg.CRV()->sectorName(sp.second)
        << " " << sp.minsum << " " << sp.maxsum << " (" << sp.nlines << " lines)" << std::endl;
    }
    return ovls.empty() && copl.empty() && nbadsum == 0;
  }

  KinKal::VEC3 KinKalGeomCheck::toGDML(VEC3 const& detpos) const {
    GeomHandle<DetectorSystem> det;
    auto const m = det->toMu2e(CLHEP::Hep3Vector(detpos.X(),detpos.Y(),detpos.Z()));
    auto const& off = config_.gdmlOffset();
    return VEC3(0.1*(m.x()+off[0]),0.1*(m.y()+off[1]),0.1*(m.z()+off[2]));
  }

  std::string KinKalGeomCheck::g4Material(VEC3 const& detpos, std::string* volname) const {
    auto const g = toGDML(detpos);
    TGeoNode* node = gGeoManager->FindNode(g.X(),g.Y(),g.Z());
    if(!node) return "";
    if(volname) *volname = node->GetVolume()->GetName();
    std::string name = node->GetVolume()->GetMaterial()->GetName();
    return name.substr(0,name.find("0x"));
  }

  bool KinKalGeomCheck::checkG4(KinKalGeom const& kkg, KKGeomCheck::Elements const& elems) const {
    TGeoManager::SetVerboseLevel(0);
    if(!TGeoManager::Import(config_.gdmlFile().c_str()) || !gGeoManager)
      throw cet::exception("KINKALGEOM") << "KinKalGeomCheck: cannot import " << config_.gdmlFile() << "\n";
    GeomHandle<KKMaterial> kkmat;
    // KinKal densities (DetMaterial: g/mm^3) in g/cm^3, per element
    std::vector<double> rho(elems.size(),0.0);
    for(size_t ie = 0; ie < elems.size(); ++ie){
      auto const mat = elems[ie].kind == KKGeomCheck::Element::crvSector ? kkmat->CRVMaterial() : kkmat->material(elems[ie].material);
      if(!mat) throw cet::exception("KINKALGEOM") << "KinKalGeomCheck: no KinKal material '" << elems[ie].material << "'\n";
      rho[ie] = 1000.0*mat->density();
    }
    auto isCRV = [](std::string const& vol){ return vol.compare(0,3,"CRS") == 0; };
    // CRV mother volumes (CRSmother, CRSmotherLayer, CRSscintLayer) are air: count only the solid parts of the CRV
    auto isCRVSolid = [&isCRV](std::string const& vol, std::string const& mat){ return isCRV(vol) && mat != "G4_AIR"; };
    auto isCRVElement = [](KKGeomCheck::Element const& el){
      return el.kind == KKGeomCheck::Element::crvSector || el.name.compare(0,14,"CRV_StrongBack") == 0; };
    std::mt19937 rng(config_.seed());
    bool ok = true;

    // (a) material at points inside each element
    std::cout << "KinKalGeomCheck: G4 material at points inside each element (KinKal density, matching fraction, most common G4 materials)" << std::endl;
    unsigned nbadel = 0;
    for(size_t ie = 0; ie < elems.size(); ++ie){
      auto const& el = elems[ie];
      std::map<std::string,unsigned> counts;
      unsigned nmatch = 0;
      for(unsigned is = 0; is < config_.nMaterialSample(); ++is){
        std::string vol;
        auto const mat = g4Material(el.sample(rng),&vol);
        bool match = false;
        if(el.kind == KKGeomCheck::Element::crvSector) match = isCRVSolid(vol,mat);
        else if(mat == el.material) match = true;
        else { auto range = aliases_.equal_range(el.material); for(auto it = range.first; it != range.second; ++it) if(it->second == mat) match = true; }
        if(match) ++nmatch;
        ++counts[el.kind == KKGeomCheck::Element::crvSector && match ? std::string("CRS*") : mat];
      }
      double const frac = nmatch/double(config_.nMaterialSample());
      bool const good = frac >= config_.minMaterialFraction();
      if(!good) ++nbadel;
      std::vector<std::pair<unsigned,std::string>> top;
      for(auto const& [m,n] : counts) top.emplace_back(n,m);
      std::sort(top.rbegin(),top.rend());
      std::cout << (good ? "  " : "  BAD ") << std::setw(22) << std::left << el.name << std::right << " " << std::setw(16) << el.material
        << " rho " << std::setprecision(4) << rho[ie] << " match " << std::setprecision(3) << frac << " :";
      for(size_t it = 0; it < std::min<size_t>(3,top.size()); ++it) std::cout << " " << top[it].second << " " << top[it].first/double(config_.nMaterialSample());
      std::cout << std::endl;
    }
    if(nbadel > 0) ok = false;

    // (b) scintillator on each CRV layer's center plane: KinKal layer spans against G4
    if(kkg.CRV()){
      std::uniform_real_distribution<double> unit(-1.0,1.0);
      // a point is in KinKal scintillator if any sector has a layer there: neighbouring sectors' layers share the overlap bands
      auto inLayer = [&kkg](VEC3 const& p){
        for(auto const& sect : kkg.CRV()->sectors()){
          auto const& r = *sect.sector_;
          VEC3 const rel = p - r.center();
          if(std::fabs(rel.Dot(r.uDirection())) > r.uHalfLength()) continue;
          double const w = rel.Dot(r.normal()), v = rel.Dot(r.vDirection());
          for(auto const& layer : sect.layers_){
            if(std::fabs(w - layer.w_) > 1.0) continue;
            for(auto const& span : layer.v_) if(v >= span.first && v <= span.second) return true;
          }
        }
        return false;
      };
      std::cout << "KinKalGeomCheck: CRV layer spans against G4 on each layer center plane (fraction of points that agree)" << std::endl;
      unsigned nbadsect = 0;
      for(auto const& sector : kkg.CRV()->sectors()){
        auto const& rect = *sector.sector_;
        unsigned nagree = 0, ntot = 0;
        for(auto const& layer : sector.layers_){
          for(unsigned is = 0; is < config_.nMaterialSample(); ++is){
            VEC3 const p = rect.center() + 0.98*rect.uHalfLength()*unit(rng)*rect.uDirection() + rect.vHalfLength()*unit(rng)*rect.vDirection()
              + layer.w_*rect.normal();
            bool const kk = inLayer(p);
            if(kk == (g4Material(p) == config_.crvScintillator())) ++nagree;
            ++ntot;
          }
        }
        if(ntot == 0) continue;
        double const frac = nagree/double(ntot);
        bool const good = frac >= config_.minLayerAgreement();
        if(!good) ++nbadsect;
        std::cout << (good ? "  " : "  BAD ") << std::setw(8) << std::left << sector.sname_ << std::right << " layers " << sector.layers_.size()
          << " agree " << std::setprecision(4) << frac << std::endl;
      }
      if(nbadsect > 0) ok = false;
    }

    // (c) column density along lines from the tracker to the CRV sectors, outside the DS bore
    if(!kkg.CRV() || kkg.CRV()->sectors().empty() || !kkg.DS() || !kkg.tracker()) return ok;
    auto const& sectors = kkg.CRV()->sectors();
    auto const& trk = kkg.tracker()->outer();
    auto const& bore = kkg.DS()->inner();
    std::uniform_real_distribution<double> flat(-1.0,1.0);
    std::uniform_int_distribution<size_t> pick(0,sectors.size()-1);
    struct SectorSum { unsigned n = 0, nbad = 0; double g4p = 0, kkp = 0, g4c = 0, kkc = 0; };
    std::map<std::string,SectorSum> bysector;
    std::map<std::string,double> uncovered; // G4 column (g/cm^2) outside every KinKal element, per G4 volume
    unsigned nbad = 0;
    double const coverStep = config_.coverageStep();
    for(unsigned iray = 0; iray < config_.nRays(); ++iray){
      // start: uniform in the tracker volume; end: the far face of a random point of a random sector
      double const r = trk.radius()*std::sqrt(0.5*(1.0+flat(rng))), phi = M_PI*flat(rng);
      VEC3 const start = trk.center() + VEC3(r*std::cos(phi),r*std::sin(phi),trk.halfLength()*flat(rng));
      auto const& sector = sectors[pick(rng)];
      auto const& rect = *sector.sector_;
      VEC3 const end = rect.center() + (rect.uHalfLength()*flat(rng))*rect.uDirection() + (rect.vHalfLength()*flat(rng))*rect.vDirection()
        + sector.whw_*rect.normal();
      double const send = (end-start).R();
      VEC3 const dir = (end-start).Unit();
      // allowed segments: [0,send] outside the bore
      std::vector<std::pair<double,double>> allowed;
      double b0, b1;
      if(KKGeomCheck::cylinderInterval(bore.center(),bore.axis(),bore.radius(),bore.halfLength(),start,dir,0.0,send,b0,b1)){
        if(b0 > 0.0) allowed.emplace_back(0.0,b0);
        if(b1 < send) allowed.emplace_back(b1,send);
      } else allowed.emplace_back(0.0,send);
      // KinKal column
      double kkp = 0.0, kkc = 0.0;
      for(size_t ie = 0; ie < elems.size(); ++ie){
        double len = 0.0;
        for(auto const& seg : allowed) len += elems[ie].pathLength(start,dir,seg.first,seg.second);
        // a CRV sector adds the fraction of its layers this line crosses, as the extrapolation does
        if(len > 0.0 && elems[ie].sector >= 0){
          auto const& sect = sectors[elems[ie].sector];
          double const dn = dir.Dot(sect.sector_->normal());
          if(std::fabs(dn) > 1.0e-9) len *= sect.layerFraction(start + ((sect.sector_->center()-start).Dot(sect.sector_->normal())/dn)*dir,dir);
        }
        (isCRVElement(elems[ie]) ? kkc : kkp) += 0.1*len*rho[ie];
      }
      // G4 column: walk the GDML geometry along the line
      double g4p = 0.0, g4c = 0.0;
      auto const g0 = toGDML(start);
      gGeoManager->InitTrack(g0.X(),g0.Y(),g0.Z(),dir.X(),dir.Y(),dir.Z());
      double s = 0.0;
      for(unsigned istep = 0; istep < 100000 && s < send; ++istep){
        TGeoNode* node = gGeoManager->GetCurrentNode();
        if(!node || gGeoManager->IsOutside()) break;
        double const density = node->GetVolume()->GetMaterial()->GetDensity();
        std::string const vol = node->GetVolume()->GetName();
        gGeoManager->FindNextBoundaryAndStep(0.1*(send-s)+1.0e-4);
        double const step = std::max(10.0*gGeoManager->GetStep(),1.0e-6);
        if(density >= config_.minDensity()){
          for(auto const& seg : allowed){
            double const a = std::max(s,seg.first), b = std::min(std::min(s+step,send),seg.second);
            if(b <= a) continue;
            (isCRV(vol) ? g4c : g4p) += 0.1*(b-a)*density;
            // coverage: sample the step and tally the part outside every KinKal element
            unsigned const npt = std::max(1u,unsigned(std::ceil((b-a)/coverStep)));
            unsigned nout = 0;
            for(unsigned ip = 0; ip < npt; ++ip){
              VEC3 const p = start + (a + (ip+0.5)*(b-a)/npt)*dir;
              if(std::none_of(elems.begin(),elems.end(),[&p](auto const& el){ return el.inside(p); })) ++nout;
            }
            if(nout > 0) uncovered[vol.substr(0,vol.find("0x"))] += 0.1*(b-a)*density*nout/npt;
          }
        }
        s += step;
      }
      double const g4 = g4p + g4c, kk = kkp + kkc;
      bool const bad = std::fabs(kk-g4) > std::max(config_.rayAbsTolerance(),config_.rayRelTolerance()*g4);
      auto& sum = bysector[sector.sname_];
      ++sum.n; sum.g4p += g4p; sum.kkp += kkp; sum.g4c += g4c; sum.kkc += kkc;
      if(bad){
        ++nbad; ++sum.nbad;
        if(config_.debug() > 1){
          GeomHandle<DetectorSystem> det;
          auto const m = det->toMu2e(CLHEP::Hep3Vector(start.X(),start.Y(),start.Z()));
          std::cout << "  BAD ray to " << sector.sname_ << " start (" << m.x() << "," << m.y() << "," << m.z() << ") dir " << dir
            << " G4 " << g4p << "+" << g4c << " KinKal " << kkp << "+" << kkc << " g/cm^2" << std::endl;
        }
      }
    }
    std::cout << "KinKalGeomCheck: tracker-to-CRV column density outside the DS bore, mean g/cm^2 (passive + CRV)" << std::endl;
    for(auto const& [name,sum] : bysector){
      std::cout << "  " << std::setw(8) << std::left << name << std::right << " rays " << std::setw(5) << sum.n
        << " G4 " << std::setw(7) << sum.g4p/sum.n << " + " << std::setw(7) << sum.g4c/sum.n
        << "  KinKal " << std::setw(7) << sum.kkp/sum.n << " + " << std::setw(7) << sum.kkc/sum.n
        << "  bad " << sum.nbad/double(sum.n) << std::endl;
    }
    double const badfrac = nbad/double(std::max(1u,config_.nRays()));
    std::cout << "KinKalGeomCheck: " << nbad << " of " << config_.nRays() << " lines outside tolerance (fraction " << badfrac << ")" << std::endl;
    if(badfrac > config_.maxBadRayFraction()) ok = false;
    std::vector<std::pair<double,std::string>> unc;
    for(auto const& [vol,col] : uncovered) unc.emplace_back(col/config_.nRays(),vol);
    std::sort(unc.rbegin(),unc.rend());
    std::cout << "KinKalGeomCheck: G4 volumes crossed outside every KinKalGeom element (mean g/cm^2 per line)" << std::endl;
    for(size_t iu = 0; iu < std::min<size_t>(25,unc.size()); ++iu) std::cout << "  " << std::setw(40) << std::left << unc[iu].second << std::right << " " << unc[iu].first << std::endl;
    return ok;
  }
}

DEFINE_ART_MODULE(mu2e::KinKalGeomCheck)

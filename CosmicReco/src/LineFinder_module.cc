//
//
// Original author D. Brown and G. Tassielli
//

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA.hh"

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "TH2F.h"

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <utility>
#include <functional>
#include <float.h>
#include <vector>
#include <map>

namespace mu2e{

  class LineFinder : public art::EDProducer {
    public:
      struct Config{
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<int> diag{Name("diag"), Comment("Create diag histograms"),0};
        fhicl::Atom<int> minpeak{Name("minPeak"), Comment("Minimum hits in accumulator peak")};
        fhicl::Atom<float> maxDOCA{Name("maxDOCA"), Comment("Largest DOCA that is considered on the track (mm)")};
        fhicl::Atom<float> t0offset{Name("t0offset"), Comment("T0 offset")};
        fhicl::Atom<int> nsteps{Name("NSteps"), Comment("Number of steps per straw")};
        fhicl::Atom<int> ntsteps{Name("NTSteps"), Comment("Number of transverse steps per straw")};
        fhicl::Atom<float> stepsize{Name("StepSize"), Comment("Size of each step in fraction of res")};
        fhicl::Atom<unsigned> nmax{Name("MaxPairs"), Comment("Max pairs to try")};
        fhicl::Atom<art::InputTag> chToken{Name("ComboHitCollection"),Comment("tag for straw hit collection")};
        fhicl::Atom<art::InputTag> tcToken{Name("TimeClusterCollection"),Comment("tag for time cluster collection")};
        fhicl::Atom<art::InputTag> ccToken{Name("CaloClusterCollection"),Comment("tag for calo cluster collection"), ""};
        fhicl::Atom<bool> addCaloClusters{Name("AddCaloClusters"), Comment("Whether to connect calo clusters to the line seeds"), false};
      };
      typedef art::EDProducer::Table<Config> Parameters;
      explicit LineFinder(const Parameters& conf);
      virtual ~LineFinder(){};
      virtual void produce(art::Event& event ) override;
      virtual void beginRun(art::Run& run) override;

    private:

      Config _conf;

      //config parameters:
      int _diag;
      int  _minPeak;
      float _maxDOCA;
      float _t0offset;
      int _Nsteps, _Ntsteps;
      float _stepSize;
      unsigned _nmax;
      art::InputTag  _chToken;
      art::InputTag  _tcToken;
      art::InputTag  _ccToken;
      bool _addCaloClusters;

      ProditionsHandle<Tracker> _alignedTracker_h;
      const CaloClusterCollection* _cccol;
      art::Handle<CaloClusterCollection> _cccolH;
      const Calorimeter* _calorimeter;

      int findLine(const ComboHitCollection& shC, std::vector<StrawHitIndex> const& shiv, art::Event const& event, CosmicTrackSeed& tseed);
      void addClusterToSeed(CosmicTrackSeed& tseed);
};


 LineFinder::LineFinder(const Parameters& conf) :
   art::EDProducer(conf),
           _diag (conf().diag()),
        _minPeak (conf().minpeak()),
        _maxDOCA (conf().maxDOCA()),
        _t0offset (conf().t0offset()),
        _Nsteps (conf().nsteps()),
        _Ntsteps (conf().ntsteps()),
        _stepSize (conf().stepsize()),
        _nmax (conf().nmax()),
        _chToken (conf().chToken()),
        _tcToken (conf().tcToken()),
        _ccToken (conf().ccToken()),
        _addCaloClusters (conf().addCaloClusters()),
        _cccol(nullptr),
        _calorimeter(nullptr)
{
  consumes<ComboHitCollection>(_chToken);
  consumes<TimeClusterCollection>(_tcToken);
  if(_addCaloClusters) {
    consumes<CaloClusterCollection>(_ccToken);
  }
  produces<CosmicTrackSeedCollection>();

 }

void LineFinder::beginRun(art::Run& run) {
  _calorimeter = GeomHandle<Calorimeter>().get();
}

void LineFinder::produce(art::Event& event ) {

  auto const& chH = event.getValidHandle<ComboHitCollection>(_chToken);
  const ComboHitCollection& chcol(*chH);
  auto  const& tcH = event.getValidHandle<TimeClusterCollection>(_tcToken);
  const TimeClusterCollection& tccol(*tcH);

  if(_addCaloClusters) {
    event.getByLabel(_ccToken, _cccolH);
    if(!_cccolH.isValid()) {
      throw cet::exception("RECO") << "[LineFinder::produce] Could not retrieve calo cluster collection with tag " << _ccToken.encode();
    }
    _cccol = _cccolH.product();
  }


  std::unique_ptr<CosmicTrackSeedCollection> seed_col(new CosmicTrackSeedCollection());
  if(_diag > 2) std::cout << "[LineFinder:: " << __func__ << "] Event " << event.id()
                         << " has " << tccol.size() << " time clusters\n";

  for (size_t index=0;index< tccol.size();++index) {
    const auto& tclust = tccol[index];
    std::vector<StrawHitIndex> shiv;
    // get collection at straw level
    auto chcp = chcol.fillStrawHitIndices(tclust.hits(),shiv,StrawIdMask::uniquestraw);
    auto chc = *chcp;

    CosmicTrackSeed tseed;
    tseed._straw_chits.setAsSubset(chH,StrawIdMask::uniquestraw);

    tseed._timeCluster = art::Ptr<TimeCluster>(tcH,index);
    tseed._track.converged = true;

    if(_diag > 2) std::cout << "  Time cluster " << index << " has " << tclust.hits().size() << " hits"
                           << " and " << chc.size() << " straw hits\n";
    int seedSize = findLine(chc, shiv, event, tseed);
    if (_diag > 1)
      std::cout << "[LineFinder::" << __func__ <<"] Line seed N(hits) = " << seedSize << std::endl;

    if (seedSize >= _minPeak){
      if (_diag > 0)
        std::cout << "LineFinder: found line (N(hits) = " << seedSize << ")"
                  << " " << tseed._track.FitParams.T0
                  << " " << tseed._track.FitParams.A0
                  << " " << tseed._track.FitParams.B0
                  << " " << tseed._track.FitParams.A1
                  << " " << tseed._track.FitParams.B1
                  << std::endl;
      if(_addCaloClusters) addClusterToSeed(tseed);
      tseed._status.merge(TrkFitFlag::Straight);
      tseed._status.merge(TrkFitFlag::hitsOK);
      tseed._status.merge(TrkFitFlag::helixOK);
      tseed._status.merge(TrkFitFlag::helixConverged);
      tseed._track.MinuitParams.cov = std::vector<double>(15, 0);
      seed_col->push_back(tseed);
    }
  }

  event.put(std::move(seed_col));
}

int LineFinder::findLine(const ComboHitCollection& shC, std::vector<StrawHitIndex> const& shiv, art::Event const& event, CosmicTrackSeed& tseed){

  //mu2e::GeomHandle<mu2e::Tracker> th;
  auto tracker = _alignedTracker_h.getPtr(event.id());//th.get();
  int bestcount = 0;
  double bestll = 0;

  CLHEP::Hep3Vector seedDir(0,0,0);
  CLHEP::Hep3Vector seedInt(0,0,0);

  bool found_all = false;
  unsigned n = 0;
  // lets get the best pairwise vector
  for (size_t i=0;i<shiv.size();i++){
    size_t iloc = shiv[i];
    if (found_all)
      break;
    Straw const& strawi = tracker->getStraw(shC[iloc].strawId());
    for (size_t j=shiv.size()-1;j>i;j--){
      size_t jloc = shiv[j];
      if (found_all)
        break;
      n += 1;
      Straw const& strawj = tracker->getStraw(shC[jloc].strawId());
      for (int is=-1*_Nsteps;is<_Nsteps+1;is++){
        CLHEP::Hep3Vector ipos = shC[iloc].posCLHEP() + strawi.getDirection()*shC[iloc].wireRes()*_stepSize*is;
        for (int js=-1*_Nsteps;js<_Nsteps+1;js++){
          CLHEP::Hep3Vector jpos = shC[jloc].posCLHEP() + strawj.getDirection()*shC[jloc].wireRes()*_stepSize*js;

          CLHEP::Hep3Vector newdir = (jpos-ipos).unit();
          CLHEP::Hep3Vector icross = (jpos-ipos).cross(strawi.getDirection()).unit();
          CLHEP::Hep3Vector jcross = (jpos-ipos).cross(strawj.getDirection()).unit();
          for (int its=-1*_Ntsteps;its<_Ntsteps+1;its++){
            ipos += icross*2.5*its;
            for (int jts=-1*_Ntsteps;jts<_Ntsteps+1;jts++){
              jpos += jcross*2.5*its;

              // now loop over all hits and see how many are in this track
              int count = 0;
              double ll = 0;
              for (size_t k=0;k<shiv.size();k++){
                size_t kloc = shiv[k];
                Straw const& strawk = tracker->getStraw(shC[kloc].strawId());
                TwoLinePCA pca( strawk.getMidPoint(), strawk.getDirection(),
                    ipos, newdir);
                double dist = (pca.point1()-strawk.getMidPoint()).dot(strawk.getDirection());
                if (pca.dca() < _maxDOCA && fabs(dist) < strawk.halfLength()){
                  count += 1;
                  ll += pow(dist-shC[kloc].wireDist(),2)/shC[kloc].wireVar();
                }
              }
              if (count > bestcount || (count == bestcount && ll < bestll)){
  //              besti = i;
  //              bestj = j;
                bestcount = count;
                bestll = ll;
                seedDir = newdir.unit();
                if (seedDir.y() > 0) seedDir *= -1;
                if (seedDir.y() != 0)
                  seedInt = ipos - newdir*ipos.y()/newdir.y();
              }
            }
          }
        }
      }
      if (bestcount == (int) shiv.size())
        found_all = true;
      if (n > _nmax){
        i = shiv.size();j = 1;
      }
    }
  }

  // get pos and direction into Z alignment
  if (seedDir.y() != 0){
    seedDir /= -1*seedDir.y();
    seedInt -= seedDir*seedInt.y()/seedDir.y();
  }

  double avg_t0 = 0;
  int good_hits = 0;
  for (size_t k=0;k<shiv.size();k++){
    size_t kloc = shiv[k];
    Straw const& strawk = tracker->getStraw(shC[kloc].strawId());
    TwoLinePCA pca( strawk.getMidPoint(), strawk.getDirection(),
        seedInt, seedDir);
    double dist = (pca.point1()-strawk.getMidPoint()).dot(strawk.getDirection());
    if (pca.dca() < _maxDOCA && fabs(dist) < strawk.halfLength()){
      double traj_time = ((pca.point2() - seedInt).dot(seedDir.unit()))/299.9;
      double hit_t0 = shC[kloc].time() - shC[kloc].driftTime() - shC[kloc].propTime() - traj_time;
      avg_t0 += hit_t0;
      good_hits++;
      ComboHit combohit;
      combohit.init(shC[kloc],kloc);
      tseed._straw_chits.push_back(std::move(combohit));
    }
  }
  avg_t0 /= good_hits;

  tseed._t0._t0 = avg_t0-_t0offset;

  tseed._track.FitParams.T0 = tseed._t0._t0;
  tseed._track.FitParams.A0 = seedInt.x();
  tseed._track.FitParams.B0 = seedInt.z();
  tseed._track.FitParams.A1 = seedDir.x();
  tseed._track.FitParams.B1 = seedDir.z();
  tseed._track.MinuitParams.T0 = tseed._t0._t0;
  tseed._track.MinuitParams.A0 = seedInt.x();
  tseed._track.MinuitParams.B0 = seedInt.z();
  tseed._track.MinuitParams.A1 = seedDir.x();
  tseed._track.MinuitParams.B1 = seedDir.z();
  XYZVectorF X(1,0,0);
  XYZVectorF Y(0,1,0);
  XYZVectorF Z(0,0,1);
  TrackAxes XYZ(X,Y,Z);
  tseed._track.InitCoordSystem = XYZ;
  tseed._track.FitCoordSystem = XYZ;
  XYZVectorF xyzint(seedInt);
  XYZVectorF xyzdir(seedDir);
  TrackEquation XYZTrack(xyzint,xyzdir);
  tseed._track.SetFitEquation(XYZTrack);
  tseed._track.SetMinuitEquation(XYZTrack);

  return good_hits;
}

void LineFinder::addClusterToSeed(CosmicTrackSeed& tseed) {
  if(!_cccol) return;

  // Get the track seed information
  const XYZVectorF seedInt = tseed._track.Pos0();
  const XYZVectorF seedDir = tseed._track.Dir();
  if(std::abs(seedDir.z()) < 1e-3) return; // avoid near-vertical tracks
  const CLHEP::Hep3Vector seedDirCLHEP(seedDir.x(), seedDir.y(), seedDir.z());
  const CLHEP::Hep3Vector seedIntCLHEP(seedInt.x(), seedInt.y(), seedInt.z());
  const double t0 = tseed._t0._t0;

  // Check if a calo cluster can be matched to the line seed
  const CaloCluster* matchedCluster = nullptr;
  int matchedClusterIndex = -1;
  double min_t0_diff = 1.e10;
  constexpr double t0_match_window = 50.; // ns, loose time window for matching
  for(size_t i_cl = 0; i_cl < _cccol->size(); ++i_cl) {
    const CaloCluster& cl = _cccol->at(i_cl);
    const CLHEP::Hep3Vector cl_pos = _calorimeter->geomUtil().mu2eToTracker(_calorimeter->geomUtil().diskToMu2e(cl.diskID(), cl.cog3Vector()));
    const double cl_time = cl.time();
    // Get the expected time at the calo disk based on the line seed
    const double speed = CLHEP::c_light; // assume speed of the particle is close to the speed of light
    const double dr = std::abs((cl_pos - seedIntCLHEP).dot(seedDirCLHEP.unit()));
    const double t_flight = dr / speed;
    // check either trajectory
    const double t0_diff_1 = std::abs((t0 + t_flight) - cl_time);
    const double t0_diff_2 = std::abs((t0 - t_flight) - cl_time);
    const double t0_diff = std::min(t0_diff_1, t0_diff_2);
    if(t0_diff < min_t0_diff && t0_diff < t0_match_window) {
      min_t0_diff = t0_diff;
      matchedCluster = &cl;
      matchedClusterIndex = i_cl;
    }
  }

  if(matchedCluster) {
    tseed._caloCluster = art::Ptr<CaloCluster>(_cccolH, matchedClusterIndex);
    if(_diag > 0) {
      std::cout << "LineFinder: matched calo cluster with energy " << matchedCluster->energyDep() << " MeV"
                << " time " << matchedCluster->time() << " ns to line seed with t0 = "
                << tseed._t0._t0 << " ns" << std::endl;
    }
  }
}

}//end mu2e namespace
using mu2e::LineFinder;
DEFINE_ART_MODULE(LineFinder)

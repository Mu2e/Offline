//
// TODO: explain purpose
// Original author S. Middleton and H. Jafree
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
#include "Offline/TrackerGeom/inc/Tracker.hh"

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh" 

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

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

  class PhotonConversionFinder : public art::EDProducer {
    public:
      struct Config{
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<int> diag{Name("diag"), Comment("Create diag histograms"),0};
        fhicl::Atom<art::InputTag> chToken{Name("ComboHitCollection"),Comment("tag for straw hit collection")};
        fhicl::Atom<art::InputTag> tcToken{Name("TimeClusterCollection"),Comment("tag for time cluster collection")};
      };
      typedef art::EDProducer::Table<Config> Parameters;
      explicit PhotonConversionFinder(const Parameters& conf);
      virtual ~PhotonConversionFinder(){};
      virtual void beginRun(art::Run& ) override;
      virtual void beginJob() override;
      virtual void produce(art::Event& event ) override;

    private:

      Config _conf;

      //config parameters:
      int _diag;
      art::InputTag  _chToken;
      art::InputTag  _tcToken;

  };


 PhotonConversionFinder::PhotonConversionFinder(const Parameters& conf) :
  art::EDProducer(conf),
  _diag (conf().diag()),
  _chToken (conf().chToken()),
  _tcToken (conf().tcToken())
  {
    consumes<ComboHitCollection>(_chToken);
    consumes<TimeClusterCollection>(_tcToken);
    produces<HelixSeedCollection>();
  }
 
  void PhotonConversionFinder::beginRun(art::Run& )  {//TODO 
  }

  void PhotonConversionFinder::beginJob() { //TODO - can add TTree and THistograms here if required
  }

  void PhotonConversionFinder::produce(art::Event& event ) {
    /*auto const& chH = event.getValidHandle<ComboHitCollection>(_chToken);
    const ComboHitCollection& chcol(*chH);
    auto  const& tcH = event.getValidHandle<TimeClusterCollection>(_tcToken);
    const TimeClusterCollection& tccol(*tcH);*/
    std::unique_ptr<HelixSeedCollection> seed_col(new HelixSeedCollection());
    //TODO - this is where your algorithm will be...
    event.put(std::move(seed_col)); //TODO - one per event
  }


}//end mu2e namespace
using mu2e::PhotonConversionFinder;
DEFINE_ART_MODULE(PhotonConversionFinder);

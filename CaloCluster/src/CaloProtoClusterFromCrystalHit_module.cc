//
// The clustering is performed in several stages, using the folowing data structures:
//  - bidimensional array that have a list of energy deposits for each crystal id (prefilter low energy deposits to speed up things)
//  - a map containing the potential seeds (the map is ordered by most energetic to lowest energetic hits)

// The clustering proceeds in three steps: form energetic proto-clusters, look at split-offs and form final clusters. The first two steps are done by CaloProtoClusterFromCrystalHit
// and produce proto-clusters. The last step is done by MakeCaloCluster, and procudes a cluster

// 1. Energetic proto-clusters (clusters above some threshold energy)
//    - start from the most energetic seed (over some threshold)
//    - form a proto-cluster by adding all simply connected cluster to the seed (simply connected = any two hits in a cluster can be joined
//      by a continuous path of clusters in the crystal)
//    - mark the correpsonding hits as used, and update the seed list
//
// 2. Split-offs: some clusters might have low-energy split-offs, and we need to find them
//    - filter the remaining unassigned hits to retain only those compatible with the time of the main clusters (there are a lot of background low energy deposits
//      incompatible with any clsuter, we don't want them)
//    - update the seed lists and form all remaining clusters as before
//
// 3. Cluster formation
//    - merge proto-clusters and split-offs to form clusters ifg they are "close" enough
//
// Note: 1. One can try to filter the hits first, before producing energetic proto-clusters (use two iterators on the time ordered crystal list, see commented code at the end).
//          The performance difference is very small compared to this implementation, but more obscure, so for clarity, I kept this one.
//       2. I tried a bunch of other optimizations but the performance increase is small and the code harder to read, so I left it as is
//
// Original author: B. Echenard
//


// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "cetlib_except/exception.h"

#include "CaloCluster/inc/ClusterFinder.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"

// Other includes.
#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <memory>




namespace mu2e {


  class CaloProtoClusterFromCrystalHit : public art::EDProducer {


  public:

    typedef std::vector<const CaloCrystalHit*>  CaloCrystalVec;
    typedef std::list<const CaloCrystalHit*>    CaloCrystalList;


    explicit CaloProtoClusterFromCrystalHit(fhicl::ParameterSet const& pset) :
      art::EDProducer{pset},
      caloCrystalToken_{consumes<CaloCrystalHitCollection>(pset.get<std::string>("caloCrystalModuleLabel"))},
      producerNameMain_(pset.get<std::string>("mainClusterCollName")),
      producerNameSplit_(pset.get<std::string>("splitClusterCollName")),
      EminSeed_(pset.get<double>("EminSeed")),
      EnoiseCut_(pset.get<double>("EnoiseCut")),
      ExpandCut_(pset.get<double>("ExpandCut")),
      timeCut_(pset.get<double>("timeCut")),
      deltaTime_(pset.get<double>("deltaTime")),
      diagLevel_(pset.get<int>("diagLevel",0)),
      messageCategory_("CLUSTER")
    {
      produces<CaloProtoClusterCollection>(producerNameMain_);
      produces<CaloProtoClusterCollection>(producerNameSplit_);
    }

    void produce(art::Event& e) override;

  private:

    art::ProductToken<CaloCrystalHitCollection> const caloCrystalToken_;
    std::string       producerNameMain_;
    std::string       producerNameSplit_;
    double            EminSeed_;
    double            EnoiseCut_;
    double            ExpandCut_;
    double            timeCut_;
    double            deltaTime_;
    int               diagLevel_;
    const std::string messageCategory_;

    void makeProtoClusters(CaloProtoClusterCollection& caloProtoClustersMain,
                           CaloProtoClusterCollection& caloProtoClustersSplit,
                           const art::Handle<CaloCrystalHitCollection>& CaloCrystalHitsHandle);

    void filterByTime(CaloCrystalList& liste,
                      const std::vector<double>& clusterTime,
                      std::list<const CaloCrystalHit*> &seedList);

    void fillCluster(CaloProtoClusterCollection& caloProtoClustersColl, const CaloCrystalList& clusterList,
                     const art::Handle<CaloCrystalHitCollection>& CaloCrystalHitsHandle);

    void dump(const std::vector<CaloCrystalList>& caloIdHitMap);

  };


  void CaloProtoClusterFromCrystalHit::produce(art::Event& event)
  {
    // Check that calorimeter geometry description exists
    art::ServiceHandle<GeometryService> geom;
    if( !(geom->hasElement<Calorimeter>()) ) return;

    // Get handles to calorimeter crystal hits
    art::Handle<CaloCrystalHitCollection> CaloCrystalHitsHandle;
    bool const success = event.getByToken(caloCrystalToken_, CaloCrystalHitsHandle);
    if (!success) return;

    // Create a new CaloCluster collection and fill it
    auto caloProtoClustersMain = std::make_unique<CaloProtoClusterCollection>();
    auto caloProtoClustersSplit = std::make_unique<CaloProtoClusterCollection>();
    makeProtoClusters(*caloProtoClustersMain,*caloProtoClustersSplit,CaloCrystalHitsHandle);

    event.put(std::move(caloProtoClustersMain),  producerNameMain_);
    event.put(std::move(caloProtoClustersSplit), producerNameSplit_);
  }


  //----------------------------------------------------------------------------------------------------------
  void CaloProtoClusterFromCrystalHit::makeProtoClusters(CaloProtoClusterCollection& caloProtoClustersMain,
                                                         CaloProtoClusterCollection& caloProtoClustersSplit,
                                                         const art::Handle<CaloCrystalHitCollection> & CaloCrystalHitsHandle)
  {
    const Calorimeter& cal = *(GeomHandle<Calorimeter>());
    const CaloCrystalHitCollection& CaloCrystalHits(*CaloCrystalHitsHandle);
    if (CaloCrystalHits.empty()) return;


    //declare and fill the hash map crystal_id -> list of CaloHits
    std::vector<CaloCrystalList>      mainClusterList, splitClusterList,caloIdHitMap(cal.nCrystal());
    std::list<const CaloCrystalHit*>  seedList;
    std::vector<double>               clusterTime;



    //fill data structures
    for (const auto& hit : CaloCrystalHits)
      {
        if (hit.energyDep() < EnoiseCut_ || hit.time() < timeCut_) continue;
        caloIdHitMap[hit.id()].push_back(&hit);
        seedList.push_back(&hit);
      }

    seedList.sort([](const CaloCrystalHit* a, const CaloCrystalHit* b) {return a->energyDep() > b->energyDep();});


    //produce main clusters
    while( !seedList.empty() )
      {
        const CaloCrystalHit* crystalSeed = *seedList.begin();
        if (crystalSeed->energyDep() < EminSeed_) break;

        ClusterFinder finder(cal,crystalSeed,deltaTime_, ExpandCut_);
        finder.formCluster(caloIdHitMap);

        mainClusterList.push_back(finder.clusterList());
        clusterTime.push_back(crystalSeed->time());

        for (const auto& hit: finder.clusterList()) seedList.remove(hit);
      }


    //filter unneeded hits
    for (unsigned i=0; i < caloIdHitMap.size(); ++i ) filterByTime(caloIdHitMap[i], clusterTime, seedList);


    //produce split-offs clusters
    while(!seedList.empty())
      {
        const CaloCrystalHit* crystalSeed = *seedList.begin();
        ClusterFinder finder(cal,crystalSeed,deltaTime_, ExpandCut_);

        finder.formCluster(caloIdHitMap);
        splitClusterList.push_back(finder.clusterList());

        for (const auto& hit: finder.clusterList()) seedList.remove(hit);
      }




    //save the main and split clusters
    for (auto cluster : mainClusterList)  fillCluster(caloProtoClustersMain,cluster,CaloCrystalHitsHandle);
    for (auto cluster : splitClusterList) fillCluster(caloProtoClustersSplit,cluster,CaloCrystalHitsHandle);



    //sort these guys
    std::sort(caloProtoClustersMain.begin(),  caloProtoClustersMain.end(), [](const CaloProtoCluster& a, const CaloProtoCluster& b) {return a.time() < b.time();});
    std::sort(caloProtoClustersSplit.begin(), caloProtoClustersSplit.end(),[](const CaloProtoCluster& a, const CaloProtoCluster& b) {return a.time() < b.time();});

  }





  //----------------------------------------------------------------------------------------------------------
  void CaloProtoClusterFromCrystalHit::fillCluster(CaloProtoClusterCollection& caloProtoClustersColl, const CaloCrystalList& clusterPtrList,
                                                   const art::Handle<CaloCrystalHitCollection>& CaloCrystalHitsHandle)
  {

    const CaloCrystalHitCollection& CaloCrystalHits(*CaloCrystalHitsHandle);
    const CaloCrystalHit* caloCrystalHitBase = &CaloCrystalHits.front();

    std::vector<art::Ptr<CaloCrystalHit>> caloCrystalHitsPtrVector;
    double totalEnergy(0),totalEnergyErr(0);
    //double timeW(0),timeWtot(0);

    for (auto clusterPrt : clusterPtrList)
      {
        //double weight = 1.0/clusterPrt->timeErr()/clusterPrt->timeErr();
        //timeW    += weight*clusterPrt->time();
        //timeWtot += weight;

        totalEnergy    += clusterPrt->energyDep();
        totalEnergyErr += clusterPrt->energyDepErr()*clusterPrt->energyDepErr();

        size_t idx = (clusterPrt - caloCrystalHitBase);
        caloCrystalHitsPtrVector.push_back( art::Ptr<CaloCrystalHit>(CaloCrystalHitsHandle,idx) );
      }

    totalEnergyErr = sqrt(totalEnergyErr);
    double time    = (*clusterPtrList.begin())->time();
    double timeErr = (*clusterPtrList.begin())->timeErr();
    //double time    = timeW/timeWtot;
    //double timeErr = 1.0/sqrt(timeWtot);


    caloProtoClustersColl.emplace_back(CaloProtoCluster(time,timeErr,totalEnergy,totalEnergyErr,caloCrystalHitsPtrVector,false));


    if (diagLevel_ > 1)
      {
        std::cout<<"This cluster contains "<<clusterPtrList.size()<<" crystals, id= ";
        for (auto clusterPrt : clusterPtrList) std::cout<<clusterPrt->id()<<" ";
        std::cout<<" with energy="<<totalEnergy<<" and time="<<time<<std::endl;;
      }
  }




  //----------------------------------------------------------------------------------------------------------
  void CaloProtoClusterFromCrystalHit::filterByTime(CaloCrystalList& liste, const std::vector<double>& clusterTime, std::list<const CaloCrystalHit*> &seedList)
  {
    for (auto it = liste.begin(); it != liste.end(); ++it)
      {
        const CaloCrystalHit* hit = *it;

        auto itTime = clusterTime.begin();
        while (itTime != clusterTime.end())
          {
            if ( (*itTime - hit->time()) < deltaTime_) break;
            ++itTime;
          }

        if (itTime == clusterTime.end() ) {seedList.remove(hit); liste.erase(it); if (it != liste.begin()) --it;}
      }
  }



  //----------------------------------------------------------------------------------------------------------
  void CaloProtoClusterFromCrystalHit::dump(const std::vector<CaloCrystalList>& caloIdHitMap)
  {
    for (unsigned int i=0; i<caloIdHitMap.size(); ++i)
      {
        if (caloIdHitMap[i].size()>0) std::cout<<"ProtoCluster crystal "<<i<<" has size "<<caloIdHitMap[i].size()<<std::endl;
        for (auto j :  caloIdHitMap[i]) std::cout<<i<<" "<<j->id()<<" "<<j->energyDep()<<" "<<j->time()<<std::endl;
      }
  }


}

DEFINE_ART_MODULE(mu2e::CaloProtoClusterFromCrystalHit);


/*

// Just in case
// This is a snippet of code to include only the hits compatible with the time of potential seeds in the map.
// The gain in performance is small, and I find it obscures the code, so I left the old version


//fast forward iterator until first crystal in time
std::vector<CaloCrystalHit>::const_iterator allCrystal = CaloCrystalHits.begin();
while (allCrystal->time() < timeCut_ && allCrystal != CaloCrystalHits.end()) ++allCrystal;
if (allCrystal == CaloCrystalHits.end()) return;


//fast forward iterator until first high energy crystal
std::vector<CaloCrystalHit>::const_iterator highCrystal = allCrystal;
while (highCrystal->energyDep() < EminSeed_ && highCrystal != CaloCrystalHits.end()) ++highCrystal;
if (highCrystal == CaloCrystalHits.end() ) return;


double maxTime = highCrystal->time() + deltaTime_Plus;
double minTime = highCrystal->time() - deltaTime_Minus;

while( allCrystal != CaloCrystalHits.end() )
{
if (allCrystal->time() > maxTime && highCrystal != CaloCrystalHits.end())
{
do ++highCrystal; while ( highCrystal != CaloCrystalHits.end() && highCrystal->energyDep() < EminSeed_);
if (highCrystal == CaloCrystalHits.end()) --highCrystal;
maxTime = highCrystal->time() + deltaTime_Plus;
minTime = highCrystal->time() - deltaTime_Minus;
continue;
}

if (allCrystal->energyDep() > EnoiseCut_ && allCrystal->time() > minTime )
{
const CaloCrystalHit* hit = &(*allCrystal);
caloIdHitMap[hit->id()].push_back(hit);
seedList.add(hit);
}
++allCrystal;
}


//--------------------------------------------------------------------------------------------------------------------
namespace
{
struct caloSeedCompare {
bool operator() (mu2e::const CaloCrystalHit* a, mu2e::const CaloCrystalHit* b) const
{
if (std::abs(a->energyDep() - b->energyDep()) > 1e-6) return a->energyDep() > b->energyDep();
if (a->id() != b->id()) return a->id() > b->id();
return a->time() > b->time();
}
};
}
*/

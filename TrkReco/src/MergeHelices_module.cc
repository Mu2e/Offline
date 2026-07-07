//
//  Merge HelixSeeds from different pattern recognition paths.  Modeled
//  on MergeHelixFinder.  P. Murat, D. Brown (LBNL)
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
// mu2e data products
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
// utilities
#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums4.hh"
// C++
#include <vector>
#include <memory>
#include <iostream>
#include <forward_list>
#include <string>

namespace mu2e {
  class MergeHelices : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<int>             debug        { Name("debugLevel"),            Comment("Debug Level"), 0};
        fhicl::Atom<unsigned>        deltanh      { Name("deltanh"),               Comment("difference in the active StrawHit counts")};
        fhicl::Atom<float>           scaleXY      { Name("scaleXY"),               Comment("scaling factor to get chi2XY/ndof distribution peak at 1")};
        fhicl::Atom<float>           scaleZPhi    { Name("scaleZPhi"),             Comment("scaling factor to get chi2ZPhi/ndof distribution peak at 1")};
        fhicl::Atom<bool>            selectbest   { Name("SelectBest"),            Comment("Select best overlapping helices for output"), true};
        fhicl::Atom<bool>            usecalo      { Name("UseCalo"),               Comment("Use CaloCluster info in comparison"), true};
        fhicl::Atom<unsigned>        minnover     { Name("MinNHitOverlap"),        Comment("Minimum number of common hits to consider helices to be 'the same'")};
        fhicl::Atom<float>           minoverfrac  { Name("MinHitOverlapFraction"), Comment("Minimum fraction of common hits to consider helices to be 'the same'")};
        fhicl::Sequence<std::string> BadHitFlags  { Name("BadHitFlags"),           Comment("HelixHit flag bits to exclude from counting"),std::vector<std::string>{"Outlier"}};
        fhicl::Sequence<std::string> HelixFinders { Name("HelixFinders"),          Comment("HelixSeed producers to merge")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit MergeHelices(const Parameters& conf);
      void produce(art::Event& evt) override;
    private:
      enum HelixComp{unique=-1,first=0,second=1};
      int _debug;
      unsigned _deltanh;
      float _scaleXY,_scaleZPhi;
      unsigned _minnover;
      float _minoverfrac;
      bool _selectbest, _usecalo;
      std::vector<std::string> _hfs;
      StrawHitFlag _badhit;
      // helper functions
      typedef std::vector<StrawHitIndex> SHIV;
      HelixComp compareHelices(art::Event const& evt, HelixSeed const& h1, HelixSeed const& h2);
      void countHits(art::Event const& evt, HelixSeed const& h1, HelixSeed const& h2, unsigned& nh1, unsigned& nh2, unsigned& nover);
      unsigned countOverlaps(SHIV const& s1, SHIV const& s2);
      void findchisq(art::Event const& evt, HelixSeed const& h2, float& chixy, float& chizphi) const;
  };

  MergeHelices::MergeHelices(const Parameters& config) : art::EDProducer{config},
    _debug(config().debug()),
    _deltanh(config().deltanh()),
    _scaleXY(config().scaleXY()),
    _scaleZPhi(config().scaleZPhi()),
    _minnover(config().minnover()),
    _minoverfrac(config().minoverfrac()),
    _selectbest(config().selectbest()),
    _usecalo(config().usecalo()),
    _hfs(config().HelixFinders()),
    _badhit(config().BadHitFlags())
    {
      consumesMany<HelixSeedCollection>    ();
      produces<HelixSeedCollection>    ();
      produces<TimeClusterCollection>    ();
    }

  void MergeHelices::produce(art::Event& event) {
    // create output
    std::unique_ptr<HelixSeedCollection> mhels(new HelixSeedCollection);
    std::unique_ptr<TimeClusterCollection> tcs(new TimeClusterCollection);
    // needed for creating Ptrs
    auto TimeClusterCollectionPID = event.getProductID<TimeClusterCollection>();
    auto TimeClusterCollectionGetter = event.productGetter(TimeClusterCollectionPID);
    // loop over helix products and flatten the helix collections into a single collection
    std::list<const HelixSeed*> hseeds;
    for(auto const& hf : _hfs) {
      art::InputTag hsct(hf);
      auto hsch = event.getValidHandle<HelixSeedCollection>(hsct);
      auto const& hsc = *hsch;
      for(auto const&  hs : hsc) {
        hseeds.insert(hseeds.end(),&hs);
      }
    }
    // now loop over all combinations
    for(auto ihel = hseeds.begin(); ihel != hseeds.end();) {
      auto jhel = ihel; jhel++;
      while( jhel != hseeds.end()) {
        // compare the helices
        auto hcomp = compareHelices(event, **ihel, **jhel);
        if(hcomp == unique)
          jhel++; // both helices are unique: simply advance the iterator to keep both
        else if(hcomp == first)
          jhel = hseeds.erase(jhel); // first helix is 'better'; remove the second
        else if(hcomp == second){ // second helix is better; remove the first and restart the loop
          ihel = hseeds.erase(ihel);
          break;
        }
      }
      // only advance the outer loop if we exhausted the inner one
      if(jhel == hseeds.end())
        ihel++;
    }
    // now hseeds contains only pointers to unique and best helices: use them to
    // create the output, which must be deep-copy, including the time cluster
    for(auto ihel = hseeds.begin(); ihel != hseeds.end(); ihel++){
      // copy the time cluster first
      tcs->push_back(*(*ihel)->timeCluster());
      // create a Ptr to the new TimeCluster
      auto tcptr = art::Ptr<TimeCluster>(TimeClusterCollectionPID,tcs->size()-1,TimeClusterCollectionGetter);
      // copy the Helix Seed and point its TimeCluster to the copy
      HelixSeed hs(**ihel);
      hs._timeCluster = tcptr;
      float chixy(0),chizphi(0);
      // Calculate the XY and ZPhi chisquared of the helix
      findchisq(event,**ihel,chixy,chizphi);
      // Replace the helix seed chisquared with the new values
      hs._helix._chi2dXY = chixy;
      hs._helix._chi2dZPhi = chizphi;
      hs._timeCluster = tcptr;
      mhels->push_back(hs);
    }
    event.put(std::move(mhels));
    event.put(std::move(tcs));
  }

  MergeHelices::HelixComp MergeHelices::compareHelices(art::Event const& evt, HelixSeed const& h1, HelixSeed const& h2) {
    HelixComp retval(unique);
    unsigned nh1, nh2, nover;
    // count the StrawHit overlap between the helices
    countHits(evt,h1,h2, nh1, nh2, nover);
    unsigned minh = std::min(nh1, nh2);
    float chih1xy(0),chih1zphi(0),chih2xy(0),chih2zphi(0);
    findchisq(evt,h1,chih1xy,chih1zphi);
    // Calculate the chi-sq of the helices
    findchisq(evt,h2,chih2xy,chih2zphi);
    // overlapping helices: decide which is best
    if(nover >= _minnover && nover/float(minh) > _minoverfrac) {
      if(h1.caloCluster().isNonnull() && h2.caloCluster().isNull())
        retval = first;
      // Pick the one with a CaloCluster first
      else if( h2.caloCluster().isNonnull() && h1.caloCluster().isNull())
        retval = second;
      // then compare active StrawHit counts and if difference of the StrawHit counts greater than deltanh
      else if((nh1 > nh2) && (nh1-nh2) > _deltanh)
        retval = first;
      else if((nh2 > nh1) && (nh2-nh1) > _deltanh)
        retval = second;
      // finally compare chisquared: sum xy and fz
      else if(chih1xy+chih1zphi  < chih2xy+chih2zphi)
        retval = first;
      else
        retval = second;
    }
    return retval;
  }

  void MergeHelices::countHits(art::Event const& evt, HelixSeed const& h1, HelixSeed const& h2, unsigned& nh1, unsigned& nh2, unsigned& nover) {
    nh1 = nh2 = nover = 0;
    SHIV shiv1, shiv2;
    for(size_t ihit=0;ihit < h1.hits().size(); ihit++){
      auto const& hh = h1.hits()[ihit];
      if(!hh.flag().hasAnyProperty(_badhit))
        h1.hits().fillStrawHitIndices(ihit,shiv1);
    }
    for(size_t ihit=0;ihit < h2.hits().size(); ihit++){
      auto const& hh = h2.hits()[ihit];
      if(!hh.flag().hasAnyProperty(_badhit))
        h2.hits().fillStrawHitIndices(ihit,shiv2);
    }
    nh1 = shiv1.size();
    nh2 = shiv2.size();
    nover = countOverlaps(shiv1,shiv2);
  }

  void MergeHelices::findchisq(art::Event const& evt, HelixSeed const& h1,float& chixy, float& chizphi) const{
    ::LsqSums4 sxy;
    ::LsqSums4 szphi;
    for(auto const& hh : h1.hits()) {
      if(!hh.flag().hasAnyProperty(_badhit)){
        // Calculate the transverse and z-phi weights of the hits
        float dx = hh.pos().x() - h1.helix().center().x();
        float dy = hh.pos().y() - h1.helix().center().y();
        float dxn = dx*hh.vDir().x()+dy*hh.vDir().y();
        float costh2 = dxn*dxn/(dx*dx+dy*dy);
        float sinth2 = 1-costh2;
        float e2xy = hh.uVar()*sinth2+hh.vVar()*costh2;
        float wtxy = 1./e2xy;
        float e2zphi = hh.uVar()*costh2+hh.vVar()*sinth2;
        float wtzphi = h1.helix().radius()*h1.helix().radius()/e2zphi;
        wtxy *= _scaleXY;
        wtzphi *= _scaleZPhi;
        sxy.addPoint(hh.pos().x(),hh.pos().y(),wtxy);
        szphi.addPoint(hh.pos().z(),hh.helixPhi(),wtzphi);
      }
    }
    chixy = sxy.chi2DofCircle();
    chizphi = szphi.chi2DofLine();
  }

  unsigned MergeHelices::countOverlaps(SHIV const& shiv1, SHIV const& shiv2) {
    unsigned retval(0);
    for(auto shi1 : shiv1)
      for(auto shi2 : shiv2)
        if(shi1 == shi2)
          retval++;
    return retval;
  }
}
DEFINE_ART_MODULE(mu2e::MergeHelices)

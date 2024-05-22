// Containerize StrawDigiBundles, with a extra functionality for resolving
// situations where multiple digis have overlapping digitization windows
//
// Ed Callaghan, 2023

#ifndef TrackerMC_StrawDigiBundleCollection_hh
#define TrackerMC_StrawDigiBundleCollection_hh

// stl
#include <map>
#include <memory>
#include <vector>

// art
#include <art/Framework/Principal/Handle.h>

// cetlib_except
#include <cetlib_except/exception.h>

// mu2e
#include <Offline/MCDataProducts/inc/StrawDigiMC.hh>
#include <Offline/TrackerConditions/inc/StrawElectronics.hh>
#include <Offline/TrackerMC/inc/StrawDigiBundle.hh>

namespace mu2e{
  using SDBC_iterator = std::vector<StrawDigiBundle>::iterator;
  using SDBC_const_iterator = std::vector<StrawDigiBundle>::const_iterator;

  class StrawDigiBundleCollection{
    public:
      // forwarded calls to underlying container
      size_t size() const;
      SDBC_iterator begin();
      SDBC_const_iterator begin() const;
      SDBC_iterator end();
      SDBC_const_iterator end() const;
      StrawDigiBundle& operator[](size_t i);
      const StrawDigiBundle& operator[](size_t i) const;

      // convenience methods for accepting new StrawDigiBundles
      // from different source situations
      template<template<typename> class Pointer>
      void Append(const Pointer<StrawDigiCollection>& digis,
                  const Pointer<StrawDigiADCWaveformCollection>& adcss,
                  const Pointer<StrawDigiMCCollection>& mcs);
      template<template<typename> class Pointer>
      void Append(const Pointer<StrawDigiCollection>& digis,
                  const Pointer<StrawDigiADCWaveformCollection>& adcss);
      void Append(const StrawDigiBundle bundle);

      // utility functions
      void FillStrawDigis(StrawDigiCollection& rv);
      void FillStrawDigiADCWaveforms(StrawDigiADCWaveformCollection& rv);
      void FillStrawDigiMCs(StrawDigiMCCollection& rv);
      StrawDigiBundleCollection operator+= (const StrawDigiBundleCollection& other);

      // accessors, as objects
      StrawDigiCollection GetStrawDigis();
      StrawDigiADCWaveformCollection GetStrawDigiADCWaveforms();
      StrawDigiMCCollection GetStrawDigiMCs();

      // accessors, as pointers
      std::unique_ptr<StrawDigiCollection> GetStrawDigiPtrs();
      std::unique_ptr<StrawDigiADCWaveformCollection> GetStrawDigiADCWaveformPtrs();
      std::unique_ptr<StrawDigiMCCollection> GetStrawDigiMCPtrs();

      // identify sets of  digis with overlapping digitization windows, and
      // reduce each such set to a single digi, representing their "sum"
      StrawDigiBundleCollection ResolveCollisions(const StrawElectronics& conditions);
    protected:
      std::vector<StrawDigiBundle> bundles;
      StrawDigiBundleCollection ResolveCollision(StrawDigiBundleCollection& collided);
    private:
      /**/
  };

  // forward size query to underlying container
  size_t StrawDigiBundleCollection::size() const{
    auto rv = this->bundles.size();
    return rv;
  }

  // forward iterator access to underlying bundles
  SDBC_iterator StrawDigiBundleCollection::begin(){
    auto rv = this->bundles.begin();
    return rv;
  }

  SDBC_const_iterator StrawDigiBundleCollection::begin() const{
    auto rv = this->bundles.begin();
    return rv;
  }

  SDBC_iterator StrawDigiBundleCollection::end(){
    auto rv = this->bundles.end();
    return rv;
  }

  SDBC_const_iterator StrawDigiBundleCollection::end() const{
    auto rv = this->bundles.end();
    return rv;
  }

  // convenience methods for accepting new StrawDigiBundles
  // from different source situations
  template<template<typename> class Pointer>
  void StrawDigiBundleCollection::Append(const Pointer<StrawDigiCollection>& digis,
                                         const Pointer<StrawDigiADCWaveformCollection>& adcss,
                                         const Pointer<StrawDigiMCCollection>& mcs){
      if ((digis->size() != adcss->size()) || (digis->size() != mcs->size())){
        std::string msg = "Attempting to append unsynchronized StrawDigi triplets: lengths = " + std::to_string(digis->size()) + ", " + std::to_string(adcss->size()) + ", " + std::to_string(mcs->size());
        throw cet::exception("TRIPLET SIZE MISMATCH") << msg << std::endl;
      }
      for (size_t i = 0 ; i < digis->size() ; i++){
        auto digi = digis->at(i);
        auto adcs = adcss->at(i);
        auto mc   = mcs->at(i);
        this->bundles.emplace_back(digi, adcs, mc);
      }
  }

  template<template<typename> class Pointer>
  void StrawDigiBundleCollection::Append(const Pointer<StrawDigiCollection>& digis,
                                         const Pointer<StrawDigiADCWaveformCollection>& adcss){
      if (digis->size() != adcss->size()){
        std::string msg = "Attempting to append unsynchronized StrawDigi doublets: lengths = " + std::to_string(digis->size()) + ", " + std::to_string(adcss->size());
        throw cet::exception("DOUBLET SIZE MISMATCH") << msg << std::endl;
      }
      for (size_t i = 0 ; i < digis->size() ; i++){
        const auto& digi = digis->at(i);
        const auto& adcs = adcss->at(i);
        StrawDigiBundle bundle(digi, adcs);
        this->bundles.emplace_back(digi, adcs);
        continue;
      }
  }
}

#endif

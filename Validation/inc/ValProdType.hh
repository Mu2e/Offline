#ifndef ValProdType_HH_
#define ValProdType_HH_


//
// This class looks for products of a certain type in the event
// and makes and fills a set of validation histograms for each
// instance
//

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include <vector>
#include <memory>


namespace mu2e {

  template <class T, class V>  class ValProdType {

  public:

    explicit ValProdType() {};
    int analyze( art::Event const&  event  );

  private:

    art::ServiceHandle<art::TFileService> _tfs;

    // the object containing the histograms is V
    // this is a pointer to that object
    typedef std::shared_ptr<V> Val;
    // need a vector - one set of histograms for each product instance
    typedef std::vector<Val> vVal;
    // the vector of pointers to histogram containers, V,  for 
    // products of type T
    vVal _val;

  };

}



template <class T, class V>
int mu2e::ValProdType<T,V>::analyze(art::Event const& event){

  // get all instances of products of type T
  std::vector<art::Handle< T >> vah;
  event.getManyByType(vah);


  std::string name;
  // loop over the list of instances of products of this type
  for (auto const & ah : vah) {
    const art::Provenance* prov = ah.provenance();

    std::string inst = prov->productInstanceName();
    if(inst.size()==0) inst="noName";
    name = prov->friendlyClassName()+"_"+
      prov->moduleLabel()+"_"+ inst;
    if(name.find("mu2e::",0)==0) name.erase(0,6);

    // see if this instance of this product is already in our list 
    // of products being histogrammed.  If not, add it to the list
    std::shared_ptr< V > prd = nullptr;
    for (auto ptr : _val ) {
      // if this instance of this product found in the list
      if(ptr->name().compare(name)==0)  prd = ptr;
    }
    // if not in the list, create a new set of histograms 
    // for this product
    if ( prd == nullptr ) {
      prd = std::make_shared<V>(name);
      // create histograms
      art::TFileDirectory tfdir = _tfs->mkdir(name);
      prd->declare(tfdir);
      // add it to the list of products being histogrammed
      _val.push_back( prd );
    }
    // histogram this event for this product
    prd->fill(*ah);
  }
  return 0;
}


# endif

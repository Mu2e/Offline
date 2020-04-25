#ifndef MVATools_HH
#define MVATools_HH

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "DataProducts/inc/MVAMask.hh"

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <vector>
#include <string>


namespace mu2e 
{

  class MVATools
  {
     public:
       struct Config
       {
	  fhicl::Atom<std::string> weights{ fhicl::Name("MVAWeights"), fhicl::Comment("MVA Weights xml file")};
       };       
       
       explicit MVATools(fhicl::ParameterSet const&);
       explicit MVATools(const Config& conf);

       virtual ~MVATools();
       void     initMVA();
       float    evalMVA(const std::vector<float>&,  const MVAMask& vmask=0xffffffff) const;
       float    evalMVA(const std::vector<double>&, const MVAMask& vmask=0xffffffff) const;
       void     showMVA() const;
       
       const std::vector<std::string>& titles() const { return title_;}     
       const std::vector<std::string>& labels() const { return label_;}     
 
 
    private:       
       enum   aType {null, tanh, sigmoid, relu};
       
       void   getGen(xercesc::DOMDocument* xmlDoc);
       void   getOpts(xercesc::DOMDocument* xmlDoc);
       void   getNorm(xercesc::DOMDocument* xmlDoc);
       void   getWgts(xercesc::DOMDocument* xmlDoc);
       float  activation(float arg) const;

       mutable std::vector<float> x_;
       mutable std::vector<float> y_;
       mutable std::vector<float> fv_;       
       std::vector<float>         wgts_;
       std::vector<unsigned>      links_;
       unsigned                   maxNeurons_;
       aType                      activeType_;
       bool                       oldMVA_;
       bool                       isNorm_;       
       std::vector<float>         voffset_;
       std::vector<float>         vscale_;
       std::vector<std::string>   title_;
       std::vector<std::string>   label_;
       std::string                activationTypeString_;
       std::string                mvaWgtsFile_;
  };
}
#endif

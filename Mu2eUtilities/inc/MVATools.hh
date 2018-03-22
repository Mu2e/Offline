#ifndef MVATools_HH
#define MVATools_HH

// framework
#include "fhiclcpp/ParameterSet.h"

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
       
       explicit MVATools(fhicl::ParameterSet const&);
       virtual ~MVATools();
       void     initMVA();
       float    evalMVA(const std::vector<float>&) const;
       float    evalMVA(const std::vector<double>&) const;
       void     showMVA()const;
     
 
    private:
       
       enum   aType {null, tanh,sigmoid,relu};
       
       void   getGen(xercesc::DOMDocument* xmlDoc);
       void   getOpts(xercesc::DOMDocument* xmlDoc);
       void   getNorm(xercesc::DOMDocument* xmlDoc);
       void   getWgts(xercesc::DOMDocument* xmlDoc);
       float activation(float arg) const;

       std::string mvaWgtsFile_;
       std::vector<std::vector<float>> wgts_;
       std::vector<float> voffset_;
       std::vector<float> vscale_;
       std::vector<std::string> title_;
       std::vector<std::string> label_;
       aType activeType_;
       bool oldMVA_;
       bool isNorm_;
       
       
       std::vector<unsigned> layerToNeurons_;
       std::vector<unsigned> synapsessPerLayer_;       
       
       // local workspace variables
       mutable std::vector<float> x_;
       mutable std::vector<float> y_;
  };
}
#endif

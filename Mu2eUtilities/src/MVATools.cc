#include <exception>

#include <cmath>
#include <stdlib.h>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "DataProducts/inc/MVAMask.hh"
#include "cetlib_except/exception.h"
// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "xercesc/dom/DOMDocument.hpp"
#include "xercesc/dom/DOMElement.hpp"
#include "xercesc/dom/DOMNamedNodeMap.hpp"
#include "xercesc/dom/DOMNode.hpp"
#include "xercesc/dom/DOMNodeList.hpp"
#include "xercesc/dom/DOMXPathResult.hpp"
#include "xercesc/parsers/AbstractDOMParser.hpp"
#include "xercesc/util/XMLException.hpp"
#include "xercesc/util/XMLString.hpp"
#include "xercesc/util/XercesDefs.hpp"
#include "xercesc/util/Xerces_autoconf_config.hpp"

#include "Mu2eUtilities/inc/MVATools.hh"

using namespace xercesc;

namespace mu2e
{

  MVATools::MVATools(const Config& config) :
    x_(),
    y_(),
    fv_(),
    wgts_(),
    maxNeurons_(0),
    activeType_(aType::null),
    oldMVA_(false),
    isNorm_(false),
    voffset_(),
    vscale_(),
    title_(),
    label_(),
    activationTypeString_("none"),
    mvaWgtsFile_()
  {
     ConfigFileLookupPolicy configFile;
     std::string weights = config.weights();
     mvaWgtsFile_ = configFile(weights);
  }

  MVATools::MVATools(fhicl::ParameterSet const& pset) :
    x_(),
    y_(),
    fv_(),
    wgts_(),
    maxNeurons_(0),
    activeType_(aType::null),
    oldMVA_(false),
    isNorm_(false),
    voffset_(),
    vscale_(),
    title_(),
    label_(),
    activationTypeString_("none"),
    mvaWgtsFile_()
  {
     ConfigFileLookupPolicy configFile;
     std::string weights = pset.get<std::string>("MVAWeights");
     mvaWgtsFile_ = configFile(weights);
  }

  MVATools::MVATools(const std::string& xmlfilename) :
    x_(),
    y_(),
    fv_(),
    wgts_(), 
    maxNeurons_(0), 
    activeType_(aType::null),
    oldMVA_(false),
    isNorm_(false),
    voffset_(), 
    vscale_(), 
    title_(), 
    label_(),
    activationTypeString_("none"),
    mvaWgtsFile_() { 

    ConfigFileLookupPolicy configFile;
    mvaWgtsFile_ = configFile(xmlfilename);
  }


  MVATools::~MVATools() {
    XMLPlatformUtils::Terminate();
  }

  void MVATools::initMVA()
  {
    if (wgts_.size()>0) throw cet::exception("RECO")<<"mu2e::MVATools: already initialized" << std::endl;

    xercesc::DOMDocument* xmlDoc = getXmlDoc();
    getGen(xmlDoc);
    getOpts(xmlDoc);
    getNorm(xmlDoc);
    getWgts(xmlDoc);

    xmlDoc->release();
  }

  void MVATools::getGen(xercesc::DOMDocument* xmlDoc)
  {
      XMLCh* ATT_GENERAL = XMLString::transcode("GeneralInfo");
      XMLCh* ATT_INFO = XMLString::transcode("Info");
      XMLCh* TAG_NAME = XMLString::transcode("name");
      XMLCh* TAG_VALUE = XMLString::transcode("value");

      XMLCh *xpathStr = XMLString::transcode("/MethodSetup/GeneralInfo");
      DOMXPathResult* xpathRes = xmlDoc->evaluate(xpathStr,xmlDoc->getDocumentElement(),NULL,
                                                  DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE, NULL);
      XMLString::release(&xpathStr) ;
      DOMElement* optElem = dynamic_cast<DOMElement* >(xpathRes->getNodeValue());
      xpathRes->release();


      DOMNodeList* children = optElem->getElementsByTagName(ATT_INFO);
      for (XMLSize_t ix = 0 ; ix < children->getLength() ; ++ix )
      {
          DOMNode* childNode = children->item(ix) ;
          DOMNamedNodeMap* attrs = childNode->getAttributes();

          std::string label,value;
          for( XMLSize_t ia = 0 ; ia < attrs->getLength() ; ++ia )
          {
	     DOMNode* attr = attrs->item(ia);
	     char* attValue = XMLString::transcode(attr->getNodeValue());
	     if (XMLString::equals(TAG_NAME,attr->getNodeName()) )     label = std::string(attValue);
	     if (XMLString::equals(TAG_VALUE,attr->getNodeName()) )    value = std::string(attValue);
             XMLString::release(&attValue);
          }
	  if (label.find("TMVA Release") != std::string::npos)
	  {
 	     std::string code = value.substr(value.find("[")+1,value.find("]")-value.find("[")-1);
	     int iversion = atoi(code.c_str());
             if (iversion < 262657) oldMVA_ = true;
	  }
      }

      XMLString::release(&ATT_GENERAL);
      XMLString::release(&ATT_INFO);
      XMLString::release(&TAG_NAME);
      XMLString::release(&TAG_VALUE);
  }

  void MVATools::getOpts(xercesc::DOMDocument* xmlDoc)
  {
      XMLCh* ATT_OPTION = XMLString::transcode("Option");
      XMLCh* TAG_NAME = XMLString::transcode("name");
      XMLCh* TAG_MODIFIED = XMLString::transcode("modified");

      XMLCh *xpathStr = XMLString::transcode("/MethodSetup/Options");
      DOMXPathResult* xpathRes = xmlDoc->evaluate(xpathStr,xmlDoc->getDocumentElement(),NULL,
                                                  DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE, NULL);
      XMLString::release(&xpathStr) ;
      DOMElement* optElem = dynamic_cast<DOMElement* >(xpathRes->getNodeValue());
      xpathRes->release();

      DOMNodeList* children = optElem->getElementsByTagName(ATT_OPTION);
      for (XMLSize_t ix = 0 ; ix < children->getLength() ; ++ix )
      {
          DOMNode* childNode = children->item(ix) ;
          DOMNamedNodeMap* attrs = childNode->getAttributes();
          char* value = XMLString::transcode(childNode->getFirstChild()->getNodeValue());

          std::string label, modified;
          for( XMLSize_t ia = 0 ; ia < attrs->getLength() ; ++ia )
          {
	     DOMNode* attr = attrs->item(ia);
	     char* attValue = XMLString::transcode(attr->getNodeValue());
	     if (XMLString::equals(TAG_NAME,attr->getNodeName()) ) label = std::string(attValue);
	     if (XMLString::equals(TAG_MODIFIED,attr->getNodeName()) ) modified = std::string(attValue);
             XMLString::release(&attValue);
          }
          if (label.find("NeuronType") != std::string::npos)
	  {
	     std::string val(value);
	     if (val.find("tanh") != std::string::npos)    activeType_ = aType::tanh;
	     if (val.find("sigmoid") != std::string::npos) activeType_ = aType::sigmoid;
	     if (val.find("ReLU") != std::string::npos)    activeType_ = aType::relu;
             activationTypeString_ = val;
	  }
          if (label.find("VarTransform") != std::string::npos)
	  {
	     std::string val(value);
	     if (modified.find("Yes") != std::string::npos) isNorm_ = true;
	     if (isNorm_ && val.find("N") == std::string::npos)
               throw cet::exception("RECO")<<"mu2e::MVATools: unknown normalization mode" << std::endl;
	  }
      }

      if (activeType_ == aType::null) throw cet::exception("RECO")<<"mu2e::MVATools: unknown activation function" << std::endl;

      XMLString::release(&ATT_OPTION);
      XMLString::release(&TAG_NAME);
      XMLString::release(&TAG_MODIFIED);
  }

  void MVATools::getNorm(xercesc::DOMDocument* xmlDoc)
  {
      XMLCh* TAG_VARIABLE = XMLString::transcode("Variable");
      XMLCh* ATT_MIN = XMLString::transcode("Min");
      XMLCh* ATT_MAX = XMLString::transcode("Max");
      XMLCh* ATT_TITLE = XMLString::transcode("Title");
      XMLCh* ATT_LABEL = XMLString::transcode("Label");

      XMLCh *xpathStr = XMLString::transcode("/MethodSetup/Variables");
      DOMXPathResult* xpathRes = xmlDoc->evaluate(xpathStr,xmlDoc->getDocumentElement(),NULL,
                                                  DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE, NULL);
      XMLString::release(&xpathStr);
      DOMElement* varElem = dynamic_cast<DOMElement* >(xpathRes->getNodeValue()) ;
      xpathRes->release();

      DOMNodeList* children = varElem->getElementsByTagName(TAG_VARIABLE);
      for( XMLSize_t ix = 0 ; ix < children->getLength() ; ++ix )
      {
         float vmin(std::numeric_limits<float>::max());
         float vmax(0.0);
         std::string label,title;

         DOMNode* childNode = children->item( ix ) ;
         DOMNamedNodeMap* attrs = childNode->getAttributes();
         for( XMLSize_t ia = 0 ; ia < attrs->getLength() ; ++ia )
         {
	    DOMNode* attr = attrs->item(ia);
	    char* attValue = XMLString::transcode(attr->getNodeValue());

	    if (XMLString::equals(ATT_MIN,attr->getNodeName()) )      vmin = strtof(attValue,NULL);
	    if (XMLString::equals(ATT_MAX,attr->getNodeName()) )      vmax = strtof(attValue,NULL);
	    if (XMLString::equals(ATT_TITLE,attr->getNodeName()) )    title = std::string(attValue);
	    if (XMLString::equals(ATT_LABEL,attr->getNodeName()) )    label = std::string(attValue);

            XMLString::release( &attValue ) ;
         }
         voffset_.push_back(vmin);
         vscale_.push_back(2.0/(vmax-vmin));
         title_.push_back(title);
         label_.push_back(label);
      }

      XMLString::release(&TAG_VARIABLE);
      XMLString::release(&ATT_MIN);
      XMLString::release(&ATT_MAX);
      XMLString::release(&ATT_LABEL);
      XMLString::release(&ATT_TITLE);
  }


  xercesc::DOMDocument* MVATools::getXmlDoc() {

    try
      {
	XMLPlatformUtils::Initialize();
      } 
    catch (XMLException& e)
      {
	char* message = XMLString::transcode( e.getMessage() ) ;
	throw cet::exception("RECO")<<"mu2e::MVATools: XML initialization error: " <<  message << std::endl;
	XMLString::release( &message ) ;
      }
  
    XercesDOMParser* parser = new XercesDOMParser();
    parser->setValidationScheme(XercesDOMParser::Val_Never);
    parser->setDoNamespaces(false);
    parser->setDoSchema(false);
    parser->setLoadExternalDTD( false );
  
    XMLCh *xmlFile = XMLString::transcode(mvaWgtsFile_.c_str());
    parser->parse(xmlFile);
    XMLString::release( &xmlFile ) ;
  
    xercesc::DOMDocument* xmlDoc = parser->adoptDocument() ; // adopt the document so that the parser no longer owns it...
    delete parser; // ...can then delete the parser and solve a memory leak

    if (!xmlDoc) {
      throw cet::exception("RECO") << "mu2e::MVATools could not create xmlDoc with filename " << mvaWgtsFile_.c_str() << std::endl;
    }

    return xmlDoc;
  }


  void MVATools::getWgts(xercesc::DOMDocument* xmlDoc)
  {
      XMLCh* ATT_NSYNAPSES = XMLString::transcode("NSynapses");
      XMLCh* ATT_INDEX = XMLString::transcode("Index");

      XMLCh *xpathStr = XMLString::transcode("/MethodSetup/Weights/Layout");
      DOMXPathResult* xpathRes = xmlDoc->evaluate(xpathStr,xmlDoc->getDocumentElement(),NULL,
                                                  DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE, NULL);
      XMLString::release(&xpathStr) ;
      xpathRes->release();


      xpathStr = XMLString::transcode("/MethodSetup/Weights/Layout/Layer/Neuron");
      xpathRes = xmlDoc->evaluate(xpathStr,xmlDoc->getDocumentElement(), NULL,
                                  DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE, NULL);
      XMLString::release(&xpathStr);

      int iCurrentLayer(-1);
      unsigned nNeurons(0);
      std::vector<std::vector<float> > wtemp;
      for (XMLSize_t i=0;i<xpathRes->getSnapshotLength();i++)
      {
          xpathRes->snapshotItem(i);

          DOMNode*    nNeuron = xpathRes->getNodeValue();
          DOMElement* neuronElement = dynamic_cast<DOMElement* >(nNeuron);
          DOMNode*    parentNode = nNeuron->getParentNode();
          DOMElement* parentElement = dynamic_cast<DOMElement* >(parentNode);

          char* layerIndex = XMLString::transcode(parentElement->getAttribute(ATT_INDEX));
          int   iLayer     = atoi(layerIndex);
          char* nSynapses  = XMLString::transcode(neuronElement->getAttribute(ATT_NSYNAPSES));
          unsigned iSynapses  = atoi(nSynapses);

          if (iLayer != iCurrentLayer)
	  {
              if (iCurrentLayer > -1) links_.push_back(nNeurons);
	      iCurrentLayer=iLayer;
              nNeurons=0;
	  }

	  if (nNeuron->getFirstChild()==NULL) continue; //output neuron has no children

	  char* cw = XMLString::transcode(nNeuron->getFirstChild()->getNodeValue());
	  std::string sw(cw),temp;
          std::stringstream ss(sw);
          std::vector<float> w;
	  while(ss >> temp) w.push_back(strtof(temp.c_str(),NULL));
	  if (w.size() != iSynapses) throw cet::exception("RECO")<<"mu2e::MVATools: internal error" << std::endl;
	  wtemp.push_back(std::move(w));

	  XMLString::release(&cw);
	  ++nNeurons;
      }

      // the weights are given in a back-propagation style, continuous weights are from a neuron of the previous layer to all neurons of the next layer
      // for evaluation, forward-propagation style is better, contiguous weights from all neurons of the previous layer to a single neuron of the next layer
      // we also need the number of synapses / layer, which is given by the number of neurons in the next layer -1 to remove bias neuron
      //
      std::vector<unsigned> layerToNeurons(1,0),synapsessPerLayer;
      for (const auto& val : links_)    layerToNeurons.push_back(val+layerToNeurons.back());
      for (size_t i=1;i<links_.size();++i) synapsessPerLayer.push_back(links_[i]-1);
      synapsessPerLayer.push_back(1);

      for (unsigned iLayer=1;iLayer < layerToNeurons.size(); ++iLayer)
      {
	  for (unsigned iOut=0;iOut<synapsessPerLayer[iLayer-1];++iOut)
	  {
	     std::vector<float> temp;
	     for (unsigned iIn=layerToNeurons[iLayer-1];iIn<layerToNeurons[iLayer];++iIn) wgts_.push_back(wtemp[iIn][iOut]);
	  }
      }

      maxNeurons_ = *std::max_element(links_.begin(),links_.end());
      x_ = std::vector<float>(maxNeurons_,0.0f);
      y_ = std::vector<float>(maxNeurons_,0.0f);

      XMLString::release(&ATT_INDEX);
      XMLString::release(&ATT_NSYNAPSES);
  }
  
  void MVATools::getCalib(std::map<float, float>& effCalib) {

    xercesc::DOMDocument* xmlDoc = getXmlDoc();

    XMLCh* TAG_CALIBRATION = XMLString::transcode("Calib");
    XMLCh* ATT_INDEX = XMLString::transcode("Index");
    XMLCh* ATT_EFF = XMLString::transcode("CalibVal");
    XMLCh* ATT_CUT = XMLString::transcode("Val");
    
    XMLCh *xpathStr = XMLString::transcode("/MethodSetup/Calibration");
    DOMXPathResult* xpathRes = xmlDoc->evaluate(xpathStr,xmlDoc->getDocumentElement(),NULL,
						DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE, NULL);
    XMLString::release(&xpathStr);
    DOMElement* varElem = dynamic_cast<DOMElement* >(xpathRes->getNodeValue()) ;
    if (!varElem) {
      throw cet::exception("MVATools") << "No calibration for this MVA (" << mvaWgtsFile_ << ")" << std::endl;
    }
    xpathRes->release();
    
    DOMNodeList* children = varElem->getElementsByTagName(TAG_CALIBRATION);
    for( XMLSize_t ix = 0 ; ix < children->getLength() ; ++ix ) {       
      float eff(0.0);
      float cut(0.0);
      
      DOMNode* childNode = children->item( ix ) ;
      DOMNamedNodeMap* attrs = childNode->getAttributes();
      for( XMLSize_t ia = 0 ; ia < attrs->getLength() ; ++ia ) {
	DOMNode* attr = attrs->item(ia);
	char* attValue = XMLString::transcode(attr->getNodeValue());
	
	if (XMLString::equals(ATT_EFF,attr->getNodeName()) )      eff = strtof(attValue,NULL);
	if (XMLString::equals(ATT_CUT,attr->getNodeName()) )      cut = strtof(attValue,NULL);
	
	XMLString::release( &attValue ) ;
      }

      effCalib.insert(std::pair<float, float>(eff, cut));      
    }
  
    XMLString::release(&TAG_CALIBRATION);
    XMLString::release(&ATT_INDEX);
    XMLString::release(&ATT_EFF);
    XMLString::release(&ATT_CUT);
    xmlDoc->release();
  }


  float MVATools::evalMVA(const std::vector<double >& v, const MVAMask& mask) const
  {
     for (size_t i=0;i<v.size();++i)  fv_[i] = static_cast<float>(v[i]);
     return evalMVA(fv_,mask);
  }

  float MVATools::evalMVA(const std::vector<float>& v, const MVAMask& mask) const
  {

      // Normalize the input data and add the bias node, skip masked values
      size_t ival(0);
      for (size_t ivar=0; ivar < v.size(); ivar++)
      {
         if ( mask & (1<<ivar) )
         {
	    x_[ival]= isNorm_ ? (v[ivar]-voffset_[ival])*vscale_[ival] - 1.0 : v[ivar];
	    ++ival;
         }
      }
      x_[ival] = 1.0;

      if (ival != links_[0]-1)
	throw cet::exception("RECO")<<"mu2e::MVATools: mismatch input dimension (ival = " << ival << ") and network architecture (links_[0]-1 = " << links_[0]-1 << ")" << std::endl;


      //perform feed forward calculation up to the last hidden layer
      unsigned idxWeight(0);
      for (unsigned k=0;k<links_.size()-1;++k)
      {
          //the number of synpases is given by the number of neurons in the next layer -1 (do not count bias neuron!)
          for (unsigned j=0;j<links_[k+1]-1;++j)
          {
	     y_[j]=0.0f;
	     for (unsigned i=0;i<links_[k];++i) y_[j] += wgts_[i+idxWeight]*x_[i];
             y_[j] = activation(y_[j]);
             idxWeight += links_[k];
          }
          x_.swap(y_);
          x_[links_[k+1]-1] = 1.0f; //add bias neuron
      }

      //calculate output neuron value
      float yf(0.0);
      for (unsigned i=0;i<links_.back();++i) yf += wgts_[i+idxWeight]*x_[i];

      if (oldMVA_) return yf;
      return  1.0/(1.0+expf(-yf));
  }




  float MVATools::activation(float arg) const
  {
     if (activeType_== aType::tanh)
     {
       if (oldMVA_) return std::tanh(arg);
       if (arg > 4.97) return 1.0;
       if (arg < -4.97) return -1.0;
       float arg2 = arg * arg;
       float a = arg * (135135.0f + arg2 * (17325.0f + arg2 * (378.0f + arg2)));
       float b = 135135.0f + arg2 * (62370.0f + arg2 * (3150.0f + arg2 * 28.0f));
       return a/b;
     }
     if (activeType_== aType::sigmoid) return 1.0/(1.0+expf(-arg));
     if (activeType_== aType::relu) return std::max(0.0f,arg);

     return -999.0;
  }


  void MVATools::showMVA() const
  {
      std::cout << "MVA weights from file:" <<     mvaWgtsFile_ << std::endl;;
      std::cout << "MVA NLayers: " << links_.size()-1 << std::endl;;
      std::cout << "MVA NVars: " << title_.size() << std::endl;;
      std::cout << "MVA Activation type: " << activationTypeString_ << std::endl;;

      std::cout.setf(std::ios::scientific);
      std::cout.precision(7);

      const std::string stars1(12,'*');
      const std::string label1 = " MVA Normalization ";
      std::cout << stars1 << label1 << stars1 << std::endl;;
      for (size_t i = 0; i <label_.size(); ++i)
        std::cout << "Var " << i <<" "<< label_[i] << " " << title_[i] << ": min=" << voffset_[i] << " max=" << 2.0/vscale_[i]+voffset_[i] << std::endl;;

      const std::string morestars1(24+label1.size(),'*');
      std::cout << morestars1 << std::endl;;

      const std::string stars2(23,'*');
      const std::string label2 = " MVA Weights ";
      std::cout << stars2 << label2 << stars2 << std::endl;


      int idx(0);
      for (unsigned k=0;k<links_.size()-1;++k)
      {
        std::cout<<"Layer "<<k<<std::endl;
        for (unsigned j=0;j<links_[k+1]-1;++j)
        {
	  std::cout<<"Synapses 1.."<<links_[k]<<" of current layer to synapse "<<j<<" of next layer"<<std::endl;
	  for (unsigned i=idx;i<idx+links_[k];++i)std::cout<<wgts_[i]<<" ";
	  std::cout<<std::endl;
          idx += links_[k];
        }
      }

      std::cout<<"Layer "<<links_.size()-1<<std::endl;
      std::cout<<"Synapses 1.."<<links_.back()<<" of current layer to synapse 0 of next layer"<<std::endl;
      for (unsigned i=idx;i<idx+links_.back();++i)std::cout<<wgts_[i]<<" ";
      std::cout<<std::endl;

      const std::string morestars2(46+label2.size(),'*');
      std::cout << morestars2 << std::endl;
  }

}


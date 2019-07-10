#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Mu2eUtilities/inc/MVATools.hh"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>
#include <limits>

using namespace xercesc;

namespace mu2e 
{


MVATools::MVATools(const Config& config) : 
  mvaWgtsFile_(), 
  wgts_(), 
  x_(), 
  y_(),
  fv_(), 
  layerToNeurons_(), 
  synapsessPerLayer_(),
  voffset_(), 
  vscale_(), 
  title_(), 
  label_(),
  activeType_(aType::null),
  activationTypeString_("none"),
  oldMVA_(false),
  isNorm_(false)
{
   ConfigFileLookupPolicy configFile;
   std::string weights = config.weights();
   mvaWgtsFile_ = configFile(weights);
}

MVATools::MVATools(fhicl::ParameterSet const& pset) : 
  mvaWgtsFile_(), 
  wgts_(), 
  x_(), 
  y_(),
  fv_(), 
  layerToNeurons_(), 
  synapsessPerLayer_(),
  voffset_(), 
  vscale_(), 
  title_(), 
  label_(),
  activeType_(aType::null),
  activationTypeString_("none"),
  oldMVA_(false),
  isNorm_(false)
{
   ConfigFileLookupPolicy configFile;
   std::string weights = pset.get<std::string>("MVAWeights");
   mvaWgtsFile_ = configFile(weights);
}

MVATools::~MVATools() {}



void MVATools::initMVA()
{
    if (wgts_.size()>0) throw cet::exception("RECO")<<"mu2e::MVATools: already initialized" << std::endl;

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
 
    xercesc::DOMDocument* xmlDoc = parser->getDocument() ;

    getGen(xmlDoc);
    getOpts(xmlDoc);
    getNorm(xmlDoc);
    getWgts(xmlDoc);
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
	   if (val.find("tanh") != std::string::npos) activeType_ = aType::tanh;  
	   if (val.find("sigmoid") != std::string::npos) activeType_ = aType::sigmoid;  
	   if (val.find("ReLU") != std::string::npos) activeType_ = aType::relu;
           activationTypeString_ = val;  	   
	}      
        if (label.find("VarTransform") != std::string::npos)
	{
	   std::string val(value);
	   if (modified.find("Yes") != std::string::npos) isNorm_=true;
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

        DOMNode* nNeuron = xpathRes->getNodeValue();
        DOMElement* neuronElement = dynamic_cast<DOMElement* >(nNeuron);
        DOMNode* parentNode = nNeuron->getParentNode();
        DOMElement* parentElement = dynamic_cast<DOMElement* >(parentNode);

        char* layerIndex = XMLString::transcode(parentElement->getAttribute(ATT_INDEX));
        int   iLayer     = atoi(layerIndex);
        char* nSynapses  = XMLString::transcode(neuronElement->getAttribute(ATT_NSYNAPSES));
        unsigned iSynapses  = atoi(nSynapses);

        if (iLayer != iCurrentLayer) 
	{
	    layerToNeurons_.push_back(nNeurons); 
	    if (iSynapses) synapsessPerLayer_.push_back(iSynapses); 
	    iCurrentLayer=iLayer;
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

    for (unsigned iLayer=1;iLayer < layerToNeurons_.size(); ++iLayer)
    {    
	for (unsigned iOut=0;iOut<synapsessPerLayer_[iLayer-1];++iOut)
	{        
	   std::vector<float> temp;
	   for (unsigned iIn=layerToNeurons_[iLayer-1];iIn<layerToNeurons_[iLayer];++iIn) 
              temp.push_back(wtemp[iIn][iOut]);                    
	   wgts_.push_back(std::move(temp));       
	}
    }
    
    unsigned maxSynapses=*std::max_element(layerToNeurons_.begin(),layerToNeurons_.end());
    for (unsigned i=0;i<layerToNeurons_[1];++i) x_.push_back(1);
    for (unsigned i=0;i<=maxSynapses;++i) y_.push_back(0);
    
    for (unsigned i=0;i<layerToNeurons_[1]-1;++i) fv_.push_back(0);
    
    
    XMLString::release(&ATT_INDEX);
    XMLString::release(&ATT_NSYNAPSES);    
}




float MVATools::evalMVA(const std::vector<double >& v,MVAMask mask) const 
{
   for (size_t i=0;i<v.size();++i)  fv_[i] = static_cast<float>(v[i]);   
   return evalMVA(fv_,mask);
}

float MVATools::evalMVA(const std::vector<float>& v,MVAMask mask) const 
{
    // Normalize the bias node is 
    // skip variables not masked
    size_t ival(0);
    unsigned vsize = voffset_.size();
    for (size_t ivar=0; ivar < v.size(); ivar++){
      if( (mask&(1<<ivar))){
	x_[ival]= isNorm_ ? (v[ivar]-voffset_[ival])*vscale_[ival] - 1.0 : v[ivar];
	ival++;
      }
    }
    x_[vsize] = 1.0;
        
    //forward propagation of internal layers
    unsigned idxWeight(0);
    for (unsigned k=0;k<synapsessPerLayer_.size()-1;++k)
    {          
      for (unsigned j=0;j<synapsessPerLayer_[k];++j)
      {
	y_[j]=0.0;
	for (unsigned i=0;i<wgts_[idxWeight].size();++i) y_[j] += wgts_[idxWeight][i]*x_[i];
        y_[j] = activation(y_[j]);
        ++idxWeight;
      }      
            
      x_.swap(y_); //faster than assignment
      x_[synapsessPerLayer_[k]] = 1.0; //add bias neuron
    }   
    
    //output layer
    float y(0.0);
    for (unsigned i=0;i<wgts_[idxWeight].size();++i) y += wgts_[idxWeight][i]*x_[i];

    if (ival != layerToNeurons_[1]-1) 
      throw cet::exception("RECO")<<"mu2e::MVATools: mismatch input dimension and network architecture" << std::endl;

  if (oldMVA_) return y;
    return  1.0/(1.0+expf(-y));
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
    std::cout << "MVA NLayers: " << synapsessPerLayer_.size() << std::endl;;
    std::cout << "MVA NVars: " << title_.size() << std::endl;;
    std::cout << "MVA Activation type: " << activationTypeString_ << std::endl;;

    std::cout.setf(std::ios::scientific);
    std::cout.precision(7);

    const std::string stars1(12,'*');
    const std::string label1 = " MVA Normalization ";
    std::cout << stars1 << label1 << stars1 << std::endl;;
    for (size_t i = 0; i <label_.size(); ++i)
      std::cout << "Var " << i << label_[i] << " " << title_[i] << ": min=" << voffset_[i] << " max=" << 2.0/vscale_[i]+voffset_[i] << std::endl;;
    
    const std::string morestars1(24+label1.size(),'*');
    std::cout << morestars1 << std::endl;;

    const std::string stars2(23,'*');
    unsigned idxWeight(0);
    const std::string label2 = " MVA Weights ";
    std::cout << stars2 << label2 << stars2 << std::endl;
    
    for (unsigned k=0;k<synapsessPerLayer_.size();++k)
    {          
      std::cout<<" Layer : "<<k<<std::endl;
      for (unsigned j=0;j<synapsessPerLayer_[k];++j)
      {
	std::cout<<"Synapses 1.."<<wgts_[idxWeight].size()<<" of current layer to synapse "<<j<<" of next layer"<<std::endl;
        for (unsigned i=0;i<wgts_[idxWeight].size();++i) std::cout <<wgts_[idxWeight][i]<<" ";
        std::cout<<std::endl;
        ++idxWeight;
      }      
    }   
    
    const std::string morestars2(46+label2.size(),'*');
    std::cout << morestars2 << std::endl;
}



}


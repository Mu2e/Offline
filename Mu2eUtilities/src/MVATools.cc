
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Mu2eUtilities/inc/MVATools.hh"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"

// Xerces XML Parser
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>

#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>

using namespace std;
using namespace xercesc;

namespace mu2e 
{
  MVATools::MVATools(fhicl::ParameterSet const& pset){
    // location-independent files
    ConfigFileLookupPolicy configFile;
    string weights = pset.get<std::string>("MVAWeights");
    _mvaWgtsFile = configFile(weights);
  }
  MVATools::~MVATools()
  {}

  // Initialize Xerces and the MVA
  void MVATools::initMVA(){
    // Initialize xerces
    try{
      XMLPlatformUtils::Initialize();
    }catch( XMLException& e ){
      char* message = XMLString::transcode( e.getMessage() ) ;
      throw cet::exception("RECO")<<"mu2e::MVATools: XML initialization error: " <<  message << endl;
      XMLString::release( &message ) ;
    }

    // Create the parser
    XercesDOMParser* parser = new XercesDOMParser();
    parser->setValidationScheme(XercesDOMParser::Val_Never);
    parser->setDoNamespaces(false);
    parser->setDoSchema(false);
    parser->setLoadExternalDTD( false );
  
    XMLCh *xmlFile = XMLString::transcode(_mvaWgtsFile.c_str());
    parser->parse(xmlFile);
    XMLString::release( &xmlFile ) ;

    _xmlDoc = parser->getDocument() ;

    // get the normalization for the mva weights
    getNorm();
    for(vector<double> norm : _wn){
      _wnr2.push_back(0.5*(norm[1]-norm[0]));
    }
    // get the weights
    getWgts();
  }

  // Find normalization ranges
  void MVATools::getNorm() {
    XMLCh* TAG_CLASS = XMLString::transcode("Class");
    XMLCh* TAG_VARIABLE = XMLString::transcode("Variable");
    XMLCh* TAG_VARINDEX = XMLString::transcode("VarIndex");
    XMLCh* TAG_CLASSINDEX = XMLString::transcode("ClassIndex");
    XMLCh* VAL_2 = XMLString::transcode("2");
    XMLCh* TAG_RANGE = XMLString::transcode("Range");
    XMLCh* ATT_NVARS = XMLString::transcode("NVar");
    XMLCh* ATT_INDEX = XMLString::transcode("Index");
    XMLCh* ATT_MIN = XMLString::transcode("Min");
    XMLCh* ATT_MAX = XMLString::transcode("Max");
    XMLCh* ATT_TITLE = XMLString::transcode("Title");
    XMLCh* ATT_LABEL = XMLString::transcode("Label");

    // Get the number of ANN(MLP) input variables
    XMLCh *xpathStr = XMLString::transcode("/MethodSetup/Variables");
    DOMXPathResult* xpathRes = _xmlDoc->evaluate(xpathStr,_xmlDoc->getDocumentElement(),NULL,
                               DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE, NULL);
    XMLString::release( &xpathStr ) ;

    DOMElement* varElem = dynamic_cast<DOMElement* >(xpathRes->getNodeValue()) ;

    char* attValue = XMLString::transcode(varElem->getAttribute(ATT_NVARS));
    int nVars = atoi(attValue);
//    XMLString::release( &attValue );
//    cout << "found " << nVars << " variables" << endl; 
    _wn.resize(nVars);
    _x.resize(nVars);
    _y.resize(nVars);
    _title.resize(nVars);
    _label.resize(nVars);

    DOMNodeList* children = varElem->getElementsByTagName(TAG_VARIABLE) ;
//    cout << "found " << children->getLength() << " children" << endl;
    for( XMLSize_t ix = 0 ; ix < children->getLength() ; ++ix ){
      DOMNode* childNode = children->item( ix ) ;
      DOMNamedNodeMap* attrs = childNode->getAttributes();
      int ivar=0;
      double vmin=std::numeric_limits<double>::max();
      double vmax=0.;
      string label,title;
      for( XMLSize_t ia = 0 ; ia < attrs->getLength() ; ++ia ){
	DOMNode* attr = attrs->item(ia);
	char* attValue = XMLString::transcode(attr->getNodeValue());
	if(XMLString::equals(TAG_VARINDEX,attr->getNodeName()) ){
	  ivar = atoi(attValue);
	}else if(XMLString::equals(ATT_MIN,attr->getNodeName()) ){
	  vmin = strtod(attValue,NULL);
	}else if(XMLString::equals(ATT_MAX,attr->getNodeName()) ){
	  vmax = strtod(attValue,NULL);
	}else if(XMLString::equals(ATT_TITLE,attr->getNodeName()) ){
	  title = string(attValue);
	}else if(XMLString::equals(ATT_LABEL,attr->getNodeName()) ){
	  label = string(attValue);
	}
        XMLString::release( &attValue ) ;
      }
      _wn[ivar].push_back(vmin);
      _wn[ivar].push_back(vmax);
      _title[ivar] = title;
      _label[ivar] = label;
    }
    xpathRes->release();

    XMLString::release(&TAG_CLASS);
    XMLString::release(&TAG_CLASSINDEX);
    XMLString::release(&TAG_RANGE);
    XMLString::release(&VAL_2);
    XMLString::release(&ATT_NVARS);
    XMLString::release(&TAG_VARINDEX);
    XMLString::release(&ATT_INDEX);
    XMLString::release(&ATT_MIN);
    XMLString::release(&ATT_MAX);
    XMLString::release(&ATT_LABEL);
    XMLString::release(&ATT_TITLE);
  }

  // Find the layer weights
  void MVATools::getWgts(){
    XMLCh* ATT_INDEX = XMLString::transcode("Index");
    XMLCh* ATT_NLYRS = XMLString::transcode("NLayers");

    // Get the number of ANN(MLP) layers
    XMLCh *xpathStr = XMLString::transcode("/MethodSetup/Weights/Layout");
    DOMXPathResult* xpathRes = _xmlDoc->evaluate(xpathStr,_xmlDoc->getDocumentElement(),NULL,
                               DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE, NULL);
    XMLString::release( &xpathStr ) ;

    DOMElement* layoutElem = dynamic_cast<DOMElement* >(xpathRes->getNodeValue()) ;
    xpathRes->release();

    char* attValue = XMLString::transcode(layoutElem->getAttribute(ATT_NLYRS));
    int nLyrs = atoi(attValue);
    XMLString::release( &attValue ) ;

    int iLastLyr = nLyrs - 1;
    _wgts.resize(iLastLyr);

    // Now get the weights for each layer
    xpathStr = XMLString::transcode("/MethodSetup/Weights/Layout/Layer/Neuron");
    xpathRes = _xmlDoc->evaluate(xpathStr,_xmlDoc->getDocumentElement(), NULL,
                                              DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE, NULL);
    XMLString::release( &xpathStr ) ;

    for(XMLSize_t i=0;i<xpathRes->getSnapshotLength();i++) {
      xpathRes->snapshotItem(i);
      DOMNode* nNeuron = xpathRes->getNodeValue();

      DOMNode* parentNode = nNeuron->getParentNode();
      DOMElement* parentElement = dynamic_cast<DOMElement* >(parentNode);

      char* layerIndex = XMLString::transcode(parentElement->getAttribute(ATT_INDEX));
      int ilayer = atoi(layerIndex);

      vector<double> wv;
      if(ilayer<iLastLyr){
        char* cw = XMLString::transcode(nNeuron->getFirstChild()->getNodeValue());
        string sw(cw);
        stringstream ss (sw);
        string temp;
        while(ss>>temp){
          wv.push_back(strtod(temp.c_str(),NULL));
        }
        _wgts[ilayer].push_back(wv);
        XMLString::release(&cw);
      }
    }

    _twgts.resize(iLastLyr);

    // Get transpose of 2d vector in each layer
    for(vector<vector<vector<double> > >::size_type i = 0; i != _wgts.size(); ++i){
      for(vector<vector<double> >::size_type col=0; col != _wgts[i][0].size(); ++col){
        vector<double>vrow;
        for(vector<double>vtmp: _wgts[i] ){
          vrow.push_back(vtmp[col]);
        }
        _twgts[i].push_back(vrow);
      }
    }
    xpathRes->release();
    XMLString::release(&ATT_NLYRS);
    XMLString::release(&ATT_INDEX);
  }

  double MVATools::evalMVA(vector<double>const &v) const {

    // Normalize
    for(vector<double>::size_type i = 0; i != _wn.size(); ++i){
      _x[i] = ((v[i]-_wn[i][0])/_wnr2[i]) - 1.;
    }

    // do internal layers
    size_t klast = _twgts.size()-1;
    for(size_t k = 0; k != klast; ++k){
      for(vector<vector<double> >::size_type j = 0; j != _twgts[k].size(); ++j){
        _y[j]=0.;
        size_t ilast = _twgts[k][j].size()-1;
        for(size_t i = 0; i != ilast; ++i){
          _y[j] += _x[i]*_twgts[k][j][i];
        }
        _y[j]+=_twgts[k][j][ilast];
        _y[j]=1./(1.+exp(-_y[j]));
      }
      _x=_y;
    }

    // do output
    double mva=0.;
    size_t ilast = _twgts[klast][0].size()-1;
    for(size_t i=0;i<ilast;i++){
      mva+=_y[i]*_twgts[klast][0][i];
    }
    mva+=_twgts[klast][0][ilast];
    return mva;
    
  }

  void MVATools::showMVA()const{
    
    cout << "MVA weights from file:" <<     _mvaWgtsFile << endl;
    cout << "MVA NLayers: " << _wgts.size() << endl;
    cout << "MVA NVars: " << _title.size() << endl;

    ios_base::fmtflags origflags = cout.flags();
    streamsize prec = cout.precision();
    streamsize width = cout.width(); 
    cout.setf(ios::scientific);
    cout.precision(7);

    const string stars1(12,'*');
    const string label1 = " MVA Normalization ";
    cout << stars1 << label1 << stars1 << endl;
    for(vector< vector<double> >::size_type i = 0; i != _wn.size(); ++i){
      cout << "Var " << i << _label[i] << " " << _title[i] << ": min=" << _wn[i][0] << " max=" << _wn[i][1] << endl;
    }
    const string morestars1(24+label1.size(),'*');
    cout << morestars1 << endl;

    const string stars2(23,'*');
    const string label2 = " MVA Weights ";
    cout << stars2 << label2 << stars2 << endl;
    cout << "From file: " << _mvaWgtsFile << endl;
    for(vector <vector< vector<double> > >::size_type k = 0; k != _twgts.size(); ++k){
      cout << "Layer " << k << ":" << endl;
      for(vector< vector<double> >::size_type j = 0; j != _twgts[k].size(); j++) {
        for(vector<double>::size_type i = 0; i != _twgts[k][j].size(); i++) {
          cout << _twgts[k][j][i] << " ";
        }
        cout << endl;
      }
    }
    const string morestars2(46+label2.size(),'*');
    cout << morestars2 << endl;

    cout.flags(origflags);
    cout.precision(prec);
    cout.width(width);
  }
}

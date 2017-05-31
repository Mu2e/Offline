#ifndef MVATools_HH
#define MVATools_HH

// framework
#include "fhiclcpp/ParameterSet.h"

// Xerces XML Parser
#include <xercesc/dom/DOM.hpp>

#include <vector>
#include <string>

namespace mu2e 
{
  class MVATools
  {
      friend class FMVATools;     
    
  public:
    explicit MVATools(fhicl::ParameterSet const&);
    virtual ~MVATools();
    void initMVA();
    double evalMVA(std::vector<double> const&);
    void showMVA()const;
  protected:
    void getNorm();
    void getWgts();
  private:
    std::string _mvaWgtsFile;
    xercesc::DOMDocument* _xmlDoc;
    std::vector <std::vector<double> >_wn;
    std::vector<double>_wnr2;
    std::vector <std::vector<std::vector<double> > > _wgts;
    std::vector <std::vector<std::vector<double> > > _twgts;
    std::vector<double>_x;
    std::vector<double>_y;
    std::vector<std::string> _title, _label;
  };
}
#endif

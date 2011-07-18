//
// Class which holds (and is able to display) information of objects displayed by the event display. It is used as one of the base classes of each shape, e.g. TPolyLine3DTrack, etc.
//
// $Id: ComponentInfo.h,v 1.10 2011/07/18 21:51:04 greenc Exp $
// $Author: greenc $
// $Date: 2011/07/18 21:51:04 $
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_dict_classes_ComponentInfo_h
#define EventDisplay_src_dict_classes_ComponentInfo_h


#include <TText.h>
#include <TPad.h>
#include <iostream>
#include <sstream>
#include <vector>
#ifndef __CINT__
#include "boost/shared_ptr.hpp"
#include "cpp0x/type_traits"
#endif

namespace mu2e_eventdisplay
{

  class ComponentInfo
  {
#ifndef __CINT__
    public:
#endif
    ComponentInfo(const ComponentInfo &);
    ComponentInfo& operator=(const ComponentInfo &);

    private:
#ifndef __CINT__
    std::vector<boost::shared_ptr<TText> > _text;
    boost::shared_ptr<std::string>         _name;
    //will be expanded later for more advanced info, e.g. histograms, etc.
    //this will be only a base class, and the specifics will be in inherited classes

    public:
    ComponentInfo()
    {
      _name=boost::shared_ptr<std::string>(new std::string);
      for(int i=0; i<5; i++)
      {
        boost::shared_ptr<TText> newLine(new TText(0.0,0.0,NULL));
        newLine->SetTextColor(1);
        newLine->SetTextSizePixels(60);
        _text.push_back(newLine);
      }
    }


    ComponentInfo(const boost::shared_ptr<ComponentInfo> c)  //no reference, since shared_ptr
    {
      _name=c->getName();
      std::vector<boost::shared_ptr<TText> >::const_iterator iter;
      for(iter=c->getText().begin(); iter!=c->getText().end(); iter++)
      {
        _text.push_back(*iter);
      }
    }

    virtual ~ComponentInfo() {} //TTexts will get deleted automatically.

    const boost::shared_ptr<std::string> getName() const {return _name;}

    void setName(const char *newName) {_name->assign(newName);}

    const std::vector<boost::shared_ptr<TText> > &getText() const {return _text;}

    void setText(const unsigned int lineNumber, const char *newText)
    {
      if(lineNumber<0 || lineNumber>=_text.size()) return; //TODO throw exception
      _text[lineNumber]->SetTitle(newText);
    }

    template <typename NUMERIC_TYPE>
      // This template will not instantiate for non-arithmetic types (cf SFINAE).
      // std::enable_if expression below will give return type void if instantiated.
      typename std::enable_if<std::is_arithmetic<NUMERIC_TYPE>::value>::type 
      // const on pass-by-value is meaningless in an argument list.
      expandLine(unsigned int lineNumber, NUMERIC_TYPE newNumber, const char *unit)
    {
      if(lineNumber<0 || lineNumber>=_text.size()) return;  //TODO throw exception
      const char *oldLine = _text[lineNumber]->GetTitle();
      if(oldLine!=NULL)
      {
        std::ostringstream newLine;
        // No manipulators is equivalent to %g for floating point types, %i for integral.
        // Add eg std::setfill('0') and/or std::setw(n) from <iomanip> for padding.
        newLine<<oldLine<<", "<<newNumber<<unit;
        _text[lineNumber]->SetTitle(newLine.str().c_str());
      }
    }

    void removeLine(const unsigned int lineNumber)
    {
      if(lineNumber<0 || lineNumber>=_text.size()) return;  //TODO throw exception
      _text[lineNumber]->SetTitle("");
    }


    void getExpectedSize(unsigned int &width, unsigned int &height) const
    {
      width=0; height=0;
      int i=1;
      std::vector<boost::shared_ptr<TText> >::const_iterator iter;
      for(iter=_text.begin(); iter!=_text.end(); iter++, i++)
      {
        if((*iter)->GetTitle()!=NULL) 
        {
          if(strlen((*iter)->GetTitle())>0)
          {
            unsigned int w,h;
            (*iter)->GetBoundingBox(w,h);
            if(w+20>width) width=w+20;
            height=i*20;
          }
        }
      }
    }

    void showInfo(const unsigned int width, const unsigned int height) const
    {
      int i=1;
      std::vector<boost::shared_ptr<TText> >::const_iterator iter;
      for(iter=_text.begin(); iter!=_text.end(); iter++, i++)
      {
        if((*iter)->GetTitle()!=NULL) 
        {
          double x=20.0/width;
          double y=1.0-i*20.0/height;
          (*iter)->SetX(x);
          (*iter)->SetY(y);
          (*iter)->Draw();
        }
      }
    }

#endif
    ClassDef(ComponentInfo,0);
  };

}
#endif /* EventDisplay_src_dict_classes_ComponentInfo_h */

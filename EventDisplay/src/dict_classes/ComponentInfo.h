//
// Class which holds (and is able to display) information of objects displayed by the event display. It is used as one of the base classes of each shape, e.g. TPolyLine3DTrack, etc. 
//
// $Id: ComponentInfo.h,v 1.4 2011/03/05 05:06:09 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/03/05 05:06:09 $
//
// Original author Ralf Ehrlich
//

#ifndef COMPONENTINFO_H
#define COMPONENTINFO_H

#include <TText.h>
#include <iostream>
#include <vector>
#ifndef __CINT__
#include "boost/shared_ptr.hpp"
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
        boost::shared_ptr<TText> newLine(new TText(0.02,0.9-0.1*i,NULL));
        newLine->SetTextColor(1);
        newLine->SetTextSize(0.05);
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

    void expandLine(const unsigned int lineNumber, const double newNumber, const char *unit)
    {
      if(lineNumber<0 || lineNumber>=_text.size()) return;  //TODO throw exception
      const char *oldLine = _text[lineNumber]->GetTitle();
      if(oldLine!=NULL)
      {
        char newLine[strlen(oldLine)+30];
        sprintf(newLine,"%s, %g%s",oldLine, newNumber, unit);
        _text[lineNumber]->SetTitle(newLine);
      }
    }

    void expandLine(const unsigned int lineNumber, const int newNumber, const char *unit)
    {
      if(lineNumber<0 || lineNumber>=_text.size()) return;  //TODO throw exception
      const char *oldLine = _text[lineNumber]->GetTitle();
      if(oldLine!=NULL)
      {
        char newLine[strlen(oldLine)+30];
        sprintf(newLine,"%s, %i%s",oldLine, newNumber, unit);
        _text[lineNumber]->SetTitle(newLine);
      }
    }

    void removeLine(const unsigned int lineNumber)
    {
      if(lineNumber<0 || lineNumber>=_text.size()) return;  //TODO throw exception
      _text[lineNumber]->SetTitle("");
    }


    void showInfo() const
    {
      std::vector<boost::shared_ptr<TText> >::const_iterator iter;
      for(iter=_text.begin(); iter!=_text.end(); iter++)
      {
        if((*iter)->GetTitle()!=NULL) (*iter)->Draw();
      }
    }

#endif
    ClassDef(ComponentInfo,0);
  };

}
#endif

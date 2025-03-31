//
// Class which holds (and is able to display) information of objects displayed by the event display. It is used as one of the base classes of each shape, e.g. TPolyLine3DTrack, etc.
//
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_dict_classes_ComponentInfo_h
#define EventDisplay_src_dict_classes_ComponentInfo_h


#include <TStyle.h>
#include <TText.h>
#include <TAxis.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TList.h>
#include <TF1.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
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

    protected:
#ifndef __CINT__
    std::vector<boost::shared_ptr<TText> >    _text;
    std::vector<boost::shared_ptr<TObject> >  _hist;
    boost::shared_ptr<std::string>            _name;

    public:
    ComponentInfo()
    {
      _name=boost::shared_ptr<std::string>(new std::string);
      for(int i=0; i<7; i++)
      {
        char const* tmp = nullptr;
        boost::shared_ptr<TText> newLine(new TText(0.0,0.0,tmp));
        newLine->SetTextColor(1);
        newLine->SetTextSize(0.07);
        _text.push_back(newLine);
      }
    }

    virtual ~ComponentInfo() {} //TTexts will get deleted automatically.

    const boost::shared_ptr<std::string> getName() const {return _name;}

    void setName(const char *newName) {_name->assign(newName);}

    std::vector<boost::shared_ptr<TObject> > &getHistVector() {return _hist;}

    const std::vector<boost::shared_ptr<TText> > &getText() const {return _text;}

    void setText(const unsigned int lineNumber, const char *newText)
    {
      if(lineNumber>=_text.size()) return; //TODO throw exception
      _text[lineNumber]->SetTitle(newText);
    }

    void expandLine(unsigned int lineNumber, const char *text)
    {
      if(lineNumber>=_text.size()) return;  //TODO throw exception
      const char *oldLine = _text[lineNumber]->GetTitle();
      if(oldLine!=nullptr)
      {
        std::ostringstream newLine;
        newLine<<oldLine<<"  "<<text;
        _text[lineNumber]->SetTitle(newLine.str().c_str());
      }
    }

    void removeLine(const unsigned int lineNumber)
    {
      if(lineNumber>=_text.size()) return;  //TODO throw exception
      _text[lineNumber]->SetTitle("");
    }


    void getExpectedSize(unsigned int &width, unsigned int &height) const
    {
      width=0; height=0;
      int i=1;
      std::vector<boost::shared_ptr<TText> >::const_iterator iter;
      for(iter=_text.begin(); iter!=_text.end(); iter++, i++)
      {
        if((*iter)->GetTitle()!=nullptr)
        {
          if(strlen((*iter)->GetTitle())>0)
          {
            unsigned int w,h;
            (*iter)->GetBoundingBox(w,h);
            w*=1.5;   //don't know why this is necessary
            if(w+40>width) width=w+40;
            height+=h+5;
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
        if((*iter)->GetTitle()!=nullptr)
        {
          double x=20.0/width;
          double y=1.0-i*20.0/height;
          (*iter)->SetX(x);
          (*iter)->SetY(y);
          (*iter)->Draw();
        }
      }
    }

    void showHist(unsigned int i) const
    {
      if(i<_hist.size())
      {
        gStyle->SetOptTitle(0);
        TGraph      *g = dynamic_cast<TGraph*>(_hist[i].get());
        TMultiGraph *m = dynamic_cast<TMultiGraph*>(_hist[i].get());
        TH1         *h = dynamic_cast<TH1*>(_hist[i].get());
        if(g)
        {
          g->GetXaxis()->SetLabelSize(0.05);
          g->GetXaxis()->SetTitleSize(0.05);
          g->GetXaxis()->SetTitleOffset(0.8);
          g->GetYaxis()->SetLabelSize(0.05);
          g->GetYaxis()->SetTitleSize(0.05);
          g->GetYaxis()->SetTitleOffset(0.8);
          g->Draw("ap");
        }
        else if(h)
        {                                    //TODO may need other settings
          h->GetXaxis()->SetLabelSize(0.05);
          h->GetXaxis()->SetTitleSize(0.05);
          h->GetXaxis()->SetTitleOffset(0.8);
          h->GetYaxis()->SetLabelSize(0.05);
          h->GetYaxis()->SetTitleSize(0.05);
          h->GetYaxis()->SetTitleOffset(0.8);
          h->Draw("ap");
        }
        else if(m)
        {
          m->GetXaxis()->SetTitle("t [ns]");
          m->GetYaxis()->SetTitle("ADC");
          m->GetXaxis()->SetLabelSize(0.05);
          m->GetXaxis()->SetTitleSize(0.05);
          m->GetXaxis()->SetTitleOffset(0.8);
          m->GetYaxis()->SetLabelSize(0.05);
          m->GetYaxis()->SetTitleSize(0.05);
          m->GetYaxis()->SetTitleOffset(0.8);
          m->Draw("a");

          TList *functionList = m->GetListOfFunctions();
          if(functionList->GetSize()>0)
          {
            //1st function is the pedestal that gets added to each reco pulse fitted function
            TF1 *fpedestal = dynamic_cast<TF1*>(functionList->At(0));
            for(int iFunction=1; iFunction<functionList->GetSize(); ++iFunction)
            {
              TF1 *f = dynamic_cast<TF1*>(functionList->At(iFunction));
              TString functionString;
              functionString.Append(fpedestal->GetExpFormula());
              if(iFunction>0) functionString.Append("+");
              functionString.Append(f->GetExpFormula());

              TF1 *functionSum = new TF1("recoFunctionSum",functionString.Data());
              functionSum->DrawF1(m->GetXaxis()->GetXmin(),m->GetXaxis()->GetXmax(),"csame");
              functionSum->SetLineWidth(2);
              functionSum->SetLineColor(2);
              if(f->GetLineStyle()==2) functionSum->SetLineStyle(2);
            }
          }

          const std::string multigraphName = m->GetName();
          if(multigraphName.compare(0,8,"Waveform")==0)
          {
            TText *t[3];
            t[0] = new TText(.12,.95,_name->c_str());
            t[1] = new TText(.12,.90,multigraphName.c_str());
            t[2] = new TText(.12,.85,"No Reco pulses");
            int sipm=atoi(&multigraphName.back());
            for(size_t itext=0; itext<_text.size(); itext++)
            {
              const std::string s = _text[itext]->GetTitle();
              size_t spos = s.find(Form("SiPM%i",sipm));
              if(spos!=std::string::npos && spos+5<s.size())
              {
                t[2]->SetTitle(Form("Reco pulse(s) %s",s.substr(spos+5).c_str()));
                break;
              }
            }
            for(int j=0; j<3; j++)
            {
              t[j]->SetNDC();
              t[j]->SetTextSize(0.05);
              t[j]->SetTextColor(1);
              t[j]->Draw();
            }
          }
        }
        else _hist[i]->Draw();
      }
    }
#endif
    ClassDef(ComponentInfo,0);
  };

}
#endif /* EventDisplay_src_dict_classes_ComponentInfo_h */

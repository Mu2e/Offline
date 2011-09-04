//
// Class which manages a root file which stores a TTree with one branch which hold TObjArrays. These TObjArrays hold all TPolyLines, TPolyLine3Ds, and TTexts of each event.
//
// $Id: RootFileManager.h,v 1.1 2011/09/04 04:43:34 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2011/09/04 04:43:34 $
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_RootFileManager_h
#define EventDisplay_src_RootFileManager_h

#include <TFile.h>
#include <TTree.h>
#include <TPad.h>
#include <TPolyLine3D.h>
#include <TPolyLine.h>
#include <TBox.h>
#include <TText.h>

namespace mu2e_eventdisplay
{

class RootFileManager
{
  bool                         _active;
  TDirectory*                  _directory;
  TObjArray*                   _objArray; 
  boost::shared_ptr<TFile>     _file;
  boost::shared_ptr<TTree>     _tree;

  RootFileManager(const RootFileManager &);
  RootFileManager& operator=(const RootFileManager &);

  public:

  RootFileManager() : _active(false) {}

  void setFile(const char *c)
  {
    if(_active) return;
    _active=true;
    TDirectory *tmpDirectory=gDirectory;
    _file=boost::shared_ptr<TFile>(new TFile(c,"RECREATE"));
    _directory=gDirectory;
    
    _objArray = NULL;
    _tree = boost::shared_ptr<TTree>(new TTree("Tree","Tree"));
    _tree->Branch("Branch","TObjArray",&_objArray);

    gDirectory=tmpDirectory;
  }

  bool isActive() {return _active;}

  void storeEvent()
  {
    if(_active)
    {
      _objArray = new TObjArray();
      _objArray->SetOwner();
      TObjLink *lnk = gPad->GetListOfPrimitives()->FirstLink();
      while(lnk)
      {
        TObject *obj = lnk->GetObject();
        if(obj->InheritsFrom(TPolyLine3D::Class()))
        {
          TPolyLine3D *p = new TPolyLine3D(*dynamic_cast<TPolyLine3D*>(obj));
          if(p) _objArray->Add(p);
        }
        if(obj->InheritsFrom(TPolyLine::Class()))
        {
          TPolyLine *p = new TPolyLine(*dynamic_cast<TPolyLine*>(obj));
          if(p) _objArray->Add(p);
        }
        if(obj->InheritsFrom(TBox::Class()))
        {
          TBox *p = new TBox(*dynamic_cast<TBox*>(obj));
          if(p) _objArray->Add(p);
        }
        if(obj->InheritsFrom(TText::Class()))
        {
          TText *p = new TText(*dynamic_cast<TText*>(obj));
          if(p) _objArray->Add(p);
        }
        lnk = lnk->Next();
      }
      _tree->Fill();
      _objArray->Clear();
      delete _objArray;
    }
  }

  void write()
  {
    if(_active)
    {
      TDirectory *tmpDirectory=gDirectory;
      gDirectory=_directory;
      _tree->Write();
      _file->Close();
      gDirectory=tmpDirectory;
    }
  }

};

}

#endif /* EventDisplay_src_RootFileManager_h */

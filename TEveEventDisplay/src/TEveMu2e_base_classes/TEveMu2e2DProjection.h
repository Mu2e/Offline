#ifndef TEveMu2e2DProjection_h
#define TEveMu2e2DProjection_h

#include <TObject.h>
#include <TEveProjectionManager.h>
#include <TEveViewer.h>
#include <TEveScene.h>

namespace mu2e {

  class TEveMu2e2DProjection: public TEveProjectionManager { 
    public:
      #ifndef __CINT__
      TEveMu2e2DProjection(){};
      virtual ~TEveMu2e2DProjection(){};
      #endif
      TEveViewer *fXYView;
      TEveViewer *fRZView;
      TEveProjectionManager *fXYMgr;
      TEveProjectionManager *fRZMgr;
      TEveScene *fDetXYScene;
      TEveScene *fDetRZScene;
      TEveScene *fEvtXYScene;
      TEveScene *fEvtRZScene;
      ClassDef(TEveMu2e2DProjection, 0);
  };
}
#endif




//
// Class which manages the "Save As" dialog boxes by providing the right file types, and checking the extension of the returned file name.
//
// $Id: SaveDialogManager.h,v 1.4 2012/02/17 21:43:05 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2012/02/17 21:43:05 $
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_SaveDialogManager_h
#define EventDisplay_src_SaveDialogManager_h

#include <TGFileDialog.h>
#include <string>

class SaveDialogManager
{
  public:
  static bool singleImage(std::string &f) //no Root files (they need to be stored in TTrees, which allow multiple "images")
  {
    const char *fileTypes[]={"GIF files","*.gif",
                             "PNG files","*.png",
//                             "JPEG files","*.jpg",  //images have a very low quality
//                             "TIFF files","*.tiff", //crashes
                             "PDF files","*.pdf",
                             "PS files","*.ps",
                             "EPS files","*.eps",
                             "XPM files","*.xpm",
                             "SVG files","*.svg",
//                             "XML files","*.xpm",   //requires a streamer for the ComponentInfoContainer
//                             "C files","*.C",       //opening the files requires loading all mu2e libraries
//                             "C++ files","*.cxx",   //opening the files requires loading all mu2e libraries
                             0,0};
    bool to_return=dialogBox(f,fileTypes);
    return(to_return);
  }

  static bool rootTree(std::string &f) //only Root TTrees
  {
    const char *fileTypes[]={"ROOT files","*.root",
                             0,0};
    bool to_return=dialogBox(f,fileTypes);
    return(to_return);
  }

  static bool animatedImage(std::string &f, bool &isRootFile) //only gif and Root TTrees 
                                                               //(which store multiple frames of one event)
  {
    const char *fileTypes[]={"GIF files","*.gif",
                             "ROOT files","*.root",
                             0,0};
    bool to_return=dialogBox(f,fileTypes);
    if(f.length()>=5)
    {
      if(f.compare(f.length()-5, 5, ".root")) isRootFile=false; else isRootFile=true;
    }
    return(to_return);
  }

  private:
  static bool dialogBox(std::string &f, const char **fileTypes)
  {
    TGFileInfo fileInfo;
    fileInfo.fFileTypes = fileTypes;
    new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDSave, &fileInfo); //ROOT takes care of deleting this
    if(!fileInfo.fFilename) return(false);

    f.assign(fileInfo.fFilename);
    const char *fileExtension=fileTypes[fileInfo.fFileTypeIdx+1];
    int  fileExtensionLength=strlen(fileExtension)-1;
    if(f.compare(f.length()-fileExtensionLength, fileExtensionLength, fileExtension+1))
    {
      f.append(fileExtension+1);
    }
    return(true);
  }
};

#endif

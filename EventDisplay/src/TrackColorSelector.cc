using namespace std;
#include "EventDisplay/src/TrackColorSelector.h"

namespace mu2e_eventdisplay
{

void TrackColorSelector::setupTrackLegend()
{
  for(int i=0; i<30; i++)
  {
    _legendParticleLine[i]=new TPolyLine();
    _legendParticleLine[i]->SetPoint(0, 0.6,0.45-i*0.05);
    _legendParticleLine[i]->SetPoint(1, 0.7,0.45-i*0.05);
    _legendParticleGroup[i]=new TText(0.6,0.44-i*0.05,"");
    _legendParticleGroup[i]->SetTextColor(_whiteBackground?kBlack:kWhite);
    _legendParticleGroup[i]->SetTextSize(0.025);
    _legendParticleText[i]=new TText(0.72,0.44-i*0.05,"");
    _legendParticleText[i]->SetTextSize(0.025);
  }
}

void TrackColorSelector::drawGroupName(int i, const std::string &name)
{
  _legendParticleGroup[i]->SetTitle(name.c_str());
  _legendParticleGroup[i]->Draw();
}

void TrackColorSelector::drawLineAndText(int i, int color, const std::string &name)
{
  _legendParticleText[i]->SetTitle(name.c_str());
  _legendParticleText[i]->SetTextColor(color);
  _legendParticleText[i]->Draw();
  _legendParticleLine[i]->SetLineColor(color);
  _legendParticleLine[i]->Draw();
}

void TrackColorSelector::drawTrackLegend(TText **legendParticleGroup, TText **legendParticleText, TPolyLine **legendParticleLine)
{
  _legendParticleGroup=legendParticleGroup;
  _legendParticleText=legendParticleText;
  _legendParticleLine=legendParticleLine;

  setupTrackLegend();
  unsigned int n=_trackInfos->size();
  switch(n)
  {
      case 1: drawGroupName(0,_trackInfos->at(0).entryText);
              drawLineAndText(1, 2, "e+, e-");
              drawLineAndText(2, 3, "mu+, mu-");
              drawLineAndText(3, 4, "gamma");
              drawLineAndText(4, 5, "n0");
              drawLineAndText(5, 6, "neutrinos");
              drawLineAndText(6, 28, "other particles");
              break;
      case 2: drawGroupName(0,_trackInfos->at(0).entryText);
              drawLineAndText(1, 2, "e+, e-");
              drawLineAndText(2, 3, "mu+, mu-");
              drawLineAndText(3, 4, "gamma");
              drawLineAndText(4, 5, "n0");
              drawLineAndText(5, 28, "other particles");
              drawGroupName(6,_trackInfos->at(1).entryText);
              drawLineAndText(7, 7, "e+, e-");
              drawLineAndText(8, 6, "mu+, mu-");
              drawLineAndText(9, 8, "gamma");
              drawLineAndText(10, 9, "n0");
              drawLineAndText(11, 46, "other particles");
              break;
      case 3: drawGroupName(0,_trackInfos->at(0).entryText);
              drawLineAndText(1, 2, "e+, e-");
              drawLineAndText(2, 3, "mu+, mu-");
              drawLineAndText(3, 28, "other particles");
              drawGroupName(4,_trackInfos->at(1).entryText);
              drawLineAndText(5, 4, "e+, e-");
              drawLineAndText(6, 5, "mu+, mu-");
              drawLineAndText(7, 46, "other particles");
              drawGroupName(8,_trackInfos->at(2).entryText);
              drawLineAndText(9, 7, "e+, e-");
              drawLineAndText(10, 6, "mu+, mu-");
              drawLineAndText(11, 9, "other particles");
              break;
      case 4: drawGroupName(0,_trackInfos->at(0).entryText);
              drawLineAndText(1, 2, "e+, e-");
              drawLineAndText(2, 28, "other particles");
              drawGroupName(3,_trackInfos->at(1).entryText);
              drawLineAndText(4, 5, "e+, e-");
              drawLineAndText(5, 46, "other particles");
              drawGroupName(6,_trackInfos->at(2).entryText);
              drawLineAndText(7, 3, "e+, e-");
              drawLineAndText(8, 8, "other particles");
              drawGroupName(9,_trackInfos->at(3).entryText);
              drawLineAndText(10, 4, "e+, e-");
              drawLineAndText(11, 9, "other particles");
              break;
     default: for(unsigned int i=0; i<8 && i<n; i++)
              { 
                drawLineAndText(i, i+2, _trackInfos->at(i).entryText);
              }
              if(n<8) drawLineAndText(n, 28, "other tracks");
              else drawLineAndText(8, 28, "other tracks");
  };
}

int TrackColorSelector::getColor(boost::shared_ptr<Track> track)
{
  int particleid=track->getParticleId();
  int trackclass=track->getTrackClass(); 
  int trackclassindex=track->getTrackClassIndex(); 
  int color=_whiteBackground?28:25;

  unsigned int n=_trackInfos->size();
  switch(n)
  {
      case 1: switch(particleid)
              {
                case   11:
                case  -11: color=2; break;   //e+,e-
                case   13:
                case  -13: color=3; break;   //mu+,mu-
                case   22: color=4; break;   //gamma
                case 2112: color=5; break;   //n0
                case   12:
                case  -12:
                case   14:
                case  -14:
                case   16:
                case  -16: color=6; break;   //neutrinos
                default  : color=_whiteBackground?28:25;
              };
              break;
      case 2: if((trackclass==_trackInfos->at(0).classID) && 
                 (trackclassindex==_trackInfos->at(0).index))
              {
                switch(particleid)
                {
                  case   11:
                  case  -11: color=2; break;   //e+,e-
                  case   13:
                  case  -13: color=3; break;   //mu+,mu-
                  case   22: color=4; break;   //gamma
                  case 2112: color=5; break;   //n0
                  default  : color=_whiteBackground?28:25;
                };
              }
              if((trackclass==_trackInfos->at(1).classID) && 
                 (trackclassindex==_trackInfos->at(1).index))
              {
                switch(particleid)
                {
                  case   11:
                  case  -11: color=7; break;   //e+,e-
                  case   13:
                  case  -13: color=6; break;   //mu+,mu-
                  case   22: color=8; break;  //gamma
                  case 2112: color=9; break;  //n0
                  default  : color=46;
                };
              }
              break;
      case 3: if((trackclass==_trackInfos->at(0).classID) && 
                 (trackclassindex==_trackInfos->at(0).index))
              {
                switch(particleid)
                {
                  case   11:
                  case  -11: color=2; break;   //e+,e-
                  case   13:
                  case  -13: color=3; break;   //mu+,mu-
                  default  : color=_whiteBackground?28:25;
                };
              }
              if((trackclass==_trackInfos->at(1).classID) && 
                 (trackclassindex==_trackInfos->at(1).index))
              {
                switch(particleid)
                {
                  case   11:
                  case  -11: color=4; break;   //e+,e-
                  case   13:
                  case  -13: color=5; break;   //mu+,mu-
                  default  : color=46;
                };
              }
              if((trackclass==_trackInfos->at(2).classID) && 
                 (trackclassindex==_trackInfos->at(2).index))
              {
                switch(particleid)
                {
                  case   11:
                  case  -11: color=7; break;   //e+,e-
                  case   13:
                  case  -13: color=6; break;   //mu+,mu-
                  default  : color=9;
                };
              }
              break;
      case 4: if((trackclass==_trackInfos->at(0).classID) && 
                 (trackclassindex==_trackInfos->at(0).index))
              {
                switch(particleid)
                {
                  case   11:
                  case  -11: color=2; break;   //e+,e-
                  default  : color=_whiteBackground?28:25;
                };
              }
              if((trackclass==_trackInfos->at(1).classID) && 
                 (trackclassindex==_trackInfos->at(1).index))
              {
                switch(particleid)
                {
                  case   11:
                  case  -11: color=5; break;   //e+,e-
                  default  : color=46;
                };
              }
              if((trackclass==_trackInfos->at(2).classID) && 
                 (trackclassindex==_trackInfos->at(2).index))
              {
                switch(particleid)
                {
                  case   11:
                  case  -11: color=3; break;   //e+,e-
                  default  : color=8;
                };
              }
              if((trackclass==_trackInfos->at(3).classID) && 
                 (trackclassindex==_trackInfos->at(3).index))
              {
                switch(particleid)
                {
                  case   11:
                  case  -11: color=4; break;   //e+,e-
                  default  : color=9;
                };
              }
              break;
     default: color=_whiteBackground?28:25;
              for(unsigned int i=0; i<8 && i<n; i++)
              {
                if((trackclass==_trackInfos->at(i).classID) && 
                   (trackclassindex==_trackInfos->at(i).index)) color=i+2;
              }
  };
  return color;
}

}

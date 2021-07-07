//
// Class assigns colors to tracks, depending on the number of track groups selected.
//
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_TrackColorSelector_h
#define EventDisplay_src_TrackColorSelector_h

#include "Offline/EventDisplay/src/ContentSelector.h"
#include "Offline/EventDisplay/src/Track.h"
#include <TText.h>
#include <TPolyLine.h>

namespace mu2e_eventdisplay
{

class TrackColorSelector
{
  TrackColorSelector();
  TrackColorSelector(const TrackColorSelector &);
  TrackColorSelector& operator=(const TrackColorSelector &);

  std::vector<ContentSelector::trackInfoStruct> const *_trackInfos;
  TText **_legendParticleGroup;
  TText **_legendParticleText;
  TPolyLine **_legendParticleLine;
  bool _whiteBackground;

  void setupTrackLegend();
  void drawGroupName(int i, const std::string &name);
  void drawLineAndText(int i, int color, const std::string &name);

  public:
  TrackColorSelector(std::vector<ContentSelector::trackInfoStruct> const *trackInfos, bool whiteBackground) : _trackInfos(trackInfos), _whiteBackground(whiteBackground) {}
  void drawTrackLegend(TText **legendParticleGroup, TText **legendParticleText, TPolyLine **legendParticleLine);
  int  getColor(boost::shared_ptr<Track> track);
};

}
#endif /* EventDisplay_src_TrackColorSelector_h */


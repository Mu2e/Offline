//
// Class assigns colors to tracks, depending on the number of track groups selected.
//
// $Id: TrackColorSelector.h,v 1.1 2011/08/14 06:31:19 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2011/08/14 06:31:19 $
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_TrackColorSelector_h
#define EventDisplay_src_TrackColorSelector_h

#include "ContentSelector.h"
#include "Track.h"
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

  void setupTrackLegend();
  void drawGroupName(int i, const std::string &name);
  void drawLineAndText(int i, int color, const std::string &name);

  public:
  TrackColorSelector(std::vector<ContentSelector::trackInfoStruct> const *trackInfos) : _trackInfos(trackInfos) {}
  void drawTrackLegend(TText **legendParticleGroup, TText **legendParticleText, TPolyLine **legendParticleLine);
  int  getColor(boost::shared_ptr<Track> track);
};

}
#endif /* EventDisplay_src_TrackColorSelector_h */


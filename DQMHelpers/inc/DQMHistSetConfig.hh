#ifndef DQMHelpers_inc_DQMHistSetConfig_hh
#define DQMHelpers_inc_DQMHistSetConfig_hh
// Validated FHiCL for DQMHistSet::Config: which copies of each histogram a job
// keeps, and which it publishes. Offline modules nest it with fhicl::Table;
// online modules build the same table from their ParameterSet, so the two
// cannot drift, and fhiclcpp rejects an unknown or mistyped key in both.
//
//   hists : {
//     annotateTitles : true      # append the range tag to every copy's title
//     stampMetadata  : true      # add the dqmSegment TNamed to every copy
//     liveSeries     : false     # book the client's graphs (online only)
//     rules : [
//       { match    : "h1_channels"        # glob over the dir-qualified path
//         modes    : [ "job", "window" ]
//         window   : { span : 50000  unit : "ewt"  subdivisions : 10 }
//         liveName : "h1_channelsLastEwt"
//         publish  : true                 # hand this hist's copies to the consumer
//         group    : "occupancy" }        # omit for one entry per copy
//     ]
//   }
//
// Binning is not here and never will be; see DQMAxis.hh. Nor is there a way to
// drop a histogram: FHiCL chooses copies, not the histogram set.
//
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/DQMHistSet.hh"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include <string>
#include <vector>

namespace mu2e {

struct DQMSubRunFhicl {
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;
  fhicl::Atom<int> keep{Name("keep"), Comment("archived subruns to keep; -1 = all"), 0};
  fhicl::Atom<bool> persist{Name("persist"), Comment("write the archived copies"), true};
  fhicl::Atom<bool> persistLive{
      Name("persistLive"), Comment("write the in-progress copy"), false};
};

struct DQMWindowFhicl {
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;
  fhicl::Atom<unsigned long long> span{
      Name("span"), Comment("window length in `unit`s"), 50000ull};
  fhicl::Atom<std::string> unit{
      Name("unit"), Comment("event | ewt | subrun"), "event"};
  fhicl::Atom<int> subdivisions{
      Name("subdivisions"), Comment("ring depth; 1 = disjoint blocks"), 10};
  fhicl::Atom<int> keep{Name("keep"), Comment("archived spans to keep"), 0};
  fhicl::Atom<bool> persist{Name("persist"), Comment("write the archived copies"), false};
  fhicl::Atom<bool> persistLive{
      Name("persistLive"), Comment("write the rolling copy"), false};
};

struct DQMRuleFhicl {
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;
  fhicl::Atom<std::string> match{
      Name("match"), Comment("glob over the directory-qualified histogram path"), "*"};
  fhicl::Sequence<std::string> modes{
      Name("modes"), Comment("any of job, subrun, window"),
      std::vector<std::string>{"job"}};
  fhicl::Atom<bool> jobPersist{
      Name("jobPersist"), Comment("write the job copy"), true};
  fhicl::OptionalTable<DQMSubRunFhicl> subrun{Name("subrun"), Comment("subrun copies")};
  fhicl::OptionalTable<DQMWindowFhicl> window{Name("window"), Comment("rolling window")};
  fhicl::Atom<std::string> liveName{
      Name("liveName"), Comment("name of the rolling copy"), ""};
  fhicl::Atom<bool> publish{
      Name("publish"), Comment("hand this histogram's copies to the consumer"), false};
  fhicl::Atom<std::string> group{
      Name("group"), Comment("label to collect published copies under"), ""};
  fhicl::Atom<std::string> archiveGroup{
      Name("archiveGroup"), Comment("label for archived copies; empty = group"), ""};
};

struct DQMHistSetFhicl {
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;
  fhicl::Atom<bool> annotateTitles{
      Name("annotateTitles"), Comment("append the range tag to every copy's title"), true};
  fhicl::Atom<bool> stampMetadata{
      Name("stampMetadata"), Comment("add the dqmSegment stamp to every copy"), true};
  fhicl::Atom<std::string> subrunDir{Name("subrunDir"), Comment("subdirectory"), "bySubrun"};
  fhicl::Atom<std::string> segmentDir{Name("segmentDir"), Comment("subdirectory"), "segments"};
  fhicl::Atom<bool> liveSeries{
      Name("liveSeries"), Comment("book DQMSeries graphs (online only; not hadd-safe)"),
      false};
  //absent means no rule at all, which is the job-only default
  fhicl::OptionalSequence<fhicl::Table<DQMRuleFhicl>> rules{
      Name("rules"), Comment("per-histogram rules; omit for job-only")};
};

// The configuration of a client with no knobs of its own. A client that has
// some declares its own struct with `hists` beside them.
struct DQMClientFhicl {
  fhicl::Table<DQMHistSetFhicl> hists{
      fhicl::Name("hists"),
      fhicl::Comment("copies and publishing; see Offline/DQMHelpers/README.md")};
};

// Throws cet::exception on an unknown mode or window unit.
DQMHistSet::Config toConfig(const DQMHistSetFhicl& c);
// Online: validates a ParameterSet against DQMHistSetFhicl first.
DQMHistSet::Config toConfig(const fhicl::ParameterSet& ps);

} // namespace mu2e

#endif /* DQMHelpers_inc_DQMHistSetConfig_hh */

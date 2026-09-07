#ifndef DQMHelpers_inc_DQMSegmentationConfig_hh
#define DQMHelpers_inc_DQMSegmentationConfig_hh
// FHiCL parsing for DQMSegmentation::Config.
//
// One parser, used by both the offline analyzers (which hand it an
// OptionalDelegatedParameter's ParameterSet) and the otsdaq online modules
// (which hand it ps.get<fhicl::ParameterSet>("segmentation")), so the two
// cannot drift. Unknown keys throw rather than being silently ignored -- a
// mistyped rule that quietly does nothing is worse online than a job that
// refuses to start.
//
// Grammar:
//
//   segmentation : {
//     annotateTitles : true      # append the range tag to every copy's title
//     stampMetadata  : true      # add the dqmSegment TNamed to every copy
//     subrunDir      : "bySubrun"
//     segmentDir     : "segments"
//     rules : [
//       { match   : "h1_channels"        # glob over the dir-qualified path
//         enabled : true                 # false: do not book it at all
//         modes   : [ "job", "window" ]
//         job     : { persist : true }
//         window  : { span : 50000  unit : "ewt"  subdivisions : 10
//                     keep : 4  persist : false  persistLive : false }
//         liveName : "h1_channelsLastEwt"
//         publish  : true            # hand this hist's copies to the consumer
//         group    : "timing_feb"    # collect them under one label; omit for
//       }                            # one entry per copy, keyed on its name
//     ]
//   }
//
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/DQMSegmentation.hh"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/OptionalDelegatedParameter.h"

namespace mu2e {

// Throws cet::exception on an unknown key, an unknown mode, or a bad unit.
DQMSegmentation::Config parseSegmentation(const fhicl::ParameterSet& ps);

// Module-facing overload: an absent block gives the job-only default. The block
// is delegated rather than validated field by field because the rule list is
// variable length, and because one parser shared with the online modules is
// worth more than two schemas that can disagree.
DQMSegmentation::Config parseSegmentation(const fhicl::OptionalDelegatedParameter& p);

// For a module that owns more than one helper. Each helper's registry matches
// rule globs against paths relative to its own directory, so one shared block
// cannot tell two helpers' identically named histograms apart -- `nEvents`
// exists in all three CRV helpers. `specific` wins outright when present
// (it replaces `shared`, it does not merge with it); otherwise `shared`
// applies. Absent both, the job-only default.
DQMSegmentation::Config parseSegmentation(const fhicl::OptionalDelegatedParameter& specific,
                                          const fhicl::OptionalDelegatedParameter& shared);

} // namespace mu2e

#endif /* DQMHelpers_inc_DQMSegmentationConfig_hh */

// Lifecycle shared by every subdetector DQM client.
//
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/DQMClient.hh"

namespace mu2e {

namespace {
DQMHistSet::Config clientConfig(DQMHistSet::Config c)
{
  c.strictBooking = true;
  c.autoNEvents = true;
  return c;
}
} // namespace

DQMClient::DQMClient(const std::string& name, int binningVersion,
                     const DQMHistSet::Config& hists) :
    name_(name), binningVersion_(binningVersion), hists_(clientConfig(hists))
{
  hists_.SetVersion(name_, binningVersion_);
  diag_.SetName(name_);
}

void DQMClient::Book(art::TFileDirectory dir)
{
  if (booked_) {
    return;
  }
  booked_ = true;
  hists_.Book(dir);
  series_.Book(dir, hists_.config().liveSeries);
  book();
  hists_.FreezeBooking();
}

void DQMClient::BeginSubRun(int run, int subrun)
{
  run_ = run;
  subrun_ = subrun;
  hists_.BeginSubRun(run, subrun);
}

void DQMClient::EndSubRun()
{
  endSubRun();
  hists_.EndSubRun();
}

void DQMClient::EndJob()
{
  if (!booked_ || ended_) {
    return;
  }
  ended_ = true;
  endJob();
  hists_.Finalize();
  series_.Persist();
  diag_.Report();
}

void DQMClient::ResetForNewRun()
{
  resetForNewRun();
  hists_.ResetContents();
  series_.ResetContents();
}

void DQMClient::beginEvent(std::optional<uint64_t> clock)
{
  ++nEvents_;
  hists_.Advance(nEvents_, clock);
}

} // namespace mu2e

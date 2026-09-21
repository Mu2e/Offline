// Capped time series for the online monitor.
//
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/DQMSeries.hh"

#include <utility>

namespace mu2e {

DQMSeries::DQMSeries(std::unique_ptr<TGraph> graph, std::size_t maxPoints) :
    graph_(std::move(graph)), maxPoints_(maxPoints)
{}

void DQMSeries::trim()
{
  while (maxPoints_ > 0 && static_cast<std::size_t>(graph_->GetN()) > maxPoints_) {
    graph_->RemovePoint(0);
  }
}

void DQMSeries::Add(double x, double y)
{
  if (graph_ == nullptr) {
    return;
  }
  graph_->SetPoint(graph_->GetN(), x, y);
  trim();
  last_ = y;
}

void DQMSeries::AddIfChanged(double x, double y)
{
  if (last_.has_value() && *last_ == y) {
    return;
  }
  Add(x, y);
}

void DQMSeries::Step(double x, double y)
{
  if (graph_ == nullptr) {
    return;
  }
  const int n = graph_->GetN();
  if (!last_.has_value() || n == 0) {
    Add(x, y);
    return;
  }
  if (*last_ != y) {
    graph_->SetPoint(n, x, *last_);
    Add(x, y);
    return;
  }
  //unchanged: keep the step's start point and move its end to x
  if (n >= 2 && graph_->GetY()[n - 2] == y) {
    graph_->SetPoint(n - 1, x, y);
  } else {
    Add(x, y);
  }
}

void DQMSeries::Clear()
{
  if (graph_ != nullptr) {
    graph_->Set(0);
  }
  last_.reset();
}

void DQMSeriesSet::Book(art::TFileDirectory dir, bool enabled)
{
  dir_ = dir;
  enabled_ = enabled;
}

DQMSeries& DQMSeriesSet::book(const std::string& path, const std::string& title,
                              std::size_t maxPoints)
{
  if (!enabled_ || !dir_) {
    return inert_;
  }
  auto existing = series_.find(path);
  if (existing != series_.end()) {
    return *existing->second;
  }
  const std::size_t slash = path.rfind('/');
  const std::string name = slash == std::string::npos ? path : path.substr(slash + 1);
  auto g = std::make_unique<TGraph>();
  g->SetName(name.c_str());
  g->SetTitle(title.c_str());
  graphs_.push_back(g.get());
  order_.push_back(path);
  auto owner = std::make_unique<DQMSeries>(std::move(g),
                                           maxPoints > 0 ? maxPoints : kDefaultMaxPoints);
  return *series_.emplace(path, std::move(owner)).first->second;
}

DQMSeries& DQMSeriesSet::get(const std::string& path)
{
  auto it = series_.find(path);
  return it == series_.end() ? inert_ : *it->second;
}

void DQMSeriesSet::ResetContents()
{
  for (auto& [path, series] : series_) {
    series->Clear();
  }
}

art::TFileDirectory& DQMSeriesSet::directoryFor(const std::string& dirPath)
{
  if (dirPath.empty()) {
    return *dir_;
  }
  auto it = subdirs_.find(dirPath);
  if (it != subdirs_.end()) {
    return it->second;
  }
  const std::size_t slash = dirPath.rfind('/');
  art::TFileDirectory& parent =
      slash == std::string::npos ? *dir_ : directoryFor(dirPath.substr(0, slash));
  const std::string leaf = slash == std::string::npos ? dirPath : dirPath.substr(slash + 1);
  return subdirs_.emplace(dirPath, parent.mkdir(leaf)).first->second;
}

void DQMSeriesSet::Persist()
{
  if (!dir_ || persisted_) {
    return;
  }
  persisted_ = true;
  for (const std::string& path : order_) {
    const TGraph* g = series_.at(path)->graph();
    const std::size_t slash = path.rfind('/');
    art::TFileDirectory& d =
        directoryFor(slash == std::string::npos ? "" : path.substr(0, slash));
    if (g->GetN() <= 0) {
      d.makeAndRegister<TGraph>(g->GetName(), g->GetTitle());
    } else {
      d.makeAndRegister<TGraph>(g->GetName(), g->GetTitle(), g->GetN(), g->GetX(),
                                g->GetY());
    }
  }
}

} // namespace mu2e

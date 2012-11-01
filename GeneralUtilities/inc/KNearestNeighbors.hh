// Find k nearest neighbors for each point in an input container
// using the given metric.
//
// BEWARE: the current implementation is a naive O(N^2) algorithm.
// There are much faster ways of doing this.
//
// $Id: KNearestNeighbors.hh,v 1.1 2012/11/01 23:32:36 gandr Exp $
// $Author: gandr $
// $Date: 2012/11/01 23:32:36 $
//
// Original author Andrei Gaponenko

#ifndef GeneralUtilities_KNearestNeighbors_hh
#define GeneralUtilities_KNearestNeighbors_hh

#include <queue>
#include <vector>

namespace mu2e {

  template<class Point>
  class KNearestNeighbors {
  public:

    template<class Distance>
    KNearestNeighbors(unsigned k, // number of neighbors to compute
                      const std::vector<Point>& group, // input points
                      const Distance& dist);

    // size() == group.size()
    std::size_t size() const { return pp_.size(); }

    //----------------
    struct Entry {
      Point point;
      double distance;

      Entry(Point p, double d) : point(p), distance(d) {}

      // Larger distances go to the top or the queue and get popped out
      bool operator<(const Entry& b) const {
        return distance < b.distance;
      }
    };

    typedef std::vector<Entry> Points;

    //----------------
    // Primary interface of the class: access neighbors of a point.
    // ipoint is index of the point in the original group container

    //FIXME: const Points& operator[](Points::size_type ipoint) const { return pp_[ipoint]; }
    const Points& operator[](unsigned ipoint) const { return pp_[ipoint]; }

  private:
    std::vector<Points> pp_;
  };

  //----------------------------------------------------------------
  template<class Point> template<class Distance>
  KNearestNeighbors<Point>::KNearestNeighbors(unsigned k,
                                              const std::vector<Point>& group,
                                              const Distance& dist)
    : pp_(group.size())
  {
    typedef std::priority_queue<Entry> Neighbors;
    typedef std::vector<Neighbors> ImplPoints;
    ImplPoints data(group.size());

    // There are better way to compute K nearest neighbors, but
    // brute forcing it is good enough for current Mu2e use case.
    for(unsigned i=0; i<group.size(); ++i) {
      for(unsigned j=i+1; j<group.size(); ++j) {
        const double r = dist(group[i], group[j]);
        data[i].push(Entry(group[j], r));
        if(k < data[i].size()) {
          data[i].pop();
        }
        data[j].push(Entry(group[i], r));
        if(k < data[j].size()) {
          data[j].pop();
        }
      }
    }

    // copy results to the final containers
    for(unsigned i=0; i<data.size(); ++i) {
      while(!data[i].empty()) {
        pp_[i].push_back(data[i].top());
        data[i].pop();
      }
    }
  } // KNearestNeighbors()

  //----------------------------------------------------------------

} // namespace mu2e

#endif /* GeneralUtilities_KNearestNeighbors_hh */

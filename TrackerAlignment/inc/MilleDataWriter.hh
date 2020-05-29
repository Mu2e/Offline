// Ryunosuke O'Neil
// roneil@fnal.gov
// ryunoneil@gmail.com
// http://github.com/ryuwd

// Writes 'mille data' with support for doubles and gzip compression.

#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include <boost/iostreams/categories.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

using namespace boost::iostreams;

template <typename WORDTYPE> class MilleDataWriter {
  static_assert(std::is_same<WORDTYPE, double>::value || std::is_same<WORDTYPE, float>::value,
                "Can only write floats or doubles.");

  std::string filename;
  std::ofstream file_stream;
  filtering_ostream gz_fstream;

  std::vector<WORDTYPE> track_buf;
  std::vector<int> label_buf;

public:
  MilleDataWriter(const std::string& file, bool gzip = true, int buf_size = 1000) :
      filename(file), file_stream(file, std::ios_base::out | std::ios_base::binary) {
    if (gzip) {
      gz_fstream.push(gzip_compressor());
    }
    gz_fstream.push(file_stream);

    track_buf.reserve(buf_size);
    label_buf.reserve(buf_size);
  }

  void pushHit(std::vector<WORDTYPE> const& local_derivatives,
               std::vector<WORDTYPE> const& global_derivatives,
               std::vector<int> const& global_labels, WORDTYPE const& measurement,
               WORDTYPE const& meas_error) {

    if (global_labels.size() != global_labels.size()) {
      std::cerr << "MilleDataWriter: size mismatch between number of global labels and global "
                   "derivatives for this hit"
                << std::endl;
      return;
    }

    if (meas_error < 0) {
      std::cerr << "MilleDataWriter: "
                << "Bad measurement" << std::endl;
      return;
    }
    push(measurement, 0);

    std::vector<int> local_labels(local_derivatives.size());
    std::iota(local_labels.begin(), local_labels.end(), 1);
    push(local_derivatives, local_labels);

    push(meas_error, 0);
    push(global_derivatives, global_labels);
  }

  void flushTrack() {
    const int words = n_words();
    gz_fstream.write(reinterpret_cast<const char*>(&words), sizeof(words));
    std::copy(track_buf.begin(), track_buf.end(), std::ostream_iterator<char>(gz_fstream));
    std::copy(label_buf.begin(), label_buf.end(), std::ostream_iterator<char>(gz_fstream));

    clear();
  }

  void clear() {
    label_buf.clear();
    track_buf.clear();
  }

private:
  void push(std::vector<WORDTYPE> const& data, std::vector<int> const& labels) {
    label_buf.insert(label_buf.end(), labels.begin(), labels.end());
    track_buf.insert(track_buf.end(), data.begin(), data.end());
  }

  void push(WORDTYPE const& data, int const& label) {
    label_buf.emplace_back(label);
    track_buf.emplace_back(data);
  }

  int n_words() { 
    int result = label_buf.size() + track_buf.size();

    // if we are storing doubles, we pass a negative word number
    // to indicate to readC(..) (readc.c) it should read doubles, not floats.
    if (std::is_same<WORDTYPE, double>::value) {
      result = -result;
    }

    return result; 
  }
};

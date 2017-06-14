#ifndef FMVATools_HH
#define FMVATools_HH

// framework
#include "Mu2eUtilities/inc/MVATools.hh"
#include "fhiclcpp/ParameterSet.h"

#include <string>
#include <vector>

namespace mu2e {
using value_t = double;
using values_t = std::vector<value_t>;
using layer_weights_t = std::vector<values_t>;
using layers_weights_t = std::vector<layer_weights_t>;
using inp_params_t = values_t;
using sizes_t = std::vector<size_t>;

constexpr size_t max_synapses = 10;

class FMVATools {
  struct FlatWeights {
    ///[ {layer0,neuron0,weight0},...{layerN,neuronN,weightN}]
    values_t _weights;
    sizes_t _dims;

    void pack(layers_weights_t const& /*layers*/);
    void pack(values_t const& /*weights*/, values_t const& /*norms*/);
  };

 public:
  explicit FMVATools(fhicl::ParameterSet const& pars) : _x(), _y(), _first_layer(), _inner_layers(), _mva_tools(pars) {}

  virtual ~FMVATools() {}

  void initMVA();
  void showMVA() const;

  value_t evalMVA(inp_params_t const& /*params*/) const;

 private:
  std::array<value_t, max_synapses> _x;
  std::array<value_t, max_synapses> _y;
  FlatWeights _first_layer;
  FlatWeights _inner_layers;

  MVATools _mva_tools;
};
}
#endif

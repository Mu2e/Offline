#include "Mu2eUtilities/inc/FMVATools.hh"
using namespace std;
namespace mu2e 
{    
  /*
    auto printer=[](auto const& values, std::string name) {
      std::cout << name << " size=" << values.size() << " values> ";
      for(auto const& value:values ) std::cout << value << ", ";
      std::cout <<endl;
    };
  */
    
  void FMVATools::FlatWeights::pack(layers_weights_t const &layers){
    _dims.reserve(20);
    _weights.reserve(100);
    
    _dims.push_back(layers.size());//number of layers
    
    for(auto const& neurons:layers){
      if(neurons.size()==0)
	throw cet::exception("RECO")<<"mu2e::FMVATools: XML initialization error: Layers must have neuron(s)." << endl;
      
      _dims.push_back(neurons.size()); //number of neurons per layer
      _dims.push_back(neurons[0].size()); //number of weights for each neuron of this layer      
      for(auto const& weights:neurons) {
	if( _dims.back()!=weights.size())
	  throw cet::exception("RECO")<<"mu2e::FMVATools: XML initialization error: All neurons of each layer must have the same number of weights." << endl;
	  
	for(auto const& weight:weights) _weights.push_back(weight);	
      }
    }
  }

  void FMVATools::FlatWeights::pack(values_t const &weights, values_t const& norms){
    if(weights.size()!=norms.size())
      	throw cet::exception("RECO")<<"mu2e::FMVATools: XML initialization error: Invalid first layer params." << endl;

    _dims.reserve(20);
    _weights.reserve(100);
    
    _dims.push_back(0);//number of layers
    _dims.push_back(weights.size());//number of input pairs
    _dims.push_back(2);//number of input pairs

    for (size_t i=0; i< weights.size();++i){
      _weights.push_back(weights[i]);
      _weights.push_back(norms[i]);
    }
  }

  void FMVATools::initMVA(){
    _mva_tools.initMVA();
    
    auto weights=values_t();
    
    for(auto const& a:_mva_tools._wn)
      weights.push_back(a[0]);

    auto norms=values_t();
    
    for(auto const& a:_mva_tools._wnr2)
      norms.push_back(1.0/a);
    
    _first_layer.pack(weights,norms);    
    _inner_layers.pack(_mva_tools._twgts);
    
    _mva_tools.showMVA();
  }

     
  value_t FMVATools::evalMVA(inp_params_t const& v){
    // Normalize
    const auto nweightpairs=_first_layer._dims[1];
    auto s= size_t(0);
   
    for(size_t i = 0; i < nweightpairs; ++i, s+=2){
      _x[i] = -_first_layer._weights[s]*_first_layer._weights[s+1] -1.0 + v[i]*_first_layer._weights[s+1];
    } 

    // do internal layers
    const auto nlayers = _inner_layers._dims[0]-1;    
    s=0;
    
    for(size_t k = 0; k < nlayers; ++k){
      const auto nneurons = _inner_layers._dims[2*k+1];
      const auto nweights = _inner_layers._dims[2*k+2]-1;      
      for(size_t j = 0; j <nneurons; ++j){	
	_y[j]=0.0;	
	for(size_t i = 0; i < nweights; ++i){
	  _y[j] +=  _x[i]*_inner_layers._weights[s++];
        }
         _y[j] +=  _inner_layers._weights[s++];
	 
        _y[j]=1./(1.+exp((float)-_y[j]));//  effect < 1e-6, saves 450 instructions
      }
      _x=_y;
    }

    // do output    
    const auto nweights = _inner_layers._dims[2*nlayers+2]-1;      
        
    value_t mva=0.0;

    for(size_t i = 0; i < nweights; ++i){
      mva +=  _y[i]*_inner_layers._weights[s++];
    }

    mva +=  _inner_layers._weights[s++];

    return mva;
  }  
}

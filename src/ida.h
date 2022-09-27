#ifndef __IDA_H__
#define __IDA_H__

//-------------------------------------------------------------------------------------------------
class InverseDistanceAlgorithm {
public:
  InverseDistanceAlgorithm (Real c, const Node* const* d, int n) : _coeff(c), _nodes(d), _size(n) {}
  //
  inline
  Real calculate (const Point& xyz, const std::vector<Real>& data) const
    {
      libmesh_assert(this->_size == data.size());
      //
      std::vector<Real> phi;
      this->init(xyz, phi);
      //
      Real out = 0.0;
      for (int i=0; i<this->_size; i++)
        out += phi[i] * data[i];
      return out;
    }

private:
  // private member function that computes the weights of the IDA algorithm
  inline
  void init (const Point& xyz, std::vector<Real>& phi) const
    {
      // evaluate the weights
      std::vector<Real> weights(this->_size);
      for (int i=0; i<this->_size; i++)
        {
          const Real ds = Point(*(this->_nodes[i])-xyz).norm();
          weights[i] = pow(ds, -this->_coeff);
        }
      const Real weights_sum = accumulate(weights.begin(), weights.end(), 0.0);
      //
      phi.clear();
      for (int i=0; i<this->_size; i++)
        phi.push_back(weights[i] / weights_sum);
    }

private:
  // exponent (positive-valued) of the inverse algorithm kernel
  Real _coeff;
  // pointer to an array of local node pointers
  const Node* const * _nodes;
  int _size;
};
//-------------------------------------------------------------------------------------------------

#endif // __IDA_H__

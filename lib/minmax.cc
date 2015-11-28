#include "minmax.hh"

using namespace sbb;

/* ********************************************************************************************* *
 * Implementation of MaximumDensityObj
 * ********************************************************************************************* */
MaximumDensityObj::MaximumDensityObj(const std::vector<RandomVariableObj *> &variables)
  : DensityObj(), _densities()
{
  _densities.reserve(variables.size());
  for (size_t i=0; i<variables.size(); i++) {
    _densities.push_back(variables[i]->density());
  }
}

MaximumDensityObj::~MaximumDensityObj() {
  // pass...
}

void
MaximumDensityObj::mark() {
  if (isMarked()) { return; }
  DensityObj::mark();
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->mark();
  }
}

void
MaximumDensityObj::eval(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  Eigen::VectorXd tmp(out.size());
  Eigen::MatrixXd pdfs(out.size(), _densities.size());
  Eigen::MatrixXd cdfs(out.size(), _densities.size());
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->eval(Tmin, Tmax, tmp);    pdfs.col(i) = tmp;
    _densities[i]->evalCDF(Tmin, Tmax, tmp); cdfs.col(i) = tmp;
  }

  out.setZero();
  for (size_t i=0; i<_densities.size(); i++) {
    tmp.setOnes();
    for (size_t j=0; j<_densities.size(); j++) {
      if (i == j) {
        tmp.array() *= pdfs.col(j).array();
      } else {
        tmp.array() *= cdfs.col(j).array();
      }
    }
    out += tmp;
  }
}

void
MaximumDensityObj::evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const {
  Eigen::VectorXd tmp(out.size());

  out.setOnes();
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->evalCDF(Tmin, Tmax, tmp);
    out.array() *= tmp.array();
  }
}

void
MaximumDensityObj::sample(Eigen::VectorXd &out) const {
  Eigen::VectorXd tmp;
  Eigen::MatrixXd samples(out.size(), _densities.size());
  for (size_t i=0; i<_densities.size(); i++) {
    _densities[i]->sample(tmp); samples.col(i) = tmp;
  }
  for (size_t i=0; i<out.size(); i++) {
    out[i] = samples.row(i).maxCoeff();
  }
}


/* ********************************************************************************************* *
 * Implementation of MaximumObj
 * ********************************************************************************************* */
MaximumObj::MaximumObj(RandomVariableObj *a, RandomVariableObj *b)
  : RandomVariableObj(), _variables(), _density(0)
{
  if (MaximumObj *max_a = dynamic_cast<MaximumObj *>(a)) {
    _variables = max_a->variables();
  } else {
    _variables.push_back(a);
  }
  if (MaximumObj *max_b = dynamic_cast<MaximumObj *>(b)) {
    for (size_t i=0; i<max_b->variables().size(); i++) {
      _variables.push_back(max_b->variables()[i]);
    }
  } else {
    _variables.push_back(b);
  }

  // Collect dependencies
  for (size_t i=0; i<_variables.size(); i++) {
    // Add implicit dependencies
    _dependencies.insert(_variables[i]->dependencies().begin(),
                         _variables[i]->dependencies().end());
    // add variable itself
    _dependencies.insert(_variables[i]);
  }

  _density = new MaximumDensityObj(_variables);
}

MaximumObj::MaximumObj(const std::vector<RandomVariableObj *> &variables)
  : RandomVariableObj(), _variables(variables), _density(0)
{
  // Collect dependencies
  for (size_t i=0; i<_variables.size(); i++) {
    // Add implicit dependencies
    _dependencies.insert(_variables[i]->dependencies().begin(),
                         _variables[i]->dependencies().end());
    // add variable itself
    _dependencies.insert(_variables[i]);
  }

  _density = new MaximumDensityObj(_variables);
}

MaximumObj::~MaximumObj() {
  // pass...
}

void
MaximumObj::mark() {
  if (isMarked()) { return; }
  RandomVariableObj::mark();
  for (size_t i=0; i<_variables.size(); i++) {
    _variables[i]->mark();
  }
  _density->mark();
}

DensityObj *
MaximumObj::density() {
  return _density;
}


/* ********************************************************************************************* *
 * Implementation of Maximum container
 * ********************************************************************************************* */
Maximum::Maximum(MaximumObj *obj)
  : RandomVariable(obj), _maximum(obj)
{
  // pass...
}

Maximum::Maximum(const RandomVariable &a, const RandomVariable &b)
  : RandomVariable(new MaximumObj(*a, *b)), _maximum(static_cast<MaximumObj *>(_randomVariable))
{
  // pass...
}

Maximum::Maximum(const Maximum &other)
  : RandomVariable(other), _maximum(other._maximum)
{
  // pass...
}

Maximum &
Maximum::operator =(const Maximum &other) {
  RandomVariable::operator =(other);
  _maximum = other._maximum;
  return *this;
}

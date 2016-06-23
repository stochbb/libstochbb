#include "derivedensity.hh"
#include "operators.hh"
#include "randomvariable.hh"
#include "compound.hh"
#include "reduction.hh"
#include "affinetrafo.hh"
#include "chain.hh"
#include "conditional.hh"
#include "minmax.hh"
#include "mixture.hh"

#include <list>
#include <map>

using namespace stochbb;

// Declaration of helper functions to derived the density of "derived" random variables.
void processCompoundVariable(const Compound &var, std::map<Var, Density> &vartable) throw (Error);
void processAffineTrafo(const AffineTrafo &var, std::map<Var, Density> &vartable) throw (Error);
void processSum(const Chain &var, std::map<Var, Density> &vartable) throw (Error);
void processConditional(const Conditional &var, std::map<Var, Density> &vartable) throw (Error);
void processCondSum(const CondSum &var, std::map<Var, Density> &vartable) throw (Error);
void processMaximum(const Maximum &var, std::map<Var, Density> &vartable) throw (Error);
void processMinimum(const Minimum &var, std::map<Var, Density> &vartable) throw (Error);
void processMixture(const Mixture &var, std::map<Var, Density> &vartable) throw (Error);

// tiny helper function to sort vectors of densities
inline int density_compare(const Density &a, const Density &b) {
  return a.compare(b);
}


/* ********************************************************************************************* *
 * convolve densities
 * ********************************************************************************************* */
Density
stochbb::convolve(const std::vector<Density> &densities, double scale, double shift) {
  // copy vector
  std::vector<Density> dens(densities);
  // Sort densities w.r.t type and parameters
  std::sort(dens.begin(), dens.end(), density_compare);

  // Get reduction rules
  ConvolutionReductions &rules = ConvolutionReductions::get();

  // Try to combine some of the densities
  std::vector<Density>::iterator last = dens.begin();
  std::vector<Density>::iterator current = dens.begin(); current++;
  while (current != dens.end()) {
    if (ConvolutionReductionRule *rule = rules.find(*last, *current)) {
      // If densities can be combined -> combine & replace last density
      *last = rule->apply(*last, *current);
      // erase combined density
      current = dens.erase(current);
    } else {
      // If densities cannot be combined -> advance iterators
      last++; current++;
    }
  }

  if (1 == dens.size()) {
    // If only one density is left -> unpack
    return dens.back();
  }

  // Otherwise construct convolution density from list of densities
  return new ConvolutionDensityObj(dens);
}


/* ********************************************************************************************* *
 * Implementation of deriveDensity
 * ********************************************************************************************* */
Density
stochbb::deriveDensity(const Var &var) throw (Error) {
  // Check if var is not null
  if (var.isNull()) {
    TypeError err;
    err << "Cannot derive PDF from an empty variable.";
    throw err;
  }

  std::list<Var> queue;
  std::map<Var, Density> vartable;
  std::set<Var> queued;
  queue.push_front(var);
  queued.insert(var);

  // Process variables
  while (queue.size()) {
    // Get the first variable from queue
    Var var = queue.front();

    // If the variable has been processed
    if (vartable.count(var)) {
      // -> remove from queue
      queue.pop_front();
      continue;
    } else if (var.is<AtomicVar>()) {
      // Atomic random variables have a density assigned
      vartable[var] = var.density();
    } else if (var.is<DerivedVar>()) {
      // cast to derived var
      DerivedVar derived(var.as<DerivedVar>());
      // a flag if all variables has been processed
      bool ready = true;

      // First process all RVs, a derived RV depends on
      for (size_t i=0; i<derived.numVariables(); i++) {
        // If the variable has not been processed yet
        if (! vartable.count(derived.variable(i))) {
          // -> add var to queue if not already there
          if (! queued.count(derived.variable(i)))
            queue.push_front(derived.variable(i));
          // Mark var as not ready yet
          ready = false;
        }
      }

      // If all variable have been processed, dispatch by type
      if (ready) {
        // Dispatch by type
        if (var.is<Compound>()) {
          processCompoundVariable(var.as<Compound>(), vartable);
        } else if (var.is<AffineTrafo>()) {
          processAffineTrafo(var.as<AffineTrafo>(), vartable);
        } else if (var.is<Chain>()) {
          processSum(var.as<Chain>(), vartable);
        } else if (var.is<Conditional>()) {
          processConditional(var.as<Conditional>(), vartable);
        } else if (var.is<CondSum>()) {
          processCondSum(var.as<CondSum>(), vartable);
        } else if (var.is<Maximum>()) {
          processMaximum(var.as<Maximum>(), vartable);
        } else if (var.is<Minimum>()) {
          processMinimum(var.as<Minimum>(), vartable);
        } else if (var.is<Mixture>()) {
          processMixture(var.as<Mixture>(), vartable);
        } else {
          TypeError err;
          err << "Internal error: Unknown derived variable type " << var
              << ": Density derivation is not implemented yet for this type.";
        }
      }
    } else {
      TypeError err;
      err << "Cannot derive density from variable " << var
          << ": It is neither an atomic nor derived random variable.";
      throw err;
    }
  }

  return vartable[var];
}


void
processCompoundVariable(const Compound &var, std::map<Var, Density> &vartable) throw (Error) {
  // First assemble vector of parameter variables
  std::vector<Var> paramvars; paramvars.reserve(var.numVariables());
  for (size_t i=0; i<var.numVariables(); i++) {
    paramvars.push_back(var.variable(i));
  }
  // check if all parameter variables are mutually independent
  if (! independent(paramvars)) {
    AssumptionError err;
    err << "Cannot derive density for compound random variable " << var
        << ": Parameter variables are not mutually independent!";
    throw err;
  }

  // Assemble parameter density vector
  std::vector<Density> params; params.reserve(var.numVariables());
  for (size_t i=0; i<var.numVariables(); i++) {
    params.push_back(vartable[var.variable(i)]);
  }

  // Construct density
  Density dens(new CompoundDensityObj(var.distribution(), params));

  // Apply reduction rules
  while (CompoundReductionRule *rule = CompoundReductions::get().find(dens)) {
    dens = rule->apply(dens);
  }

  // store result in vartable
  vartable[var] = dens;
}

void
processAffineTrafo(const AffineTrafo &var, std::map<Var, Density> &vartable) throw (Error) {
  // Pretty staight forward...
  // ... get the dens of the var ...
  Density dens = vartable[var.variable(0)];
  // ... and transform it
  vartable[var] = dens.affine(var.scale(), var.shift());
}

void
processSum(const Chain &var, std::map<Var, Density> &vartable) throw (Error) {
  // Get variables to be summed & densities
  std::vector<Var> variables; variables.reserve(var.numVariables());
  std::vector<Density> densities; densities.reserve(var.numVariables());
  for (size_t i=0; i<var.numVariables(); i++) {
    variables.push_back(var.variable(i));
    densities.push_back(vartable[variables.back()]);
  }

  // Check for independence
  if (! independent(variables)) {
    AssumptionError err;
    err << "Cannot create sum: Variables not independent.";
    throw err;
  }

  // Assemble convolution density
  vartable[var] = convolve(densities);
}

void
processConditional(const Conditional &var, std::map<Var, Density> &vartable) throw (Error) {
  // Get variables
  Var X1 = var.variable(0), X2 = var.variable(1), Y1 = var.variable(2), Y2 = var.variable(3);
  // Check for independence (Y1 and Y2 are allowed to be dependent RVs).
  if (! independent(std::vector<Var> {X1, X2, Y1})) {
    AssumptionError err;
    err << "Cannot instantiate conditional (X1<X2) ? Y1 : Y2."
           " Variables X1, X2, Y1 are not mutually independent.";
    throw err;
  }
  if (! independent(std::vector<Var> {X1, X2, Y2})) {
    AssumptionError err;
    err << "Cannot instantiate conditional (X1<X2) ? Y1 : Y2."
           " Variables X1, X2, Y2 are not mutually independent.";
    throw err;
  }

  // Store density
  vartable[var] = new ConditionalDensityObj(X1, X2, Y1, Y2);
}

void
processCondSum(const CondSum &var, std::map<Var, Density> &vartable) throw (Error) {
  // Get variables
  Var X1 = var.variable(0), X2 = var.variable(1), Y1 = var.variable(2), Y2 = var.variable(3);
  // Check for independence (Y1 and Y2 are allowed to be dependent RVs).
  if (! independent(std::vector<Var> {X1, X2, Y1})) {
    AssumptionError err;
    err << "Cannot instantiate conditional sum (X1<X2) ? Y1+X1 : Y2+X2."
           " Variables X1, X2, Y1 are not mutually independent.";
    throw err;
  }
  if (! independent(std::vector<Var> {X1, X2, Y2})) {
    AssumptionError err;
    err << "Cannot instantiate conditional sum (X1<X2) ? Y1+X1 : Y2+X2."
           " Variables X1, X2, Y2 are not mutually independent.";
    throw err;
  }

  // Store density
  vartable[var] = new CondSumDensityObj(X1, X2, Y1, Y2);
}

void
processMaximum(const Maximum &var, std::map<Var, Density> &vartable) throw (Error) {
  // Get variables
  std::vector<Var> variables; variables.reserve(var.numVariables());
  for (size_t i=0; i<var.numVariables(); i++) {
    variables.push_back(var.variable(i));
  }
  // check for independence
  if (! independent(variables)) {
    AssumptionError err;
    err << "Cannot derive maximum density: Variables are not mutually independent.";
    throw err;
  }
  vartable[var] = new MaximumDensityObj(variables);
}

void
processMinimum(const Minimum &var, std::map<Var, Density> &vartable) throw(Error) {
  // Get variables
  std::vector<Var> variables; variables.reserve(var.numVariables());
  for (size_t i=0; i<var.numVariables(); i++) {
    variables.push_back(var.variable(i));
  }
  // check for independence
  if (! independent(variables)) {
    AssumptionError err;
    err << "Cannot derive minimum density: Variables are not mutually independent.";
    throw err;
  }
  vartable[var] = new MinimumDensityObj(variables);
}

void
processMixture(const Mixture &var, std::map<Var, Density> &vartable) throw(Error) {
  // Get variables
  std::vector<Var> variables; variables.reserve(var.numVariables());
  std::vector<double> weights; weights.reserve(var.numVariables());
  for (size_t i=0; i<var.numVariables(); i++) {
    variables.push_back(var.variable(i));
    weights.push_back(var.weight(i));
  }
  // check for independence
  if (! independent(variables)) {
    AssumptionError err;
    err << "Cannot derive mixture density: Variables are not mutually independent.";
    throw err;
  }
  vartable[var] = new MixtureDensityObj(weights, variables);
}

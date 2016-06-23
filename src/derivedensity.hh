#ifndef PDFASSEMBLER_HH
#define PDFASSEMBLER_HH

#include "api.hh"

namespace stochbb {

/** Returns the convoultion of the given densities.
 * This function first tries to perform the convolutions analytically and
 * resorts to the numerical convolution if necessary.
 * @ingroup density */
Density convolve(const std::vector<Density> &densities, double scale=1, double shift=0);

/** Derives the density for the given variable.
 * The derivation may fail with an exception if some assumptions (usually independence assuptions)
 * are not met, that are required for the numeric or analytic derivation of the density. */
Density deriveDensity(const Var &var) throw (Error);

}

#endif // PDFASSEMBLER_HH

#ifndef KINEMATICS_H
#define KINEMATICS_H

/**
 * Expected 1/beta for a particle of mass [GeV/c^2] and momentum p [GeV/c].
 * 1/beta = E/p = sqrt(m^2 + p^2) / p.
 * Returns -1 if p <= 0.
 */
double OneOverBetaExpected(double mass, double p);

/**
 * Delta(1/beta) = 1/beta_measured - 1/beta_expected(mass, p).
 * Returns the signed difference; use fabs() for a cut on |Delta(1/beta)|.
 * If betaMeasured <= 0 or p <= 0, returns 999.0 so caller can skip.
 */
double DeltaOneOverBeta(double betaMeasured, double mass, double p);

#endif

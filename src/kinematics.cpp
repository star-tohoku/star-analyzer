#include "kinematics.h"
#include "TMath.h"

double OneOverBetaExpected(double mass, double p) {
  if (p <= 0.) return -1.;
  double E = TMath::Sqrt(mass * mass + p * p);
  return E / p;
}

double DeltaOneOverBeta(double betaMeasured, double mass, double p) {
  const double kInvalid = 999.0;
  if (betaMeasured <= 1e-6 || p <= 1e-6) return kInvalid;
  double oneOverBetaMeas = 1.0 / betaMeasured;
  double oneOverBetaExp = OneOverBetaExpected(mass, p);
  if (oneOverBetaExp < 0.) return kInvalid;
  return oneOverBetaMeas - oneOverBetaExp;
}

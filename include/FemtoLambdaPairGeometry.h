#ifndef FEMTO_LAMBDA_PAIR_GEOMETRY_H
#define FEMTO_LAMBDA_PAIR_GEOMETRY_H
#include "TVector2.h"
#include <cmath>
#include <map>
#include <stdexcept>
namespace femto_lambda_pair_geometry {
struct TrackPoint {
  double phi,eta;
  bool valid;
  TrackPoint():phi(0.),eta(0.),valid(false){}
};
// Track-array indices are local to an event. Each point uses its own event field.
struct EventGeometry {
  int eventIndex;
  std::map<int,TrackPoint> points;
  EventGeometry():eventIndex(-1){}
  void Clear(){eventIndex=-1;points.clear();}
};
inline const TrackPoint& Find(const EventGeometry* event,int expectedEvent,int index) {
  if(!event||expectedEvent<1||index<0||event->eventIndex!=expectedEvent)
    throw std::runtime_error("Lambda pair geometry has a missing/wrong event identity");
  const std::map<int,TrackPoint>::const_iterator point=event->points.find(index);
  if(point==event->points.end())throw std::runtime_error("Lambda pair geometry is missing a required saved track");
  return point->second;
}
inline bool Difference(const TrackPoint& a,const TrackPoint& b,double& phi,double& eta) {
  if(!a.valid||!b.valid||!std::isfinite(a.phi)||!std::isfinite(a.eta)||
     !std::isfinite(b.phi)||!std::isfinite(b.eta))return false;
  const double rawPhi=a.phi-b.phi;
  if(!std::isfinite(rawPhi))return false;
  phi=TVector2::Phi_mpi_pi(rawPhi);eta=a.eta-b.eta;
  return std::isfinite(phi)&&std::isfinite(eta);
}
inline bool InsideRectangle(double phi,double eta,double phiMax,double etaMax) {
  return std::isfinite(phi)&&std::isfinite(eta)&&std::fabs(phi)<phiMax&&std::fabs(eta)<etaMax;
}
}
#endif

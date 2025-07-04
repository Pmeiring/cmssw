#ifndef DataFormats_TrackingRecHitSoA_test_alpaka_Hits_test_h
#define DataFormats_TrackingRecHitSoA_test_alpaka_Hits_test_h

#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::testTrackingRecHitSoA {

  void runKernels(::reco::TrackingRecHitView& hits, ::reco::HitModuleSoAView& mods, Queue& queue);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::testTrackingRecHitSoA

#endif  // DataFormats_TrackingRecHitSoA_test_alpaka_Hits_test_h

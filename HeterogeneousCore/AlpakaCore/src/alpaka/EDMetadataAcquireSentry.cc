#include "HeterogeneousCore/AlpakaCore/interface/alpaka/chooseDevice.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDMetadataAcquireSentry.h"
#include "HeterogeneousCore/AlpakaCore/interface/EventCache.h"
#include "HeterogeneousCore/AlpakaCore/interface/QueueCache.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace detail {
    EDMetadataAcquireSentry::EDMetadataAcquireSentry(edm::StreamID streamID,
                                                     edm::WaitingTaskWithArenaHolder holder,
                                                     bool synchronize)
        : EDMetadataAcquireSentry(detail::chooseDevice(streamID), std::move(holder), synchronize) {}

    EDMetadataAcquireSentry::EDMetadataAcquireSentry(Device const& device,
                                                     edm::WaitingTaskWithArenaHolder holder,
                                                     bool synchronize)
        : waitingTaskHolder_(std::move(holder)), synchronize_(synchronize) {
#ifdef ALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED
      // all synchronous backends
      metadata_ = std::make_shared<EDMetadata>(cms::alpakatools::getQueueCache<Queue>().get(device));
#else
      // all asynchronous backends
      metadata_ = std::make_shared<EDMetadata>(cms::alpakatools::getQueueCache<Queue>().get(device),
                                               cms::alpakatools::getEventCache<Event>().get(device));
#endif
    }

#ifndef ALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED
    // all asynchronous backends
    std::shared_ptr<EDMetadata> EDMetadataAcquireSentry::finish() {
      if (synchronize_) {
        alpaka::wait(metadata_->queue());
      } else {
        metadata_->enqueueCallback(std::move(waitingTaskHolder_));
      }
      return std::move(metadata_);
    }
#endif
  }  // namespace detail
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

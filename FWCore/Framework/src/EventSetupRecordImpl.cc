// -*- C++ -*-
//
// Package:     Framework
// Class  :     EventSetupRecordImpl
//
// Implementation:
//     <Notes on implementation>
//
// Author:      Chris Jones
// Created:     Sat Mar 26 18:06:32 EST 2005
//

// system include files
#include <algorithm>
#include <cassert>
#include <iterator>
#include <sstream>
#include <utility>

// user include files
#include "FWCore/Framework/interface/EventSetupRecordImpl.h"
#include "FWCore/Framework/interface/ESProductResolver.h"
#include "FWCore/Framework/interface/ComponentDescription.h"
#include "FWCore/ServiceRegistry/interface/ESParentContext.h"

#include "FWCore/Utilities/interface/ConvertException.h"
#include "FWCore/Utilities/interface/Exception.h"

namespace edm {
  namespace eventsetup {

    EventSetupRecordImpl::EventSetupRecordImpl(EventSetupRecordKey const& iKey,
                                               ActivityRegistry const* activityRegistry,
                                               unsigned int iovIndex)
        : validity_(),
          key_(iKey),
          activityRegistry_(activityRegistry),
          cacheIdentifier_(1),  // initial value meaningless, gets overwritten before use
          iovIndex_(iovIndex),
          isAvailable_(true),
          validityModificationUnderway_(false) {}

    EventSetupRecordImpl::EventSetupRecordImpl(EventSetupRecordImpl&& source)
        : validity_{source.validity_},
          key_{source.key_},
          keysForResolvers_{std::move(source.keysForResolvers_)},
          resolvers_(std::move(source.resolvers_)),
          activityRegistry_{source.activityRegistry_},
          cacheIdentifier_{source.cacheIdentifier_},
          iovIndex_{source.iovIndex_},
          isAvailable_{source.isAvailable_.load()},
          validityModificationUnderway_{source.validityModificationUnderway_.load()} {}

    EventSetupRecordImpl& EventSetupRecordImpl::operator=(EventSetupRecordImpl&& rhs) {
      validity_ = rhs.validity_;
      key_ = rhs.key_;
      keysForResolvers_ = std::move(rhs.keysForResolvers_);
      resolvers_ = std::move(rhs.resolvers_);
      activityRegistry_ = rhs.activityRegistry_;
      cacheIdentifier_ = rhs.cacheIdentifier_;
      iovIndex_ = rhs.iovIndex_;
      isAvailable_.store(rhs.isAvailable_.load());
      validityModificationUnderway_.store(validityModificationUnderway_.load());
      return *this;
    }

    ValidityInterval EventSetupRecordImpl::validityInterval() const {
      bool expected = false;
      while (not validityModificationUnderway_.compare_exchange_strong(expected, true)) {
        expected = false;
      }
      ValidityInterval temp = validity_;
      validityModificationUnderway_ = false;
      return temp;
    }

    void EventSetupRecordImpl::initializeForNewIOV(unsigned long long iCacheIdentifier,
                                                   ValidityInterval const& iValidityInterval,
                                                   bool hasFinder) {
      cacheIdentifier_ = iCacheIdentifier;
      validity_ = iValidityInterval;
      if (hasFinder) {
        for (auto& productResolver : resolvers_) {
          productResolver->initializeForNewIOV();
        }
      }
    }

    void EventSetupRecordImpl::setSafely(const ValidityInterval& iInterval) const {
      bool expected = false;
      while (not validityModificationUnderway_.compare_exchange_strong(expected, true)) {
        expected = false;
      }
      validity_ = iInterval;
      validityModificationUnderway_ = false;
    }

    void EventSetupRecordImpl::getESProducers(std::vector<ComponentDescription const*>& esproducers) const {
      esproducers.clear();
      esproducers.reserve(resolvers_.size());
      for (auto const& iData : resolvers_) {
        ComponentDescription const* componentDescription = iData->providerDescription();
        if (!componentDescription->isLooper_ && !componentDescription->isSource_) {
          esproducers.push_back(componentDescription);
        }
      }
    }

    std::vector<ComponentDescription const*> EventSetupRecordImpl::componentsForRegisteredDataKeys() const {
      std::vector<ComponentDescription const*> ret;
      ret.reserve(resolvers_.size());
      for (auto const& resolver : resolvers_) {
        ret.push_back(resolver->providerDescription());
      }
      return ret;
    }

    std::vector<unsigned int> EventSetupRecordImpl::produceMethodIDsForRegisteredDataKeys() const {
      std::vector<unsigned int> ret;
      ret.reserve(resolvers_.size());
      for (auto const& resolver : resolvers_) {
        ret.push_back(resolver->produceMethodID());
      }
      return ret;
    }

    bool EventSetupRecordImpl::add(const DataKey& iKey, ESProductResolver* iResolver) {
      const ESProductResolver* resolver = find(iKey);
      if (nullptr != resolver) {
        //
        // we already know the field exist, so do not need to check against end()
        //

        // POLICY: If a Producer and a Source both claim to deliver the same data, the
        //  Producer 'trumps' the Source. If two modules of the same type claim to deliver the
        //  same data, this is an error unless the configuration specifically states which one
        //  is to be chosen.  A Looper trumps both a Producer and a Source.

        assert(resolver->providerDescription());
        assert(iResolver->providerDescription());
        if (iResolver->providerDescription()->isLooper_) {
          resolvers_[std::distance(keysForResolvers_.begin(),
                                   std::lower_bound(keysForResolvers_.begin(), keysForResolvers_.end(), iKey))] =
              iResolver;
          return true;
        }

        if (resolver->providerDescription()->isSource_ == iResolver->providerDescription()->isSource_) {
          //should lookup to see if there is a specified 'chosen' one and only if not, throw the exception
          throw cms::Exception("EventSetupConflict")
              << "two EventSetup " << (resolver->providerDescription()->isSource_ ? "Sources" : "Producers")
              << " want to deliver type=\"" << iKey.type().name() << "\" label=\"" << iKey.name().value() << "\"\n"
              << " from record " << key().type().name() << ". The two providers are \n"
              << "1) type=\"" << resolver->providerDescription()->type_ << "\" label=\""
              << resolver->providerDescription()->label_ << "\"\n"
              << "2) type=\"" << iResolver->providerDescription()->type_ << "\" label=\""
              << iResolver->providerDescription()->label_ << "\"\n"
              << "Please either\n   remove one of these "
              << (resolver->providerDescription()->isSource_ ? "Sources" : "Producers")
              << "\n   or find a way of configuring one of them so it does not deliver this data"
              << "\n   or use an es_prefer statement in the configuration to choose one.";
        } else if (resolver->providerDescription()->isSource_) {
          resolvers_[std::distance(keysForResolvers_.begin(),
                                   std::lower_bound(keysForResolvers_.begin(), keysForResolvers_.end(), iKey))] =
              iResolver;
        } else {
          return false;
        }
      } else {
        auto lb = std::lower_bound(keysForResolvers_.begin(), keysForResolvers_.end(), iKey);
        auto index = std::distance(keysForResolvers_.begin(), lb);
        keysForResolvers_.insert(lb, iKey);
        resolvers_.insert(resolvers_.begin() + index, iResolver);
      }
      return true;
    }

    void EventSetupRecordImpl::clearResolvers() {
      keysForResolvers_.clear();
      resolvers_.clear();
    }

    void EventSetupRecordImpl::invalidateResolvers() {
      for (auto& productResolver : resolvers_) {
        productResolver->invalidate();
      }
    }

    void EventSetupRecordImpl::resetIfTransientInResolvers() {
      for (auto& productResolver : resolvers_) {
        productResolver->resetIfTransient();
      }
    }

    void const* EventSetupRecordImpl::getFromResolverAfterPrefetch(ESResolverIndex iResolverIndex,
                                                                   bool iTransientAccessOnly,
                                                                   ComponentDescription const*& iDesc,
                                                                   DataKey const*& oGottenKey) const {
      const ESProductResolver* resolver = resolvers_[iResolverIndex.value()];
      assert(nullptr != resolver);
      iDesc = resolver->providerDescription();

      auto const& key = keysForResolvers_[iResolverIndex.value()];
      oGottenKey = &key;

      void const* hold = nullptr;
      try {
        convertException::wrap([&]() { hold = resolver->getAfterPrefetch(*this, key, iTransientAccessOnly); });
      } catch (cms::Exception& e) {
        addTraceInfoToCmsException(e, key.name().value(), resolver->providerDescription(), key);
        throw;
      }
      return hold;
    }

    const ESProductResolver* EventSetupRecordImpl::find(const DataKey& iKey) const {
      auto lb = std::lower_bound(keysForResolvers_.begin(), keysForResolvers_.end(), iKey);
      if ((lb == keysForResolvers_.end()) or (*lb != iKey)) {
        return nullptr;
      }
      return resolvers_[std::distance(keysForResolvers_.begin(), lb)].get();
    }

    void EventSetupRecordImpl::prefetchAsync(WaitingTaskHolder iTask,
                                             ESResolverIndex iResolverIndex,
                                             EventSetupImpl const* iEventSetupImpl,
                                             ServiceToken const& iToken,
                                             ESParentContext iParent) const noexcept {
      if UNLIKELY (iResolverIndex == ESResolverIndex::moduleLabelDoesNotMatch() ||
                   iResolverIndex == ESResolverIndex::noResolverConfigured()) {
        return;
      }

      const ESProductResolver* resolver = resolvers_[iResolverIndex.value()];
      if (nullptr != resolver) {
        auto const& key = keysForResolvers_[iResolverIndex.value()];
        resolver->prefetchAsync(iTask, *this, key, iEventSetupImpl, iToken, iParent);
      }
    }

    bool EventSetupRecordImpl::wasGotten(const DataKey& aKey) const {
      const ESProductResolver* resolver = find(aKey);
      if (nullptr != resolver) {
        return resolver->cacheIsValid();
      }
      return false;
    }

    edm::eventsetup::ComponentDescription const* EventSetupRecordImpl::providerDescription(const DataKey& aKey) const {
      const ESProductResolver* resolver = find(aKey);
      if (nullptr != resolver) {
        return resolver->providerDescription();
      }
      return nullptr;
    }

    void EventSetupRecordImpl::fillRegisteredDataKeys(std::vector<DataKey>& oToFill) const {
      oToFill = keysForResolvers_;
    }

    void EventSetupRecordImpl::addTraceInfoToCmsException(cms::Exception& iException,
                                                          const char* iName,
                                                          const ComponentDescription* iDescription,
                                                          const DataKey& iKey) const {
      std::ostringstream ost;
      ost << "Using EventSetup component " << iDescription->type_ << "/'" << iDescription->label_ << "' to make data "
          << iKey.type().name() << "/'" << iName << "' in record " << this->key().type().name();
      iException.addContext(ost.str());
    }

  }  // namespace eventsetup
}  // namespace edm

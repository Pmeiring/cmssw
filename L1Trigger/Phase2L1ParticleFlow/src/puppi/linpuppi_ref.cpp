#include "L1Trigger/Phase2L1ParticleFlow/interface/puppi/linpuppi_ref.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/puppi/linpuppi_bits.h"
#include <cmath>
#include <algorithm>

#include "L1Trigger/Phase2L1ParticleFlow/interface/common/bitonic_hybrid_sort_ref.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/common/bitonic_sort_ref.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/dbgPrintf.h"

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/allowedValues.h"
#include "FWCore/Utilities/interface/transform.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#endif

#include "L1Trigger/Phase2L1ParticleFlow/interface/NNVtxAssoc.h"

using namespace l1ct;
using namespace linpuppi;

l1ct::LinPuppiEmulator::LinPuppiEmulator(unsigned int nTrack,
                                         unsigned int nIn,
                                         unsigned int nOut,
                                         unsigned int nVtx,
                                         unsigned int dR2Min,
                                         unsigned int dR2Max,
                                         unsigned int iptMax,
                                         unsigned int dzCut,
                                         glbeta_t etaCut,
                                         double ptSlopeNe_0,
                                         double ptSlopeNe_1,
                                         double ptSlopePh_0,
                                         double ptSlopePh_1,
                                         double ptZeroNe_0,
                                         double ptZeroNe_1,
                                         double ptZeroPh_0,
                                         double ptZeroPh_1,
                                         double alphaSlope_0,
                                         double alphaSlope_1,
                                         double alphaZero_0,
                                         double alphaZero_1,
                                         double alphaCrop_0,
                                         double alphaCrop_1,
                                         double priorNe_0,
                                         double priorNe_1,
                                         double priorPh_0,
                                         double priorPh_1,
                                         pt_t ptCut_0,
                                         pt_t ptCut_1,
                                         bool useMLAssociation,
                                         const double associationThreshold,
                                         std::string associationGraphPath,
                                         std::vector<double> associationNetworkZ0binning,
                                         std::vector<double> associationNetworkEtaBounds,
                                         std::vector<double> associationNetworkZ0ResBins,
                                         unsigned int nFinalSort,
                                         SortAlgo finalSortAlgo)
    : nTrack_(nTrack),
      nIn_(nIn),
      nOut_(nOut),
      nVtx_(nVtx),
      dR2Min_(dR2Min),
      dR2Max_(dR2Max),
      iptMax_(iptMax),
      dzCut_(dzCut),
      absEtaBins_(1, etaCut),
      ptSlopeNe_(2),
      ptSlopePh_(2),
      ptZeroNe_(2),
      ptZeroPh_(2),
      alphaSlope_(2),
      alphaZero_(2),
      alphaCrop_(2),
      priorNe_(2),
      priorPh_(2),
      ptCut_(2),
      useMLAssociation_(useMLAssociation),
      nFinalSort_(nFinalSort ? nFinalSort : nOut),
      finalSortAlgo_(finalSortAlgo),
      debug_(false),
      fakePuppi_(false) {
  ptSlopeNe_[0] = ptSlopeNe_0;
  ptSlopeNe_[1] = ptSlopeNe_1;
  ptSlopePh_[0] = ptSlopePh_0;
  ptSlopePh_[1] = ptSlopePh_1;
  ptZeroNe_[0] = ptZeroNe_0;
  ptZeroNe_[1] = ptZeroNe_1;
  ptZeroPh_[0] = ptZeroPh_0;
  ptZeroPh_[1] = ptZeroPh_1;
  alphaSlope_[0] = alphaSlope_0;
  alphaSlope_[1] = alphaSlope_1;
  alphaZero_[0] = alphaZero_0;
  alphaZero_[1] = alphaZero_1;
  alphaCrop_[0] = alphaCrop_0;
  alphaCrop_[1] = alphaCrop_1;
  priorNe_[0] = priorNe_0;
  priorNe_[1] = priorNe_1;
  priorPh_[0] = priorPh_0;
  priorPh_[1] = priorPh_1;
  ptCut_[0] = ptCut_0;
  ptCut_[1] = ptCut_1;

  if (useMLAssociation_ and withinCMSSW_) {
    nnVtxAssoc_ = std::make_unique<NNVtxAssoc>(NNVtxAssoc(associationGraphPath,
                                                          associationThreshold,
                                                          associationNetworkZ0binning,
                                                          associationNetworkEtaBounds,
                                                          associationNetworkZ0ResBins));
  }
}

#ifdef CMSSW_GIT_HASH
l1ct::LinPuppiEmulator::LinPuppiEmulator(const edm::ParameterSet &iConfig)
    : nTrack_(iConfig.getParameter<uint32_t>("nTrack")),
      nIn_(iConfig.getParameter<uint32_t>("nIn")),
      nOut_(iConfig.getParameter<uint32_t>("nOut")),
      nVtx_(iConfig.getParameter<uint32_t>("nVtx")),
      dR2Min_(l1ct::Scales::makeDR2FromFloatDR(iConfig.getParameter<double>("drMin"))),
      dR2Max_(l1ct::Scales::makeDR2FromFloatDR(iConfig.getParameter<double>("dr"))),
      iptMax_(l1ct::Scales::intPt(l1ct::Scales::makePtFromFloat(iConfig.getParameter<double>("ptMax")))),
      dzCut_(l1ct::Scales::makeZ0(iConfig.getParameter<double>("dZ"))),
      absEtaBins_(
          edm::vector_transform(iConfig.getParameter<std::vector<double>>("absEtaCuts"), l1ct::Scales::makeGlbEta)),
      ptSlopeNe_(iConfig.getParameter<std::vector<double>>("ptSlopes")),
      ptSlopePh_(iConfig.getParameter<std::vector<double>>("ptSlopesPhoton")),
      ptZeroNe_(iConfig.getParameter<std::vector<double>>("ptZeros")),
      ptZeroPh_(iConfig.getParameter<std::vector<double>>("ptZerosPhoton")),
      alphaSlope_(iConfig.getParameter<std::vector<double>>("alphaSlopes")),
      alphaZero_(iConfig.getParameter<std::vector<double>>("alphaZeros")),
      alphaCrop_(iConfig.getParameter<std::vector<double>>("alphaCrop")),
      priorNe_(iConfig.getParameter<std::vector<double>>("priors")),
      priorPh_(iConfig.getParameter<std::vector<double>>("priorsPhoton")),
      ptCut_(edm::vector_transform(iConfig.getParameter<std::vector<double>>("ptCut"), l1ct::Scales::makePtFromFloat)),
      useMLAssociation_(iConfig.getParameter<bool>("useMLAssociation")),
      nFinalSort_(iConfig.getParameter<uint32_t>("nFinalSort")),
      debug_(iConfig.getUntrackedParameter<bool>("debug", false)),
      fakePuppi_(iConfig.getParameter<bool>("fakePuppi")) {
  if (absEtaBins_.size() + 1 != ptSlopeNe_.size())
    throw cms::Exception("Configuration", "size mismatch for ptSlopes parameter");
  if (absEtaBins_.size() + 1 != ptSlopePh_.size())
    throw cms::Exception("Configuration", "size mismatch for ptSlopesPhoton parameter");
  if (absEtaBins_.size() + 1 != ptZeroPh_.size())
    throw cms::Exception("Configuration", "size mismatch for ptZeros parameter");
  if (absEtaBins_.size() + 1 != ptZeroNe_.size())
    throw cms::Exception("Configuration", "size mismatch for ptZerosPhotons parameter");
  if (absEtaBins_.size() + 1 != priorPh_.size())
    throw cms::Exception("Configuration", "size mismatch for priors parameter");
  if (absEtaBins_.size() + 1 != priorNe_.size())
    throw cms::Exception("Configuration", "size mismatch for priorsPhotons parameter");
  if (absEtaBins_.size() + 1 != alphaSlope_.size())
    throw cms::Exception("Configuration", "size mismatch for alphaSlope parameter");
  if (absEtaBins_.size() + 1 != alphaZero_.size())
    throw cms::Exception("Configuration", "size mismatch for alphaZero parameter");
  if (absEtaBins_.size() + 1 != alphaCrop_.size())
    throw cms::Exception("Configuration", "size mismatch for alphaCrop parameter");
  if (absEtaBins_.size() + 1 != ptCut_.size())
    throw cms::Exception("Configuration", "size mismatch for ptCut parameter");
  if (useMLAssociation_) {
    edm::ParameterSet nnVtxAssocPSet_ = iConfig.getParameter<edm::ParameterSet>("NNVtxAssociation");
    edm::FileInPath associationGraphPathFIP =
        edm::FileInPath(nnVtxAssocPSet_.getParameter<std::string>("associationGraph"));

    nnVtxAssoc_ = std::make_unique<NNVtxAssoc>(
        NNVtxAssoc(associationGraphPathFIP.fullPath(),
                   nnVtxAssocPSet_.getParameter<double>("associationThreshold"),
                   nnVtxAssocPSet_.getParameter<std::vector<double>>("associationNetworkZ0binning"),
                   nnVtxAssocPSet_.getParameter<std::vector<double>>("associationNetworkEtaBounds"),
                   nnVtxAssocPSet_.getParameter<std::vector<double>>("associationNetworkZ0ResBins")));
  }
  const std::string &sortAlgo = iConfig.getParameter<std::string>("finalSortAlgo");
  if (sortAlgo == "Insertion")
    finalSortAlgo_ = SortAlgo::Insertion;
  else if (sortAlgo == "BitonicRUFL")
    finalSortAlgo_ = SortAlgo::BitonicRUFL;
  else if (sortAlgo == "BitonicHLS")
    finalSortAlgo_ = SortAlgo::BitonicHLS;
  else if (sortAlgo == "BitonicVHDL")
    finalSortAlgo_ = SortAlgo::BitonicVHDL;
  else if (sortAlgo == "Hybrid")
    finalSortAlgo_ = SortAlgo::Hybrid;
  else if (sortAlgo == "FoldedHybrid")
    finalSortAlgo_ = SortAlgo::FoldedHybrid;
  else
    throw cms::Exception("Configuration", "unsupported finalSortAlgo '" + sortAlgo + "'");
}

edm::ParameterSetDescription l1ct::LinPuppiEmulator::getParameterSetDescription() {
  edm::ParameterSetDescription description;
  description.add<uint32_t>("nTrack");
  description.add<uint32_t>("nIn");
  description.add<uint32_t>("nOut");
  description.add<uint32_t>("nVtx", 1);
  description.add<double>("dZ");
  description.add<double>("dr");
  description.add<double>("drMin");
  description.add<double>("ptMax");
  description.add<std::vector<double>>("absEtaCuts");
  description.add<std::vector<double>>("ptCut");
  description.add<bool>("useMLAssociation");
  description.add<edm::ParameterSetDescription>("NNVtxAssociation", NNVtxAssoc::getParameterSetDescription());
  description.add<std::vector<double>>("ptSlopes");
  description.add<std::vector<double>>("ptSlopesPhoton");
  description.add<std::vector<double>>("ptZeros");
  description.add<std::vector<double>>("ptZerosPhoton");
  description.add<std::vector<double>>("alphaSlopes");
  description.add<std::vector<double>>("alphaZeros");
  description.add<std::vector<double>>("alphaCrop");
  description.add<std::vector<double>>("priors");
  description.add<std::vector<double>>("priorsPhoton");
  description.add<uint32_t>("nFinalSort");
  description.ifValue(edm::ParameterDescription<std::string>("finalSortAlgo", "Insertion", true),
                      edm::allowedValues<std::string>(
                          "Insertion", "BitonicRUFL", "BitonicHLS", "Hybrid", "FoldedHybrid", "BitonicVHDL"));
  description.add<bool>("fakePuppi", false);
  description.addUntracked<bool>("debug", false);
  return description;
}
#endif

void l1ct::LinPuppiEmulator::puppisort_and_crop_ref(unsigned int nOutMax,
                                                    const std::vector<l1ct::PuppiObjEmu> &in,
                                                    std::vector<l1ct::PuppiObjEmu> &out,
                                                    SortAlgo sortAlgo) {
  const unsigned int nOut = std::min<unsigned int>(nOutMax, in.size());
  out.resize(nOut);
  for (unsigned int iout = 0; iout < nOut; ++iout) {
    out[iout].clear();
  }

  if (sortAlgo == SortAlgo::Insertion) {
    for (unsigned int it = 0, nIn = in.size(); it < nIn; ++it) {
      for (int iout = int(nOut) - 1; iout >= 0; --iout) {
        if (out[iout].hwPt <= in[it].hwPt) {
          if (iout == 0 || out[iout - 1].hwPt > in[it].hwPt) {
            out[iout] = in[it];
          } else {
            out[iout] = out[iout - 1];
          }
        }
      }
    }
  } else if (sortAlgo == SortAlgo::BitonicRUFL) {
    bitonic_sort_and_crop_ref(in.size(), nOut, &in[0], &out[0]);
  } else if (sortAlgo == SortAlgo::BitonicHLS || sortAlgo == SortAlgo::Hybrid) {
    hybrid_bitonic_sort_and_crop_ref(in.size(), nOut, &in[0], &out[0], sortAlgo == SortAlgo::Hybrid);
  } else if (sortAlgo == SortAlgo::FoldedHybrid) {
    folded_hybrid_bitonic_sort_and_crop_ref(in.size(), nOut, &in[0], &out[0], true);
  } else if (sortAlgo == SortAlgo::BitonicVHDL) {
    // The VHDL version always takes power-of-2 inputs
    // (Nominally it produces the same size output, though things may get optimized away in the implementation)

    // find the po2 that's bigger than the input size
    unsigned int nextpo2 = 1;
    while (nextpo2 < in.size()) {
      nextpo2 <<= 1;
    }
    std::vector<l1ct::PuppiObjEmu> inPadded(nextpo2);
    for (unsigned int i = 0; i < in.size(); i++) {
      inPadded[i] = in[i];
    }
    for (unsigned int i = in.size(); i < nextpo2; i++) {
      inPadded[i].clear();
    }
    std::vector<l1ct::PuppiObjEmu> outPadded(nextpo2);
    bitonic_sort_and_crop_ref(nextpo2, nextpo2, inPadded.data(), outPadded.data());
    for (unsigned int i = 0; i < nOut; i++) {
      out[i] = outPadded[i];
    }
  }
}

void l1ct::LinPuppiEmulator::linpuppi_chs_ref(const PFRegionEmu &region,
                                              const std::vector<PVObjEmu> &pv,
                                              const std::vector<PFChargedObjEmu> &pfch /*[nTrack]*/,
                                              std::vector<PuppiObjEmu> &outallch /*[nTrack]*/) const {
  const unsigned int nTrack = std::min<unsigned int>(nTrack_, pfch.size());
  const unsigned int nVtx = std::min<unsigned int>(nVtx_, pv.size());
  outallch.resize(nTrack);
  for (unsigned int i = 0; i < nTrack; ++i) {
    int pZ0 = pfch[i].hwZ0;
    int z0diff = -99999;
    bool pass_network = false;
    for (unsigned int j = 0; j < nVtx; ++j) {
      int pZ0Diff = pZ0 - pv[j].hwZ0;
      if (std::abs(z0diff) > std::abs(pZ0Diff))
        z0diff = pZ0Diff;
      if (useMLAssociation_ and withinCMSSW_ &&
          nnVtxAssoc_->TTTrackNetworkSelector<const l1ct::PFChargedObjEmu>(region, pfch[i], pv[j]) == 1)
        pass_network = true;
    }
    bool accept = pfch[i].hwPt != 0;
    if (!fakePuppi_ && !useMLAssociation_)
      accept = accept && region.isFiducial(pfch[i]) && (std::abs(z0diff) <= int(dzCut_) || pfch[i].hwId.isMuon());
    if (!fakePuppi_ && useMLAssociation_)
      accept = accept && region.isFiducial(pfch[i]) && (pass_network || pfch[i].hwId.isMuon());
    if (accept) {
      outallch[i].fill(region, pfch[i]);
      if (fakePuppi_) {                           // overwrite Dxy & TkQuality with debug information
        outallch[i].setHwDxy(dxy_t(pv[0].hwZ0));  ///hack to get this to work
        outallch[i].setHwTkQuality(region.isFiducial(pfch[i]) ? 1 : 0);
      }
      if (debug_ && pfch[i].hwPt > 0)
        dbgPrintf("ref candidate %02u pt %7.2f pid %1d   vz %+6d  dz %+6d (cut %5d), fid %1d -> pass, packed %s\n",
                  i,
                  pfch[i].floatPt(),
                  pfch[i].intId(),
                  int(pfch[i].hwZ0),
                  z0diff,
                  dzCut_,
                  region.isFiducial(pfch[i]),
                  outallch[i].pack().to_string(16).c_str());
    } else {
      outallch[i].clear();
      if (debug_ && pfch[i].hwPt > 0)
        dbgPrintf("ref candidate %02u pt %7.2f pid %1d   vz %+6d  dz %+6d (cut %5d), fid %1d -> fail\n",
                  i,
                  pfch[i].floatPt(),
                  pfch[i].intId(),
                  int(pfch[i].hwZ0),
                  z0diff,
                  dzCut_,
                  region.isFiducial(pfch[i]));
    }
  }
}

unsigned int l1ct::LinPuppiEmulator::find_ieta(const PFRegionEmu &region, eta_t eta) const {
  int n = absEtaBins_.size();
  glbeta_t abseta = region.hwGlbEta(eta);
  if (abseta < 0)
    abseta = -abseta;
  for (int i = 0; i < n; ++i) {
    if (abseta <= absEtaBins_[i])
      return i;
  }
  return n;
}

std::pair<pt_t, puppiWgt_t> l1ct::LinPuppiEmulator::sum2puppiPt_ref(
    sumTerm_t sum, pt_t pt, unsigned int ieta, bool isEM, int icand) const {
  const alpha_t logOffset = std::log2(linpuppi::PT2DR2_LSB) - SUM_BITSHIFT;
  const ptSlope_t ptSlopeNe = ptSlopeNe_[ieta];
  const ptSlope_t ptSlopePh = ptSlopePh_[ieta];
  const l1ct::pt_t ptZeroNe = ptZeroNe_[ieta];
  const l1ct::pt_t ptZeroPh = ptZeroPh_[ieta];
  const x2_t alphaCrop = alphaCrop_[ieta];
  // we put a log(2) in all alphaSlopes here since we compute alpha as log2(sum) instead of ln(sum)
  const alphaSlope_t alphaSlopeNe = alphaSlope_[ieta] * std::log(2.);
  const alphaSlope_t alphaSlopePh = alphaSlope_[ieta] * std::log(2.);
  const alpha_t alphaZeroNe = alphaZero_[ieta] / std::log(2.);
  const alpha_t alphaZeroPh = alphaZero_[ieta] / std::log(2.);
  const x2_t priorNe = priorNe_[ieta];
  const x2_t priorPh = priorPh_[ieta];

  // emulate computing
  //    alpha = log2(sum)
  //    x2a = alphaSlope*(alpha - alphaZero)
  // we use a 10-bit LUT for the log2(sum), and to save firmware resources we
  // also pack the computation of x2a in the same LUT.
  // In this emulator, however, we also compute the separately alpha, as it is
  // useful for debugging and comparison to the floating-point emulator.
  const int log2lut_bits = 10;
  alpha_t alpha = 0;
  uint64_t logarg = sum.bits_to_uint64();
  if (logarg > 0) {
    alpha = logOffset;
    while (logarg >= (1 << log2lut_bits)) {
      logarg = logarg >> 1;
      alpha += 1;
    }
    alpha += alpha_t(std::log2(float(logarg)));
  }
  alpha_t alphaZero = (isEM ? alphaZeroPh : alphaZeroNe);
  alphaSlope_t alphaSlope = (isEM ? alphaSlopePh : alphaSlopeNe);
  x2_t x2a = -alphaSlope * alphaZero;
  logarg = sum.bits_to_uint64();
  if (logarg > 0) {
    x2a += alphaSlope * logOffset;
    while (logarg >= (1 << log2lut_bits)) {
      logarg = logarg >> 1;
      x2a += alphaSlope;
    }
    x2a += alphaSlope * alpha_t(std::log2(float(logarg)));
  }
  /* // debug printout, can be useful to study ranges and precisions for LUTs and coefficients
    dbgPrintf("ref:  x2a(sum = %10.5f, raw = %9lu): logarg = %9lu, sumterm = %.6f, table[logarg] = %d = %.6f, ret pre-crop = %d = %.6f\n", 
          sum.to_float(), sum.bits_to_uint64(), logarg, 
          (alphaSlope * (logOffset - alphaZero)).to_float(),
          (alphaSlope * alpha_t(std::log2(float(logarg)))).V.to_int(), 
          (alphaSlope * alpha_t(std::log2(float(logarg)))).to_float(), 
          x2a.V.to_int(), x2a.to_float());
  */
  x2a = std::min(std::max<x2_t>(x2a, -alphaCrop), alphaCrop);

  l1ct::pt_t ptZero = (isEM ? ptZeroPh : ptZeroNe);
  ptSlope_t ptSlope = (isEM ? ptSlopePh : ptSlopeNe);
  x2_t x2pt = ptSlope * (pt - ptZero);

  x2_t prior = (isEM ? priorPh : priorNe);

  x2_t x2 = x2a + x2pt - prior;

  puppiWgt_t weight = 1.0 / (1.0 + std::exp(-x2.to_float()));
  typedef ap_fixed<pt_t::width, pt_t::iwidth, AP_RND, AP_SAT> pt_rounding_t;
  pt_t ptPuppi = pt_rounding_t(pt * weight);

  if (debug_)
    dbgPrintf(
        "ref candidate %02d pt %7.2f  em %1d  ieta %1d: sum %10.4f alpha %+7.2f   x2a %+7.3f  x2pt %+7.3f   x2 %+7.3f  "
        "--> weight %.4f  puppi pt %7.2f\n",
        icand,
        Scales::floatPt(pt),
        int(isEM),
        ieta,
        sum.to_float(),
        alpha.to_float() * std::log(2.),
        x2a.to_float(),
        x2pt.to_float(),
        x2.to_float(),
        weight.to_float(),
        Scales::floatPt(ptPuppi));

  return std::make_pair(ptPuppi, weight);
}

void l1ct::LinPuppiEmulator::fwdlinpuppi_ref(const PFRegionEmu &region,
                                             const std::vector<HadCaloObjEmu> &caloin /*[nIn]*/,
                                             std::vector<PuppiObjEmu> &outallne_nocut /*[nIn]*/,
                                             std::vector<PuppiObjEmu> &outallne /*[nIn]*/,
                                             std::vector<PuppiObjEmu> &outselne /*[nOut]*/) const {
  const unsigned int nIn = std::min<unsigned int>(nIn_, caloin.size());
  const unsigned int PTMAX2 = (iptMax_ * iptMax_);

  outallne_nocut.resize(nIn);
  outallne.resize(nIn);
  for (unsigned int in = 0; in < nIn; ++in) {
    outallne_nocut[in].clear();
    outallne[in].clear();
    if (caloin[in].hwPt == 0)
      continue;
    sumTerm_t sum = 0;  // (pt^2)/(dr2)
    for (unsigned int it = 0; it < nIn; ++it) {
      if (it == in || caloin[it].hwPt == 0)
        continue;
      unsigned int dr2 = dr2_int(caloin[it].hwEta, caloin[it].hwPhi, caloin[in].hwEta, caloin[in].hwPhi);
      if (dr2 <= dR2Max_) {                                             // if dr is inside puppi cone
        unsigned int dr2short = (dr2 >= dR2Min_ ? dr2 : dR2Min_) >> 5;  // reduce precision to make divide LUT cheaper
        unsigned int pt = caloin[it].intPt(), pt2 = pt * pt;
        dr2inv_t dr2inv = 1.0f / float(dr2short);
        sumTerm_t term = std::min(pt2 >> 5, PTMAX2 >> 5) * dr2inv;
        sum += term;
      }
    }
    unsigned int ieta = find_ieta(region, caloin[in].hwEta);
    std::pair<pt_t, puppiWgt_t> ptAndW = sum2puppiPt_ref(sum, caloin[in].hwPt, ieta, caloin[in].hwIsEM(), in);

    outallne_nocut[in].fill(region, caloin[in], ptAndW.first, ptAndW.second);
    if (region.isFiducial(caloin[in]) && outallne_nocut[in].hwPt >= ptCut_[ieta]) {
      outallne[in] = outallne_nocut[in];
    }
  }
  puppisort_and_crop_ref(nOut_, outallne, outselne);
}

void l1ct::LinPuppiEmulator::linpuppi_ref(const PFRegionEmu &region,
                                          const std::vector<TkObjEmu> &track /*[nTrack]*/,
                                          const std::vector<PVObjEmu> &pv, /*[nVtx]*/
                                          const std::vector<PFNeutralObjEmu> &pfallne /*[nIn]*/,
                                          std::vector<PuppiObjEmu> &outallne_nocut /*[nIn]*/,
                                          std::vector<PuppiObjEmu> &outallne /*[nIn]*/,
                                          std::vector<PuppiObjEmu> &outselne /*[nOut]*/) const {
  const unsigned int nIn = std::min<unsigned>(nIn_, pfallne.size());
  const unsigned int nTrack = std::min<unsigned int>(nTrack_, track.size());
  const unsigned int PTMAX2 = (iptMax_ * iptMax_);

  outallne_nocut.resize(nIn);
  outallne.resize(nIn);
  for (unsigned int in = 0; in < nIn; ++in) {
    outallne_nocut[in].clear();
    outallne[in].clear();
    if (pfallne[in].hwPt == 0)
      continue;
    sumTerm_t sum = 0;  // (pt^2)/(dr2)
    for (unsigned int it = 0; it < nTrack; ++it) {
      if (track[it].hwPt == 0)
        continue;

      int pZMin = 99999;
      bool pass_network = false;
      for (unsigned int v = 0; v < nVtx_; ++v) {
        if (v < pv.size()) {
          int ppZMin = std::abs(int(track[it].hwZ0 - pv[v].hwZ0));
          if (pZMin > ppZMin)
            pZMin = ppZMin;
          if (useMLAssociation_ and withinCMSSW_ &&
              nnVtxAssoc_->TTTrackNetworkSelector<const l1ct::TkObjEmu>(region, track[it], pv[v]) == 1)
            pass_network = true;
        }
      }
      if (useMLAssociation_ && !pass_network)
        continue;
      if (!useMLAssociation_ && std::abs(pZMin) > int(dzCut_))
        continue;
      unsigned int dr2 = dr2_int(pfallne[in].hwEta, pfallne[in].hwPhi, track[it].hwEta, track[it].hwPhi);
      if (dr2 <= dR2Max_) {                                             // if dr is inside puppi cone
        unsigned int dr2short = (dr2 >= dR2Min_ ? dr2 : dR2Min_) >> 5;  // reduce precision to make divide LUT cheaper
        unsigned int pt = track[it].intPt(), pt2 = pt * pt;
        dr2inv_t dr2inv = 1.0f / float(dr2short);
        sumTerm_t term = std::min(pt2 >> 5, PTMAX2 >> 5) * dr2inv;
        /* // printout useful for comparing internals steps of computation to the ones done in the firmware code
        dbgPrintf("cand pT %5.1f eta %+6.2f ref term %2d %2d: dr = %8d  pt2_shift = %8lu  term = %12lu  apterm = %12.6f\n",
                  pfallne[it].floatPt(),
                  region.floatGlbEtaOf(pfallne[it]),
                  in,
                  it,
                  dr2,
                  std::min<uint64_t>(pt2 >> 5, PTMAX2 >> 5),
                  term.to_float());
        */
        sum += term;
      }
    }

    unsigned int ieta = find_ieta(region, pfallne[in].hwEta);
    bool isEM = (pfallne[in].hwId.isPhoton());
    std::pair<pt_t, puppiWgt_t> ptAndW = sum2puppiPt_ref(sum, pfallne[in].hwPt, ieta, isEM, in);
    if (!fakePuppi_) {
      outallne_nocut[in].fill(region, pfallne[in], ptAndW.first, ptAndW.second);
      if (region.isFiducial(pfallne[in]) && outallne_nocut[in].hwPt >= ptCut_[ieta]) {
        outallne[in] = outallne_nocut[in];
      }
    } else {  // fakePuppi: keep the full candidate, but set the Puppi weight and some debug info into it
      outallne_nocut[in].fill(region, pfallne[in], pfallne[in].hwPt, ptAndW.second);
      outallne_nocut[in].hwData[9] = region.isFiducial(pfallne[in]);
      outallne_nocut[in].hwData(20, 10) = ptAndW.first(10, 0);
      outallne[in] = outallne_nocut[in];
    }
    if (debug_ && pfallne[in].hwPt > 0 && outallne_nocut[in].hwPt > 0) {
      dbgPrintf("ref candidate %02u pt %7.2f  -> puppi pt %7.2f, fiducial %1d, packed %s\n",
                in,
                pfallne[in].floatPt(),
                outallne_nocut[in].floatPt(),
                int(region.isFiducial(pfallne[in])),
                outallne_nocut[in].pack().to_string(16).c_str());
    }
  }
  puppisort_and_crop_ref(nOut_, outallne, outselne);
}

std::pair<float, float> l1ct::LinPuppiEmulator::sum2puppiPt_flt(
    float sum, float pt, unsigned int ieta, bool isEM, int icand) const {
  float alphaZero = alphaZero_[ieta], alphaSlope = alphaSlope_[ieta], alphaCrop = alphaCrop_[ieta];
  float alpha = sum > 0 ? std::log(sum) : -9e9;
  float x2a = std::min(std::max(alphaSlope * (alpha - alphaZero), -alphaCrop), alphaCrop);

  float ptZero = (isEM ? ptZeroPh_[ieta] : ptZeroNe_[ieta]);
  float ptSlope = (isEM ? ptSlopePh_[ieta] : ptSlopeNe_[ieta]);
  float x2pt = ptSlope * (pt - ptZero);

  float prior = (isEM ? priorPh_[ieta] : priorNe_[ieta]);

  float x2 = x2a + x2pt - prior;

  float weight = 1.0 / (1.0 + std::exp(-x2));

  float puppiPt = pt * weight;
  if (debug_)
    dbgPrintf(
        "flt candidate %02d pt %7.2f  em %1d  ieta %1d: alpha %+7.2f   x2a         %+7.3f  x2pt         %+7.3f   x2    "
        "     %+7.3f  --> weight        %.4f  puppi pt %7.2f\n",
        icand,
        pt,
        int(isEM),
        ieta,
        std::max(alpha, -99.99f),
        x2a,
        x2pt,
        x2,
        weight,
        puppiPt);

  return std::make_pair(puppiPt, weight);
}

void l1ct::LinPuppiEmulator::fwdlinpuppi_flt(const PFRegionEmu &region,
                                             const std::vector<HadCaloObjEmu> &caloin /*[nIn]*/,
                                             std::vector<PuppiObjEmu> &outallne_nocut /*[nIn]*/,
                                             std::vector<PuppiObjEmu> &outallne /*[nIn]*/,
                                             std::vector<PuppiObjEmu> &outselne /*[nOut]*/) const {
  const unsigned int nIn = std::min<unsigned int>(nIn_, caloin.size());
  const float f_ptMax = Scales::floatPt(Scales::makePt(iptMax_));

  outallne_nocut.resize(nIn);
  outallne.resize(nIn);
  for (unsigned int in = 0; in < nIn; ++in) {
    outallne_nocut[in].clear();
    outallne[in].clear();
    if (caloin[in].hwPt == 0)
      continue;
    float sum = 0;
    for (unsigned int it = 0; it < nIn; ++it) {
      if (it == in || caloin[it].hwPt == 0)
        continue;
      unsigned int dr2 = dr2_int(caloin[it].hwEta, caloin[it].hwPhi, caloin[in].hwEta, caloin[in].hwPhi);
      if (dr2 <= dR2Max_) {  // if dr is inside puppi cone
        sum += std::pow(std::min(caloin[it].floatPt(), f_ptMax), 2) / (std::max(dr2, dR2Min_) * DR2_LSB);
      }
    }

    unsigned int ieta = find_ieta(region, caloin[in].hwEta);
    std::pair<float, float> ptAndW = sum2puppiPt_flt(sum, caloin[in].floatPt(), ieta, caloin[in].hwIsEM(), in);
    outallne_nocut[in].fill(region, caloin[in], Scales::makePtFromFloat(ptAndW.first), l1ct::puppiWgt_t(ptAndW.second));
    if (region.isFiducial(caloin[in]) && outallne_nocut[in].hwPt >= ptCut_[ieta]) {
      outallne[in] = outallne_nocut[in];
    }
  }

  puppisort_and_crop_ref(nOut_, outallne, outselne);
}

void l1ct::LinPuppiEmulator::linpuppi_flt(const PFRegionEmu &region,
                                          const std::vector<TkObjEmu> &track /*[nTrack]*/,
                                          const std::vector<PVObjEmu> &pv,
                                          const std::vector<PFNeutralObjEmu> &pfallne /*[nIn]*/,
                                          std::vector<PuppiObjEmu> &outallne_nocut /*[nIn]*/,
                                          std::vector<PuppiObjEmu> &outallne /*[nIn]*/,
                                          std::vector<PuppiObjEmu> &outselne /*[nOut]*/) const {
  const unsigned int nIn = std::min<unsigned>(nIn_, pfallne.size());
  const unsigned int nTrack = std::min<unsigned int>(nTrack_, track.size());
  const float f_ptMax = Scales::floatPt(Scales::makePt(iptMax_));

  outallne_nocut.resize(nIn);
  outallne.resize(nIn);
  for (unsigned int in = 0; in < nIn; ++in) {
    outallne_nocut[in].clear();
    outallne[in].clear();
    if (pfallne[in].hwPt == 0)
      continue;
    float sum = 0;
    for (unsigned int it = 0; it < nTrack; ++it) {
      if (track[it].hwPt == 0)
        continue;

      int pZMin = 99999;
      bool pass_network = false;
      for (unsigned int v = 0, nVtx = std::min<unsigned int>(nVtx_, pv.size()); v < nVtx; ++v) {
        int ppZMin = std::abs(int(track[it].hwZ0 - pv[v].hwZ0));
        if (pZMin > ppZMin)
          pZMin = ppZMin;
        if (useMLAssociation_ and withinCMSSW_ &&
            nnVtxAssoc_->TTTrackNetworkSelector<const l1ct::TkObjEmu>(region, track[it], pv[v]) == 1)
          pass_network = true;
      }
      if (useMLAssociation_ && !pass_network)
        continue;
      if (!useMLAssociation_ && std::abs(pZMin) > int(dzCut_))
        continue;
      unsigned int dr2 = dr2_int(
          pfallne[in].hwEta, pfallne[in].hwPhi, track[it].hwEta, track[it].hwPhi);  // if dr is inside puppi cone
      if (dr2 <= dR2Max_) {
        sum += std::pow(std::min(track[it].floatPt(), f_ptMax), 2) / (std::max(dr2, dR2Min_) * DR2_LSB);
      }
    }
    unsigned int ieta = find_ieta(region, pfallne[in].hwEta);
    bool isEM = pfallne[in].hwId.isPhoton();
    std::pair<float, float> ptAndW = sum2puppiPt_flt(sum, pfallne[in].floatPt(), ieta, isEM, in);
    outallne_nocut[in].fill(
        region, pfallne[in], Scales::makePtFromFloat(ptAndW.first), l1ct::puppiWgt_t(ptAndW.second));
    if (region.isFiducial(pfallne[in]) && outallne_nocut[in].hwPt >= ptCut_[ieta]) {
      outallne[in] = outallne_nocut[in];
    }
  }
  puppisort_and_crop_ref(nOut_, outallne, outselne);
}

void l1ct::LinPuppiEmulator::run(const PFInputRegion &in,
                                 const std::vector<l1ct::PVObjEmu> &pvs,
                                 OutputRegion &out) const {
  if (debug_) {
    dbgPrintf("\nWill run LinPuppi in region eta %+5.2f, phi %+5.2f, pv0 int Z %+d\n",
              in.region.floatEtaCenter(),
              in.region.floatPhiCenter(),
              pvs.front().hwZ0.to_int());
  }
  if (std::abs(in.region.floatEtaCenter()) < 2.5) {  // within tracker
    std::vector<PuppiObjEmu> outallch, outallne_nocut, outallne, outselne;
    linpuppi_chs_ref(in.region, pvs, out.pfcharged, outallch);
    linpuppi_ref(in.region, in.track, pvs, out.pfneutral, outallne_nocut, outallne, outselne);
    // ensure proper sizes of the vectors, to get accurate sorting wrt firmware
    const std::vector<PuppiObjEmu> &ne = (nOut_ == nIn_ ? outallne : outselne);
    unsigned int nch = outallch.size(), nne = ne.size(), i;
    outallch.resize(nTrack_ + nOut_);
    for (i = nch; i < nTrack_; ++i)
      outallch[i].clear();
    for (unsigned int j = 0; j < nne; ++i, ++j)
      outallch[i] = ne[j];
    for (; i < nTrack_ + nOut_; ++i)
      outallch[i].clear();
    puppisort_and_crop_ref(nFinalSort_, outallch, out.puppi, finalSortAlgo_);
    // trim if needed
    while (!out.puppi.empty() && out.puppi.back().hwPt == 0)
      out.puppi.pop_back();
    out.puppi.shrink_to_fit();
  } else {  // forward
    std::vector<PuppiObjEmu> outallne_nocut, outallne;
    fwdlinpuppi_ref(in.region, in.hadcalo, outallne_nocut, outallne, out.puppi);
  }
}

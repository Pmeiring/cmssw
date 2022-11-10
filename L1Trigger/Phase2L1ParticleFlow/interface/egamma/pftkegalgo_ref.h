#ifndef PFTKEGALGO_REF_H
#define PFTKEGALGO_REF_H

#include "DataFormats/L1TParticleFlow/interface/layer1_emulator.h"
#include "DataFormats/L1TParticleFlow/interface/egamma.h"
#include "DataFormats/L1TParticleFlow/interface/pf.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/conifer.h"

namespace edm {
  class ParameterSet;
}

namespace l1ct {

  struct PFTkEGAlgoEmuConfig {
    unsigned int nTRACK;
    unsigned int nTRACK_EGIN;
    unsigned int nEMCALO_EGIN;
    unsigned int nEM_EGOUT;

    bool filterHwQuality;
    bool doBremRecovery;
    bool writeBeforeBremRecovery;
    int caloHwQual;
    bool doEndcapHwQual;
    float emClusterPtMin;  // GeV
    float dEtaMaxBrem;
    float dPhiMaxBrem;

    std::vector<double> absEtaBoundaries;
    std::vector<double> dEtaValues;
    std::vector<double> dPhiValues;
    float trkQualityPtMin;  // GeV
    bool writeEgSta;

    struct IsoParameters {
      IsoParameters(const edm::ParameterSet &);
      IsoParameters(float tkQualityPtMin, float dZ, float dRMin, float dRMax)
          : tkQualityPtMin(Scales::makePtFromFloat(tkQualityPtMin)),
            dZ(Scales::makeZ0(dZ)),
            dRMin2(Scales::makeDR2FromFloatDR(dRMin)),
            dRMax2(Scales::makeDR2FromFloatDR(dRMax)) {}
      pt_t tkQualityPtMin;
      ap_int<z0_t::width + 1> dZ;
      int dRMin2;
      int dRMax2;
    };

    IsoParameters tkIsoParams_tkEle;
    IsoParameters tkIsoParams_tkEm;
    IsoParameters pfIsoParams_tkEle;
    IsoParameters pfIsoParams_tkEm;
    bool doTkIso;
    bool doPfIso;
    EGIsoEleObjEmu::IsoType hwIsoTypeTkEle;
    EGIsoObjEmu::IsoType hwIsoTypeTkEm;

    bool doCompositeTkEle;
    struct CompIDParameters {
      CompIDParameters(const edm::ParameterSet &);
      CompIDParameters(double hoeMin_highpt, double hoeMax_highpt, double tkptMin_highpt, double tkptMax_highpt, double srrtotMin_highpt, double srrtotMax_highpt, double detaMin_highpt, double detaMax_highpt, double dptMin_highpt, double dptMax_highpt, double meanzMin_highpt, double meanzMax_highpt, double dphiMin_highpt, double dphiMax_highpt, double tkchi2Min_highpt, double tkchi2Max_highpt, double tkz0Min_highpt, double tkz0Max_highpt, double tknstubsMin_highpt, double tknstubsMax_highpt, double BDTcut_wp97p5, double BDTcut_wp95p0,
                       double coreShowerLengthMin_lowpt, double coreShowerLengthMax_lowpt, double tkptMin_lowpt, double tkptMax_lowpt, double srrtotMin_lowpt, double srrtotMax_lowpt, double detaMin_lowpt, double detaMax_lowpt, double dptMin_lowpt, double dptMax_lowpt, double meanzMin_lowpt, double meanzMax_lowpt, double dphiMin_lowpt, double dphiMax_lowpt, double tkchi2Min_lowpt, double tkchi2Max_lowpt, double tkz0Min_lowpt, double tkz0Max_lowpt, double tknstubsMin_lowpt, double tknstubsMax_lowpt)
          : hoeMin_highpt(hoeMin_highpt),           hoeMax_highpt(hoeMax_highpt),
            tkptMin_highpt(tkptMin_highpt),         tkptMax_highpt(tkptMax_highpt),
            srrtotMin_highpt(srrtotMin_highpt),     srrtotMax_highpt(srrtotMax_highpt),
            detaMin_highpt(detaMin_highpt),         detaMax_highpt(detaMax_highpt),
            dptMin_highpt(dptMin_highpt),           dptMax_highpt(dptMax_highpt),
            meanzMin_highpt(meanzMin_highpt),       meanzMax_highpt(meanzMax_highpt),
            dphiMin_highpt(dphiMin_highpt),         dphiMax_highpt(dphiMax_highpt),
            tkchi2Min_highpt(tkchi2Min_highpt),     tkchi2Max_highpt(tkchi2Max_highpt),
            tkz0Min_highpt(tkz0Min_highpt),         tkz0Max_highpt(tkz0Max_highpt),
            tknstubsMin_highpt(tknstubsMin_highpt), tknstubsMax_highpt(tknstubsMax_highpt),
            BDTcut_wp97p5(BDTcut_wp97p5),           BDTcut_wp95p0(BDTcut_wp95p0),
            coreShowerLengthMin_lowpt(coreShowerLengthMin_lowpt),coreShowerLengthMax_lowpt(coreShowerLengthMax_lowpt),
            tkptMin_lowpt(tkptMin_lowpt),           tkptMax_lowpt(tkptMax_lowpt),
            srrtotMin_lowpt(srrtotMin_lowpt),       srrtotMax_lowpt(srrtotMax_lowpt),
            detaMin_lowpt(detaMin_lowpt),           detaMax_lowpt(detaMax_lowpt),
            dptMin_lowpt(dptMin_lowpt),             dptMax_lowpt(dptMax_lowpt),
            meanzMin_lowpt(meanzMin_lowpt),         meanzMax_lowpt(meanzMax_lowpt),
            dphiMin_lowpt(dphiMin_lowpt),           dphiMax_lowpt(dphiMax_lowpt),
            tkchi2Min_lowpt(tkchi2Min_lowpt),       tkchi2Max_lowpt(tkchi2Max_lowpt),
            tkz0Min_lowpt(tkz0Min_lowpt),           tkz0Max_lowpt(tkz0Max_lowpt),
            tknstubsMin_lowpt(tknstubsMin_lowpt),   tknstubsMax_lowpt(tknstubsMax_lowpt){}
      double hoeMin_highpt, hoeMax_highpt; // Parameters for the high pT model
      double tkptMin_highpt, tkptMax_highpt;
      double srrtotMin_highpt, srrtotMax_highpt;
      double detaMin_highpt, detaMax_highpt;
      double dptMin_highpt, dptMax_highpt;
      double meanzMin_highpt, meanzMax_highpt;
      double dphiMin_highpt, dphiMax_highpt;
      double tkchi2Min_highpt, tkchi2Max_highpt;
      double tkz0Min_highpt, tkz0Max_highpt;
      double tknstubsMin_highpt, tknstubsMax_highpt;
      double BDTcut_wp97p5, BDTcut_wp95p0;
      double coreShowerLengthMin_lowpt, coreShowerLengthMax_lowpt;  // Similar for the low pT model
      double tkptMin_lowpt, tkptMax_lowpt;
      double srrtotMin_lowpt, srrtotMax_lowpt;
      double detaMin_lowpt, detaMax_lowpt;
      double dptMin_lowpt, dptMax_lowpt;
      double meanzMin_lowpt, meanzMax_lowpt;
      double dphiMin_lowpt, dphiMax_lowpt;
      double tkchi2Min_lowpt, tkchi2Max_lowpt;
      double tkz0Min_lowpt, tkz0Max_lowpt;
      double tknstubsMin_lowpt, tknstubsMax_lowpt;
    };

    CompIDParameters myCompIDparams;

    int debug = 0;

    PFTkEGAlgoEmuConfig(const edm::ParameterSet &iConfig);
    PFTkEGAlgoEmuConfig(unsigned int nTrack,
                        unsigned int nTrack_in,
                        unsigned int nEmCalo_in,
                        unsigned int nEmOut,
                        bool filterHwQuality,
                        bool doBremRecovery,
                        bool writeBeforeBremRecovery = false,
                        int caloHwQual = 4,
                        bool doEndcapHwQual = false,
                        float emClusterPtMin = 2.,
                        float dEtaMaxBrem = 0.02,
                        float dPhiMaxBrem = 0.1,
                        const std::vector<double> &absEtaBoundaries = {0.0, 1.5},
                        const std::vector<double> &dEtaValues = {0.015, 0.01},
                        const std::vector<double> &dPhiValues = {0.07, 0.07},
                        float trkQualityPtMin = 10.,
                        bool writeEgSta = false,
                        const IsoParameters &tkIsoParams_tkEle = {2., 0.6, 0.03, 0.2},
                        const IsoParameters &tkIsoParams_tkEm = {2., 0.6, 0.07, 0.3},
                        const IsoParameters &pfIsoParams_tkEle = {1., 0.6, 0.03, 0.2},
                        const IsoParameters &pfIsoParams_tkEm = {1., 0.6, 0.07, 0.3},
                        bool doTkIso = true,
                        bool doPfIso = false,
                        EGIsoEleObjEmu::IsoType hwIsoTypeTkEle = EGIsoEleObjEmu::IsoType::TkIso,
                        EGIsoObjEmu::IsoType hwIsoTypeTkEm = EGIsoObjEmu::IsoType::TkIsoPV,
                        bool doCompositeTkEle = false,
                        const CompIDParameters &myCompIDparams = {-1.0, 1566.547607421875, 1.9501149654388428, 11102.0048828125, 0.0, 0.01274710614234209, -0.24224889278411865, 0.23079538345336914, 0.010325592942535877, 184.92538452148438, 325.0653991699219, 499.6089782714844, -6.281332015991211, 6.280326843261719, 0.024048099294304848, 1258.37158203125, -14.94140625, 14.94140625, 4.0, 6.0, 0.5406244, 0.9693441,
                                                                  1.0, 26.0, 1.9501149654388428, 11102.0048828125, 0.0, 0.013055052608251572, -0.24158358573913574, 0.2402184009552002, 0.07829763740301132, 1558.2740478515625, 322.10272216796875, 504.7358093261719, -6.280882835388184, 6.28144645690918, 0.00784884113818407, 1252.6024169921875, -14.94140625, 14.94140625, 4.0, 6.0},
                        int debug = 0)

        : nTRACK(nTrack),
          nTRACK_EGIN(nTrack_in),
          nEMCALO_EGIN(nEmCalo_in),
          nEM_EGOUT(nEmOut),
          filterHwQuality(filterHwQuality),
          doBremRecovery(doBremRecovery),
          writeBeforeBremRecovery(writeBeforeBremRecovery),
          caloHwQual(caloHwQual),
          doEndcapHwQual(doEndcapHwQual),
          emClusterPtMin(emClusterPtMin),
          dEtaMaxBrem(dEtaMaxBrem),
          dPhiMaxBrem(dPhiMaxBrem),
          //absEtaBoundaries(std::move(absEtaBoundaries)),
          //dEtaValues(std::move(dEtaValues)),
          //dPhiValues(std::move(dPhiValues)),
          absEtaBoundaries(absEtaBoundaries),
          dEtaValues(dEtaValues),
          dPhiValues(dPhiValues),
          trkQualityPtMin(trkQualityPtMin),
          writeEgSta(writeEgSta),
          tkIsoParams_tkEle(tkIsoParams_tkEle),
          tkIsoParams_tkEm(tkIsoParams_tkEm),
          pfIsoParams_tkEle(pfIsoParams_tkEle),
          pfIsoParams_tkEm(pfIsoParams_tkEm),
          doTkIso(doTkIso),
          doPfIso(doPfIso),
          hwIsoTypeTkEle(hwIsoTypeTkEle),
          hwIsoTypeTkEm(hwIsoTypeTkEm),
          doCompositeTkEle(doCompositeTkEle),
          myCompIDparams(myCompIDparams),
          debug(debug) {}
  };

  class PFTkEGAlgoEmulator {
  public:

    PFTkEGAlgoEmulator(const PFTkEGAlgoEmuConfig &config) : cfg(config), debug_(cfg.debug) {}


    virtual ~PFTkEGAlgoEmulator() {}

    void toFirmware(const PFInputRegion &in, PFRegion &region, EmCaloObj calo[/*nCALO*/], TkObj track[/*nTRACK*/]) const;
    void toFirmware(const OutputRegion &out, EGIsoObj out_egphs[], EGIsoEleObj out_egeles[]) const;
    void toFirmware(const PFInputRegion &in,
                    const l1ct::PVObjEmu &pvin,
                    PFRegion &region,
                    TkObj track[/*nTRACK*/],
                    PVObj &pv) const;

    void run(const PFInputRegion &in, OutputRegion &out) const;
    void runIso(const PFInputRegion &in, const std::vector<l1ct::PVObjEmu> &pvs, OutputRegion &out) const;

    void setDebug(int verbose) { debug_ = verbose; }

    bool writeEgSta() const { return cfg.writeEgSta; }

  private:


    void link_emCalo2emCalo(const std::vector<EmCaloObjEmu> &emcalo, std::vector<int> &emCalo2emCalo) const;

    void link_emCalo2tk_elliptic(const PFRegionEmu &r,
                                const std::vector<EmCaloObjEmu> &emcalo,
                                const std::vector<TkObjEmu> &track,
                                std::vector<int> &emCalo2tk) const;

    void link_emCalo2tk_composite(const PFRegionEmu &r,
                                const std::vector<EmCaloObjEmu> &emcalo,
                                const std::vector<TkObjEmu> &track,
                                std::vector<int> &emCalo2tk,
                                std::vector<float> &emCaloTkBdtScore_highpt,
                                std::vector<float> &emCaloTkBdtScore_lowpt) const;

    struct CompositeCandidate {
      unsigned int cluster_idx;
      unsigned int track_idx;
      double dpt; // For sorting
    };

    float compute_composite_score_highpt(CompositeCandidate &cand,
                                  const std::vector<EmCaloObjEmu> &emcalo,
                                  const std::vector<TkObjEmu> &track,
                                  const PFTkEGAlgoEmuConfig::CompIDParameters &params) const;

    float compute_composite_score_lowpt(CompositeCandidate &cand,
                                  const std::vector<EmCaloObjEmu> &emcalo,
                                  const std::vector<TkObjEmu> &track,
                                  const PFTkEGAlgoEmuConfig::CompIDParameters &params) const;

    //FIXME: still needed
    float deltaPhi(float phi1, float phi2) const;

    void sel_emCalo(unsigned int nmax_sel,
                    const std::vector<EmCaloObjEmu> &emcalo,
                    std::vector<EmCaloObjEmu> &emcalo_sel) const;

    void eg_algo(const PFRegionEmu &region,
                 const std::vector<EmCaloObjEmu> &emcalo,
                 const std::vector<TkObjEmu> &track,
                 const std::vector<int> &emCalo2emCalo,
                 const std::vector<int> &emCalo2tk,
                 const std::vector<float>&emCaloTkBdtScore_highpt,
                 const std::vector<float>&emCaloTkBdtScore_lowpt,
                 std::vector<EGObjEmu> &egstas,
                 std::vector<EGIsoObjEmu> &egobjs,
                 std::vector<EGIsoEleObjEmu> &egeleobjs) const;

    void addEgObjsToPF(std::vector<EGObjEmu> &egstas,
                       std::vector<EGIsoObjEmu> &egobjs,
                       std::vector<EGIsoEleObjEmu> &egeleobjs,
                       const std::vector<EmCaloObjEmu> &emcalo,
                       const std::vector<TkObjEmu> &track,
                       const int calo_idx,
                       const unsigned int hwQual,
                       const pt_t ptCorr,
                       const int tk_idx,
                       const float bdtScore_highpt,
                       const float bdtScore_lowpt,
                       const std::vector<unsigned int> &components = {}) const;

    EGObjEmu &addEGStaToPF(std::vector<EGObjEmu> &egobjs,
                           const EmCaloObjEmu &calo,
                           const unsigned int hwQual,
                           const pt_t ptCorr,
                           const std::vector<unsigned int> &components) const;

    EGIsoObjEmu &addEGIsoToPF(std::vector<EGIsoObjEmu> &egobjs,
                              const EmCaloObjEmu &calo,
                              const unsigned int hwQual,
                              const pt_t ptCorr) const;

    EGIsoEleObjEmu &addEGIsoEleToPF(std::vector<EGIsoEleObjEmu> &egobjs,
                                    const EmCaloObjEmu &calo,
                                    const TkObjEmu &track,
                                    const unsigned int hwQual,
                                    const pt_t ptCorr,
                                    const float bdtScore_highpt,
                                    const float bdtScore_lowpt) const;

    // FIXME: reimplemented from PFAlgoEmulatorBase
    template <typename T>
    void ptsort_ref(int nIn, int nOut, const std::vector<T> &in, std::vector<T> &out) const {
      // // std based sort
      // out = in;
      // std::sort(out.begin(), out.end(),  std::greater<T>());
      // out.resize(nOut);
      out.resize(nOut);
      for (int iout = 0; iout < nOut; ++iout) {
        out[iout].clear();
      }
      for (int it = 0; it < nIn; ++it) {
        for (int iout = 0; iout < nOut; ++iout) {
          if (in[it].hwPt >= out[iout].hwPt) {
            for (int i2 = nOut - 1; i2 > iout; --i2) {
              out[i2] = out[i2 - 1];
            }
            out[iout] = in[it];
            break;
          }
        }
      }
    }

    template <typename T>
    int deltaR2(const T &charged, const EGIsoObjEmu &egphoton) const {
      // NOTE: we compare Tk at vertex against the calo variable...
      return dr2_int(charged.hwVtxEta(), charged.hwVtxPhi(), egphoton.hwEta, egphoton.hwPhi);
    }

    template <typename T>
    int deltaR2(const T &charged, const EGIsoEleObjEmu &egele) const {
      return dr2_int(charged.hwVtxEta(), charged.hwVtxPhi(), egele.hwVtxEta(), egele.hwVtxPhi());
    }

    int deltaR2(const PFNeutralObjEmu &neutral, const EGIsoObjEmu &egphoton) const {
      return dr2_int(neutral.hwEta, neutral.hwPhi, egphoton.hwEta, egphoton.hwPhi);
    }

    int deltaR2(const PFNeutralObjEmu &neutral, const EGIsoEleObjEmu &egele) const {
      // NOTE: we compare Tk at vertex against the calo variable...
      return dr2_int(neutral.hwEta, neutral.hwPhi, egele.hwVtxEta(), egele.hwVtxPhi());
    }

    template <typename T>
    ap_int<z0_t::width + 1> deltaZ0(const T &charged, const EGIsoObjEmu &egphoton, z0_t z0) const {
      ap_int<z0_t::width + 1> delta = charged.hwZ0 - z0;
      if (delta < 0)
        delta = -delta;
      return delta;
    }

    template <typename T>
    ap_int<z0_t::width + 1> deltaZ0(const T &charged, const EGIsoEleObjEmu &egele, z0_t z0) const {
      ap_int<z0_t::width + 1> delta = charged.hwZ0 - egele.hwZ0;
      if (delta < 0)
        delta = -delta;
      return delta;
    }

    template <typename TCH, typename TEG>
    void compute_sumPt(iso_t &sumPt,
                       iso_t &sumPtPV,
                       const std::vector<TCH> &objects,
                       unsigned int nMaxObj,
                       const TEG &egobj,
                       const PFTkEGAlgoEmuConfig::IsoParameters &params,
                       z0_t z0) const {
      for (unsigned int itk = 0; itk < std::min<unsigned>(objects.size(), nMaxObj); ++itk) {
        const auto &obj = objects[itk];

        if (obj.hwPt < params.tkQualityPtMin)
          continue;

        int dR2 = deltaR2(obj, egobj);

        if (dR2 > params.dRMin2 && dR2 < params.dRMax2) {
          sumPt += obj.hwPt;
          if (deltaZ0(obj, egobj, z0) < params.dZ) {
            sumPtPV += obj.hwPt;
          }
        }
      }
    }

    template <typename TEG>
    void compute_sumPt(iso_t &sumPt,
                       iso_t &sumPtPV,
                       const std::vector<PFNeutralObjEmu> &objects,
                       unsigned int nMaxObj,
                       const TEG &egobj,
                       const PFTkEGAlgoEmuConfig::IsoParameters &params,
                       z0_t z0) const {
      for (unsigned int itk = 0; itk < std::min<unsigned>(objects.size(), nMaxObj); ++itk) {
        const auto &obj = objects[itk];

        if (obj.hwPt < params.tkQualityPtMin)
          continue;

        int dR2 = deltaR2(obj, egobj);

        if (dR2 > params.dRMin2 && dR2 < params.dRMax2) {
          sumPt += obj.hwPt;
          // PF neutrals are not constrained by PV (since their Z0 is 0 by design)
          sumPtPV += obj.hwPt;
        }
      }
    }

    void compute_isolation(std::vector<EGIsoObjEmu> &egobjs,
                           const std::vector<TkObjEmu> &objects,
                           const PFTkEGAlgoEmuConfig::IsoParameters &params,
                           z0_t z0) const;
    void compute_isolation(std::vector<EGIsoEleObjEmu> &egobjs,
                           const std::vector<TkObjEmu> &objects,
                           const PFTkEGAlgoEmuConfig::IsoParameters &params,
                           z0_t z0) const;
    void compute_isolation(std::vector<EGIsoObjEmu> &egobjs,
                           const std::vector<PFChargedObjEmu> &charged,
                           const std::vector<PFNeutralObjEmu> &neutrals,
                           const PFTkEGAlgoEmuConfig::IsoParameters &params,
                           z0_t z0) const;
    void compute_isolation(std::vector<EGIsoEleObjEmu> &egobjs,
                           const std::vector<PFChargedObjEmu> &charged,
                           const std::vector<PFNeutralObjEmu> &neutrals,
                           const PFTkEGAlgoEmuConfig::IsoParameters &params,
                           z0_t z0) const;

    PFTkEGAlgoEmuConfig cfg;
    std::unique_ptr<conifer::BDT<ap_fixed<22,3,AP_RND_CONV,AP_SAT>,ap_fixed<22,3,AP_RND_CONV,AP_SAT>,0>> composite_bdt_highpt_;
    std::unique_ptr<conifer::BDT<ap_fixed<22,3,AP_RND_CONV,AP_SAT>,ap_fixed<22,3,AP_RND_CONV,AP_SAT>,0>> composite_bdt_lowpt_;
    int debug_;
  };
}  // namespace l1ct

#endif

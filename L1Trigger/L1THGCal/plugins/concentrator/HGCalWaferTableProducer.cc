#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "L1Trigger/L1THGCal/interface/HGCalProcessorBase.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerTools.h"
#include "L1Trigger/L1THGCal/interface/concentrator/AEinputUtil.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

class HGCalWaferTableProducer : public edm::stream::EDProducer<>{
public:
  explicit HGCalWaferTableProducer(const edm::ParameterSet&);
  void beginRun(const edm::Run&, const edm::EventSetup&) override;
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  double calibrate(double, int, unsigned);
  void simhits(const edm::Event& e,
               std::unordered_map<uint32_t, double>& simhits_ee,
               std::unordered_map<uint32_t, double>& simhits_fh,
               std::unordered_map<uint32_t, double>& simhits_bh);
  float getsimenergy(const unsigned& i, 
        const std::unordered_map<uint32_t, double>& simhits_ee, 
        const std::unordered_map<uint32_t, double>& simhits_fh,
        const std::unordered_map<uint32_t, double>& simhits_bh);

  edm::ESHandle<HGCalTriggerGeometryBase> triggerGeometry_;
  const edm::ESGetToken<HGCalTriggerGeometryBase, CaloGeometryRecord> triggerGeomToken_;
  edm::EDGetTokenT<l1t::HGCalTriggerCellBxCollection> trigger_cells_token_;
  edm::EDGetTokenT<edm::PCaloHitContainer> simhits_ee_token_, simhits_fh_token_, simhits_bh_token_;
  HGCalTriggerTools triggerTools_;
  AEinputUtil aeInputUtil_;

  double keV2fC_;
  std::vector<double> fcPerMip_;
  std::vector<double> layerWeights_;
  std::vector<double> thicknessCorrections_;

  static constexpr int nInputs_ = 64;

  int wafer_n_;
  std::vector<uint32_t> wafer_id_;
  std::vector<int> wafer_side_;
  std::vector<int> wafer_layer_;
  std::vector<int> wafer_subdet_;
  std::vector<int> wafer_waferu_;
  std::vector<int> wafer_waferv_;
  std::vector<int> wafer_wafertype_;
  std::vector<float> wafer_mipPt_;
  std::vector<float> wafer_energy_;
  std::vector<float> wafer_simenergy_;
  std::vector<float> wafer_eta_;
  std::vector<float> wafer_phi_;
  std::vector<float> wafer_x_;
  std::vector<float> wafer_y_;
  std::vector<float> wafer_z_;
  std::vector<uint64_t> wafer_sumCALQ_;
  std::vector<std::vector<uint32_t>> wafer_norm_{nInputs_};
  std::vector<std::vector<uint32_t>> wafer_ADC_{nInputs_};
  std::vector<std::vector<uint32_t>> wafer_CALQ_{nInputs_};
  std::vector<std::vector<float>> wafer_AEin_{nInputs_};

};

HGCalWaferTableProducer::HGCalWaferTableProducer(const edm::ParameterSet& conf)
  : triggerGeomToken_(esConsumes<HGCalTriggerGeometryBase, CaloGeometryRecord, edm::Transition::BeginRun>()),
    trigger_cells_token_(consumes<l1t::HGCalTriggerCellBxCollection>(conf.getParameter<edm::InputTag>("TriggerCells"))),
    simhits_ee_token_(consumes<edm::PCaloHitContainer>(conf.getParameter<edm::InputTag>("eeSimHits"))),
    simhits_fh_token_(consumes<edm::PCaloHitContainer>(conf.getParameter<edm::InputTag>("fhSimHits"))),
    simhits_bh_token_(consumes<edm::PCaloHitContainer>(conf.getParameter<edm::InputTag>("bhSimHits"))),
    aeInputUtil_(conf.getParameter<unsigned>("bitsPerADC"),
                 conf.getParameter<unsigned>("bitsPerNorm"),
                 conf.getParameter<unsigned>("bitsPerCALQ"),
                 conf.getParameter<unsigned>("bitsPerInput"),
                 conf.getParameter<bool>("useModuleFactor"),
                 conf.getParameter<bool>("bitShiftNormalize"),
                 conf.getParameter<bool>("useTransverseADC"),
                 conf.getParameter<bool>("normByMax")),
    keV2fC_(conf.getParameter<double>("keV2fC")),
    fcPerMip_(conf.getParameter<std::vector<double>>("fcPerMip")),
    layerWeights_(conf.getParameter<std::vector<double>>("layerWeights")),
    thicknessCorrections_(conf.getParameter<std::vector<double>>("thicknessCorrections"))
{
  produces<nanoaod::FlatTable>("wafers");
}

void HGCalWaferTableProducer::beginRun(const edm::Run& /*run*/, const edm::EventSetup& es) {
  triggerGeometry_ = es.getHandle(triggerGeomToken_);
}


float HGCalWaferTableProducer::getsimenergy(const unsigned& id, 
        const std::unordered_map<uint32_t, double>& simhits_ee, 
        const std::unordered_map<uint32_t, double>& simhits_fh,
        const std::unordered_map<uint32_t, double>& simhits_bh) {
    double energy = 0;
    unsigned layer = triggerTools_.layerWithOffset(id);
    // search for simhit for all the cells inside the trigger cell
    const auto& cells = triggerTools_.getTriggerGeometry()->getCellsFromTriggerCell(id);
    for (uint32_t c_id : cells) {
      int thickness = triggerTools_.thicknessIndex(c_id);
      if (triggerTools_.isEm(id)) {
        auto itr = simhits_ee.find(c_id);
        if (itr != simhits_ee.end()){
          energy += calibrate(itr->second, thickness, layer);
        }
      } else if (triggerTools_.isSilicon(id)) {
        auto itr = simhits_fh.find(c_id);
        if (itr != simhits_fh.end()){
          energy += calibrate(itr->second, thickness, layer);
        }
      } else {
        auto itr = simhits_bh.find(c_id);
        if (itr != simhits_bh.end()){
          energy += itr->second;
        }
      }
    }
    return energy;
}

double HGCalWaferTableProducer::calibrate(double energy, int thickness, unsigned layer) {
  double fcPerMip = fcPerMip_[thickness];
  double thicknessCorrection = thicknessCorrections_[thickness];
  double layerWeight = layerWeights_[layer];
  double TeV2GeV = 1.e3;
  return energy * keV2fC_ / fcPerMip * layerWeight * TeV2GeV / thicknessCorrection;
}

void HGCalWaferTableProducer::simhits(const edm::Event& e,
                                                std::unordered_map<uint32_t, double>& simhits_ee,
                                                std::unordered_map<uint32_t, double>& simhits_fh,
                                                std::unordered_map<uint32_t, double>& simhits_bh) {
  edm::Handle<edm::PCaloHitContainer> ee_simhits_h;
  e.getByToken(simhits_ee_token_, ee_simhits_h);
  const edm::PCaloHitContainer& ee_simhits = *ee_simhits_h;
  edm::Handle<edm::PCaloHitContainer> fh_simhits_h;
  e.getByToken(simhits_fh_token_, fh_simhits_h);
  const edm::PCaloHitContainer& fh_simhits = *fh_simhits_h;
  edm::Handle<edm::PCaloHitContainer> bh_simhits_h;
  e.getByToken(simhits_bh_token_, bh_simhits_h);
  const edm::PCaloHitContainer& bh_simhits = *bh_simhits_h;

  // EE
  for (const auto& simhit : ee_simhits) {
    DetId id = simhit.id();
    if (id.rawId() == 0)
      continue;
    auto itr_insert = simhits_ee.emplace(id, 0.);
    itr_insert.first->second += simhit.energy();
  }
  //  FH
  for (const auto& simhit : fh_simhits) {
    DetId id = simhit.id();
    if (id.rawId() == 0)
      continue;
    auto itr_insert = simhits_fh.emplace(id, 0.);
    itr_insert.first->second += simhit.energy();
  }
  //  BH
  for (const auto& simhit : bh_simhits) {
    DetId id = simhit.id();
    if (id.rawId() == 0)
      continue;
    auto itr_insert = simhits_bh.emplace(id, 0.);
    itr_insert.first->second += simhit.energy();
  }
}

void HGCalWaferTableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // 1) Retrieve trigger cells
  edm::Handle<l1t::HGCalTriggerCellBxCollection> trigger_cells_h;
  iEvent.getByToken(trigger_cells_token_, trigger_cells_h);
  const l1t::HGCalTriggerCellBxCollection& trigger_cells = *trigger_cells_h;

  // 2) Sim hit association
  std::unordered_map<uint32_t, double> simhits_ee;
  std::unordered_map<uint32_t, double> simhits_fh;
  std::unordered_map<uint32_t, double> simhits_bh;
  std::unordered_map<uint32_t, bool> simhit_has_rechit;
  simhits(iEvent, simhits_ee, simhits_fh, simhits_bh);

  // 3) Access geometry
  triggerTools_.setGeometry(triggerGeometry_.product());
  aeInputUtil_.setGeometry(triggerGeometry_.product());

  // 4) Group trigger cells into wafers
  std::unordered_map<uint32_t, std::vector<l1t::HGCalTriggerCell>> wafers;
  for(const auto& tc : trigger_cells){
    if (triggerTools_.isScintillator(tc.detId())) continue;
    uint32_t wafer_id = triggerTools_.getTriggerGeometry()->getModuleFromTriggerCell(tc.detId());
    wafers[wafer_id].push_back(tc);
  }

  // 5) Compute wafer features
  for(const auto& wafer : wafers){
    ++wafer_n_;

    DetId id0(wafer.second[0].detId());

    // Wafer IDs, types, coordinates
    wafer_id_.emplace_back(wafer.first);
    wafer_side_.emplace_back(triggerTools_.zside(id0));
    wafer_layer_.emplace_back(triggerTools_.layerWithOffset(id0));
    if(id0.det() == DetId::HGCalTrigger){
      HGCalTriggerDetId idtrg0(id0);
      wafer_subdet_.emplace_back(idtrg0.subdet());
      wafer_waferu_.emplace_back(idtrg0.waferU());
      wafer_waferv_.emplace_back(idtrg0.waferV());
      wafer_wafertype_.emplace_back(idtrg0.type());
    } else if(id0.det() == DetId::HGCalHSc){
      HGCScintillatorDetId idsc0(id0);
      wafer_subdet_.emplace_back(idsc0.subdet());
      wafer_waferu_.emplace_back(-999);
      wafer_waferv_.emplace_back(-999);
      wafer_wafertype_.emplace_back(idsc0.type());
    } else {
      throw cms::Exception("HGCalTriggerNtupleWafers") << "Unknown detid " << id0.det() << std::endl;
    }

    // Wafer energies
    float mipPt=0, energy=0, simenergy=0;
    for(const auto& tc : wafer.second){
      mipPt += tc.mipPt();
      energy += tc.energy();
      simenergy += getsimenergy(tc.detId(), simhits_ee, simhits_fh, simhits_bh);
    }
    wafer_mipPt_.emplace_back(mipPt);
    wafer_energy_.emplace_back(energy);
    wafer_simenergy_.emplace_back(simenergy);

    // Wafer positions
    auto pos = triggerTools_.getTriggerGeometry()->getModulePosition(wafer.first);
    wafer_eta_.emplace_back(pos.eta());
    wafer_phi_.emplace_back(pos.phi());
    wafer_x_.emplace_back(pos.x());
    wafer_y_.emplace_back(pos.y());
    wafer_z_.emplace_back(pos.z());

    // Deposits within the wafers
    aeInputUtil_.run(wafer.second);
    wafer_sumCALQ_.emplace_back(aeInputUtil_.getModSum());
    for(unsigned i=0; i<nInputs_; ++i){
      wafer_ADC_[i].emplace_back(aeInputUtil_.getADC(i));
      wafer_norm_[i].emplace_back(aeInputUtil_.getNorm(i));
      wafer_CALQ_[i].emplace_back(aeInputUtil_.getCALQ(i));
      wafer_AEin_[i].emplace_back(aeInputUtil_.getInput(i)/aeInputUtil_.getInputNorm());
    }
  }

  // 6) Build FlatTable
  auto tab = std::make_unique<nanoaod::FlatTable>(wafer_id_.size(), "L1THGCAL_wafer", false, false);
  tab->addColumn<int>("id", wafer_id_, "");
  tab->addColumn<int>("zside", wafer_side_, "");
  tab->addColumn<int>("layer", wafer_layer_, "");
  tab->addColumn<int>("subdet", wafer_subdet_, "");
  tab->addColumn<int>("waferu", wafer_waferu_, "");
  tab->addColumn<int>("waferv", wafer_waferv_, "");
  tab->addColumn<int>("wafertype", wafer_wafertype_, "");
  tab->addColumn<float>("mipPt", wafer_mipPt_, "");
  tab->addColumn<float>("energy", wafer_energy_, "");
  tab->addColumn<float>("simenergy", wafer_simenergy_, "");
  tab->addColumn<float>("eta", wafer_eta_, "");
  tab->addColumn<float>("phi", wafer_phi_, "");
  tab->addColumn<float>("x", wafer_x_, "");
  tab->addColumn<float>("y", wafer_y_, "");
  tab->addColumn<float>("z", wafer_z_, "");
  tab->addColumn<uint64_t>("sumCALQ", wafer_sumCALQ_, "");
  for (unsigned i = 0; i < nInputs_; ++i) {
    tab->addColumn<float>(Form("ADC_%u", i), wafer_ADC_[i], "");
    tab->addColumn<float>(Form("norm_%u", i), wafer_norm_[i], "");
    tab->addColumn<float>(Form("CALQ_%u", i), wafer_CALQ_[i], "");
    tab->addColumn<float>(Form("AEin_%u", i), wafer_AEin_[i], "");
  }

  iEvent.put(std::move(tab), "wafers");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HGCalWaferTableProducer);

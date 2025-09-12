import FWCore.ParameterSet.Config as cms

from L1Trigger.L1THGCalUtilities.hgcalTriggerNtuples_cfi import *

L1THGCalTriggerNtuples = cms.Sequence(l1tHGCalTriggerNtuplizer)

from PhysicsTools.NanoAOD.common_cff import Var, ExtVar
def LazyVar(expr, valtype, doc=None, precision=-1):
    return Var(expr, valtype, doc, precision, lazyEval=True)

hgcClusterTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag('l1tHGCalBackEndLayer2Producer:HGCalBackendLayer2Processor3DClustering'),
    doc = cms.string("HGCal 3D clusters"),
    cut  = cms.string("pt > 1"),
    name = cms.string("HGCal3DCl"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table
    variables = cms.PSet(
        pt  = LazyVar("pt",  float,precision=8),
        phi = LazyVar("phi", float,precision=8),
        eta  = LazyVar("eta", float,precision=8),
        nTcs = LazyVar("constituents.size",  int,precision=8),
        ptEm = LazyVar("iPt('EM')",  float,precision=8),
        hwQual = LazyVar("hwQual",  int,precision=8),
        showerlength = LazyVar("showerLength", int),
        coreshowerlength = LazyVar("coreShowerLength", int),
        firstlayer = LazyVar("firstLayer", int),
        maxlayer = LazyVar("maxLayer", int),
        seetot = LazyVar("sigmaEtaEtaTot", float),
        seemax = LazyVar("sigmaEtaEtaMax", float),
        spptot = LazyVar("sigmaPhiPhiTot", float),
        sppmax = LazyVar("sigmaPhiPhiMax", float),
        szz = LazyVar("sigmaZZ", float),
        srrtot = LazyVar("sigmaRRTot", float),
        srrmax = LazyVar("sigmaRRMax", float),
        srrmean = LazyVar("sigmaRRMean", float),
        emaxe = LazyVar("eMax/energy", float),
        hoe = LazyVar("hOverE", float),
        meanz = LazyVar("abs(zBarycenter)", float),
        layer10 = LazyVar("layer10percent", float),
        layer50 = LazyVar("layer50percent", float),
        layer90 = LazyVar("layer90percent", float),
        ntc67 = LazyVar("triggerCells67percent", float),
        ntc90 = LazyVar("triggerCells90percent", float),
        varRR = LazyVar("varRR", float),
        varZZ = LazyVar("varZZ", float),
        varEtaEta = LazyVar("varEtaEta", float),
        varPhiPhi = LazyVar("varPhiPhi", float),
        first1layers = LazyVar("first1layers", float),
        first3layers = LazyVar("first3layers", float),
        first5layers = LazyVar("first5layers", float),
        firstHcal1layers = LazyVar("firstHcal1layers", float),
        firstHcal3layers = LazyVar("firstHcal3layers", float),
        firstHcal5layers = LazyVar("firstHcal5layers", float),
        last1layers = LazyVar("last1layers", float),
        last3layers = LazyVar("last3layers", float),
        last5layers = LazyVar("last5layers", float),
        emax1layers = LazyVar("emax1layers", float),
        emax3layers = LazyVar("emax3layers", float),
        emax5layers = LazyVar("emax5layers", float),
        eot = LazyVar("eot", float),
        ebm0 = LazyVar("ebm0", int),
        ebm1 = LazyVar("ebm1", int),
        hbm = LazyVar("hbm", int),
        )
)

from L1Trigger.Phase2L1ParticleFlow.l1tPFClustersFromHGC3DClusters_cfi import l1tPFClustersFromHGC3DClusters
from L1Trigger.L1THGCal.egammaIdentification import egamma_identification_histomax
hgcClusterExtTable = cms.EDProducer("L1HGC3DclTableProducer",
                src = cms.InputTag('l1tHGCalBackEndLayer2Producer:HGCalBackendLayer2Processor3DClustering'),
                cut  = cms.string("pt > 1"),
                name = cms.string("HGCal3DCl"),
                emVsPionID=l1tPFClustersFromHGC3DClusters.emVsPionID,
                emVsPUID=l1tPFClustersFromHGC3DClusters.emVsPUID,
                EGIdentification = egamma_identification_histomax.clone(),
            )

## HGCAL objects
p2HGCTablesTask = cms.Task(
    hgcClusterTable,
    hgcClusterExtTable,
)

hgcPh2NanoSequence = cms.Sequence(p2HGCTablesTask)
# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: --python_filename private-RunIIAutumn18NanoAODv7_NANO.py --no_exec -n -1 --step NANO --eventcontent NANOAODSIM --conditions 102X_upgrade2018_realistic_v21 --era Run2_2018,run2_nanoAOD_102Xv1 --mc --customise_commands=process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)));process.MessageLogger.cerr.FwkReport.reportEvery=100 --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAODSIM --filein file:input.root --fileout file:output.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('NANO',eras.Run2_2018,eras.run2_nanoAOD_102Xv1)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:input.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('--python_filename nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:output.root'),
    outputCommands = process.NANOAODSIMEventContent.outputCommands
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v21', '')

# Path and EndPath definitions
process.nanoAOD_step = cms.Path(process.nanoSequenceMC)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step,process.NANOAODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeMC 

#call to customisation function nanoAOD_customizeMC imported from PhysicsTools.NanoAOD.nano_cff
process = nanoAOD_customizeMC(process)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions

# Customisation from command line

process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)));process.MessageLogger.cerr.FwkReport.reportEvery=100
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

## private customization for SOS Displaced production
process.finalMuons.cut = "isPFMuon"
#
SOSDisp_Muon_extras = cms.PSet(
    nPixelLayers = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('number of pixel hits'),
        expr = cms.string('?track.isNonnull?innerTrack().hitPattern().pixelLayersWithMeasurement():0'),
        mcOnly = cms.bool(False),
        precision = cms.int32(-1),
        type = cms.string('int')
    ),
    TMOneStationTight = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('trk track matched with at least one muon station'),
        expr = cms.string("isGood(\'TMOneStationTight\')"),
        mcOnly = cms.bool(False),
        precision = cms.int32(-1),
        type = cms.string('bool')
    )
)

SOSDisp_GenPart_extras = cms.PSet(
    vtx_x = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('GenPart Vertex X-Coord'),
        expr = cms.string('vx'),
        mcOnly = cms.bool(True),
        precision = cms.int32(-1),
        type = cms.string('float')
    ),
    vtx_y = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('GenPart Vertex Y-Coord'),
        expr = cms.string('vy'),
        mcOnly = cms.bool(True),
        precision = cms.int32(-1),
        type = cms.string('float')
    ),
    vtx_z = cms.PSet(
        compression = cms.string('none'),
        doc = cms.string('GenPart Vertex Z-Coord'),
        expr = cms.string('vz'),
        mcOnly = cms.bool(True),
        precision = cms.int32(-1),
        type = cms.string('float')
    )
)

process.muonTable.variables = cms.PSet(process.muonTable.variables, SOSDisp_Muon_extras)
process.genParticleTable.variables = cms.PSet(process.genParticleTable.variables, SOSDisp_GenPart_extras)

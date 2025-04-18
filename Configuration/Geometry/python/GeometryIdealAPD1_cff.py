import FWCore.ParameterSet.Config as cms

#
# Geometry master configuration
#
# Ideal geometry, needed for simulation
from Geometry.CMSCommonData.cmsIdealGeometryAPD1XML_cfi import *
from Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cff import *
# Reconstruction geometry services
#  Tracking Geometry
from Geometry.CommonTopologies.globalTrackingGeometry_cfi import *

#Tracker
from RecoTracker.GeometryESProducer.TrackerRecoGeometryESProducer_cfi import *
from Geometry.TrackerGeometryBuilder.TrackerAdditionalParametersPerDet_cfi import *

#Muon
from Geometry.MuonNumbering.muonGeometryConstants_cff import *
from Geometry.MuonNumbering.muonNumberingInitialization_cfi import *
from Geometry.MuonNumbering.muonOffsetESProducer_cff import *
from RecoMuon.DetLayers.muonDetLayerGeometry_cfi import *

#  Alignment
from Geometry.TrackerGeometryBuilder.idealForDigiTrackerGeometry_cff import *
from Geometry.CSCGeometryBuilder.idealForDigiCscGeometry_cff import *
from Geometry.DTGeometryBuilder.idealForDigiDtGeometry_cff import *

#  Calorimeters
from Geometry.CaloEventSetup.CaloTopology_cfi import *
from Geometry.CaloEventSetup.CaloGeometry_cff import *
from Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi import *
from Geometry.EcalMapping.EcalMapping_cfi import *
from Geometry.EcalMapping.EcalMappingRecord_cfi import *
from Geometry.EcalCommonData.ecalSimulationParameters_cff import *
from Geometry.HcalCommonData.hcalDDConstants_cff import *
from Geometry.ForwardGeometry.zdcTopologyEP_cfi import *

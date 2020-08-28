import FWCore.ParameterSet.Config as cms

import PhysicsTools.PatAlgos.tools.helpers as configtools
from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet
from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
#0;136;0cfrom PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC, miniAOD_customizeAllData
# lower the pt threshold
def lowerTauPt(process,postfix='',tauPt=8, jetPt=5):
    from FWCore.ParameterSet.MassReplace import massSearchReplaceParam
    massSearchReplaceParam(getattr(process,'PATTauSequence'+postfix),'minJetPt',14,jetPt)
    getattr(process,'selectedPatTaus'+postfix).cut = cms.string("pt > {} && tauID(\'decayModeFindingNewDMs\')> 0.5".format(tauPt))

def addCustomization(process,options,**kwargs):
    doMM = kwargs.pop('doMM',False)
    doMT = kwargs.pop('doMT',False)
    #doSlimming = kwargs.pop('doSlimming',False)

    process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
    patAlgosToolsTask = configtools.getPatAlgosToolsTask(process)

    process.PATTauSequence = cms.Sequence(process.PFTau+process.makePatTaus+process.selectedPatTaus)
    #process.recoTauAK4PFJets08Region.pfCandAssocMapSrc = cms.InputTag("")


    #########################
    ### Muon Cleaned Taus ###
    #########################
    process.recoMuonsForJetCleaning = cms.EDFilter('MuonRefSelector',
        src = cms.InputTag('muons'),
        cut = cms.string('pt > 3.0 && isPFMuon && (isGlobalMuon || isTrackerMuon)'),
    )

    process.ak4PFJetsMuonCleaned = cms.EDProducer(
        'MuonCleanedJetProducer',
        jetSrc = cms.InputTag("ak4PFJets"),
        muonSrc = cms.InputTag("recoMuonsForJetCleaning"),
        pfCandSrc = cms.InputTag("particleFlow"),
        pfCandCollection=cms.InputTag("particleFlow"),     
    )
    
    process.muonCleanedHPSPFTausTask = cms.Task(
        process.recoMuonsForJetCleaning,
        process.ak4PFJetsMuonCleaned
    )
    patAlgosToolsTask.add(process.muonCleanedHPSPFTausTask)

    jetSrc = 'ak4PFJetsMuonCleaned'


    pfAssocMap = cms.InputTag('ak4PFJetsMuonCleaned','pfCandAssocMapForIsolation')
    process.PATTauSequenceMuonCleaned = cloneProcessingSnippet(process,process.PATTauSequence, 'MuonCleaned', addToTask = True)
    massSearchReplaceAnyInputTag(process.PATTauSequenceMuonCleaned,cms.InputTag('ak4PFJets'),cms.InputTag(jetSrc))

    process.recoTauAK4PFJets08RegionMuonCleaned.pfCandAssocMapSrc = pfAssocMap
    process.slimmedTausMuonCleaned = process.slimmedTausNoDeepIDs.clone(src = cms.InputTag('selectedPatTausMuonCleaned'))
    patAlgosToolsTask.add(process.slimmedTausMuonCleaned)
    lowerTauPt(process,'MuonCleaned')
    
    if options.isMC:
        process.tauGenJetsMuonCleaned.GenParticles = "prunedGenParticles"
        process.patTausMuonCleaned.embedGenMatch = False
    else:
        from PhysicsTools.PatAlgos.tools.coreTools import _removeMCMatchingForPATObject
        attrsToDelete = []
        postfix = ''
        print "removing MC dependencies for tausMuonCleaned"
        _removeMCMatchingForPATObject(process, 'tauMatch', 'patTausMuonCleaned',postfix)
        ## remove mc extra configs for taus
        tauProducer = getattr(process,'patTausMuonCleaned'+postfix)
        tauProducer.addGenJetMatch   = False
        tauProducer.embedGenJetMatch = False
        attrsToDelete += [tauProducer.genJetMatch.getModuleLabel()]
        tauProducer.genJetMatch      = ''
        attrsToDelete += ['tauGenJetsMuonCleaned'+postfix]
        attrsToDelete += ['tauGenJetsSelectorAllHadronsMuonCleaned'+postfix]
        for attr in attrsToDelete:
            if hasattr(process,attr): delattr(process,attr)

    process.MuonCollectionClean = cms.EDProducer("MuonCollectionClean",
        MuTag   = cms.InputTag("slimmedMuons"),
        TauTag  = cms.InputTag("slimmedTausMuonCleaned"),
        ParticleFlowCandTag = cms.InputTag("packedPFCandidates"),
        muonID  = cms.string('loose'),
        ptCut   = cms.double(3.0),
    )

    patAlgosToolsTask.add(process.MuonCollectionClean)

    # add tauID
    updatedTauName = "slimmedTausNewIDMuonCleaned"
    #import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
    import BoostedDiTau.Skimmer.runTauIdMVA as tauIdConfig
    tauIdEmbedderMuonCleaned = tauIdConfig.TauIDEmbedder(process, cms, debug = True,
                                                         updatedTauName = updatedTauName,
                                                         tauSrc = 'slimmedTausMuonCleaned',
                                                         toKeep = ["2017v2",
                                                                   "deepTau2017v2p1", #deepTau TauIDs
                                                         ])
    tauIdEmbedderMuonCleaned.runTauID()
    
    process.tauIdMuonCleaned = cms.Task(
        process.rerunMvaIsolationTaskMuonCleaned,
        process.slimmedTausNewIDMuonCleaned,
    )

    patAlgosToolsTask.add(process.tauIdMuonCleaned)

    #############################
    ### Electron cleaned taus ###
    #############################
    process.recoElectronsForJetCleaning = cms.EDFilter('ElectronFilter',
                                               vertex = cms.InputTag("offlinePrimaryVerticesWithBS"),
                                               Rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                               electrons = cms.InputTag("gedGsfElectrons"),
                                               conv = cms.InputTag("conversions"),
                                               BM = cms.InputTag("offlineBeamSpot"),
                                               Tracks = cms.InputTag("electronGsfTracks"),
                                               #Passcount =cms.uint32(1),
                                               )
    
    process.ak4PFJetsElectronCleaned = cms.EDProducer(
        'ElectronCleanedJetProducer',
        jetSrc = cms.InputTag("ak4PFJets"),
        electronSrc = cms.InputTag("recoElectronsForJetCleaning","LooseElectronRef"),
        pfCandSrc = cms.InputTag("particleFlow"),
        pfCandCollection=cms.InputTag("particleFlow"),
    )

##     process.ak4PFJetsElectronCleanedFlat = cms.EDProducer(
##         'ElectronCleanedJetProducer',
##         jetSrc = cms.InputTag("ak4PFJets"),
##         electronSrc = cms.InputTag("recoElectronsForJetCleaning","LooseElectronRefFlat"),
##         pfCandSrc = cms.InputTag("particleFlow"),
##         pfCandCollection=cms.InputTag("particleFlow"),
##         )

    process.electronCleanedHPSPFTausTask = cms.Task(
        process.recoElectronsForJetCleaning,
        process.ak4PFJetsElectronCleaned
        #process.ak4PFJetsElectronCleanedFlat
    )
    
    patAlgosToolsTask.add(process.electronCleanedHPSPFTausTask)
    
    jetSrc = 'ak4PFJetsElectronCleaned'
    
    pfAssocMaps = cms.InputTag('ak4PFJetsElectronCleaned','pfCandAssocMapForIsolation')
    process.PATTauSequenceElectronCleaned = cloneProcessingSnippet(process,process.PATTauSequence, 'ElectronCleaned', addToTask = True)
    massSearchReplaceAnyInputTag(process.PATTauSequenceElectronCleaned,cms.InputTag('ak4PFJets'),cms.InputTag(jetSrc))
    
    process.recoTauAK4PFJets08RegionElectronCleaned.pfCandAssocMapSrc = pfAssocMaps
    process.slimmedTausElectronCleaned = process.slimmedTausNoDeepIDs.clone(src = cms.InputTag('selectedPatTausElectronCleaned'))
    patAlgosToolsTask.add(process.slimmedTausElectronCleaned)

    lowerTauPt(process,'ElectronCleaned')

    process.ElectronCollectionClean = cms.EDProducer("ElectronCollectionClean",
        EleTag      = cms.InputTag("slimmedElectrons"),
        TauTag      = cms.InputTag("slimmedTausElectronCleaned"),
        ParticleFlowCandTag = cms.InputTag("packedPFCandidates"),
        rhoTag      = cms.InputTag("fixedGridRhoAll"),
        electronID  = cms.string('loose'),
        ptCut       = cms.double(7.0),
        passRelIso  = cms.bool(False),
        effAreasConfigFile = cms.FileInPath("BoostedDiTau/Skimmer/data/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
    )

    patAlgosToolsTask.add(process.ElectronCollectionClean)

    # add tauID
    updatedTauName = "slimmedTausNewIDElectronCleaned"
    #import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
    import BoostedDiTau.Skimmer.runTauIdMVA as tauIdConfig
    tauIdEmbedderElectronCleaned = tauIdConfig.TauIDEmbedder(process, cms, debug = True,
                                                         updatedTauName = updatedTauName,
                                                         tauSrc = 'slimmedTausElectronCleaned',
                                                         toKeep = ["2017v2",
                                                                   "deepTau2017v2p1", #deepTau TauIDs
                                                         ])
    tauIdEmbedderElectronCleaned.runTauID()
    
    process.tauIdElectronCleaned = cms.Task(
        process.rerunMvaIsolationTaskElectronCleaned,
        process.slimmedTausNewIDElectronCleaned,
    )

    patAlgosToolsTask.add(process.tauIdElectronCleaned)

    

##     jetSrcF = 'ak4PFJetsElectronCleanedFlat'
##     
##     pfAssocMapsFlat = cms.InputTag('ak4PFJetsElectronCleanedFlat','pfCandAssocMapForIsolation')
##     process.PATTauSequenceElectronCleanedFlat = cloneProcessingSnippet(process,process.PATTauSequence, 'ElectronCleanedFlat', addToTask = True)
##     massSearchReplaceAnyInputTag(process.PATTauSequenceElectronCleanedFlat,cms.InputTag('ak4PFJets'),cms.InputTag(jetSrcF))
##     
##     process.recoTauAK4PFJets08RegionElectronCleanedFlat.pfCandAssocMapSrc = pfAssocMapsFlat
##     process.slimmedTausElectronCleanedFlat = process.slimmedTaus.clone(src = cms.InputTag('selectedPatTausElectronCleanedFlat'))
##     patAlgosToolsTask.add(process.slimmedTausElectronCleanedFlat)
##     
##     lowerTauPt(process,'ElectronCleanedFlat')
    
    if options.isMC:
        process.tauGenJetsElectronCleaned.GenParticles = "prunedGenParticles"
        process.patTausElectronCleaned.embedGenMatch = False
    else:
        from PhysicsTools.PatAlgos.tools.coreTools import _removeMCMatchingForPATObject
        attrsToDelete = []
        postfix = ''
        print "removing MC dependencies for tausElectronCleaned"
        _removeMCMatchingForPATObject(process, 'tauMatch', 'patTausElectronCleaned',postfix)
        ## remove mc extra configs for taus
        tauProducer = getattr(process,'patTausElectronCleaned'+postfix)
        tauProducer.addGenJetMatch   = False
        tauProducer.embedGenJetMatch = False
        attrsToDelete += [tauProducer.genJetMatch.getModuleLabel()]
        tauProducer.genJetMatch      = ''
        attrsToDelete += ['tauGenJetsElectronCleaned'+postfix]
        attrsToDelete += ['tauGenJetsSelectorAllHadronsElectronCleaned'+postfix]
        for attr in attrsToDelete:
            if hasattr(process,attr): delattr(process,attr)

    
    #############################
    ### lower pt for nonclean ###
    #############################
    #process.combinatoricRecoTaus.minJetPt = cms.double(8.0)
    #process.recoTauAK4PFJets08Region.minJetPt = cms.double(8.0)
    #process.ak4PFJetsLegacyHPSPiZeros.minJetPt = cms.double(8.0)
    #process.ak4PFJetsRecoTauChargedHadrons.minJetPt = cms.double(8.0)
    #process.selectedPatTaus.cut = cms.string("pt > 8. && tauID(\'decayModeFindingNewDMs\')> 0.5")
    lowerTauPt(process)
    
    ############################
    ### lower pt for boosted ###
    ############################
    process.ca8PFJetsCHSprunedForBoostedTaus.jetPtMin = cms.double(20.0)
    process.ca8PFJetsCHSprunedForBoostedTaus.subjetPtMin = cms.double(8.0)
    #process.combinatoricRecoTausBoosted.minJetPt = cms.double(8.0)
    #process.recoTauAK4PFJets08RegionBoosted.minJetPt = cms.double(8.0)
    #process.selectedPatTausBoosted.cut = cms.string("pt > 8. && tauID(\'decayModeFindingNewDMs\')> 0.5")
    lowerTauPt(process,'Boosted')

    #######################################
    ### Update btagging for cleaned jet ###
    #######################################
    
    from RecoJets.JetAssociationProducers.j2tParametersVX_cfi import j2tParametersVX
    process.ak4PFJetsMuonCleanedTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
        j2tParametersVX,
        jets = cms.InputTag("ak4PFJetsMuonCleaned")
    )
    process.patJetMuonCleanedCharge = cms.EDProducer("JetChargeProducer",
        src = cms.InputTag("ak4PFJetsMuonCleanedTracksAssociatorAtVertex"),
        var = cms.string('Pt'),
        exp = cms.double(1.0)
    )
    
    from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
    addJetCollection(process, postfix   = "", labelName = 'MuonCleaned', jetSource = cms.InputTag('ak4PFJetsMuonCleaned'),
                    jetCorrections = ('AK4PF', ['L2Relative', 'L3Absolute'], ''),
                    algo= 'AK', rParam = 0.4, btagDiscriminators = map(lambda x: x.value() ,process.patJets.discriminatorSources)
    )
    
    if options.isMC: process.patJetGenJetMatchMuonCleaned.matched = 'slimmedGenJets'
    process.patJetsMuonCleaned.jetChargeSource = cms.InputTag("patJetMuonCleanedCharge")
    
    process.slimmedJetsMuonCleaned = process.slimmedJets.clone(src = cms.InputTag("selectedPatJetsMuonCleaned"))


    #################
    ### Skim Path ###
    #################
    #process.main_path = cms.Path()
    #process.main_path_mm = cms.Path()
    #process.main_path_em = cms.Path()
    #process.main_path_mt = cms.Path()
    #process.main_path_et = cms.Path()
    #process.main_path_tt = cms.Path()
    
    #process.main_path_mt *= process.analysisMuonsNoIsoMTCount
    #process.main_path_mt *= process.analysisTausMuonCleaned
    #process.main_path_mt *= process.analysisTausMuonCleanedCount
    #process.main_path_ut *= process.analysisTausUnCleaned
    #process.main_path_ut *= process.analysisTausUnCleanedCount
    #process.main_path_et *= process.analysisTausElectronCleaned
    #process.main_path_et *= process.analysisTausElectronCleanedCount
    #process.z_tau_eff_path *= process.analysisTaus
    #process.z_tau_eff_path *= process.analysisTausCount
    #process.z_tau_eff_muclean_path *= process.analysisTausMuonCleaned
    #process.z_tau_eff_muclean_path *= process.analysisTausMuonCleanedCount

    # and for mt require third muon
    
    ############################
    ### Tau Eff requirements ###
    ############################
##     process.mumuZTauEff = cms.EDProducer("CandViewShallowCloneCombiner",
##         decay = cms.string("{0} {1}".format('slimmedMuons','analysisTaus')),
##         checkCharge = cms.bool(False),
##         cut   = cms.string("30<mass<210 && deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi)>0.5"),
##     )
##     process.mumuZCountTauEff = cms.EDFilter("PATCandViewCountFilter",
##          minNumber = cms.uint32(1),
##          maxNumber = cms.uint32(999),
##          src = cms.InputTag('mumuZTauEff'),
##     )
##     process.mumuZMuonCleanedTauEff = cms.EDProducer("CandViewShallowCloneCombiner",
##         decay = cms.string("{0} {1}".format('slimmedMuons','analysisTausMuonCleaned')),
##         checkCharge = cms.bool(False),
##         cut   = cms.string("30<mass<210 && deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi)>0.5"),
##     )
##     process.mumuZMuonCleanedCountTauEff = cms.EDFilter("PATCandViewCountFilter",
##          minNumber = cms.uint32(1),
##          maxNumber = cms.uint32(999),
##          src = cms.InputTag('mumuZMuonCleanedTauEff'),
##     )
##     process.z_tau_eff_path *= process.mumuZTauEff
##     process.z_tau_eff_path *= process.mumuZCountTauEff
##     process.z_tau_eff_muclean_path *= process.mumuZMuonCleanedTauEff
##     process.z_tau_eff_muclean_path *= process.mumuZMuonCleanedCountTauEff
    
    #################
    ### Finish up ###
    #################
    # add to schedule
    #process.schedule.append(process.main_path)
    #process.schedule.append(process.main_path_em)
    #process.schedule.append(process.main_path_et)
    #process.schedule.append(process.main_path_mt)
    #process.schedule.append(process.main_path_ut)
    #process.schedule.append(process.main_path_tt)
    #process.schedule.append(process.z_path)
    #process.schedule.append(process.z_tau_eff_path)
    #process.schedule.append(process.z_tau_eff_muclean_path)

    # lumi summary
    if not options.isMC:
        process.TFileService = cms.Service("TFileService",
                                           fileName = cms.string(options.outputFile.split('.root')[0]+'_lumi.root'),
    )
    
        process.lumiTree = cms.EDAnalyzer("LumiTree",
                                          genEventInfo = cms.InputTag("generator"),
                                          lheEventProduct = cms.InputTag("externalLHEProducer"),
        )
        process.lumi_step = cms.Path(process.lumiTree)
        process.schedule.append(process.lumi_step)
    
        process.lumiSummary = cms.EDProducer("LumiSummaryProducer",
                                             genEventInfo = cms.InputTag("generator"),
                                             lheEventProduct = cms.InputTag("externalLHEProducer"),
        )
        process.lumiSummary_step = cms.Path(process.lumiSummary)
        process.schedule.append(process.lumiSummary_step)
    
    if not options.isMC:
        process.MINIAODSIMEventContent.outputCommands += [
            'drop *_ctppsLocalTrackLiteProducer_*_*', # Don't know what this is, but it prevents running in older releases
        ]
    
    # additional skims
    if doMM:
        process.MINIAODSIMoutputZSKIM = process.MINIAODSIMoutput.clone(
            SelectEvents = cms.untracked.PSet(
                SelectEvents = cms.vstring('z_path'),
            ),
            fileName = cms.untracked.string(options.outputFile.split('.root')[0]+'_zskim.root'),
        )
        process.MINIAODSIMoutputZSKIM_step = cms.EndPath(process.MINIAODSIMoutputZSKIM)
        process.schedule.append(process.MINIAODSIMoutputZSKIM_step)
    
    if doMT:
        process.MINIAODSIMoutputZMUTAUSKIM = process.MINIAODSIMoutput.clone(
            SelectEvents = cms.untracked.PSet(
                SelectEvents = cms.vstring('z_tau_eff_path','z_tau_eff_muclean_path'),
            ),
            fileName = cms.untracked.string(options.outputFile.split('.root')[0]+'_zmutauskim.root'),
        )
        process.MINIAODSIMoutputZMUTAUSKIM_step = cms.EndPath(process.MINIAODSIMoutputZMUTAUSKIM)
        process.schedule.append(process.MINIAODSIMoutputZMUTAUSKIM_step)


    if not options.isStandard:
        process.MINIAODSIMEventContent.outputCommands += [
            'keep *_slimmedTausMuonCleaned_*_*',
            'keep *_slimmedTausElectronCleaned_*_*',
            'keep *_ak4PFJetsElectronCleaned_*_*',
            'keep *_ak4PFJetsMuonCleaned_*_*',
            #'keep *_ak4PFJetsElectronCleanedFlat_*_*',
            #'keep *_ak4PFJets_*_*',
        ]

    #process.MINIAODoutput.SelectEvents = cms.untracked.PSet(
    #    SelectEvents = cms.vstring('main_path'),
    #    #SelectEvents = cms.vstring('main_path_em','main_path_et','main_path_mt','main_path_tt'),
    #)

    if options.isReMINIAOD:
        process.MINIAODSIMEventContent.outputCommands += [
            'keep *_slimmedTausNewIDMuonCleaned_*_*',
            'keep *_slimmedTausNewIDElectronCleaned_*_*',
        ]
    
    if options.doSlimming:
    
##        ###############
##        ### Trigger ###
##        ###############
##        process.HLT =cms.EDFilter("HLTHighLevel",
##             TriggerResultsTag = cms.InputTag("TriggerResults","","RECO"),
##             #HLTPaths = cms.vstring("HLT_IsoMu27_v*", "HLT_IsoTkMu27_v*"), # 2017
##             HLTPaths = cms.vstring("HLT_IsoMu24_v*", "HLT_IsoTkMu24_v*", "HLT_IsoMu27_v*", "HLT_IsoTkMu27_v*"),
##             eventSetupPathsKey = cms.string(''),
##             andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
##             throw = cms.bool(False) # throw exception on unknown path names
##        )
##        #process.main_path *= process.HLT
        
    
        ###############
        ### Muon ID ###
        ###############
        process.selectedPATMuons = cms.EDFilter('PATMuonBaseLineSelection',
            muons = cms.InputTag("slimmedMuons"),
            vertex = cms.InputTag("offlineSlimmedPrimaryVertices")
        )
        
    
        ###################
        ### Electron ID ###
        ###################
        process.selectedPATElectrons = cms.EDFilter('PATElectronBaseLineSelection',
            electrons = cms.InputTag("slimmedElectrons"),
            vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
            Rho = cms.InputTag("fixedGridRhoFastjetAll"),
            conv = cms.InputTag("conversions"),
            BM = cms.InputTag("offlineBeamSpot"),
            isDebug = cms.untracked.bool(False)
        )
        
        ##############
        ### Jet ID ###
        ##############
        process.selectedPATJetsNoID = cms.EDFilter('PATJetSelector',
            src = cms.InputTag('slimmedJets'),
            cut = cms.string('pt > 20.0 && abs(eta)<2.5')
        )
    
        from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
        process.selectedPATJets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                                filterParams = cms.PSet(
                                                    version = cms.string('WINTER17'),
                                                    quality = cms.string('TIGHTLEPVETO')
                                                ),
                                                src = cms.InputTag("selectedPATJetsNoID"),
                                                filter = cms.bool(True)
        )
    
        process.PATJetBaseLineSelectionTask = cms.Task(
            process.selectedPATJetsNoID,
            process.selectedPATJets
        )
        
##        #########################
##        ### Trigger Threshold ###
##        #########################
##        #process.triggerMuon = cms.EDFilter('PATMuonSelector',
##        process.triggerMuon = cms.EDFilter('MuonSelector',
##            #src = cms.InputTag('secondMuon'),
##            src = cms.InputTag('analysisMuonsNoIso'),
##            #cut = cms.string('pt > 27.0'),
##            cut = cms.string('pt > 24.0'),
##        )
##        process.triggerMuonCount = cms.EDFilter("PATCandViewCountFilter",
##             minNumber = cms.uint32(1),
##             maxNumber = cms.uint32(999),
##             src = cms.InputTag('triggerMuon'),
##        )
        
        ########################
        ### Tau requirements ###
        ########################
        process.selectedPATTaus = cms.EDFilter('PATTauSelector',
            src = cms.InputTag('NewTauIDsEmbedded'),
#            src = cms.InputTag('slimmedTaus'),
            cut = cms.string('pt > 10.0 && abs(eta)<2.3 && tauID(\'decayModeFinding\')> 0.5'),
#            cut = cms.string('pt > 10.0 && abs(eta)<2.3 && tauID(\'decayModeFinding\')> 0.5 && tauID(\'byIsolationMVArun2v1DBoldDMwLTraw\')>-0.5'),
        )
    
        process.selectedPATTausMuonCleaned = cms.EDFilter('PATTauSelector',
            src = cms.InputTag('NewTauIDsEmbeddedMuonCleaned'),
#            src = cms.InputTag('slimmedTausMuonCleaned'),
            cut = cms.string('pt > 10.0 && abs(eta)<2.3 && tauID(\'decayModeFinding\')> 0.5'),
#            cut = cms.string('pt > 10.0 && abs(eta)<2.3 && tauID(\'decayModeFinding\')> 0.5 && tauID(\'byIsolationMVArun2v1DBoldDMwLTraw\')>-0.5'),
        )
    
        process.selectedPATTausElectronCleaned = cms.EDFilter('PATTauSelector',
            src = cms.InputTag('NewTauIDsEmbeddedElectronCleaned'),
#            src = cms.InputTag('slimmedTausElectronCleaned'),
            cut = cms.string('pt > 10.0 && abs(eta)<2.3 && tauID(\'decayModeFinding\')> 0.5'),
#            cut = cms.string('pt > 10.0 && abs(eta)<2.3 && tauID(\'decayModeFinding\')> 0.5 && tauID(\'byIsolationMVArun2v1DBoldDMwLTraw\')>-0.5'),
        )

##         process.selectedPATTausElectronCleanedFlat = cms.EDFilter('PATTauSelector',
##             src = cms.InputTag('NewTauIDsEmbeddedElectronCleanedFlat'),
#            src = cms.InputTag('slimmedTausElectronCleaned'),
##             cut = cms.string('pt > 10.0 && abs(eta)<2.3 && tauID(\'decayModeFinding\')> 0.5'),
#            cut = cms.string('pt > 10.0 && abs(eta)<2.3 && tauID(\'decayModeFinding\')> 0.5 && tauID(\'byIsolationMVArun2v1DBoldDMwLTraw\')>-0.5'),
##         )
    
        process.selectedPATTausBoosted = cms.EDFilter('PATTauSelector',
            src = cms.InputTag('slimmedTausElectronCleaned'),
            cut = cms.string('pt > 10.0 && abs(eta)<2.3 && tauID(\'decayModeFinding\')> 0.5'),
#            cut = cms.string('pt > 10.0 && abs(eta)<2.3 && tauID(\'decayModeFinding\')> 0.5 && tauID(\'byIsolationMVArun2v1DBoldDMwLTraw\')>-0.5'),
        )
    
        patAlgosToolsTask.add(process.selectedPATMuons)
        patAlgosToolsTask.add(process.selectedPATElectrons)
        patAlgosToolsTask.add(process.PATJetBaseLineSelectionTask)
        patAlgosToolsTask.add(process.selectedPATTaus)
        patAlgosToolsTask.add(process.selectedPATTausMuonCleaned)
        patAlgosToolsTask.add(process.selectedPATTausElectronCleaned)
        patAlgosToolsTask.add(process.selectedPATTausElectronCleanedFlat)
        patAlgosToolsTask.add(process.selectedPATTausBoosted)
    
        # additional changes to standard MiniAOD content
        process.MINIAODSIMEventContent.outputCommands = cms.untracked.vstring(
            'drop *',
            #'keep *_slimmedTausMuonCleaned_*_*',
            #'keep *_slimmedTausElectronCleaned_*_*',
            'keep *GsfElectron*_reducedEgamma_*_*',
            'keep *SuperCluster*_reducedEgamma_*_*',
            'keep *_fixedGridRhoFastjetAll_*_*',
            'keep *_TriggerResults_*_*',
            'keep recoGenParticles_prunedGenParticles_*_*',
            'keep *_generator_*_*',
            'keep *_slimmedGenJets_*_*',
            'keep *_selectedPATMuons_*_*',
            'keep *_selectedPATElectrons_*_*',
            'keep *_selectedPATJets_*_*',
            'keep *_selectedPATTaus_*_*',
            'keep *_selectedPATTausMuonCleaned_*_*',
            'keep *_selectedPATTausElectronCleaned_*_*',
            #'keep *_selectedPATTausElectronCleanedFlat_*_*',
            'keep *_selectedPATTausBoosted_*_*',
            'keep *_slimmedMETs_*_*',
            #'keep patPackedCandidates_packedPFCandidates_*_*',
            #'keep *_offlineSlimmedPrimaryVertices_*_*',
            #'keep *_ak4PFJetsElectronCleaned_*_*',
            #'keep *_ak4PFJetsMuonCleaned_*_*',
            #'keep *_ak4PFJets_*_*',
        )

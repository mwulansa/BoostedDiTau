import os
import ROOT
import sys
import copy
import re
from array import array
import math

from argparse import ArgumentParser
parser = ArgumentParser()

parser.add_argument('--files', type=str,
                  dest='files',
                  help='Input files')

parser.add_argument('--maxFiles', type=int, 
                  dest='maxFiles',
                  default=-1,
                  help='Maximum number of files to use')

parser.add_argument('--outname', type=str,
                  default='outplots.root',
                  dest='outname',
                  help='Name of output file')

parser.add_argument('--verbose', action='store_true',
                  default=False,
                  dest='verbose',
                  help='Print debugging info')

parser.add_argument('--correctJets', type=str,
                  default=None,
                  dest='correctJets',
                  help='Apply latest jet corrections. Specify JEC name based on https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC.')

parser.add_argument('--isData', action='store_true',
                  default=False,
                  dest='isData',
                  help='Run on data')


parser.add_argument('--smearJets', action='store_true',
                  default=False,
                  dest='smearJets',
                  help='Smear jet energy.')


parser.add_argument('--maxevents', type=int,
                  default=-1,
                  dest='maxevents',
                  help='Number of events to run. -1 is all events')

parser.add_argument('--maxjets', type=int,
                  default=999,
                  dest='maxjets',
                  help='Number of jets to plot. To plot all jets, set to a big number like 999')

parser.add_argument('--minAK8Pt', type=float,
                  default=200.,
                  dest='minAK8Pt',
                  help='Minimum PT for AK8 jets')

parser.add_argument('--maxAK8Rapidity', type=float,
                  default=2.4,
                  dest='maxAK8Rapidity',
                  help='Maximum AK8 rapidity')


parser.add_argument('--minAK4Pt', type=float,
                  default=0,
                  dest='minAK4Pt',
                  help='Minimum PT for AK4 jets')

parser.add_argument('--maxAK4Rapidity', type=float,
                  default=5.4,
                  dest='maxAK4Rapidity',
                  help='Maximum AK4 rapidity')

parser.add_argument('--xrootd', type=str,
                  #default='cmsxrootd.fnal.gov',
                  default = "xrootd-cms.infn.it",
                  dest='xrootd',
                  help='xrootd redirect string')

parser.add_argument('--matchPdgIdAK4', type=str, nargs=2,
                  dest='matchPdgIdAK4',
                  help='Perform truth matching of AK4 jets (specify PDG ID and DeltaR)')

parser.add_argument('--matchPdgIdAK8', type=str, nargs=2,
                  dest='matchPdgIdAK8',
                  help='Perform truth matching of AK8 jets (specify PDG ID and DeltaR)')


args = parser.parse_args()

#############################

from DataFormats.FWLite import Events, Handle

jethandle0  = Handle ("std::vector<pat::Jet>")
jetlabel0 = ("slimmedJets")

rhoHandle = Handle ("double")
rhoLabel = ("fixedGridRhoAll")

pvHandle = Handle("std::vector<reco::Vertex>")
pvLabel = ("offlineSlimmedPrimaryVertices")

#############
#if args.correctJets:

vPar = ROOT.vector(ROOT.JetCorrectorParameters)()
vPar.push_back( ROOT.JetCorrectorParameters('Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt'))
vPar.push_back( ROOT.JetCorrectorParameters('Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt'))
vPar.push_back( ROOT.JetCorrectorParameters('Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt'))
jec = ROOT.FactorizedJetCorrector( vPar )
jecUnc = ROOT.JetCorrectionUncertainty('Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt')
##########################3
####################################

f = ROOT.TFile(args.outname, "RECREATE")
f.cd()

h_nevents = ROOT.TH1F("h_nevents", "NEvents", 1, -0.5, 0.5)
h_ptAK4       = ROOT.TH1F("h_ptAK4", "Jet p_{T};p_{T} (GeV)", 2000, 0, 2000)
h_ptUpAK4     = ROOT.TH1F("h_ptUpAK4", "JEC Up Jet p_{T};p_{T} (GeV)", 2000, 0, 2000)
h_ptDownAK4   = ROOT.TH1F("h_ptDownAK4", "JEC Down Jet p_{T};p_{T} (GeV)", 2000, 0, 2000)
h_ptUncorrAK4 = ROOT.TH1F("h_ptUncorrAK4", "UnCorrected Jet p_{T};p_{T} (GeV)", 2000, 0, 2000)
h_JECValueAK4 = ROOT.TH1F("h_JECValueAK4", "Value of JEC for AK4 Jet", 100, 0.8, 1.1)
h_etaAK4      = ROOT.TH1F("h_etaAK4", "AK4 Jet #eta;#eta", 120, -6, 6)
h_yAK4        = ROOT.TH1F("h_yAK4", "AK4 Jet Rapidity;y", 120, -6, 6)
h_phiAK4      = ROOT.TH1F("h_phiAK4", "AK4 Jet #phi;#phi (radians)",100,-ROOT.Math.Pi(),ROOT.Math.Pi())
h_mAK4        = ROOT.TH1F("h_mAK4", "AK4 Jet Mass;Mass (GeV)", 100, 0, 1000)
h_areaAK4     = ROOT.TH1F("h_areaAK4", "AK4 Jet Area;Area", 250, 0, 5.0)

# h_ptAK4Gen   = ROOT.TH1F("h_ptAK4Gen", "AK4Gen Jet p_{T};p_{T} (GeV)", 300, 0, 3000)
# h_etaAK4Gen  = ROOT.TH1F("h_etaAK4Gen", "AK4Gen Jet #eta;#eta", 120, -6, 6)
# h_yAK4Gen    = ROOT.TH1F("h_yAK4Gen", "AK4Gen Jet Rapidity;y", 120, -6, 6)
# h_phiAK4Gen  = ROOT.TH1F("h_phiAK4Gen", "AK4Gen Jet #phi;#phi (radians)",100,-ROOT.Math.Pi(),ROOT.Math.Pi())
# h_mAK4Gen    = ROOT.TH1F("h_mAK4Gen", "AK4Gen Jet Mass;Mass (GeV)", 100, 0, 1000)
# h_areaAK4Gen = ROOT.TH1F("h_areaAK4Gen", "AK4Gen Jet Area;Area", 250, 0, 5.0)

varTree = ROOT.TTree("varTree", "varTree")
ak4pt           = array('f', [-1.])
ak4mass         = array('f', [-1.])

varTree.Branch('ak4pt', ak4pt, 'ak4pt/F')
varTree.Branch('ak4mass', ak4mass, 'ak4mass/F')

###########################################

filelist = file(args.files)
filesraw = filelist.readlines()
files = []
nevents = 0
print len(filesraw)
for i, ifile in enumerate(filesraw):
    if args.maxFiles >= 0:
#        print ifile
        files.append(ifile.rstrip())
#        print files
        if i >= args.maxFiles:
            break
    if len( ifile ) > 2 and ifile[:6] == "/store": 
        s = 'root://' + args.xrootd + '/' + ifile.rstrip()
        files.append( s )
        print 'Added ' + s
    elif len(ifile) > 2 and ifile[:4] == "/afs":
        files.append(ifile.rstrip())

for ifile in files :
    nevents = 0
    if args.maxevents > 0 and nevents > args.maxevents :
        break
    print 'Processing file ' + ifile
    events = Events(ifile)

    # loop over events in this file
    i = 0
    for event in events:
        if args.maxevents > 0 and nevents > args.maxevents :
            break
        i += 1
        nevents += 1
        h_nevents.Fill(0)

        if i % 1000 == 0 :
            print '    ---> Event ' + str(i)

            ######################################

            # get rho and vertices for JEC
        event.getByLabel (rhoLabel, rhoHandle)
        event.getByLabel (pvLabel, pvHandle)

        rhoValue = rhoHandle.product()
        pvs = pvHandle.product()
        
        # if doMatchingAK4 or doMatchingAK8:
        #     event.getByLabel(genParticlesLabel, genParticlesHandle)
        #     genParticles = genParticlesHandle.product()

        #     if doMatchingAK4:
        #         genPartsAK4 = []
        #         for genParticle in genParticles:
        #             if abs(genParticle.pdgId()) == matchAK4PdgId:
        #                 genPartsAK4.append(genParticle)

        if args.verbose :
            print '------ AK4 jets ------'
        # use getByLabel, just like in cmsRun
        event.getByLabel (jetlabel0, jethandle0)
        # get the product
        jets0 = jethandle0.product()
        # loop over jets and fill hists
        ijet = 0

        selected_jets=[]
        selected_uncorrjets=[]

        for jet in jets0 :
            if ijet >= args.maxjets :
                break

            if jet.pt() > args.minAK4Pt and abs(jet.rapidity()) < args.maxAK4Rapidity :

                uncorrJet = copy.copy( jet.correctedP4(0) ) # For some reason, in python this is interfering with jet.genJet() in strange ways without the copy.copy

                if uncorrJet.E() < 0.00001 :
                    print 'Very strange. Uncorrected jet E = ' + str( uncorrJet.E()) + ', but Corrected jet E = ' + str( jet.energy() )
                    continue

                # Apply loose jet ID to uncorrected jet  https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
                # the folling is valid only for eta < 2.4. see twiki for full definition
                # nhf = jet.neutralHadronEnergy() / uncorrJet.E()
                # nef = jet.neutralEmEnergy() / uncorrJet.E()
                # chf = jet.chargedHadronEnergy() / uncorrJet.E()
                # cef = jet.chargedEmEnergy() / uncorrJet.E()
                # nconstituents = jet.numberOfDaughters()
                # nch = jet.chargedMultiplicity()

                NHF  = jet.neutralHadronEnergyFraction()
                NEMF = jet.neutralEmEnergyFraction()
                CHF  = jet.chargedHadronEnergyFraction()
                MUF  = jet.muonEnergyFraction()
                CEMF = jet.chargedEmEnergyFraction()
                NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity()
                NumNeutralParticles =jet.neutralMultiplicity()
                CHM      = jet.chargedMultiplicity()

                # goodJet = \
                #   abs(jet.eta())<2.5 and \
                #   MUF < 0.8 and \
                #   CEMF < 0.9 and \
                #   (NHF < 0.90 and NEMF < 0.90 and NumConst>1) and \
                #   ((abs(jet.eta())<=2.4 and CHF > 0.00 and CHM > 0.00 and CEMF < 0.99) \
                #   or abs(jet.eta())>2.4) and abs(jet.eta())<=2.7

                # if not goodJet :
                #     continue

                if MUF > 0.8: continue
                if CEMF > 0.9: continue
                if (NHF<0.90 and NEMF<0.90 and NumConst>1) and ((abs(uncorrJet.eta())<=2.4 and CHF>0 and CHM>0 and CEMF<0.99) or abs(uncorrJet.eta())>2.4) and abs(uncorrJet.eta())<=2.7:

                    selected_uncorrjets+=[uncorrJet]
                    selected_jets+=[jet]

        if len(selected_jets)>0 :selected_jets.sort(key=lambda x: x.pt(), reverse=True)
        if len(selected_uncorrjets)>0 :selected_jets.sort(key=lambda x: x.pt(), reverse=True)

        corr = 1.0
        corrUp = 1.0
        corrDn = 1.0

        #if args.correctJets : 
        jec.setJetEta( selected_uncorrjets[0].eta() )
        jec.setJetPt ( selected_uncorrjets[0].pt() )
        jec.setJetE  ( selected_uncorrjets[0].energy() )
        jec.setJetA  ( selected_jets[0].jetArea() )
        jec.setRho   ( rhoValue[0] )
        jec.setNPV   ( len(pvs) )
        icorr = jec.getCorrection()
        corr *= icorr
        corrUp *= icorr
        corrDn *= icorr

        #JEC Uncertainty
        jecUnc.setJetEta( selected_uncorrjets[0].eta() )
        jecUnc.setJetPhi( selected_uncorrjets[0].phi() )
        jecUnc.setJetPt( corr * selected_uncorrjets[0].pt() )
        corrUp += jecUnc.getUncertainty(1)
        jecUnc.setJetEta( selected_uncorrjets[0].eta() )
        jecUnc.setJetPhi( selected_uncorrjets[0].phi() )
        jecUnc.setJetPt( corr * selected_uncorrjets[0].pt() )
        corrDn -= jecUnc.getUncertainty(0)

################

#                    if (corr * uncorrJet.pt()) > 20 and abs(jet.eta())<2.5 :
        h_ptAK4.Fill( corr * selected_uncorrjets[0].pt() )
        h_JECValueAK4.Fill( corr )
        h_ptUncorrAK4.Fill( selected_uncorrjets[0].pt() )
        h_ptDownAK4.Fill( corrDn * selected_uncorrjets[0].pt() )
        h_ptUpAK4.Fill( corrUp * selected_uncorrjets[0].pt() )   

        h_yAK4.Fill( selected_jets[0].y() )
        h_phiAK4.Fill( selected_jets[0].phi() )
        h_mAK4.Fill( selected_jets[0].mass() )
        h_areaAK4.Fill( selected_jets[0].jetArea() )

#            ak4mass[0] = jet.mass()
#            ak4pt[0] = jet.pt()
#            varTree.Fill()

            # genJet = jet.genJet()
            # if genJet != None :
            #     h_ptAK4Gen.Fill( genJet.pt() )
            #     h_etaAK4Gen.Fill( genJet.eta() )
            #     h_yAK4Gen.Fill( genJet.y() )
            #     h_phiAK4Gen.Fill( genJet.phi() )
            #     h_mAK4Gen.Fill( genJet.mass() )
            #     h_areaAK4Gen.Fill( genJet.jetArea() )
            # if args.verbose == True : 
            #     print ("Jet {0:4.0f}, orig pt = {1:10.2f}, eta = {2:6.2f}, phi = {3:6.2f}, m = {4:6.2f}, " +
            #            "nda = {5:3.0f}, vtxmass = {6:6.2f}, area = {7:6.2f}, corr = {8:6.3f} +{9:6.3f} -{10:6.3f} ").format(
            #         ijet, jet.pt(), jet.eta(), jet.phi(), jet.mass(), jet.numberOfDaughters(), jet.userFloat('vtxMass'),
            #         jet.jetArea(), corr, abs(corrUp - corr), abs(corr - corrDn)
            #         ),
            #     if genJet != None :
            #         print (", gen pt = {0:6.2f}").format( genJet.pt() )
            #     else :
            #         print ''

f.cd()
f.Write()
f.Close()

import math
from ROOT import TMath

def dEtaInSeed(ele):
    if ele.superCluster().isNonnull() and ele.superCluster().seed().isNonnull():
        return ele.deltaEtaSuperClusterTrackAtVtx()-ele.superCluster().eta()+ele.superCluster().seed().eta()
    else:
        return sys.float_info.max

def GsfEleEInverseMinusPInverse(ele):
    ecal_energy_inverse = 1.0/ele.ecalEnergy()
    eSCoverP = ele.eSuperClusterOverP()
    return abs(1.0 - eSCoverP)*ecal_energy_inverse

def GsfEleMissingHitsCut(ele):
    mHits = ele.gsfTrack().hitPattern().numberOfAllHits(1)
    return mHits

def GsfEleConversionVetoCut(ele, bs, convs):
    if not bs==None or not convs==None:
        return not hasMatchedConversion(ele,convs,bs.position())
    else:
        print "Couldn't find a necessary collection, returning true!"
        return True

def hasMatchedConversion(ele,convs,bsp):
    for conv in convs:
        if not matchesConversion(ele, conv): continue
        if not isGoodConversion(conv, bsp): continue
        return True
    return False

def matchesConversion(ele, conv):
    allowCkfMatch=True
    for track in conv.tracks():
        print ele.gsfTrack().id().id(), track.id().id(), ele.closestCtfTrackRef().id().id()
        if ele.gsfTrack().isNonnull() and ele.gsfTrack().id().id()==track.id().id() and ele.gsfTrack().key() == track.key():
            return True
        elif allowCkfMatch and ele.closestCtfTrackRef().isNonnull() and ele.closestCtfTrackRef().id().id()==track.id().id() and ele.closestCtfTrackRef().key() == track.key():
            return True
    return False

def isGoodConversion(conv, bsp):
    #print mom, vtx.chi2()
    probMin=1e-6
    lxyMin=2.0
    nHitsBeforeVtxMax=1
    vtx=conv.conversionVertex()
    if vtx==None: return False
    if TMath.Prob(vtx.chi2(), int(vtx.ndof()))<probMin: return False

    mom=conv.refittedPairMomentum()
    dbsx = vtx.x() - bsp.x()
    dbsy = vtx.y() - bsp.y()
    lxy = (mom.x()*dbsx + mom.y()*dbsy)/mom.rho()

    if lxy<lxyMin: return False

    for i in range(conv.nHitsBeforeVtx()):
        if i>nHitsBeforeVtxMax: return False
        
    return True

def GsfEleEffAreaPFIsoCut(ele, rho):
    Etamin = [0.,1.000,1.4790,2.000,2.2000,2.3000,2.4000]
    Etamax = [1.0000,1.4790,2.0000,2.2000,2.3000,2.4000,5.000]
    EAs = [0.1703,0.1715,0.1213,0.1230,0.1635,0.1937,0.2393]

    for i in range(len(EAs)):
        if abs(ele.eta())>=Etamin[i] and abs(ele.eta())<Etamax[i]:
            effA=EAs[i]

    pfIso=ele.pfIsolationVariables()

    chad = pfIso.sumChargedHadronPt
    nhad = pfIso.sumNeutralHadronEt
    pho = pfIso.sumPhotonEt

    #print abs(ele.eta()), effA
    iso = chad + max(0.0, nhad + pho-effA*rho)
        
    iso /= ele.pt()

    return iso

def muonIsoCut(muon):
    return (muon.pfIsolationR04().sumChargedHadronPt+max(0,muon.pfIsolationR04().sumPhotonEt+muon.pfIsolationR04().sumNeutralHadronEt-0.5*muon.pfIsolationR04().sumPUPt))/muon.pt()

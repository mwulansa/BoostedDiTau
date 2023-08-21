import os, sys, math
import ROOT

# Cross section calc----------------

pi = math.pi

c_t = 1
m_t = 1.77686
fa = 1000

m_a = [10,30,50]
xsec_gt = [1.0256092e-07, 3.133914e-07, 5.072265e-07]

Br = [0, 6.8 , 9.7 , 5.7 , 6.2 , 3.0 , 3.0 , 9.7 , 0.88 , 1.9 , 1.8 , 2.1 , 2.9] #so Model 1 is M[_][1]

xsec_M = [[],[],[]]

for j in range (3):
    for i in range (13):
        Gamma_tt = ((c_t**2*m_a[j]/(8*pi))*m_t**2*math.sqrt(1-4*m_t**2/m_a[j]**2))/fa**2
        xsec_g_TCP = xsec_gt[j] / Gamma_tt
        xsec_M_list = xsec_g_TCP * Br[i]
        xsec_M[j].append(xsec_M_list)
        print(xsec_M[j][i])

xsec_TCP10 = xsec_M[0][1]
xsec_TCP30 = xsec_M[1][1]
xsec_TCP50 = xsec_M[2][1]

print("xsec_TCP10", xsec_TCP10)
print("xsec_TCP30", xsec_TCP30)
print("xsec_TCP50", xsec_TCP50)

#----------------------------------

xsecs={
    'TCP_m10':{'':xsec_TCP10},
    'TCP_m30':{'':xsec_TCP30},
    'TCP_m50':{'':xsec_TCP50},
    'DYJetsToLL':{
        'M-1To5_HT-600toInf':16.72*(1.23/1.375),
        'M-1To5_HT-400to600':65.9*(1.23/1.375),
        'M-1To5_HT-200to400':789.8*(1.23/1.375),
        'M-1To5_HT-150to200':1124.0*(1.23/1.375), 
#        'M-1To5_HT-600toInf':16.72*1.23,
#        'M-1To5_HT-400to600':65.9*123,
#        'M-1To5_HT-200to400':789.8*1.23,
#        'M-1To5_HT-150to200':1124*1.23, 
        'M-5to50_HT-70to100':301.0*1.23,
        'M-5to50_HT-100to200':224.4*1.23, 
        'M-5to50_HT-200to400':37.87*1.23, 
        'M-5to50_HT-400to600':3.628*1.23, 
        'M-5to50_HT-600toInf':1.107*1.23, 
        'M-50_HT-70to100':169.9*1.23,
        'M-50_HT-100to200':147.4*1.23, 
        'M-50_HT-200to400':40.99*1.23, 
        'M-50_HT-400to600':5.678*1.23, 
        'M-50_HT-600to800':1.367*1.23, 
        'M-50_HT-800to1200':0.6304*1.23, 
        'M-50_HT-1200to2500':0.1514*1.23, 
        'M-50_HT-2500toInf':0.003565*1.23},
    'DYJetsToQQ':{
        'HT180':1208},
    'TTJets':{
        'Dilept':87.315},
    'ST':{
        's-channel_4f_leptonDecays':3.36,
        't-channel_antitop_4f_inclusiveDecays':26.38,
        't-channel_top_4f_inclusiveDecays':44.33,
        'tW_antitop_5f_inclusiveDecays':35.85,
        'tW_top_5f_inclusiveDecays':35.85},
    'Diboson':{
        'WW':118.7,
        'WZ':47.13,
        'ZZ':16.523},
    'QCD':{
        'HT50to100':246300000,
        'HT100to200':27990000,
        'HT200to300':1712000,
        'HT300to500':347700,
        'HT500to700':32100,
        'HT700to1000':6831,
        'HT1000to1500':1207,
        'HT1500to2000':119.9,
        'HT2000toInf':25.24},
    'WJetsToQQ':{
        'HT-600ToInf':99.65},
    'ZJetsToQQ':{
        'HT600toInf':581.9}
}

def weightBackgroundHists(hists, files, version, reco, var):
    for sample in list(xsecs):
        for mass in list(xsecs[sample]):
            if gen==True:
                filename="h_Gen_"+sample+"_"+mass+"_"+version+".root"
            else:
                filename="h_"+sample+"_"+mass+"_"+version+".root"          
            histname="h"+reco+"_"+var 
            print(filename, histname)
            fil=ROOT.TFile(filename, 'r')
            if fil.IsZombie(): continue
            files+=[fil]
            hist=fil.Get(histname) 
            nEvt=fil.Get("hNEvent").GetBinContent(2)
            xsec=xsecs[sample][mass] 
            if list(xsecs[sample]).index(mass)==0:  
                h = hist.Clone(sample+"_"+histname)
                print(h)
                if scaling == 'xsection':
                    h.Scale((xsec/nEvt)*L)
                else: 
                    scale = 1/(h.Integral())
                    h.Scale(scale)
                hists[sample] = h 
            else: 
                if scaling == 'xsection': 
                    hist.Scale((xsec/nEvt)*L)
                else:
                    hist.Scale(scale)
                hists[sample].Add(hist)
            
hists = {}
files = []

version = 'vBT-7'

reco ='BTau'

L = 35900

#-------------

gen = False
scaling = 'xsection' #'xsection' else: 'unity'

#-------------

Sample = ['DYJetsToLL','TCP_m10','QCD','TTJets','ZJetsToQQ'] #sample names

histlist = ['Metcut_M_nA_veto','dR_M_METPt_nA_veto','Trig_M_nA_veto','Baseline_M_nA_veto'] #histograms to plot

for var in histlist:
    weightBackgroundHists(hists, files, version, reco, var)

    if scaling == 'xsection':
         out = ROOT.TFile("h_"+reco+"_"+var+"_"+version+".root",'recreate')
    else: 
         out = ROOT.TFile("h_"+reco+"_"+var+"_"+version+"_unity.root",'recreate')

    out.cd()

    for name in Sample:
#        hists[name].Rebin(5)                                                                                                           

        if name=='DYJetsToLL':
            hists[name].SetFillColor(ROOT.kRed)
        if name == 'TCP':
            hists[name].SetFillStyle(3335)
        elif name=='TTJets':
            hists[name].SetFillColor(ROOT.kOrange)
        elif name=='ST':
            hists[name].SetFillColor(ROOT.kMagenta)
        elif name == 'ZJetsToQQ':
            hists[name].SetFillColor(ROOT.kBlue)
        elif name == 'QCD':
            hists[name].SetFillColor(ROOT.kGray+3)
        elif name=='WJetsToQQ':
            hists[name].SetFillColor(ROOT.kGreen+3)

        hists[name].Write()
        print(name)

    out.Close()

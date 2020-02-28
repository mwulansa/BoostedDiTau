import ROOT
import math

def SetStyle1D():
    ROOT.gROOT.ForceStyle()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadBottomMargin(0.11)
    ROOT.gStyle.SetPadTopMargin(0.09)
    ROOT.gStyle.SetPadRightMargin(0.07)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    return

def SetStyle2D():
    ROOT.gROOT.ForceStyle()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPadLeftMargin(0.1)
    ROOT.gStyle.SetPadBottomMargin(0.15)
    ROOT.gStyle.SetPadTopMargin(0.05)
    ROOT.gStyle.SetPadRightMargin(0.2)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    return

def SetHistStyle(h,color,lineStyle,lineWidth,markerSize):
    h.SetLineColor(color)
    h.SetLineStyle(lineStyle)
    h.SetLineWidth(lineWidth)
    h.SetMarkerSize(markerSize)
    return

def SetupAxis(h, whichAxis, Ndivisions, labelSize, labelOffset, title, titleSize, titleOffset, pSet=None):
    # initial Values
    if Ndivisions=="":
        Ndivisions=510
    if labelSize=="Default":
        labelSize=0.05
    if labelOffset=="Default":
        labelOffset=0.005
    if titleSize=="Default":
        titleSize=0.06
    if titleOffset=="Default":
        titleOffset=1.2
    if whichAxis=="X":
        h.GetXaxis().SetNdivisions(Ndivisions)
        h.GetXaxis().SetLabelSize(labelSize)
        h.GetXaxis().SetLabelOffset(labelOffset)
        h.GetXaxis().SetTitle(title)
        h.GetXaxis().SetTitleSize(titleSize)
        h.GetXaxis().SetTitleOffset(titleOffset)
    if whichAxis=="Y":
        h.GetYaxis().SetNdivisions(Ndivisions)
        h.GetYaxis().SetLabelSize(labelSize)
        h.GetYaxis().SetLabelOffset(labelOffset)
        h.GetYaxis().SetTitle(title)
        h.GetYaxis().SetTitleSize(titleSize)
        h.GetYaxis().SetTitleOffset(titleOffset)
    if whichAxis=="Z":
        h.GetZaxis().SetNdivisions(Ndivisions)
        h.GetZaxis().SetLabelSize(labelSize)
        h.GetZaxis().SetLabelOffset(labelOffset)
        h.GetZaxis().SetTitle(title)
        h.GetZaxis().SetTitleSize(titleSize)
        h.GetZaxis().SetTitleOffset(titleOffset)
    return
        
def SetupLegend(l, testSize):
    l.SetTextSize(0.04)
    l.SetFillStyle(0)
    return

def SetupRatioPlotPads():
    c = ROOT.TCanvas("Comparisons", "Comparisons", 0, 0, 650, 800)
    c.Divide(1,1)
    c.cd()
    padDown=ROOT.TPad("","",0,0,1,0.3)
    padDown.SetTopMargin(0.08)
    padDown.SetBottomMargin(0.25)
    #padDown.Draw()
    c.cd()
    padUp=ROOT.TPad("","",0,0.3,1,1)
    padUp.SetBottomMargin(0.001)
    padUp.SetTopMargin(0.07)
    #padUp.Draw()
    return c, padDown, padUp

RatioPlotTopXaxisSetup=["X", "", 0, 0, "", 0, 0]
RatioPlotTopYaxisSetup=["Y", "", 0.04, 0.005, "N", 0.05, 1.2]
RatioPlotBottomXaxisSetup=["X", "", 0.04/(3./6.5), 0.008, "Energy [GeV]", 0.05/(3./6.5), 1.]
RatioPlotBottomYaxisSetup=["Y", 505, 0.04/(3./6.5), 0.008, "Mixing/PreMixing", 0.05/(3./6.5), 0.5]

def SetRatioPlotTick(h):
    h.GetXaxis().SetTickLength(0.07)
    h.GetYaxis().SetTickLength(0.04)

def sumHistTimesBinCenter(h):
    nbins=h.GetNbinsX()
    sum=0
    for nbin in range(nbins):
        sum+=h.GetBinContent(nbin+1)*h.GetBinCenter(nbin+1)
    return sum

def EToET(E, eta):
    ET=E*math.sin(2*math.atan(math.exp(-eta)))
    return ET

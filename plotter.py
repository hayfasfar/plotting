from style import *
import sys
import style
import ROOT
import math
import os 
import json
import yaml

# path to processed nanoAOD ntuples
ntuple_path = "/vols/cms/hsfar/nanoAOD_friends/01Apr20/"
lumi = 35.88

def find_xsec(path, xsecs):
    for key, val in xsecs.items():
        if key in path:     
            return val

# Read in event yields and cross-sections for normalisation
with open("eventyields.json") as json_file:
    yields = json.load(json_file)

with open("xsecs.yaml") as yaml_file:
    xsecs = yaml.load(yaml_file, Loader=yaml.FullLoader)

# Define categories by cuts for plotting

categories = {}
#categories["mu1mu2jet_DY_1T1L_CR_"] = "(ntightMuons == 1)*(nlooseMuons == 1)*(dimuon_mass > 85. and dimuon_mass <= 110.)*(MET_pt < 50.)"

#categories["mu1mu2jet_DY_2tights_CR_"] = "(ntightMuons == 2)*(dimuon_mass > 85.)*(nselectedJets>=2)"

#categories["mu1mu2jet_Wjets_CR_"] = "(ntightMuons == 1)*(nlooseMuons == 1)*(dimuon_mass > 110. or ( dimuon_mass < 80. and dimuon_mass > 15.))*(Jet_Muon_MET_Mu_mT >=60.)*(MET_pt < 100.)*(nselectedJets<6)*(Jet_Muon_ht <=180.)*(Jet_Muon_minPhi >1.)*(lepJet_deltaR>=1.)"
categories["mu1mu2jet_TT_CR_"] ="(ntightMuon == 1)*(nlooseMuons == 1)*(nselectedJets>=8)"

#categories["mu1mu2jet_QCD_CR_"] = "(ntightMuons == 1)*(nlooseMuons == 1) * (dimuon_mass < 10.) *(Jet_Muon_minPhi<0.8) * (MET_pt < 50.) " 

# This class is responsible for making the histogram and plotting it for a given variable
class Variable:
    def __init__(self, varexp, name, nbins, xmin, xmax, jesUp=None, jesDown=None, jerUp=None, jerDown=None, sysUnc=True, logy=True):
        self.varexp = varexp
        self.args = (varexp, varexp, nbins, xmin, xmax,)
        self.stack = ROOT.THStack(varexp, varexp)
        self.sumMC = ROOT.TH1F(name, name, nbins, xmin, xmax)
	self.sumUncUp = ROOT.TH1F(name+"Up", name+"Up", nbins, xmin, xmax)
	self.sumUncDown = ROOT.TH1F(name+"Down",name+"Down", nbins, xmin, xmax)
        self.signals = []
        self.data = None
        self.name = name
        self.logy = logy
        self.leg = makeLegend(0.5,0.6,0.89,0.89)
        self.leg.SetTextSize(self.leg.GetTextSize()*0.8)
        self.xmin = xmin
        self.xmax = xmax
	self.jesUp = jesUp 
	self.jesDown = jesDown
	self.jerUp = jerUp
	self.jerDown = jerDown 
	self.sysUnc = sysUnc
    def Add(self, hist, title, isSignal=False, isData=False):
        hist[-1].SetDirectory(0)
        if isSignal:
            self.signals.append(hist[-1])
            hist.SetLineStyle(len(self.signals))
            self.leg.AddEntry(hist[-1], title, "l")
        elif isData:
            self.data = hist[-1]
            self.leg.AddEntry(hist[-1], title, "p")
        else:
	    
	    self.stack.Add(hist[-1])
            self.sumMC.Add(hist[-1])
            self.leg.AddEntry(hist[-1], title, "f")

	    if self.sysUnc  :
        	hist[-3].SetDirectory(0)
		self.sumUncUp.Add(hist[-3])
        	hist[-2].SetDirectory(0)
		self.sumUncDown.Add(hist[-2])
	    
    def Draw(self, suffix, opt):
        canvas = makeCanvas(name=self.varexp)
        upperPad = ROOT.TPad("upperPad", "upperPad", 0, 0.33, 1, 1)
        lowerPad = ROOT.TPad("lowerPad", "lowerPad", 0, 0, 1, 0.33)
        upperPad.SetBottomMargin(0.00001)
        upperPad.SetBorderMode(0)
        upperPad.SetTopMargin(0.15)
        lowerPad.SetTopMargin(0.00001)
        lowerPad.SetBottomMargin(0.4)
        lowerPad.SetBorderMode(0)
        canvas.SetBottomMargin(0.2)
        canvas.SetTopMargin(0.1)
        upperPad.Draw()
        lowerPad.Draw()
        upperPad.cd()

	## Draw MC histograms which are stacked together in self.stack. 

        self.stack.Draw(opt)
        self.stack.SetMinimum(1)
        if self.logy:
            self.stack.SetMaximum(self.stack.GetMaximum()*1000)
            upperPad.SetLogy()
        else:
            self.stack.SetMaximum(self.stack.GetMaximum()*1.6)

	## Draw Signal if it is signal MC comparison study. otherwise the array self.signals will be empty.
        for signal in self.signals:
            signal.SetFillStyle(0)
            signal.Draw("HIST SAME")

	### now draw data if it is CR Data study. 
        self.data.Draw("P SAME")
        self.hist_ratio = self.data.Clone("ratio histogram")
        self.hist_ratio.Divide(self.sumMC)
	## add SysUnc here
 
	## Clone data first. you need to divide data by observable SystUp. 
	#self.histSysUp_ratio = self.data.Clone("ratio histogram")
	#self.histSysUpClone = self.sumUncUp.Clone("copy up")
	#self.histSysUp_ratio.Divide(self.sumUncUp)
	self.histSysUp_ratio = self.sumUncUp.Clone("ratio plot")
	self.histSysUp_ratio.Divide(self.sumMC) 
	
	self.histSysDown_ratio = self.sumUncDown.Clone("ratio plot")	
	self.histSysDown_ratio.Divide(self.sumMC)
	#self.histSysDown_ratio = self.data.Clone()
	#self.histSysDown_ratio.Divide(self.sumUncDown)
	lowerPad.cd()
        self.hist_ratio.SetMinimum(0.8)
        self.hist_ratio.SetMaximum(1.2)
        self.hist_ratio.GetXaxis().SetTitle(self.name)
        self.hist_ratio.GetXaxis().SetTitleOffset(2.5)
        self.hist_ratio.Draw("P")
	self.histSysUp_ratio.SetFillStyle(0)
	self.histSysUp_ratio.SetLineColor(4)
	#self.histSysUp_ratio.SetFillColorAlpha(ROOT.kBlue, 0.35)
	self.histSysUp_ratio.Draw("HIST SAME")
	self.histSysDown_ratio.SetFillStyle(0)
	self.histSysDown_ratio.SetLineColor(2)
	#self.histSysDown_ratio.SetFillColorAlpha(ROOT.kBlue, 0.35)
	self.histSysDown_ratio.Draw("HIST SAME")
	rootObj = []
        line = ROOT.TLine(self.xmin, 1, self.xmax, 1)
        line.Draw("SAME")
	for ibin in range(self.sumMC.GetNbinsX()):
            c = self.sumMC.GetBinCenter(ibin+1)
            w = self.sumMC.GetBinWidth(ibin+1)
            m = self.sumMC.GetBinContent(ibin+1)
            if m > 0.0:
                h = min(self.sumMC.GetBinError(ibin+1)/m, 0.399)
                box = ROOT.TBox(c-0.5*w, 1-h, c+0.5*w, 1+h)
                box.SetFillStyle(3345)
                box.SetLineColor(ROOT.kGray+1)
                box.SetFillColor(ROOT.kGray)
                rootObj.append(box)
                box.Draw("SameF")
                box2 = ROOT.TBox(c-0.5*w, 1-h, c+0.5*w, 1+h)
                box2.SetFillStyle(0)
                box2.SetLineColor(ROOT.kGray+1)
                box2.SetFillColor(ROOT.kGray)
                rootObj.append(box2)
                box2.Draw("SameL")
        canvas.cd()
        self.leg.Draw("SAME")
        makeCMSText(0.13, 0.97,additionalText="Simulation Preliminary")
        makeText(0.13, 0.78, 0.4, 0.78, "nselectedJets>= 8")
        makeLumiText(0.7, 0.97)
        canvas.SaveAs("plots/"+suffix+self.varexp+"nselectedJets8.pdf")

# This class prepares a given sample by scaling to int. luminosity
class Sample:
    def __init__(self, name, paths, isMC=True):
        self.paths = paths
        self.name = name
        self.file_list = ROOT.std.vector('string')()
        self.sum_weight = 0
        self.isMC = isMC
        for path in self.paths:
	    #print "path.join is : ", os.path.join(ntuple_path, path)  
            for f in os.listdir(os.path.join(ntuple_path, path)):
                self.file_list.push_back(os.path.join(ntuple_path, path, f))
            if self.isMC:
                self.sum_weight += yields[path]["weighted"]
        self.rdf = ROOT.RDataFrame("Friends", self.file_list)
        if self.isMC:
            self.rdf = self.rdf.Define("weightLumi", "genweight*puweight*tightMuon_weight_iso_nominal*tightMuon_weight_id_nominal*%s*1000.0*%s/%s" %(lumi, find_xsec(path, xsecs), self.sum_weight))
        else:
            self.rdf = self.rdf.Define("weightLumi", "1")
        for category, weight in categories.items():
            self.rdf = self.rdf.Define(category, "weightLumi*%s" %(weight))
        self.hists = []
        #print("RDF "+name+ " has entries: "+str(self.rdf.Count().GetValue()))

# A process is a combination of several "Samples" which are all added up internally
class Process:
    def __init__(self, name, title, color):
        self.name = name
        self.title = title
        self.color = color	
        self.hists = []
        self.rdfs = []

    def add(self, *args):
        for arg in args:
            self.rdfs.append(arg.rdf)


    def Histo1D(self, args, varexp , jesUp, jesDown, jerUp, jerDown, sysUnc, weight):
        for i, rdf in enumerate(self.rdfs):
	    var = varexp.split("[")[0]
            if i == 0:
		hist = rdf.Define("var",varexp).Histo1D(args, "var", weight)
                #hist = rdf.Histo1D(args, varexp, weight)
		if sysUnc :
			histjesUp = rdf.Histo1D(args, jesUp, weight)
			histjesDown = rdf.Histo1D(args, jesDown , weight)
			histjerUp = rdf.Histo1D(args, jerUp , weight)
			histjerDown = rdf.Histo1D(args, jerDown , weight)
            else:
                #tmp_hist = rdf.Histo1D(args, varexp, weight)
		tmp_hist = rdf.Define("var" , varexp).Histo1D(args, "var", weight)
                hist.Add(tmp_hist.GetValue())
		if sysUnc: 
			tmp_histjesUp = rdf.Histo1D(args, jesUp, weight)
			tmp_histjesDown = rdf.Histo1D(args, jesDown , weight)
			tmp_histjerUp = rdf.Histo1D(args, jerUp, weight)
			tmp_histjerDown = rdf.Histo1D(args, jerDown , weight)
	    		histjesUp.Add(tmp_hist.GetValue())
			histjesDown.Add(tmp_hist.GetValue())
	    		histjerUp.Add(tmp_hist.GetValue())
			histjerDown.Add(tmp_hist.GetValue())
	if sysUnc : 
	        histTotUp = hist.Clone() 
		histTotDown = hist.Clone() 			
		for i in range(1,hist.GetNbinsX()+1) :
			bin_nominal = hist.GetBinContent(i) 
			bin_jesUp = histjesUp.GetBinContent(i) 
			bin_jesDown = histjesDown.GetBinContent(i)
			bin_jerUp = histjerUp.GetBinContent(i) 
			bin_jerDown = histjerDown.GetBinContent(i)

			

			eps_jesUp = bin_jesUp - bin_nominal
			eps_jesDown = bin_jesDown - bin_nominal  
			eps_jerUp =  bin_jerUp - bin_nominal
                	eps_jerDown = bin_jerDown - bin_nominal
			epJesAvr = (eps_jesUp + eps_jesDown)/2.
			epJerAvr = (eps_jesDown + eps_jerDown)/2.

			sysTotalAvr =  math.sqrt(epJesAvr*epJesAvr + eps_jesDown*eps_jesDown)
			histTotUp.SetBinContent(i , bin_nominal  + sysTotalAvr)
                	histTotDown.SetBinContent(i , bin_nominal  - sysTotalAvr)

		
	if sysUnc :
		histTotUp.SetLineWidth(2)
		histTotDown.SetLineWidth(2) 
		self.hists.append(histTotUp.Clone())
		self.hists.append(histTotDown.Clone())
        hist.GetXaxis().SetTitle(args[1])
        hist.SetLineColor(ROOT.TColor.GetColor(colorscale(self.color, 0.8)))
        hist.SetFillColor(ROOT.TColor.GetColor(colorscale(self.color, 1)))
        hist.SetLineWidth(2)
        self.hists.append(hist.Clone())
	return self.hists

# self.hists return a TH1F histogram nselectedJets.

w0jets = Sample("w0jets", ["WToLNu_0J_13TeV-amcatnloFXFX-pythia8-2016"])
w1jets = Sample("w1jets", ["WToLNu_1J_13TeV-amcatnloFXFX-pythia8-2016"])
w2jets = Sample("w2jets", ["WToLNu_2J_13TeV-amcatnloFXFX-pythia8-ext4-2016"])
dy10to50 = Sample("dy10to50",  ["DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-2016"])
dy50 = Sample("dy50",  ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-ext2-2016"])
ttsemilep = Sample("ttsemilep", ["TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8-2016"])
TTo2L2Nu = Sample("TTo2L2Nu", ["TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8-2016"])
qcd_15to20 = Sample("qcd_15to20", ["QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"]) 
qcd_20to30 = Sample("qcd_20to30", ["QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_30to50 = Sample("qcd_30to50", ["QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_50to80 = Sample("qcd_50to80", ["QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_80to120 = Sample("qcd_80to120", ["QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_120to170 = Sample("qcd_120to170", ["QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_170to300 = Sample("qcd_170to300", ["QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia-2016"])
qcd_170to300 = Sample("qcd_170to300", ["QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia-2016"])
qcd_300to470 = Sample("qcd_300to470", ["QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_470to600 = Sample("qcd_470to600", ["QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_600to800 = Sample("qcd_600to800", ["QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_800to1000 = Sample("qcd_800to1000", ["QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
qcd_1000toInf = Sample("qcd_1000toInf", ["QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8-2016"])
hnlM4_V0p00183575597507 = Sample("HNL2", ["HeavyNeutrino_lljj_M-4_V-0_00183575597507_mu_Dirac_Moriond17_aug2018_miniAODv3-2016"])
hnlM8_V0p000415932686862 = Sample("HNL", ["HeavyNeutrino_lljj_M-8_V-0_000415932686862_mu_Dirac_Moriond17_aug2018_miniAODv3-2016"])

runb = Sample("Run2016B", ["SingleMuon_Run2016B_ver2"], isMC=False)
runc = Sample("Run2016C", ["SingleMuon_Run2016C"], isMC=False)
rund = Sample("Run2016D", ["SingleMuon_Run2016D"], isMC=False)
rune = Sample("Run2016E", ["SingleMuon_Run2016E"], isMC=False)
runf = Sample("Run2016F", ["SingleMuon_Run2016F"], isMC=False)
rung = Sample("Run2016G", ["SingleMuon_Run2016G"], isMC=False)
runh = Sample("Run2016H", ["SingleMuon_Run2016H"], isMC=False)

wjets = Process("W+Jets", "W+jets", "#388e3c")
wjets.add(w0jets, w1jets, w2jets)
dyjets = Process("DY+Jets", "DY+Jets", "#1976d2")
dyjets.add(dy10to50, dy50)
tt = Process("ttbar", "t#bar{t}",  "#ef5350")
tt.add(ttsemilep, TTo2L2Nu)
hnl1 = Process("HNL1", "m_{N} = 8 GeV, |V_{#mu}|^{2} = 4.2#times10^{-7}", "#087858")
hnl1.add(hnlM8_V0p000415932686862)
hnl2 = Process("HNL2", "m_{N} = 4 GeV, |V_{#mu}|^{2} = 3.4#times10^{-6}", "#087858")
hnl2.add(hnlM4_V0p00183575597507)
qcd = Process("qcd", "QCD", "#bdbdbd")
qcd.add(qcd_15to20, qcd_20to30, qcd_30to50, qcd_50to80, qcd_80to120, qcd_120to170, qcd_170to300, qcd_300to470, qcd_470to600, qcd_600to800, qcd_800to1000, qcd_1000toInf)

data = Process("data", "data", "#000000")
data.add(runb, runc, rund, rune, runf, rung, runh)
processes = [ wjets, qcd , tt , dyjets, data]
#processes = [ qcd , tt , data]

# variable is 1 line in the variable.yaml file. 
# suffix is for example "mu1mu2jet_TT_CR_" 
#variable.varexp is : 
#variable.args should be the histogram arguments. 
 
for suffix, weight in categories.items():
   print(weight)
   with open("variables.yaml") as yaml_file:
        variables = yaml.load(yaml_file, Loader=yaml.FullLoader)
   for index , variable in enumerate(variables):
     if index == int(sys.argv[1]) : 
        print(variable)
        variable = Variable(*variable)
        print "variable args is : ", variable.args, " variable.varexp is :  ", variable.varexp
        for process in processes:
            print(process.name)
            if "HNL" in process.name:
                isSignal = True
            else:
                isSignal = False
            if process.name == "data":
                isData = True
		variable.sysUnc = False
            else:
                isData = False
	    print "suffix is :", suffix
	    #histo1D should return 3 histograms. 
	    variable.Add(process.Histo1D(variable.args, variable.varexp,variable.jesUp, variable.jesDown, variable.jerUp, variable.jerDown,variable.sysUnc, suffix), process.title, isSignal=isSignal, isData=isData)
            #variable.Add(process.Histo1D(variable.args, variable.varexp, suffix), process.title, isSignal=isSignal, isData=isData)
        variable.Draw(suffix, "hist")

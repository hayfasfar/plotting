from style import *
import sys
import style
import ROOT
import math
import os 
import json
import yaml

# path to processed nanoAOD ntuples
ntuple_path = "/home/hep/hsfar/private/nanoTools/CMSSW_10_2_18/src/nanoAOD_friends_ptJets30/HNL/18Mar20"
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
categories["mu1mu2jet_TT_CR_"] ="(ntightMuons == 1)*(nlooseMuons == 1)*(MET_pt > 100)"

#categories["mu1mu2jet_QCD_CR_"] = "(ntightMuons == 1)*(nlooseMuons == 1) * (dimuon_mass < 10.) *(Jet_Muon_minPhi<0.8) * (MET_pt < 50.) " 

# This class is responsible for making the histogram and plotting it for a given variable
class Variable:
    def __init__(self, varexp, name, nbins, xmin, xmax, jerUp=None, jerDown=None, jesUp=None, jesDown=None, sysUnc=True, logy=True):
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
	self.jerUp = jerUp
	self.jerDown = jerDown 
	self.jesUp = jesUp 
	self.jesDown = jesDown
	self.sysUnc = True
    def Add(self, hist, title, isSignal=False, isData=False):
        hist.SetDirectory(0)
        if isSignal:
            self.signals.append(hist)
            hist.SetLineStyle(len(self.signals))
            self.leg.AddEntry(hist, title, "l")
        elif isData:
            self.data = hist
            self.leg.AddEntry(hist, title, "p")
        else:
            self.stack.Add(hist)
            self.sumMC.Add(hist)
            self.leg.AddEntry(hist, title, "f")
	    #if sysUnc : 
	    
    def Draw(self, suffix, opt):
        print ("plotting "+self.varexp)
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

        self.stack.Draw(opt)
        self.stack.SetMinimum(1)
        if self.logy:
            self.stack.SetMaximum(self.stack.GetMaximum()*1000)
            upperPad.SetLogy()
        else:
            self.stack.SetMaximum(self.stack.GetMaximum()*1.6)
        for signal in self.signals:
            signal.SetFillStyle(0)
            signal.Draw("HIST SAME")
        self.data.Draw("P SAME")
        self.hist_ratio = self.data.Clone("ratio histogram")
        self.hist_ratio.Divide(self.sumMC) 
        lowerPad.cd()
        self.hist_ratio.SetMinimum(0.5)
        self.hist_ratio.SetMaximum(1.5)
        self.hist_ratio.GetXaxis().SetTitle(self.name)
        self.hist_ratio.GetXaxis().SetTitleOffset(2.5)
        self.hist_ratio.Draw("P")

        line = ROOT.TLine(self.xmin, 1, self.xmax, 1)
        line.Draw("SAME")
        canvas.cd()
        self.leg.Draw("SAME")
        makeCMSText(0.13, 0.97,additionalText="Simulation Preliminary")
        makeText(0.13, 0.78, 0.4, 0.78, "jets pt > 30. GeV")
        makeLumiText(0.7, 0.97)
        canvas.SaveAs("plots/"+suffix+self.varexp+".pdf")

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
            self.rdf = self.rdf.Define("weightLumi", "genweight*puweight*tightMuons_weight_iso_nominal*tightMuons_weight_id_nominal*%s*1000.0*%s/%s" %(lumi, find_xsec(path, xsecs), self.sum_weight))
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
	#self.histTotalSystUp = ROOT.TH1F(args[name], args.name, args.nbins, args.xmin, args.xmax)
	self.histTotalSystUp_0 = ROOT.TH1F(args[0], args[1], args[2],args[3],args[4])
	self.histTotalSystUp = ROOT.TH1F(args[0], args[1], args[2],args[3],args[4])
        for i, rdf in enumerate(self.rdfs):
	    var = varexp.split("[")[0]
            if i == 0:
		hist = rdf.Define("var",varexp).Histo1D(args, "var", weight)
                #hist = rdf.Histo1D(args, varexp, weight)
		if sysUnc : 
			histjesUp = rdf.Histo1D(args, jesUp , weight)
			histjerUp = rdf.Histo1D(args, jerUp, weight)
			for i in hist.GetNbinsX() : 
				alpha = hist.GetBinContent(i) 
				beta = histjesUp.GetBinContent(i) 
				gamma = histjerUp.GetBinContent(i) 
				epsilon1 = beta - alpha 
				epsilon2 = gamma - alpha 
				systTotalUp = math.sqrt(epsilon1*epsilon1 + epsilon2*epsilon2)	
				histTotalSystUp_0.SetBinContent(i , alpha + systTotalUp)
	    		histTotalSystUp = histTotalSysUp_0.Clone(histTotalSystUp_0)
            else:
                #tmp_hist = rdf.Histo1D(args, varexp, weight)
		tmp_hist = rdf.Define("var" , varexp).Histo1D(args, "var", weight)
                hist.Add(tmp_hist.GetValue())
		
		

        hist.GetXaxis().SetTitle(args[1])
        hist.SetLineColor(ROOT.TColor.GetColor(colorscale(self.color, 0.8)))
        hist.SetFillColor(ROOT.TColor.GetColor(colorscale(self.color, 1)))
        hist.SetLineWidth(2)
        self.hists.append(hist.Clone())
#	print "histo1D type is : ", type(self.hists[-1])
        return self.hists[-1]

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
            else:
                isData = False
	    print "suffix is :", suffix
	    #Variable.Add should take the name of the uncertainties variables. 
	    #histo1D should return 3 histograms. 
	    variable.Add(process.Histo1D(variable.args, variable.varexp,variable.jesUp, variable.jerUp, variable.jesDown, variable.jerDown,variable.sysUnc, suffix), process.title, isSignal=isSignal, isData=isData)
            #variable.Add(process.Histo1D(variable.args, variable.varexp, suffix), process.title, isSignal=isSignal, isData=isData)
        variable.Draw(suffix, "hist")

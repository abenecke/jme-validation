#! /usr/bin/env python
from utils import *
import ctypes
import re

class CompareSamples():
    def __init__(self, samples, histsname, refsample=None, year='2022', outputPath = None, pdfextraname='', tausel = ['noJetSel']):
        self.samples = samples
        self.tausel = tausel
        self.histsname = histsname
        self.refsample = refsample
        if self.refsample is None:
            self.refsample = list(self.samples.keys())[0]
        self.year = year
        if outputPath is None:
            from ModuleRunner import GenericPath
            outputPath = GenericPath().jme_path
        self.outputPath = os.path.join(outputPath,'pdfs/CompareSamples')
        os.system('mkdir -p '+os.path.join(self.outputPath))
        self.pdfextraname = pdfextraname
        TDR.extraText   = 'Simulation'
        TDR.extraText2  = 'Preliminary'
        TDR.cms_lumi = TDR.commonScheme['legend'][self.year]
        TDR.cms_energy = TDR.commonScheme['energy'][self.year]
        self.infos = {
            'color': {
                'NOMINAL':      rt.kBlack,
                'Tau4GeV':  rt.kOrange+1,
                'Tau':  rt.kOrange+1,
                'HPSTau':  rt.kOrange+4,
                'Tau10GeV': rt.kGreen+2,
                'FromPV2Tau0GeV': rt.kGreen+2,
                'FromPV2Tau4GeV': rt.kRed+1,
                'noCandRemoval': rt.kBlue+1,
                 'CHS': rt.kBlue+1,
                'FromPV2Tau4GeVHPSTauvsHPSTau': rt.kOrange+4,
                'FromPV2Tau4GeVHPSTauvsPUPPI' :  rt.kRed,
                'NOMINALHPSTauvsCHS' : rt.kBlue+1,
                'FromPV2Tau4GeVHPSTauPNETvsPUPPI' : rt.kMagenta,
                'FromPV2Tau4GeVRecoveredTauvsPUPPI' : rt.kGreen+2,
                'FromPV2Tau4GeVPUPPIvsHPSTau' : rt.kRed+2,
                'FromPV2Tau4GeVHPSTauDeepTauvsPUPPI': rt.kYellow+2,
                },
            'color_default': [rt.kBlack, rt.kGreen+2, rt.kRed+1, rt.kOrange+1]
        }
        self.eta_bins = ['0p0to1p3','1p3to2p4','2p4to2p7','2p7to3p0','3p0to5p2']
    
    def LoadInputs(self):
        self.files = OrderedDict()
        self.hists = OrderedDict()
        for sname, fname in self.samples.items():
            # print("sample",sname)
            self.files[sname] = rt.TFile(fname)
            for hname in self.histsname:
                # print("histname",hname)
                h = self.files[sname].Get(hname)
                h.SetDirectory(0)
                h.Rebin(2)
                h.Scale(1./h.Integral())
                self.hists[sname+hname] = h
            if not any([x in self.pdfextraname for x in ['QCD', 'DY']]): continue
            # print("pdfext",self.pdfextraname)
            fname = '/'.join(fname.split('/')[:-2]+['pdfs/'+self.pdfextraname[1:]+'/resp_eff_pur_plots.root'])
            self.files[sname+'resp_eff_pur'] = rt.TFile(fname)

            for eta_bin in self.eta_bins:
                # print("eta_bin",eta_bin)
                h = self.files[sname+'resp_eff_pur'].Get(f'effPurity_eff_pteta_eta{eta_bin}')
                h.SetDirectory(0)
                self.hists[sname+'eff'+eta_bin] = h
                
                h = self.files[sname+'resp_eff_pur'].Get(f'effPurity_purity_pteta_eta{eta_bin}')
                h.SetDirectory(0)
                self.hists[sname+'purity'+eta_bin] = h
                
                h = self.files[sname+'resp_eff_pur'].Get(f'dijet_AK4responseresponse{eta_bin}jer')
                self.hists[sname+'jer'+eta_bin] = h
                
                h = self.files[sname+'resp_eff_pur'].Get(f'dijet_AK4responseresponse{eta_bin}jes')
                self.hists[sname+'jes'+eta_bin] = h

                if "2p7to3p0" in eta_bin or "3p0to5p2" in eta_bin: continue
                if "Tau" not in self.pdfextraname: continue
                for sel_tag in self.tausel:
                    # print("sel_tag",sel_tag)
                    list_histos = [el.GetName() for el in self.files[sname+'resp_eff_pur'].GetListOfKeys() if sel_tag in el.GetName() and eta_bin in el.GetName()]
                    
                    for el in list_histos:
                        if "status" in el and "status" not in sel_tag: continue

                        h = self.files[sname+'resp_eff_pur'].Get(el)
                        h.SetDirectory(0)
                        # print(sname+el)
                        self.hists[sname+el] = h


                    # h = self.files[sname+'resp_eff_pur'].Get(f'{sel_tag}_taueff_leadingtau0p2_TAU_pteta_eta{eta_bin}')
                    # h.SetDirectory(0)
                    # self.hists[sname+f'{sel_tag}tau_leadingtau0p2_{eta_bin}'] = h

                    # if "NOMINAL" in sname:
                    #     h = self.files[sname+'resp_eff_pur'].Get(f'{sel_tag}CHS_taueff_leadingtau0p2_TAU_pteta_eta{eta_bin}')
                    #     h.SetDirectory(0)
                    #     self.hists[f'CHS{sel_tag}tau_leadingtau0p2_{eta_bin}'] = h
                    #     h = self.files[sname+'resp_eff_pur'].Get(f'{sel_tag}Tau_taueff_leadingtau0p2_TAU_pteta_eta{eta_bin}')
                    #     h.SetDirectory(0)
                    #     self.hists[f'Tau{sel_tag}tau_leadingtau0p2_{eta_bin}'] = h
                    #     print("HPS taus",f'{sel_tag}HPSTau_taueff_leadingtau0p2_TAU_pteta_eta{eta_bin}')
                    #     h = self.files[sname+'resp_eff_pur'].Get(f'{sel_tag}HPSTau_taueff_leadingtau0p2_TAU_pteta_eta{eta_bin}')
                    #     h.SetDirectory(0)
                    #     self.hists[f'HPSTau{sel_tag}tau_leadingtau0p2_{eta_bin}'] = h
    
                
    def Close(self):
        for f_ in self.files.values():
            f_.Close()
   
    def Plot(self, histname):
        log_list = ['nJets','pt']
        h_ref = self.hists[self.refsample+histname]
        h_ref_ratio = self.hists[self.refsample+histname].Clone('h_ref_ratio')
        x_min = h_ref.GetBinLowEdge(1)
        x_max = h_ref.GetBinLowEdge(h_ref.GetNbinsX())
        y_max = h_ref.GetMaximum() *3 if not any(el in histname for el in log_list) else h_ref.GetMaximum() * 100
        self.canv = tdrDiCanvas(histname, x_min, x_max, 1e-04, y_max, 0.4, 1.6, histname, 'A.U.', f'Var./{self.refsample}')
        # self.canv = tdrDiCanvas(histname, x_min, y_max, 1, 1.6*h_ref.GetMaximum(), 0.5, 1.5, histname, 'A.U.', 'Ratio')
        # self.canv.cd(1).SetLogy(True)
        self.leg  = tdrLeg(0.40,0.90-(len(self.samples)+1)*0.040,0.90,0.90)
        for s_ind, sname in enumerate(self.samples.keys()):
            h = self.hists[sname+histname]
            h.Scale(h_ref.Integral())
            self.canv.cd(1)
            color = self.infos['color'][sname] if sname in self.infos['color'] else self.infos['color_default'][s_ind]
            tdrDraw(h, 'P', mcolor=color, fstyle=0)
            self.leg.AddEntry(h, sname, 'lp')

            self.canv.cd(2)
            self.hists[sname+histname+'ratio'] = h.Clone(sname+histname+'ratio')
            self.hists[sname+histname+'ratio'].Divide(h_ref_ratio)
            tdrDraw(self.hists[sname+histname+'ratio'], 'P', mcolor=color, fstyle=0)
        if any(el in histname for el in log_list):
            self.canv.cd(1)
            rt.gPad.SetLogy()
        self.canv.SaveAs(os.path.join(self.outputPath, f'{histname}{self.pdfextraname}.pdf'))
        self.canv.Close()
    
    def PlotEff(self, eta_bin, quant = 'eff',sel_tag = ''):
        if not any([x in self.pdfextraname for x in ['QCD', 'DY']]): return
        taustatus = '99'
        if "tau_leadingtau0p2_" in quant:
            # self.samples["CHS"]="CHS"            
            # self.samples["Tau"]="Tau"
            # self.samples["HPSTau"]="HPSTau"
            taustatus = 'inclusive' if not re.findall(r'status\d*',quant) else re.findall(r'status\d*',quant)[0].replace("status","")
        

        h_ref = self.hists[self.refsample+quant+eta_bin]
        h_ref_ratio = self.hists[self.refsample+quant+eta_bin].Clone('h_ref_ratio')
        eta_min, eta_max = eta_bin.replace('p','.').split('to')
        TDR.extraText3 = []
        TDR.extraText3.append(f'{eta_min} < |#eta| < {eta_max}')
        if '99' not in taustatus: TDR.extraText3.append(f'Tau.status = {taustatus}')
        Y_min,Y_max = (0.8,1.4) if not 'jer' in quant else (0.,0.3)
        X_min,X_max = (20,3500) 
        ratioy_min, ratioy_max = (0.8, 1.1)
        if "tau" in quant: 
            Y_min,Y_max = (0.5,1.3)
            X_min,X_max = (20,100) 
            ratioy_min,ratioy_max = (0.85,2) 
        self.canv = tdrDiCanvas(eta_bin, X_min,X_max, Y_min,Y_max, ratioy_min, ratioy_max, 'p_{T}^{gen}',quant if not "tau" in quant else "Eff.", f'Var./{self.refsample}')
        # self.canv = tdrDiCanvas(eta_bin, x_min, y_max, 1e-04, 1, 0.5, 1.5, eta_bin, 'A.U.', 'Ratio')
        # self.canv.cd(1).SetLogy(True)
        # if not "tau" in quant:
        #     self.canv.cd(1).SetLogx(True)
        #     self.canv.cd(2).SetLogx(True)
        FixXAsisPartition(self.canv.cd(2), shift=0.77, textsize=0.11, bins=[30,300,3000])
        self.canv.cd(1)
        self.leg  = tdrLeg(0.40,0.90-(len(self.samples)+1)*0.040,0.90,0.90)
        for s_ind, sname in enumerate(self.samples.keys()):
            h = self.hists[sname+quant+eta_bin]
            self.canv.cd(1)
            color = self.infos['color'][sname] if sname in self.infos['color'] else self.infos['color_default'][s_ind]
            tdrDraw(h, 'P', mcolor=color, fstyle=0)
            self.leg.AddEntry(h, sname, 'lp')
            self.canv.cd(2)
            if 'je' in quant:
                self.hists[sname+eta_bin+'ratio'] = h.Clone(sname+eta_bin+'ratio')
                for xbin in range(1,h_ref_ratio.GetN()):
                    x1,y1,x2,y2 = ctypes.c_double(),ctypes.c_double(),ctypes.c_double(),ctypes.c_double()
                    h_ref_ratio.GetPoint(xbin,x1,y1)
                    self.hists[sname+eta_bin+'ratio'].GetPoint(xbin,x2,y2)
                    self.hists[sname+eta_bin+'ratio'].SetPoint(xbin,x1.value,y2.value/y1.value if y1.value else 0.)
            else:
                self.hists[sname+eta_bin+'ratio'] = h.Clone(sname+eta_bin+'ratio')
                self.hists[sname+eta_bin+'ratio'].Divide(h_ref_ratio)

            tdrDraw(self.hists[sname+eta_bin+'ratio'], 'P', mcolor=color, fstyle=0)
        self.canv.SaveAs(os.path.join(self.outputPath, f'{quant}_{eta_bin}{self.pdfextraname}.pdf'))
        self.canv.Close()
        # if "tau_leadingtau0p2_" in quant:
        #     del self.samples['CHS']
        #     del self.samples['Tau']
        #     del self.samples['HPSTau']


    def PlotHistos(self,hists, refname, eta_bin):
        ratiohists = OrderedDict()
        taustatus = '99'
        if "tau_leadingtau0p2_" in refname:
            taustatus = 'inclusive' if not re.findall(r'status\d*',refname) else re.findall(r'status\d*',refname)[0].replace("status","")
        

        h_ref = hists[refname]
        h_ref_ratio = hists[refname].Clone('h_ref_ratio')
        eta_min, eta_max = eta_bin.replace('p','.').split('to')
        TDR.extraText3 = []
        TDR.extraText3.append(f'{eta_min} < |#eta| < {eta_max}')
        if '99' not in taustatus: TDR.extraText3.append(f'Tau.status = {taustatus}')
        Y_min,Y_max = (0.8,1.4) if not 'jer' in refname else (0.,0.3)
        X_min,X_max = (20,3500) 
        ratioy_min, ratioy_max = (0.8, 1.1)
        if "tau" in refname: 
            Y_min,Y_max = (0.7,1.3)
            X_min,X_max = (20,100) 
            ratioy_min,ratioy_max = (0.65,1.1) 
        self.canv = tdrDiCanvas(eta_bin, X_min,X_max, Y_min,Y_max, ratioy_min, ratioy_max, 'p_{T}^{gen}',refname if not "tau" in refname else "Eff.", f'Var./NOMINAL')
        FixXAsisPartition(self.canv.cd(2), shift=0.63, textsize=0.11, bins=[30,300,3000])
        self.canv.cd(1)
        self.leg  = tdrLeg(0.40,0.90-(len(hists)+1)*0.040,0.90,0.90)
        for ind, name in enumerate(hists.keys()):
            h = hists[name]
            self.canv.cd(1)
            colorname = name.replace(f"_taueff_leadingtau0p2_TAU_pteta_eta{eta_bin}","").replace("customSel","")
            # print(colorname, colorname in self.infos['color'])
            color = self.infos['color'][colorname] if colorname in self.infos['color'] else self.infos['color_default'][ind]
            tdrDraw(h, 'P', mcolor=color, fstyle=0)
            self.leg.AddEntry(h, name.replace("customSel","").replace(f"_taueff_leadingtau0p2_TAU_pteta_eta{eta_bin}",''), 'lp')
            self.canv.cd(2)
            ratiohists[name+'ratio'] = h.Clone(name+'ratio')
            ratiohists[name+'ratio'].Divide(h_ref_ratio)

            tdrDraw(ratiohists[name+'ratio'], 'P', mcolor=color, fstyle=0)
        self.canv.SaveAs(os.path.join(self.outputPath, f'customSel_{refname}_{eta_bin}{self.pdfextraname}.pdf'))
        self.canv.Close()

  
    def PlotAll(self):
        self.LoadInputs()
        for hname in self.histsname:
            self.Plot(histname=hname)

        for eta_bin in self.eta_bins:
            if "Tau" not in self.pdfextraname: 
                self.PlotEff(eta_bin=eta_bin, quant='eff')
                self.PlotEff(eta_bin=eta_bin, quant='purity')
                self.PlotEff(eta_bin=eta_bin, quant='jer')
                self.PlotEff(eta_bin=eta_bin, quant='jes')
            else:
                if "2p7to3p0" in eta_bin or "3p0to5p2" in eta_bin: continue
                for sel_tag in self.tausel:
                    if "customSel" not in sel_tag: 
                        self.PlotEff(eta_bin=eta_bin, quant=f'{sel_tag}TAU_taueff_leadingtau0p2_TAU_pteta_eta')
                    else:
                        tmp_dict = OrderedDict()
                        for name in self.hists:
                            if sel_tag not in name: continue
                            if eta_bin not in name: continue
                            if "status" in name: continue
                            if "NOMINAL" in name and "PUPPI" in name: continue
                            if "NOMINAL" not in name and "PUPPI" not in name: continue
                            tmp_dict[name] = self.hists[name]

                        self.PlotHistos(tmp_dict, f"NOMINALcustomSelHPSTauvsHPSTau_taueff_leadingtau0p2_TAU_pteta_eta{eta_bin}", eta_bin)
        self.Close()




def main():

    variable = ['eta','pt']
    types = ['AK4Jets','Jet1','Jet2','AK8Jets','FatJet1','clak8jet_FatJet1']

    hnames = [f'{typ}_{var}' for typ in types for var in variable]
    hnames+=["AK4Jets_nJets","AK8Jets_nJets"]
    hnames+=[f'FatJet1_{var}' for var in ['nConst','tau21','tau32','mSD']]
    hnames+=[f'FatJet1_{var}' for var in ['particleNetWithMass_QCD','particleNetWithMass_TvsQCD','particleNetWithMass_H4qvsQCD','particleNetWithMass_HbbvsQCD','particleNetWithMass_WvsQCD']]
    hnames+=[f'Jet1_{var}' for var in ['btagDeepFlavB','btagDeepFlavCvB','btagDeepFlavCvL','btagDeepFlavQG','btagPNetB','btagPNetCvB','btagPNetCvL','btagPNetQvG','PNetRegPtRawCorr','PNetRegPtRawCorrNeutrino','PNetRegPtRawRes']]
    hnames+=[f'Jet1_{var}' for var in ['nConst','nConst_eta2p0to3p0']]

    selections = ['noSel']
    histsname = [f'{sel}_{hname}' for sel in selections for hname in hnames]
    histsname += ["dijet_taueff_leadingtau0p2_jet1_PNetTauvsJet","dijet_taueff_leadingtau0p2_jet_PNetTauvsJet"]

    
    samples = {}

    # names = ['NOMINAL', 'Tau4GeV','Tau10GeV','FromPV2Tau4GeV','noCandRemoval']
    names = ['NOMINAL','FromPV2Tau4GeV']
    for name in names:
        samples[name] = f'outputs/QCDModule/QCD_2022_G_Summer22EE{name}_Summer22EERun3_V0_MC_Summer22EERun3_RunF_V0_DATA/results/QCD_Flat.root'
    
    pdfextraname = '_QCD'

    # MP = CompareSamples(samples=samples, histsname=histsname, pdfextraname=pdfextraname, outputPath = 'outputs/QCDModule/',refsample = 'NOMINAL').PlotAll()


    ###### DYTau
    histsname =["noJetSel_taueff_leadingtau0p2_matchedjet_PNetTauvsJet_0","noJetSel_taueff_leadingtau0p2_jet1_PNetTauvsJet","noJetSel_taueff_leadingtau0p2_jet_PNetTauvsJet"]
    # names+=["FromPV2Tau0GeV"]
    for name in names:
        samples[name] = f'outputs/DYModule/DY_2022_G_Summer22EE{name}_Summer22EERun3_V0_MC_Summer22EERun3_RunF_V0_DATA/results/DYTau.root'
    
    pdfextraname = '_DYTau'

    # MP = CompareSamples(samples=samples, histsname=histsname, pdfextraname=pdfextraname, outputPath = 'outputs/DYModule/',refsample = 'NOMINAL',tausel = ["noJetSel","noJetSeltaustatus0","noJetSeltaustatus1","noJetSeltaustatus2","noJetSeltaustatus10","noJetSeltaustatus11","noJetSeltaustatus15"]).PlotAll()

    MP = CompareSamples(samples=samples, histsname=histsname, pdfextraname=pdfextraname, outputPath = 'outputs/DYModule/',refsample = 'NOMINAL',tausel = ["customSel"]).PlotAll() #,"customSeltaustatus0","customSeltaustatus1","customSeltaustatus2","customSeltaustatus10","customSeltaustatus11","customSeltaustatus15"]).PlotAll()

    


if __name__ == '__main__':
    main()

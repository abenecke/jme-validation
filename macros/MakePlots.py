#! /usr/bin/env python

import os
import re
from math import sqrt
from collections import OrderedDict
from tdrstyle_JME import *
import tdrstyle_JME as TDR
rt.gROOT.SetBatch(rt.kTRUE)
rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptFit(0)


def ominus(a,b):
    return sqrt(max(a*a-b*b,0.0))

def oplus(a,b):
    return sqrt(a*a+b*b)

def getMedianError(hist):
    if hist.GetEffectiveEntries()!=0:
        return 1.253 * hist.GetRMS() / sqrt(hist.GetEffectiveEntries())
    else: return 0

def MakeRatioHistograms(h1, h2, name):
    hratio = h1.Clone(name)
    for bin in range(1,h1.GetNbinsX()+1):
        r, dr = (0,0)
        num = h1.GetBinContent(bin)
        den = h2.GetBinContent(bin) 
        if den!=0 and num!=0:
            r = num/den
            dr = r*oplus(h1.GetBinError(bin)/num,h2.GetBinError(bin)/den)
        hratio.SetBinContent(bin, r)
        hratio.SetBinError(bin, dr)
    return hratio

def Confidence(h1, median, confLevel = 0.95):
    ix = h1.GetXaxis().FindBin(median)
    ixlow = ix
    ixhigh = ix
    nb = h1.GetNbinsX()
    ntot = h1.Integral()
    nsum = h1.GetBinContent(ix)
    width = h1.GetBinWidth(ix)
    if ntot==0: return 0
    while (nsum < confLevel * ntot):
        nlow = h1.GetBinContent(ixlow-1) if ixlow>0 else 0
        nhigh = h1.GetBinContent(ixhigh+1) if ixhigh<nb else 0
        if (nsum+max(nlow,nhigh) < confLevel * ntot):
            if (nlow>=nhigh and ixlow>0):
                nsum += nlow
                ixlow -=1
                width += h1.GetBinWidth(ixlow)
            elif ixhigh<nb:
                nsum += nhigh
                ixhigh+=1
                width += h1.GetBinWidth(ixhigh)
            else: raise ValueError('BOOM')
        else:
            if (nlow>nhigh):
                width += h1.GetBinWidth(ixlow-1) * (confLevel * ntot - nsum) / nlow
            else:
                width += h1.GetBinWidth(ixhigh+1) * (confLevel * ntot - nsum) / nhigh
            nsum = ntot
    return width


class MakePlots():
    def __init__(self, year, path, fname, pdfextraname='', campaign = '',jec_algo=('','')):
        # self.eta_bins = ['0p0to5p2', '0p0to1p3','1p3to2p4','2p4to2p7','2p7to3p0','3p0to5p2']
        self.eta_bins = ['0p0to1p3','1p3to2p4','2p4to2p7','2p7to3p0','3p0to5p2']
        # self.eta_bins = ['0p0to1p3']
        self.pt_bins = ['15to17','17to20','20to23','23to27','27to30','30to35','35to40','40to45','45to57','57to72','72to90','90to120','120to150','150to200','200to300','300to400','400to550','550to750','750to1000','1000to1500','1500to2000','2000to2500','2500to3000','3000to3500','3500to4000','4000to4500','4500to5000']
        # self.pt_bins = ['150to200']

        self.types = {
            'effPurity': [('allrecojets',  ['reco','gen','unmatchedgen','unmatchedreco' ]), ('puritymatched', ['reco']), ('effmatched', ['gen'])]
        }
        # self.types = {}
        self.response_names = ['dijet','noLepton','noSel','Zmasscut']
        self.rawresponse_names = ['dijet','Zmasscut']
        
        self.year = year
        self.campaign = campaign
        self.fname = fname
        self.inputPath = path
        self.outputPath = os.path.join(self.inputPath,'pdfs', fname.replace("results/",""))
        for fold in self.types.keys():
            os.system('mkdir -p '+os.path.join(self.outputPath,fold))
        self.pdfextraname = pdfextraname
        self.jec_algo_AK4,self.jec_algo_AK8 = jec_algo


    
    def LoadInputs(self):
        self.files = OrderedDict()
        self.hists = OrderedDict()
        self.graphs = OrderedDict()
        self.files[self.fname] = rt.TFile(os.path.join(self.inputPath,f'{self.fname}.root'))
        f_ = self.files[self.fname]
        # remove response names, that do not exist
        lista = [ el.GetName().split("_")[0] for el in f_.GetListOfKeys() if "response" in el.GetName() and "rawresponse" not in el.GetName()]
        rawInFile = [ el.GetName().split("_")[0] for el in f_.GetListOfKeys() if "rawresponse" in el.GetName()]

        self.response_names = [el for el in self.response_names if el in lista]
        self.rawresponse_names = [el for el in self.rawresponse_names if el in rawInFile]

        self.quant   = array('d',[0.5])
        self.quant_y   = array('d',[0.5])
        fAK4 = open("dummy.txt",'w')
        fAK8 = open("dummy.txt",'w')
        if len(self.rawresponse_names):
            fdirectory = "textFiles/"+self.campaign + self.year+"_USER_"+"MC/"
            os.system('mkdir -p '+fdirectory)
            print("Creating MCtruth form the raw response in folder: ",fdirectory)
            fAK4 = open(fdirectory+self.campaign + self.year+"_USER_"+"MC_L2Relative_AK4PFPuppi.txt","w")
            fAK4.write("{2 JetEta JetPt 1 Rho [0] Correction L2Relative}\n")
            fAK8 = open(fdirectory+self.campaign + self.year+"_USER_"+"MC_L2Relative_AK8PFPuppi.txt","w")
            fAK8.write("{2 JetEta JetPt 1 Rho [0] Correction L2Relative}\n")


        for response_name in self.response_names:
            self.GetResponse(f_, fAK4, response_name+"_AK4response")
            self.GetResponse(f_, fAK8, response_name+"_AK8response")

        for response_name in self.rawresponse_names:
            if(len(self.rawresponse_names)>1):
                print("You can not put multiple corrections in one file")
                break
            else:
                self.GetResponse(f_, fAK4, response_name+"_AK4rawresponse", response_name = "rawresponse")
                self.GetResponse(f_, fAK8, response_name+"_AK8rawresponse", response_name = "rawresponse")
                fAK4.close()
                fAK8.close()

        self.GetEffPurity(f_,"eff","effPurity")
        self.GetEffPurity(f_,"purity","effPurity")

        #read in tau histos and create eff
        #2D hist of pT and eta
        self.tauefflist = list(set([ re.sub(r'_\d','',el.GetName().replace("_denum","").replace("_num","").replace("_pteta","").replace("_TAU","")) for el in f_.GetListOfKeys() if "TAU" in el.GetName()]))
        for el in self.tauefflist:
            self.GetEffPurity(f_,"TAU",selection_name=el, quant = "pteta")

                
    def Close(self):
        self.files['save'] = rt.TFile(os.path.join(self.outputPath, 'resp_eff_pur_plots.root'), 'Recreate')
        for name, graph in self.graphs.items():
            graph.Write(name)
        for name, hist in self.hists.items():
            hist.Write(name)

        for f_ in self.files.values():
            f_.Close()


    def CreateCanvas(self, canvName='', zoom=True, nEntries=7):
        if 'canv' in self.__dict__: self.canv.Close()
        XMin, XMax = (15, 7500)
        YMin, YMax = (0.9,1.1) if zoom else (0.0,2.)
        if 'resolution' in canvName.lower(): YMin, YMax = (0.0,0.5)
        if 'tau' in canvName.lower(): 
            YMin, YMax = (0.7,1.1)
            XMin, XMax = (15, 200)
        if 'response' in canvName.lower(): xName, yName = ('Gen jet p_{T} [GeV]', 'Response')
        if 'resolution' in canvName.lower(): xName, yName = ('Gen jet p_{T} [GeV]', 'Resolution')
        elif 'effmatched' in canvName.lower(): xName, yName = ('Gen jet p_{T} [GeV]', 'Efficiency')
        elif 'puritymatched' in canvName.lower(): xName, yName = ('Reco jet p_{T} [GeV]', 'Purity')
        elif 'TauEff' in canvName.lower(): xName, yName = ('Tau p_{T} [GeV]', 'Tau Efficiency')
        else: xName, yName = ('', '')
        TDR.extraText   = 'Simulation'
        TDR.extraText2  = 'Preliminary'
        TDR.cms_lumi = TDR.commonScheme['legend'][self.year]
        TDR.cms_energy = TDR.commonScheme['energy'][self.year]
        TDR.extraText3 = []
        algo = 'AK4' if not 'ak8' in canvName.lower() else 'AK8'
        PU = self.jec_algo_AK4.upper() if not 'ak8' in canvName.lower() else self.jec_algo_AK8.upper()
        TDR.extraText3.append(algo+' '+PU)
        self.canv = tdrCanvas(canvName, XMin, XMax, YMin, YMax, xName, yName, square=kSquare, isExtraSpace=True)
        self.leg  = tdrLeg(0.60,0.90-(nEntries+1)*0.040,0.90,0.90)
        FixXAsisPartition(self.canv)
        self.canv.SetLogx(True if not "tau" in canvName.lower() else False)
        self.lines= {}
        for y in [0]:
            self.lines[y] = rt.TLine(XMin, y, XMax, y)
            tdrDrawLine(self.lines[y], lcolor=rt.kBlack, lstyle=rt.kDashed, lwidth=1)
    
   
    def GetResponse(self, f_,  f, selection_name="dijet",response_name = "response"):

        orig_hname = f'{selection_name}_ratio'
        orig_hist = f_.Get(orig_hname)
        for z_bin in range(1, orig_hist.GetNbinsZ()+1):
            pts, jes, jer, pts_err, jes_err, jer_err = ([],[],[],[],[], [])
            eta_bin = self.eta_bins[z_bin-1]
            for y_bin in range(1,orig_hist.GetNbinsY()+1):
                pt = self.pt_bins[y_bin-1]
                hname = f'{selection_name}_{response_name}_eta{eta_bin}_pt{pt}'
                # print("hname",hname)

                hist = rt.TH1D(hname, hname, orig_hist.GetNbinsX(), orig_hist.GetXaxis().GetXmin(), orig_hist.GetXaxis().GetXmax())

                # Calculate the product for each X bin
                for x_bin in range(1, orig_hist.GetNbinsX() + 1):
                    value_y = orig_hist.GetBinContent(x_bin, y_bin, z_bin)
                    error_y = orig_hist.GetBinError(x_bin, y_bin, z_bin)

                    # Add the product to the corresponding bin in the current TH1 histogram
                    hist.SetBinContent(x_bin, value_y)
                    hist.SetBinError(x_bin, error_y)


                hist.GetQuantiles(1,self.quant_y,self.quant)
                pt_min, pt_max = pt.split('to')
                pt_min, pt_max = int(pt_min), int(pt_max)
                pts.append((pt_max+pt_min)/2)
                pts_err.append((pt_max-pt_min)/2)
                jes.append(self.quant_y[0])
                jer.append(Confidence(hist, hist.GetMean(), confLevel = 0.87)/(2*1.514))
                jes_err.append(getMedianError(hist))
                jer_err.append(hist.GetRMSError())
                if "raw" in response_name: 
                    if pt_min == 15: pt_min = 0
                    eta_min, eta_max = eta_bin.replace("p",".").split("to")
                    line_plus = f"%s    %s    %i    %i    %i   %i  %i    %.2f \n"%(eta_min, eta_max, pt_min, pt_max, 3,0,200, 1./self.quant_y[0])
                    line_minus = f"%s    %s    %i    %i    %i   %i  %i    %.2f \n"%("-"+eta_min,"-" + eta_max, pt_min, pt_max, 3,0,200, 1./self.quant_y[0])
                    f.write(line_plus)
                    f.write(line_minus)
                
            # print("jes",jes)
            # print("jer",jer)
            self.graphs[selection_name+response_name+eta_bin+'jes'] = rt.TGraphErrors(len(pts), array('d',pts), array('d',jes), array('d',pts_err), array('d',jes_err))
            self.graphs[selection_name+response_name+eta_bin+'jer'] = rt.TGraphErrors(len(pts), array('d',pts), array('d',jer), array('d',pts_err), array('d',jer_err))


    def GetEffPurity(self, f_, effPurity = "eff", selection_name="dijet",quant="pteta"):
        den_hname = f'{selection_name}_{effPurity}_{quant}_denum'
        den_2Dhist = f_.Get(den_hname)
        num_hname = f'{selection_name}_{effPurity}_{quant}_num'
        num_2Dhist = f_.Get(num_hname)

        for y_bin in range(1, num_2Dhist.GetNbinsY()+1):
            eta_bin = self.eta_bins[y_bin-1]
            if "TAU" in effPurity and ("2p7to3p0" in eta_bin or "3p0to5p2" in eta_bin): continue
            hname = f'{selection_name}_{effPurity}_{quant}_eta{eta_bin}'
            # print("hname",hname)

            edges = array("d")
            for x_bin in range(1, num_2Dhist.GetNbinsX() + 1):
                edges.append(num_2Dhist.GetXaxis().GetBinLowEdge(x_bin))
            edges.append(num_2Dhist.GetXaxis().GetXmax())
            # print("edges",edges)

            num_hist = rt.TH1D(hname+"num", hname+"num", len(edges)-1, edges)
            den_hist = rt.TH1D(hname+"den", hname+"den", len(edges)-1, edges)
            self.hists[hname] = num_hist.Clone(hname)

            # Calculate the product for each X bin
            for x_bin in range(1, num_2Dhist.GetNbinsX() + 1):
                # print("eta",eta_bin,"x_bin",x_bin,"num",num_2Dhist.GetBinContent(x_bin,y_bin), "den", den_2Dhist.GetBinContent(x_bin,y_bin))
                # print("num_2Dhist.GetXaxis",num_2Dhist.GetXaxis().GetBinLowEdge(x_bin),"num_hist.GetXaxis",num_hist.GetXaxis().GetBinLowEdge(x_bin))
                num_hist.SetBinContent(x_bin, num_2Dhist.GetBinContent(x_bin,y_bin))
                den_hist.SetBinContent(x_bin, den_2Dhist.GetBinContent(x_bin,y_bin))

                num_hist.SetBinError(x_bin, num_2Dhist.GetBinError(x_bin,y_bin))
                den_hist.SetBinError(x_bin, den_2Dhist.GetBinError(x_bin,y_bin))

            
            self.hists[hname].Divide(num_hist, den_hist, 1, 1, 'B')
        


    def PlotResponse(self, selection_name = "dijet", response_name ='response'):
        infos = {
            '0p0to1p3' : {'legend': '0.0 < |#eta| < 1.3', 'color':  rt.kGreen+3,  'marker': rt.kFullSquare,       'msize':  1.0},
            '1p3to2p4' : {'legend': '1.3 < |#eta| < 2.4', 'color':  rt.kBlue-4,   'marker': rt.kFullStar,         'msize':  1.3},
            '2p4to2p7' : {'legend': '2.4 < |#eta| < 2.7', 'color':  rt.kAzure+2,  'marker': rt.kFullCrossX,       'msize':  1.3},
            '2p7to3p0' : {'legend': '2.7 < |#eta| < 3.0', 'color':  rt.kRed+1,    'marker': rt.kFullTriangleUp,   'msize':  1.0},
            '3p0to5p2' : {'legend': '3.0 < |#eta| < 5.2', 'color':  rt.kOrange+1, 'marker': rt.kFullTriangleDown, 'msize':  1.0},
            '0p0to5p2' : {'legend': '0.0 < |#eta| < 5.2', 'color':  rt.kBlack,    'marker': rt.kFullCircle,       'msize':  1.0},
        }

        self.CreateCanvas(canvName=response_name+"_"+selection_name, zoom=True if not 'raw' in response_name else False)
        for name, graph in self.graphs.items():
            if not 'jes' in name: continue
            if not selection_name in name: continue
            if ("raw" in response_name and "raw" not in name) or ("raw" not in response_name and "raw" in name): continue
            info = infos[name.replace('jes','').replace(selection_name,'').replace(response_name,'').replace("_","")]
            tdrDraw(graph, 'P', msize=info['msize'], marker=info['marker'], mcolor=info['color'])
            self.leg.AddEntry(graph, info['legend'], 'lp')
        self.canv.SaveAs(os.path.join(self.outputPath, response_name+'_'+selection_name + self.pdfextraname + '.pdf'))

        self.CreateCanvas(canvName='resolution'+"_"+selection_name, zoom=True)
        for name, graph in self.graphs.items():
            if not 'jer' in name: continue
            if not selection_name in name: continue
            if ("raw" in response_name and "raw" not in name) or ("raw" not in response_name and "raw" in name): continue
            info = infos[name.replace('jer','').replace(selection_name,'').replace(response_name,'').replace('_','')]
            tdrDraw(graph, 'P', msize=info['msize'], marker=info['marker'], mcolor=info['color'])
            self.leg.AddEntry(graph, info['legend'], 'lp')
        
        self.canv.SaveAs(os.path.join(self.outputPath,('raw' if 'raw' in response_name else '')+ 'resolution'+'_'+selection_name + self.pdfextraname + '.pdf'))

    def PlotEffPurity(self, effPurity="eff",selection_name="effPurity",quant="pteta"):
        for eta_bin in self.eta_bins:
            if "TAU" in effPurity and ("2p7to3p0" in eta_bin or "3p0to5p2" in eta_bin): continue
            hname = f'{selection_name}_{effPurity}_{quant}_eta{eta_bin}'
            self.CreateCanvas(canvName=hname, zoom=False)
            hist = self.hists[hname]
            tdrDraw(hist, 'P', fstyle=0)
            # tdrDraw(h, opt, marker=rt.kFullCircle, msize=1.0, mcolor=rt.kBlack, lstyle=rt.kSolid, lcolor=-1, fstyle=1001, fcolor=rt.kYellow+1, alpha=-1):
            self.canv.SaveAs(os.path.join(self.outputPath,  hname + self.pdfextraname + '.pdf'))
    
    def PlotAll(self):
        self.LoadInputs()

        for response_name in self.response_names:
            self.PlotResponse(selection_name=response_name+"_AK4response", response_name = "response")
            self.PlotResponse(selection_name=response_name+"_AK8response", response_name = "response")
        for response_name in self.rawresponse_names:
            self.PlotResponse(selection_name=response_name+"_AK4rawresponse", response_name = "rawresponse")
            self.PlotResponse(selection_name=response_name+"_AK8rawresponse", response_name = "rawresponse")

        self.PlotEffPurity()
        self.PlotEffPurity(effPurity="purity")
        for el in self.tauefflist:
            self.PlotEffPurity(effPurity="TAU",selection_name = el)
        self.Close()



def main():
    MP = MakePlots().PlotAll()

if __name__ == '__main__':
    main()

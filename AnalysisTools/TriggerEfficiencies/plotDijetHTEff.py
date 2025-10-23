#!/usr/bin/env python3
'''
DESCRIPTION:


'''
#================
# Import modules
#================
from argparse import ArgumentParser
import ROOT
import math 
import time

def getCanvas():
    d = ROOT.TCanvas("", "", 800, 700)
    d.SetLeftMargin(0.12)
    d.SetRightMargin(0.15)
    d.SetLeftMargin(0.13)
    return d

def AddPrivateWorkText(setx=0.21, sety=0.905):
    tex = ROOT.TLatex(0.,0.,'Private Work'); #'HLT tutorial');
    tex.SetNDC();
    tex.SetX(setx);
    tex.SetY(sety);
    tex.SetTextFont(53);
    tex.SetTextSize(28);
    tex.SetLineWidth(2)
    return tex

def AddCMSText(setx=0.205, sety=0.905):
    texcms = ROOT.TLatex(0.,0., 'CMS');
    texcms.SetNDC();
    texcms.SetTextAlign(31);
    texcms.SetX(setx);
    texcms.SetY(sety);
    texcms.SetTextFont(63);
    texcms.SetLineWidth(2);
    texcms.SetTextSize(30);
    return texcms

def createLegend():
    legend = ROOT.TLegend(0.58, 0.20, 0.82, 0.35)
    legend.SetFillColor(0)
    legend.SetFillStyle(0);
    legend.SetBorderSize(0);
    return legend

def SetStyle(h, color, marker_style):
    h.SetLineColor(color)
    h.SetMarkerColor(color)
    h.SetMarkerStyle(marker_style)
    return h

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTextFont(42)

def SetStyle(h, COLOR):
    h.SetMarkerStyle(21)
    h.SetMarkerColor(COLOR)
    h.SetLineColor(COLOR)
    return h

colors = {0: ROOT.kBlack,
          1: ROOT.kBlue,
          2: ROOT.kGreen+1,
          3: ROOT.kRed+1,
          4: ROOT.kOrange-3,
          5: ROOT.kMagenta+2,
          6: ROOT.kTeal+3,
          }
          

def main(args):

    f = ROOT.TFile(args.rfile, "READ")
    # Post Processor histDirName
    fdir = f.GetDirectory("DijethtTrigAnalyzerNanoAOD")
    statOption = ROOT.TEfficiency.kFCP

    variables  = ["minv_1", "minv_2", "minv_3", "minv_4"] #, "pv"] After pass ref trigger
    triggers   = ["passed"] #HLT or DST (only the las tone now)

    for var in variables:
        c = getCanvas()
        leg = createLegend()
        # Pass Ref Trigger & Offline Selection
        den = fdir.Get(f'h_{var}_all')

        nums = {}
        effs = {}
        for j, trg in enumerate(triggers):
            # Pass Ref Trigger, Signal Trigger & Offline Selection
            nums[trg] = fdir.Get(f'h_{var}_{trg}')
            effs[trg] = ROOT.TEfficiency(nums[trg], den)
            effs[trg].SetStatisticOption(statOption)
            effs[trg] = SetStyle(effs[trg], colors[j])
            if j == 0:
                effs[trg].Draw()
            else:
                effs[trg].Draw("same")
            #leg.AddEntry(effs[trg], "OR)
            leg.AddEntry(effs[trg], trg.replace("passed", "DST_PFScouting_JetHT_v"))
                
        c.Modified()
        c.Update()
        #effs["passedL1"].GetPaintedGraph().GetYaxis().SetRangeUser(0.0, 1.2)
        leg.Draw("same")
        
        # Styling stuff
        tex_cms = AddCMSText()
        tex_cms.Draw("same")
        
        private = AddPrivateWorkText()
        private.Draw("same")
        
        header = ROOT.TLatex()
        header.SetTextSize(0.04)
        header.DrawLatexNDC(0.56, 0.905, "2024I, #sqrt{s} = 13.6 TeV")
        #header.DrawLatexNDC(0.57, 0.905, "2024I, #sqrt{s} = 13 TeV")
        
        c.Update()
        c.Modified()
        for fs in args.formats:
            savename = f'plots/TrgEffs_DijetHT_Minv_AN_Ortogonal_{var}{fs}'
            c.SaveAs(savename)
            
    '''
    # Plot 2D:
    c2D = getCanvas()

    den = fdir.Get('h_pfht_vs_pv_all')
    num = fdir.Get('h_pfht_vs_pv_passed')
    
    eff2D = ROOT.TEfficiency(num, den)
    eff2D.Draw("COLZ")
    c2D.Modified()
    c2D.Update()
    tex_cms = AddCMSText()
    tex_cms.Draw("same")
    private = AddPrivateWorkText()
    private.Draw("same")
    header = ROOT.TLatex()
    header.SetTextSize(0.04)
    header.DrawLatexNDC(0.57, 0.905, "2023, #sqrt{s} = 13.6 TeV")
    c2D.Update()
    c2D.Modified()
    for fs in args.formats:
        c2D.SaveAs("plots/Eff2D_PFHTtrigger_PFHTvsPV%s" % (fs))
    '''
if __name__ == "__main__":

    VERBOSE       = True
    YEAR          = "2024"
    # Post Processor histFileName
    TRGROOTFILE   = "histos_DijetHTTrigNanoAOD.root"
    FORMATS       = ['.png', '.pdf']

    parser = ArgumentParser(description="Derive the trigger scale factors")
    parser.add_argument("-v", "--verbose", dest="verbose", default=VERBOSE, action="store_true", help="Verbose mode for debugging purposes [default: %s]" % (VERBOSE))
    parser.add_argument("--rfile", dest="rfile", type=str, action="store", default=TRGROOTFILE, help="ROOT file containing the denominators and numerators [default: %s]" % (TRGROOTFILE))
    parser.add_argument("--year", dest="year", action="store", default=YEAR, help="Process year")
    parser.add_argument("--formats", dest="formats", default=FORMATS, action="store", help="Formats to save histograms")

    args = parser.parse_args()
    main(args)

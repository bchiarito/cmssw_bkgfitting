from __future__ import print_function
import math
import ROOT
import sys
import os
import argparse
import array
import fitting_utils as util
import constants as VALS

# command line options
parser = argparse.ArgumentParser(description="")
parser.add_argument("input", metavar="INPUT", help="input root file")
parser.add_argument("--testBin", default=None, help="specify bin to test")
parser.add_argument("--testRegion", default=None, choices = ["iso_sym", "iso_asym", "noniso_sym", "noniso_asym"], help="specify region to test")
parser.add_argument("--name", default="plots", help="create name for plots pdf")
parser.add_argument("--useUnscaledTight", default=False, action="store_true", help="used unscaled tight in preference to scaled")
parser.add_argument("--show", default=False, action="store_true", help="")
args = parser.parse_args()

# constants
egamma_rootfile = "summed_egamma.root"

# plotting style
leg_x1, leg_x2, leg_y1, leg_y2 = 0.7, 0.60, 0.89, 0.89
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)
ROOT.gStyle.SetLegendFillColor(ROOT.TColor.GetColorTransparent(ROOT.kWhite, 0.01));
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

# select regions
regions = ["iso_sym", "iso_asym", "noniso_sym", "noniso_asym"]
if args.testRegion: regions = [args.testRegion]
eta_regions = ["barrel", "endcap"]
photon_regions = ["tight", "loose"]
if args.testBin is not None: test_bin = (args.testBin).split(" ")

# init
os.chdir(args.input)
infile1 = ROOT.TFile(egamma_rootfile)
c1 = ROOT.TCanvas("c1", "c1", 800, 600)
c1.Print(args.name + ".pdf[")

# run
loose_fit_files = []
scaled_tight_files = []
tight_temp_files = []
tight_deg_files = []
tight_chi2_files = []
tight_poly_files = []
file_counter_loose = 0
file_counter_tight = 0
file_counter_temp = 0
file_counter_deg = 0
file_counter_chi2 = 0
file_counter_poly = 0
for region in regions:  # loop through twoprong sideband regions
    if args.testBin is not None: 
        if not region == test_bin[0]: continue
    for i in range(len(VALS.PT_EDGES)):  # loop through pt bins for a fixed twoprong sideband
        if args.testBin is not None: 
            if not VALS.PT_EDGES[i] == int(test_bin[2]): continue
        for eta_reg in eta_regions:  # loop through eta regions for fixed pt-bin and fixed twoprong sideband
            if args.testBin is not None: 
                if not eta_reg == test_bin[1]: continue
           
            # Generate correct plots names to access from summed histogram files
            egamma_loose_plots = "plots/twoprong_masspi0_" + region + "_" + eta_reg

            if i == len(VALS.PT_EDGES) - 1: egamma_loose_plots += "_" + str(VALS.PT_EDGES[i]) + "+"
            else: egamma_loose_plots += "_" + str(VALS.PT_EDGES[i]) + "_" + str(VALS.PT_EDGES[i+1])
            if i == len(VALS.PT_EDGES) - 1: hist_name = region + "_" + eta_reg + "_" + str(VALS.PT_EDGES[i]) + "+"
            else: hist_name = region + "_" + eta_reg + "_" + str(VALS.PT_EDGES[i]) + "_" + str(VALS.PT_EDGES[i+1]) 

            # Reference name of the histogram created in the backend 
            egamma_loose_plots += "_loose"

            h_egamma_loose = infile1.Get(egamma_loose_plots)
            h_egamma_loose.SetBinErrorOption(ROOT.TH1.kPoisson)
            h_egamma_loose.SetLineColor(ROOT.kBlack)
            h_egamma_loose.SetTitle("")
            h_egamma_loose.GetXaxis().SetTitle("")
            h_egamma_loose.GetYaxis().SetTitle("")
           
            # Input loose templates
            if not os.path.exists('loose_fit_hists'): 
                print("ERROR: Must create directory for loose templates titled loose_fit_hists")
                exit()
            else:
                loose_fit_files.append(ROOT.TFile("loose_fit_hists/" + hist_name + "_loose.root"))
                loose_fit_as_hist = loose_fit_files[file_counter_loose].Get(hist_name+"_loose")
                loose_fit_as_hist.SetLineColor(ROOT.kBlack)
                loose_fit_as_hist.SetTitle("")
                loose_fit_as_hist.GetXaxis().SetTitle("")
                loose_fit_as_hist.GetYaxis().SetTitle("")
                file_counter_loose += 1

            if args.useUnscaledTight or region == "iso_sym":
                # Get the histograms from the input file
                h_egamma_tight = infile1.Get("plots/twoprong_masspi0_"+hist_name+"_tight")
                h_egamma_tight.SetLineColor(ROOT.kBlack)
                h_egamma_tight.SetStats(0)
                h_egamma_tight.SetTitle("")
                h_egamma_tight.GetXaxis().SetTitle("")
                h_egamma_tight.GetYaxis().SetTitle("")
            else:
                # Input scaled tight data
                if not os.path.exists('scaled_tight_hists'):
                    print("ERROR: Must create directory for scaled data titled scaled_tight_hists")
                    exit()
                else:
                    scaled_tight_files.append(ROOT.TFile("scaled_tight_hists/" + hist_name + "_tight.root"))
                    h_egamma_tight = scaled_tight_files[file_counter_tight].Get(hist_name+"_tight")
                    h_egamma_tight.SetLineColor(ROOT.kBlack)
                    h_egamma_tight.SetStats(0)
                    h_egamma_tight.SetTitle("")
                    h_egamma_tight.GetXaxis().SetTitle("")
                    h_egamma_tight.GetYaxis().SetTitle("")
                    file_counter_tight += 1
            
            # Set Poisson errors for tight histogram
            h_egamma_tight.SetBinErrorOption(ROOT.TH1.kPoisson)

            # Input tight template
            if args.useUnscaledTight:
                if not os.path.exists('tight_templates/templates_noscaling'):
                    print("ERROR: Must create directory for scaled data as follows: tight_templates/templates_noscaling")
                    exit()
                else:
                    tight_temp_files.append(ROOT.TFile("tight_templates/templates_noscaling/" + hist_name + "_tight_temp.root"))
                    tight_fit_as_hist = tight_temp_files[file_counter_temp].Get(hist_name+"_tight_temp")
                    file_counter_temp += 1
            else:
                if not os.path.exists('tight_templates/templates'):
                    print("ERROR: Must create directory for scaled data as follows: tight_templates/templates")
                    exit()
                else:
                    tight_temp_files.append(ROOT.TFile("tight_templates/templates/" + hist_name + "_tight_temp.root"))
                    tight_fit_as_hist = tight_temp_files[file_counter_temp].Get(hist_name+"_tight_temp")
                    file_counter_temp += 1

            # Input tight degrees
            if not os.path.exists('tight_templates/degrees'):
                print("ERROR: Must create directory for scaled data Bernstein degrees as follows: tight_templates/degrees")
                exit()
            else:
                tight_deg_files.append(ROOT.TFile("tight_templates/degrees/" + hist_name + "_tight_temp_deg.root"))
                tight_deg_hist = tight_deg_files[file_counter_deg].Get(hist_name+"_tight_temp_deg")
                file_counter_deg += 1
                bern_deg = int(tight_deg_hist.GetBinContent(1))

            # Input chi2 values
            if not os.path.exists('tight_templates/chi2s'):
                print("ERROR: Must create directory for scaled data Bernstein degrees as follows: tight_templates/chi2s")
                exit()
            else:
                tight_chi2_files.append(ROOT.TFile("tight_templates/chi2s/" + hist_name + "_tight_poly_chi2.root"))
                tight_chi2_hist = tight_chi2_files[file_counter_chi2].Get(hist_name+"_tight_poly_chi2")
                file_counter_chi2 += 1
                chi2_val = round(tight_chi2_hist.GetBinContent(1), 4)

            # Input bern poly TF1 (to extract params)
            if not os.path.exists('tight_templates/polys'):
                print("ERROR: Must create directory for scaled data Bernstein degrees as follows: tight_templates/polys")
                exit()
            else:
                tight_poly_files.append(ROOT.TFile("tight_templates/polys/" + hist_name + "_tight_poly.root"))
                bern_poly = tight_poly_files[file_counter_poly].Get(hist_name+"_tight_poly") 
                file_counter_poly += 1

            if bern_deg == 0: just_poly = ROOT.TF1("bern_polynomial", util.bern_deg0, 0, 25, 1)
            if bern_deg == 1: just_poly = ROOT.TF1("bern_polynomial", util.bern_deg1, 0, 25, 2)
            if bern_deg == 2: just_poly = ROOT.TF1("bern_polynomial", util.bern_deg2, 0, 25, 3)
            if bern_deg == 3: just_poly = ROOT.TF1("bern_polynomial", util.bern_deg3, 0, 25, 4)
            if bern_deg == 4: just_poly = ROOT.TF1("bern_polynomial", util.bern_deg4, 0, 25, 5)
            if bern_deg == 5: just_poly = ROOT.TF1("bern_polynomial", util.bern_deg5, 0, 25, 6)
            if bern_deg == 6: just_poly = ROOT.TF1("bern_polynomial", util.bern_deg6, 0, 25, 7)
            if bern_deg == 7: just_poly = ROOT.TF1("bern_polynomial", util.bern_deg7, 0, 25, 8)
            
            # Set bernstein poly params
            for j in range(bern_deg+1): just_poly.SetParameter(j, bern_poly.GetParameter(j))

            fitted_func = util.HistogramToFunction(loose_fit_as_hist)
            fitted_func_times_constant, _, _ = util.MultiplyWithPolyToTF1(fitted_func, 0, poly=0)
            #fit_result = h_egamma_tight.Fit(fitted_func_times_constant, '0L')
            tight_fit_w_constant = util.TemplateToHistogram(fitted_func_times_constant, 1000, 0, 50)

            #just_poly = util.ExtractPolyFromTightFit(func_with_poly, poly=3)
            tlast_bin = h_egamma_tight.GetNbinsX() + 1
            for b in reversed(range(h_egamma_tight.GetNbinsX())):
              if h_egamma_tight.GetBinContent(b+1)==0: continue
              else:
                tlast_bin = b + 1
                break
            rightmost_tightdata = h_egamma_loose.GetBinLowEdge(tlast_bin+1)
            #just_poly.SetRange(0,rightmost_tightdata)

            h_loose_pull_num = h_egamma_loose.Clone()
            h_loose_pull_num.Reset()
            h_loose_pull = h_egamma_loose.Clone()
            h_loose_pull.Reset()
            h_loose_pull_num.Add(h_egamma_loose, loose_fit_as_hist, 1, -1)  # Numerator of pull hist is data - fit

            for j in range(h_loose_pull_num.GetNbinsX()): 
                if h_egamma_loose.GetBinContent(j+1) == 0: 
                    err = h_egamma_loose.GetBinErrorUp(j+1)
                else: 
                    if loose_fit_as_hist.GetBinContent(j+1) > h_egamma_loose.GetBinContent(j+1):
                        err = h_egamma_loose.GetBinErrorUp(j+1)
                    else:
                        err = h_egamma_loose.GetBinErrorLow(j+1)
                h_loose_pull.SetBinContent(j+1, h_loose_pull_num.GetBinContent(j+1)/err)
                h_loose_pull.SetBinError(j+1, 1)
            
            h_tight_pull_num = h_egamma_tight.Clone()
            h_tight_pull_num.Reset()
            h_tight_pull = h_egamma_tight.Clone()
            h_tight_pull.Reset()
            h_tight_pull_num.Add(h_egamma_tight, tight_fit_as_hist, 1, -1)  # Numerator of pull hist is data - fit
            bin_bin_error = 0
            
            h_tight_pull_error = h_tight_pull.Clone() # to visualize binbin error
            h_tight_pull_error.Reset()

            for j in range(h_tight_pull_num.GetNbinsX()): 
                if h_egamma_tight.GetBinContent(j+1) == 0:
                    err = h_egamma_tight.GetBinErrorUp(j+1)
                else: 
                    if tight_fit_as_hist.GetBinContent(j+1) > h_egamma_tight.GetBinContent(j+1):
                        err = h_egamma_tight.GetBinErrorUp(j+1)
                    else:
                        err = h_egamma_tight.GetBinErrorLow(j+1)
                h_tight_pull.SetBinContent(j+1, h_tight_pull_num.GetBinContent(j+1)/err)
                h_tight_pull.SetBinError(j+1, 1)
                h_tight_pull_error.SetBinContent(j+1, 0)
                h_tight_pull_error.SetBinError(j+1, (tight_fit_as_hist.GetBinContent(j+1)*bin_bin_error)/err)

            # pull for tight fit with constant
            h_tight_pullc = h_egamma_tight.Clone()
            h_tight_pullc.Reset()
            h_tight_pullc_num = h_egamma_tight.Clone()
            h_tight_pullc_num.Reset()
            h_tight_pullc_num.Add(h_egamma_tight, tight_fit_w_constant, 1, -1)  # Numerator of pull hist is data - fit
            for j in range(h_tight_pullc_num.GetNbinsX()): 
                if h_egamma_tight.GetBinContent(j+1) == 0:
                    err = h_egamma_tight.GetBinErrorUp(j+1)
                else: 
                    if tight_fit_w_constant.GetBinContent(j+1) > h_egamma_tight.GetBinContent(j+1):
                        err = h_egamma_tight.GetBinErrorUp(j+1)
                    else:
                        err = h_egamma_tight.GetBinErrorLow(j+1)
                h_tight_pullc.SetBinContent(j+1, h_tight_pullc_num.GetBinContent(j+1)/err)
                h_tight_pullc.SetBinError(j+1, 0) # no error bar
            
            # Legend creation
            legend1 = ROOT.TLegend(0.45, 0.75, 0.88, 0.87)
            legend1.AddEntry(h_egamma_loose, "Loose Photon Sideband", "l")
            legend1.AddEntry(loose_fit_as_hist, "Landau + Exponential Fit", "l")
            legend1.SetTextSize(0.03)

            legend2 = ROOT.TLegend(0.53, 0.72, 0.88, 0.87)
            if not region == "iso_sym":
                legend2.AddEntry(h_egamma_tight, "Full Selection", "l")
                legend2.AddEntry(tight_fit_as_hist, 'Bernstein, Degree ' + str(int(bern_deg)), 'l')
                legend2.AddEntry('', 'Chi2/Ndof: ' + str(chi2_val), '')
                legend2.SetTextSize(0.03)
           
            # Create Title 
            title1 = "p_{T} " + str(VALS.PT_EDGES[i])
            if i == len(VALS.PT_EDGES)-1: title1 += "+ GeV" 
            else: title1 += "-" + str(VALS.PT_EDGES[i+1]) + " GeV"
            if eta_reg == "barrel": eta1 = "Barrel TP"
            else: eta1 = "Endcap TP"
           
            title2 = "TP Control Region: "
            if region == "iso_sym": title2 += "Isolated Symmetric"
            elif region == "iso_asym": title2 += "Isolated Asymmetric"
            elif region == "noniso_sym": title2 += "Nonisolated Symmetric"
            else: title2 += "Nonisolated Asymmetric"

            # Draw plots
            # TOP-LEFT PANEL
            c1.cd()
            pad1 = ROOT.TPad('pad1', 'pad1', 0, 0.3, 0.5, 1)
            pad1.SetLeftMargin(0.15)
            pad1.SetBottomMargin(0.1)
            pad1.Draw()
            pad1.cd()
            # Draw loose hist and fit template
            h_egamma_loose.GetXaxis().SetTitle("")
            h_egamma_loose.GetXaxis().SetLabelSize(0)
            h_egamma_loose.GetYaxis().SetTitleSize(0.05)
            h_egamma_loose.SetTitleOffset(1, "Y")
            h_egamma_loose.GetYaxis().SetTitle("Entries / 50 MeV")
            h_egamma_loose.SetMaximum()
            h_egamma_loose.SetMinimum(0.1)
            h_egamma_loose.Draw("e")
            loose_fit_as_hist.SetLineColor(ROOT.kRed+2)
            loose_fit_as_hist.Draw("same")
            ROOT.gPad.SetLogy()
            # Dynamically set the x-axis range 
            h_egamma_loose.GetXaxis().SetRangeUser(0, util.findMaxXVal(loose_fit_as_hist)*3/5)
            h_egamma_loose.GetYaxis().SetRangeUser(0.1, h_egamma_loose.GetMaximum()*3.5)
            legend1.Draw("same")
            # Draw CMS logo
            l1_cms = ROOT.TLatex(0.165, 0.84, "CMS")
            l1_cms.SetNDC(ROOT.kTRUE)
            l1_cms.Draw("same")
            # Draw lumi at top of plot
            l1_lumi = ROOT.TLatex(0.6, 0.915, "137.6 fb^{-1} (13 TeV)")
            l1_lumi.SetNDC(ROOT.kTRUE)
            l1_lumi.SetTextSize(0.04)
            l1_lumi.SetTextFont(42)
            l1_lumi.Draw("same")
            l1_loose = ROOT.TLatex(0.15, 0.915, "Loose Photon")
            l1_loose.SetNDC(ROOT.kTRUE)
            l1_loose.SetTextSize(0.04)
            l1_loose.SetTextFont(42)
            l1_loose.Draw("same")
            # Specify pt bin in loose plot
            title1_bin = ROOT.TLatex(0.60, 0.71, title1)
            title1_bin.SetTextSize(0.04)
            title1_bin.SetNDC(ROOT.kTRUE)
            title1_bin.Draw("same")
            title1_eta = ROOT.TLatex(0.60, 0.66, eta1)
            title1_eta.SetTextSize(0.04)
            title1_eta.SetNDC(ROOT.kTRUE)
            title1_eta.Draw("same")
            control_reg = ROOT.TLatex(0.05, 0.96, title2)  # 0.05, 0.96
            control_reg.SetTextSize(0.042)
            control_reg.SetNDC(ROOT.kTRUE)
            control_reg.Draw("")
            ROOT.gPad.Update()
            
            # TOP-RIGHT PANEL
            c1.cd()
            pad2 = ROOT.TPad('pad2', 'pad2', 0.47, 0.3, 1, 1)
            pad2.Draw()
            pad2.cd()
            if not region == "iso_sym":
                ROOT.gPad.SetLogy()
                h_egamma_tight.GetXaxis().SetLabelSize(0)
                h_egamma_tight.GetYaxis().SetTitleSize(0.05)
                h_egamma_tight.SetTitleOffset(1, "Y")
                h_egamma_tight.GetYaxis().SetTitle("Entries / 50 MeV")
                h_egamma_tight.GetXaxis().SetRangeUser(0, util.findMaxXVal(loose_fit_as_hist)*3/5)
                h_egamma_tight.SetMinimum(0.1)
                h_egamma_tight.Draw("e")
                tight_fit_w_constant.SetLineColor(ROOT.kBlue)
                tight_fit_as_hist.SetLineColor(ROOT.kRed+2)
                tight_fit_as_hist.SetLineWidth(1)
                tight_fit_as_hist.Draw("same hist")
                h_egamma_tight.Draw("e same")
                h_egamma_tight.GetYaxis().SetRangeUser(0.1, h_egamma_tight.GetMaximum()*2.3)
                l2 = ROOT.TLatex(0.12, 0.84, "CMS")
                l2.SetNDC(ROOT.kTRUE)
                l2.Draw("same")
                l2_lumi = ROOT.TLatex(0.615, 0.915, "137.6 fb^{-1} (13 TeV)")
                l2_lumi.SetNDC(ROOT.kTRUE)
                l2_lumi.SetTextSize(0.04)
                l2_lumi.SetTextFont(42)
                l2_lumi.Draw("same")
                legend2.Draw("same")
                l2_tight = ROOT.TLatex(0.105, 0.915, "Tight Photon")
                l2_tight.SetNDC(ROOT.kTRUE)
                l2_tight.SetTextSize(0.04)
                l2_tight.SetTextFont(42)
                l2_tight.Draw("same")
                overlay = ROOT.TPad("overlay","",0, 0.06, 1, 0.5)
                overlay.SetFillStyle(4000)
                overlay.SetFillColor(0)
                overlay.SetFrameFillStyle(4000)
                overlay.SetFrameLineWidth(0)
                overlay.Draw()
                overlay.cd()
                empty = ROOT.TH1F(util.getname('empty'), '', 100, 0, 50)
                empty.SetLineColor(ROOT.kRed)
                empty.GetXaxis().SetRangeUser(0, util.findMaxXVal(loose_fit_as_hist)*3/5)
                empty.GetYaxis().SetRangeUser(min(0, just_poly.GetMinimum()), just_poly.GetMaximum())
                empty.Draw('AH')
                just_poly.SetRange(0, 50)
                just_poly.SetTitle("")
                just_poly.Draw("AI L same")
                ROOT.gPad.Update()
                rightaxis = ROOT.TGaxis(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUymax(), 510, "L+")
                rightaxis.SetLineColor(ROOT.kRed);
                rightaxis.SetLabelColor(ROOT.kRed);
                rightaxis.Draw()
                ROOT.gPad.Update()
            
            # BOTTOM-LEFT PANEL
            c1.cd()
            pad3 = ROOT.TPad('pad3', 'pad3', 0, 0, 0.5, 0.33)
            pad3.SetLeftMargin(0.15)
            pad3.SetTopMargin(0.03)
            pad3.SetBottomMargin(0.23)
            pad3.Draw()
            pad3.cd()
            h_loose_pull.SetTitleOffset(0.5, "Y")
            h_loose_pull.GetYaxis().SetTitleSize(0.1)
            h_loose_pull.GetYaxis().SetLabelSize(0.06)
            h_loose_pull.GetYaxis().SetTitle("Pull")
            h_loose_pull.SetLineColor(ROOT.kBlack)
            h_loose_pull.Draw('pe')
            h_loose_pull.GetXaxis().SetLabelSize(0.07)
            h_loose_pull.GetXaxis().SetTitleSize(0.1)
            h_loose_pull.GetXaxis().SetTitle("m_{TP} (GeV)")
            h_loose_pull.SetMarkerStyle(21)
            h_loose_pull.SetMarkerSize(0.25)
            h_loose_pull.SetStats(0)
            h_loose_pull.GetXaxis().SetRangeUser(0, util.findMaxXVal(loose_fit_as_hist)*3/5)
            loose_pull_abs_max = max(abs(h_loose_pull.GetMinimum()), h_loose_pull.GetMaximum()) 
            tight_pull_abs_max = max(abs(h_tight_pull.GetMinimum()), h_tight_pull.GetMaximum())
            abs_pull_max = max(loose_pull_abs_max, tight_pull_abs_max)
            h_loose_pull.GetYaxis().SetRangeUser(-abs_pull_max-2, abs_pull_max+2)
            line3 = ROOT.TLine(0.2, 0.6, 0.9, 0.6)
            line3.SetNDC(ROOT.kTRUE)
            line3.SetLineWidth(1)
            line3.SetLineColorAlpha(ROOT.kBlack, 0.5)
            line3.Draw("same")
            
            # BOTTOM-RIGHT PANEL
            c1.cd()
            pad4 = ROOT.TPad('pad4', 'pad4', 0.47, 0, 1, 0.33)
            pad4.SetTopMargin(0.03)
            pad4.SetBottomMargin(0.23)
            pad4.Draw()
            pad4.cd()
            if not region == "iso_sym":
                h_tight_pull.SetTitleOffset(0.5, "Y")
                h_tight_pull.GetYaxis().SetTitleSize(0.1)
                h_tight_pull.GetYaxis().SetTitle("Pull")
                h_tight_pull.GetYaxis().SetLabelSize(0.06)
                h_tight_pull.GetXaxis().SetLabelSize(0.07)
                h_tight_pull.GetXaxis().SetTitleSize(0.1)
                h_tight_pull.GetXaxis().SetTitle("m_{TP} (GeV)")
                h_tight_pull.Draw('pe')
                #h_tight_pull_error.SetLineColor(ROOT.kGray+2)
                #h_tight_pull_error.SetFillColor(ROOT.kGray+2)
                #h_tight_pull_error.Draw('same e2')
                h_tight_pull.SetLineColor(ROOT.kBlack)
                h_tight_pull.Draw('pe same')
                h_tight_pull.SetMarkerStyle(21)
                h_tight_pull.SetMarkerSize(0.25)
                h_tight_pull.SetStats(0)
                h_tight_pull.GetXaxis().SetRangeUser(0, util.findMaxXVal(loose_fit_as_hist)*3/5)
                h_tight_pull.GetYaxis().SetRangeUser(-abs_pull_max-2, abs_pull_max+2)
                """
                h_tight_pullc.SetMarkerColor(ROOT.kBlue)
                h_tight_pullc.SetMarkerStyle(8)
                h_tight_pullc.SetMarkerSize(0.25)
                h_tight_pullc.Draw('pe same')
                """
                line4 = ROOT.TLine(0.2, 0.6, 0.9, 0.6)
                line4.SetNDC(ROOT.kTRUE)
                line4.SetLineWidth(1)
                line4.SetLineColorAlpha(ROOT.kBlack, 0.5)
                line4.Draw("same")
            ROOT.gPad.Update()
            c1.Print(args.name + ".pdf")
            loose_fit_files[file_counter_loose-1].Close()
            tight_temp_files[file_counter_temp-1].Close()
            tight_deg_files[file_counter_deg-1].Close()
            tight_chi2_files[file_counter_deg-1].Close()

if args.show: input("Finished. Press Enter.")
c1.Print(args.name + ".pdf]")

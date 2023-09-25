from __future__ import print_function
import math
import ROOT
import sys
import os
import argparse
import array
import fitting_utils as util
import scipy
import scipy.stats as stats


def findMaxXVal(hist):
    max_x_val = 0
    for i in range(hist.GetMaximumBin(), hist.GetNbinsX()):
        if hist.GetBinContent(i+1) < 1e-1:
            max_x_val = 0.1*(i+1)
            break
    return max_x_val


def binConverter(test_bin):
    bin_list = test_bin.split(" ")
    return bin_list


# command line options
parser = argparse.ArgumentParser(description="")
parser.add_argument("input", metavar="INPUT", help="input root file")

# plot specification
parser.add_argument("--test", default=False, action="store_true", help="create test plots")
parser.add_argument("--testBin", default=None, help="specify bin to test")
parser.add_argument("--name", default="plots", help="create name for plots pdf")

# parse args
args = parser.parse_args()
infile1 = ROOT.TFile(sys.argv[1])

# other config
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)
ROOT.gStyle.SetLegendFillColor(ROOT.TColor.GetColorTransparent(ROOT.kWhite, 0.01));
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
leg_x1, leg_x2, leg_y1, leg_y2 = 0.7, 0.60, 0.89, 0.89

c1 = ROOT.TCanvas("c1", "c1", 800, 600)
#if not args.sanity: ROOT.TPad.Divide(c1, 1, 2)
c1.Print(args.name + ".pdf[")

eta_regions = ["barrel", "endcap"]
regions = ["iso_sym", "iso_asym", "noniso_sym", "noniso_asym"]

if args.testBin is None: test_regions = ["iso_sym"] 
elif "noniso_asym" in args.testBin: test_regions = ["noniso_asym"]
elif "noniso_sym" in args.testBin: test_regions = ["noniso_sym"]
elif "iso_asym" in args.testBin: test_regions = ["iso_asym"]
elif "iso_sym" in args.testBin: test_regions = ["iso_sym"]

if args.test: regions = test_regions

bins = [20,40,60,80,100,140,180,220,300,380]

if args.testBin is not None: test_bin = binConverter(args.testBin)
chi2_pvalues = []
loose_fit_files = []
scaled_tight_files = []
tight_temp_files = []
tight_deg_files = []
file_counter_loose = 0
file_counter_tight = 0
file_counter_temp = 0
file_counter_deg = 0
for region in regions:  # loop through twoprong sideband regions
    if args.testBin is not None: 
        if not region == test_bin[0]: continue
    for i in range(len(bins)):  # loop through pt bins for a fixed twoprong sideband
        if args.testBin is not None: 
            if not bins[i] == int(test_bin[2]): continue
        for eta_reg in eta_regions:  # loop through eta regions for fixed pt-bin and fixed twoprong sideband
            if args.testBin is not None: 
                if not eta_reg == test_bin[1]: continue

            # Generate correct plots names to access from summed histogram files
            egamma_loose_plots = "plots/twoprong_masspi0_" + region + "_" + eta_reg

            if i == len(bins) - 1: egamma_loose_plots += "_" + str(bins[i]) + "+"
            else: egamma_loose_plots += "_" + str(bins[i]) + "_" + str(bins[i+1])

            # Reference name of the histogram created in the backend 
            egamma_loose_plots += "_loose"

            # Get the histograms from the input file
            h_egamma_loose = infile1.Get(egamma_loose_plots)
            h_egamma_loose.SetBinErrorOption(ROOT.TH1.kPoisson)
            h_egamma_loose.SetLineColor(ROOT.kBlack)
            h_egamma_loose.SetTitle("")
            h_egamma_loose.GetXaxis().SetTitle("")
            h_egamma_loose.GetYaxis().SetTitle("")

            if i == len(bins) - 1: hist_name = region + "_" + eta_reg + "_" + str(bins[i]) + "+"
            else: hist_name = region + "_" + eta_reg + "_" + str(bins[i]) + "_" + str(bins[i+1]) 
           
            # Input loose templates
            if not os.path.exists('loose_fit_hists'): 
                print("ERROR: Must create directory for loose templates titled loose_fit_hists")
            else:
                loose_fit_files.append(ROOT.TFile("loose_fit_hists/" + hist_name + "_loose.root"))
                loose_fit_as_hist = loose_fit_files[file_counter_loose].Get(hist_name+"_loose")
                loose_fit_as_hist.SetLineColor(ROOT.kBlack)
                loose_fit_as_hist.SetTitle("")
                loose_fit_as_hist.GetXaxis().SetTitle("")
                loose_fit_as_hist.GetYaxis().SetTitle("")
                file_counter_loose += 1

            # Input scaled tight data
            if not os.path.exists('scaled_tight_hists'):
                print("ERROR: Must create directory for scaled data titled scaled_tight_hists")
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
            if not os.path.exists('tight_templates/templates'):
                print("ERROR: Must create directory for scaled data as follows: tight_templates/templates")
            else:
                tight_temp_files.append(ROOT.TFile("tight_templates/templates/" + hist_name + "_tight_temp.root"))
                tight_fit_as_hist = tight_temp_files[file_counter_temp].Get(hist_name+"_tight_temp")
                file_counter_temp += 1

            # Input tight template
            if not os.path.exists('tight_templates'):
                print("ERROR: Must create directory for scaled data Bernstein degrees as follows: tight_templates/degrees")
            else:
                tight_deg_files.append(ROOT.TFile("tight_templates/degrees/" + hist_name + "_tight_temp_deg.root"))
                tight_deg_hist = tight_deg_files[file_counter_deg].Get(hist_name+"_tight_temp_deg")
                file_counter_deg += 1
            
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
            legend1 = ROOT.TLegend(0.53, 0.77, 0.88, 0.87)
            legend1.SetX1NDC(ROOT.kTRUE)
            legend1.AddEntry(h_egamma_loose, "Loose Photon Sideband", "l")
            legend1.AddEntry(loose_fit_as_hist, "Landau + Exponential Fit", "l")
            legend1.SetTextSize(0.024)

            legend2 = ROOT.TLegend(0.58, 0.77, 0.93, 0.87)
            if not region == "iso_sym":
                legend2.AddEntry(h_egamma_tight, "Full Selection", "l")
                legend2.AddEntry(tight_fit_as_hist, 'Bernstein, Degree ' + str(int(tight_deg_hist.GetBinContent(1))), 'l')
                legend2.SetTextSize(0.024)
           
            # Create Title 
            title1 = "TwoProng p_{T} " + str(bins[i])
            if i == len(bins)-1: title1 += "+, " 
            else: title1 += "-" + str(bins[i+1]) + ", "
            if eta_reg == "barrel": title1 += "Barrel Photon"
            else: title1 += "Endcap Photon"
           
            title2 = "TwoProng Control Region: "
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
            h_egamma_loose.GetXaxis().SetRangeUser(0, findMaxXVal(loose_fit_as_hist)*3/5)
            legend1.Draw("same")
            # Draw CMS logo
            l1_cms = ROOT.TLatex(0.165, 0.84, "CMS")
            l1_cms.SetNDC(ROOT.kTRUE)
            l1_cms.Draw("same")
            # Draw lumi at top of plot
            l1_lumi = ROOT.TLatex(0.63, 0.915, "137.6 fb^{-1} (13 TeV)")
            l1_lumi.SetNDC(ROOT.kTRUE)
            l1_lumi.SetTextSize(0.035)
            l1_lumi.SetTextFont(42)
            l1_lumi.Draw("same")
            # Draw plot title (first half)
            title1_l = ROOT.TLatex(0.2, 0.96, title1)
            title1_l.SetTextSize(0.04)
            title1_l.SetNDC(ROOT.kTRUE)
            title1_l.Draw("same")
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
                h_egamma_tight.GetXaxis().SetRangeUser(0, findMaxXVal(loose_fit_as_hist)*3/5)
                h_egamma_tight.SetMinimum(0.1)
                h_egamma_tight.Draw("e")
                tight_fit_w_constant.SetLineColor(ROOT.kBlue)
                tight_fit_as_hist.SetLineColor(ROOT.kRed+2)
                tight_fit_as_hist.SetLineWidth(1)
                tight_fit_as_hist.Draw("same hist")
                h_egamma_tight.Draw("e same")
                h_egamma_tight.GetYaxis().SetRangeUser(0.1, h_egamma_tight.GetMaximum()+50)
                l2 = ROOT.TLatex(0.12, 0.84, "CMS")
                l2.SetNDC(ROOT.kTRUE)
                l2.Draw("same")
                l2_lumi = ROOT.TLatex(0.645, 0.915, "137.6 fb^{-1} (13 TeV)")
                l2_lumi.SetNDC(ROOT.kTRUE)
                l2_lumi.SetTextSize(0.035)
                l2_lumi.SetTextFont(42)
                l2_lumi.Draw("same")
                legend2.Draw("same")
                title2_l = ROOT.TLatex(0.05, 0.96, title2)
                title2_l.SetTextSize(0.04)
                title2_l.SetNDC(ROOT.kTRUE)
                title2_l.Draw("")
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
            h_loose_pull.GetXaxis().SetTitle("m_{TwoProng} (GeV)")
            h_loose_pull.SetMarkerStyle(21)
            h_loose_pull.SetMarkerSize(0.25)
            h_loose_pull.SetStats(0)
            h_loose_pull.GetXaxis().SetRangeUser(0, findMaxXVal(loose_fit_as_hist)*3/5)
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
                h_tight_pull.GetXaxis().SetTitle("m_{TwoProng} (GeV)")
                h_tight_pull.Draw('pe')
                #h_tight_pull_error.SetLineColor(ROOT.kGray+2)
                #h_tight_pull_error.SetFillColor(ROOT.kGray+2)
                #h_tight_pull_error.Draw('same e2')
                h_tight_pull.SetLineColor(ROOT.kBlack)
                h_tight_pull.Draw('pe same')
                h_tight_pull.SetMarkerStyle(21)
                h_tight_pull.SetMarkerSize(0.25)
                h_tight_pull.SetStats(0)
                h_tight_pull.GetXaxis().SetRangeUser(0, findMaxXVal(loose_fit_as_hist)*3/5)
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
            if args.testBin is not None: input()
            c1.Print(args.name + ".pdf")

c1.Print(args.name + ".pdf]")





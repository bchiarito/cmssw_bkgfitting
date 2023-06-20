from __future__ import print_function
import math
import ROOT
import sys
import os
import argparse
import array
import fitting_utils as util


def count_nonzero_bins(hist, bound):
  count = 0
  for i in range(hist.GetNbinsX()):
    if hist.GetBinLowEdge(i+1) < bound: continue
    if not hist.GetBinContent(i+1) == 0: count += 1
  return count


def RSS(func, hist, bound, integral=False):
  rss = 0
  by_bin = []
  for i in range(hist.GetNbinsX()):
    if hist.GetBinContent(i+1) == 0: continue
    if hist.GetBinLowEdge(i+1) < bound: continue
    if not integral:
      val = (hist.GetBinContent(i+1) - func.Eval(hist.GetBinCenter(i+1)))**2
    else:
      val = ( hist.GetBinContent(i+1) - (func.Integral(hist.GetBinLowEdge(i+1), hist.GetBinLowEdge(i+1) + hist.GetBinWidth(i+1)))/hist.GetBinWidth(i+1) )**2
    by_bin.append(val)
    rss += val
  return rss, by_bin


def binConverter(test_bin):
    bin_list = test_bin.split(" ")
    return bin_list


# command line options
parser = argparse.ArgumentParser(description="")
parser.add_argument("input", metavar="INPUT", help="input root file")

# plot specification
parser.add_argument("--sanity", "-s", default=False, action="store_true", help="create sanity plots")
parser.add_argument("--test", default=False, action="store_true", help="create test plots")
parser.add_argument("--testBin", default=None, help="specify bin to test")
parser.add_argument("--ratio", default=False, action="store_true", help="create ratio plots instead of fit plots")
parser.add_argument("--fit", default=False, action="store_true", help="create fit plots only")
parser.add_argument("--name", default="plots", help="create name for plots pdf")

# parse args
args = parser.parse_args()
infile1 = ROOT.TFile(sys.argv[1])

# other config
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)
ROOT.gStyle.SetLegendFillColor(ROOT.TColor.GetColorTransparent(ROOT.kWhite, 0.01));
ROOT.gStyle.SetLegendBorderSize(0)
leg_x1, leg_x2, leg_y1, leg_y2 = 0.7, 0.60, 0.89, 0.89

c1 = ROOT.TCanvas("c1", "c1", 800, 600)
#if not args.sanity: ROOT.TPad.Divide(c1, 1, 2)
c1.Print(args.name + ".pdf[")

# pi0: masspi0 plots for all eta regions, barrel, and endcap; pi0_bins: pt-binned masspi0 plots in barrel and endcap; overlay; pt-binned plots with overlayed ratios for each twoprong region
sanity_plots = ["sieie", "pfRelIso03_chg", "hadTow"]  
main_plots = ["pi0_bins"]
test_plots = ["pi0_bins"]
if args.sanity: plots = sanity_plots
elif args.test: plots = test_plots
else: plots = main_plots

eta_regions = ["all", "barrel", "endcap"]
regions = ["iso_sym", "iso_asym", "noniso_sym", "noniso_asym"]
test_regions = ["noniso_sym"]
if args.test: regions = test_regions
photon_regions = ["tight", "loose"]
bins = [20,40,60,70,80,100,120,140,160,180,200,240,300,380,460]


for item in plots:
    if item == "pfRelIso03_chg" or item == "sieie" or item == "hoe" or item == "hadTow":  # sanity plots
        for region in photon_regions:
            iso_sym = "photon_" + region + "_" + item + "_iso_sym"
            iso_asym = "photon_" + region + "_" + item + "_iso_asym"
            noniso_sym = "photon_" + region + "_" + item + "_noniso_sym"
            noniso_asym = "photon_" + region + "_" + item + "_noniso_asym"

            h_iso_sym = infile1.Get("plots/" + iso_sym)
            h_iso_asym = infile1.Get("plots/" + iso_asym)
            h_noniso_sym = infile1.Get("plots/" + noniso_sym)
            h_noniso_asym = infile1.Get("plots/" + noniso_asym)

            h_iso_sym.SetLineColor(ROOT.kRed)
            h_iso_asym.SetLineColor(ROOT.kBlack)
            h_noniso_sym.SetLineColor(ROOT.kGreen+2)
            h_noniso_asym.SetLineColor(ROOT.kBlue)

            if not h_iso_sym.Integral() == 0: h_iso_sym.Scale(1.0/h_iso_sym.Integral())
            if not h_noniso_sym.Integral() == 0: h_noniso_sym.Scale(1.0/h_noniso_sym.Integral())
            if not h_iso_asym.Integral() == 0: h_iso_asym.Scale(1.0/h_iso_asym.Integral())
            if not h_noniso_asym.Integral() == 0: h_noniso_asym.Scale(1.0/h_noniso_asym.Integral())
            
            # Legend creation
            legend = ROOT.TLegend(leg_x1, leg_x2, leg_y1, leg_y2) 
            legend.AddEntry(h_noniso_asym, "NonIso_Asym TwoProng", "l")
            legend.AddEntry(h_iso_sym, "Iso_Sym TwoProng", "l")
            legend.AddEntry(h_iso_asym, "Iso_Asym TwoProng", "l")
            legend.AddEntry(h_noniso_sym, "NonIso_Sym TwoProng", "l")

            title = region + " Photon"
            h_iso_sym.SetTitle(title)

            # Draw plots
            c1.cd(1)
            h_iso_sym.Draw("e")
            h_noniso_asym.Draw("samee")
            h_iso_asym.Draw("samee")
            h_noniso_sym.Draw("samee")
            ROOT.gPad.SetLogy()
            if item == "eta": ROOT.gPad.SetGridx(1)
            legend.Draw("same")

            if item == "sieie": h_iso_sym.GetXaxis().SetRangeUser(0, 0.1)
            elif item == "hoe" and region == "tight": h_iso_sym.GetXaxis().SetRangeUser(0, 0.4)
            elif item == "hoe" and region == "loose": h_iso_sym.GetXaxis().SetRangeUser(0, 2)
            elif item == "pfRelIso03_chg" and region == "tight": h_iso_sym.GetXaxis().SetRangeUser(0, 0.1)
            elif item == "pfRelIso03_chg" and region == "loose": h_iso_sym.GetXaxis().SetRangeUser(0, 0.2)
            elif item == "hadTow" and region == "tight": h_iso_sym.GetXaxis().SetRangeUser(0, 0.1) 
            elif item == "hadTow" and region == "loose": h_iso_sym.GetXaxis().SetRangeUser(0, 0.2) 
            c1.Print(args.name + ".pdf")
    elif item == "pi0":  # un-pt-binned massPi0 plots
        ROOT.TPad.Divide(c1, 1, 2)
        for region in regions:  # loop through twoprong regions
            for eta_reg in eta_regions:  # loop through eta regions for a fixed twoprong sideband
                if not eta_reg == "barrel" and not eta_reg == "endcap":
                    h_egamma_tight = infile1.Get("plots/twoprong_masspi0_" + region + "_tight")
                    h_egamma_loose = infile1.Get("plots/twoprong_masspi0_" + region + "_loose")
                else:
                    h_egamma_tight = infile1.Get("plots/twoprong_masspi0_" + region + "_" + eta_reg + "_tight")
                    h_egamma_loose = infile1.Get("plots/twoprong_masspi0_" + region + "_" + eta_reg + "_loose")
                
                h_egamma_tight.SetLineColor(ROOT.kBlack)
                h_egamma_loose.SetLineColor(ROOT.kGreen+2)
                h_egamma_loose.SetFillColor(ROOT.kGreen+2)
                 
                if not h_egamma_tight.Integral() == 0: h_egamma_tight.Scale(1.0/h_egamma_tight.Integral())
                if not h_egamma_loose.Integral() == 0: h_egamma_loose.Scale(1.0/h_egamma_loose.Integral())
                
                h_ratio = h_egamma_tight.Clone()
                h_ratio.Reset()
                h_ratio.SetLineColor(ROOT.kBlack)
                h_ratio.Divide(h_egamma_tight, h_egamma_loose)
               
                # Create title for plot 
                title = region + " Twoprong"
                if eta_reg == "barrel": title += ", Barrel"
                elif eta_reg == "endcap": title += ", Endcap"
               
                # Legend creation
                legend = ROOT.TLegend(leg_x1, leg_x2, leg_y1, leg_y2)
                legend.AddEntry(h_egamma_tight, "Tight Photon, " + str(h_egamma_tight.GetEntries()), "l")
                legend.AddEntry(h_egamma_loose, "Loose Photon, " + str(h_egamma_loose.GetEntries()), "f")

                h_egamma_tight.SetTitle(title)
                h_ratio.SetTitle("Tight / Loose Photon")  
                
                # Draw plots
                c1.cd(1)
                h_egamma_tight.Draw("")
                mc_stack = ROOT.THStack('hs', 'hs')
                mc_stack.Add(h_egamma_loose)
                mc_stack.Draw("hist same")
                h_egamma_tight.Draw("samee")
                ROOT.gPad.SetLogy()
                h_egamma_tight.GetXaxis().SetRangeUser(0, 26)
                legend.Draw("same")
                
                c1.cd(2)
                h_ratio.Draw("e")
                h_ratio.GetXaxis().SetRangeUser(0, 26)
                h_ratio.GetYaxis().SetRangeUser(-2, 4)
                h_ratio.SetStats(0)
                ROOT.gPad.SetGridy(1)
                ROOT.gPad.Update()
                c1.Print(args.name + ".pdf")    
    elif item == "pi0_bins":
        if args.ratio: ROOT.TPad.Divide(c1, 1, 2)
        if args.testBin is not None: test_bin = binConverter(args.testBin)
        for region in regions:  # loop through twoprong sideband regions
            if args.testBin is not None: 
                if not region == test_bin[0]: continue
            for i in range(len(bins)):  # loop through pt bins for a fixed twoprong sideband
                if args.testBin is not None: 
                    if not bins[i] == int(test_bin[2]): continue
                for eta_reg in eta_regions:  # loop through eta regions for fixed pt-bin and fixed twoprong sideband
                    if args.testBin is not None: 
                        if not eta_reg == test_bin[1]: continue
                    if not eta_reg == "barrel" and not eta_reg == "endcap": continue  # no pt-bin plots for barrel and endcap combined, so skip this case

                    # Generate correct plots names to access from summed histogram files
                    egamma_tight_plots = "plots/twoprong_masspi0_" + region + "_" + eta_reg
                    egamma_loose_plots = "plots/twoprong_masspi0_" + region + "_" + eta_reg
                
                    if i == len(bins) - 1:
                        egamma_tight_plots += "_" + str(bins[i]) + "+"
                        egamma_loose_plots += "_" + str(bins[i]) + "+"
                    else:
                        egamma_tight_plots += "_" + str(bins[i]) + "_" + str(bins[i+1])
                        egamma_loose_plots += "_" + str(bins[i]) + "_" + str(bins[i+1]) 
                    
                    # Reference name of the histogram created in the backend 
                    egamma_tight_plots += "_tight"
                    egamma_loose_plots += "_loose"
                    
                    # Get the histograms from the input file
                    h_egamma_tight = infile1.Get(egamma_tight_plots)
                    h_egamma_loose = infile1.Get(egamma_loose_plots)
                    
                    # Configure display options
                    h_egamma_tight.SetLineColor(ROOT.kBlack)
                    h_egamma_loose.SetLineColor(ROOT.kBlack)

                    if args.ratio:  # create unfitted ratio plots between tight and loose photons
                        if not h_egamma_tight.Integral() == 0: h_egamma_tight.Scale(1.0/h_egamma_tight.Integral())
                        if not h_egamma_loose.Integral() == 0: h_egamma_loose.Scale(1.0/h_egamma_loose.Integral())
                         
                        h_ratio = h_egamma_tight.Clone()
                        h_ratio.Reset()
                        h_ratio.SetLineColor(ROOT.kBlack)
                        h_ratio.Divide(h_egamma_tight, h_egamma_loose)
                        h_ratio.SetTitle("Tight / Loose")
                        
                        # Create legend
                        legend1 = ROOT.TLegend(0.65, 0.45, 0.9, 0.6)
                        legend1.AddEntry(h_egamma_loose, "Loose Photon, " + str(h_egamma_loose.GetEntries()), "l")
                        legend1.AddEntry(h_egamma_tight, "Tight Photon, " + str(h_egamma_tight.GetEntries()), "l")
                        
                        # Create title for plot 
                        title = region + " Twoprong"
                        if eta_reg == "barrel": title += ", Barrel"
                        elif eta_reg == "endcap": title += ", Endcap"
                        if i == len(bins) - 1: title += ", pt > " + str(bins[i])
                        else: title += ", " + str(bins[i]) + " < pt < " + str(bins[i+1])
                        
                        c1.cd(1)
                        h_egamma_loose.SetTitle(title)
                        h_egamma_loose.SetMaximum()
                        h_egamma_loose.Draw("e")
                        h_egamma_tight.Draw("samee")
                        ROOT.gPad.SetLogy()
                        if bins[i] < 80: h_egamma_loose.GetXaxis().SetRangeUser(0, 5)
                        elif bins[i] < 120: h_egamma_loose.GetXaxis().SetRangeUser(0, 10)
                        elif bins[i] < 200: h_egamma_loose.GetXaxis().SetRangeUser(0, 15)
                        elif bins[i] < 380: h_egamma_loose.GetXaxis().SetRangeUser(0, 20)
                        else: h_egamma_loose.GetXaxis().SetRangeUser(0, 26)
                        legend1.Draw("same")

                        c1.cd(2)
                        h_ratio.Draw("e")
                        if bins[i] < 80: h_ratio.GetXaxis().SetRangeUser(0, 5)
                        elif bins[i] < 120: h_ratio.GetXaxis().SetRangeUser(0, 10)
                        elif bins[i] < 200: h_ratio.GetXaxis().SetRangeUser(0, 15)
                        elif bins[i] < 380: h_ratio.GetXaxis().SetRangeUser(0, 20)
                        else: h_ratio.GetXaxis().SetRangeUser(0, 26)
                        h_ratio.GetYaxis().SetRangeUser(-2, 4)
                        h_ratio.SetStats(0)
                        ROOT.gPad.SetGridy(1)
                        ROOT.gPad.Update()
                        c1.Print(args.name + ".pdf")
                    else:
                        if i == len(bins) - 1: print("############### PT BIN: " + str(bins[i]) + "+, " + eta_reg.upper() + " ###############")
                        else: print("############### PT BIN: " + str(bins[i]) + "-" + str(bins[i+1]) + ", " + eta_reg.upper() + " ###############")
                        
                        # varibles for f-test
                        fits = []
                        histo = ''
                        num_param = []
                        by_bins = []
                        bound1s = []

                        ### new idea, fit only rising and falling ###

                        # determine left and right bounds
                        # fit landaus from first_bin to left_bin
                        # fit exps from right_bin to last_bin
                        ENTRIES_CUTOFF = 1000
                        first_bin = 4
                        left_bin = 0
                        for b in range(h_egamma_loose.GetNbinsX()):
                          if h_egamma_loose.GetBinContent(b+1)<ENTRIES_CUTOFF: continue
                          else:
                            left_bin = b + 1
                            break
                        right_bin = h_egamma_loose.GetNbinsX()+1
                        for b in reversed(range(h_egamma_loose.GetNbinsX())):
                          if h_egamma_loose.GetBinContent(b+1)<ENTRIES_CUTOFF: continue
                          else:
                            right_bin = b + 1
                            break
                        last_bin = h_egamma_loose.GetNbinsX() + 1
                        for b in reversed(range(h_egamma_loose.GetNbinsX())):
                          if h_egamma_loose.GetBinContent(b+1)==0: continue
                          else:
                            last_bin = b + 1
                            break
                        # convert to x-ranges
                        first = h_egamma_loose.GetBinLowEdge(first_bin)
                        left = h_egamma_loose.GetBinLowEdge(left_bin+1)
                        right = h_egamma_loose.GetBinLowEdge(right_bin)
                        last = h_egamma_loose.GetBinLowEdge(last_bin+1)
                        
                        landau_range = 2
                        exp_range = 3
                        
                        # loop through combinations of 1,2 landaus and 1,2,3 exponentials
                        for j in range(landau_range):
                            for k in range(exp_range):
                                nLandau = j+1
                                nExp = k+1
                                if not nLandau == 1 or not nExp == 2: continue 
                                
                                landau_guess = None
                                exp_guess = None
                                
                                if nLandau == 1 and nExp == 2:
                                    pass
                                if nLandau == 2 and nExp == 3:
                                    landau_guess = [37840, 0.7493, 0.1121, 0.555, 1.371, 0.3455]
                                    exp_guess = [2e+05, -2.082, 1.555, -1.551, 1.5,  -1]
                                
                                func_rising, fitresult_rising = util.fit_hist(h_egamma_loose, 'landau', first, left, N=nLandau, initial_guesses=landau_guess)
                                rising_fit_as_hist = util.TemplateToHistogram(func_rising, 1000, 0, 50)
                                h_egamma_loose.Draw()
                                c1.Update()
                                stats1 = h_egamma_loose.GetListOfFunctions().FindObject("stats").Clone("stats1")
                                c1.Clear()
                                c1.Update()
                                stats1.SetY1NDC(.5)
                                stats1.SetY2NDC(.7)

                                func_falling, fitresult_falling = util.fit_hist(h_egamma_loose, 'exp', right, last, N=nExp, initial_guesses=exp_guess)
                                falling_fit_as_hist = util.TemplateToHistogram(func_falling, 1000, 0, 50)
                                h_egamma_loose.Draw()
                                c1.Update()
                                stats2 = h_egamma_loose.GetListOfFunctions().FindObject("stats").Clone("stats2")
                                c1.Clear()
                                c1.Update()

                                # create overall fitted histogram as: rising - bulk - falling
                                loose_fit_as_hist = h_egamma_loose.Clone()
                                loose_fit_as_hist.Reset()
                                for b in range(h_egamma_loose.GetNbinsX()):
                                    if b < first_bin:
                                        loose_fit_as_hist.SetBinContent(b+1, 0)
                                    elif b < left_bin:
                                        loose_fit_as_hist.SetBinContent(b+1, rising_fit_as_hist.GetBinContent(b+1))
                                    elif b <= right_bin:
                                        loose_fit_as_hist.SetBinContent(b+1, h_egamma_loose.GetBinContent(b+1))
                                    elif b <= last_bin:
                                        loose_fit_as_hist.SetBinContent(b+1, falling_fit_as_hist.GetBinContent(b+1))
                                    else:
                                        loose_fit_as_hist.SetBinContent(b+1, 0)
                        
                                fitted_func = util.HistogramToFunction(loose_fit_as_hist)
                                fitted_func_times_constant, _ = util.MultiplyWithPolyToTF1(fitted_func, 0, cheb=0)
                                h_egamma_tight.Fit(fitted_func_times_constant, '0L') 
                                tight_fit_w_constant = util.TemplateToHistogram(fitted_func_times_constant, 1000, 0, 50)

                                func_with_poly, _ = util.MultiplyWithPolyToTF1(fitted_func, 3, cheb=0)
                                h_egamma_tight.Fit(func_with_poly, '0L') 
                                tight_fit_as_hist = util.TemplateToHistogram(func_with_poly, 1000, 0, 50)

                                h_loose_pull_num = h_egamma_loose.Clone()
                                h_loose_pull_num.Reset()
                                h_loose_pull = h_egamma_loose.Clone()
                                h_loose_pull.Reset()
                                h_loose_pull_num.Add(h_egamma_loose, loose_fit_as_hist, 1, -1)  # Numerator of pull hist is data - fit

                                for j in range(h_loose_pull_num.GetNbinsX()): 
                                    if h_egamma_loose.GetBinContent(j+1) == 0: sqrt_err = 1.8
                                    else: sqrt_err = h_egamma_loose.GetBinError(j+1)
                                    h_loose_pull.SetBinContent(j+1, h_loose_pull_num.GetBinContent(j+1)/sqrt_err)
                                    h_loose_pull.SetBinError(j+1, 1)
                                
                                h_tight_pull_num = h_egamma_tight.Clone()
                                h_tight_pull_num.Reset()
                                h_tight_pull = h_egamma_tight.Clone()
                                h_tight_pull.Reset()
                                h_tight_pull_num.Add(h_egamma_tight, tight_fit_as_hist, 1, -1)  # Numerator of pull hist is data - fit
                                
                                for j in range(h_tight_pull_num.GetNbinsX()): 
                                    if h_egamma_tight.GetBinContent(j+1) == 0: sqrt_err = 1.8
                                    else: sqrt_err = h_egamma_tight.GetBinError(j+1)
                                    h_tight_pull.SetBinContent(j+1, h_tight_pull_num.GetBinContent(j+1)/sqrt_err)
                                    h_tight_pull.SetBinError(j+1, 1)

                                # pull for tight fit with constant
                                h_tight_pullc_num = h_egamma_tight.Clone()
                                h_tight_pullc_num.Reset()
                                h_tight_pullc = h_egamma_tight.Clone()
                                h_tight_pullc.Reset()
                                h_tight_pullc_num.Add(h_egamma_tight, tight_fit_w_constant, 1, -1)  # Numerator of pull hist is data - fit
                                
                                for j in range(h_tight_pullc_num.GetNbinsX()): 
                                    if h_egamma_tight.GetBinContent(j+1) == 0: sqrt_err = 1.8
                                    else: sqrt_err = h_egamma_tight.GetBinError(j+1)
                                    h_tight_pullc.SetBinContent(j+1, h_tight_pullc_num.GetBinContent(j+1)/sqrt_err)
                                    h_tight_pullc.SetBinError(j+1, 0) # no error bar
                                
                                # Create title for plot 
                                title = region + " Twoprong"
                                if eta_reg == "barrel": title += ", Barrel"
                                elif eta_reg == "endcap": title += ", Endcap"
                                if i == len(bins) - 1: title += ", pt > " + str(bins[i])
                                else: title += ", " + str(bins[i]) + " < pt < " + str(bins[i+1]) 
                                if nExp == 1: title += ", 1 exp"
                                elif nExp == 2: title += ", 2 exp"
                                elif nExp == 3: title += ", 3 exp"

                                # Legend creation
                                legend1 = ROOT.TLegend(0.65, 0.35, 0.9, 0.45)
                                legend1.AddEntry(h_egamma_loose, "Loose Photon, " + str(h_egamma_loose.GetEntries()), "l")
                                #if not args.ratio: legend1.AddEntry(0, "Chi2/NDF: " + str(chi2 / ndf), "")
                                legend2 = ROOT.TLegend(0.65, 0.35, 0.9, 0.45)
                                legend2.AddEntry(h_egamma_tight, "Tight Photon, " + str(h_egamma_tight.GetEntries()), "l")
                                legend2.AddEntry(tight_fit_as_hist, "Fitted Tight", "f")
                                
                                # Draw plots
                                if args.fit: c1.cd(1)
                                else:
                                  c1.cd()
                                  pad1 = ROOT.TPad('pad1', 'pad1', 0, 0.3, 0.5, 1)
                                  pad1.Draw()
                                  pad1.cd()
                                
                                h_egamma_loose.SetTitle(title)
                                h_egamma_loose.SetMaximum()
                                h_egamma_loose.Draw("e")
                                loose_fit_as_hist.SetLineColor(ROOT.kRed+1)
                                loose_fit_as_hist.Draw("same")
                                # f2.Draw("same")
                                stats1.Draw()
                                stats2.Draw()
                                ROOT.gPad.SetLogy()
                                if bins[i] < 80: h_egamma_loose.GetXaxis().SetRangeUser(0, 5)
                                elif bins[i] < 120: h_egamma_loose.GetXaxis().SetRangeUser(0, 10)
                                elif bins[i] < 200: h_egamma_loose.GetXaxis().SetRangeUser(0, 15)
                                elif bins[i] < 380: h_egamma_loose.GetXaxis().SetRangeUser(0, 20)
                                else: h_egamma_loose.GetXaxis().SetRangeUser(0, 26)
                                legend1.Draw("same")
                                ROOT.gPad.Update()
                                
                                if not args.fit:
                                    if not region == "iso_sym":
                                        c1.cd()
                                        pad2 = ROOT.TPad('pad2', 'pad2', 0.5, 0.3, 1, 1)
                                        pad2.Draw()
                                        pad2.cd()
                                        h_egamma_tight.Draw("e")
                                        tight_fit_w_constant.SetLineColor(ROOT.kBlue)
                                        tight_fit_as_hist.SetLineColor(ROOT.kRed)
                                        tight_fit_as_hist.SetLineWidth(1)
                                        tight_fit_as_hist_errorbars = tight_fit_as_hist.Clone()
                                        tight_fit_as_hist_errorbars.SetFillColor(ROOT.kRed+2)
                                        tight_fit_as_hist_errorbars.Draw("same e2")
                                        tight_fit_as_hist.Draw("same hist")
                                        tight_fit_w_constant.Draw('same')
                                        h_egamma_tight.Draw("e same")
                                        ROOT.gPad.SetLogy()
                                        if bins[i] < 80: h_egamma_tight.GetXaxis().SetRangeUser(0, 5)
                                        elif bins[i] < 120: h_egamma_tight.GetXaxis().SetRangeUser(0, 10)
                                        elif bins[i] < 200: h_egamma_tight.GetXaxis().SetRangeUser(0, 15)
                                        elif bins[i] < 380: h_egamma_tight.GetXaxis().SetRangeUser(0, 20)
                                        else: h_egamma_tight.GetXaxis().SetRangeUser(0, 26)
                                        legend2.Draw("same")
                                
                                if not args.fit:
                                    c1.cd()
                                    pad3 = ROOT.TPad('pad3', 'pad3', 0, 0, 0.5, 0.3)
                                    pad3.Draw()
                                    pad3.cd()
                                    h_loose_pull.SetTitle("(Loose - Fit) / Error")
                                    h_loose_pull.SetLineColor(ROOT.kBlack)
                                    h_loose_pull.Draw('pe')
                                    h_loose_pull.SetMarkerStyle(8)
                                    h_loose_pull.SetMarkerSize(0.25)
                                    h_loose_pull.GetYaxis().SetRangeUser(-15, 15)
                                    h_loose_pull.SetStats(0)
                                    if bins[i] < 80: h_loose_pull.GetXaxis().SetRangeUser(0, 5)
                                    elif bins[i] < 120: h_loose_pull.GetXaxis().SetRangeUser(0, 10)
                                    elif bins[i] < 200: h_loose_pull.GetXaxis().SetRangeUser(0, 15)
                                    elif bins[i] < 380: h_loose_pull.GetXaxis().SetRangeUser(0, 20)
                                    else: h_loose_pull.GetXaxis().SetRangeUser(0, 26)
                                    
                                    if not region == "iso_sym":
                                        c1.cd()
                                        pad4 = ROOT.TPad('pad4', 'pad4', 0.5, 0, 1, 0.3)
                                        pad4.Draw()
                                        pad4.cd()
                                        h_tight_pull.SetTitle("(Tight - Fit) / Error")
                                        h_tight_pull.SetLineColor(ROOT.kBlack)
                                        h_tight_pull.Draw('pe')
                                        h_tight_pull.SetMarkerStyle(8)
                                        h_tight_pull.SetMarkerSize(0.25)
                                        h_tight_pull.GetYaxis().SetRangeUser(-15, 15)
                                        h_tight_pull.SetStats(0)
                                        if bins[i] < 80: h_tight_pull.GetXaxis().SetRangeUser(0, 5)
                                        elif bins[i] < 120: h_tight_pull.GetXaxis().SetRangeUser(0, 10)
                                        elif bins[i] < 200: h_tight_pull.GetXaxis().SetRangeUser(0, 15)
                                        elif bins[i] < 380: h_tight_pull.GetXaxis().SetRangeUser(0, 20)
                                        else: h_tight_pull.GetXaxis().SetRangeUser(0, 26)
                                        h_tight_pullc.SetMarkerColor(ROOT.kBlue)
                                        h_tight_pullc.SetMarkerStyle(8)
                                        h_tight_pullc.SetMarkerSize(0.25)
                                        h_tight_pullc.Draw('pe same')
                                
                                ROOT.gPad.Update()
                                if args.testBin is not None: raw_input()
                                c1.Print(args.name + ".pdf")
                        # after loop on fits, do f-test
                        """
                        lower_bound = min(*bound1s)
                        print(bound1s)
                        rss1, by_bin = RSS(fits[0], saved_loose_histo, lower_bound); by_bins.append(by_bin)
                        rss2, by_bin = RSS(fits[1], saved_loose_histo, lower_bound); by_bins.append(by_bin)
                        rss3, by_bin = RSS(fits[2], saved_loose_histo, lower_bound); by_bins.append(by_bin)
                        rss4, by_bin = RSS(fits[3], saved_loose_histo, lower_bound); by_bins.append(by_bin)
                        p1 = num_param[0]
                        p2 = num_param[1]
                        p3 = num_param[2]
                        p4 = num_param[3]
                        n = count_nonzero_bins(h_egamma_loose, lower_bound)
                        F21 = ((rss1 - rss2)/(p2 - p1)) / (rss2/(n - p2))
                        F31 = ((rss1 - rss3)/(p3 - p1)) / (rss3/(n - p3))
                        F32 = ((rss2 - rss3)/(p3 - p2)) / (rss3/(n - p3))
                        F41 = ((rss1 - rss4)/(p4 - p1)) / (rss4/(n - p4))
                        F42 = ((rss2 - rss4)/(p4 - p2)) / (rss4/(n - p4))
                        F43 = ((rss3 - rss4)/(p4 - p3)) / (rss4/(n - p4))
                        print str(p1)+" "+str(p2)+" "+str(p3)+" "+str(p4)
                        print str(rss1)+" "+str(rss2)+" "+str(rss3)+" "+str(rss4)
                        print str(n)
                        print str(lower_bound)
                        print ""
                        for b in range(len(by_bins[0])):
                          s = str(b)+": "
                          for array in by_bins:
                            s = s+" "+"{:.3}".format(array[b])
                          print s
                        print ""
                        print "F21: "+ str(F21)
                        print "  ({}, {}) degrees of freedom".format(p2-p1, n-p2)
                        print "F31: "+ str(F21)
                        print "  ({}, {}) degrees of freedom".format(p3-p1, n-p3)
                        print "F32: "+ str(F32)
                        print "  ({}, {}) degrees of freedom".format(p3-p2, n-p3)
                        print "F41: "+ str(F41)
                        print "  ({}, {}) degrees of freedom".format(p4-p1, n-p4)
                        print "F42: "+ str(F42)
                        print "  ({}, {}) degrees of freedom".format(p4-p2, n-p4)
                        print "F43: "+ str(F43)
                        print "  ({}, {}) degrees of freedom".format(p4-p3, n-p4)
                        raw_input()
                        """
                            
c1.Print(args.name + ".pdf]")
infile1.Close()
if args.testBin is not None: raw_input()

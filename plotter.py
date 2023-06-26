from __future__ import print_function
import math
import ROOT
import sys
import os
import argparse
import array
import fitting_utils as util

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
parser.add_argument("--ftest", default=None, help="format: '<CHEB_TYPE> <MAXDEGREE>', default is no f-test")

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
test_regions = ["iso_asym"]
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
                        exp_range = 4

                        # loop through combinations of 1,2 landaus and 1,2,3 exponentials
                        for k in range(landau_range):
                            for l in range(exp_range):
                                nLandau = k+1
                                nExp = l+1
                                
                                old_method = False
                                landau_guess = None
                                exp_guess = None
                                
                                if region == "noniso_sym":
                                    if not nExp == 3 or not nLandau == 2: continue
                                if region == "iso_asym":
                                    if not nExp == 3 or not nLandau == 1: continue

                                if region == "iso_sym" or region == "iso_asym": 
                                    old_method = True
                                    guesses = None

                                # NonIso-Sym Guesses
                                if region == "noniso_sym" and eta_reg == "barrel" and bins[i] == 20:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [3.491e+06, 0.8838, 0.145, 0.2689, 1.171, 0.2137] 
                                        exp_guess = [8.593e+05, -4.278, 1.595, -7.263, 0.2, -9.997]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "endcap" and bins[i] == 20:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [1.705e+04, 1.108, 0.1949, 0.5, 1.108, 0.1949] 
                                        exp_guess = [4182, -0.9553, 1.4, -5.89, 0.5, -1]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "barrel" and bins[i] == 40:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [1.705e+04, 1.33, 0.205, 0.3737, 1.085, 0.1943]
                                        exp_guess = [3.636e+05, -2.809, 2.176, -4.8, 0.6, -1.807]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "endcap" and bins[i] == 40:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [1.123e+05, 0.875, 0.1524, 0.5311, 0.4702, 0.4837]
                                        exp_guess = [4.549e+04, -2.09, 1.915, -3.572, 0.6484, -6.328]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "endcap" and bins[i] == 60:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [1.021e+04, 0.7237, 0.107, 0.2622, 1.169, 0.2502]
                                        exp_guess = [5970, -0.9294, 1.702, -2.435, 1.262, -5.012]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "barrel" and bins[i] == 70:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [1.857e+04, 0.682, 0.09099, 0.5348, 1.346, 0.2239]
                                        exp_guess = [1.349e+04, -1.075, 2.269, -2.444, 1.235, -5.598]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "endcap" and bins[i] == 70:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [1.163e+04, 0.7851, 0.1258, 0.1535, 1.196, 0.2558]
                                        exp_guess = [1.045e+04, -1.111, 1.645, -2.263, 1.974, -6.166]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "barrel" and bins[i] == 80:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [1.826e+04, 0.6239, 0.07162, 0.5, 1.336, 0.2698]
                                        exp_guess = [1.506e+04, -0.9609, 2.4, -2.24, 2.2, -0.4279]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "endcap" and bins[i] == 80:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [1.826e+04, 0.6239, 0.07162, 0.6, 1.336, 0.2698]
                                        exp_guess = [1.506e+04, -0.9609, 2.4, -2.24, 2.2, -0.4279]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "barrel" and bins[i] == 100:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [1.705e+04, 0.6496, 0.0804, 0.7, 0.6496, 0.0804] 
                                        exp_guess = [1.755e+04, -0.92, 3.125, -1.706, 2, -1]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "endcap" and bins[i] == 100:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [1.705e+04, 0.6496, 0.2, 0.7, 0.6496, 0.2]
                                        exp_guess = [1.755e+04, -0.92, 2.4, -1.706, 2.6, -1]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "barrel" and bins[i] == 120:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [1.705e+04, 0.6496, 0.3, 0.7, 0.6496, 0.3] 
                                        exp_guess = [1.755e+04, -0.92, 3.125, -1.706, 2, -1]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "endcap" and bins[i] == 120:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [1.705e+04, 0.6496, 0.3, 0.7, 0.6496, 0.3] 
                                        exp_guess = [1.755e+04, -0.92, 3.125, -1.706, 2, -1]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "barrel" and bins[i] == 140:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [8972, 0.8638, 0.1507, 0.8, 1.079, 0.2261]
                                        exp_guess = [9448, -1.041, 3.2, -1.423, 3, -3]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "endcap" and bins[i] == 140:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [1.761e+04, 0.7322, 0.1067, 0.4, 1.079, 0.2261]
                                        exp_guess = [9448, -1.041, 3.2, -1.423, 3, -3]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "barrel" and bins[i] == 160:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [9535, 0.6861, 0.08948, 0.1568, 1.522, 0.3291]
                                        exp_guess = [5.971e+04, -1.224, 4, -3.118, 4, -4.477]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "endcap" and bins[i] == 160:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [7468, 0.9087, 0.1663, 0.1634, 1.12, 0.2401]
                                        exp_guess = [1.107e+04, -1.035, 3, -1.572, 2, -7.523]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "barrel" and bins[i] == 180:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [1.183e+04, 0.7977, 0.1293, 0.1502, 1.249, 0.27]
                                        exp_guess = [4.739e+04, -1.172, 3.5, -3.605, 5.5, -3.677]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "endcap" and bins[i] == 180:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [7031, 1.004, 0.1965, 0.2383, 1.071, 0.2191]
                                        exp_guess = [1.742e+04, -1.198, 3.5, -3.296, 2.5, -1.366]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "barrel" and bins[i] == 200:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [1.204e+04, 0.7742, 0.1196, 0.4986, 1.483, 0.3501]
                                        exp_guess = [4.491e+04, -1.059, 5.5, -3.574, 4.5, -2.234]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "endcap" and bins[i] == 200:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [1.204e+04, 0.7742, 0.1196, 0.8, 1.483, 0.3501]
                                        exp_guess = [4.491e+04, -1.059, 5.5, -3.574, 4.5, -2.234]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "barrel" and bins[i] == 300:
                                    if nLandau == 2 and nExp == 3:
                                        landau_guess = [8442, 1.149, 0.2453, 0.8, 1.115, 0.2326]
                                        exp_guess = [3.14e+04, -1.023, 6.75, -0.5922, 1.25, -1]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "endcap" and bins[i] == 300:
                                    nExp += 1
                                    if nLandau == 2 and nExp == 4:
                                        old_method = True
                                        guesses = [1.706e+04, 1.529, 0.3766, 1.706e+04, 1.529, 0.3766, -0.01528, -0.7213, -1.017, -0.5931, 0.5, 0.8, 0.5, 1.27, 3.6]  
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "barrel" and bins[i] == 380:
                                    nExp += 1
                                    if nLandau == 2 and nExp == 4:
                                        old_method = True
                                        guesses = [4848, 1.63, 0.4089, 4848, 1.63, 0.4089, -3.553e-13, -0.6753, -0.8161, -0.5092, 0.5, 1.404, 0.5423, 1.478, 4.575]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "endcap" and bins[i] == 380:
                                    nExp += 1
                                    if nLandau == 2 and nExp == 4:
                                        old_method = True
                                        guesses = [890.4, 1.723, 0.4403, 890.4, 1.723, 0.4403, -1.441, -0.9773, -0.5518, -0.2026, 0.5, 4.1, 0.5, 2.8, 7]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "barrel" and bins[i] == 460:
                                    nExp += 1
                                    if nLandau == 2 and nExp == 4:
                                        old_method = True
                                        guesses = [2487, 1.731, 0.4342, 2487, 1.731, 0.4342, -0.08184, -0.6556, -0.6027, -0.3448, 0.5, 1.187, 0.5, 3.168, 5.617]
                                    else: continue
                                if region == "noniso_sym" and eta_reg == "endcap" and bins[i] == 460:
                                    nExp += 1
                                    if nLandau == 2 and nExp == 4:
                                        old_method = True
                                        guesses = [301.3, 1.842, 0.482, 301.3, 1.842, 0.482, -1.031, -0.5918, -1.591, -0.2211, 0.5, 3.275, 0.8257, 6, 4]
                                    else: continue
                                
                                if old_method:
                                    N = str(nLandau) + str(nExp)
                                    func_full, fitresult_full = util.fit_hist(h_egamma_loose, 'full', 0, 50, int(N), initial_guesses=guesses)
                                    loose_fit_as_hist = util.TemplateToHistogram(func_full, 1000, 0, 50)
                                else:
                                    func_rising, fitresult_rising = util.fit_hist(h_egamma_loose, 'landau', first, left, N=nLandau, initial_guesses=landau_guess)
                                    rising_fit_as_hist = util.TemplateToHistogram(func_rising, 1000, 0, 50)
                                    h_egamma_loose.Draw()
                                    c1.Update()
                                    stats1 = h_egamma_loose.GetListOfFunctions().FindObject("stats").Clone("stats1")
                                    c1.Clear()
                                    c1.Update()
                                    stats1.SetY1NDC(.4)
                                    stats1.SetY2NDC(.6)

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

                                FTEST = True if args.ftest else False
                                if not FTEST:
                                    func_with_poly, _ = util.MultiplyWithPolyToTF1(fitted_func, 2, cheb=1)
                                    h_egamma_tight.Fit(func_with_poly, '0L') 
                                    tight_fit_as_hist = util.TemplateToHistogram(func_with_poly, 1000, 0, 50)

                                if FTEST:
                                    parse = args.ftest.split()
                                    NUM_DEGREES = int(parse[1])
                                    CHEB_TYPE = int(parse[0])
                                    fitfuncs = []
                                    fitresults = []
                                    statboxes = []
                                    for degree in range(NUM_DEGREES+1):
                                        func_with_poly, _ = util.MultiplyWithPolyToTF1(fitted_func, degree, cheb=CHEB_TYPE)
                                        fitresult = h_egamma_tight.Fit(func_with_poly, '0SL')
                                        tight_fit_as_hist = util.TemplateToHistogram(func_with_poly, 1000, 0, 50)

                                        fitfuncs.append(func_with_poly)
                                        fitresults.append(fitresult)
                                        h_egamma_tight.Draw()
                                        c1.Update()
                                        statboxes.append(h_egamma_tight.GetListOfFunctions().FindObject("stats").Clone("stat"+str(degree)))
                                        c1.Clear()
                                        c1.Update()
                                    h_egamma_tight.SetStats(0)

                                    best_d = 0
                                    for d2 in range(NUM_DEGREES+1):
                                        for d1 in range(d2):
                                            if not d1 == best_d: continue
                                            decision, ftest, target, dof1, dof2 = util.ftest(
                                                h_egamma_tight, fitfuncs[d2], fitresults[d2], fitfuncs[d1], fitresults[d1])
                                            print(d2, '>', d1, decision)
                                            print('  F={} target={}'.format(ftest, target), '({}, {}) dof'.format(dof1, dof2))
                                            print('')
                                            if decision: best_d = d2
                                    print('Best: ', best_d)
                                    func_with_poly = fitfuncs[best_d]
                                    tight_fit_as_hist = util.TemplateToHistogram(func_with_poly, 1000, 0, 50)
                                    tight_stat = statboxes[best_d]

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
                                if not old_method:
                                    if nLandau == 1: title += ", 1 land"
                                    elif nLandau == 2: title += ", 2 land"
                                    if nExp == 1: title += ", 1 exp"
                                    elif nExp == 2: title += ", 2 exp"
                                    elif nExp == 3: title += ", 3 exp"
                                    elif nExp == 4: title += ", 4 exp"
                                if old_method:
                                    title += ", full fit (" + str(nLandau) + " land " + str(nExp) + " exp)"

                                # Legend creation
                                if old_method: legend1 = ROOT.TLegend(0.35, 0.78, 0.6, 0.88)
                                else: legend1 = ROOT.TLegend(0.62, 0.27, 0.9, 0.37)
                                legend1.AddEntry(h_egamma_loose, "Loose Photon, " + str(h_egamma_loose.GetEntries()), "l")
                                #if not args.ratio: legend1.AddEntry(0, "Chi2/NDF: " + str(chi2 / ndf), "")
                                legend2 = ROOT.TLegend(0.65, 0.25, 0.9, 0.35)
                                legend2.AddEntry(h_egamma_tight, "Tight Photon, " + str(h_egamma_tight.GetEntries()), "l")
                                legend2.AddEntry(tight_fit_as_hist, "Fitted Tight (degree "+str(func_with_poly.GetNpar()-1)+")", "l")
                                legend2.AddEntry(tight_fit_w_constant, "Constant fit, C = {:.5}".format(fitted_func_times_constant.GetParameter(0)), "l")
                                
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
                                if not(left == 0 and right == 50):
                                    try:
                                        stats1.Draw()
                                        stats2.Draw()
                                    except NameError:
                                        pass
                                ROOT.gPad.SetLogy()
                                if bins[i] < 60: h_egamma_loose.GetXaxis().SetRangeUser(0, 5)
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
                                        if FTEST: tight_stat.Draw()
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
                                        if bins[i] < 60: h_egamma_tight.GetXaxis().SetRangeUser(0, 5)
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
                                    h_loose_pull.GetYaxis().SetRangeUser(-10, 10)
                                    h_loose_pull.SetStats(0)
                                    if bins[i] < 60: h_loose_pull.GetXaxis().SetRangeUser(0, 5)
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
                                        h_tight_pull.GetYaxis().SetRangeUser(-10, 10)
                                        h_tight_pull.SetStats(0)
                                        if bins[i] < 60: h_tight_pull.GetXaxis().SetRangeUser(0, 5)
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
                            
c1.Print(args.name + ".pdf]")
infile1.Close()

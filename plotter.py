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
parser.add_argument("--integral", default=False, action="store_true", help="add I to tight fit")
parser.add_argument("--printFtest", "--printftest", default=False, action="store_true", help="for a fixed test bin, create a pdf of all possible ftest fits")

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

if args.testBin is None: test_regions = ["iso_asym"] 
elif "noniso_asym" in args.testBin: test_regions = ["noniso_asym"]
elif "noniso_sym" in args.testBin: test_regions = ["noniso_sym"]
elif "iso_asym" in args.testBin: test_regions = ["iso_asym"]
elif "iso_sym" in args.testBin: test_regions = ["iso_sym"]

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
            if args.printFtest and args.testBin is None:
                print("EMPTY PDF: Must have --testBin option when using --printFtest")
                break
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
                    
                    # Set Poisson errors for tight histogram
                    h_egamma_tight.SetBinErrorOption(ROOT.TH1.kPoisson)

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
                        if i == 0: print("====================== " + region.upper() + " =====================")
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
                        
                        old_method = False
                        landau_guess = None
                        exp_guess = None
                        
                        nLandau = 1
                        nExp = 3
                        
                        # Rarely any "bulk" region, i.e. few bins with more than 1000 entries
                        if region == "iso_sym" or region == "iso_asym": 
                            old_method = True

                        # Iso-Sym Guesses
                        if region == "iso_sym":
                            if eta_reg == "barrel":
                                if bins[i] == 20:
                                    nExp -= 1
                                    guesses = [1313, 0.7032, 0.1068, -5.47, -10, 1.014, 0.6314]
                                if bins[i] == 40:
                                    nExp -= 2
                                    guesses = [1523, 0.9313, 0.1856, -3.766, 1.024]
                                if bins[i] == 60:
                                    nExp -= 2
                                    guesses = [785.9, 1.117, 0.2567, -3.498, 1.184]
                                if bins[i] == 70:
                                    nExp -= 1
                                    guesses = [759.8, 1.11, 0.2469, -0.6545, -3.363, 1.025, 0.3268]
                                if bins[i] == 80:
                                    nExp -= 1
                                    guesses = [1648, 1.205, 0.2847, -1.647, -3.407, 1.275, 0.321]
                                if bins[i] == 100:
                                    guesses = [1703, 1.241, 0.2967, -1.365, -2.817, -0.5, 1.375, 0.3328, 1.4] 
                                if bins[i] == 120:
                                    nExp -= 1
                                    guesses = [1776, 1.28, 0.3121, -1.098, -2.588, 1.425, 0.3672] 
                                if bins[i] == 140:
                                    nExp -= 1
                                    guesses = [1887, 1.314, 0.3242, -1.142, -2.368, 1.674, 0.2515]
                                if bins[i] == 160:
                                    nExp -= 1
                                    guesses = [2088, 1.366, 0.3423, -1.125, -2.096, 1.625, 0.4]
                                if bins[i] == 180:
                                    guesses = [2380, 1.361, 0.3352, -1.189, -2.023, -0.05, 1.775, 0.4189, 1.5]
                                if bins[i] == 200:
                                    guesses = [5458, 1.394, 0.3427, -1.012, -1.784, -0.5, 2, 2, 1]
                                if bins[i] == 240:
                                    old_method = False
                                    landau_guess = [5458, 1.394, 0.3427]
                                    exp_guess = [5458, -1.012, 1.678, -1.784, 0.5775, -1]
                                if bins[i] == 300:
                                    guesses = [4131, 1.451, 0.3552, -0.0001324, -1.054, -1.36, 1.369, 0.3256, 1.168]
                                if bins[i] == 380:
                                    guesses = [1559, 1.54, 0.3837, -0.08148, -0.882, -1.123, 1.591, 0.2366, 0.337]
                                if bins[i] == 460:
                                    guesses = [1143, 1.676, 0.4336, -0.3109, -0.9949, -0.5755, 1.561, 0.574, 4.59]
                            elif eta_reg == "endcap":
                                if bins[i] == 20:
                                    nExp -= 1
                                    guesses = [2198, 0.6681, 0.09261, -5.133, -10, 0.854, 0.721]
                                if bins[i] == 40:
                                    nExp -= 2 
                                    guesses = [1944, 0.8865, 0.1719, -3.938, 1.001]
                                if bins[i] == 60:
                                    nExp -= 1 
                                    guesses = [961.9, 1.039, 0.219, -0.06295, -3.632, 0.925, 0.2]
                                if bins[i] == 70:
                                    nExp -= 2
                                    guesses = [885.3, 1.086, 0.2314, -3.771, 1.3]
                                if bins[i] == 80:
                                    nExp -= 1
                                    guesses = [1590, 1.129, 0.2523, -1.411, -3.112, 1.251, 0.2239]
                                if bins[i] == 100:
                                    nExp -= 1
                                    guesses = [1480, 1.228, 0.2869, -1.143, -2.683, 1.387, 0.2375]
                                if bins[i] == 120:
                                    nExp -= 1
                                    guesses = [1444, 1.277, 0.3012, -0.837, -2.624, 1.56, 0.2]
                                if bins[i] == 140:
                                    nExp -= 1
                                    guesses = [1424, 1.296, 0.3049, -0.458, -2.143, 1.348, 0.4153]
                                if bins[i] == 160:
                                    nExp -= 1
                                    guesses = [1480, 1.321, 0.3095, -0.5116, -2.012, 1.375, 0.513]
                                if bins[i] == 180:
                                    nExp -= 1
                                    guesses = [1559, 1.436, 0.3578, -0.8864, -2.058, 1.525, 0.6868]
                                if bins[i] == 200:
                                    nExp -= 1
                                    guesses = [3150, 1.473, 0.3669, -1.522, -0.5, 2, 2.5]
                                if bins[i] == 240:
                                    guesses = [3040, 1.432, 0.3386, -0.1704, -1.026, -1.618, 1.405, 0.4994, 0.6806]
                                if bins[i] == 300:
                                    guesses = [1596, 1.523, 0.3731, -0.07532, -0.9641, -1.38, 1.405, 0.4691, 1.187]
                                if bins[i] == 380:
                                    nExp -= 1
                                    guesses = [475.8, 1.65, 0.4243, -0.5718, -1.198, 1.825, 0.8539]
                                if bins[i] == 460:
                                    guesses = [195.5, 1.389, 0.2903, -0.0001153, -0.7601, -0.9159, 1.305, 0.8198, 0.2123]
                        
                        # Iso-Asym Guesses
                        if region == "iso_asym":
                            if eta_reg == "barrel":
                                if bins[i] == 20:
                                    nExp -= 1
                                    guesses = [246.5, 0.8312, 0.1213, -4.7, -10, 0.875, 0.4371]
                                if bins[i] == 40:
                                    nExp -= 2
                                    guesses = [422.5, 0.9891, 0.1768, -4.297, 1.116]
                                if bins[i] == 60:
                                    nExp -= 1
                                    guesses = [224.5, 1.061, 0.193, -2.913, -5.832, 1.141, 1.184]
                                if bins[i] == 70:
                                    nExp -= 1
                                    guesses = [266.6, 1.228, 0.266, -1.351, -5.152, 1.029, 0.8208]
                                if bins[i] == 80:
                                    nExp -= 1
                                    guesses = [532.3, 1.19, 0.2404, -0.1541, -2.411, 1.075, 0.2]
                                if bins[i] == 100:
                                    nExp -= 1
                                    guesses = [585.1, 1.328, 0.2896, -1.456, -2.425, 1.362, 0.2857]
                                if bins[i] == 120:
                                    guesses = [683.6, 1.396, 0.3051, -1.74, -2.069, -4.83, 1.462, 0.4043, 1.386]
                                if bins[i] == 140:
                                    guesses = [772.1, 1.426, 0.3213, -0.09985, -2.09, -10, 1.379, 0.3031, 2.843]
                                if bins[i] == 160:
                                    nExp -= 1
                                    guesses = [955.5, 1.547, 0.3559, -1.125, -2.009, 1.725, 0.25]
                                if bins[i] == 180:
                                    nExp -= 1
                                    guesses = [1236, 1.713, 0.4136, -1.181, -1.847, 1.674, 0.3422]
                                if bins[i] == 200:
                                    guesses = [3034, 1.763, 0.431, -0.7333, -1.646, -1.62, 1.813, 0.2465, 1.083]
                                if bins[i] == 240:
                                    guesses = [3956, 1.879, 0.4745, -1.248, -1.578, -0.7454, 2.209, 0.2557, 3.601]
                                if bins[i] == 300:
                                    guesses = [1.706e+04, 1.529, 0.3766, -0.01528, -0.7213, -1.017, 4.5, 4, 2] 
                                if bins[i] == 380:
                                    guesses = [1096, 2.116, 0.5494, -0.4469, -1.293, -0.8872, 2.428, 0.4568, 2.857]
                                if bins[i] == 460:
                                    guesses = [793.4, 2.157, 0.5492, -0.01715, -0.9391, -0.5238, 2.047, 0.8956, 5.37]
                            elif eta_reg == "endcap":
                                if bins[i] == 20:
                                    nExp -= 1
                                    guesses = [299.9, 0.7635, 0.1048, -7.313, -10, 0.9439, 0.2901]
                                if bins[i] == 40:
                                    nExp -= 2
                                    guesses = [546.6, 0.8803, 0.1328, -4.339, 1.089]
                                if bins[i] == 60:
                                    nExp -= 2
                                    guesses = [346.3, 1.128, 0.2239, -3.488, 0.9861]
                                if bins[i] == 70:
                                    nExp -= 1
                                    guesses = [290.3, 1.008, 0.1734, -3.002, -9.928, 1.331, 1.363]
                                if bins[i] == 80:
                                    guesses = [787.1, 1.41, 0.3193, -0.9958, -2.985, -2.563, 1.049, 0.3013, 0.9245]
                                if bins[i] == 100:
                                    nExp -= 1
                                    guesses = [678.8, 1.319, 0.2715, -2.841e-10, -2.461, 1.08, 0.2679]
                                if bins[i] == 120:
                                    guesses = [646.3, 1.374, 0.2959, -1.153, -2.426, -9.998, 1.394, 0.2735, 2.457]
                                if bins[i] == 140:
                                    nExp -= 1
                                    guesses = [713.7, 1.365, 0.2823, -2.554e-11, -2.152, 1.299, 0.3368]
                                if bins[i] == 160:
                                    nExp -= 1
                                    guesses = [851.9, 1.539, 0.3443, -1.333, -2.108, 1.625, 0.391]
                                if bins[i] == 180:
                                    guesses = [1041, 1.725, 0.4141, -1.863, -1.999, -0.9901, 1.807, 0.3725, 2.195]
                                if bins[i] == 200:
                                    guesses = [2319, 1.77, 0.4112, -1.939, -1.37, -1.744, 2.033, 1.962, 1.919]
                                if bins[i] == 240:
                                    guesses = [2679, 2.01, 0.505, -1.64, -1.099, -1.108, 2.146, 2.514, 2.416]
                                if bins[i] == 300:
                                    guesses = [1606, 1.895, 0.4265, -0.7, -1.483, -1, 2, 0.5, 2.5]
                                if bins[i] == 380:
                                    guesses = [500, 2.301, 0.6017, -1.264, -1.592, -0.8502, 2.804, 0.971, 0.9496]
                                if bins[i] == 460:
                                    nExp -= 1
                                    guesses = [270.3, 2.457, 0.6263, -0.9771, -0.5, 3.075, 5]

                        # NonIso-Sym Guesses
                        if region == "noniso_sym":
                            if eta_reg == "barrel":
                                if bins[i] == 20:
                                    landau_guess = [3.491e+06, 0.8838, 0.145] 
                                    exp_guess = [8.593e+05, -4.278, 1.595, -7.263, 0.3, -8]
                                if bins[i] == 40:
                                    landau_guess = [1.705e+04, 1.33, 0.205]
                                    exp_guess = [3.636e+05, -2.809, 2.176, -4.8, 0.6, -1.807]
                                if bins[i] == 60:
                                    nExp -= 1
                                    landau_guess = [2.212e+04, 0.6783, 0.09052]
                                    exp_guess = [4.467e+05, -2.789, 2.974, -5.139]
                                if bins[i] == 70:
                                    landau_guess = [1.857e+04, 0.682, 0.09099]
                                    exp_guess = [1.349e+04, -1.075, 2.269, -2.444, 1.235, -5.598]
                                if bins[i] == 80:
                                    landau_guess = [1.826e+04, 0.6239, 0.07162]
                                    exp_guess = [1.506e+04, -0.9609, 2.4, -2.24, 2.2, -0.4279]
                                if bins[i] == 100:
                                    landau_guess = [1.705e+04, 0.6496, 0.0804] 
                                    exp_guess = [1.755e+04, -0.92, 3.125, -1.706, 2, -1]
                                if bins[i] == 120:
                                    landau_guess = [1.705e+04, 0.6496, 0.3] 
                                    exp_guess = [1.755e+04, -0.92, 3.125, -1.706, 2, -1]
                                if bins[i] == 140:
                                    landau_guess = [8972, 0.8638, 0.1507]
                                    exp_guess = [9448, -1.041, 3.2, -1.423, 3, -3]
                                if bins[i] == 160:
                                    landau_guess = [9535, 0.6861, 0.08948]
                                    exp_guess = [5.971e+04, -1.224, 4, -3.118, 4, -4.477]
                                if bins[i] == 180:
                                    landau_guess = [1.183e+04, 0.7977, 0.1293]
                                    exp_guess = [4.739e+04, -1.172, 3.5, -3.605, 5.5, -3.677]
                                if bins[i] == 200:
                                    nExp -= 1
                                    landau_guess = [1.204e+04, 0.7742, 0.1196]
                                    exp_guess = [2.194e+05, -1.407, 5.087, -0.8335]
                                if bins[i] == 240:
                                    landau_guess = [8717, 0.8303, 0.1381]
                                    exp_guess = [9.701e+04, -1.193, 6.173, -0.6467, 6.452, -2.239]
                                if bins[i] == 300:
                                    landau_guess = [1145, 1.149, 0.2453]
                                    exp_guess = [3.14e+04, -1.023, 6.75, -0.5922, 1.25, -1]
                                if bins[i] == 380:
                                    old_method = True
                                    guesses = [4848, 1.63, 0.4089, -3.553e-13, -0.6753, -0.8161, 0.5, 1.404, 0.5423]
                                if bins[i] == 460:
                                    old_method = True
                                    guesses = [2487, 1.731, 0.4342, -0.08184, -0.6556, -0.6027, 0.5, 1.187, 0.5]
                            elif eta_reg == "endcap":
                                if bins[i] == 20:
                                    landau_guess = [1.705e+04, 1.108, 0.1949] 
                                    exp_guess = [4182, -0.9553, 1.4, -5.89, 0.5, -1]
                                if bins[i] == 40:
                                    landau_guess = [1.123e+05, 0.875, 0.1524]
                                    exp_guess = [4.549e+04, -2.09, 1.915, -3.572, 0.6484, -6.328]
                                if bins[i] == 60:
                                    landau_guess = [1.021e+04, 0.7237, 0.107]
                                    exp_guess = [5970, -0.9294, 1.702, -2.435, 1.262, -5.012]
                                if bins[i] == 70:
                                    landau_guess = [1.163e+04, 0.7851, 0.1258]
                                    exp_guess = [1.045e+04, -1.111, 1.645, -2.263, 1.974, -6.166]
                                if bins[i] == 80:
                                    landau_guess = [1.826e+04, 0.6239, 0.07162]
                                    exp_guess = [1.506e+04, -0.9609, 2.4, -2.24, 2.2, -0.4279]
                                if bins[i] == 100:
                                    landau_guess = [1.705e+04, 0.6496, 0.2]
                                    exp_guess = [1.755e+04, -0.92, 2.4, -1.706, 2.6, -1]
                                if bins[i] == 120:
                                    landau_guess = [1.705e+04, 0.6496, 0.3] 
                                    exp_guess = [1.755e+04, -0.92, 3.125, -1.706, 2, -1]
                                if bins[i] == 140:
                                    landau_guess = [1.761e+04, 0.7322, 0.1067]
                                    exp_guess = [9448, -1.041, 3.2, -1.423, 3, -3]
                                if bins[i] == 160:
                                    landau_guess = [7468, 0.9087, 0.1663]
                                    exp_guess = [1.107e+04, -1.035, 3, -1.572, 2, -7.523]
                                if bins[i] == 180:
                                    landau_guess = [7031, 1.004, 0.1965]
                                    exp_guess = [1.742e+04, -1.198, 3.5, -3.296, 2.5, -1.366]
                                if bins[i] == 200:
                                    landau_guess = [1.204e+04, 0.7742, 0.1196]
                                    exp_guess = [4.491e+04, -1.059, 5.5, -3.574, 4.5, -2.234]
                                if bins[i] == 240:
                                    landau_guess = [7937, 1.269, 0.2811]
                                    exp_guess = [4.457e+04, -2, 5, -1, 1.8, -0.5]
                                if bins[i] == 300:
                                    old_method = True
                                    guesses = [1.706e+04, 1.529, 0.3766, -0.01528, -0.7213, -1.017, 0.5, 0.8, 0.5]  
                                if bins[i] == 380:
                                    old_method = True
                                    guesses = [890.4, 1.723, 0.4403, -1.441, -0.9773, -0.2518, 2.5, 2.5, 5]
                                if bins[i] == 460:
                                    old_method = True
                                    guesses = [301.3, 1.842, 0.482, -1.031, -0.5918, -1.591, 0.5, 3.275, 0.8257]
                               
                        # NonIso_Asym Guesses   
                        if region == "noniso_asym":
                            if eta_reg == "barrel":
                                if bins[i] == 20:
                                    landau_guess = [9.069e+04, 0.9345, 0.1463]
                                    exp_guess = [4.834e+04, -3.378, 1.206, -8.486, 0.7373, -4.555]
                                if bins[i] == 40:
                                    nExp -= 1
                                    landau_guess = [2.459e+04, 0.8719, 0.1261]
                                    exp_guess = [9328, -1.29, 1.599, -5.368]
                                if bins[i] == 60:
                                    landau_guess = [1.021e+04, 0.7237, 0.107]
                                    exp_guess = [5970, -0.9294, 1.702, -2.435, 1.262, -5.012]
                                if bins[i] == 70:
                                    nExp -= 1
                                    landau_guess = [1.262e+04, 1.14, 0.2139]
                                    exp_guess = [7.614e+04, -2.433, 2.075, -3.571]
                                if bins[i] == 80:
                                    landau_guess = [1.474e+04, 1.084, 0.1992]
                                    exp_guess = [7428, -0.8801, 2.325, -2.899, 1.044, -2.326]
                                if bins[i] == 100:
                                    landau_guess = [1.705e+04, 0.6496, 0.0804] 
                                    exp_guess = [1.986e+04, -0.8552, 2.8, -2.188, 1.2, -1.68]
                                if bins[i] == 120:
                                    landau_guess = [1.705e+04, 0.6496, 0.3] 
                                    exp_guess = [1.755e+04, -0.92, 3.125, -1.706, 2, -1]
                                if bins[i] == 140:
                                    landau_guess = [8211, 1.16, 0.2245] 
                                    exp_guess = [3.881e+04, -1.409, 3.225, -1.699, 2.5, -1.078]
                                if bins[i] == 160:
                                    landau_guess = [8663, 1.297, 0.2693]
                                    exp_guess = [2.436e+04, -1.231, 2.659, -1.467, 3.5, -1]
                                if bins[i] == 180:
                                    landau_guess = [1.183e+04, 0.7977, 0.1293]
                                    exp_guess = [4.739e+04, -1.172, 3.5, -3.605, 5.5, -3.677]
                                if bins[i] == 200:
                                    landau_guess = [1.204e+04, 0.7742, 0.1196]
                                    exp_guess = [4.491e+04, -1.059, 3, -1, 3, -2.234]
                                if bins[i] == 240:
                                    landau_guess = [8442, 1.149, 0.2453]
                                    exp_guess = [3.14e+04, -1.023, 5, -0.5922, 3, -2]
                                if bins[i] == 300:
                                    nExp -= 1
                                    landau_guess = [5089, 2.157, 0.5678]
                                    exp_guess = [1.148e+04, -1.051, 6.125, -0.5868]
                                if bins[i] == 380:
                                    old_method = True
                                    guesses = [1430, 2.107, 0.5122, -0.04034, -1.359, -0.5, 3, 3, 2]
                                if bins[i] == 460:
                                    nExp -= 1
                                    old_method = True
                                    guesses = [775.1, 2.767, 0.7776, -0.7585, -0.3751, 3.079, 5.396]
                            elif eta_reg == "endcap":
                                if bins[i] == 20:
                                    nExp -= 1
                                    landau_guess = [6045, 0.7751, 0.1028]
                                    exp_guess = [7331, -3.713, 1.18, -7.362]
                                if bins[i] == 40:
                                    landau_guess = [1e+04, 0.9233, 0.1527]
                                    exp_guess = [3937, -0.4774, 0.9104, -2.799, 0.7293, -4.712]
                                if bins[i] == 60:
                                    landau_guess = [1.021e+04, 0.7237, 0.107]
                                    exp_guess = [5970, -0.9294, 1.702, -2.435, 1.262, -5.012]
                                if bins[i] == 70:
                                    nExp -= 1
                                    landau_guess = [5725, 1.131, 0.2179]
                                    exp_guess = [9440, -1.938, 1.9, -2.955]
                                if bins[i] == 80:
                                    landau_guess = [1.826e+04, 0.6239, 0.07162]
                                    exp_guess = [1.506e+04, -0.9609, 2.4, -2.24, 2.2, -0.4279]
                                if bins[i] == 100:
                                    landau_guess = [1.705e+04, 0.6496, 0.2]
                                    exp_guess = [1.755e+04, -0.92, 2.4, -1.706, 2.6, -1]
                                if bins[i] == 120:
                                    landau_guess = [1.705e+04, 0.6496, 0.3]
                                    exp_guess = [1.755e+04, -0.92, 3.125, -1.706, 2, -1]
                                if bins[i] == 140:
                                    landau_guess = [6416, 1.397, 0.3019]
                                    exp_guess = [8019, -1.101, 1.875, -1.559, 3.2, -1.032]
                                if bins[i] == 160:
                                    landau_guess = [5895, 1.577, 0.3635]
                                    exp_guess = [8294, -1.077, 1.691, -1.543, 1.484, -1.361]
                                if bins[i] == 180:
                                    landau_guess = [7031, 1.004, 0.1965]
                                    exp_guess = [1.742e+04, -1.198, 3, -3.296, 3, -1.366]
                                if bins[i] == 200:
                                    nExp -= 1
                                    landau_guess = [5881, 1.694, 0.3963]
                                    exp_guess = [2.347e+04, -1.451, 3.478, -1.168]
                                if bins[i] == 240:
                                    nExp -= 1
                                    old_method = True
                                    guesses = [4657, 2.046, 0.5207, -1.225, -0.7, 2.5, 3.5]
                                if bins[i] == 300:
                                    old_method = True
                                    guesses = [1470, 1.687, 0.3483, -2, -1, -1, 3, 2, 1.5]
                                if bins[i] == 380:
                                    old_method = True
                                    guesses = [1430, 2.107, 0.5122, -0.04034, -1.359, -0.5, 3, 3, 2]
                                if bins[i] == 460:
                                    old_method = True
                                    guesses = [345.4, 2.06, 0.485, -2, -1, -0.5, 3, 2, 1]
                      
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
                        fitted_func_times_constant, _, _ = util.MultiplyWithPolyToTF1(fitted_func, 0, poly=0)
                        h_egamma_tight.Fit(fitted_func_times_constant, '0L' if not args.integral else '0LI')
                        tight_fit_w_constant = util.TemplateToHistogram(fitted_func_times_constant, 1000, 0, 50)

                        FTEST = True if args.ftest else False
                        NUM_PLOTS = 1
                        if not FTEST:
                            POLY_TYPE = 3
                            DEGREE = 1
                            func_with_poly, _, _ = util.MultiplyWithPolyToTF1(fitted_func, DEGREE, poly=POLY_TYPE)
                            h_egamma_tight.Fit(func_with_poly, '0L' if not args.integral else '0LI')
                            tight_fit_as_hist = util.TemplateToHistogram(func_with_poly, 1000, 0, 50)
                        
                        if FTEST:
                            parse = args.ftest.split()
                            NUM_DEGREES = int(parse[1])
                            if args.printFtest and args.testBin is not None: NUM_PLOTS = NUM_DEGREES+1
                            POLY_TYPE = int(parse[0])
                            fitfuncs = []
                            fitresults = []
                            statboxes = []
                            for degree in range(NUM_DEGREES+1):
                                func_with_poly, _, _ = util.MultiplyWithPolyToTF1(fitted_func, degree, poly=POLY_TYPE)
                                fitresult = h_egamma_tight.Fit(func_with_poly, '0SL' if not args.integral else '0SLI')
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
                                    decision, ftest, target, rss1, rss2, dof1, dof2 = util.ftest(
                                        h_egamma_tight, fitfuncs[d2], fitresults[d2], fitfuncs[d1], fitresults[d1])
                                    print(d2, '>', d1, decision)
                                    print('  F={} target={}'.format(ftest, target), '({}, {}) dof'.format(dof1, dof2))
                                    if decision: best_d = d2
                            print('Best: ', best_d)
                            func_with_poly = fitfuncs[best_d]
                            tight_fit_as_hist = util.TemplateToHistogram(func_with_poly, 1000, 0, 50)
                            tight_stat = statboxes[best_d]
                        
                        for plot in range(NUM_PLOTS):
                            if args.printFtest and args.testBin:
                                func_with_poly = fitfuncs[plot]
                                tight_fit_as_hist = util.TemplateToHistogram(func_with_poly, 1000, 0, 50)
                                tight_stat = statboxes[plot]

                            just_poly = util.ExtractPolyFromTightFit(func_with_poly, poly=POLY_TYPE)
                            tlast_bin = h_egamma_tight.GetNbinsX() + 1
                            for b in reversed(range(h_egamma_tight.GetNbinsX())):
                              if h_egamma_tight.GetBinContent(b+1)==0: continue
                              else:
                                tlast_bin = b + 1
                                break
                            rightmost_tightdata = h_egamma_loose.GetBinLowEdge(tlast_bin+1)
                            just_poly.SetRange(0,rightmost_tightdata)

                            # determine bin-by-bin error
                            STEP_SIZE = 0.001
                            integral = False
                            hist = h_egamma_tight
                            fit = func_with_poly
                            #ndof = util.count_nonzero_bins(hist) - fit.GetNpar()
                            ndof = tlast_bin - fit.GetNpar()
                            chi2, by_bin = util.RSS(fit, hist, error=0, integral=integral, chi2=True)
                            chi2_ndof = chi2 / ndof
                            error = 0.00
                            #print('initial')
                            #print('chi2', chi2, 'ndof', ndof, 'chi2_ndof', chi2_ndof)
                            while chi2_ndof > 1.0:
                                error += STEP_SIZE
                                #print('trying', error)
                                chi2, _ = util.RSS(fit, hist, error=error, integral=integral, chi2=True)
                                chi2_ndof = chi2 / ndof
                                #print('  ', 'chi2', chi2, 'chi2/ndof', chi2_ndof)
                            bin_bin_error = error

                            '''
                            print('example:')
                            chi2 = 0
                            for j in range(h_egamma_tight.GetNbinsX()):
                                if j > tlast_bin: continue
                                edge = h_egamma_tight.GetBinLowEdge(j+1)
                                hist = h_egamma_tight.GetBinContent(j+1)
                                fit = tight_fit_as_hist.GetBinContent(j+1)
                                if h_egamma_tight.GetBinContent(j+1) == 0:
                                    err = h_egamma_tight.GetBinErrorUp(j+1)
                                else: 
                                    if tight_fit_as_hist.GetBinContent(j+1) > h_egamma_tight.GetBinContent(j+1):
                                        err = h_egamma_tight.GetBinErrorUp(j+1)
                                    else: err = h_egamma_tight.GetBinErrorLow(j+1)
                                pull = (hist - fit)/err
                                pull2 = pull**2
                                chi2 += pull2
                                print('{:<2} bin {:<4} hist {:<8.4} fit {:<12.4} err {:<8.4} pull {:<12.4} pull^2 {:<12.4}'.format(j, edge, hist, fit, err, pull, pull2))
                                #print('{:<2} bin {:<4} hist {:<8.4} fit {:<12.4} err {:<8.4} pull {:<12.4} pull^2 {:<12.4} func {}'.format(j, edge, hist, fit, err, pull, pull2, by_bin[j][4]))
                            print('chi2', chi2, 'ndof', ndof, 'chi2/ndof', chi2/ndof)
                            '''

                            h_loose_pull_num = h_egamma_loose.Clone()
                            h_loose_pull_num.Reset()
                            h_loose_pull = h_egamma_loose.Clone()
                            h_loose_pull.Reset()
                            h_loose_pull_num.Add(h_egamma_loose, loose_fit_as_hist, 1, -1)  # Numerator of pull hist is data - fit

                            for j in range(h_loose_pull_num.GetNbinsX()): 
                                if h_egamma_loose.GetBinContent(j+1) == 0: err = 1.8
                                else: err = h_egamma_loose.GetBinError(j+1)
                                h_loose_pull.SetBinContent(j+1, h_loose_pull_num.GetBinContent(j+1)/err)
                                h_loose_pull.SetBinError(j+1, 1)
                            
                            h_tight_pull_num = h_egamma_tight.Clone()
                            h_tight_pull_num.Reset()
                            h_tight_pull = h_egamma_tight.Clone()
                            h_tight_pull.Reset()
                            h_tight_pull_num.Add(h_egamma_tight, tight_fit_as_hist, 1, -1)  # Numerator of pull hist is data - fit
                            
                            h_tight_pull_error = h_tight_pull.Clone() # to visualize binbin error
                            h_tight_pull_error.Reset()

                            for j in range(h_tight_pull_num.GetNbinsX()): 
                                if h_egamma_tight.GetBinContent(j+1) == 0:
                                    err = h_egamma_tight.GetBinErrorUp(j+1)
                                else: 
                                    if tight_fit_as_hist.GetBinContent(j+1) > h_egamma_tight.GetBinContent(j+1):
                                        err = h_egamma_tight.GetBinErrorUp(j+1)
                                    else: err = h_egamma_tight.GetBinErrorLow(j+1)
                                h_tight_pull.SetBinContent(j+1, h_tight_pull_num.GetBinContent(j+1)/err)
                                h_tight_pull.SetBinError(j+1, 1)
                                h_tight_pull_error.SetBinContent(j+1, 0)
                                h_tight_pull_error.SetBinError(j+1, (tight_fit_as_hist.GetBinContent(j+1)*bin_bin_error)/err)

                            # pull for tight fit with constant
                            h_tight_pullc_num = h_egamma_tight.Clone()
                            h_tight_pullc_num.Reset()
                            h_tight_pullc = h_egamma_tight.Clone()
                            h_tight_pullc.Reset()
                            h_tight_pullc_num.Add(h_egamma_tight, tight_fit_w_constant, 1, -1)  # Numerator of pull hist is data - fit
                            
                            for j in range(h_tight_pullc_num.GetNbinsX()): 
                                if h_egamma_tight.GetBinContent(j+1) == 0: err = 1.8
                                else: 
                                    if tight_fit_as_hist.GetBinContent(j+1) > h_egamma_tight.GetBinContent(j+1): err = h_egamma_tight.GetBinErrorUp(j+1)
                                    else: err = h_egamma_tight.GetBinErrorLow(j+1)
                                h_tight_pullc.SetBinContent(j+1, h_tight_pullc_num.GetBinContent(j+1)/err)
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
                            legend2 = ROOT.TLegend(0.29, 0.70, 0.62, 0.89)
                            legend2.AddEntry(h_egamma_tight, "Tight Photon, " + str(h_egamma_tight.GetEntries()), "l")
                            if FTEST: legend2.AddEntry(tight_fit_as_hist, "Fit w f-test (Degree "+str(func_with_poly.GetNpar()-1)+")", "l")
                            else: legend2.AddEntry(tight_fit_as_hist, "Fit (Degree "+str(func_with_poly.GetNpar()-1)+")", "l")
                            legend2.AddEntry(tight_fit_w_constant, "Constant fit, C = {:.4}".format(fitted_func_times_constant.GetParameter(0)), "l")
                            legend2.AddEntry('', 'Chi2/Ndof: {:.3f}'.format(chi2_ndof), '')
                            legend2.AddEntry('', 'Bin Error: {:.1%}'.format(bin_bin_error), '')
                            
                            # Draw plots
                            if args.fit: c1.cd(1)
                            else:
                              c1.cd()
                              pad1 = ROOT.TPad('pad1', 'pad1', 0, 0.3, 0.5, 1)
                              pad1.Draw()
                              pad1.cd()
                            
                            h_egamma_loose.SetTitle(title)
                            h_egamma_loose.SetMaximum()
                            h_egamma_loose.SetMinimum(0.1)
                            h_egamma_loose.Draw("e")
                            loose_fit_as_hist.SetLineColor(ROOT.kRed+1)
                            loose_fit_as_hist.Draw("same")
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
                                    ROOT.gPad.SetLogy()
                                    if bins[i] < 60: h_egamma_tight.GetXaxis().SetRangeUser(0, 5)
                                    elif bins[i] < 120: h_egamma_tight.GetXaxis().SetRangeUser(0, 10)
                                    elif bins[i] < 200: h_egamma_tight.GetXaxis().SetRangeUser(0, 15)
                                    elif bins[i] < 380: h_egamma_tight.GetXaxis().SetRangeUser(0, 20)
                                    else: h_egamma_tight.GetXaxis().SetRangeUser(0, 26)
                                    h_egamma_tight.SetMinimum(0.1)
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
                                    ROOT.gPad.Update()
                                    
                                    legend2.Draw("same")
                                    overlay = ROOT.TPad("overlay","",0, 0.06, 1, 0.5)
                                    overlay.SetFillStyle(4000)
                                    overlay.SetFillColor(0)
                                    overlay.SetFrameFillStyle(4000)
                                    overlay.SetFrameLineWidth(0)
                                    overlay.Draw()
                                    overlay.cd()
                                    empty = ROOT.TH1F(util.getname('empty'), '', 100, 0, 50)
                                    empty.SetLineColor(ROOT.kRed)
                                    if bins[i] < 60: empty.GetXaxis().SetRangeUser(0, 5)
                                    elif bins[i] < 120: empty.GetXaxis().SetRangeUser(0, 10)
                                    elif bins[i] < 200: empty.GetXaxis().SetRangeUser(0, 15)
                                    elif bins[i] < 380: empty.GetXaxis().SetRangeUser(0, 20)
                                    else: empty.GetXaxis().SetRangeUser(0, 26)
                                    empty.GetYaxis().SetRangeUser(just_poly.GetMinimum(), just_poly.GetMaximum())
                                    empty.Draw('AH')
                                    just_poly.SetTitle("")
                                    just_poly.Draw("AI L same")
                                    
                                    if not args.fit:
                                        if not region == "iso_sym":
                                            c1.cd()
                                            pad2 = ROOT.TPad('pad2', 'pad2', 0.5, 0.3, 1, 1)
                                            pad2.Draw()
                                            pad2.cd()
                                            ROOT.gPad.SetLogy()
                                            if bins[i] < 60: h_egamma_tight.GetXaxis().SetRangeUser(0, 5)
                                            elif bins[i] < 120: h_egamma_tight.GetXaxis().SetRangeUser(0, 10)
                                            elif bins[i] < 200: h_egamma_tight.GetXaxis().SetRangeUser(0, 15)
                                            elif bins[i] < 380: h_egamma_tight.GetXaxis().SetRangeUser(0, 20)
                                            else: h_egamma_tight.GetXaxis().SetRangeUser(0, 26)
                                            h_egamma_tight.SetMinimum(0.1)
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
                                            ROOT.gPad.Update()
                                            legend2.Draw("same")
                                            overlay = ROOT.TPad("overlay","",0, 0.06, 1, 0.5)
                                            overlay.SetFillStyle(4000)
                                            overlay.SetFillColor(0)
                                            overlay.SetFrameFillStyle(4000)
                                            overlay.SetFrameLineWidth(0)
                                            overlay.Draw()
                                            overlay.cd()
                                            empty = ROOT.TH1F(util.getname('empty'), '', 100, 0, 50)
                                            empty.SetLineColor(ROOT.kRed)
                                            if bins[i] < 60: empty.GetXaxis().SetRangeUser(0, 5)
                                            elif bins[i] < 120: empty.GetXaxis().SetRangeUser(0, 10)
                                            elif bins[i] < 200: empty.GetXaxis().SetRangeUser(0, 15)
                                            elif bins[i] < 380: empty.GetXaxis().SetRangeUser(0, 20)
                                            else: empty.GetXaxis().SetRangeUser(0, 26)
                                            empty.GetYaxis().SetRangeUser(min(0, just_poly.GetMinimum()), just_poly.GetMaximum())
                                            empty.Draw('AH')
                                            just_poly.SetTitle("")
                                            just_poly.SetRange(0, 50)
                                            just_poly.Draw("AI L same")
                                            ROOT.gPad.Update()
                                            rightaxis = ROOT.TGaxis(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUymax(), 510, "L+")
                                            rightaxis.SetLineColor(ROOT.kRed);
                                            rightaxis.SetLabelColor(ROOT.kRed);
                                            rightaxis.Draw()
                                            ROOT.gPad.Update()
                                            #topaxis = ROOT.TGaxis(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUxmin(), ROOT.gPad.GetUxmax(), 510, "+L")
                                            #topaxis.SetLineColor(ROOT.kRed);
                                            #topaxis.SetLabelColor(ROOT.kRed);
                                            #topaxis.Draw()
                                    
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
                                            h_tight_pull.Draw('pe')
                                            h_tight_pull_error.SetLineColor(ROOT.kGray+2)
                                            h_tight_pull_error.SetFillColor(ROOT.kGray+2)
                                            h_tight_pull_error.Draw('same e2')
                                            h_tight_pull.SetTitle("(Tight - Fit) / Error")
                                            h_tight_pull.SetLineColor(ROOT.kBlack)
                                            h_tight_pull.Draw('pe same')
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
                                    rightaxis = ROOT.TGaxis(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUymax(), 510, "L+")
                                    rightaxis.SetLineColor(ROOT.kRed);
                                    rightaxis.SetLabelColor(ROOT.kRed);
                                    rightaxis.Draw()
                                    ROOT.gPad.Update()
                                    #topaxis = ROOT.TGaxis(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUxmin(), ROOT.gPad.GetUxmax(), 510, "+L")
                                    #topaxis.SetLineColor(ROOT.kRed);
                                    #topaxis.SetLabelColor(ROOT.kRed);
                                    #topaxis.Draw()
                            
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
                                    h_tight_pull.Draw('pe')
                                    h_tight_pull_error.SetLineColor(ROOT.kGray+2)
                                    h_tight_pull_error.SetFillColor(ROOT.kGray+2)
                                    h_tight_pull_error.Draw('same e2')
                                    h_tight_pull.SetTitle("(Tight - Fit) / Error")
                                    h_tight_pull.SetLineColor(ROOT.kBlack)
                                    h_tight_pull.Draw('pe same')
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
                            if args.testBin is not None and not args.printFtest: input()
                            c1.Print(args.name + ".pdf")
                                
    c1.Print(args.name + ".pdf]")
    infile1.Close()

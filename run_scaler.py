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
import itertools

def binConverter(test_bin):
    bin_list = test_bin.split(" ")
    return bin_list

# "scale" a control-region tight photon distribution to its corresponding signal distribution
def removeEntries(bkg_hist, sig_hist):
    if bkg_hist.Integral() == 0: return False
   
    p = float(sig_hist.Integral()) / bkg_hist.Integral() # probability of removing entry
    if p >= 1: return False
    
    for i in range(1, bkg_hist.GetNbinsX()+1):
        N = round(bkg_hist.GetBinContent(i))
        if N < 100:
            bkg_hist.SetBinContent(i, round(rand.Binomial(N, p)))
        else: 
            content = rand.Gaus(N*p, (N*p*(1-p))**0.5)
            if content < 0: bkg_hist.SetBinContent(i, 0) 
            else: bkg_hist.SetBinContent(i, round(content))

# command line options
parser = argparse.ArgumentParser(description="")
parser.add_argument("input", metavar="INPUT", help="input root file")
parser.add_argument("--testBin", default=None, help="specify bin to test")
# seed?
# cutoff?
args = parser.parse_args()

# Constants
directory = "scaled_tight_hists"
bins = [20,40,60,80,100,140,180,220,300,380]
eta_regions = ["barrel", "endcap"]
regions = ["iso_sym", "iso_asym", "noniso_sym", "noniso_asym"]
rand = ROOT.TRandom3()
rand.SetSeed(12445)

# init
infile1 = ROOT.TFile(args.input)
if args.testBin is None: test_regions = ["noniso_sym"] 
elif "noniso_asym" in args.testBin: test_regions = ["noniso_asym"]
elif "noniso_sym" in args.testBin: test_regions = ["noniso_sym"]
elif "iso_asym" in args.testBin: test_regions = ["iso_asym"]
elif "iso_sym" in args.testBin: test_regions = ["iso_sym"]
if args.testBin is not None: test_bin = binConverter(args.testBin)
if not os.path.exists(directory): os.mkdir(directory)
os.chdir(directory)
    
"""
for region in regions:  # loop through twoprong sideband regions
    if args.testBin is not None: 
        if not region == test_bin[0]: continue
    print("at begin region loop:", region)
    for i in range(len(bins)):  # loop through pt bins for a fixed twoprong sideband
        if args.testBin is not None: 
            if not bins[i] == int(test_bin[2]): continue
        print("  at begin bin loop:", region, bins[i])
        for eta_reg in eta_regions:  # loop through eta regions for fixed pt-bin and fixed twoprong sideband
            if args.testBin is not None: 
                if not eta_reg == test_bin[1]: continue
            if not eta_reg == "barrel" and not eta_reg == "endcap": continue  # no pt-bin plots for barrel and endcap combined, so skip this case
            
            print("    at begin eta loop:", region, bins[i], eta_reg)
"""
scaled_tight_files = []
file_counter_tight = 0
for region, i, eta_reg in itertools.product(regions, range(len(bins)), eta_regions):
    if args.testBin:
        if region != test_bin[0]: continue
        if eta_reg != test_bin[1]: continue
        if bins[i] != int(test_bin[2]): continue

    egamma_tight_plots = "plots/twoprong_masspi0_" + region + "_" + eta_reg
    if i == len(bins)-1: egamma_tight_plots += "_" + str(bins[i]) + "+"
    else: egamma_tight_plots += "_" + str(bins[i]) + "_" + str(bins[i+1])
    egamma_tight_plots += "_tight"

    h_egamma_tight = infile1.Get(egamma_tight_plots)

    # Scale
    if not region == "iso_sym": removeEntries(h_egamma_tight, infile1.Get(egamma_tight_plots.replace(region, "iso_sym")))

    # Save
    if i == len(bins)-1: title = region + "_" + eta_reg + "_" + str(bins[i]) + "+_tight"
    else: title = region + "_" + eta_reg + "_" + str(bins[i]) + "_" + str(bins[i+1]) + "_tight" 
    outfile = ROOT.TFile(title + ".root", "RECREATE")
    outfile.cd()
    tight_hist = ROOT.TH1F(title, title, 1000, 0, 50) 
    for b in range(h_egamma_tight.GetNbinsX()):
        tight_hist.SetBinContent(b+1,h_egamma_tight.GetBinContent(b+1))
    tight_hist.SetName(title)
    tight_hist.Write()
    outfile.Close()

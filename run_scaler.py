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

# command line options
parser = argparse.ArgumentParser(description="")
parser.add_argument("input", metavar="INPUT", help="input directory file with summed_egamma.root")
parser.add_argument("--testBin", default=None, help="specify bin to test")
parser.add_argument("--phislice", default=False, action="store_true", help="")
args = parser.parse_args()

# Constants
directory_name = "scaled_tight_hists"
phislice_directory_name = "scaled_phislice_tight_hists"
bins = [20,40,60,80,100,140,180,220,300,380]
phi_bins = [0, 500, 1000, 1500, 2000]
eta_regions = ["barrel", "endcap"]
regions = ["iso_sym", "iso_asym", "noniso_sym", "noniso_asym"]
SEED = 12445

# init
os.chdir(args.input)
infile1 = ROOT.TFile('summed_egamma.root')
if args.testBin is None: test_regions = ["noniso_sym"] 
elif "noniso_asym" in args.testBin: test_regions = ["noniso_asym"]
elif "noniso_sym" in args.testBin: test_regions = ["noniso_sym"]
elif "iso_asym" in args.testBin: test_regions = ["iso_asym"]
elif "iso_sym" in args.testBin: test_regions = ["iso_sym"]
if args.testBin is not None: test_bin = (args.testBin).split(" ")
if not os.path.exists(directory_name): os.mkdir(directory_name)
os.chdir(directory_name)

# pt-binned only histos
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
    if not region == "iso_sym": util.removeEntries(h_egamma_tight, infile1.Get(egamma_tight_plots.replace(region, "iso_sym")), seed=SEED)

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

    egamma_tight_plots = "plots/twoprong_masspi0_phi" + region + "_" + eta_reg
    if i == len(bins)-1: egamma_tight_plots += "_" + str(bins[i]) + "+"
    else: egamma_tight_plots += "_" + str(bins[i]) + "_" + str(bins[i+1])
    egamma_tight_plots += "_tight"

if not args.phislice: exit()
os.chdir('../')
if not os.path.exists(phislice_directory_name): os.mkdir(phislice_directory_name)
os.chdir(phislice_directory_name)

# pt-and-phi-binned histos
for region, i, eta_reg, k in itertools.product(regions, range(len(bins)), eta_regions, range(len(phi_bins))):
    if args.testBin:
        if region != test_bin[0]: continue
        if eta_reg != test_bin[1]: continue
        if bins[i] != int(test_bin[2]): continue
    histogram_prefix_mass = "twoprong_masspi0_"
    phi_low = phi_bins[k]
    phi_high = phi_bins[k+1] if k != len(phi_bins)-1 else "Inf"
    control_region = region + "_" + eta_reg
    pt_low = bins[i]
    pt_high = bins[i+1] if i != len(bins)-1 else "Inf"
    tight_plot_name = 'plots/{}phi{}-{}_{}_pt{}-{}_tight'.format(histogram_prefix_mass, phi_low, phi_high, control_region, pt_low, pt_high)
    h_egamma_tight = infile1.Get(tight_plot_name)
    h_egamma_signal_region = infile1.Get(tight_plot_name.replace(region, "iso_sym"))

    # Scale
    if not region == "iso_sym": util.removeEntries(h_egamma_tight, h_egamma_signal_region)

    # Save
    save_name = '{}phi{}-{}_{}_pt{}-{}_tight'.format(histogram_prefix_mass, phi_low, phi_high, control_region, pt_low, pt_high)
    outfile = ROOT.TFile(save_name + ".root", "RECREATE")
    outfile.cd()
    save_hist = h_egamma_tight.Clone()
    save_hist.SetName(save_name)
    save_hist.Reset()
    for b in range(h_egamma_tight.GetNbinsX()):
        save_hist.SetBinContent(b+1,h_egamma_tight.GetBinContent(b+1))
    save_hist.Write()
    outfile.Close()

from __future__ import print_function
import math
import ROOT
import sys
import os
import argparse
import array
import fitting_utils as util
import constants as VALS
import scipy
import scipy.stats as stats
import itertools

# command line options
parser = argparse.ArgumentParser(description="")
parser.add_argument("input", metavar="INPUT", help="input directory file with summed_egamma.root")
parser.add_argument("--testBin", default=None, help="specify bin to test")
parser.add_argument("--nophislice", default=False, action="store_true", help="")
parser.add_argument("--scaleTo", "-s", default="same", choices=["same", "overall"], help="")
args = parser.parse_args()

# Constants
SEED = 12445
directory_name = "scaled_tight_hists"
phislice_directory_name = "scaled_phislice_tight_hists"
eta_regions = ["barrel", "endcap"]
regions = ["iso_sym", "iso_asym", "noniso_sym", "noniso_asym"]

# init
os.chdir(args.input)
infile1 = ROOT.TFile('summed_egamma.root')
rand = ROOT.TRandom3()
rand.SetSeed(SEED)
if args.testBin is not None: test_bin = (args.testBin).split(" ")
if not os.path.exists(directory_name): os.mkdir(directory_name)
os.chdir(directory_name)

# pt-binned only histos
for region, i, eta_reg in itertools.product(regions, range(len(VALS.PT_EDGES)), eta_regions):
    if region == "iso_sym": continue
    if args.testBin:
        if region != test_bin[0]: continue
        if eta_reg != test_bin[1]: continue
        if VALS.PT_EDGES[i] != int(test_bin[2]): continue
    egamma_tight_plots = "plots/twoprong_masspi0_" + region + "_" + eta_reg
    if i == len(VALS.PT_EDGES)-1: egamma_tight_plots += "_" + str(VALS.PT_EDGES[i]) + "+"
    else: egamma_tight_plots += "_" + str(VALS.PT_EDGES[i]) + "_" + str(VALS.PT_EDGES[i+1])
    egamma_tight_plots += "_tight"
    h_egamma_tight = infile1.Get(egamma_tight_plots)

    # Scale
    target = infile1.Get(egamma_tight_plots.replace(region, "iso_sym")).Integral()
    util.removeEntries(h_egamma_tight, target, rand=rand)

    # Save
    if i == len(VALS.PT_EDGES)-1: title = region + "_" + eta_reg + "_" + str(VALS.PT_EDGES[i]) + "+_tight"
    else: title = region + "_" + eta_reg + "_" + str(VALS.PT_EDGES[i]) + "_" + str(VALS.PT_EDGES[i+1]) + "_tight" 
    outfile = ROOT.TFile(title + ".root", "RECREATE")
    outfile.cd()
    tight_hist = ROOT.TH1F(title, title, 1000, 0, 50) 
    for b in range(h_egamma_tight.GetNbinsX()):
        tight_hist.SetBinContent(b+1,h_egamma_tight.GetBinContent(b+1))
    tight_hist.SetName(title)
    tight_hist.Write()
    outfile.Close()

    egamma_tight_plots = "plots/twoprong_masspi0_phi" + region + "_" + eta_reg
    if i == len(VALS.PT_EDGES)-1: egamma_tight_plots += "_" + str(VALS.PT_EDGES[i]) + "+"
    else: egamma_tight_plots += "_" + str(VALS.PT_EDGES[i]) + "_" + str(VALS.PT_EDGES[i+1])
    egamma_tight_plots += "_tight"

if args.nophislice: exit()
os.chdir('../')
if not os.path.exists(phislice_directory_name): os.mkdir(phislice_directory_name)
os.chdir(phislice_directory_name)

if args.scaleTo == "overall":
    full_integrals = {}
    for region, eta_reg in itertools.product(regions, eta_regions):
        total_integral = 0
        for i, k in itertools.product(range(len(VALS.PT_EDGES)), range(len(VALS.PHI_EDGES))):
            histogram_prefix_mass = "twoprong_masspi0_"
            phi_low = VALS.PHI_EDGES[k]
            phi_high = VALS.PHI_EDGES[k+1] if k != len(VALS.PHI_EDGES)-1 else "Inf"
            control_region = region + "_" + eta_reg
            pt_low = VALS.PT_EDGES[i]
            pt_high = VALS.PT_EDGES[i+1] if i != len(VALS.PT_EDGES)-1 else "Inf"
            tight_plot_name = 'plots/{}phi{}-{}_{}_pt{}-{}_tight'.format(histogram_prefix_mass, phi_low, phi_high, control_region, pt_low, pt_high)
            h_egamma_tight = infile1.Get(tight_plot_name)
            total_integral += h_egamma_tight.Integral()
        full_integrals[region + "_" + eta_reg] = total_integral
    print(full_integrals)

# pt-and-phi-binned histos
for region, i, eta_reg, k in itertools.product(regions, range(len(VALS.PT_EDGES)), eta_regions, range(len(VALS.PHI_EDGES))):
    if region == "iso_sym": continue
    if args.testBin:
        if region != test_bin[0]: continue
        if eta_reg != test_bin[1]: continue
        if VALS.PT_EDGES[i] != int(test_bin[2]): continue
    histogram_prefix_mass = "twoprong_masspi0_"
    phi_low = VALS.PHI_EDGES[k]
    phi_high = VALS.PHI_EDGES[k+1] if k != len(VALS.PHI_EDGES)-1 else "Inf"
    control_region = region + "_" + eta_reg
    pt_low = VALS.PT_EDGES[i]
    pt_high = VALS.PT_EDGES[i+1] if i != len(VALS.PT_EDGES)-1 else "Inf"
    tight_plot_name = 'plots/{}phi{}-{}_{}_pt{}-{}_tight'.format(histogram_prefix_mass, phi_low, phi_high, control_region, pt_low, pt_high)
    h_egamma_tight = infile1.Get(tight_plot_name)
    h_egamma_signal_region = infile1.Get(tight_plot_name.replace(region, "iso_sym"))

    # Scale
    if args.scaleTo == "overall":
        signal_region = "iso_sym_" + eta_reg
        target = h_egamma_tight.Integral() * (full_integrals[signal_region] / full_integrals[control_region])
    if args.scaleTo == "same": target = h_egamma_signal_region.Integral()
    util.removeEntries(h_egamma_tight, target)

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

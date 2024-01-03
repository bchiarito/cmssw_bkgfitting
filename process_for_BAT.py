import uproot
import numpy as np
import sys
import os
import argparse

parser = argparse.ArgumentParser("")
parser.add_argument("input_dir", help="")
parser.add_argument("--old", action='store_true', help="don't include phi-binned tight data")
args = parser.parse_args()

# Constants
pt_bins = [20,40,60,80,100,140,180,220,300,380]
phi_bins = [0, 500, 1000, 1500, 2000]
control_region = "noniso_asym_barrel"
signal_region = "iso_sym_barrel"
egamma_filename = "summed_egamma.root"
mc_filename = "summed_gjets.root"
tight_data_basename = "tight_data_"
tight_data_old_basename = "tight_data_old_"
loose_fit_basename = "loose_fit_"
bin_centers_filename = "bin_centers.txt"
pt_spectrum_basename = "pt_spectrum_"
histogram_prefix_mass = "twoprong_masspi0_"
histogram_prefix_pt = "twoprong_pt_"
loose_dir = "loose_fit_hists"
scaled_tight_dir = "scaled_tight_hists"

rootdir = os.path.normpath(args.input_dir)
file = uproot.open(rootdir+"/"+egamma_filename)
file_mc = uproot.open(rootdir+"/"+mc_filename)

for i in range(len(pt_bins)):
    try:
        pt_low = pt_bins[i]
        pt_high = pt_bins[i+1] if i != len(pt_bins)-1 else "Inf"

        tight_filename = rootdir+'/'+scaled_tight_dir+"/{}_{}_{}_tight.root".format(control_region, pt_low, pt_high)
        loose_filename = rootdir+'/'+loose_dir+"/{}_{}_{}_loose.root".format(control_region, pt_low, pt_high)
        tight_filename = tight_filename.replace("_Inf", "+")
        loose_filename = loose_filename.replace("_Inf", "+")

        tight_hist_name = "{}_{}_{}_tight".format(control_region, pt_low, pt_high)
        loose_hist_name = "{}_{}_{}_loose".format(control_region, pt_low, pt_high)
        tight_hist_name = tight_hist_name.replace("_Inf", "+")
        loose_hist_name = loose_hist_name.replace("_Inf", "+")
        file1 = uproot.open(tight_filename)
        file2 = uproot.open(loose_filename)

        if args.scaled:
            tight_hist = file1[tight_hist_name]
        else:
            tight_hist_name = "plots/{}".format(histogram_prefix_mass) + tight_hist_name
            tight_hist = file[tight_hist_name]

        loose_fit_hist = file2[loose_hist_name]

        tight_data, edges = tight_hist.to_numpy()
        loose_fit, _ = loose_fit_hist.to_numpy()
        centers = np.array([(edges[j] + edges[j])/2.0 for j in range(len(edges)-1)])

        np.savetxt(rootdir+'/{}{}.txt'.format(tight_data_old_basename, i+1), tight_data)
        np.savetxt(rootdir+'/{}{}.txt'.format(loose_fit_basename, i+1), loose_fit)
        if i==0: np.savetxt(rootdir+'/{}'.format(bin_centers_filename), centers)
    except FileNotFoundError as err:
        print("Failed for pt bin: ", pt_bins[i])
        continue

if args.old:
    mc_fraction_hist_name = "plots/{}".format(histogram_prefix_pt)+signal_region.replace("_barrel", "").replace("_endcap", "")+"_tight"
    mc_fraction_hist = file[mc_fraction_hist_name]
    mc_fraction, _ = mc_fraction_hist.to_numpy()
    np.savetxt(rootdir+'/{}.txt'.format(pt_spectrum_basename), mc_fraction)
    exit()

for i in range(len(pt_bins)):
    pt_low = pt_bins[i]
    pt_high = pt_bins[i+1] if i != len(pt_bins)-1 else "Inf"
    for k in range(len(phi_bins)):
        phi_low = phi_bins[k]
        phi_high = phi_bins[k+1] if k != len(phi_bins)-1 else "Inf"
        try:
            tight_hist = file['plots/{}phi{}-{}_{}_pt{}-{}_tight'.format(histogram_prefix_mass, phi_low, phi_high, control_region, pt_low, pt_high)]
            tight_data, _ = tight_hist.to_numpy()
            np.savetxt(rootdir+'/{}phi{}_pt{}.txt'.format(tight_data_basename, k+1, i+1), tight_data)
        except uproot.KeyInFileError:
            print("Could not find", 'plots/{}_phi{}-{}_{}_pt{}-{}_tight'.format(histogram_prefix_mass, phi_low, phi_high, control_region, pt_low, pt_high))

    tight_hist = file['plots/{}phi{}_{}_pt{}-{}_tight'.format(histogram_prefix_mass, 'All', control_region, pt_low, pt_high)]
    tight_data, _ = tight_hist.to_numpy()
    np.savetxt(rootdir+'/{}phi{}_pt{}.txt'.format(tight_data_basename, 'All', i+1), tight_data)

for k in range(len(phi_bins)):
    phi_low = phi_bins[k]
    phi_high = phi_bins[k+1] if k != len(phi_bins)-1 else "Inf"
    mc_fraction_hist_name = "plots/twoprong_pt_phi{}-{}_".format(phi_low, phi_high)+signal_region+"_ptX-X_tight"
    mc_fraction_hist = file_mc[mc_fraction_hist_name]
    mc_fraction, _ = mc_fraction_hist.to_numpy()
    np.savetxt(rootdir+'/{}{}.txt'.format(pt_spectrum_basename, k+1), mc_fraction)

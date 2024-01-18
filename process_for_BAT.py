import uproot
import numpy as np
import sys
import os
import argparse
import subprocess
import itertools

parser = argparse.ArgumentParser("")
parser.add_argument("input_dir", help="")
parser.add_argument("--clean", action="store_true", default=False, help="")
args = parser.parse_args()

# Constants
pt_bins = [20,40,60,80,100,140,180,220,300,380]
phi_bins = [0, 500, 1000, 1500, 2000]
regions = ["iso_sym", "iso_asym", "noniso_sym", "noniso_asym"]
eta_regions = ["barrel", "endcap"]
signal_regions = ["iso_sym_barrel", "iso_sym_endcap"]
DEBUG = False

# Names
egamma_filename = "summed_egamma.root"
mc_filename = "summed_gjets.root"
bat_dir_base = "text_format_"
loose_fit_basename = "loose_fit_"
bin_centers_filename = "bin_centers.txt"
pt_spectrum_basename = "pt_spectrum_"
histogram_prefix_mass = "twoprong_masspi0_"
histogram_prefix_pt = "twoprong_pt_"
loose_dir = "loose_fit_hists"
scaled_tight_dir = "scaled_tight_hists"
scaled_phislice_tight_dir = "scaled_phislice_tight_hists"
tight_data_scaled_basename = "tight_data_"
tight_data_unscaled_basename = "tight_data_unscaled_"
tight_data_old_scaled_basename = "tight_data_old_"
tight_data_old_unscaled_basename = "tight_data_unscaled_old_"

# Init
rootdir = os.path.normpath(args.input_dir)
os.chdir(rootdir)
file = uproot.open(egamma_filename)
file_mc = uproot.open(mc_filename)

control_regions = [el[0]+'_'+el[1] for el in itertools.product(regions, eta_regions)]
for control_region in control_regions:
    if DEBUG: print(control_region)
    bat_dir = bat_dir_base + control_region

    # Don't run
    if args.clean:
        subprocess.run("rm -r {}".format(bat_dir), shell=True)
        continue

    if control_region in signal_regions: continue
    if not os.path.exists(bat_dir): os.mkdir(bat_dir)
    os.chdir(bat_dir)

    # pt spectrum
    for k in range(len(phi_bins)):
        phi_low = phi_bins[k]
        phi_high = phi_bins[k+1] if k != len(phi_bins)-1 else "Inf"
        if 'barrel' in control_region: region_for_pt_spectrum = "iso_sym_barrel"
        if 'endcap' in control_region: region_for_pt_spectrum = "iso_sym_endcap"
        mc_fraction_hist_name = "plots/twoprong_pt_phi{}-{}_".format(phi_low, phi_high)+region_for_pt_spectrum+"_ptX-X_tight"
        mc_fraction_hist = file_mc[mc_fraction_hist_name]
        mc_fraction, _ = mc_fraction_hist.to_numpy()
        np.savetxt('{}{}.txt'.format(pt_spectrum_basename, k+1), mc_fraction)

    # loose fits and pt-binned tight data
    for i in range(len(pt_bins)):
        pt_low = pt_bins[i]
        pt_high = pt_bins[i+1] if i != len(pt_bins)-1 else "Inf"

        tight_filename = '../'+scaled_tight_dir+"/{}_{}_{}_tight.root".format(control_region, pt_low, pt_high)
        loose_filename = '../'+loose_dir+"/{}_{}_{}_loose.root".format(control_region, pt_low, pt_high)
        tight_filename = tight_filename.replace("_Inf", "+")
        loose_filename = loose_filename.replace("_Inf", "+")

        tight_hist_name = "{}_{}_{}_tight".format(control_region, pt_low, pt_high)
        loose_hist_name = "{}_{}_{}_loose".format(control_region, pt_low, pt_high)
        tight_hist_name = tight_hist_name.replace("_Inf", "+")
        loose_hist_name = loose_hist_name.replace("_Inf", "+")

        scaled_tight_file = uproot.open(tight_filename)
        loose_file = uproot.open(loose_filename)

        loose_fit_hist = loose_file[loose_hist_name]
        tight_hist_scaled = scaled_tight_file[tight_hist_name]

        tight_hist_name = "plots/{}".format(histogram_prefix_mass) + tight_hist_name
        tight_hist_unscaled = file[tight_hist_name]

        loose_fit, edges = loose_fit_hist.to_numpy()

        tight_data_scaled, _ = tight_hist_scaled.to_numpy()
        tight_data_unscaled, _ = tight_hist_unscaled.to_numpy()

        if i==0:
            centers = np.array([(edges[j] + edges[j])/2.0 for j in range(len(edges)-1)])
            np.savetxt('{}'.format(bin_centers_filename), centers)
        np.savetxt('{}{}.txt'.format(loose_fit_basename, i+1), loose_fit)
        np.savetxt('{}{}.txt'.format(tight_data_old_scaled_basename, i+1), tight_data_scaled)
        np.savetxt('{}{}.txt'.format(tight_data_old_unscaled_basename, i+1), tight_data_unscaled)

    # phi+pt-binned tight data
    for i in range(len(pt_bins)):
        pt_low = pt_bins[i]
        pt_high = pt_bins[i+1] if i != len(pt_bins)-1 else "Inf"
        for k in range(len(phi_bins)):
            phi_low = phi_bins[k]
            phi_high = phi_bins[k+1] if k != len(phi_bins)-1 else "Inf"
            tight_hist_unscaled = file['plots/{}phi{}-{}_{}_pt{}-{}_tight'.format(histogram_prefix_mass, phi_low, phi_high, control_region, pt_low, pt_high)]
            tight_data_unscaled, _ = tight_hist_unscaled.to_numpy()
            np.savetxt('{}phi{}_pt{}.txt'.format(tight_data_unscaled_basename, k+1, i+1), tight_data_unscaled)
        
            scaled_name = '{}phi{}-{}_{}_pt{}-{}_tight'.format(histogram_prefix_mass, phi_low, phi_high, control_region, pt_low, pt_high)
            scaled_tight_file = uproot.open('../'+scaled_phislice_tight_dir+'/'+scaled_name+'.root')
            tight_hist_scaled = scaled_tight_file[scaled_name]
            tight_data_scaled, _ = tight_hist_scaled.to_numpy()
            np.savetxt('{}phi{}_pt{}.txt'.format(tight_data_scaled_basename, k+1, i+1), tight_data_scaled)

            tight_data_unscaled, _ = tight_hist_unscaled.to_numpy()
            np.savetxt('{}phi{}_pt{}.txt'.format(tight_data_unscaled_basename, k+1, i+1), tight_data_unscaled)

        tight_hist_unscaled = file['plots/{}phi{}_{}_pt{}-{}_tight'.format(histogram_prefix_mass, 'All', control_region, pt_low, pt_high)]
        tight_data_unscaled, _ = tight_hist_unscaled.to_numpy()
        np.savetxt('{}phi{}_pt{}.txt'.format(tight_data_unscaled_basename, 'All', i+1), tight_data_unscaled)
    os.chdir("../")

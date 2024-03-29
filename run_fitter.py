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
#import scipy.stats as stats

# command line options
parser = argparse.ArgumentParser(description="")
parser.add_argument("input", metavar="INPUT", help="input root file")
parser.add_argument("--testRegion", default=None, choices = ["iso_sym", "iso_asym", "noniso_sym", "noniso_asym"], help="specify region to test")
parser.add_argument("--testBin", default=None, help="specify bin to test")
parser.add_argument("--name", default="plots", help="create name for plots pdf")
parser.add_argument("--show", default=False, action="store_true", help="don't close after finished running (to view canvas)")
run_args = parser.add_argument_group("run options")
run_args.add_argument("--createLooseFits", default=False, action="store_true", help="will recreate loose fits even if present")
run_args.add_argument("--useUnscaledTight", default=False, action="store_true", help="use unscaled tight data in preference to scaled")
run_args.add_argument("--ftest", default="3 4", help="change ftest, format: '<CHEB_TYPE> <MAXDEGREE>', default is cheby degree 4")
run_args.add_argument("--integral", default=False, action="store_true", help="add I to tight fit")
plot_args = parser.add_argument_group("plotting options")
plot_args.add_argument("--checkPull", default=False, action="store_true", help="print on legend if there are four consecutive pull bins greater than 1.5 sigma")
plot_args.add_argument("--specifyFtestDegree", "--fdeg", default=None, help="specify which degree ftest should pick for visualization purposes")
plot_args.add_argument("--printFtest", "--printftest", default=False, action="store_true", help="for a fixed test bin, create a pdf of all possible ftest fits")
plot_args.add_argument("--onlyLoose", default=False, action="store_true", help="create loose fit plots only")
plot_args.add_argument("--sanity", "-s", default=False, action="store_true", help="create sanity plots")

# parse args
args = parser.parse_args()

# constants
egamma_rootfile = 'summed_egamma.root'

# plot style
leg_x1, leg_x2, leg_y1, leg_y2 = 0.7, 0.60, 0.89, 0.89
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)
ROOT.gStyle.SetLegendFillColor(ROOT.TColor.GetColorTransparent(ROOT.kWhite, 0.01));
ROOT.gStyle.SetLegendBorderSize(0)

# select regions
# pi0: masspi0 plots for all eta regions, barrel, and endcap
# pi0_VALS.PT_EDGES: pt-binned masspi0 plots in barrel and endcap; 
# overlay; pt-binned plots with overlayed ratios for each twoprong region
if args.sanity: plots = ["sieie", "pfRelIso03_chg", "hadTow"]  
else: plots = ["pi0_VALS.PT_EDGES"]
regions = ["iso_sym", "iso_asym", "noniso_sym", "noniso_asym"]
if args.testRegion: regions = [args.testRegion]
eta_regions = ["barrel", "endcap"]
photon_regions = ["tight", "loose"]

# init
os.chdir(args.input)
infile1 = ROOT.TFile(egamma_rootfile)
c1 = ROOT.TCanvas("c1", "c1", 800, 600)
c1.Print(args.name + ".pdf[")

# run
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
    elif item == "pi0":  # un-binned plots (needs to be updated)
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

    elif item == "pi0_VALS.PT_EDGES":  # binned plots (PRIMARILY USED)
        if args.testBin is not None: test_bin = (args.testBin).split(" ")
        chi2_pvalues = []
        if not os.path.exists('loose_fit_hists'): os.mkdir("loose_fit_hists")
        if not os.path.exists('tight_templates'): os.mkdir("tight_templates")
        os.chdir('tight_templates')
        if args.useUnscaledTight:
            if not os.path.exists('templates_noscaling'): os.mkdir("templates_noscaling")
        else:
            if not os.path.exists('templates'): os.mkdir("templates")
        if not os.path.exists('degrees'): os.mkdir("degrees")
        if not os.path.exists('chi2s'): os.mkdir("chi2s")
        os.chdir('../')
        if not os.path.exists('tight_templates'): os.mkdir("tight_templates")
        os.chdir('tight_templates')
        if not os.path.exists('polys'): os.mkdir('polys')
        os.chdir('../')
            
        loose_fit_files = []
        scaled_tight_files = []
        file_counter_loose = 0
        file_counter_tight = 0
        for region in regions:  # loop through twoprong sideband regions
            if args.printFtest and args.testBin is None:
                print("EMPTY PDF: Must have --testBin option when using --printFtest")
                break
            if args.testBin is not None: 
                if not region == test_bin[0]: continue
            for i in range(len(VALS.PT_EDGES)):  # loop through pt VALS.PT_EDGES for a fixed twoprong sideband
                if args.testBin is not None: 
                    if not VALS.PT_EDGES[i] == int(test_bin[2]): continue
                for eta_reg in eta_regions:  # loop through eta regions for fixed pt-bin and fixed twoprong sideband
                    if args.testBin is not None: 
                        if not eta_reg == test_bin[1]: continue
                    if not eta_reg == "barrel" and not eta_reg == "endcap": continue  # no pt-bin plots for barrel and endcap combined, so skip this case
                    
                    # Generate correct plots names to access from summed histogram files
                    egamma_tight_plots = "plots/twoprong_masspi0_" + region + "_" + eta_reg
                    egamma_loose_plots = "plots/twoprong_masspi0_" + region + "_" + eta_reg
                    
                    # this must follow the naming convention of the input histograms
                    if i == len(VALS.PT_EDGES) - 1:
                        egamma_tight_plots += "_" + str(VALS.PT_EDGES[i]) + "+"
                        egamma_loose_plots += "_" + str(VALS.PT_EDGES[i]) + "+"
                    else:
                        egamma_tight_plots += "_" + str(VALS.PT_EDGES[i]) + "_" + str(VALS.PT_EDGES[i+1])
                        egamma_loose_plots += "_" + str(VALS.PT_EDGES[i]) + "_" + str(VALS.PT_EDGES[i+1]) 
                    
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
                    
                    if i == len(VALS.PT_EDGES) - 1: hist_name = region + "_" + eta_reg + "_" + str(VALS.PT_EDGES[i]) + "+"
                    else: hist_name = region + "_" + eta_reg + "_" + str(VALS.PT_EDGES[i]) + "_" + str(VALS.PT_EDGES[i+1]) 
                    
                    # scaled tight data
                    if not args.useUnscaledTight and region != "iso_sym":
                        scaled_tight_files.append(ROOT.TFile("scaled_tight_hists/" + hist_name + "_tight.root"))
                        h_scaled_tight = scaled_tight_files[file_counter_tight].Get(hist_name+"_tight")
                        h_egamma_tight.Reset()
                        for b in range(h_scaled_tight.GetNbinsX()): 
                            h_egamma_tight.SetBinContent(b+1,h_scaled_tight.GetBinContent(b+1))
                        file_counter_tight += 1

                    # FITTING
                    if i == 0 and eta_reg == "barrel": print("====================== " + region.upper() + " =====================")
                    if i == len(VALS.PT_EDGES) - 1: print("############### " + region.upper() + " " + eta_reg.upper() + " " + str(VALS.PT_EDGES[i]) + "+ ###############")
                    else: print("############### " + region.upper() + " " + eta_reg.upper() + " " + str(VALS.PT_EDGES[i]) + "-" + str(VALS.PT_EDGES[i+1]) + " ###############")
                    
                    ### new idea, fit only rising and falling ###
                    # determine left and right bounds
                    # fit landaus from first_bin to left_bin
                    # fit exps from right_bin to last_bin
                    ENTRIES_CUTOFF = 1000  # adjust based on bulk region
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
                    first = h_egamma_loose.GetBinLowEdge(first_bin)  # first bin with data
                    left = h_egamma_loose.GetBinLowEdge(left_bin+1)  # first bin exceeding ENTRIES_CUTOFF
                    right = h_egamma_loose.GetBinLowEdge(right_bin)  # last bin exceeding ENTRIES_CUTOFF
                    last = h_egamma_loose.GetBinLowEdge(last_bin+1)  # last bin with data
                    
                    fit_init = util.lookup_fit_guesses(region, eta_reg, VALS.PT_EDGES[i])
                    old_method = fit_init['old_method']
                    guesses = fit_init['guesses']
                    nLandau = fit_init['nLandau']
                    landau_guess = fit_init['landau_guess']
                    nExp = fit_init['nExp']
                    exp_guess = fit_init['exp_guess']

                    # Loose spectrum fitting
                    os.chdir("loose_fit_hists")
                    if i == len(VALS.PT_EDGES) - 1: title = region + "_" + eta_reg + "_" + str(VALS.PT_EDGES[i]) + "+_loose"
                    else: title = region + "_" + eta_reg + "_" + str(VALS.PT_EDGES[i]) + "_" + str(VALS.PT_EDGES[i+1]) + "_loose" 
                    if not os.path.exists(title + ".root") or args.createLooseFits: # create new loose fits
                        if old_method:
                            N = str(nLandau) + str(nExp)
                            func_full, fitresult_full = util.fit_hist(h_egamma_loose, 'full', 0, 50, int(N), initial_guesses=guesses)
                            loose_fit_as_hist = util.TemplateToHistogram(func_full, 1000, 0, 50)  # the histogram bin definition must align with the input loose and tight histograms
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
                                if b < left_bin:
                                    loose_fit_as_hist.SetBinContent(b+1, rising_fit_as_hist.GetBinContent(b+1))
                                elif b <= right_bin:
                                    loose_fit_as_hist.SetBinContent(b+1, h_egamma_loose.GetBinContent(b+1)) 
                                else:
                                    loose_fit_as_hist.SetBinContent(b+1, falling_fit_as_hist.GetBinContent(b+1))
                        
                        # Save the loose fits in a separate file
                        #if i == len(VALS.PT_EDGES) - 1: title = region + "_" + eta_reg + "_" + str(VALS.PT_EDGES[i]) + "+_loose"
                        #else: title = region + "_" + eta_reg + "_" + str(VALS.PT_EDGES[i]) + "_" + str(VALS.PT_EDGES[i+1]) + "_loose" 
                        outfile = ROOT.TFile(title + ".root", "RECREATE")
                        outfile.cd()
                        loose_hist = ROOT.TH1F(title, title, 1000, 0, 50) 
                        for b in range(loose_fit_as_hist.GetNbinsX()):
                            loose_hist.SetBinContent(b+1,loose_fit_as_hist.GetBinContent(b+1))
                        loose_hist.SetName(title)
                        loose_hist.Write()
                        outfile.Close()
                        os.chdir("../")
                    else: # read from already present fits
                        os.chdir("../")
                        
                    # Input loose fit templates saved previously for fitting to the tight data
                    loose_fit_files.append(ROOT.TFile("loose_fit_hists/" + hist_name + "_loose.root"))
                    loose_fit_as_hist = loose_fit_files[file_counter_loose].Get(hist_name+"_loose")
                    file_counter_loose += 1

                    fitted_func = util.HistogramToFunction(loose_fit_as_hist)
                    fitted_func_times_constant, _, _ = util.MultiplyWithPolyToTF1(fitted_func, 0, poly=0)  # this part can appear a bit buggy in the plots for some reason
                    fit_result = h_egamma_tight.Fit(fitted_func_times_constant, '0L' if not args.integral else '0LI')
                    tight_fit_w_constant = util.TemplateToHistogram(fitted_func_times_constant, 1000, 0, 50)  

                    # Decide whether an F-test should be used to pick the best polynomial degree
                    FTEST = True
                    NUM_PLOTS = 1  # this variable is used so that the --printFtest option works
                    if not FTEST:
                        POLY_TYPE = 3
                        DEGREE = 0 
                        func_with_poly, func_with_ploy_py, _ = util.MultiplyWithPolyToTF1(fitted_func, DEGREE, poly=POLY_TYPE)
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
                            func_with_poly, func_with_poly_py, _ = util.MultiplyWithPolyToTF1(fitted_func, degree, poly=POLY_TYPE)
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
                        
                        if args.specifyFtestDegree is not None: best_d = int(args.specifyFtestDegree)
                        print('Best: ', best_d)
                        func_with_poly = fitfuncs[best_d]
                        tight_fit_as_hist = util.TemplateToHistogram(func_with_poly, 1000, 0, 50)
                        tight_stat = statboxes[best_d]
                
                    # Save the tight templates to be plotted separately
                    if not args.useUnscaledTight: os.chdir("tight_templates/templates/")
                    else: os.chdir("tight_templates/templates_noscaling/")

                    # Save the loose fits in a separate file
                    if i == len(VALS.PT_EDGES) - 1: title = region + "_" + eta_reg + "_" + str(VALS.PT_EDGES[i]) + "+_tight_temp"
                    else: title = region + "_" + eta_reg + "_" + str(VALS.PT_EDGES[i]) + "_" + str(VALS.PT_EDGES[i+1]) + "_tight_temp" 
                    outfile1 = ROOT.TFile(title + ".root", "RECREATE")
                    outfile1.cd()
                    tight_fit_hist = ROOT.TH1F(title, title, 1000, 0, 50) 
                    for b in range(tight_fit_as_hist.GetNbinsX()): tight_fit_hist.SetBinContent(b+1,tight_fit_as_hist.GetBinContent(b+1))
                    tight_fit_hist.SetName(title)
                    tight_fit_hist.Write()
                    outfile1.Close()
                    os.chdir("../degrees/")
                    outfile2 = ROOT.TFile(title+"_deg.root", "RECREATE")
                    outfile2.cd()
                    tight_fit_deg = ROOT.TH1F(title+"_deg", title+"_deg", 10, 0, 10)
                    for b in range(tight_fit_deg.GetNbinsX()): tight_fit_deg.SetBinContent(b+1,best_d)
                    tight_fit_deg.SetName(title+"_deg")
                    tight_fit_deg.Write()
                    outfile2.Close()
                    os.chdir("../../")

                    if i == len(VALS.PT_EDGES) - 1: hist_name = region + "_" + eta_reg + "_" + str(VALS.PT_EDGES[i]) + "+"
                    else: hist_name = region + "_" + eta_reg + "_" + str(VALS.PT_EDGES[i]) + "_" + str(VALS.PT_EDGES[i+1]) 
                    
                    # The plotting done here are for rough visualization purposes
                    # For finalized plots, a separate plotting script should be used
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

                        os.chdir("tight_templates/polys/")
                        # Save the loose fits in a separate file
                        if i == len(VALS.PT_EDGES) - 1: title = region + "_" + eta_reg + "_" + str(VALS.PT_EDGES[i]) + "+_tight_poly"
                        else: title = region + "_" + eta_reg + "_" + str(VALS.PT_EDGES[i]) + "_" + str(VALS.PT_EDGES[i+1]) + "_tight_poly" 
                        outfile1 = ROOT.TFile(title + ".root", "RECREATE")
                        outfile1.cd()
                        bern_poly = just_poly.Clone()
                        bern_poly.SetName(title)
                        bern_poly.Write()
                        outfile1.Close()
                        os.chdir("../../")
                    
                        # determine bin-by-bin error
                        STEP_SIZE = 0.005
                        integral = False
                        hist = h_egamma_tight
                        fit = func_with_poly
                        #ndof = util.count_nonzero_VALS.PT_EDGES(hist) - fit.GetNpar()
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

                        chi2_mod, mod_bins = util.RSS(fit, hist, error=0, integral=integral, chi2=True, cutoff=5)
                        num_bins = len(mod_bins)
                        if not num_bins == 0:
                            chi2_mod_ndof = chi2_mod / (num_bins-fit.GetNpar())
                            chi2_pvalues.append(scipy.stats.chi2.sf(chi2_mod, num_bins-fit.GetNpar()))
                        else:
                            chi2_mod_ndof = chi2_mod
                            chi2_pvalues.append(-1)
                        #print(chi2_mod, len(num_bins), chi2_mod_ndof)
                        
                        os.chdir("tight_templates/chi2s/")
                        outfile3 = ROOT.TFile(title+"_chi2.root", "RECREATE")
                        outfile3.cd()
                        tight_fit_chi2 = ROOT.TH1D(title+"_chi2", title+"_chi2", 1000, 0, 100)
                        for b in range(tight_fit_chi2.GetNbinsX()): tight_fit_chi2.SetBinContent(b+1, chi2_mod_ndof)
                        tight_fit_chi2.SetName(title+"_chi2")
                        tight_fit_chi2.Write()
                        outfile2.Close()
                        os.chdir("../../")

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
                        
                        # Create title for plot 
                        title = region + " Twoprong"
                        if eta_reg == "barrel": title += ", Barrel"
                        elif eta_reg == "endcap": title += ", Endcap"
                        if i == len(VALS.PT_EDGES) - 1: title += ", pt > " + str(VALS.PT_EDGES[i])
                        else: title += ", " + str(VALS.PT_EDGES[i]) + " < pt < " + str(VALS.PT_EDGES[i+1]) 
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
                        legend1 = ROOT.TLegend(0.35, 0.78, 0.65, 0.9)
                        legend1.AddEntry(h_egamma_loose, "Loose Photon, " + str(round(h_egamma_loose.GetEntries())), "l")
                        if region == "iso_sym": legend1.AddEntry(h_egamma_tight, "Tight Photon, " + str(round(h_egamma_tight.GetEntries())), "l")
                        legend2 = ROOT.TLegend(0.29, 0.70, 0.62, 0.89)
                        legend2.AddEntry(h_egamma_tight, "Tight Photon, " + str(h_egamma_tight.GetEntries()), "l")
                        if POLY_TYPE == 0: legend2.AddEntry('', 'Polynomial', '')
                        if POLY_TYPE == 1: legend2.AddEntry('', 'Chebyshev 1st kind', '')
                        if POLY_TYPE == 2: legend2.AddEntry('', 'Chebyshev 2nd kind', '')
                        if POLY_TYPE == 3: legend2.AddEntry('', 'Bernstein', '')
                        func_n_par = func_with_poly.GetNpar()
                        if POLY_TYPE == 3 and func_n_par == 1: fit_degree = 0
                        elif POLY_TYPE == 3 and func_n_par > 1: fit_degree = func_n_par - 1
                        elif POLY_TYPE < 3: fit_degree = func_n_par - 1
                        if FTEST: legend2.AddEntry(tight_fit_as_hist, "Fit w f-test (Degree "+str(fit_degree)+")", "l")
                        else: legend2.AddEntry(tight_fit_as_hist, "Fit (Degree "+str(fit_degree)+")", "l")
                        legend2.AddEntry(tight_fit_w_constant, "Constant fit p0 = {:.4}".format(fitted_func_times_constant.GetParameter(0)), "l")
                        legend2.AddEntry('', 'Chi2/Ndof: {:.3f}'.format(chi2_ndof), '')
                        legend2.AddEntry('', 'Chi2_mod/Ndof: {:.3f}'.format(chi2_mod_ndof), '')
                        legend2.AddEntry('', 'Bin Error: {:.1%}'.format(bin_bin_error), '')
                        if args.checkPull:
                            legend2.AddEntry('', '3 Consecutive Bins > 2 sigma: ' + str(utilcheckPullHist(h_tight_pull, 3, 2)), '')
                            legend2.AddEntry('', '4 Consecutive Bins > 1.5 sigma: ' + str(util.checkPullHist(h_tight_pull, 4, 1.5)), '')
                            legend2.AddEntry('', '5 Consecutive Bins > 1 sigma: ' + str(util.checkPullHist(h_tight_pull, 5, 1)), '')

                        # Draw plots
                        if args.onlyLoose: c1.cd(1)
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
                        if VALS.PT_EDGES[i] < 60: h_egamma_loose.GetXaxis().SetRangeUser(0, 5)
                        elif VALS.PT_EDGES[i] < 120: h_egamma_loose.GetXaxis().SetRangeUser(0, 10)
                        elif VALS.PT_EDGES[i] < 200: h_egamma_loose.GetXaxis().SetRangeUser(0, 15)
                        elif VALS.PT_EDGES[i] < 380: h_egamma_loose.GetXaxis().SetRangeUser(0, 20)
                        else: h_egamma_loose.GetXaxis().SetRangeUser(0, 26)
                        legend1.Draw("same")
                        ROOT.gPad.Update()

                        if not args.onlyLoose:
                            if not region == "iso_sym":
                                c1.cd()
                                pad2 = ROOT.TPad('pad2', 'pad2', 0.5, 0.3, 1, 1)
                                pad2.Draw()
                                pad2.cd()
                                ROOT.gPad.SetLogy()
                                if VALS.PT_EDGES[i] < 60: h_egamma_tight.GetXaxis().SetRangeUser(0, 5)
                                elif VALS.PT_EDGES[i] < 120: h_egamma_tight.GetXaxis().SetRangeUser(0, 10)
                                elif VALS.PT_EDGES[i] < 200: h_egamma_tight.GetXaxis().SetRangeUser(0, 15)
                                elif VALS.PT_EDGES[i] < 380: h_egamma_tight.GetXaxis().SetRangeUser(0, 20)
                                else: h_egamma_tight.GetXaxis().SetRangeUser(0, 26)
                                h_egamma_tight.SetMinimum(0.1)
                                h_egamma_tight.Draw("e")
                                if FTEST: tight_stat.Draw()
                                tight_fit_w_constant.SetLineColor(ROOT.kBlue)
                                tight_fit_as_hist.SetLineColor(ROOT.kRed)
                                tight_fit_as_hist.SetLineWidth(1)
                                #tight_fit_as_hist_errorbars = tight_fit_as_hist.Clone()
                                #tight_fit_as_hist_errorbars.SetFillColor(ROOT.kRed+2)
                                #tight_fit_as_hist_errorbars.Draw("same e2")
                                tight_fit_as_hist.Draw("same hist")
                                tight_fit_w_constant.Draw('same')
                                h_egamma_tight.Draw("e same")
                                h_egamma_tight.GetYaxis().SetRangeUser(0.1, h_egamma_tight.GetMaximum()+50)
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
                                if VALS.PT_EDGES[i] < 60: empty.GetXaxis().SetRangeUser(0, 5)
                                elif VALS.PT_EDGES[i] < 120: empty.GetXaxis().SetRangeUser(0, 10)
                                elif VALS.PT_EDGES[i] < 200: empty.GetXaxis().SetRangeUser(0, 15)
                                elif VALS.PT_EDGES[i] < 380: empty.GetXaxis().SetRangeUser(0, 20)
                                else: empty.GetXaxis().SetRangeUser(0, 26)
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
                                #topaxis = ROOT.TGaxis(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUxmin(), ROOT.gPad.GetUxmax(), 510, "+L")
                                #topaxis.SetLineColor(ROOT.kRed);
                                #topaxis.SetLabelColor(ROOT.kRed);
                                #topaxis.Draw()

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
                            if VALS.PT_EDGES[i] < 60: h_loose_pull.GetXaxis().SetRangeUser(0, 5)
                            elif VALS.PT_EDGES[i] < 120: h_loose_pull.GetXaxis().SetRangeUser(0, 10)
                            elif VALS.PT_EDGES[i] < 200: h_loose_pull.GetXaxis().SetRangeUser(0, 15)
                            elif VALS.PT_EDGES[i] < 380: h_loose_pull.GetXaxis().SetRangeUser(0, 20)
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
                                if VALS.PT_EDGES[i] < 60: h_tight_pull.GetXaxis().SetRangeUser(0, 5)
                                elif VALS.PT_EDGES[i] < 120: h_tight_pull.GetXaxis().SetRangeUser(0, 10)
                                elif VALS.PT_EDGES[i] < 200: h_tight_pull.GetXaxis().SetRangeUser(0, 15)
                                elif VALS.PT_EDGES[i] < 380: h_tight_pull.GetXaxis().SetRangeUser(0, 20)
                                else: h_tight_pull.GetXaxis().SetRangeUser(0, 26)
                                h_tight_pullc.SetMarkerColor(ROOT.kBlue)
                                h_tight_pullc.SetMarkerStyle(8)
                                h_tight_pullc.SetMarkerSize(0.25)
                                h_tight_pullc.Draw('pe same')
                                ROOT.gPad.Update()

                        ROOT.gPad.Update()
                        c1.Print(args.name + ".pdf")

    if args.show: input("Finished. Press Enter.")
    c1.Print(args.name + ".pdf]")
    infile1.Close()

#!/bin/env python

import sys
sys.path.append("/home/users/phchang/public_html/phys/rooutil")
import plotmaker
import ROOT as r
from rooutil import plottery_wrapper as p
from plottery import plottery as ply
import math

r.gStyle.SetPalette( r.kInvertedDarkBodyRadiator )

###################################################################################################
# List of root files with histograms after running ScanChain.C
###################################################################################################
tfile_data   = r.TFile("outputs/fakerate_data.root")
tfile_wj     = r.TFile("outputs/fakerate_wj.root")
tfile_dy     = r.TFile("outputs/fakerate_dy.root")
tfile_ttbar  = r.TFile("outputs/fakerate_ttbar.root")
tfile_vv     = r.TFile("outputs/fakerate_vv.root")
tfile_qcd_mu = r.TFile("outputs/fakerate_qcd_mu.root")
tfile_qcd_el = r.TFile("outputs/fakerate_qcd_el.root")

# ROOT is a bit weird and need to keep the TH1's in a global variable to make it persist.
# Otherwise, the histograms get destroyed upon exiting the function.
# Could use SetDirectory(0) like in C++, but PyROOT is a bit finicky with this function.
# It's easier to deal with this with a global list instance.
hists = []

###################################################################################################
#
#
# Functions to access the histograms in the root file
#
#
###################################################################################################
def tightmu_name(varname, region, prefix=""): return prefix+"histo_{}_{}_mu".format(varname, region)
def loosemu_name(varname, region, prefix=""): return prefix+"histo_{}_{}_loose_mu".format(varname, region)
def tightel_name(varname, region, prefix=""): return prefix+"histo_{}_{}_el".format(varname, region)
def looseel_name(varname, region, prefix=""): return prefix+"histo_{}_{}_loose_el".format(varname, region)

def tightmu_file(varname, region, f, n, prefix=""): hists.append(p.move_overflow(f.Get(tightmu_name(varname, region, prefix)).Clone(n))); return hists[-1]
def loosemu_file(varname, region, f, n, prefix=""): hists.append(p.move_overflow(f.Get(loosemu_name(varname, region, prefix)).Clone(n))); return hists[-1]
def tightel_file(varname, region, f, n, prefix=""): hists.append(p.move_overflow(f.Get(tightel_name(varname, region, prefix)).Clone(n))); return hists[-1]
def looseel_file(varname, region, f, n, prefix=""): hists.append(p.move_overflow(f.Get(looseel_name(varname, region, prefix)).Clone(n))); return hists[-1]

def tightmu_data(varname, region, SYSTVAR=""): return tightmu_file(varname, region, tfile_data, "Data", prefix=SYSTVAR)
def loosemu_data(varname, region, SYSTVAR=""): return loosemu_file(varname, region, tfile_data, "Data", prefix=SYSTVAR)
def tightel_data(varname, region, SYSTVAR=""): return tightel_file(varname, region, tfile_data, "Data", prefix=SYSTVAR)
def looseel_data(varname, region, SYSTVAR=""): return looseel_file(varname, region, tfile_data, "Data", prefix=SYSTVAR)

def tightmu_wj(varname, region, SYSTVAR=""): return tightmu_file(varname, region, tfile_wj, "W", prefix=SYSTVAR)
def loosemu_wj(varname, region, SYSTVAR=""): return loosemu_file(varname, region, tfile_wj, "W", prefix=SYSTVAR)
def tightel_wj(varname, region, SYSTVAR=""): return tightel_file(varname, region, tfile_wj, "W", prefix=SYSTVAR)
def looseel_wj(varname, region, SYSTVAR=""): return looseel_file(varname, region, tfile_wj, "W", prefix=SYSTVAR)

def tightmu_dy(varname, region, SYSTVAR=""): return tightmu_file(varname, region, tfile_dy, "Z", prefix=SYSTVAR)
def loosemu_dy(varname, region, SYSTVAR=""): return loosemu_file(varname, region, tfile_dy, "Z", prefix=SYSTVAR)
def tightel_dy(varname, region, SYSTVAR=""): return tightel_file(varname, region, tfile_dy, "Z", prefix=SYSTVAR)
def looseel_dy(varname, region, SYSTVAR=""): return looseel_file(varname, region, tfile_dy, "Z", prefix=SYSTVAR)

def tightmu_ttbar(varname, region, SYSTVAR=""): return tightmu_file(varname, region, tfile_ttbar, "Top", prefix=SYSTVAR)
def loosemu_ttbar(varname, region, SYSTVAR=""): return loosemu_file(varname, region, tfile_ttbar, "Top", prefix=SYSTVAR)
def tightel_ttbar(varname, region, SYSTVAR=""): return tightel_file(varname, region, tfile_ttbar, "Top", prefix=SYSTVAR)
def looseel_ttbar(varname, region, SYSTVAR=""): return looseel_file(varname, region, tfile_ttbar, "Top", prefix=SYSTVAR)

def tightmu_vv(varname, region, SYSTVAR=""): return tightmu_file(varname, region, tfile_vv, "VV", prefix=SYSTVAR)
def loosemu_vv(varname, region, SYSTVAR=""): return loosemu_file(varname, region, tfile_vv, "VV", prefix=SYSTVAR)
def tightel_vv(varname, region, SYSTVAR=""): return tightel_file(varname, region, tfile_vv, "VV", prefix=SYSTVAR)
def looseel_vv(varname, region, SYSTVAR=""): return looseel_file(varname, region, tfile_vv, "VV", prefix=SYSTVAR)

def tightmu_qcd(varname, region, SYSTVAR=""): return tightmu_file(varname, region, tfile_qcd_mu, "QCD", prefix=SYSTVAR)
def loosemu_qcd(varname, region, SYSTVAR=""): return loosemu_file(varname, region, tfile_qcd_mu, "QCD", prefix=SYSTVAR)
def tightel_qcd(varname, region, SYSTVAR=""): return tightel_file(varname, region, tfile_qcd_mu, "QCD", prefix=SYSTVAR)
def looseel_qcd(varname, region, SYSTVAR=""): return looseel_file(varname, region, tfile_qcd_mu, "QCD", prefix=SYSTVAR)

def tightel_qcd(varname, region, SYSTVAR=""): return tightel_file(varname, region, tfile_qcd_el, "QCD", prefix=SYSTVAR)
def looseel_qcd(varname, region, SYSTVAR=""): return looseel_file(varname, region, tfile_qcd_el, "QCD", prefix=SYSTVAR)
def tightel_qcd(varname, region, SYSTVAR=""): return tightel_file(varname, region, tfile_qcd_el, "QCD", prefix=SYSTVAR)
def looseel_qcd(varname, region, SYSTVAR=""): return looseel_file(varname, region, tfile_qcd_el, "QCD", prefix=SYSTVAR)

###################################################################################################
#
#
# Function that plots the variable with "Data" and "Bkgs". Also computes "Data - Bkg" and returns Ewk subtracted data histogram.
#
#
###################################################################################################
# 1D version
def plot(varname, region, t, option={}, addqcd=False, nfs=[], syst=""):
    fulloptions = option
    bkgs = ["wj", "dy", "ttbar", "vv"]
    if addqcd:
        bkgs.append("qcd")
    for bkg in bkgs:
        cmd = "h_{b} = p.apply_nf_w_error({lep}_{b}(varname, region, syst), nfs)".format(lep=t, b=bkg)
        exec cmd
    exec "h_data = {lep}_data(varname, region, syst)".format(lep=t)
    h_totalbkg = h_wj.Clone("bkg")
    h_totalbkg.Reset()
    for bkg in bkgs:
        exec "h_totalbkg.Add(h_{b})".format(b=bkg)
    maxbkg = h_totalbkg.GetMaximum()
    minbkg = h_totalbkg.GetMaximum()
    maxdata = h_data.GetMaximum()
    mindata = h_data.GetMaximum()
    ymax = max(maxbkg, maxdata)
    ymin = min(minbkg, mindata)
    yrange = [ymin, ymax]
    if "yaxis_log" in option and option["yaxis_log"]:
        yrange = [ymin / 10000., ymax * 100.]
    if varname.find("varbin") != -1:
        option["divide_by_bin_width"] = True
    option["yaxis_range"] = yrange
    h_ratio = h_data.Clone("ratio")
    h_ratio.Add(h_totalbkg, -1)
    if addqcd:
        plotcmd = "p.plot_hist( sigs = [], bgs = [h_wj, h_dy, h_ttbar, h_vv, h_qcd], data = h_data, colors = [2001, 2003, 7004, 7005, 2005], options = fulloptions)".format(lep=t)
    else:
        plotcmd = "p.plot_hist( sigs = [], bgs = [h_wj, h_dy, h_ttbar, h_vv], data = h_data, colors = [2001, 2003, 7004, 7005], options = fulloptions)".format(lep=t)
    exec plotcmd
    return h_ratio

# 2D version
def plot2d(varname, region, t, option={}, addqcd=False, nfs=[], syst=""):
    fulloptions = option
    bkgs = ["wj", "dy", "ttbar", "vv"]
    if addqcd:
        bkgs.append("qcd")
    for bkg in bkgs:
        cmd = "h_{b} = p.apply_nf_w_error_2d({lep}_{b}(varname, region, syst), nfs)".format(lep=t, b=bkg)
        exec cmd
    exec "h_data = {lep}_data(varname, region, syst)".format(lep=t)
    h_totalbkg = h_wj.Clone("bkg")
    h_totalbkg.Reset()
    for bkg in bkgs:
        exec "h_totalbkg.Add(h_{b})".format(b=bkg)
    #if addqcd:
    #    plotcmd = "p.plot_hist( sigs = [], bgs = [h_wj, h_dy, h_ttbar, h_vv, h_qcd], data = h_data, colors = [2001, 2003, 7004, 7005, 2005], options = fulloptions)".format(lep=t)
    #else:
    #    plotcmd = "p.plot_hist( sigs = [], bgs = [h_wj, h_dy, h_ttbar, h_vv], data = h_data, colors = [2001, 2003, 7004, 7005], options = fulloptions)".format(lep=t)
    #exec plotcmd
    h_ratio = h_data.Clone("ratio")
    h_ratio.Add(h_totalbkg, -1)
    return h_ratio

###################################################################################################
#
#
# Function to compute normalization factors for the EWK subtraction.
#
#
###################################################################################################
# "CR" : the traditional CR of, high MT-window, high MET, NF in bins of conecorrpt
def nfs_from_CR_mu(syst):
    h_wj    = tightmu_wj   ("conecorrptvarbin", "CR", syst)
    h_dy    = tightmu_dy   ("conecorrptvarbin", "CR", syst)
    h_ttbar = tightmu_ttbar("conecorrptvarbin", "CR", syst)
    h_vv    = tightmu_vv   ("conecorrptvarbin", "CR", syst)
    h_ratio = tightmu_data ("conecorrptvarbin", "CR", syst).Clone("ratio")
    h_bkg   = h_wj.Clone   ("ratio")
    h_bkg.Add(h_dy)
    h_bkg.Add(h_ttbar)
    h_bkg.Add(h_vv)
    h_ratio.Divide(h_bkg)
    nfs = [ [ h_ratio.GetBinContent(i), h_ratio.GetBinError(i) ] for i in xrange(1, h_ratio.GetNbinsX()+1) ]
    return nfs

def nfs_from_CR_el(syst):
    h_wj    = tightel_wj   ("conecorrptvarbin", "CR", syst)
    h_dy    = tightel_dy   ("conecorrptvarbin", "CR", syst)
    h_ttbar = tightel_ttbar("conecorrptvarbin", "CR", syst)
    h_vv    = tightel_vv   ("conecorrptvarbin", "CR", syst)
    h_ratio = tightel_data ("conecorrptvarbin", "CR", syst).Clone("ratio")
    h_bkg   = h_wj.Clone   ("ratio")
    h_bkg.Add(h_dy)
    h_bkg.Add(h_ttbar)
    h_bkg.Add(h_vv)
    h_ratio.Divide(h_bkg)
    nfs = [ [ h_ratio.GetBinContent(i), h_ratio.GetBinError(i) ] for i in xrange(1, h_ratio.GetNbinsX()+1) ]
    return nfs

# "CR3" : MET < 20. MT > 60, single NF
def nfs_from_CR3_mu(syst):
    h_wj    = tightmu_wj   ("mt", "CR3", syst)
    h_dy    = tightmu_dy   ("mt", "CR3", syst)
    h_ttbar = tightmu_ttbar("mt", "CR3", syst)
    h_vv    = tightmu_vv   ("mt", "CR3", syst)
    h_ratio = tightmu_data ("mt", "CR3", syst).Clone("ratio")
    h_ratio.Rebin(50)
    h_bkg   = h_wj.Clone   ("ratio")
    h_bkg.Add(h_dy)
    h_bkg.Add(h_ttbar)
    h_bkg.Add(h_vv)
    h_bkg.Rebin(50)
    h_ratio.Divide(h_bkg)
    print h_ratio.GetBinContent(1)
    print h_ratio.GetBinError(1)
    nfs = [ [ h_ratio.GetBinContent(i), h_ratio.GetBinError(i) ] for i in xrange(1, h_ratio.GetNbinsX()+1) ]
    return nfs

def nfs_from_CR3_el(syst):
    h_wj    = tightel_wj   ("mt", "CR3", syst)
    h_dy    = tightel_dy   ("mt", "CR3", syst)
    h_ttbar = tightel_ttbar("mt", "CR3", syst)
    h_vv    = tightel_vv   ("mt", "CR3", syst)
    h_ratio = tightel_data ("mt", "CR3", syst).Clone("ratio")
    h_ratio.Rebin(50)
    h_bkg   = h_wj.Clone   ("ratio")
    h_bkg.Add(h_dy)
    h_bkg.Add(h_ttbar)
    h_bkg.Add(h_vv)
    h_bkg.Rebin(50)
    h_ratio.Divide(h_bkg)
    print h_ratio.GetBinContent(1)
    print h_ratio.GetBinError(1)
    nfs = [ [ h_ratio.GetBinContent(i), h_ratio.GetBinError(i) ] for i in xrange(1, h_ratio.GetNbinsX()+1) ]
    return nfs

# "CR2" : MET < 20, MT (80, 120), single NF
def nfs_from_CR2_mu(syst):
    h_wj    = tightmu_wj   ("mt", "CR2", syst)
    h_dy    = tightmu_dy   ("mt", "CR2", syst)
    h_ttbar = tightmu_ttbar("mt", "CR2", syst)
    h_vv    = tightmu_vv   ("mt", "CR2", syst)
    h_ratio = tightmu_data ("mt", "CR2", syst).Clone("ratio")
    h_ratio.Rebin(50)
    h_bkg   = h_wj.Clone   ("ratio")
    h_bkg.Add(h_dy)
    h_bkg.Add(h_ttbar)
    h_bkg.Add(h_vv)
    h_bkg.Rebin(50)
    h_ratio.Divide(h_bkg)
    print h_ratio.GetBinContent(1)
    print h_ratio.GetBinError(1)
    nfs = [ [ h_ratio.GetBinContent(i), h_ratio.GetBinError(i) ] for i in xrange(1, h_ratio.GetNbinsX()+1) ]
    return nfs

def nfs_from_CR2_el(syst):
    h_wj    = tightel_wj   ("mt", "CR2", syst)
    h_dy    = tightel_dy   ("mt", "CR2", syst)
    h_ttbar = tightel_ttbar("mt", "CR2", syst)
    h_vv    = tightel_vv   ("mt", "CR2", syst)
    h_ratio = tightel_data ("mt", "CR2", syst).Clone("ratio")
    h_ratio.Rebin(50)
    h_bkg   = h_wj.Clone   ("ratio")
    h_bkg.Add(h_dy)
    h_bkg.Add(h_ttbar)
    h_bkg.Add(h_vv)
    h_bkg.Rebin(50)
    h_ratio.Divide(h_bkg)
    print h_ratio.GetBinContent(1)
    print h_ratio.GetBinError(1)
    nfs = [ [ h_ratio.GetBinContent(i), h_ratio.GetBinError(i) ] for i in xrange(1, h_ratio.GetNbinsX()+1) ]
    return nfs

###################################################################################################
#
#
# Functions to return the final fakerate histograms.
#
#
###################################################################################################
def fakerate_1d_mu_data_hist(syst=""):
    h_mu_tight = plot("conecorrptvarbin", "MR", "tightmu", {"output_name": "frplots/plot_tightmu.png", "yaxis_log":True}, False, nfs_from_CR2_mu(syst)).Clone("tight")
    h_mu_loose = plot("conecorrptvarbin", "MR", "loosemu", {"output_name": "frplots/plot_loosemu.png", "yaxis_log":True}, False, nfs_from_CR2_mu(syst)).Clone("loose")
    h_mu_tight.Divide(h_mu_loose)
    return h_mu_tight

def fakerate_1d_el_data_hist(syst=""):
    h_el_tight = plot("conecorrptvarbin", "MR", "tightel", {"output_name": "frplots/plot_tightel.png", "yaxis_log":True}, False, nfs_from_CR2_el(syst)).Clone("tight")
    h_el_loose = plot("conecorrptvarbin", "MR", "looseel", {"output_name": "frplots/plot_looseel.png", "yaxis_log":True}, False, nfs_from_CR2_el(syst)).Clone("loose")
    h_el_tight.Divide(h_el_loose)
    return h_el_tight

def fakerate_1d_mu_qcd_hist(syst=""):
    h_mu_tight_qcd = tightmu_qcd("conecorrptvarbin", "MR")
    h_mu_loose_qcd = loosemu_qcd("conecorrptvarbin", "MR")
    h_mu_tight_qcd.Divide(h_mu_loose_qcd)
    return h_mu_tight_qcd

def fakerate_1d_el_qcd_hist(syst=""):
    h_el_tight_qcd = tightel_qcd("conecorrptvarbin", "MR")
    h_el_loose_qcd = looseel_qcd("conecorrptvarbin", "MR")
    h_el_tight_qcd.Divide(h_el_loose_qcd)
    return h_el_tight_qcd

def fakerate_2d_mu_data_hist(syst=""):
    h_mu_tight = plot2d("conecorrpt_v_eta_muvarbin", "MR", "tightmu", {"output_name": "frplots/plot2d_tightmu.png"}, False, nfs_from_CR2_mu(syst)).Clone("tight")
    h_mu_loose = plot2d("conecorrpt_v_eta_muvarbin", "MR", "loosemu", {"output_name": "frplots/plot2d_loosemu.png"}, False, nfs_from_CR2_mu(syst)).Clone("loose")
    h_mu_tight.Divide(h_mu_loose)
    return h_mu_tight

def fakerate_2d_el_data_hist(syst=""):
    h_el_tight = plot2d("conecorrpt_v_eta_elvarbin", "MR", "tightel", {"output_name": "frplots/plot2d_tightel.png"}, False, nfs_from_CR2_el(syst)).Clone("tight")
    h_el_loose = plot2d("conecorrpt_v_eta_elvarbin", "MR", "looseel", {"output_name": "frplots/plot2d_looseel.png"}, False, nfs_from_CR2_el(syst)).Clone("loose")
    h_el_tight.Divide(h_el_loose)
    return h_el_tight

def fakerate_2d_mu_qcd_hist(syst=""):
    h_mu_tight_qcd = tightmu_qcd("conecorrpt_v_eta_muvarbin", "MR")
    h_mu_loose_qcd = loosemu_qcd("conecorrpt_v_eta_muvarbin", "MR")
    h_mu_tight_qcd.Divide(h_mu_loose_qcd)
    return h_mu_tight_qcd

def fakerate_2d_el_qcd_hist(syst=""):
    h_el_tight_qcd = tightel_qcd("conecorrpt_v_eta_elvarbin", "MR")
    h_el_loose_qcd = looseel_qcd("conecorrpt_v_eta_elvarbin", "MR")
    h_el_tight_qcd.Divide(h_el_loose_qcd)
    return h_el_tight_qcd

###################################################################################################
#
#
# Utility functions
#
#
###################################################################################################
def get_full_error(nominal, systup, systdn):
    diffup = nominal.Clone("diffup")
    diffdn = nominal.Clone("diffdn")
    diffup.Add(systup, -1)
    diffdn.Add(systdn, -1)
    for i in xrange(0, nominal.GetNbinsX()+2):
        dup = abs(diffup.GetBinContent(i))
        ddn = abs(diffdn.GetBinContent(i))
        d = max(dup, ddn)
        be = nominal.GetBinError(i)
        nominal.SetBinError(i, math.sqrt(be**2 + d**2))
    return nominal

def get_full_error_2d(nominal, systup, systdn):
    diffup = nominal.Clone("diffup")
    diffdn = nominal.Clone("diffdn")
    diffup.Add(systup, -1)
    diffdn.Add(systdn, -1)
    for i in xrange(0, nominal.GetNbinsX()+2):
        for j in xrange(0, nominal.GetNbinsY()+2):
            dup = abs(diffup.GetBinContent(i, j))
            ddn = abs(diffdn.GetBinContent(i, j))
            d = max(dup, ddn)
            be = nominal.GetBinError(i, j)
            nominal.SetBinError(i, j, math.sqrt(be**2 + d**2))
    return nominal

###################################################################################################
#
#
# Finally, the scripts to actually perform the jobs.
#
#
###################################################################################################
def draw_fakerate_1d_mu():
    h_fakerate_1d_mu_data = fakerate_1d_mu_data_hist().Clone("h_fakerate_1d_mu_data")
    h_fakerate_1d_mu_data.Print("all")
    h_fakerate_1d_mu_data_syst13 = fakerate_1d_mu_data_hist("syst13_").Clone("h_fakerate_1d_mu_data")
    h_fakerate_1d_mu_data_syst13.Print("all")
    h_fakerate_1d_mu_data_syst14 = fakerate_1d_mu_data_hist("syst14_").Clone("h_fakerate_1d_mu_data")
    h_fakerate_1d_mu_data_syst14.Print("all")
    h_fakerate_1d_mu_data_fullerror = get_full_error(h_fakerate_1d_mu_data, h_fakerate_1d_mu_data_syst13, h_fakerate_1d_mu_data_syst14)
    h_fakerate_1d_mu_data_fullerror.Print("all")
    h_fakerate_1d_mu_qcd  = fakerate_1d_mu_qcd_hist().Clone("QCD")
    h_fakerate_1d_mu_qcd .Print("all")
    p.plot_hist(sigs=[], bgs=[h_fakerate_1d_mu_qcd], data=h_fakerate_1d_mu_data_fullerror, colors=[2], options={"output_name":"frplots/fakerate_mu.png", "ratio_range":[0.0, 2.0], "draw_points":True})

def draw_fakerate_1d_el():
    h_fakerate_1d_el_data = fakerate_1d_el_data_hist().Clone("h_fakerate_1d_el_data")
    h_fakerate_1d_el_data.Print("all")
    h_fakerate_1d_el_data_syst13 = fakerate_1d_el_data_hist("syst13_").Clone("h_fakerate_1d_el_data")
    h_fakerate_1d_el_data_syst13.Print("all")
    h_fakerate_1d_el_data_syst14 = fakerate_1d_el_data_hist("syst14_").Clone("h_fakerate_1d_el_data")
    h_fakerate_1d_el_data_syst14.Print("all")
    h_fakerate_1d_el_data_fullerror = get_full_error(h_fakerate_1d_el_data, h_fakerate_1d_el_data_syst13, h_fakerate_1d_el_data_syst14)
    h_fakerate_1d_el_data_fullerror.Print("all")
    h_fakerate_1d_el_qcd  = fakerate_1d_el_qcd_hist().Clone("QCD")
    h_fakerate_1d_el_qcd .Print("all")
    p.plot_hist(sigs=[], bgs=[h_fakerate_1d_el_qcd], data=h_fakerate_1d_el_data_fullerror, colors=[2], options={"output_name":"frplots/fakerate_el.png", "ratio_range":[0.0, 2.0], "draw_points":True})

def draw_fakerate_2d_mu():
    h_fakerate_2d_mu_data = fakerate_2d_mu_data_hist().Clone("h_fakerate_2d_mu_data")
    h_fakerate_2d_mu_data.Print("all")
    h_fakerate_2d_mu_data_syst13 = fakerate_2d_mu_data_hist("syst13_").Clone("h_fakerate_2d_mu_data")
    h_fakerate_2d_mu_data_syst13.Print("all")
    h_fakerate_2d_mu_data_syst14 = fakerate_2d_mu_data_hist("syst14_").Clone("h_fakerate_2d_mu_data")
    h_fakerate_2d_mu_data_syst14.Print("all")
    h_fakerate_2d_mu_data_fullerror = get_full_error_2d(h_fakerate_2d_mu_data, h_fakerate_2d_mu_data_syst13, h_fakerate_2d_mu_data_syst14)
    h_fakerate_2d_mu_data_fullerror.Print("all")
    h_fakerate_2d_mu_qcd  = fakerate_2d_mu_qcd_hist().Clone("QCD")
    h_fakerate_2d_mu_qcd .Print("all")
    max_data = h_fakerate_2d_mu_data_fullerror.GetMaximum()
    min_data = h_fakerate_2d_mu_data_fullerror.GetMinimum()
    max_qcd = h_fakerate_2d_mu_qcd.GetMaximum()
    min_qcd = h_fakerate_2d_mu_qcd.GetMinimum()
    ply.plot_hist_2d( h_fakerate_2d_mu_data_fullerror, options = { "output_name": "frplots/fakerate_2d_mu_data.png", "zaxis_range": [min_data/1.5, 1.5*max_data], "zaxis_log": False, "bin_text_smart": False, "us_flag": False, "output_ic": False, "zaxis_noexponents": True, "draw_option_2d": "textecolz", "bin_text_format": ".3f", "xaxis_log": True, "bin_text_size": 1.0, "palette_name": "radiation" })
    ply.plot_hist_2d( h_fakerate_2d_mu_qcd           , options = { "output_name": "frplots/fakerate_2d_mu_qcd.png" , "zaxis_range": [min_qcd /1.5, 1.5*max_qcd ], "zaxis_log": False, "bin_text_smart": False, "us_flag": False, "output_ic": False, "zaxis_noexponents": True, "draw_option_2d": "textecolz", "bin_text_format": ".3f", "xaxis_log": True, "bin_text_size": 1.0, "palette_name": "radiation" })

def draw_fakerate_2d_el():
    h_fakerate_2d_el_data = fakerate_2d_el_data_hist().Clone("h_fakerate_2d_el_data")
    h_fakerate_2d_el_data.Print("all")
    h_fakerate_2d_el_data_syst13 = fakerate_2d_el_data_hist("syst13_").Clone("h_fakerate_2d_el_data")
    h_fakerate_2d_el_data_syst13.Print("all")
    h_fakerate_2d_el_data_syst14 = fakerate_2d_el_data_hist("syst14_").Clone("h_fakerate_2d_el_data")
    h_fakerate_2d_el_data_syst14.Print("all")
    h_fakerate_2d_el_data_fullerror = get_full_error_2d(h_fakerate_2d_el_data, h_fakerate_2d_el_data_syst13, h_fakerate_2d_el_data_syst14)
    h_fakerate_2d_el_data_fullerror.Print("all")
    h_fakerate_2d_el_qcd  = fakerate_2d_el_qcd_hist().Clone("QCD")
    h_fakerate_2d_el_qcd .Print("all")
    max_data = h_fakerate_2d_el_data_fullerror.GetMaximum()
    min_data = h_fakerate_2d_el_data_fullerror.GetMinimum()
    max_qcd = h_fakerate_2d_el_qcd.GetMaximum()
    min_qcd = h_fakerate_2d_el_qcd.GetMinimum()
    ply.plot_hist_2d( h_fakerate_2d_el_data_fullerror, options = { "output_name": "frplots/fakerate_2d_el_data.png", "zaxis_range": [min_data/1.5, 1.5*max_data], "zaxis_log": False, "bin_text_smart": False, "us_flag": False, "output_ic": False, "zaxis_noexponents": True, "draw_option_2d": "textecolz", "bin_text_format": ".3f", "xaxis_log": True, "bin_text_size": 1.0, "palette_name": "radiation" })
    ply.plot_hist_2d( h_fakerate_2d_el_qcd           , options = { "output_name": "frplots/fakerate_2d_el_qcd.png" , "zaxis_range": [min_qcd /1.5, 1.5*max_qcd ], "zaxis_log": False, "bin_text_smart": False, "us_flag": False, "output_ic": False, "zaxis_noexponents": True, "draw_option_2d": "textecolz", "bin_text_format": ".3f", "xaxis_log": True, "bin_text_size": 1.0, "palette_name": "radiation" })

if __name__ == "__main__":
    draw_fakerate_1d_mu()
    draw_fakerate_1d_el()
    draw_fakerate_2d_mu()
    draw_fakerate_2d_el()


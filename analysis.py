#!/bin/env python

import sys
sys.path.append("/home/users/phchang/public_html/phys/rooutil")
import plotmaker
import ROOT as r
from rooutil import plottery_wrapper as p
from plottery import plottery as ply
import math

r.gStyle.SetPalette( r.kInvertedDarkBodyRadiator )

try:
    output_directory = sys.argv[1]
except:
    output_directory = "outputs"

###################################################################################################
# List of root files with histograms after running ScanChain.C
###################################################################################################
tfile_data   = r.TFile("{}/fakerate_data.root".format(output_directory))
tfile_wj     = r.TFile("{}/fakerate_wj.root".format(output_directory))
tfile_dy     = r.TFile("{}/fakerate_dy.root".format(output_directory))
tfile_ttbar  = r.TFile("{}/fakerate_ttbar.root".format(output_directory))
tfile_vv     = r.TFile("{}/fakerate_vv.root".format(output_directory))
tfile_qcd_mu = r.TFile("{}/fakerate_qcd_mu.root".format(output_directory))
tfile_qcd_el = r.TFile("{}/fakerate_qcd_el.root".format(output_directory))

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
        cmd = "h_{b} = p.apply_nf({lep}_{b}(varname, region, syst), nfs)".format(lep=t, b=bkg)
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
    yrange = [0, ymax*2.]
    if "yaxis_log" in option and option["yaxis_log"]:
        yrange = [ymin / 10000., ymax * 100.]
        option["yaxis_range"] = yrange
    if varname.find("varbin") != -1:
        option["divide_by_bin_width"] = True
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
    categs = ["wj", "dy", "ttbar", "vv", "data"]
    #if addqcd:
    #    bkgs.append("qcd")
    for bkg in bkgs:
        cmd = "h_{b} = p.apply_nf_2d({lep}_{b}(varname, region, syst), nfs)".format(lep=t, b=bkg)
        exec cmd
    # Special rebinning
    exec "h_data = {lep}_data(varname, region, syst)".format(lep=t)
    #for c in categs:
    #    exec "bc = h_{b}.GetBinContent(4, 2)".format(b=c)
    #    exec "be = h_{b}.GetBinError(4, 2)".format(b=c)
    #    exec "h_{b}.SetBinContent(4, 3, bc)".format(b=c)
    #    exec "h_{b}.SetBinError(4, 3, be)".format(b=c)
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
def nfs_from_CR_mu(syst, inclqcd):
    h_wj    = tightmu_wj   ("conecorrptvarbin", "CR", syst)
    h_dy    = tightmu_dy   ("conecorrptvarbin", "CR", syst)
    h_ttbar = tightmu_ttbar("conecorrptvarbin", "CR", syst)
    h_vv    = tightmu_vv   ("conecorrptvarbin", "CR", syst)
    h_qcd   = tightmu_qcd  ("conecorrptvarbin", "CR", "")
    h_data  = tightmu_data ("conecorrptvarbin", "CR", syst).Clone("Data")
    h_ratio = tightmu_data ("conecorrptvarbin", "CR", syst).Clone("ratio")
    h_ratio.Rebin(4)
    h_bkg   = h_wj.Clone   ("ratio")
    h_bkg.Add(h_dy)
    h_bkg.Add(h_ttbar)
    h_bkg.Add(h_vv)
    if inclqcd: h_ratio.Add(h_qcd, -1)
    h_bkg.Rebin(4)
    h_ratio.Divide(h_bkg)
    nfs = [ [ h_ratio.GetBinContent(i), h_ratio.GetBinError(i) ] for i in xrange(1, h_ratio.GetNbinsX()+1) ]
    if syst == "":
        p.plot_hist( sigs = [], bgs = [h_wj.Clone("W"), h_dy.Clone("Z"), h_ttbar.Clone("Top"), h_vv.Clone("VV")], data = h_data.Clone("Data"), colors = [2001, 2003, 7004, 7005], options = {"output_name":"frplots/plot_cr_tightmu_conecorrpt_noqcd.png", "divide_by_bin_width":True, "xaxis_log":True, "ratio_xaxis_title":"#it{p}_{T,cone-corr} [GeV]"})
        p.plot_hist( sigs = [], bgs = [h_wj.Clone("W"), h_dy.Clone("Z"), h_ttbar.Clone("Top"), h_vv.Clone("VV"), h_qcd.Clone("QCD")], data = h_data.Clone("Data"), colors = [2001, 2003, 7004, 7005, 2005], options = {"output_name":"frplots/plot_cr_tightmu_conecorrpt.png", "divide_by_bin_width":True, "xaxis_log":True, "ratio_xaxis_title":"#it{p}_{T,cone-corr} [GeV]"})
    return nfs

def nfs_from_CR_el(syst, inclqcd):
    h_wj    = tightel_wj   ("conecorrptvarbin", "CR", syst)
    h_dy    = tightel_dy   ("conecorrptvarbin", "CR", syst)
    h_ttbar = tightel_ttbar("conecorrptvarbin", "CR", syst)
    h_vv    = tightel_vv   ("conecorrptvarbin", "CR", syst)
    h_qcd   = tightel_qcd  ("conecorrptvarbin", "CR", "")
    h_data  = tightel_data ("conecorrptvarbin", "CR", syst).Clone("Data")
    h_ratio = tightel_data ("conecorrptvarbin", "CR", syst).Clone("ratio")
    h_ratio.Rebin(4)
    h_bkg   = h_wj.Clone   ("ratio")
    h_bkg.Add(h_dy)
    h_bkg.Add(h_ttbar)
    h_bkg.Add(h_vv)
    if inclqcd: h_ratio.Add(h_qcd, -1)
    h_bkg.Rebin(4)
    h_ratio.Divide(h_bkg)
    nfs = [ [ h_ratio.GetBinContent(i), h_ratio.GetBinError(i) ] for i in xrange(1, h_ratio.GetNbinsX()+1) ]
    if syst == "":
        p.plot_hist( sigs = [], bgs = [h_wj.Clone("W"), h_dy.Clone("Z"), h_ttbar.Clone("Top"), h_vv.Clone("VV")], data = h_data.Clone("Data"), colors = [2001, 2003, 7004, 7005], options = {"output_name":"frplots/plot_cr_tightel_conecorrpt_noqcd.png", "divide_by_bin_width":True, "xaxis_log":True, "ratio_xaxis_title":"#it{p}_{T,cone-corr} [GeV]"})
        p.plot_hist( sigs = [], bgs = [h_wj.Clone("W"), h_dy.Clone("Z"), h_ttbar.Clone("Top"), h_vv.Clone("VV"), h_qcd.Clone("QCD")], data = h_data.Clone("Data"), colors = [2001, 2003, 7004, 7005, 2005], options = {"output_name":"frplots/plot_cr_tightel_conecorrpt.png", "divide_by_bin_width":True, "xaxis_log":True, "ratio_xaxis_title":"#it{p}_{T,cone-corr} [GeV]"})
    return nfs

def nfs_from_CR_mu_ptbin(syst, inclqcd):
    h_wj    = tightmu_wj   ("conecorrptvarbin", "CR", syst)
    h_dy    = tightmu_dy   ("conecorrptvarbin", "CR", syst)
    h_ttbar = tightmu_ttbar("conecorrptvarbin", "CR", syst)
    h_vv    = tightmu_vv   ("conecorrptvarbin", "CR", syst)
    h_qcd   = tightmu_qcd  ("conecorrptvarbin", "CR", "")
    h_data  = tightmu_data ("conecorrptvarbin", "CR", syst).Clone("Data")
    h_ratio = tightmu_data ("conecorrptvarbin", "CR", syst).Clone("ratio")
    h_bkg   = h_wj.Clone   ("ratio")
    h_bkg.Add(h_dy)
    h_bkg.Add(h_ttbar)
    h_bkg.Add(h_vv)
    if inclqcd: h_ratio.Add(h_qcd, -1)
    h_ratio.Divide(h_bkg)
    nfs = [ [ h_ratio.GetBinContent(i), h_ratio.GetBinError(i) ] for i in xrange(1, h_ratio.GetNbinsX()+1) ]
    if syst == "":
        p.plot_hist( sigs = [], bgs = [h_wj.Clone("W"), h_dy.Clone("Z"), h_ttbar.Clone("Top"), h_vv.Clone("VV")], data = h_data.Clone("Data"), colors = [2001, 2003, 7004, 7005], options = {"output_name":"frplots/plot_cr_tightmu_conecorrpt_noqcd.png", "divide_by_bin_width":True, "xaxis_log":True, "ratio_xaxis_title":"#it{p}_{T,cone-corr} [GeV]"})
        p.plot_hist( sigs = [], bgs = [h_wj.Clone("W"), h_dy.Clone("Z"), h_ttbar.Clone("Top"), h_vv.Clone("VV"), h_qcd.Clone("QCD")], data = h_data.Clone("Data"), colors = [2001, 2003, 7004, 7005, 2005], options = {"output_name":"frplots/plot_cr_tightmu_conecorrpt.png", "divide_by_bin_width":True, "xaxis_log":True, "ratio_xaxis_title":"#it{p}_{T,cone-corr} [GeV]"})
    return nfs

def nfs_from_CR_el_ptbin(syst, inclqcd):
    h_wj    = tightel_wj   ("conecorrptvarbin", "CR", syst)
    h_dy    = tightel_dy   ("conecorrptvarbin", "CR", syst)
    h_ttbar = tightel_ttbar("conecorrptvarbin", "CR", syst)
    h_vv    = tightel_vv   ("conecorrptvarbin", "CR", syst)
    h_qcd   = tightel_qcd  ("conecorrptvarbin", "CR", "")
    h_data  = tightel_data ("conecorrptvarbin", "CR", syst).Clone("Data")
    h_ratio = tightel_data ("conecorrptvarbin", "CR", syst).Clone("ratio")
    h_bkg   = h_wj.Clone   ("ratio")
    h_bkg.Add(h_dy)
    h_bkg.Add(h_ttbar)
    h_bkg.Add(h_vv)
    if inclqcd: h_ratio.Add(h_qcd, -1)
    h_ratio.Divide(h_bkg)
    nfs = [ [ h_ratio.GetBinContent(i), h_ratio.GetBinError(i) ] for i in xrange(1, h_ratio.GetNbinsX()+1) ]
    if syst == "":
        p.plot_hist( sigs = [], bgs = [h_wj.Clone("W"), h_dy.Clone("Z"), h_ttbar.Clone("Top"), h_vv.Clone("VV")], data = h_data.Clone("Data"), colors = [2001, 2003, 7004, 7005], options = {"output_name":"frplots/plot_cr_tightel_conecorrpt_noqcd.png", "divide_by_bin_width":True, "xaxis_log":True, "ratio_xaxis_title":"#it{p}_{T,cone-corr} [GeV]"})
        p.plot_hist( sigs = [], bgs = [h_wj.Clone("W"), h_dy.Clone("Z"), h_ttbar.Clone("Top"), h_vv.Clone("VV"), h_qcd.Clone("QCD")], data = h_data.Clone("Data"), colors = [2001, 2003, 7004, 7005, 2005], options = {"output_name":"frplots/plot_cr_tightel_conecorrpt.png", "divide_by_bin_width":True, "xaxis_log":True, "ratio_xaxis_title":"#it{p}_{T,cone-corr} [GeV]"})
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
    #print h_ratio.GetBinContent(1)
    #print h_ratio.GetBinError(1)
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
    #print h_ratio.GetBinContent(1)
    #print h_ratio.GetBinError(1)
    nfs = [ [ h_ratio.GetBinContent(i), h_ratio.GetBinError(i) ] for i in xrange(1, h_ratio.GetNbinsX()+1) ]
    return nfs

# "CR2" : MET < 20, MT (80, 120), single NF
def nfs_from_CR2_mu(syst, inclqcd):
    h_wj    = tightmu_wj   ("mt", "CR2", syst)
    h_dy    = tightmu_dy   ("mt", "CR2", syst)
    h_ttbar = tightmu_ttbar("mt", "CR2", syst)
    h_vv    = tightmu_vv   ("mt", "CR2", syst)
    h_qcd   = tightmu_qcd  ("mt", "CR2", "")
    h_ratio = tightmu_data ("mt", "CR2", syst).Clone("ratio")
    h_bkg   = h_wj.Clone   ("ratio")
    h_bkg.Add(h_dy)
    h_bkg.Add(h_ttbar)
    h_bkg.Add(h_vv)
    if inclqcd: h_ratio.Add(h_qcd, -1)
    h_ratio.Rebin(50)
    h_bkg.Rebin(50)
    h_ratio.Divide(h_bkg)
    #print h_ratio.GetBinContent(1)
    #print h_ratio.GetBinError(1)
    nfs = [ [ h_ratio.GetBinContent(i), h_ratio.GetBinError(i) ] for i in xrange(1, h_ratio.GetNbinsX()+1) ]
    return nfs

def plot_conecorrpt_CR2_el(syst, inclqcd):
    h_wj    = tightel_wj   ("conecorrptvarbin", "CR2", syst)
    h_dy    = tightel_dy   ("conecorrptvarbin", "CR2", syst)
    h_ttbar = tightel_ttbar("conecorrptvarbin", "CR2", syst)
    h_vv    = tightel_vv   ("conecorrptvarbin", "CR2", syst)
    h_qcd   = tightel_qcd  ("conecorrptvarbin", "CR2", "")
    h_data  = tightel_data ("conecorrptvarbin", "CR2", syst).Clone("Data")
    h_ratio = tightel_data ("conecorrptvarbin", "CR2", syst).Clone("ratio")
    h_bkg   = h_wj.Clone   ("ratio")
    h_bkg.Add(h_dy)
    h_bkg.Add(h_ttbar)
    h_bkg.Add(h_vv)
    if inclqcd: h_bkg.Add(h_qcd)
    h_ratio.Divide(h_bkg)
    nfs = [ [ h_ratio.GetBinContent(i), h_ratio.GetBinError(i) ] for i in xrange(1, h_ratio.GetNbinsX()+1) ]
    if syst == "":
        p.plot_hist( sigs = [], bgs = [h_wj, h_dy, h_ttbar, h_vv, h_qcd], data = h_data, colors = [2001, 2003, 7004, 7005, 2005], options = {"output_name":"frplots/plot_cr2_tightel_conecorrpt.png"})
        p.plot_hist( sigs = [], bgs = [h_wj, h_dy, h_ttbar, h_vv], data = h_data, colors = [2001, 2003, 7004, 7005], options = {"output_name":"frplots/plot_cr2_tightel_conecorrpt_noqcd.png"})
    return nfs

def nfs_from_CR2_el(syst, inclqcd):
    h_wj    = tightel_wj   ("mt", "CR2", syst)
    h_dy    = tightel_dy   ("mt", "CR2", syst)
    h_ttbar = tightel_ttbar("mt", "CR2", syst)
    h_vv    = tightel_vv   ("mt", "CR2", syst)
    h_qcd   = tightel_qcd  ("mt", "CR2", "")
    h_ratio = tightel_data ("mt", "CR2", syst).Clone("ratio")
    h_bkg   = h_wj.Clone   ("ratio")
    h_bkg.Add(h_dy)
    h_bkg.Add(h_ttbar)
    h_bkg.Add(h_vv)
    if inclqcd: h_ratio.Add(h_qcd, -1)
    h_ratio.Rebin(50)
    h_bkg.Rebin(50)
    h_ratio.Divide(h_bkg)
    #print h_ratio.GetBinContent(1)
    #print h_ratio.GetBinError(1)
    nfs = [ [ h_ratio.GetBinContent(i), h_ratio.GetBinError(i) ] for i in xrange(1, h_ratio.GetNbinsX()+1) ]
    plot_conecorrpt_CR2_el(syst, inclqcd)
    return nfs

###################################################################################################
#
#
# Functions to return the final fakerate histograms.
#
#
###################################################################################################
def fakerate_1d_mu_data_hist(syst="", nfscheme=nfs_from_CR_mu, nfsinclqcd=False):
    h_mu_tight = plot("conecorrptvarbin", "MR", "tightmu", {"output_name": "frplots/{}plot_tightmu.png".format(syst), "yaxis_log":True, "no_ratio":True, "xaxis_label":"#it{p}_{T,cone-corr} [GeV]", "yaxis_label":"Events / Bin Width", "hist_disable_xerrors":False}, False, nfscheme(syst, nfsinclqcd), syst).Clone("tight")
    h_mu_loose = plot("conecorrptvarbin", "MR", "loosemu", {"output_name": "frplots/{}plot_loosemu.png".format(syst), "yaxis_log":True, "no_ratio":True, "xaxis_label":"#it{p}_{T,cone-corr} [GeV]", "yaxis_label":"Events / Bin Width", "hist_disable_xerrors":False}, False, nfscheme(syst, nfsinclqcd), syst).Clone("loose")
    h_mu_tight.Divide(h_mu_loose)
    return h_mu_tight

def fakerate_1d_el_data_hist(syst="", nfscheme=nfs_from_CR_el, nfsinclqcd=False):
    h_el_tight = plot("conecorrptvarbin", "MR", "tightel", {"output_name": "frplots/{}plot_tightel.png".format(syst), "yaxis_log":True, "no_ratio":True, "xaxis_label":"#it{p}_{T,cone-corr} [GeV]", "yaxis_label":"Events / Bin Width", "hist_disable_xerrors":False}, False, nfscheme(syst, nfsinclqcd), syst).Clone("tight")
    h_el_loose = plot("conecorrptvarbin", "MR", "looseel", {"output_name": "frplots/{}plot_looseel.png".format(syst), "yaxis_log":True, "no_ratio":True, "xaxis_label":"#it{p}_{T,cone-corr} [GeV]", "yaxis_label":"Events / Bin Width", "hist_disable_xerrors":False}, False, nfscheme(syst, nfsinclqcd), syst).Clone("loose")
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

def fakerate_2d_mu_data_hist(syst="", nfscheme=nfs_from_CR_mu, nfsinclqcd=False):
    h_mu_tight = plot2d("conecorrpt_v_eta_muvarbin", "MR", "tightmu", {"output_name": "frplots/{}plot2d_tightmu.png".format(syst)}, False, nfscheme(syst, nfsinclqcd), syst).Clone("tight")
    h_mu_loose = plot2d("conecorrpt_v_eta_muvarbin", "MR", "loosemu", {"output_name": "frplots/{}plot2d_loosemu.png".format(syst)}, False, nfscheme(syst, nfsinclqcd), syst).Clone("loose")
    h_mu_tight.Divide(h_mu_loose)
    return h_mu_tight

def fakerate_2d_el_data_hist(syst="", nfscheme=nfs_from_CR_el, nfsinclqcd=False):
    h_el_tight = plot2d("conecorrpt_v_eta_elvarbin", "MR", "tightel", {"output_name": "frplots/{}plot2d_tightel.png".format(syst)}, False, nfscheme(syst, nfsinclqcd), syst).Clone("tight")
    h_el_loose = plot2d("conecorrpt_v_eta_elvarbin", "MR", "looseel", {"output_name": "frplots/{}plot2d_looseel.png".format(syst)}, False, nfscheme(syst, nfsinclqcd), syst).Clone("loose")
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
        #d = math.sqrt(dup**2 + ddn**2)
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
def draw_fakerate_1d_mu(nfscheme=nfs_from_CR_mu, nfsinclqcd=False):
    h_fakerate_1d_mu_data = fakerate_1d_mu_data_hist(nfscheme=nfscheme, nfsinclqcd=nfsinclqcd).Clone("h_fakerate_1d_mu_data")
    #h_fakerate_1d_mu_data.Print("all")
    h_fakerate_1d_mu_data_syst13 = fakerate_1d_mu_data_hist("syst13_", nfscheme=nfscheme, nfsinclqcd=nfsinclqcd).Clone("h_fakerate_1d_mu_data")
    #h_fakerate_1d_mu_data_syst13.Print("all")
    h_fakerate_1d_mu_data_syst14 = fakerate_1d_mu_data_hist("syst14_", nfscheme=nfscheme, nfsinclqcd=nfsinclqcd).Clone("h_fakerate_1d_mu_data")
    #h_fakerate_1d_mu_data_syst14.Print("all")
    h_fakerate_1d_mu_data_fullerror = get_full_error(h_fakerate_1d_mu_data, h_fakerate_1d_mu_data_syst13, h_fakerate_1d_mu_data_syst14)
    #h_fakerate_1d_mu_data_fullerror.Print("all")
    h_fakerate_1d_mu_qcd  = fakerate_1d_mu_qcd_hist().Clone("QCD")
    #h_fakerate_1d_mu_qcd .Print("all")
    p.plot_hist(sigs=[], bgs=[h_fakerate_1d_mu_qcd.Clone("QCD")], data=h_fakerate_1d_mu_data_fullerror, colors=[2], options={"output_name":"frplots/fakerate_mu.png" , "ratio_range":[0.0, 2.0], "draw_points":True, "yaxis_label":"Fake Rate", "no_ratio":True, "xaxis_label":"#it{p}_{T,cone-corr} [GeV]", "ymax_scale":0.65})
    p.plot_hist(sigs=[], bgs=[h_fakerate_1d_mu_qcd.Clone("QCD")], data=h_fakerate_1d_mu_data_syst13, colors=[2], options={"output_name":"frplots/syst13_fakerate_mu.png", "ratio_range":[0.0, 2.0], "draw_points":True, "yaxis_label":"Fake Rate", "no_ratio":True, "xaxis_label":"#it{p}_{T,cone-corr} [GeV]", "ymax_scale":0.65})
    p.plot_hist(sigs=[], bgs=[h_fakerate_1d_mu_qcd.Clone("QCD")], data=h_fakerate_1d_mu_data_syst14, colors=[2], options={"output_name":"frplots/syst14_fakerate_mu.png", "ratio_range":[0.0, 2.0], "draw_points":True, "yaxis_label":"Fake Rate", "no_ratio":True, "xaxis_label":"#it{p}_{T,cone-corr} [GeV]", "ymax_scale":0.65})

def draw_fakerate_1d_el(nfscheme=nfs_from_CR_el, nfsinclqcd=False):
    h_fakerate_1d_el_data = fakerate_1d_el_data_hist(nfscheme=nfscheme, nfsinclqcd=nfsinclqcd).Clone("h_fakerate_1d_el_data")
    #h_fakerate_1d_el_data.Print("all")
    h_fakerate_1d_el_data_syst13 = fakerate_1d_el_data_hist("syst13_", nfscheme=nfscheme, nfsinclqcd=nfsinclqcd).Clone("h_fakerate_1d_el_data")
    #h_fakerate_1d_el_data_syst13.Print("all")
    h_fakerate_1d_el_data_syst14 = fakerate_1d_el_data_hist("syst14_", nfscheme=nfscheme, nfsinclqcd=nfsinclqcd).Clone("h_fakerate_1d_el_data")
    #h_fakerate_1d_el_data_syst14.Print("all")
    h_fakerate_1d_el_data_fullerror = get_full_error(h_fakerate_1d_el_data, h_fakerate_1d_el_data_syst13, h_fakerate_1d_el_data_syst14)
    #h_fakerate_1d_el_data_fullerror.Print("all")
    h_fakerate_1d_el_qcd  = fakerate_1d_el_qcd_hist().Clone("QCD")
    #h_fakerate_1d_el_qcd .Print("all")
    p.plot_hist(sigs=[], bgs=[h_fakerate_1d_el_qcd], data=h_fakerate_1d_el_data_fullerror, colors=[2], options={"output_name":"frplots/fakerate_el.png", "ratio_range":[0.0, 2.0], "draw_points":True, "yaxis_label":"Fake Rate", "no_ratio":True, "xaxis_label":"#it{p}_{T,cone-corr} [GeV]", "ymax_scale":1.0})

def draw_fakerate_2d_mu(nfscheme=nfs_from_CR_mu, nfsinclqcd=False):
    h_fakerate_2d_mu_data = fakerate_2d_mu_data_hist(nfscheme=nfscheme, nfsinclqcd=nfsinclqcd).Clone("h_fakerate_2d_mu_data")
    #h_fakerate_2d_mu_data.Print("all")
    h_fakerate_2d_mu_data_syst13 = fakerate_2d_mu_data_hist("syst13_", nfscheme=nfscheme, nfsinclqcd=nfsinclqcd).Clone("h_fakerate_2d_mu_data")
    #h_fakerate_2d_mu_data_syst13.Print("all")
    h_fakerate_2d_mu_data_syst14 = fakerate_2d_mu_data_hist("syst14_", nfscheme=nfscheme, nfsinclqcd=nfsinclqcd).Clone("h_fakerate_2d_mu_data")
    #h_fakerate_2d_mu_data_syst14.Print("all")
    h_fakerate_2d_mu_data_fullerror = get_full_error_2d(h_fakerate_2d_mu_data, h_fakerate_2d_mu_data_syst13, h_fakerate_2d_mu_data_syst14)
    #h_fakerate_2d_mu_data_fullerror.Print("all")
    h_fakerate_2d_mu_qcd  = fakerate_2d_mu_qcd_hist().Clone("QCD")
    #h_fakerate_2d_mu_qcd .Print("all")
    max_data = h_fakerate_2d_mu_data_fullerror.GetMaximum()
    min_data = h_fakerate_2d_mu_data_fullerror.GetMinimum()
    max_qcd = h_fakerate_2d_mu_qcd.GetMaximum()
    min_qcd = h_fakerate_2d_mu_qcd.GetMinimum()
    h_fakerate_2d_mu_data_fullerror.Print("all")
    ply.plot_hist_2d( h_fakerate_2d_mu_data_fullerror, options = { "output_name": "frplots/fakerate_2d_mu_data.png", "zaxis_range": [min_data/1.5, 1.5*max_data], "zaxis_log": False, "bin_text_smart": False, "us_flag": False, "output_ic": False, "zaxis_noexponents": True, "draw_option_2d": "textecolz", "bin_text_format": ".3f", "xaxis_log": True, "bin_text_size": 1.0, "palette_name": "radiation", "xaxis_label":"#it{p}_{T,cone-corr} [GeV]", "yaxis_label":"|#eta|", "xaxis_title_offset":1.4, "yaxis_title_offset":1.4 })
    ply.plot_hist_2d( h_fakerate_2d_mu_qcd           , options = { "output_name": "frplots/fakerate_2d_mu_qcd.png" , "zaxis_range": [min_qcd /1.5, 1.5*max_qcd ], "zaxis_log": False, "bin_text_smart": False, "us_flag": False, "output_ic": False, "zaxis_noexponents": True, "draw_option_2d": "textecolz", "bin_text_format": ".3f", "xaxis_log": True, "bin_text_size": 1.0, "palette_name": "radiation", "xaxis_label":"#it{p}_{T,cone-corr} [GeV]", "yaxis_label":"|#eta|", "xaxis_title_offset":1.4, "yaxis_title_offset":1.4 })
    ply.plot_hist_2d( h_fakerate_2d_mu_data_syst13   , options = { "output_name": "frplots/syst13_fakerate_2d_mu_data.png", "zaxis_range": [min_data/1.5, 1.5*max_data], "zaxis_log": False, "bin_text_smart": False, "us_flag": False, "output_ic": False, "zaxis_noexponents": True, "draw_option_2d": "textecolz", "bin_text_format": ".3f", "xaxis_log": True, "bin_text_size": 1.0, "palette_name": "radiation", "xaxis_label":"#it{p}_{T,cone-corr} [GeV]", "yaxis_label":"|#eta|", "xaxis_title_offset":1.4, "yaxis_title_offset":1.4 })
    ply.plot_hist_2d( h_fakerate_2d_mu_data_syst14   , options = { "output_name": "frplots/syst14_fakerate_2d_mu_data.png", "zaxis_range": [min_data/1.5, 1.5*max_data], "zaxis_log": False, "bin_text_smart": False, "us_flag": False, "output_ic": False, "zaxis_noexponents": True, "draw_option_2d": "textecolz", "bin_text_format": ".3f", "xaxis_log": True, "bin_text_size": 1.0, "palette_name": "radiation", "xaxis_label":"#it{p}_{T,cone-corr} [GeV]", "yaxis_label":"|#eta|", "xaxis_title_offset":1.4, "yaxis_title_offset":1.4 })
    return h_fakerate_2d_mu_data_fullerror, h_fakerate_2d_mu_qcd

def draw_fakerate_2d_el(nfscheme=nfs_from_CR_el, nfsinclqcd=False):
    h_fakerate_2d_el_data = fakerate_2d_el_data_hist(nfscheme=nfscheme, nfsinclqcd=nfsinclqcd).Clone("h_fakerate_2d_el_data")
    #h_fakerate_2d_el_data.Print("all")
    h_fakerate_2d_el_data_syst13 = fakerate_2d_el_data_hist("syst13_", nfscheme=nfscheme, nfsinclqcd=nfsinclqcd).Clone("h_fakerate_2d_el_data")
    #h_fakerate_2d_el_data_syst13.Print("all")
    h_fakerate_2d_el_data_syst14 = fakerate_2d_el_data_hist("syst14_", nfscheme=nfscheme, nfsinclqcd=nfsinclqcd).Clone("h_fakerate_2d_el_data")
    #h_fakerate_2d_el_data_syst14.Print("all")
    h_fakerate_2d_el_data_fullerror = get_full_error_2d(h_fakerate_2d_el_data, h_fakerate_2d_el_data_syst13, h_fakerate_2d_el_data_syst14)
    #h_fakerate_2d_el_data_fullerror.Print("all")
    h_fakerate_2d_el_qcd  = fakerate_2d_el_qcd_hist().Clone("QCD")
    #h_fakerate_2d_el_qcd .Print("all")
    max_data = h_fakerate_2d_el_data_fullerror.GetMaximum()
    min_data = h_fakerate_2d_el_data_fullerror.GetMinimum()
    max_qcd = h_fakerate_2d_el_qcd.GetMaximum()
    min_qcd = h_fakerate_2d_el_qcd.GetMinimum()
    ply.plot_hist_2d( h_fakerate_2d_el_data_fullerror, options = { "output_name": "frplots/fakerate_2d_el_data.png", "zaxis_range": [min_data/1.5, 1.5*max_data], "zaxis_log": False, "bin_text_smart": False, "us_flag": False, "output_ic": False, "zaxis_noexponents": True, "draw_option_2d": "textecolz", "bin_text_format": ".3f", "xaxis_log": True, "bin_text_size": 1.0, "palette_name": "radiation", "xaxis_label":"#it{p}_{T,cone-corr} [GeV]", "yaxis_label":"|#eta|", "xaxis_title_offset":1.4, "yaxis_title_offset":1.4 })
    ply.plot_hist_2d( h_fakerate_2d_el_qcd           , options = { "output_name": "frplots/fakerate_2d_el_qcd.png" , "zaxis_range": [min_qcd /1.5, 1.5*max_qcd ], "zaxis_log": False, "bin_text_smart": False, "us_flag": False, "output_ic": False, "zaxis_noexponents": True, "draw_option_2d": "textecolz", "bin_text_format": ".3f", "xaxis_log": True, "bin_text_size": 1.0, "palette_name": "radiation", "xaxis_label":"#it{p}_{T,cone-corr} [GeV]", "yaxis_label":"|#eta|", "xaxis_title_offset":1.4, "yaxis_title_offset":1.4 })
    return h_fakerate_2d_el_data_fullerror, h_fakerate_2d_el_qcd

###################################################################################################
#
#
# Finally, the scripts to actually perform the jobs.
#
#
###################################################################################################
def draw_ewkcr_1d_mu(): plot("mt", "CR", "tightmu", {"output_name": "frplots/plot_cr_tightmu.png", "no_ratio":True, "xaxis_label":"#it{m}_{T} [GeV]"}, True)
def draw_ewkcr_1d_el(): plot("mt", "CR", "tightel", {"output_name": "frplots/plot_cr_tightel.png", "no_ratio":True, "xaxis_label":"#it{m}_{T} [GeV]"}, True)
def draw_ewkcr2_1d_mu(): plot("mt", "CR2", "tightmu", {"output_name": "frplots/plot_cr2_tightmu.png", "no_ratio":True, "xaxis_label":"#it{m}_{T} [GeV]"}, True)
def draw_ewkcr2_1d_el(): plot("mt", "CR2", "tightel", {"output_name": "frplots/plot_cr2_tightel.png", "no_ratio":True, "xaxis_label":"#it{m}_{T} [GeV]"}, True)
def draw_ewkcr3_1d_mu(): plot("mt", "CR3", "tightmu", {"output_name": "frplots/plot_cr3_tightmu.png", "no_ratio":True, "xaxis_label":"#it{m}_{T} [GeV]"}, True)
def draw_ewkcr3_1d_el(): plot("mt", "CR3", "tightel", {"output_name": "frplots/plot_cr3_tightel.png", "no_ratio":True, "xaxis_label":"#it{m}_{T} [GeV]"}, True)
def draw_ewkcrnomt_1d_el(): plot("mt", "CR_noMT", "tightel", {"output_name": "frplots/plot_crnomt_tightel.png", "no_ratio":True, "xaxis_label":"#it{m}_{T} [GeV]", "yaxis_range":[0.,2.5e6]}, True)

###################################################################################################
#
#
# Main function
#
#
###################################################################################################

def main():
    of = r.TFile("frplots/fakerate.root", "recreate")
    draw_fakerate_1d_mu()
    draw_fakerate_1d_el(nfsinclqcd=True)
    d, q = draw_fakerate_2d_mu()
    of.cd()
    d.Clone("fakerate_mu_data").Write()
    q.Clone("fakerate_mu_qcd").Write()
    d, q = draw_fakerate_2d_el()
    of.cd()
    d.Clone("fakerate_el_data").Write()
    q.Clone("fakerate_el_qcd").Write()
    draw_ewkcr_1d_mu()
    draw_ewkcr_1d_el()
    draw_ewkcrnomt_1d_el()
    draw_ewkcr2_1d_mu()
    draw_ewkcr2_1d_el()
    draw_ewkcr3_1d_mu()
    draw_ewkcr3_1d_el()
    print nfs_from_CR_mu("",False)
    print nfs_from_CR2_mu("",False)
    print nfs_from_CR_el("",False)
    print nfs_from_CR2_el("",False)

def main2():
    of = r.TFile("frplots/fakerate.root", "recreate")
    draw_ewkcr2_1d_mu()
    print nfs_from_CR2_mu("",inclqcd=True)
    draw_fakerate_1d_mu(nfscheme=nfs_from_CR2_mu, nfsinclqcd=True)
    d, q = draw_fakerate_2d_mu(nfscheme=nfs_from_CR2_mu, nfsinclqcd=True)
    d.Clone("fakerate_mu_data").Write()
    q.Clone("fakerate_mu_qcd").Write()
    draw_ewkcr2_1d_el()
    print nfs_from_CR2_el("",inclqcd=True)
    draw_fakerate_1d_el(nfscheme=nfs_from_CR2_el, nfsinclqcd=True)
    d, q = draw_fakerate_2d_el(nfscheme=nfs_from_CR2_el, nfsinclqcd=True)
    d.Clone("fakerate_el_data").Write()
    q.Clone("fakerate_el_qcd").Write()


if __name__ == "__main__":
    main2()
    #main()
    #print nfs_from_CR_mu("", False)
    #draw_fakerate_1d_mu(nfs_from_CR_mu, nfsinclqcd=False)

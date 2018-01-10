#!/bin/env python

import sys
sys.path.append("/home/users/phchang/public_html/phys/rooutil")
import plotmaker
import ROOT as r

r.gStyle.SetPalette( r.kInvertedDarkBodyRadiator )

tfile_data   = r.TFile("outputs/fakerate_data.root")
tfile_wj     = r.TFile("outputs/fakerate_wj.root")
tfile_dy     = r.TFile("outputs/fakerate_dy.root")
tfile_ttbar  = r.TFile("outputs/fakerate_ttbar.root")
tfile_vv     = r.TFile("outputs/fakerate_vv.root")
tfile_qcd_mu = r.TFile("outputs/fakerate_qcd_mu.root")
tfile_qcd_el = r.TFile("outputs/fakerate_qcd_el.root")
#try:
#    tfile_test   = r.TFile("fakerate_DY_madgraph_0.root")
#except:
#    pass

from array import *

global prefix
prefix = ""
norebin = False
addvv = True
addqcd = False

def quadsum( numlist ):
    sys.path.append("/home/users/phchang/syncfiles/pyfiles/")
    from errors import E
    errors = [ E(0.0, x) for x in numlist ]
    return sum( errors ).err

###################################################################################################
# Rebinning histograms when needed.
# Currently only the "ptvarbin" variable and "met" is being rebinned to do some "studies".
#
def rebin(hist, name):
    return hist
    if norebin:
        return hist
    if name.find("ptvarbin") != -1:
#        bins = array('d', [10., 15., 20., 25., 30., 35., 50., 170.])
#        tmphist = hist.Rebin(7, hist.GetName(), bins)
#        bins = array('d', [10., 15., 20., 25., 30., 35., 170.])
#        tmphist = hist.Rebin(6, hist.GetName(), bins)
        bins = array('d', [10., 15., 20., 25., 30., 35., 50., 170.])
        tmphist = hist.Rebin(7, hist.GetName(), bins)
        tmphist.SetDirectory( 0 )
        return tmphist
    if name.find("histo_nvtxrewgt_met_") != -1:
        bins = array('d', [0., 20., 30., 200.])
        tmphist = hist.Rebin(3, hist.GetName(), bins)
        tmphist.SetDirectory( 0 )
        return tmphist
    return hist

###################################################################################################
# Rebinning histograms when needed.
# Currently only the "ptvarbin" variable and "met" is being rebinned to do some "studies".
#
def rebin2d(hist, name):
    return hist
    ptbins  = array('d', [10., 15., 20., 25., 30., 35., 50., 170.])
    etabins_mu = array('d', [0., 1.2, 2.1, 2.4])
    etabins_el = array('d', [0., 0.8, 1.479, 2.5])
    if norebin:
        return hist
    tmphist = None
    if name.find("ptvarbin") != -1 and name.find("mu") != -1:
        tmphist = r.TH2F( name, name, 7, ptbins, 3, etabins_mu)
        tmphist.SetDirectory( 0 )
    elif name.find("ptvarbin") != -1 and name.find("el") != -1:
        tmphist = r.TH2F( name, name, 7, ptbins, 3, etabins_el)
        tmphist.SetDirectory( 0 )
    for ix in xrange( 1, 7 ):
        for iy in xrange( 1, 4 ):
            tmphist.SetBinContent( ix, iy, hist.GetBinContent( ix, iy ) )
            tmphist.SetBinError( ix, iy, hist.GetBinError( ix, iy ) )
    tmphist.SetBinContent( 7, 1, hist.GetBinContent( 7, 1 ) + hist.GetBinContent( 8, 1 ) + hist.GetBinContent( 9, 1 ) + hist.GetBinContent( 10, 1 ) )
    tmphist.SetBinContent( 7, 2, hist.GetBinContent( 7, 2 ) + hist.GetBinContent( 8, 2 ) + hist.GetBinContent( 9, 2 ) + hist.GetBinContent( 10, 2 ) )
    tmphist.SetBinContent( 7, 3, hist.GetBinContent( 7, 3 ) + hist.GetBinContent( 8, 3 ) + hist.GetBinContent( 9, 3 ) + hist.GetBinContent( 10, 3 ) )
    tmphist.SetBinError( 7, 1, quadsum( [ hist.GetBinError( 7, 1 ) , hist.GetBinError( 8, 1 ) , hist.GetBinError( 9, 1 ) , hist.GetBinError( 10, 1 ) ] ) )
    tmphist.SetBinError( 7, 2, quadsum( [ hist.GetBinError( 7, 2 ) , hist.GetBinError( 8, 2 ) , hist.GetBinError( 9, 2 ) , hist.GetBinError( 10, 2 ) ] ) )
    tmphist.SetBinError( 7, 3, quadsum( [ hist.GetBinError( 7, 3 ) , hist.GetBinError( 8, 3 ) , hist.GetBinError( 9, 3 ) , hist.GetBinError( 10, 3 ) ] ) )
    return tmphist

###################################################################################################
# Rebinning histograms when needed.
# Currently only the "ptvarbin" variable and "met" is being rebinned to do some "studies".
#
def rebin2d_2dmu(hist, name):
    ix = 7
    val = hist.GetBinContent( ix, 3 ) + hist.GetBinContent( ix, 2 )
    err = quadsum( [ hist.GetBinError( ix, 3 ), hist.GetBinError( ix, 2 ) ] )
    hist.SetBinContent( ix, 3, val )
    hist.SetBinContent( ix, 2, val )
    hist.SetBinError( ix, 3, err )
    hist.SetBinError( ix, 2, err )
    return hist

####################################################################################################
## Rebinning histograms when needed.
## Currently only the "ptvarbin" variable and "met" is being rebinned to do some "studies".
##
#def printFakeRateFunction( ratio, name ):
#    print "double getFR_%s( float conecorrpt, 

###################################################################################################
# Plot 1D distribution (default : data, W, DY. optional : add QCD, or ttbar)
#
def draw(histname, varname, extraopt="", scale=1.00, printcontent=False):

    # Get histograms
    histo__data = tfile_data .Get(histname).Clone("Data")
    histo__wj   = tfile_wj   .Get(histname).Clone("WJets")
    histo__dy   = tfile_dy   .Get(histname).Clone("DY")
    histo__ttbar= tfile_ttbar.Get(histname).Clone("TTbar")
    if addvv: histo__vv   = tfile_vv.Get(histname).Clone("VV")
    try:
        if addqcd: histo__qcd   = tfile_qcd_mu.Get(histname).Clone("QCD")
    except:
        if addqcd: histo__qcd   = tfile_qcd_el.Get(histname).Clone("QCD")
    histo__data .SetDirectory( 0 )
    histo__wj   .SetDirectory( 0 )
    histo__dy   .SetDirectory( 0 )
    histo__ttbar.SetDirectory( 0 )
    if addvv: histo__vv   .SetDirectory( 0 )
    if addqcd: histo__qcd   .SetDirectory( 0 )
    histo__wj   .Scale( scale )
    histo__dy   .Scale( scale )
    histo__ttbar.Scale( scale )
    if addvv: histo__vv   .Scale( scale )
    if addqcd: histo__qcd   .Scale( scale )
    histo__data .SetLineColor( 1 )
    histo__wj   .SetLineColor( 7006 )
    histo__dy   .SetLineColor( 7002 )
    histo__ttbar.SetLineColor( 7004 )
    if addvv: histo__vv   .SetLineColor( 7005 )
    if addqcd: histo__qcd   .SetLineColor( 5 )
    histo__wj   .SetFillColor( 7006 )
    histo__dy   .SetFillColor( 7002 )
    histo__ttbar.SetFillColor( 7004 )
    if addvv: histo__vv   .SetFillColor( 7005 )
    if addqcd: histo__qcd   .SetFillColor( 5 )
    v_bkg_hists = r.vector("TH1*")()
    v_bkg_hists.push_back(histo__wj)
    v_bkg_hists.push_back(histo__dy)
    v_bkg_hists.push_back(histo__ttbar)
    if addvv: v_bkg_hists.push_back(histo__vv)
    if addqcd: v_bkg_hists.push_back(histo__qcd)

    #if printcontent:
    #    totalbkg = r.getTotalBkgHists(v_bkg_hists)
    #    totalbkg.Rebin(4)
    #    totalbkg.Print("all")
    #    #for hbg in v_bkg_hists:
    #    #    hbg.Print("all")
    #    histo__data.Rebin(4)
    #    histo__data.Print("all")
    #    b = totalbkg.GetBinContent(3) + totalbkg.GetBinContent(4) + totalbkg.GetBinContent(5)
    #    d = histo__data.GetBinContent(3) + histo__data.GetBinContent(4) + histo__data.GetBinContent(5)
    #    print b, d

    # Run the plot and return ratio plot
    # The plotmaker function returns the data/MC ratio TH1 created during the plotting process.
    return r.plotmaker( """
                  --yTitle N leptons
                  --xTitle %s
                  --plotOutputName plots/%s
                  --ratio_Maximum 1.4
                  --ratio_Minimum 0.6
                  --showOverflow
                  --MaximumMultiplier 1.5
                  --autoStack
                  %s
                  --ratio_DrawOpt ep
                  """%(varname, histname, extraopt) ,
                  histo__data, v_bkg_hists )
                  #--ratioPaneAtBottom

###################################################################################################
# Plot 1D distribution (default : data, W, DY. optional : add QCD, or ttbar)
#
def draw_fakerate_mu(extraopt="", scale=1.0):
    global prefix

    histname = prefix + "evt_lvl_histo_conecorrptvarbin_meas_mu"

    # Get histograms
    histo__tight = rebin( tfile_data .Get(histname).Clone("Tight"), histname )
    histo__wj    = rebin( tfile_wj   .Get(histname).Clone("WJets"), histname )
    histo__dy    = rebin( tfile_dy   .Get(histname).Clone("DY")   , histname )
    histo__ttbar = rebin( tfile_ttbar.Get(histname).Clone("TTbar"), histname )
    if addvv: histo__vv    = rebin( tfile_vv   .Get(histname).Clone("VV")   , histname )
    histo__tight.SetDirectory( 0 )
    histo__wj   .SetDirectory( 0 )
    histo__dy   .SetDirectory( 0 )
    histo__ttbar.SetDirectory( 0 )
    if addvv: histo__vv   .SetDirectory( 0 )
    histo__tight.Add( histo__wj, -scale )
    histo__tight.Add( histo__dy, -scale )
    histo__tight.Add( histo__ttbar, -scale )
    if addvv: histo__tight.Add( histo__vv, -scale )
    histo__tight.SetLineColor( 1 )
    histo__tight.SetFillColor( r.kGray + 1 )

    histname = prefix + "evt_lvl_histo_conecorrptvarbin_loose_meas_mu"

    # Get histograms
    histo__loose = rebin( tfile_data .Get(histname).Clone("Loose"), histname )
    histo__wj    = rebin( tfile_wj   .Get(histname).Clone("WJets"), histname )
    histo__dy    = rebin( tfile_dy   .Get(histname).Clone("DY")   , histname )
    histo__ttbar = rebin( tfile_ttbar.Get(histname).Clone("TTbar"), histname )
    if addvv: histo__vv    = rebin( tfile_vv   .Get(histname).Clone("VV")   , histname )
    histo__loose.SetDirectory( 0 )
    histo__wj   .SetDirectory( 0 )
    histo__dy   .SetDirectory( 0 )
    histo__ttbar.SetDirectory( 0 )
    if addvv: histo__vv   .SetDirectory( 0 )
    histo__loose.Add( histo__wj, -scale )
    histo__loose.Add( histo__dy, -scale )
    histo__loose.Add( histo__ttbar, -scale )
    if addvv: histo__loose.Add( histo__vv, -scale )
    histo__loose.SetLineColor( 1 )

    v_bkg_hists = r.vector("TH1*")()
    v_bkg_hists.push_back(histo__tight)

    #print "HERE2"
    histo__tight.Print("all")
    histo__loose.Print("all")

    # Run the plot and return ratio plot
    # The plotmaker function returns the data/MC ratio TH1 created during the plotting process.
    return r.plotmaker( """
                  --yTitle N leptons
                  --xTitle p_{T,#mu} [GeV]
                  --plotOutputName plots/%sfakerate_mu_data
                  --ratio_DrawOpt ep
                  --ratio_Maximum 0.27
                  --ratio_Minimum 0.0
                  --reverseRatio
                  --saveMainPad
                  --showOverflow
                  --divideByBinWidth
                  %s
                  --MaximumMultiplier 6
                  --autoStack
                  """%(prefix, extraopt) ,
                  histo__loose, v_bkg_hists )

###################################################################################################
# Plot 1D distribution (default : data, W, DY. optional : add QCD, or ttbar)
#
def draw_nvtxrewgt_fakerate_el(extraopt="", scale=1.0):
    global prefix

    histname = prefix + "evt_lvl_nvtxrewgt_histo_conecorrptvarbin_meas_el"

    # Get histograms
    histo__tight = rebin( tfile_data .Get(histname).Clone("Tight"), histname )
    histo__wj    = rebin( tfile_wj   .Get(histname).Clone("WJets"), histname )
    histo__dy    = rebin( tfile_dy   .Get(histname).Clone("DY")   , histname )
    histo__ttbar = rebin( tfile_ttbar.Get(histname).Clone("TTbar"), histname )
    if addvv: histo__vv    = rebin( tfile_vv   .Get(histname).Clone("VV")   , histname )
    histo__tight.SetDirectory( 0 )
    histo__wj   .SetDirectory( 0 )
    histo__dy   .SetDirectory( 0 )
    histo__ttbar.SetDirectory( 0 )
    if addvv: histo__vv   .SetDirectory( 0 )
    histo__tight.Add( histo__wj, -scale )
    histo__tight.Add( histo__dy, -scale )
    histo__tight.Add( histo__ttbar, -scale )
    if addvv: histo__tight.Add( histo__vv   , -scale )
    histo__tight.SetLineColor( 1 )
    histo__tight.SetFillColor( r.kGray + 1 )

    histname = prefix + "evt_lvl_nvtxrewgt_histo_conecorrptvarbin_loose_meas_el"

    # Get histograms
    histo__loose = rebin( tfile_data .Get(histname).Clone("Loose"), histname )
    histo__wj    = rebin( tfile_wj   .Get(histname).Clone("WJets"), histname )
    histo__dy    = rebin( tfile_dy   .Get(histname).Clone("DY")   , histname )
    histo__ttbar = rebin( tfile_ttbar.Get(histname).Clone("TTbar"), histname )
    if addvv: histo__vv    = rebin( tfile_vv   .Get(histname).Clone("VV")   , histname )
    histo__loose.SetDirectory( 0 )
    histo__wj   .SetDirectory( 0 )
    histo__dy   .SetDirectory( 0 )
    histo__ttbar.SetDirectory( 0 )
    if addvv: histo__vv   .SetDirectory( 0 )
    histo__loose.Add( histo__wj, -scale )
    histo__loose.Add( histo__dy, -scale )
    histo__loose.Add( histo__ttbar, -scale )
    if addvv: histo__loose.Add( histo__vv   , -scale )
    histo__loose.SetLineColor( 1 )

    v_bkg_hists = r.vector("TH1*")()
    v_bkg_hists.push_back(histo__tight)

    histo__loose.Print("all")
    histo__tight.Print("all")

    # Run the plot and return ratio plot
    # The plotmaker function returns the data/MC ratio TH1 created during the plotting process.
    return r.plotmaker( """
                  --yTitle N leptons
                  --xTitle p^{corr}_{T,e} [GeV]
                  --plotOutputName plots/%sfakerate_nvtxrewgt_el_data
                  --ratio_DrawOpt ep
                  --ratio_Maximum 0.37
                  --ratio_yTitle Fake Rate
                  --ratio_Minimum 0.0
                  --saveMainPad
                  --reverseRatio
                  --showOverflow
                  --divideByBinWidth
                  %s
                  --MaximumMultiplier 6
                  --autoStack
                  """%(prefix, extraopt) ,
                  histo__loose, v_bkg_hists )

###################################################################################################
# Plot 1D distribution (default : data, W, DY. optional : add QCD, or ttbar)
#
def draw_nvtxrewgt_fakerate_mu(extraopt="", scale=1.0, tightname="evt_lvl_nvtxrewgt_histo_conecorrptvarbin_meas_mu", loosename="evt_lvl_nvtxrewgt_histo_conecorrptvarbin_loose_meas_mu", outputname="fakerate_nvtxrewgt_mu_data"):
    global prefix
    print "hereHEREHREHREH", prefix

    histname = prefix + tightname

    # Get histograms
    histo__tight = rebin( tfile_data .Get(histname).Clone("Tight"), histname )
    histo__wj    = rebin( tfile_wj   .Get(histname).Clone("WJets"), histname )
    histo__dy    = rebin( tfile_dy   .Get(histname).Clone("DY")   , histname )
    histo__ttbar = rebin( tfile_ttbar.Get(histname).Clone("TTbar"), histname )
    if addvv: histo__vv    = rebin( tfile_vv   .Get(histname).Clone("VV")   , histname )
    histo__tight.SetDirectory( 0 )
    histo__wj   .SetDirectory( 0 )
    histo__dy   .SetDirectory( 0 )
    histo__ttbar.SetDirectory( 0 )
    if addvv: histo__vv   .SetDirectory( 0 )
    histo__tight.Add( histo__wj, -scale )
    histo__tight.Add( histo__dy, -scale )
    histo__tight.Add( histo__ttbar, -scale )
    if addvv: histo__tight.Add( histo__vv   , -scale )
    histo__tight.SetLineColor( 1 )
    histo__tight.SetFillColor( r.kGray + 1 )

    histname = prefix + loosename

    # Get histograms
    histo__loose = rebin( tfile_data .Get(histname).Clone("Loose"), histname )
    histo__wj    = rebin( tfile_wj   .Get(histname).Clone("WJets"), histname )
    histo__dy    = rebin( tfile_dy   .Get(histname).Clone("DY")   , histname )
    histo__ttbar = rebin( tfile_ttbar.Get(histname).Clone("TTbar"), histname )
    if addvv: histo__vv    = rebin( tfile_vv   .Get(histname).Clone("VV")   , histname )
    histo__loose.SetDirectory( 0 )
    histo__wj   .SetDirectory( 0 )
    histo__dy   .SetDirectory( 0 )
    histo__ttbar.SetDirectory( 0 )
    if addvv: histo__vv   .SetDirectory( 0 )
    histo__loose.Add( histo__wj, -scale )
    histo__loose.Add( histo__dy, -scale )
    histo__loose.Add( histo__ttbar, -scale )
    if addvv: histo__loose.Add( histo__vv   , -scale )
    histo__loose.SetLineColor( 1 )

    v_bkg_hists = r.vector("TH1*")()
    v_bkg_hists.push_back(histo__tight)

    histo__tight.Print("All")
    histo__loose.Print("All")

    # Run the plot and return ratio plot
    # The plotmaker function returns the data/MC ratio TH1 created during the plotting process.
    return r.plotmaker( """
                  --yTitle N leptons
                  --xTitle p^{corr}_{T,#mu} [GeV]
                  --plotOutputName plots/%s%s
                  --ratio_DrawOpt ep
                  --ratio_Maximum 0.27
                  --ratio_Minimum 0.0
                  --ratio_yTitle Fake Rate
                  --saveMainPad
                  --reverseRatio
                  --showOverflow
                  --divideByBinWidth
                  %s
                  --MaximumMultiplier 6
                  --autoStack
                  """%(prefix, outputname, extraopt) ,
                  histo__loose, v_bkg_hists )

###################################################################################################
# Plot 1D distribution (default : data, W, DY. optional : add QCD, or ttbar)
#
def draw_fakerate_el(extraopt="", scale=1.0):
    global prefix

    histname = prefix + "evt_lvl_histo_conecorrptvarbin_meas_el"

    # Get histograms
    histo__tight = rebin( tfile_data .Get(histname).Clone("Tight"), histname )
    histo__wj    = rebin( tfile_wj   .Get(histname).Clone("WJets"), histname )
    histo__dy    = rebin( tfile_dy   .Get(histname).Clone("DY")   , histname )
    histo__ttbar = rebin( tfile_ttbar.Get(histname).Clone("TTbar"), histname )
    if addvv: histo__vv    = rebin( tfile_vv   .Get(histname).Clone("VV")   , histname )
    histo__tight.SetDirectory( 0 )
    histo__wj   .SetDirectory( 0 )
    histo__dy   .SetDirectory( 0 )
    histo__ttbar.SetDirectory( 0 )
    if addvv: histo__vv   .SetDirectory( 0 )
    histo__tight.Add( histo__wj, -scale )
    histo__tight.Add( histo__dy, -scale )
    histo__tight.Add( histo__ttbar, -scale )
    if addvv: histo__tight.Add( histo__vv   , -scale )
    histo__tight.SetLineColor( 1 )
    histo__tight.SetFillColor( r.kGray + 1 )

    histname = prefix + "evt_lvl_histo_conecorrptvarbin_loose_meas_el"

    # Get histograms
    histo__loose = rebin( tfile_data .Get(histname).Clone("Loose"), histname )
    histo__wj    = rebin( tfile_wj   .Get(histname).Clone("WJets"), histname )
    histo__dy    = rebin( tfile_dy   .Get(histname).Clone("DY")   , histname )
    histo__ttbar = rebin( tfile_ttbar.Get(histname).Clone("TTbar"), histname )
    if addvv: histo__vv    = rebin( tfile_vv   .Get(histname).Clone("VV")   , histname )
    histo__loose.SetDirectory( 0 )
    histo__wj   .SetDirectory( 0 )
    histo__dy   .SetDirectory( 0 )
    histo__ttbar.SetDirectory( 0 )
    if addvv: histo__vv   .SetDirectory( 0 )
    histo__loose.Add( histo__wj, -scale )
    histo__loose.Add( histo__dy, -scale )
    histo__loose.Add( histo__ttbar, -scale )
    if addvv: histo__loose.Add( histo__vv   , -scale )
    histo__loose.SetLineColor( 1 )

    v_bkg_hists = r.vector("TH1*")()
    v_bkg_hists.push_back(histo__tight)

    # Run the plot and return ratio plot
    # The plotmaker function returns the data/MC ratio TH1 created during the plotting process.
    return r.plotmaker( """
                  --yTitle N leptons
                  --xTitle p_{T,e} [GeV]
                  --plotOutputName plots/%sfakerate_el_data
                  --ratio_DrawOpt ep
                  --ratio_Maximum 0.27
                  --ratio_Minimum 0.0
                  --reverseRatio
                  --showOverflow
                  --ratioPaneAtBottom
                  --divideByBinWidth
                  %s
                  --MaximumMultiplier 6
                  --autoStack
                  """%(prefix, extraopt) ,
                  histo__loose, v_bkg_hists )

###################################################################################################
# Plot 1D distribution (default : data, W, DY. optional : add QCD, or ttbar)
#
def draw_fakerate_qcd_mu(extraopt="", tightname="evt_lvl_histo_conecorrptvarbin_meas_mu", loosename="evt_lvl_histo_conecorrptvarbin_loose_meas_mu", outputname="fakerate_mu_qcd"):
    global prefix

    histname = tightname

    # Get histograms
    histo__tight = rebin( tfile_qcd_mu .Get(histname).Clone("Tight"), histname )
    histo__tight.SetDirectory( 0 )
    histo__tight.SetLineColor( 1 )
    histo__tight.SetFillColor( r.kYellow )

    histname = loosename

    # Get histograms
    histo__loose = rebin( tfile_qcd_mu .Get(histname).Clone("Loose"), histname )
    histo__loose.SetDirectory( 0 )
    histo__loose.SetLineColor( 2 )
    histo__loose.SetMarkerColor( 2 )

    v_bkg_hists = r.vector("TH1*")()
    v_bkg_hists.push_back(histo__tight)

    print "HERE1"
    histo__tight.Print("all")
    histo__loose.Print("all")

    # Run the plot and return ratio plot
    # The plotmaker function returns the data/MC ratio TH1 created during the plotting process.
    return r.plotmaker( """
                  --yTitle N leptons
                  --xTitle p^{corr}_{T,#mu} [GeV]
                  --plotOutputName plots/%s
                  --ratio_DrawOpt ep
                  --ratio_Maximum 0.27
                  --ratio_Minimum 0.0
                  --reverseRatio
                  --showOverflow
                  --divideByBinWidth
                  %s
                  --MaximumMultiplier 6
                  --autoStack
                  """%(outputname, extraopt) ,
                  histo__loose, v_bkg_hists )

###################################################################################################
# Plot 1D distribution (default : data, W, DY. optional : add QCD, or ttbar)
#
def draw_fakerate_qcd_el(extraopt=""):
    global prefix

    histname = "evt_lvl_histo_conecorrptvarbin_meas_el"

    # Get histograms
    histo__tight = rebin( tfile_qcd_el .Get(histname).Clone("Tight"), histname )
    histo__tight.SetDirectory( 0 )
    histo__tight.SetLineColor( 1 )
    histo__tight.SetFillColor( r.kYellow )

    histname = "evt_lvl_histo_conecorrptvarbin_loose_meas_el"

    # Get histograms
    histo__loose = rebin( tfile_qcd_el .Get(histname).Clone("Loose"), histname )
    histo__loose.SetDirectory( 0 )
    histo__loose.SetLineColor( 2 )
    histo__loose.SetMarkerColor( 2 )

    v_bkg_hists = r.vector("TH1*")()
    v_bkg_hists.push_back(histo__tight)

    # Run the plot and return ratio plot
    # The plotmaker function returns the data/MC ratio TH1 created during the plotting process.
    return r.plotmaker( """
                  --yTitle N leptons
                  --xTitle p^{corr}_{T,e} [GeV]
                  --ratio_DrawOpt ep
                  --ratio_Maximum 0.27
                  --ratio_Minimum 0.0
                  --reverseRatio
                  --showOverflow
                  --MaximumMultiplier 6
                  --autoStack
                  %s
                  --plotOutputName plots/fakerate_el_qcd
                  """%(extraopt) ,
                  histo__loose, v_bkg_hists )

###################################################################################################
# Plot 2D distribution (default : data, W, DY. optional : add QCD, or ttbar)
#
def draw_nvtxrewgt_fakerate_mu_2d(extraopt="", scale=1.0):
    global prefix

    histname = prefix + "evt_lvl_nvtxrewgt_histo_conecorrptvarbin_etavarbin_meas_mu"

    # Get histograms
    histo__tight = rebin2d_2dmu( tfile_data .Get(histname).Clone("Tight"), histname )
    histo__wj    = rebin2d_2dmu( tfile_wj   .Get(histname).Clone("WJets"), histname )
    histo__dy    = rebin2d_2dmu( tfile_dy   .Get(histname).Clone("DY"), histname )
    histo__ttbar = rebin2d_2dmu( tfile_ttbar.Get(histname).Clone("TTbar"), histname )
    if addvv: histo__vv    = rebin2d_2dmu( tfile_vv   .Get(histname).Clone("VV")   , histname )
    histo__tight.SetDirectory( 0 )
    histo__wj   .SetDirectory( 0 )
    histo__dy   .SetDirectory( 0 )
    histo__ttbar.SetDirectory( 0 )
    if addvv: histo__vv   .SetDirectory( 0 )
    histo__tight.Print( "all" )
    histo__wj   .Print( "all" )
    histo__dy   .Print( "all" )
    histo__ttbar.Print( "all" )
    if addvv: histo__vv   .Print( "all" )
    histo__tight.Add( histo__wj, -scale )
    histo__tight.Add( histo__dy, -scale )
    histo__tight.Add( histo__ttbar, -scale )
    if addvv: histo__tight.Add( histo__vv   , -scale )
    histo__tight.SetLineColor( 1 )
    histo__tight.SetFillColor( r.kGray + 1 )

    histname = prefix + "evt_lvl_nvtxrewgt_histo_conecorrptvarbin_etavarbin_loose_meas_mu"

    # Get histograms
    histo__loose = rebin2d_2dmu( tfile_data .Get(histname).Clone("Loose"), histname )
    histo__wj    = rebin2d_2dmu( tfile_wj   .Get(histname).Clone("WJets"), histname )
    histo__dy    = rebin2d_2dmu( tfile_dy   .Get(histname).Clone("DY"), histname )
    histo__ttbar = rebin2d_2dmu( tfile_ttbar.Get(histname).Clone("TTbar"), histname )
    if addvv: histo__vv    = rebin2d_2dmu( tfile_vv   .Get(histname).Clone("VV")   , histname )
    histo__loose.SetDirectory( 0 )
    histo__wj   .SetDirectory( 0 )
    histo__dy   .SetDirectory( 0 )
    histo__ttbar.SetDirectory( 0 )
    if addvv: histo__vv   .SetDirectory( 0 )
    histo__loose.Add( histo__wj, -scale )
    histo__loose.Add( histo__dy, -scale )
    histo__loose.Add( histo__ttbar, -scale )
    if addvv: histo__loose.Add( histo__vv   , -scale )
    histo__loose.SetLineColor( 1 )

    v_bkg_hists = r.vector("TH1*")()
    v_bkg_hists.push_back(histo__tight)

    ratio = histo__tight.Clone( "muon_fakerate_conecorrpt_v_eta" )
    ratio.SetDirectory( 0 )
    ratio.Divide( histo__loose )

    c1 = r.TCanvas( "", "", 0, 0, 800, 800 )
    c1.SetLogx()
    r.gStyle.SetPaintTextFormat( ".3f" )
    ratio.SetMaximum(0.5)
    ratio.SetMinimum(-0.15)
    ratio.SetMarkerColor(1)
    ratio.SetMarkerSize(0.3)
    ratio.SetContour(100)
    ratio.GetXaxis().SetMoreLogLabels()
    ratio.Draw( "colztexte" )
    c1.SaveAs( "plots/2dmuonfakerate.pdf" )
    return ratio

###################################################################################################
# Plot 2D distribution (default : data, W, DY. optional : add QCD, or ttbar)
#
def draw_qcd_fakerate_mu_2d(extraopt="", scale=1.0):
    global prefix

    histname = prefix + "evt_lvl_nvtxrewgt_histo_conecorrptvarbin_etavarbin_meas_mu"

    # Get histograms
    histo__tight = rebin2d( tfile_qcd_mu .Get(histname).Clone("Tight"), histname )
    histo__tight.SetDirectory( 0 )
    histo__tight.SetLineColor( 1 )
    histo__tight.SetFillColor( r.kGray + 1 )

    histname = prefix + "evt_lvl_nvtxrewgt_histo_conecorrptvarbin_etavarbin_loose_meas_mu"

    # Get histograms
    histo__loose = rebin2d( tfile_qcd_mu .Get(histname).Clone("Loose"), histname )
    histo__loose.SetDirectory( 0 )
    histo__loose.SetLineColor( 1 )

    v_bkg_hists = r.vector("TH1*")()
    v_bkg_hists.push_back(histo__tight)

    ratio = histo__tight.Clone( "muon_fakerate_conecorrpt_v_eta_qcd" )
    ratio.SetDirectory( 0 )
    ratio.Divide( histo__loose )

    c1 = r.TCanvas( "", "", 0, 0, 800, 800 )
    c1.SetLogx()
    c1.SetFillColor( -1 )
    r.gPad.Update()
    r.gPad.SetBottomMargin( 0.15 )
    r.gPad.SetLeftMargin( 0.15 )
    r.gStyle.SetPaintTextFormat( ".3f" )
    ratio.SetMaximum(0.5)
    ratio.SetMinimum(-0.15)
    ratio.SetMarkerColor(1)
    ratio.SetMarkerSize(0.8)
    ratio.SetContour(100)
    ratio.GetXaxis().SetMoreLogLabels()
    ratio.SetMarkerColor(1)
    ratio.SetMarkerSize(1.8)
    ratio.SetContour(100)
    ratio.GetXaxis().SetMoreLogLabels()
    ratio.GetXaxis().SetTitle( "p^{corr}_{T,#mu} [GeV]" )
    ratio.GetYaxis().SetTitle( "#eta" )
    ratio.GetXaxis().SetTitleSize( 0.05 )
    ratio.GetYaxis().SetTitleSize( 0.05 )
    ratio.GetXaxis().SetLabelSize( 0.05 )
    ratio.GetYaxis().SetLabelSize( 0.05 )
    ratio.GetXaxis().SetTitleOffset( 1.4 )
    ratio.GetYaxis().SetTitleOffset( 1.4 )
    ratio.Draw( "colztexte" )
    c1.SaveAs( "plots/2dmuonfakerate_qcd.pdf" )
    c1.SaveAs( "plots/2dmuonfakerate_qcd.png" )
    return ratio

###################################################################################################
# Plot 2D distribution (default : data, W, DY. optional : add QCD, or ttbar)
#
def draw_nvtxrewgt_fakerate_el_2d(extraopt="", scale=1.0):
    global prefix

    histname = prefix + "evt_lvl_nvtxrewgt_histo_conecorrptvarbin_etavarbin_meas_el"

    # Get histograms
    histo__tight = rebin2d( tfile_data .Get(histname).Clone("Tight"), histname )
    histo__wj    = rebin2d( tfile_wj   .Get(histname).Clone("WJets"), histname )
    histo__dy    = rebin2d( tfile_dy   .Get(histname).Clone("DY"), histname )
    histo__ttbar = rebin2d( tfile_ttbar.Get(histname).Clone("TTbar"), histname )
    if addvv: histo__vv    = rebin2d( tfile_vv   .Get(histname).Clone("VV")   , histname )
    histo__tight.SetDirectory( 0 )
    histo__wj   .SetDirectory( 0 )
    histo__dy   .SetDirectory( 0 )
    histo__ttbar.SetDirectory( 0 )
    if addvv: histo__vv   .SetDirectory( 0 )
    histo__tight.Add( histo__wj, -scale )
    histo__tight.Add( histo__dy, -scale )
    histo__tight.Add( histo__ttbar, -scale )
    if addvv: histo__tight.Add( histo__vv   , -scale )
    histo__tight.SetLineColor( 1 )
    histo__tight.SetFillColor( r.kGray + 1 )

    histname = prefix + "evt_lvl_nvtxrewgt_histo_conecorrptvarbin_etavarbin_loose_meas_el"

    # Get histograms
    histo__loose = rebin2d( tfile_data .Get(histname).Clone("Loose"), histname )
    histo__wj    = rebin2d( tfile_wj   .Get(histname).Clone("WJets"), histname )
    histo__dy    = rebin2d( tfile_dy   .Get(histname).Clone("DY"), histname )
    histo__ttbar = rebin2d( tfile_ttbar.Get(histname).Clone("TTbar"), histname )
    if addvv: histo__vv    = rebin2d( tfile_vv   .Get(histname).Clone("VV")   , histname )
    histo__loose.SetDirectory( 0 )
    histo__wj   .SetDirectory( 0 )
    histo__dy   .SetDirectory( 0 )
    histo__ttbar.SetDirectory( 0 )
    if addvv: histo__vv   .SetDirectory( 0 )
    histo__loose.Add( histo__wj, -scale )
    histo__loose.Add( histo__dy, -scale )
    histo__loose.Add( histo__ttbar, -scale )
    if addvv: histo__loose.Add( histo__vv   , -scale )
    histo__loose.SetLineColor( 1 )

    v_bkg_hists = r.vector("TH1*")()
    v_bkg_hists.push_back(histo__tight)

    ratio = histo__tight.Clone( "elec_fakerate_conecorrpt_v_eta" )
    ratio.SetDirectory( 0 )
    ratio.Divide( histo__loose )

    c1 = r.TCanvas( "", "", 0, 0, 800, 800 )
    c1.SetLogx()
    r.gStyle.SetPaintTextFormat( ".3f" )
    ratio.SetMaximum(0.5)
    ratio.SetMarkerColor(1)
    ratio.SetMarkerSize(0.8)
    ratio.SetContour(100)
    ratio.GetXaxis().SetMoreLogLabels()
    ratio.Draw( "colztexte" )
    c1.SaveAs( "plots/2delectronfakerate.pdf" )
    c1.SaveAs( "plots/2delectronfakerate.png" )
    return ratio

###################################################################################################
# Plot 2D distribution (default : data, W, DY. optional : add QCD, or ttbar)
#
def draw_qcd_fakerate_el_2d(extraopt="", scale=1.0):
    global prefix

    histname = prefix + "evt_lvl_nvtxrewgt_histo_conecorrptvarbin_etavarbin_meas_el"

    # Get histograms
    histo__tight = rebin2d( tfile_qcd_el .Get(histname).Clone("Tight"), histname )
    histo__tight.SetDirectory( 0 )
    histo__tight.SetLineColor( 1 )
    histo__tight.SetFillColor( r.kGray + 1 )

    histname = prefix + "evt_lvl_nvtxrewgt_histo_conecorrptvarbin_etavarbin_loose_meas_el"

    # Get histograms
    histo__loose = rebin2d( tfile_qcd_el .Get(histname).Clone("Loose"), histname )
    histo__loose.SetDirectory( 0 )
    histo__loose.SetLineColor( 1 )

    v_bkg_hists = r.vector("TH1*")()
    v_bkg_hists.push_back(histo__tight)

    ratio = histo__tight.Clone( "elec_fakerate_conecorrpt_v_eta_qcd" )
    ratio.SetDirectory( 0 )
    ratio.Divide( histo__loose )

    c1 = r.TCanvas( "", "", 0, 0, 800, 800 )
    c1.SetLogx()
    c1.SetFillColor( -1 )
    r.gPad.Update()
    r.gPad.SetBottomMargin( 0.15 )
    r.gPad.SetLeftMargin( 0.15 )
    r.gStyle.SetPaintTextFormat( ".3f" )
    ratio.SetMaximum(0.5)
    ratio.SetMarkerColor(1)
    ratio.SetMarkerSize(0.8)
    ratio.SetContour(100)
    ratio.GetXaxis().SetMoreLogLabels()
    ratio.SetMarkerColor(1)
    ratio.SetMarkerSize(1.8)
    ratio.SetContour(100)
    ratio.GetXaxis().SetMoreLogLabels()
    ratio.GetXaxis().SetTitle( "p^{corr}_{T,e} [GeV]" )
    ratio.GetYaxis().SetTitle( "#eta" )
    ratio.GetXaxis().SetTitleSize( 0.05 )
    ratio.GetYaxis().SetTitleSize( 0.05 )
    ratio.GetXaxis().SetLabelSize( 0.05 )
    ratio.GetYaxis().SetLabelSize( 0.05 )
    ratio.GetXaxis().SetTitleOffset( 1.4 )
    ratio.GetYaxis().SetTitleOffset( 1.4 )
    ratio.Draw( "colztexte" )
    c1.SaveAs( "plots/2delectronfakerate_qcd.pdf" )
    c1.SaveAs( "plots/2delectronfakerate_qcd.png" )
    return ratio

def main() :

    global prefix

    e8ips = draw("mll_e8i", "mll_e8i", "--xNbin 1 --reverseRatio" )[0].GetBinContent( 1 )
    e17ips = draw("mll_e17i", "mll_e17i", "--xNbin 1 --reverseRatio" )[0].GetBinContent( 1 )
    m8ips = draw("mll_m8i", "mll_m8i", "--xNbin 1 --reverseRatio" )[0].GetBinContent( 1 )
    m17ips = draw("mll_m17i", "mll_m17i", "--xNbin 1 --reverseRatio" )[0].GetBinContent( 1 )

    draw("mll_e8i", "mll_e8i", "--reverseRatio" )
    draw("mll_e17i", "mll_e17i", "--reverseRatio" )
    draw("mll_m8i", "mll_m8i", "--reverseRatio" )
    draw("mll_m17i", "mll_m17i", "--reverseRatio" )

    draw("histo_pt_el", "p_{T,e} [GeV]", "--saveMainPad --legendOnMainPad --Minimum 1000 --xNbin 200 --reverseRatio" )
    draw("histo_pt_mu", "p_{T,#mu} [GeV]", "--saveMainPad --legendOnMainPad --Minimum 1000 --xNbin 200 --reverseRatio" )
    draw("histo_nowgt_pt_el", "p_{T,e} [GeV]", "--xNbin 200 --saveMainPad --legendOnMainPad --Minimum 100 --reverseRatio" )
    draw("histo_nowgt_pt_mu", "p_{T,#mu} [GeV]", "--xNbin 200 --saveMainPad --legendOnMainPad --Minimum 100 --reverseRatio" )

    draw(prefix + "evt_lvl_histo_nvtx_highpt50_mu", "N_{vtx} (#mu events)", "--saveMainPad --legendOnMainPad --legendOnMainPad --xNbin 1" )[0].Print("all")
    draw(prefix + "evt_lvl_histo_nvtx_highpt50_el", "N_{vtx} (e events)", "--saveMainPad --legendOnMainPad --legendOnMainPad" )[0].Print("all")
    draw(prefix + "evt_lvl_histo_nvtx_highpt50_mu", "N_{vtx} (#mu events)", "--saveMainPad --legendOnMainPad --legendOnMainPad" )[0].Print("all")

    draw(prefix + "evt_lvl_nvtxrewgt_histo_nvtx_highpt50_el", "N_{vtx} (e events)", "--saveMainPad --legendOnMainPad --legendOnMainPad"  )[0].Print("all")
    draw(prefix + "evt_lvl_nvtxrewgt_histo_nvtx_highpt50_mu", "N_{vtx} (#mu events)", "--saveMainPad --legendOnMainPad --legendOnMainPad"  )[0].Print("all")

    f = r.TFile( "plots/fakerate_pt_v_eta.root", "recreate" )
    #sf_nonrewgt_mu = draw("evt_lvl_histo_mt_cr_mu", "m_{T,#mu} [GeV]", "--printBkg --xNbin 5" )[0].GetBinContent( 3 )
    #sf_nonrewgt_mu_syst1 = draw("syst13_evt_lvl_histo_mt_cr_mu", "m_{T,#mu} [GeV]", "--xNbin 5" )[0].GetBinContent( 3 )
    #sf_nonrewgt_mu_syst2 = draw("syst14_evt_lvl_histo_mt_cr_mu", "m_{T,#mu} [GeV]", "--xNbin 5" )[0].GetBinContent( 3 )

    #sf_nonrewgt_el = draw("evt_lvl_histo_mt_cr_el", "m_{T,e} [GeV]", "--xNbin 5" )[0].GetBinContent( 3 )
    #sf_nonrewgt_el_syst1 = draw("syst13_evt_lvl_histo_mt_cr_el", "m_{T,e} [GeV]", "--xNbin 5" )[0].GetBinContent( 3 )
    #sf_nonrewgt_el_syst2 = draw("syst14_evt_lvl_histo_mt_cr_el", "m_{T,e} [GeV]", "--xNbin 5" )[0].GetBinContent( 3 )
    #print sf_nonrewgt_mu, sf_nonrewgt_mu_syst1, sf_nonrewgt_mu_syst2
    #print sf_nonrewgt_el, sf_nonrewgt_el_syst1, sf_nonrewgt_el_syst2

#    draw("evt_lvl_nvtxrewgt_histo_nvtx_cr_mu_0", "nvtx", "" )[0].GetBinContent( 3 )
#    draw("evt_lvl_nvtxrewgt_histo_nvtx_cr_mu_1", "nvtx", "" )[0].GetBinContent( 3 )
#    draw("evt_lvl_nvtxrewgt_histo_nvtx_cr_mu_2", "nvtx", "" )[0].GetBinContent( 3 )
#    draw("evt_lvl_nvtxrewgt_histo_nvtx_cr_mu_3", "nvtx", "" )[0].GetBinContent( 3 )
#    draw("evt_lvl_nvtxrewgt_histo_nvtx_cr_mu_4", "nvtx", "" )[0].GetBinContent( 3 )
#
#    sf_mu_0 = draw("evt_lvl_nvtxrewgt_histo_mt_cr_mu_0", "m_{T,#mu} [GeV]", "--xNbin 5" )[0].GetBinContent( 3 )
#    sf_mu_1 = draw("evt_lvl_nvtxrewgt_histo_mt_cr_mu_1", "m_{T,#mu} [GeV]", "--xNbin 5" )[0].GetBinContent( 3 )
#    sf_mu_2 = draw("evt_lvl_nvtxrewgt_histo_mt_cr_mu_2", "m_{T,#mu} [GeV]", "--xNbin 5" )[0].GetBinContent( 3 )
#    sf_mu_3 = draw("evt_lvl_nvtxrewgt_histo_mt_cr_mu_3", "m_{T,#mu} [GeV]", "--xNbin 5" )[0].GetBinContent( 3 )
#    sf_mu_4 = draw("evt_lvl_nvtxrewgt_histo_mt_cr_mu_4", "m_{T,#mu} [GeV]", "--xNbin 5" )[0].GetBinContent( 3 )
#    print "sf_mu individ ", sf_mu_0, sf_mu_1, sf_mu_2, sf_mu_3, sf_mu_4

    #sf_mu = draw("evt_lvl_nvtxrewgt_histo_mt_cr_mu", "m_{T,#mu} [GeV]", "--xNbin 5" )[0].GetBinContent( 3 )
    #sf_mu_syst1 = draw("syst13_evt_lvl_nvtxrewgt_histo_mt_cr_mu", "m_{T,#mu} [GeV]", "--xNbin 5" )[0].GetBinContent( 3 )
    #sf_mu_syst2 = draw("syst14_evt_lvl_nvtxrewgt_histo_mt_cr_mu", "m_{T,#mu} [GeV]", "--xNbin 5" )[0].GetBinContent( 3 )

#    print "HEREHEREHERE2"
#    sf_mu = draw("evt_lvl_nvtxrewgt_histo_mt_cr_mu", "m_{T,#mu} [GeV]", "--xNbin 5", printcontent=True )[0].GetBinContent( 4 )
#    sf_mu_syst1 = draw("syst13_evt_lvl_nvtxrewgt_histo_mt_cr_mu", "m_{T,#mu} [GeV]", "--xNbin 5", printcontent=True )[0].GetBinContent( 4 )
#    sf_mu_syst2 = draw("syst14_evt_lvl_nvtxrewgt_histo_mt_cr_mu", "m_{T,#mu} [GeV]", "--xNbin 5", printcontent=True )[0].GetBinContent( 4 )

    print "HEREHEREHERE2"
    sf_mu = draw("evt_lvl_nvtxrewgt_histo_mt_meas2_mu", "m_{T,#mu} [GeV]", "--xNbin 5", printcontent=True )[0].GetBinContent( 3 )
    sf_mu_syst1 = draw("syst13_evt_lvl_nvtxrewgt_histo_mt_meas2_mu", "m_{T,#mu} [GeV]", "--xNbin 5", printcontent=True )[0].GetBinContent( 3 )
    sf_mu_syst2 = draw("syst14_evt_lvl_nvtxrewgt_histo_mt_meas2_mu", "m_{T,#mu} [GeV]", "--xNbin 5", printcontent=True )[0].GetBinContent( 3 )

    #print "HEREHEREHERE2"
    #sf_mu = draw("evt_lvl_nvtxrewgt_histo_mt_meas2_mu", "m_{T,#mu} [GeV]", "--xNbin 1", printcontent=True )[0].GetBinContent( 1 )
    #sf_mu_syst1 = draw("syst13_evt_lvl_nvtxrewgt_histo_mt_meas2_mu", "m_{T,#mu} [GeV]", "--xNbin 1", printcontent=True )[0].GetBinContent( 1 )
    #sf_mu_syst2 = draw("syst14_evt_lvl_nvtxrewgt_histo_mt_meas2_mu", "m_{T,#mu} [GeV]", "--xNbin 1", printcontent=True )[0].GetBinContent( 1 )

    #sf_mu = 1.10000

    #print "HEREHEREHERE2"
    #sf_mu = draw("evt_lvl_histo_mt_meas2_mu", "m_{T,#mu} [GeV]", "--xNbin 5", printcontent=True )[0].GetBinContent( 3 )
    #sf_mu_syst1 = draw("syst13_evt_lvl_histo_mt_meas2_mu", "m_{T,#mu} [GeV]", "--xNbin 5", printcontent=True )[0].GetBinContent( 3 )
    #sf_mu_syst2 = draw("syst14_evt_lvl_histo_mt_meas2_mu", "m_{T,#mu} [GeV]", "--xNbin 5", printcontent=True )[0].GetBinContent( 3 )

    #print sf_mu, sf_mu_syst2, sf_mu_syst1
    #sf_mu = 1.21891299696

    sf_el = draw("evt_lvl_nvtxrewgt_histo_mt_cr_el", "m_{T,e} [GeV]", "--xNbin 5" )[0].GetBinContent( 3 )
    sf_el_syst1 = draw("syst13_evt_lvl_nvtxrewgt_histo_mt_cr_el", "m_{T,e} [GeV]", "--xNbin 5" )[0].GetBinContent( 3 )
    sf_el_syst2 = draw("syst14_evt_lvl_nvtxrewgt_histo_mt_cr_el", "m_{T,e} [GeV]", "--xNbin 5" )[0].GetBinContent( 3 )

    draw("evt_lvl_histo_mt_cr_mu", "m_{T,#mu} [GeV]", "" )
    draw("syst13_evt_lvl_histo_mt_cr_mu", "m_{T,#mu} [GeV]", "" )
    draw("syst14_evt_lvl_histo_mt_cr_mu", "m_{T,#mu} [GeV]", "" )

    draw("evt_lvl_histo_mt_cr_el", "m_{T,e} [GeV]", "" )
    draw("syst13_evt_lvl_histo_mt_cr_el", "m_{T,e} [GeV]", "" )
    draw("syst14_evt_lvl_histo_mt_cr_el", "m_{T,e} [GeV]", "" )

    draw("evt_lvl_histo_mt_lowmet_mu", "m_{T,#mu} [GeV]", "" )
    draw("syst13_evt_lvl_histo_mt_lowmet_mu", "m_{T,#mu} [GeV]", "" )
    draw("syst14_evt_lvl_histo_mt_lowmet_mu", "m_{T,#mu} [GeV]", "" )

    draw("evt_lvl_histo_mt_lowmet_el", "m_{T,e} [GeV]", "" )
    draw("syst13_evt_lvl_histo_mt_lowmet_el", "m_{T,e} [GeV]", "" )
    draw("syst14_evt_lvl_histo_mt_lowmet_el", "m_{T,e} [GeV]", "" )

    draw("evt_lvl_nvtxrewgt_histo_pt_cr_mu", "p_{T,#mu} [GeV]", "--legendOnMainPad --saveMainPad" )
    draw("evt_lvl_nvtxrewgt_histo_pt_loose_cr_mu", "p_{T,#mu} [GeV]", "--legendOnMainPad --saveMainPad" )
    draw("evt_lvl_nvtxrewgt_histo_ptratio_loose_cr_mu", "p_{T,ratio #mu} [GeV]", "--legendOnMainPad --saveMainPad" )
    draw("evt_lvl_nvtxrewgt_histo_reliso_loose_cr_mu", "reliso", "--legendOnMainPad --saveMainPad" )
    draw("evt_lvl_nvtxrewgt_histo_mt_loose_cr_mu", "m_{T,#mu} [GeV]", "--legendOnMainPad --saveMainPad" )
    draw("evt_lvl_nvtxrewgt_histo_met_loose_cr_mu", "MET [GeV]", "--legendOnMainPad --saveMainPad" )
    draw("evt_lvl_nvtxrewgt_histo_nvtx_loose_cr_mu", "nvtx", "--legendOnMainPad --saveMainPad" )
    draw("evt_lvl_nvtxrewgt_histo_conecorrptvarbin_fine_meas_mu", "p^{corr}_{T,#mu} [GeV]" , "--data_DrawOpt ep --divideByBinWidth --saveMainPad --legendOnMainPad --Minimum 0 --Maximum 50000 --reverseRatio" , sf_mu)

    draw("evt_lvl_nvtxrewgt_histo_mt_cr_mu", "m_{T,#mu} [GeV]", "--legendOnMainPad --saveMainPad" )
    draw("syst13_evt_lvl_nvtxrewgt_histo_mt_cr_mu", "m_{T,#mu} [GeV]", "" )
    draw("syst14_evt_lvl_nvtxrewgt_histo_mt_cr_mu", "m_{T,#mu} [GeV]", "" )

    draw("evt_lvl_nvtxrewgt_histo_mt_cr_el", "m_{T,e} [GeV]", "--legendOnMainPad --saveMainPad" )
    draw("syst13_evt_lvl_nvtxrewgt_histo_mt_cr_el", "m_{T,e} [GeV]", "" )
    draw("syst14_evt_lvl_nvtxrewgt_histo_mt_cr_el", "m_{T,e} [GeV]", "" )

    draw("evt_lvl_nvtxrewgt_histo_trf_mu"       , "Extrapolation", "--printYieldsTable --legendOnMainPad --saveMainPad" )
    draw("syst13_evt_lvl_nvtxrewgt_histo_trf_mu", "Extrapolation", "--printYieldsTable --legendOnMainPad --saveMainPad" )
    draw("syst14_evt_lvl_nvtxrewgt_histo_trf_mu", "Extrapolation", "--printYieldsTable --legendOnMainPad --saveMainPad" )

    draw("evt_lvl_nvtxrewgt_histo_trf_el"       , "Extrapolation", "--printYieldsTable --legendOnMainPad --saveMainPad" )
    draw("syst13_evt_lvl_nvtxrewgt_histo_trf_el", "Extrapolation", "--printYieldsTable --legendOnMainPad --saveMainPad" )
    draw("syst14_evt_lvl_nvtxrewgt_histo_trf_el", "Extrapolation", "--printYieldsTable --legendOnMainPad --saveMainPad" )

    draw("evt_lvl_nvtxrewgt_histo_conecorrptvarbin_loose_meas_el" , "p^{corr}_{T,e} [GeV]"   , "--data_DrawOpt ep --divideByBinWidth --saveMainPad --legendOnMainPad --Minimum 100 --reverseRatio" , sf_el)
    draw("evt_lvl_nvtxrewgt_histo_conecorrptvarbin_loose_meas_mu" , "p^{corr}_{T,#mu} [GeV]" , "--data_DrawOpt ep --divideByBinWidth --saveMainPad --legendOnMainPad --Minimum 100 --reverseRatio" , sf_mu)
    draw("evt_lvl_nvtxrewgt_histo_conecorrptvarbin_meas_el"       , "p^{corr}_{T,e} [GeV]"   , "--data_DrawOpt ep --divideByBinWidth --saveMainPad --legendOnMainPad --Minimum 100 --reverseRatio" , sf_el)
    draw("evt_lvl_nvtxrewgt_histo_conecorrptvarbin_meas_mu"       , "p^{corr}_{T,#mu} [GeV]" , "--data_DrawOpt ep --divideByBinWidth --saveMainPad --legendOnMainPad --Minimum 100 --reverseRatio" , sf_mu)

    #fakerate_nvtxrewgt_mu_data = draw_nvtxrewgt_fakerate_mu( "--legendOnMainPad --data_DrawOpt ep --onlyLog --onlyLog", sf_mu, "evt_lvl_nvtxrewgt_histo_conecorrptvarbin_fine_meas_mu", "evt_lvl_nvtxrewgt_histo_conecorrptvarbin_fine_loose_meas_mu", "fakerate_nvtxrewgt_fine_mu_data")
    fakerate_nvtxrewgt_mu_data = draw_nvtxrewgt_fakerate_mu( "--legendOnMainPad --data_DrawOpt ep --onlyLog --onlyLog", sf_mu)
    #fakerate_nvtxrewgt_mu_data = draw_fakerate_mu( "--legendOnMainPad --data_DrawOpt ep --onlyLog --onlyLog", sf_mu)
    prefix = "syst13_"
    fakerate_nvtxrewgt_mu_data_syst1 = draw_nvtxrewgt_fakerate_mu("--onlyLog --plotOutputName plots/fakerate_nvtxrewgt_mu_data_syst1", sf_mu_syst1)
    #fakerate_nvtxrewgt_mu_data_syst1 = draw_fakerate_mu("--onlyLog --plotOutputName plots/fakerate_nvtxrewgt_mu_data_syst1", sf_mu_syst1)
    #fakerate_nvtxrewgt_mu_data_syst1 = draw_nvtxrewgt_fakerate_mu("--onlyLog --plotOutputName plots/fakerate_nvtxrewgt_mu_data_syst1", sf_mu * 1.0750804349)
    prefix = "syst14_"
    fakerate_nvtxrewgt_mu_data_syst2 = draw_nvtxrewgt_fakerate_mu("--onlyLog --plotOutputName plots/fakerate_nvtxrewgt_mu_data_syst2", sf_mu_syst2)
    #fakerate_nvtxrewgt_mu_data_syst2 = draw_fakerate_mu("--onlyLog --plotOutputName plots/fakerate_nvtxrewgt_mu_data_syst2", sf_mu_syst2)
    #fakerate_nvtxrewgt_mu_data_syst2 = draw_nvtxrewgt_fakerate_mu("--onlyLog --plotOutputName plots/fakerate_nvtxrewgt_mu_data_syst2", sf_mu * 0.9266846409)
    prefix = ""
    fakerate_nvtxrewgt_el_data = draw_nvtxrewgt_fakerate_el("--legendOnMainPad --data_DrawOpt ep --onlyLog --onlyLog", sf_el)
    prefix = "syst13_"
    fakerate_nvtxrewgt_el_data_syst1 = draw_nvtxrewgt_fakerate_el("--onlyLog --plotOutputName plots/fakerate_nvtxrewgt_el_data_syst1", sf_el_syst1)
    prefix = "syst14_"
    fakerate_nvtxrewgt_el_data_syst2 = draw_nvtxrewgt_fakerate_el("--onlyLog --plotOutputName plots/fakerate_nvtxrewgt_el_data_syst2", sf_el_syst2)
    prefix = ""
    fakerate_mu_qcd  = draw_fakerate_qcd_mu("--saveMainPad --legendOnMainPad --data_DrawOpt ep --onlyLog --onlyLog")
    fakerate_el_qcd  = draw_fakerate_qcd_el("--saveMainPad --legendOnMainPad --data_DrawOpt ep --onlyLog --onlyLog")

    print sf_mu, sf_mu_syst1, sf_mu_syst2

    print "without systematics"
    fakerate_nvtxrewgt_mu_data[0].Print("all")
    fakerate_nvtxrewgt_el_data[0].Print("all")
    print "without systematics"

    fakerate_nvtxrewgt_mu_data[0].SetName( "Data" )
    fakerate_nvtxrewgt_el_data[0].SetName( "Data" )
    fakerate_mu_qcd[0].SetName( "QCD MC" )
    fakerate_el_qcd[0].SetName( "QCD MC" )

    syst_data_vectors = r.vector( "TH1*" )()
    syst_data_vectors.push_back( fakerate_nvtxrewgt_el_data_syst1[0] )
    syst_data_vectors.push_back( fakerate_nvtxrewgt_el_data_syst2[0] )
    fakerate_nvtxrewgt_el_data_systsum = r.getSystByMaxDiff( fakerate_nvtxrewgt_el_data[0], syst_data_vectors )
    syst_data_vectors = r.vector( "TH1*" )()
    syst_data_vectors.push_back( fakerate_nvtxrewgt_mu_data_syst1[0] )
    syst_data_vectors.push_back( fakerate_nvtxrewgt_mu_data_syst2[0] )
    fakerate_nvtxrewgt_mu_data_systsum = r.getSystByMaxDiff( fakerate_nvtxrewgt_mu_data[0], syst_data_vectors )

    print "systematics only"
    fakerate_nvtxrewgt_mu_data_systsum.Print("all")
    fakerate_nvtxrewgt_el_data_systsum.Print("all")
    print "systematics only"

    data_vectors = r.vector( "TH1*" )()
    bkg_vectors  = r.vector( "TH1*" )()
    sig_vectors  = r.vector( "TH1*" )()
    syst_data_vectors = r.vector( "TH1*" )()
    syst_bkg_vectors  = r.vector( "TH1*" )()
    syst_sig_vectors  = r.vector( "TH1*" )()
    data_vectors.push_back( fakerate_nvtxrewgt_mu_data[0] )
    syst_data_vectors.push_back( fakerate_nvtxrewgt_mu_data_systsum )
    bkg_vectors.push_back( fakerate_mu_qcd[0] )
    syst_bkg_vectors.push_back( fakerate_mu_qcd[0] )
    r.plotmaker( """
                  --yTitle Fake rate
                  --xTitle p^{corr}_{T,#mu} [GeV]
                  --plotOutputName plots/fakerate_nvtxrewgt_mu_full
                  --ratio_DrawOpt ep
                  --ratio_Maximum 2.50
                  --ratio_Minimum 0.0
                  --Maximum 0.15
                  --printYieldsTable
                  --autoStack
                  --systByDiff
                  --onlyLin
                  --saveMainPad
                  --userhack0
                  --stack_DrawOpt hist
                  --legend_bkgDrawOpt lf
                  --legendOnMainPad
                  """,
                  data_vectors, bkg_vectors, sig_vectors, syst_data_vectors, syst_bkg_vectors, syst_sig_vectors )[0].Print("all")

    data_vectors = r.vector( "TH1*" )()
    bkg_vectors  = r.vector( "TH1*" )()
    sig_vectors  = r.vector( "TH1*" )()
    syst_data_vectors = r.vector( "TH1*" )()
    syst_bkg_vectors  = r.vector( "TH1*" )()
    syst_sig_vectors  = r.vector( "TH1*" )()
    data_vectors.push_back( fakerate_nvtxrewgt_el_data[0] )
    syst_data_vectors.push_back( fakerate_nvtxrewgt_el_data_systsum )
    bkg_vectors.push_back( fakerate_el_qcd[0] )
    syst_bkg_vectors.push_back( fakerate_el_qcd[0] )
    r.plotmaker( """
                  --yTitle Fake rate
                  --xTitle p^{corr}_{T,e} [GeV]
                  --plotOutputName plots/fakerate_nvtxrewgt_el_full
                  --ratio_DrawOpt ep
                  --ratio_Maximum 2.00
                  --ratio_Minimum 0.0
                  --printYieldsTable
                  --Maximum 0.37
                  --autoStack
                  --systByDiff
                  --onlyLin
                  --saveMainPad
                  --userhack0
                  --stack_DrawOpt hist
                  --legend_bkgDrawOpt lf
                  --legendOnMainPad
                  """,
                  data_vectors, bkg_vectors, sig_vectors, syst_data_vectors, syst_bkg_vectors, syst_sig_vectors )[0].Print("all")

    prefix = "syst13_"
    fakerate_mu_2d_sys1 = draw_nvtxrewgt_fakerate_mu_2d("", sf_mu)
    prefix = "syst14_"
    fakerate_mu_2d_sys2 = draw_nvtxrewgt_fakerate_mu_2d("", sf_mu)
    prefix = ""
    fakerate_mu_2d_nom  = draw_nvtxrewgt_fakerate_mu_2d("", sf_mu)

    v = r.vector( "TH1*" )()
    v.push_back( fakerate_mu_2d_sys1 )
    v.push_back( fakerate_mu_2d_sys2 )
    fakerate_mu_2d_totsys = r.getSyst2DByMaxDiff( fakerate_mu_2d_nom, v )
    fakerate_mu_2d = r.hist2DWithFullError( fakerate_mu_2d_nom, fakerate_mu_2d_totsys )

    r.gStyle.SetPaintTextFormat( ".3f" )

    prefix = "syst13_"
    fakerate_el_2d_sys1 = draw_nvtxrewgt_fakerate_el_2d("", sf_el)
    prefix = "syst14_"
    fakerate_el_2d_sys2 = draw_nvtxrewgt_fakerate_el_2d("", sf_el)
    prefix = ""
    fakerate_el_2d_nom  = draw_nvtxrewgt_fakerate_el_2d("", sf_el)

    v = r.vector( "TH1*" )()
    v.push_back( fakerate_el_2d_sys1 )
    v.push_back( fakerate_el_2d_sys2 )
    fakerate_el_2d_totsys = r.getSyst2DByMaxDiff( fakerate_el_2d_nom, v )
    fakerate_el_2d = r.hist2DWithFullError( fakerate_el_2d_nom, fakerate_el_2d_totsys )

    print "2d"
    fakerate_mu_2d.Print( "all" )
    fakerate_el_2d.Print( "all" )
    print "2d"

    draw_qcd_fakerate_mu_2d()
    draw_qcd_fakerate_el_2d()


    ############
    r.gStyle.SetPadBottomMargin( 0.15 )
    r.gStyle.SetPadLeftMargin( 0.15 )

    fakerate_el_2d.SetMaximum(0.9)
    fakerate_el_2d.SetMinimum(-0.15)
    fakerate_el_2d.SetMarkerColor(1)
    fakerate_el_2d.SetMarkerSize(1.8)
    fakerate_el_2d.SetContour(100)
    fakerate_el_2d.GetXaxis().SetMoreLogLabels()
    fakerate_el_2d.GetXaxis().SetTitle( "p^{corr}_{T,e} [GeV]" )
    fakerate_el_2d.GetYaxis().SetTitle( "#eta" )
    fakerate_el_2d.GetXaxis().SetTitleSize( 0.05 )
    fakerate_el_2d.GetYaxis().SetTitleSize( 0.05 )
    fakerate_el_2d.GetXaxis().SetLabelSize( 0.05 )
    fakerate_el_2d.GetYaxis().SetLabelSize( 0.05 )
    fakerate_el_2d.GetXaxis().SetTitleOffset( 1.4 )
    fakerate_el_2d.GetYaxis().SetTitleOffset( 1.4 )
    c2 = r.TCanvas( "c2", "c2", 0, 0, 800, 800 )
    c2.SetFillColor( -1 )
    r.gPad.Update()
    r.gPad.SetBottomMargin( 0.15 )
    r.gPad.SetLeftMargin( 0.15 )
    c2.SetLogx()
    fakerate_el_2d.Draw( "coltexte" )
    c2.SaveAs( "plots/2delecfakerate_full.pdf" )
    c2.SaveAs( "plots/2delecfakerate_full.png" )

    fakerate_mu_2d.SetMaximum(0.7)
    fakerate_mu_2d.SetMinimum(-0.15)
    fakerate_mu_2d.SetMarkerColor(1)
    fakerate_mu_2d.SetMarkerSize(0.8)
    fakerate_mu_2d.SetContour(100)
    fakerate_mu_2d.GetXaxis().SetMoreLogLabels()
    fakerate_mu_2d.GetYaxis().SetNdivisions(605)
    fakerate_mu_2d.GetXaxis().SetTitle( "p^{corr}_{T,#mu} [GeV]" )
    fakerate_mu_2d.GetYaxis().SetTitle( "#eta" )
    fakerate_mu_2d.GetXaxis().SetTitleSize( 0.05 )
    fakerate_mu_2d.GetYaxis().SetTitleSize( 0.05 )
    fakerate_mu_2d.GetXaxis().SetLabelSize( 0.05 )
    fakerate_mu_2d.GetYaxis().SetLabelSize( 0.05 )
    fakerate_mu_2d.GetXaxis().SetTitleOffset( 1.4 )
    fakerate_mu_2d.GetYaxis().SetTitleOffset( 1.4 )
    c1 = r.TCanvas( "c1", "c1", 0, 0, 800, 800 )
    c1.SetFillColor( -1 )
    r.gPad.Update()
    r.gPad.SetBottomMargin( 0.15 )
    r.gPad.SetLeftMargin( 0.15 )
    c1.SetLogx()
    fakerate_mu_2d.Draw( "coltexte" )
    c1.SaveAs( "plots/2dmuonfakerate_full.pdf" )
    c1.SaveAs( "plots/2dmuonfakerate_full.png" )

    fakerate_mu_2d.Write()
    fakerate_el_2d.Write()
    draw_qcd_fakerate_mu_2d("").Write()
    draw_qcd_fakerate_el_2d("").Write()

    print ""
    print "Scale factor (mu) : ", sf_mu
    print "Scale factor (el) : ", sf_el
    print ""

    print ""
    print ""
    print " Trigger prescale table in AN"
    print "============================================================"
    print "\\begin{table}[htb]"
    print "\caption{\label{tab:auxtrig:prescale} Trigger prescale value derived from selecting $Z$ boson events.}"
    print "\centering"
    print "\\begin{tabular}{|l|l|l|}"
    print "\hline"
    print "Lepton Flavor & Trigger & rederived prescale value \\\\"
    print "\hline"
    print "\multirow{2}{*}{$\mu$} & HLT\_Mu8\_TrkIsoVVL  & %d \\\\"%(m8ips)
    print "                       & HLT\_Mu17\_TrkIsoVVL & %d \\\\"%(m17ips)
    print "\hline"
    print "\multirow{2}{*}{$e$}   & HLT\_Ele8\_CaloIdL\_TrackIdL\_IsoVL\_PFJet30  & %d \\\\"%(e8ips)
    print "                       & HLT\_Ele17\_CaloIdL\_TrackIdL\_IsoVL\_PFJet30 & %d \\\\"%(e17ips)
    print "\hline"
    print "\end{tabular}"
    print "\end{table}"
    print "e8i  prescale value = ", int(e8ips)
    print "e17i prescale value = ", int(e17ips)
    print "m8i  prescale value = ", int(m8ips)
    print "m17i prescale value = ", int(m17ips)
    print "============================================================"
    print ""
    print ""

    muvals = []
    elvals = []
    for ibin in xrange(1, fakerate_nvtxrewgt_mu_data[0].GetNbinsX()+1):
        nom = fakerate_nvtxrewgt_mu_data[0].GetBinContent(ibin)
        stat = fakerate_nvtxrewgt_mu_data[0].GetBinError(ibin)
        syst = fakerate_nvtxrewgt_mu_data_systsum.GetBinContent(ibin)
        muvals.append(nom)
        muvals.append(stat)
        muvals.append(syst)
    for ibin in xrange(1, fakerate_nvtxrewgt_el_data[0].GetNbinsX()+1):
        nom = fakerate_nvtxrewgt_el_data[0].GetBinContent(ibin)
        stat = fakerate_nvtxrewgt_el_data[0].GetBinError(ibin)
        syst = fakerate_nvtxrewgt_el_data_systsum.GetBinContent(ibin)
        elvals.append(nom)
        elvals.append(stat)
        elvals.append(syst)

    print ""
    print " 1D Fake rate table for AN"
    print "============================================================"
    print "\\begin{table}[htb]"
    print "    \caption{\label{tab:fakeratept} Fake rate as a function of cone-corrected $\pt$. The fake rate is presented in the following format: (fake rate)~$\pm$~(stat.)~$\pm$~(syst.)}"
    print "    \centering"
    print "    \\begin{tabular}{|l|c|c|c|c|}"
    print "        \hline"
    print "                 & $10 < \pt^{corr} < 20 $ & $20 < \pt^{corr} < 30 $ & $30 < \pt^{corr} < 50 $ & $\pt^{corr} > 50 $ \\\\"
    print "        \hline"
    print "        $\mu$    & {:.3f}\pm{:.3f}\pm{:.3f}$     & ${:.3f}\pm{:.3f}\pm{:.3f}$     & ${:.3f}\pm{:.3f}\pm{:.3f}$     & ${:.2f}\pm{:.2f}\pm{:.2f}$        \\\\".format(*muvals)
    print "        \hline"
    print "        $e$      & {:.3f}\pm{:.3f}\pm{:.3f}$     & ${:.3f}\pm{:.3f}\pm{:.3f}$     & ${:.3f}\pm{:.3f}\pm{:.3f}$     & ${:.2f}\pm{:.2f}\pm{:.2f}$       \\\\".format(*elvals)
    print "        \hline"
    print "    \end{tabular}"
    print "\end{table}"
    print "============================================================"
    print ""
    print ""

    vals = []
    for ybin in xrange(1, fakerate_mu_2d.GetNbinsY()+1):
        val = []
        for xbin in xrange(1, fakerate_mu_2d.GetNbinsX()+1):
            val.append(fakerate_mu_2d.GetBinContent(xbin, ybin))
            val.append(fakerate_mu_2d.GetBinError(xbin, ybin))
        vals.append(val)

    print ""
    print ""
    print " Fake rate muon table for AN"
    print "============================================================"
    print "\\begin{table}[htb]"
    print "    \caption{\label{tab:fakerateptetamu} Fake rate summary table for muons. The error includes both the statistical error and the systematic error.}"
    print "    \centering"
    print "    \\begin{tabular}{|l|c|c|c|c|}"
    print "        \hline"
    print "                         & $10 < \pt^{corr} < 20 $ & $20 < \pt^{corr} < 30 $ & $30 < \pt^{corr} < 50 $ & $\pt^{corr} > 50 $ \\\\"
    print "        \hline"
    print "        $0<|\eta|<1.2$   & ${:.2f} \pm{:.2f} $         & {:.2f}5 \pm{:.2f} $         & ${:.2f} \pm{:.2f} $         & {:.2f}8\pm{:.2f}$        \\\\".format(*vals[0])
    print "        \hline"
    print "        $1.2<|\eta|<2.1$ & ${:.2f} \pm{:.2f} $         & {:.2f}9 \pm{:.2f} $         & ${:.2f} \pm{:.2f} $         & {:.2f}2\pm{:.2f}$        \\\\".format(*vals[1])
    print "        \hline"
    print "        $2.1<|\eta|<2.4$ & ${:.2f} \pm{:.2f} $         & {:.2f}2 \pm{:.2f} $         & ${:.2f} \pm{:.2f} $         & {:.2f}2\pm{:.2f}$        \\\\".format(*vals[2])
    print "        \hline"
    print "    \end{tabular}"
    print "\end{table}"
    print "============================================================"

    vals = []
    for ybin in xrange(1, fakerate_el_2d.GetNbinsY()+1):
        val = []
        for xbin in xrange(1, fakerate_el_2d.GetNbinsX()+1):
            val.append(fakerate_el_2d.GetBinContent(xbin, ybin))
            val.append(fakerate_el_2d.GetBinError(xbin, ybin))
        vals.append(val)

    print ""
    print ""
    print " Fake rate elec table for AN"
    print "============================================================"
    print "\\begin{table}[htb]"
    print "    \caption{\label{tab:fakerateptetael} Fake rate summary table for electrons. The error includes both the statistical error and the systematic error.}"
    print "    \centering"
    print "    \\begin{tabular}{|l|c|c|c|c|}"
    print "        \hline"
    print "                           & $10 < \pt^{corr} < 20 $ & $20 < \pt^{corr} < 30 $ & $30 < \pt^{corr} < 50 $ & $\pt^{corr} > 50 $ \\\\"
    print "        \hline"
    print "        $0<|\eta|<0.8$     & {:.2f} \pm{:.2f}$         & ${:.2f}\pm{:.2f}$         & ${:.2f}\pm{:.2f}$         & ${:.2f}pm{:.2f}    \\\\".format(*vals[0])
    print "        \hline"                                                                                                                                            
    print "        $0.8<|\eta|<1.479$ & {:.2f} \pm{:.2f}$         & ${:.2f}\pm{:.2f}$         & ${:.2f}\pm{:.2f}$         & ${:.2f}pm{:.2f}    \\\\".format(*vals[1])
    print "        \hline"                                                                                                                                            
    print "        $1.479<|\eta|<2.5$ & {:.2f} \pm{:.2f}$         & ${:.2f}\pm{:.2f}$         & ${:.2f}\pm{:.2f}$         & ${:.2f}pm{:.2f}    \\\\".format(*vals[2])
    print "        \hline"
    print "    \end{tabular}"
    print "\end{table}"
    print "============================================================"

    return

if __name__ == "__main__":

    prefix = ""
    main()
    #sf_mu = draw("evt_lvl_nvtxrewgt_histo_mt_cr_mu", "m_{T,#mu} [GeV]", "--xNbin 4", printcontent=True )[0].GetBinContent( 2 )

    #draw("evt_lvl_histo_met_mtwindow_mu", "test" , "--ratio_Maximum 0.02 --ratio_Minimum 0")
    #draw("evt_lvl_histo_nvtx_mtwindow_mu", "test" , "--ratio_Maximum 0.02 --ratio_Minimum 0")[0].Print("all")

    #draw("evt_lvl_histo_deltaphi_meas2_mu", "test" , "")[0].Print("all")
    #draw("evt_lvl_histo_met_meas2_mu", "test" , "")[0].Print("all")
    #draw("evt_lvl_histo_mt_meas2_mu", "test" , "")[0].Print("all")
    draw("evt_lvl_nvtxrewgt_histo_deltaphi_meas2_mu", "test" , "")[0].Print("all")
    draw("evt_lvl_nvtxrewgt_histo_deltaphi_meas2_mu", "test" , "--xNbin 2")[0].Print("all")
    draw("evt_lvl_nvtxrewgt_histo_met_meas2_mu", "test" , "")[0].Print("all")
    draw("evt_lvl_nvtxrewgt_histo_mt_meas2_mu", "test" , "")[0].Print("all")
    #draw("evt_lvl_nvtxrewgt_histo_mt_meas2_mu", "test" , "--xNbin 5")[0].Print("all")
    #draw("syst13_evt_lvl_nvtxrewgt_histo_mt_meas2_mu", "test" , "--xNbin 5")[0].Print("all")
    #draw("syst14_evt_lvl_nvtxrewgt_histo_mt_meas2_mu", "test" , "--xNbin 5")[0].Print("all")
    #draw("evt_lvl_histo_mt_meas2_mu", "test" , "--xNbin 5")[0].Print("all")
    #draw("syst13_evt_lvl_histo_mt_meas2_mu", "test" , "--xNbin 5")[0].Print("all")
    #draw("syst14_evt_lvl_histo_mt_meas2_mu", "test" , "--xNbin 5")[0].Print("all")

    #draw_fakerate_qcd_mu("--saveMainPad --legendOnMainPad --data_DrawOpt ep --onlyLog --onlyLog", "evt_lvl_nvtxrewgt_histo_conecorrptvarbin_fine_meas_mu", "evt_lvl_nvtxrewgt_histo_conecorrptvarbin_fine_loose_meas_mu", "fakerate_mu_fine_qcd")
    #fakerate_nvtxrewgt_mu_data = draw_nvtxrewgt_fakerate_mu( "--legendOnMainPad --data_DrawOpt ep --onlyLog --onlyLog", 1.06, "evt_lvl_nvtxrewgt_histo_conecorrptvarbin_meas_mu", "evt_lvl_nvtxrewgt_histo_conecorrptvarbin_loose_meas_mu", "fakerate_nvtxrewgt_mu_data")
    #fakerate_nvtxrewgt_mu_data = draw_nvtxrewgt_fakerate_mu( "--legendOnMainPad --data_DrawOpt ep --onlyLog --onlyLog", 1.05, "evt_lvl_nvtxrewgt_histo_conecorrpt_meas_mu", "evt_lvl_nvtxrewgt_histo_conecorrpt_loose_meas_mu", "fakerate_nvtxrewgt_regbin_mu_data")

    #draw("evt_lvl_nvtxrewgt_histo_conecorrptvarbin_fine_meas_mu", "p^{corr}_{T,#mu} [GeV]" , "--data_DrawOpt ep --divideByBinWidth --saveMainPad --legendOnMainPad --Minimum 0 --Maximum 50000 --reverseRatio" , 1)

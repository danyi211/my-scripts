import ROOT as r
import numpy as np
import uproot
import awkward
import concurrent.futures
import time
import copy
import os
import logging


# branches to get
branch_suffix = 'v12'
ecal_veto_epAng = 'EcalVeto_{}/epAng_'.format(branch_suffix)
ecal_veto_recoilPx = 'EcalVeto_{}/recoilPx_'.format(branch_suffix)
ecal_veto_discValue = "EcalVeto_{}/discValue_".format(branch_suffix)

branched_to_get_ldmx = []
branched_to_get_ldmx.append(ecal_veto_epAng)
branched_to_get_ldmx.append(ecal_veto_discValue)
branched_to_get_flat = ['epAng']


def DrawDeltathetaHist(table_ldmx, table_flat, hname1, hname2, tag):
    a_discValue = np.array(table_ldmx[ecal_veto_discValue])
    a_epAng = np.array(table_ldmx[ecal_veto_epAng])
    a_epAng2 = np.array(table_flat['epAng'])

    # Delta theta_ep
    a_delta_epAng = a_epAng2 - a_epAng

    # Create histogram
    h_diff = r.TH1D(hname1, '', 150, -100., 50.) 
    h_diff_BDT = r.TH1D(hname2, '', 150, -100., 50.) # with BDT

    nevents = np.size(a_epAng)
    for i in range(nevents):
        h_diff.Fill(a_delta_epAng[i])
        if a_discValue[i] > 0.99:
            h_diff_BDT.Fill(a_delta_epAng[i])
    
    # Plotting
    can = r.TCanvas("MyCan", "",700, 650)
    can.Draw()
    r.gStyle.SetNdivisions(505,"xy")
    r.gStyle.SetOptStat(0)

    pad1 = r.TPad("pad1", "pad1", 0.00, 0.33, 1.00, 1.00)
    pad2 = r.TPad("pad2", "pad2", 0.00, 0.00, 1.00, 0.33)

    pad1.SetBottomMargin(0.05)
    pad1.SetBorderMode(0)
    pad1.SetLeftMargin(0.12)
    pad1.SetTickx(1)
    pad1.SetTicky(1)
    pad1.SetLogy(1)  # pad1 logY
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.4)
    pad2.SetBorderMode(0)
    pad2.SetLeftMargin(0.12)
    pad2.SetTickx(1)
    pad2.SetTicky(1)
    pad2.SetGridy()
    r.gStyle.SetOptTitle(r.kFALSE)
    r.gStyle.SetOptStat(0)
    pad1.Draw()
    pad2.Draw()
    pad1.cd()

    # Normalize to 1
    # factor = 1.
    # h_diff.Scale(factor / h_diff.Integral())
    # h_diff_BDT.Scale(factor / h_diff_BDT.Integral())

    xtitle = "#Delta epAng [degree]"
    h_diff.GetYaxis().SetTitle("Events")
    h_diff.Draw('hist')
    h_diff_BDT.SetLineColor(2)
    h_diff_BDT.Draw('hist same')

    l = r.TLegend(0.15, 0.70, 0.4, 0.85)
    l.AddEntry(h_diff, 'Trigger', 'lp')
    l.AddEntry(h_diff_BDT, 'BDT score > 0.99', 'lp')
    l.SetBorderSize(0)
    l.SetFillStyle(0)
    l.Draw()

    text = r.TPaveText( 0.7, 0.75, 0.9, 0.9, "BRNDC" )
    text.AddText("Signal {}".format(tag))
    text.SetTextSize(0.04)
    text.SetFillColor(0)
    text.SetTextAlign(12)
    text.SetBorderSize(0)
    text.SetFillStyle(0)
    text.Draw()

    ratio0 = copy.deepcopy(h_diff)
    ratio0.Divide(h_diff)
    ratio1 = copy.deepcopy(h_diff_BDT)
    ratio1.Divide(h_diff)
    ratio1.Scale(h_diff.Integral() / h_diff_BDT.Integral())

    pad2.cd()
    ratio1.GetXaxis().SetTitle(xtitle)
    ratio1.GetXaxis().SetLabelSize(0.09)
    ratio1.GetXaxis().SetTitleSize(0.12)
    ratio1.GetXaxis().SetTickLength(0.04)
    ratio1.GetXaxis().SetTitleOffset(1.1)
    ratio1.GetXaxis().SetTitleFont(42)
    ratio1.GetYaxis().SetTitle("Trigger / BDT")
    ratio1.GetYaxis().CenterTitle(1)
    ratio1.GetYaxis().SetLabelSize(0.05)
    ratio1.GetYaxis().SetTitleSize(0.09)
    ratio1.GetYaxis().SetTitleOffset(0.6)
    ratio1.GetYaxis().SetLabelOffset(0.014)
    ratio1.GetYaxis().SetNdivisions(204,r.kFALSE)
    ratio1.GetYaxis().SetTitleOffset(0.5)
    ratio1.SetMarkerStyle(20)
    ratio1.SetStats(False)
    # ratio1.SetMaximum(2.)
    # ratio1.SetMinimum(0.0)
    # ratio1.GetYaxis().SetMaxDigits(0)
    ratio1.SetLineColor(r.kRed)
    ratio1.SetMarkerColor(r.kRed)
    ratio1.SetMarkerSize(0.8)
    ratio1.Draw("E1")

    ratio0.SetMarkerStyle(20)
    ratio0.SetStats(False)
    ratio0.SetLineColor(r.kBlue)
    ratio0.SetMarkerColor(r.kBlue)
    ratio0.SetMarkerSize(0.8)
    # ratio0.Draw("same histo")

    s = 'signal'
    can.SaveAs('anglePlots/hist_Delta_epAng_BDTcompare_ratio_{}_{}.pdf'.format(s, tag))

def main():
    # input files
    ldmx_dir = '/home/pmasterson/events/v3.0.0_tag_standard_skimmed/signal/'
    flat_dir = '/home/danyi/ldmx/tracking/BDTsamples/xinyi/signal/'
    treename_ldmx = "LDMX_Events"
    treename_flat = "EcalVeto"

    trees_ldmx_0p001 = {}
    trees_ldmx_0p01 = {}
    trees_ldmx_0p1 = {}
    trees_ldmx_1 = {}
    trees_flat_0p001 = {}
    trees_flat_0p01 = {}
    trees_flat_0p1 = {}
    trees_flat_1 = {}

    for rootfile in os.listdir(flat_dir):
        fullname_flat = flat_dir + rootfile
        fullname_original = ldmx_dir + rootfile.replace('_flatout_unsorted','')
        # Find corresbonding original rootfile
        if not os.path.exists(fullname_original):
            logging.info('Original root file {} not found.'.format(fullname_original))
        elif rootfile.find("0.001") >= 0:
            # print("found "+ rootfile)
            trees_flat_0p001[fullname_flat] = treename_flat
            trees_ldmx_0p001[fullname_original] = treename_ldmx
        elif rootfile.find("0.01") >= 0:
            # print("found "+ rootfile)
            trees_flat_0p01[fullname_flat] = treename_flat
            trees_ldmx_0p01[fullname_original] = treename_ldmx
        elif rootfile.find("0.1") >= 0:
            # print("found "+ rootfile)
            trees_flat_0p1[fullname_flat] = treename_flat
            trees_ldmx_0p1[fullname_original] = treename_ldmx
        elif rootfile.find("1.0") >= 0:
            # print("found "+ rootfile)
            trees_flat_1[fullname_flat] = treename_flat
            trees_ldmx_1[fullname_original] = treename_ldmx
    
    print("0.001 GeV: ", trees_ldmx_0p001)
    print("0.01 GeV: ", trees_ldmx_0p01)
    print("0.1 GeV: ", trees_ldmx_0p1)
    print("1 GeV: ", trees_ldmx_1)
    
    # Load data
    table_ldmx_0p001 = uproot.concatenate(trees_ldmx_0p001, branched_to_get_ldmx)
    table_ldmx_0p01 = uproot.concatenate(trees_ldmx_0p01, branched_to_get_ldmx)
    table_ldmx_0p1 = uproot.concatenate(trees_ldmx_0p1, branched_to_get_ldmx)
    table_ldmx_1 = uproot.concatenate(trees_ldmx_1, branched_to_get_ldmx)
    table_flat_0p001 = uproot.concatenate(trees_flat_0p001, branched_to_get_flat)
    table_flat_0p01 = uproot.concatenate(trees_flat_0p01, branched_to_get_flat)
    table_flat_0p1 = uproot.concatenate(trees_flat_0p1, branched_to_get_flat)
    table_flat_1 = uproot.concatenate(trees_flat_1, branched_to_get_flat)
    print("Loaded {} signal 0.001 GeV events.".format(len(table_ldmx_0p001[ecal_veto_epAng])))
    print("Loaded {} signal 0.01 GeV events.".format(len(table_ldmx_0p01[ecal_veto_epAng])))
    print("Loaded {} signal 0.1 GeV events.".format(len(table_ldmx_0p1[ecal_veto_epAng])))
    print("Loaded {} signal 1.0 GeV events.".format(len(table_ldmx_1[ecal_veto_epAng])))

    # Histograms
    DrawDeltathetaHist(table_ldmx_0p001, table_flat_0p001, 'h_Dtheta_0p001', 'h_Dtheta_0p001_wBDT', '0p001GeV')
    DrawDeltathetaHist(table_ldmx_0p01, table_flat_0p01, 'h_Dtheta_0p01', 'h_Dtheta_0p01_wBDT', '0p01GeV')
    DrawDeltathetaHist(table_ldmx_0p1, table_flat_0p1, 'h_Dtheta_0p1', 'h_Dtheta_0p1_wBDT','0p1GeV')
    DrawDeltathetaHist(table_ldmx_1, table_flat_1, 'h_Dtheta_1', 'h_Dtheta_1_wBDT', '1GeV')
    
    
if __name__=="__main__":
    main()

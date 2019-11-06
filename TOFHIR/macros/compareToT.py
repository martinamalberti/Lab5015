#! /usr/bin/env python
import os
import shutil
import glob
import math
import array
import sys
import time


import ROOT
import CMS_lumi, tdrstyle

#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetLabelSize(0.04,'X')
ROOT.gStyle.SetLabelSize(0.04,'Y')
ROOT.gStyle.SetTitleSize(0.04,'X')
ROOT.gStyle.SetTitleSize(0.04,'Y')
ROOT.gStyle.SetTitleOffset(1.1,'X')
ROOT.gStyle.SetTitleOffset(1.2,'Y')
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)

vth1s = [0, 20]
vth2s = [10, 20, 30 , 40, 50]

f = {}
g = {}
fitFun = {}

canvas = ROOT.TCanvas('c_amp_vs_tot_vs_thresholds','c_amp_vs_tot_vs_thresholds')
canvas.SetGridx()
canvas.SetGridy()
hdummy = ROOT.TH2F('hdummy','', 100, 0, 350, 100, 0, 500)
hdummy.GetXaxis().SetTitle('measured ToT (ns)')
hdummy.GetYaxis().SetTitle('amplitude (a.u.)')
hdummy.Draw()

leg = ROOT.TLegend(0.15, 0.6, 0.45, 0.89)
leg.SetBorderSize(0)

for i1, vth1 in enumerate(vth1s):
    f[vth1] = {}
    g[vth1] = {}
    fitFun[vth1] = {}
    for i2, vth2 in enumerate(vth2s):
        fname = '../output_TotNonlinearity_vth1_%d_vth2_%d.root'%(vth1, vth2)
        if ( os.path.isfile(fname) == 0 ):
            continue
        f[vth1][vth2] = ROOT.TFile.Open(fname)
        g[vth1][vth2] = f[vth1][vth2].Get('g_amp_vs_tot')
        #if (vth1 ==  0 and vth2 ==10): g[vth1][vth2].RemovePoint(0)
        #if (vth1 ==  0 and vth2 ==20): g[vth1][vth2].RemovePoint(0)
        #if (vth1 == 20 and vth2 ==20): g[vth1][vth2].RemovePoint(0)
        
        #fitFun[vth1][vth2] = g[vth1][vth2].GetFunction('fCorr')

        #expo
        #fitFun[vth1][vth2] = ROOT.TF1('fitFun_%d_%d'%(vth1, vth2),'expo',0,1000)

        # exp
        #fitFun[vth1][vth2] = ROOT.TF1('fitFun_%d_%d'%(vth1, vth2),'[0] * exp([1]*x) + [2]',0,1000)
        #fitFun[vth1][vth2].SetParameter(0, 100)
        #fitFun[vth1][vth2].SetParameter(1, 0.015)
        #fitFun[vth1][vth2].SetParameter(2, -1)
        
        #exponential passing by 0,0
        fitFun[vth1][vth2] = ROOT.TF1('fitFun_%d_%d'%(vth1, vth2),'[0] * ( exp([1]*x) - 1 )',0,1000)
        fitFun[vth1][vth2].SetParameter(1, 0.015)

        # sum of exponential passing by 0,0
        #fitFun[vth1][vth2] = ROOT.TF1('fitFun_%d_%d'%(vth1, vth2),'[0] * ( exp([1]*x) -  exp([2]*x))',0,1000)
        #fitFun[vth1][vth2].SetParameter(0, 300);
        #fitFun[vth1][vth2].SetParameter(1, 0.015);
        #fitFun[vth1][vth2].SetParameter(2, 0.001);

        # polinomial passing by 0,0
        #fitFun[vth1][vth2] = ROOT.TF1('fitFun_%d_%d'%(vth1, vth2),'[0]*x + [1]*x*x + [2]*x*x*x',0,1000)

       
        fitFun[vth1][vth2].SetRange(0,1000)
        #fitFun[vth1][vth2].SetRange(g[vth1][vth2].GetX()[0]-10, g[vth1][vth2].GetX()[5]+10 )
        fitFun[vth1][vth2].SetLineColor(800+20*i2)
        
        if (vth1!=0):
            g[vth1][vth2].SetMarkerStyle(24)
            g[vth1][vth2].SetLineStyle(2)
            fitFun[vth1][vth2].SetLineStyle(2)
        if (vth1==0):
            g[vth1][vth2].SetMarkerStyle(20)
            g[vth1][vth2].SetMarkerSize(0.8)
        g[vth1][vth2].SetMarkerColor(800+20*i2)
        g[vth1][vth2].SetLineColor(800+20*i2)

        
        g[vth1][vth2].Fit(fitFun[vth1][vth2],'R')
        g[vth1][vth2].Draw('psame')

        leg.AddEntry(g[vth1][vth2],'vth1 = %d, vth2 = %d'%(vth1, vth2), 'PL')

leg.Draw()


for i1, vth1 in enumerate(vth1s):
    for i2, vth2 in enumerate(vth2s):
        if (vth2 in fitFun[vth1].keys()):
            print vth1, vth2, 'exp1 = %.04f +/-  %.04f'%( fitFun[vth1][vth2].GetParameter(1),fitFun[vth1][vth2].GetParError(1) )
            #print vth1, vth2, 'exp2 = %.04f +/-  %.04f'%( fitFun[vth1][vth2].GetParameter(2),fitFun[vth1][vth2].GetParError(2) )

raw_input('ok?')

canvas.SaveAs('../NonLinearityToT/'+canvas.GetName()+'.pdf')
canvas.SaveAs('../NonLinearityToT/'+canvas.GetName()+'.png')



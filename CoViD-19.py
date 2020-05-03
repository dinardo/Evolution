######################################
# Program to analyse CoViD-19 spread #
# Author: Mauro E. Dinardo           #
######################################

import csv
from datetime import datetime
from random   import seed, random, gauss
from math     import sqrt, log, exp

from pyWget    import saveDataFromURL
from evolution import evolution

from ROOT import gROOT, gStyle, gApplication, TGaxis, TCanvas, TGraphErrors, TGraph, TH1D, TH2D, TLine


def SetStyle():
    gROOT.SetStyle('Plain')
    gROOT.ForceStyle()
    gStyle.SetTextFont(42)

    gStyle.SetOptTitle(0)
    gStyle.SetOptFit(1112)
    gStyle.SetOptStat(1110)

    gStyle.SetPadRightMargin(0.08)
    gStyle.SetPadTopMargin(0.11)
    gStyle.SetPadBottomMargin(0.12)

    gStyle.SetTitleFont(42,'x')
    gStyle.SetTitleFont(42,'y')
    gStyle.SetTitleFont(42,'z')

    gStyle.SetTitleOffset(1.05,'x')
    gStyle.SetTitleOffset(1.00,'y')

    gStyle.SetTitleSize(0.05,'x')
    gStyle.SetTitleSize(0.05,'y')
    gStyle.SetTitleSize(0.05,'z')

    gStyle.SetLabelFont(42,'x')
    gStyle.SetLabelFont(42,'y')
    gStyle.SetLabelFont(42,'z')

    gStyle.SetLabelSize(0.05,'x')
    gStyle.SetLabelSize(0.05,'y')
    gStyle.SetLabelSize(0.05,'z')

    TGaxis.SetMaxDigits(3)
    gStyle.SetStatY(0.9)


def readDataFromFile(fileName, column, country = '', province = ''):
    myDict = {}

    with open(fileName) as csvFile:
        data = csv.reader(csvFile, delimiter = ',')
        isFirstRow = True

        for r in data:
            if not r:
                continue

            if isFirstRow == True:
                # print 'Column names are', r
                isFirstRow = False
                firstRow = r
            elif column != -1:
                myDict[r[0]] = float(r[column])
            elif r[1] == country and r[0] == province:
                for c in range(4,len(r)):
                    myDict[datetime.strptime(firstRow[c], '%m/%d/%y').date()] = float(r[c])
                break

    return myDict


def assignErrors(yValues):
    scaleError = 3
    return [sqrt(y) * scaleError for y in yValues]


def analyzeData(country, total, active, recovered, deaths, tStart, tStop, totalPopulation, symptomaticFraction, transmissionProbability, recoveryRate, doSmearing):
    ntuple = [tStart, tStop, totalPopulation, symptomaticFraction, transmissionProbability, recoveryRate]
    farFromMax = 0.95
    growthRate = 0.13
    carryingCapacity = 4e5


    ################
    # Active cases #
    ################
    myCanvActive = TCanvas('myCanvActive_' + country, 'Active cases ' + country)

    xValues    = [i         for i in range(len(active.keys()))          if i >= tStart and i <= tStop]
    yValues    = [active[k] for i,k in enumerate(sorted(active.keys())) if i >= tStart and i <= tStop]
    erryValues = assignErrors(yValues)

    myGraphActive = TGraphErrors()
    myGraphActive.SetMarkerStyle(20)
    for i in range(len(xValues)):
        myGraphActive.SetPoint(myGraphActive.GetN(), xValues[i], yValues[i])
        myGraphActive.SetPointError(myGraphActive.GetN()-1, 0, erryValues[i])
    myGraphActive.Draw('APE1')
    myGraphActive.GetHistogram().GetXaxis().SetTitle('Time (days)')
    myGraphActive.GetHistogram().GetYaxis().SetTitle('Active cases affected by CoViD-19')

    historyActive = 0.
    for i,k in enumerate(sorted(active.keys())):
        if i < tStart:
            historyActive += active[k]
    ntuple.extend([historyActive, xValues, yValues, erryValues])
    print '==> History active cases:', historyActive

    evActive = evolution([yValues[0], carryingCapacity, recoveryRate, growthRate], tStart, tStop, totalPopulation, symptomaticFraction, transmissionProbability, historyActive)
    evActive.runOptimization(xValues, yValues, erryValues, [2], doSmearing)
    evActive.evolve(evActive.tStop, evActive.parValues, True)
    if doSmearing == True:
        evActive.smearing()
    evActiveGraphN = evActive.getGraphN()
    evActiveGraphN.Draw('PL same')
    statActive = evActive.addStats(evActive.parNames, evActive.parValues)

    print '==> Active cases, history active cases * dt, p-infected, Carrying capacity', evActive.evolve(tStop, evActive.parValues), '@', tStop, 'day'

    now = TLine(len(active)-1, 0, len(active)-1, evActive.fitFun.GetMaximum(len(active)-1))
    now.SetLineColor(4)
    now.SetLineWidth(2)
    now.Draw('same')

    willbe = TLine(evActive.fitFun.GetMaximumX(), 0, evActive.fitFun.GetMaximumX(), evActive.fitFun.GetMaximum())
    willbe.SetLineColor(6)
    willbe.SetLineWidth(2)
    willbe.Draw('same')

    myCanvActive.SetGrid()
    myCanvActive.Modified()
    myCanvActive.Update()

    myCanvActiveR0 = TCanvas('myCanvActiveR0_' + country, 'R0 ' + country)

    evActiveGraphR0 = evActive.getGraphR0(evActiveGraphN)
    evActiveGraphR0.Draw('APL')
    evActiveGraphR0.GetHistogram().GetXaxis().SetTitle('Time (days)')
    evActiveGraphR0.GetHistogram().GetYaxis().SetTitle('R')

    myCanvActiveR0.SetGrid()
    myCanvActiveR0.Modified()
    myCanvActiveR0.Update()

    myCanvActiveP = TCanvas('myCanvActiveP_' + country, 'Probability infected ' + country)

    evActiveGraphP = evActive.getGraphPinfect()
    evActiveGraphP.Draw('APL')
    evActiveGraphP.GetHistogram().GetXaxis().SetTitle('Time (days)')
    evActiveGraphP.GetHistogram().GetYaxis().SetTitle('Probability of being infected')

    myCanvActiveP.SetGrid()
    myCanvActiveP.Modified()
    myCanvActiveP.Update()


    #################
    # CONTROL PLOTS #
    #################

    ###############
    # Total cases #
    ###############
    myCanvTotal = TCanvas('myCanvTotal_' + country, 'Total cases ' + country)
    myGraphTotal = TGraphErrors()
    myGraphTotal.SetMarkerStyle(20)

    for k in sorted(total.keys()):
        myGraphTotal.SetPoint(myGraphTotal.GetN(), myGraphTotal.GetN(), total[k])
        myGraphTotal.SetPointError(myGraphTotal.GetN()-1, 0, sqrt(total[k]))

    myGraphTotal.Draw('APE1')
    myGraphTotal.GetHistogram().GetXaxis().SetTitle('Time (days)')
    myGraphTotal.GetHistogram().GetYaxis().SetTitle('Total cases affected by CoViD-19')

    myCanvTotal.SetGrid()
    myCanvTotal.Modified()
    myCanvTotal.Update()


    ##########
    # Deaths #
    ##########
    myCanvDeaths = TCanvas('myCanvDeaths_' + country, 'Deaths ' + country)
    myGraphDeaths = TGraphErrors()
    myGraphDeaths.SetMarkerStyle(20)

    for k in sorted(deaths.keys()):
        myGraphDeaths.SetPoint(myGraphDeaths.GetN(), myGraphDeaths.GetN(), deaths[k])
        myGraphDeaths.SetPointError(myGraphDeaths.GetN()-1, 0, sqrt(deaths[k]))

    myGraphDeaths.Draw('APE1')
    myGraphDeaths.GetHistogram().GetXaxis().SetTitle('Time (days)')
    myGraphDeaths.GetHistogram().GetYaxis().SetTitle('Total deaths')

    myCanvDeaths.SetGrid()
    myCanvDeaths.Modified()
    myCanvDeaths.Update()


    ########################
    # Ratio deaths / total #
    ########################
    myCanvRatio01 = TCanvas('myCanvRatio01_' + country, 'Ratio ' + country)
    myGraphRatio01 = TGraphErrors()
    myGraphRatio01.SetMarkerStyle(20)

    for k in sorted(deaths.keys()):
        myGraphRatio01.SetPoint(myGraphRatio01.GetN(), myGraphRatio01.GetN(), deaths[k]/total[k])
        myGraphRatio01.SetPointError(myGraphRatio01.GetN()-1, 0, myGraphRatio01.GetY()[myGraphRatio01.GetN()-1] * sqrt(deaths[k]/(deaths[k]*deaths[k]) + total[k]/(total[k]*total[k])))

    myGraphRatio01.Draw('APE1')
    myGraphRatio01.GetHistogram().GetXaxis().SetTitle('Time (days)')
    myGraphRatio01.GetHistogram().GetYaxis().SetTitle('Total deaths / Total cases')

    myCanvRatio01.SetGrid()
    myCanvRatio01.Modified()
    myCanvRatio01.Update()


    ############################################
    # Ratio delta(deaths + recovered) / active #
    ############################################
    myCanvRatio02 = TCanvas('myCanvRatio02_' + country, 'Ratio ' + country)
    myGraphRatio02 = TGraphErrors()
    myGraphRatio02.SetMarkerStyle(20)

    sortedKeys = sorted(deaths.keys())
    for i,k in enumerate(sortedKeys[1:]):
        myGraphRatio02.SetPoint(myGraphRatio02.GetN(), myGraphRatio02.GetN()+1, ((deaths[k] - deaths[sortedKeys[i]] + recovered[k] - recovered[sortedKeys[i]]) / active[k]) if active[k] != 0 else 0)
        numerator = deaths[k] + deaths[sortedKeys[i]] + recovered[k] + recovered[sortedKeys[i]]
        myGraphRatio02.SetPointError(myGraphRatio02.GetN()-1, 0, (myGraphRatio02.GetY()[myGraphRatio02.GetN()-1] * sqrt(numerator/pow(numerator,2) + active[k]/(active[k]*active[k]))) if active[k] != 0 else 0)

    myGraphRatio02.Draw('AP')
    myGraphRatio02.GetHistogram().GetXaxis().SetTitle('Time (days)')
    myGraphRatio02.GetHistogram().GetYaxis().SetTitle('#Delta Recovered (alive + dead) / Active cases')

    myCanvRatio02.SetGrid()
    myCanvRatio02.Modified()
    myCanvRatio02.Update()

    print '==> Doubling time:', round(log(2)/evActive.parValues[3],1), 'days'

    return [ntuple,
            myCanvTotal,   myGraphTotal,
            myCanvActive,  myGraphActive, evActive, evActiveGraphN, statActive, now, willbe, myCanvActiveR0, evActiveGraphR0, myCanvActiveP, evActiveGraphP,
            myCanvDeaths,  myGraphDeaths,
            myCanvRatio01, myGraphRatio01,
            myCanvRatio02, myGraphRatio02]


def scanParameter(ntuple, N, minVal, maxVal, doSmearing, doGlobalFit):
    growthRate = 0.13
    carryingCapacity = 4e5

    tStart = ntuple[0]
    tStop  = ntuple[1]

    totalPopulation         = ntuple[2]
    symptomaticFraction     = ntuple[3]
    transmissionProbability = ntuple[4]
    recoveryRate            = ntuple[5]
    historyActive           = ntuple[6]

    xValues    = ntuple[7]
    yValues    = ntuple[8]
    erryValues = ntuple[9]

    seed(1)

    histo = TH1D('Histo', 'Histo', 100, 0, 1)
    histo.GetXaxis().SetTitle('g')
    histo.GetYaxis().SetTitle('Entries')

    scatter = TH2D('Scatter', 'Scatter par-chi2', 100, minVal, maxVal, 100, 0, 10)
    scatter.GetXaxis().SetTitle('Transmission probability')
    scatter.GetYaxis().SetTitle('#chi^{2}/d.o.f.')

    scatterPar = TH2D('ScatterPar', 'Scatter par-g', 100, minVal, maxVal, 100, 0, 1)
    scatterPar.GetXaxis().SetTitle('Transmission probability')
    scatterPar.GetYaxis().SetTitle('g')

    # 27 - 100
    scatterC0 = TH2D('ScatterC0', 'Scatter par-C0', 100, minVal, maxVal, 100, 4e5, 8e5)
    scatterCn = TH2D('ScatterCn', 'Scatter par-Cn', 100, minVal, maxVal, 100, 0.4e6, 1.6e6)
    # 15 - 27
#    scatterC0 = TH2D('ScatterC0', 'Scatter par-C0', 100, minVal, maxVal, 100, 0.8e5, 6e5)
#    scatterCn = TH2D('ScatterCn', 'Scatter par-Cn', 100, minVal, maxVal, 100, 1e5, 8e5)
    # 9 - 15
#    scatterC0 = TH2D('ScatterC0', 'Scatter par-C0', 100, minVal, maxVal, 100, 0.8e4, 1e5)
#    scatterCn = TH2D('ScatterCn', 'Scatter par-Cn', 100, minVal, maxVal, 100, 1e4, 3e5)
    # 0 - 9
#    scatterC0 = TH2D('ScatterC0', 'Scatter par-C0', 100, minVal, maxVal, 100, 1e3, 8e4)
#    scatterCn = TH2D('ScatterCn', 'Scatter par-Cn', 100, minVal, maxVal, 100, 0.2e4, 1e5)

    scatterC0.GetXaxis().SetTitle('Transmission probability')
    scatterC0.GetYaxis().SetTitle('C_{0}')

    scatterCn.GetXaxis().SetTitle('Transmission probability')
    scatterCn.GetYaxis().SetTitle('C_{n}')

    for i in range(N):
        print '\n==> Scan number:', i
        val = minVal + (random() * (maxVal - minVal))
        transmissionProbability = val

        if doGlobalFit == False:
            evolve = evolution([yValues[0], carryingCapacity, recoveryRate, growthRate], tStart, tStop, totalPopulation, symptomaticFraction, transmissionProbability, historyActive)
            status = evolve.runOptimization(xValues, yValues, erryValues, [2], doSmearing, -1)

            if status == 3:
                histo.     Fill(evolve.parValues[3])
                scatter.   Fill(val, evolve.chi2/evolve.dof)
                scatterPar.Fill(val, evolve.parValues[3])
                scatterC0. Fill(val, evolve.parValues[1])
                scatterCn. Fill(val, evolve.evolve(tStop, evolve.parValues)[3])
        else:
            timeList = ntuple[10]

            evolve = evolution([ntuple[11], ntuple[12], recoveryRate, ntuple[13]], tStart, timeList[0], totalPopulation, symptomaticFraction, transmissionProbability)
            evolutions = [evolution([0, 0, recoveryRate, ntuple[14+i]], timeList[i], timeList[i+1], totalPopulation, symptomaticFraction, transmissionProbability) for i in range(len(timeList)-1)]

            status, parValues, parNames = evolve.runGlobalOptimization(evolutions, xValues, yValues, erryValues, [2], doSmearing, -1)

            if status == 3:
                histo.     Fill(parValues[3])
                scatter.   Fill(val, evolve.chi2/evolve.dof)
                scatterPar.Fill(val, parValues[3])
                scatterC0. Fill(val, parValues[1])

    myCanv01 = TCanvas('myCanv01','Histogram')
    histo.Draw()
    myCanv01.Modified()
    myCanv01.Update()

    myCanv02 = TCanvas('myCanv02','Scatter par-chi2')
    scatter.Draw('gcolz')
    myCanv02.Modified()
    myCanv02.Update()

    myCanv03 = TCanvas('myCanv03','Scatter par-r')
    scatterPar.Draw('gcolz')
    myCanv03.Modified()
    myCanv03.Update()

    myCanv04 = TCanvas('myCanv04','C0')
    scatterC0.Draw('gcolz')
    myCanv04.Modified()
    myCanv04.Update()

    myCanv05 = TCanvas('myCanv05','Cn')
    scatterCn.Draw('gcolz')
    myCanv05.Modified()
    myCanv05.Update()

    return [myCanv01, histo,
            myCanv02, scatter,
            myCanv03, scatterPar,
            myCanv04, scatterC0,
            myCanv05, scatterCn]


def runModel(totalPopulation, symptomaticFraction, transmissionProbability, recoveryRate, doSmearing):
    myCanvModels = TCanvas('myCanvModels','Models')
    myCanvModels.Divide(1,2)
    myGraphTmp1 = TGraph()
    myGraphTmp2 = TGraph()

    timeList  = [9, 9+6, 9+6+11, 9+6+11+43, 9+6+11+43+60, 9+6+11+43+60+30, 9+6+11+43+60+30+60, 9+6+11+43+60+30+60 +120]
    parValues = [223, 8870, recoveryRate, 0.405, 0.318, 0.239, 0.161, 0.405, 0.161, 0.405, 0.161]

    myCanvModels.cd(1)
    myGraphTmp1.SetPoint(myGraphTmp1.GetN(), 0, 0)
    myGraphTmp1.SetPoint(myGraphTmp1.GetN(), timeList[-1], 0)
    myGraphTmp1.SetMarkerStyle(1)
    myGraphTmp1.Draw()
    myGraphTmp1.GetXaxis().SetTitle('Time (days)')
    myGraphTmp1.GetYaxis().SetTitle('Active cases')

    myCanvModels.cd(2)
    myGraphTmp2.SetPoint(myGraphTmp2.GetN(), 0, 0)
    myGraphTmp2.SetPoint(myGraphTmp2.GetN(), timeList[-1], 0)
    myGraphTmp2.SetMarkerStyle(1)
    myGraphTmp2.Draw()
    myGraphTmp2.GetXaxis().SetTitle('Time (days)')
    myGraphTmp2.GetYaxis().SetTitle('R')

    evolve = evolution(parValues[0:4], 0, timeList[0], totalPopulation, symptomaticFraction, transmissionProbability)
    evolutions = [
        evolution([0,  51800, recoveryRate, parValues[4]],  timeList[0], timeList[1], totalPopulation, symptomaticFraction, transmissionProbability),
        evolution([0, 211000, recoveryRate, parValues[5]],  timeList[1], timeList[2], totalPopulation, symptomaticFraction, transmissionProbability),
        evolution([0, 419000, recoveryRate, parValues[6]],  timeList[2], timeList[3], totalPopulation, symptomaticFraction, transmissionProbability),
        evolution([0,      0, recoveryRate, parValues[7]],  timeList[3], timeList[4], totalPopulation, symptomaticFraction, transmissionProbability),
        evolution([0,      0, recoveryRate, parValues[8]],  timeList[4], timeList[5], totalPopulation, symptomaticFraction, transmissionProbability),
        evolution([0,      0, recoveryRate, parValues[9]],  timeList[5], timeList[6], totalPopulation, symptomaticFraction, transmissionProbability),
        evolution([0,      0, recoveryRate, parValues[10]], timeList[6], timeList[7], totalPopulation, symptomaticFraction, transmissionProbability)]

    graphN = evolve.evolveGlobal(evolutions, evolutions[-1].tStop, parValues, True, True)

    if doSmearing == True:
        evolve.smearing()
    graphN = evolve.getGraphN()

    graphR0 = evolve.getGraphR0(graphN)
    graphN.SetLineColor(4)
    myCanvModels.cd(1)
    graphN.Draw('L same')
    myCanvModels.cd(2)
    graphR0.Draw('L same')
    """
    evolve1 = evolution([200, 8900, recoveryRate, 0.17], 0, 1000, totalPopulation, symptomaticFraction, transmissionProbability)
    evolve1.evolve(evolve1.tStop, evolve1.parValues, True)
    graph1 = evolve1.getGraphN()
    graph1.Draw('L same')

    evolve2 = evolution([200, 8900, 0.023, 0.17], 0, 1000, 60e6, 0.3, 4.7e-3)
    evolve2.evolve(evolve2.tStop, evolve2.parValues, True)
    graph2 = evolve2.getGraphN()
    graph2.SetLineColor(4)
    graph2.Draw('L same')

    evolve3 = evolution([200, 8900, 0.023, 0.17], 0, 1000, 60e6, 0.3, 4.7e-3)
    evolve3.evolve(evolve3.tStop, evolve3.parValues, True)
    graph3 = evolve3.getGraphN()
    graph3.SetLineColor(1)
    graph3.Draw('L same')
    """
    myCanvModels.SetGrid()
    myCanvModels.Modified()
    myCanvModels.Update()

    return [myCanvModels, myGraphTmp1, myGraphTmp2,
#            graph1, graph2, graph3]
            graphN, graphR0]


def runToyMC(evolve, nEv, nToy, doSmearing):
    evolve.evolve(evolve.tStop, evolve.parValues, True)

    nBins   = int(round(evolve.tStop - evolve.tStart,1))
    xValues = [evolve.tStart + i for i in range(nBins)]
    nSigma  = 6

    histo01 = TH1D('Histo01', evolve.parNames[0], 100, 10, 40)
    histo01.GetXaxis().SetTitle('Pulls ' + evolve.parNames[0])
    histo01.GetYaxis().SetTitle('Entries')

    histo02 = TH1D('Histo02', evolve.parNames[1], 100, -nSigma, nSigma)
    histo02.GetXaxis().SetTitle('Pulls ' + evolve.parNames[1])
    histo02.GetYaxis().SetTitle('Entries')

    histo03 = TH1D('Histo03', evolve.parNames[2], 100, -nSigma, nSigma)
    histo03.GetXaxis().SetTitle('Pulls ' + evolve.parNames[2])
    histo03.GetYaxis().SetTitle('Entries')

    histo04 = TH1D('Histo04', evolve.parNames[3], 100, -nSigma, nSigma)
    histo04.GetXaxis().SetTitle('Pulls ' + evolve.parNames[3])
    histo04.GetYaxis().SetTitle('Entries')

    for i in range(nToy):
        print '\n==> Toy number:', i
        yValues = [0. for i in range(nBins)]

        for j in range(nEv):
            x = evolve.funLookUp.GetRandom()

            it = 0
            while it < nBins and xValues[it] < x:
                it += 1
            yValues[it-1] += 1

        erryValues = [sqrt(yValues[j]) for j in range(nBins)]

        evToy = evolution(evolve.parValues[:],
                          evolve.tStart,
                          evolve.tStop,
                          evolve.totalPopulation,
                          evolve.symptomaticFraction,
                          evolve.transmissionProbability,
                          evolve.historyActiveDt)
        evToy.runOptimization(xValues, yValues, erryValues, [2], doSmearing, -1)

        histo01.Fill(((evToy.parValues[0] - evolve.parValues[0]) / evToy.fitErr[0]) if evToy.fitErr[0] != 0 else 0)
        histo02.Fill(((evToy.parValues[1] - evolve.parValues[1]) / evToy.fitErr[1]) if evToy.fitErr[1] != 0 else 0)
        histo03.Fill(((evToy.parValues[2] - evolve.parValues[2]) / evToy.fitErr[2]) if evToy.fitErr[2] != 0 else 0)
        histo04.Fill(((evToy.parValues[3] - evolve.parValues[3]) / evToy.fitErr[3]) if evToy.fitErr[3] != 0 else 0)

    canv01 = TCanvas('Canv01','Pulls ' + evolve.parNames[0])
    histo01.Draw()
    histo01.Fit('gaus')
    canv01.Modified()
    canv01.Update()

    canv02 = TCanvas('Canv02','Pulls ' + evolve.parNames[1])
    histo02.Draw()
    histo02.Fit('gaus')
    canv02.Modified()
    canv02.Update()

    canv03 = TCanvas('Canv03','Pulls ' + evolve.parNames[2])
    histo03.Draw()
    histo03.Fit('gaus')
    canv03.Modified()
    canv03.Update()

    canv04 = TCanvas('Canv04','Pulls ' + evolve.parNames[3])
    histo04.Draw()
    histo04.Fit('gaus')
    canv04.Modified()
    canv04.Update()

    return [canv01, histo01,
            canv02, histo02,
            canv03, histo03,
            canv04, histo04]


def runGlobalFit(country, active, totalPopulation, symptomaticFraction, transmissionProbability, recoveryRate, doSmearing):
    tStart = 0
    tStop  = 9+6+11+37 +20
    timeList = [9, 9+6, 9+6+11, tStop]
#    tStart = 9
#    tStop  = 9+6+11+37 +20
#    timeList = [9+6, 9+6+11, tStop]

    ntuple = [tStart, tStop, totalPopulation, symptomaticFraction, transmissionProbability, recoveryRate]

    xValues    = [i         for i in range(len(active.keys()))          if i >= tStart and i <= tStop]
    yValues    = [active[k] for i,k in enumerate(sorted(active.keys())) if i >= tStart and i <= tStop]
    erryValues = assignErrors(yValues)

    ntuple.extend([0, xValues, yValues, erryValues, timeList, 223, 8870, 0.405, 0.318, 0.239, 0.161])
#    ntuple.extend([0, xValues, yValues, erryValues, timeList, 2640, 51800, 0.318, 0.239, 0.161])

    evActive = evolution([ntuple[11], ntuple[12], recoveryRate, ntuple[13]], tStart, timeList[0], totalPopulation, symptomaticFraction, transmissionProbability)
    evolutions = [evolution([0, 0, recoveryRate, ntuple[14+i]], timeList[i], timeList[i+1], totalPopulation, symptomaticFraction, transmissionProbability) for i in range(len(timeList)-1)]
    istat, parValues, parNames = evActive.runGlobalOptimization(evolutions, xValues, yValues, erryValues, [2], doSmearing)

    myCanvActive = TCanvas('myCanvActive_' + country, 'Active cases ' + country)
    myGraphActive = TGraphErrors()
    myGraphActive.SetMarkerStyle(20)

    for i in range(len(xValues)):
        myGraphActive.SetPoint(myGraphActive.GetN(), xValues[i], yValues[i])
        myGraphActive.SetPointError(myGraphActive.GetN()-1, 0, erryValues[i])
    myGraphActive.Draw('APE1')
    myGraphActive.GetHistogram().GetXaxis().SetTitle('Time (days)')
    myGraphActive.GetHistogram().GetYaxis().SetTitle('Active cases affected by CoViD-19')

    evActive.evolveGlobal(evolutions, evolutions[-1].tStop, parValues, True)
    if doSmearing == True:
        evActive.smearing()
    evActiveGraphN = evActive.getGraphN()
    evActiveGraphN.Draw('PL same')
    statActive = evActive.addStats(parNames, parValues)

    print '==> Active cases, history active cases * dt, p-infected, Carrying capacity', evActive.evolveGlobal(evolutions, tStop, parValues), '@', tStop, 'day'

    now = TLine(len(active)-1, 0, len(active)-1, evActive.fitFun.Eval(len(active)-1))
    now.SetLineColor(4)
    now.SetLineWidth(2)
    now.Draw('same')

    willbe = TLine(evActive.fitFun.GetMaximumX(), 0, evActive.fitFun.GetMaximumX(), evActive.fitFun.GetMaximum())
    willbe.SetLineColor(6)
    willbe.SetLineWidth(2)
    willbe.Draw('same')

    myCanvActive.SetGrid()
    myCanvActive.Modified()
    myCanvActive.Update()

    myCanvActiveR0 = TCanvas('myCanvActiveR0_' + country, 'R0 ' + country)

    evActiveGraphR0 = evActive.getGraphR0(evActiveGraphN)
    evActiveGraphR0.Draw('APL')
    evActiveGraphR0.GetHistogram().GetXaxis().SetTitle('Time (days)')
    evActiveGraphR0.GetHistogram().GetYaxis().SetTitle('R')

    myCanvActiveR0.SetGrid()
    myCanvActiveR0.Modified()
    myCanvActiveR0.Update()

    myCanvActiveP = TCanvas('myCanvActiveP_' + country, 'Probability infected ' + country)

    evActiveGraphP = evActive.getGraphGlobalPinfect(evolutions, parValues)
    evActiveGraphP.Draw('APL')
    evActiveGraphP.GetHistogram().GetXaxis().SetTitle('Time (days)')
    evActiveGraphP.GetHistogram().GetYaxis().SetTitle('Probability of being infected')

    myCanvActiveP.SetGrid()
    myCanvActiveP.Modified()
    myCanvActiveP.Update()

    return [ntuple,
            myCanvActive,  myGraphActive, evActive, evActiveGraphN, statActive, now, willbe, myCanvActiveR0, evActiveGraphR0, myCanvActiveP, evActiveGraphP]


######################
# Start main program #
######################
SetStyle()

graphModel = runModel(60e6, 0.3, 4.7e-3, 0.023, False)

##################################
# Read data from database: Italy #
##################################
print '=== Downloading data ==='
url      = 'https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv'
fileName = url.split('/')[-1]
saveDataFromURL(url)
print '=== Done ===\n'
active    = readDataFromFile(fileName,  6)
recovered = readDataFromFile(fileName,  9)
deaths    = readDataFromFile(fileName, 10)
total     = readDataFromFile(fileName, 11)
#graphItaly = analyzeData('Italy', total, active, recovered, deaths, 27, 100, 60e6, 0.3, 4.7e-3, 0.023, False)

#graphItaly[3].cd()
#graphModel[3].Draw('same')
#graphScan = scanParameter(graphItaly[0], 100, 0.001, 0.02, False, False)
#graphToy = runToyMC(graphItaly[5], 3340564, 400, False)

#graphGlobalFit = runGlobalFit('Italy', active, 60e6, 0.3, 5.5e-3, 0.023, True)
#graphGlobalScan = scanParameter(graphGlobalFit[0], 400, 0.001, 0.02, False, True)

"""
##################################
# Read data from database: World #
##################################
country  = 'China'
province = 'Hubei'

print '\n=== Downloading data ==='
url           = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv'
fileNameTotal = url.split('/')[-1]
saveDataFromURL(url)

url            = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv'
fileNameDeaths = url.split('/')[-1]
saveDataFromURL(url)

url               = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv'
fileNameRecovered = url.split('/')[-1]
saveDataFromURL(url)
print '=== Done ===\n'

total     = readDataFromFile(fileNameTotal,     -1, country, province)
deaths    = readDataFromFile(fileNameDeaths,    -1, country, province)
recovered = readDataFromFile(fileNameRecovered, -1, country, province)
active    = {}
for k in sorted(total.keys()):
    active[k] = total[k] - deaths[k] - recovered[k]
graphWorld = analyzeData(country, total, active, recovered, deaths, 22, 100, 60e6, 0.3, 4.8e-3, 0.031, False)
"""

raw_input('\nPress <ret> to end -> ')

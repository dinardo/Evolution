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
        data = csv.reader(csvFile, delimiter=',')
        isFirstRow = True
        for r in data:
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


def analyzeItaly(tStart, tStop, totalPopulation, symptomaticFraction, transmissionProbability, recoveryRate):
    ntuple = [tStart, tStop, totalPopulation, symptomaticFraction, transmissionProbability, recoveryRate]
    scaleError = 3
    farFromMax = 0.95
    growthRate = 0.13
    carryingCapacity = 4e5


    ###########################
    # Read data from database #
    ###########################
    print '=== Downloading data ==='
    url      = 'https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv'
    fileName = url.split('/')[-1]
    saveDataFromURL(url)
    print '=== Done ===\n'

    active    = readDataFromFile(fileName,  6)
    recovered = readDataFromFile(fileName,  9)
    deaths    = readDataFromFile(fileName, 10)
    total     = readDataFromFile(fileName, 11)


    ################
    # Active cases #
    ################
    myCanvActive = TCanvas('myCanvActive','Active cases')
    myGraphActive = TGraphErrors()
    myGraphActive.SetMarkerStyle(20)

    for k in sorted(active.keys()):
        myGraphActive.SetPoint(myGraphActive.GetN(), myGraphActive.GetN(), active[k])
        myGraphActive.SetPointError(myGraphActive.GetN()-1, 0, sqrt(active[k]) * scaleError)
    myGraphActive.SetPoint(myGraphActive.GetN(), tStop, 0)
    myGraphActive.SetPointError(myGraphActive.GetN()-1, 0, 0)
    myGraphActive.Draw('APE1')
    myGraphActive.GetHistogram().GetXaxis().SetTitle('Time (days)')
    myGraphActive.GetHistogram().GetYaxis().SetTitle('Active cases affected by CoViD-19')

    xValues    = [i                            for i in range(len(active.keys()))          if i >= tStart and i <= tStop]
    yValues    = [active[k]                    for i,k in enumerate(sorted(active.keys())) if i >= tStart and i <= tStop]
    erryValues = [sqrt(active[k]) * scaleError for i,k in enumerate(sorted(active.keys())) if i >= tStart and i <= tStop]

    historyActive = 0.
    for i,k in enumerate(sorted(active.keys())):
        if i < tStart:
            historyActive += active[k]
    ntuple.extend([historyActive, xValues, yValues, erryValues])
    print '==> History active cases:', historyActive

    evActive = evolution([yValues[0], growthRate, recoveryRate, carryingCapacity], tStart, tStop, totalPopulation, symptomaticFraction, transmissionProbability, historyActive)
    evActive.runOptimization(xValues, yValues, erryValues, [2])
    evActiveGraphN = evActive.getGraphN()
    evActiveGraphN.Draw('PL same')
    statActive = evActive.addStats()

    print '==> Active cases, history active cases * dt, p-infected, Carrying capacity', evActive.evolveActive(tStop, evActive.parValues), '@', tStop, 'day'

    now = TLine(len(active)-1, 0, len(active)-1, evActive.fitFun.GetMaximum())
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

    myCanvActiveR0 = TCanvas('myCanvActiveR0','R0')

    evActiveGraphR0 = evActive.getGraphR0(evActiveGraphN)
    evActiveGraphR0.Draw('APL')
    evActiveGraphR0.GetHistogram().GetXaxis().SetTitle('Time (days)')
    evActiveGraphR0.GetHistogram().GetYaxis().SetTitle('R_{0}')

    myCanvActiveR0.SetGrid()
    myCanvActiveR0.Modified()
    myCanvActiveR0.Update()

    myCanvActiveP = TCanvas('myCanvActiveP','Probability infected')

    evActiveGraphP = evActive.getGraphPinfect()
    evActiveGraphP.Draw('APL')
    evActiveGraphP.GetHistogram().GetXaxis().SetTitle('Time (days)')
    evActiveGraphP.GetHistogram().GetYaxis().SetTitle('Probability')

    myCanvActiveP.SetGrid()
    myCanvActiveP.Modified()
    myCanvActiveP.Update()


    #################
    # CONTROL PLOTS #
    #################

    ###############
    # Total cases #
    ###############
    myCanvTotal = TCanvas('myCanvTotal', 'Total cases')
    myGraphTotal = TGraphErrors()
    myGraphTotal.SetMarkerStyle(20)

    for k in sorted(total.keys()):
        myGraphTotal.SetPoint(myGraphTotal.GetN(), myGraphTotal.GetN(), total[k])
        myGraphTotal.SetPointError(myGraphTotal.GetN()-1, 0, sqrt(total[k]) * scaleError)

    myGraphTotal.Draw('APE1')
    myGraphTotal.GetHistogram().GetXaxis().SetTitle('Time (days)')
    myGraphTotal.GetHistogram().GetYaxis().SetTitle('Total cases affected by CoViD-19')

    myCanvTotal.SetGrid()
    myCanvTotal.Modified()
    myCanvTotal.Update()


    ##########
    # Deaths #
    ##########
    myCanvDeaths = TCanvas('myCanvDeaths','Deaths')
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
    myCanvRatio01 = TCanvas('myCanvRatio01','Ratio')
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
    myCanvRatio02 = TCanvas('myCanvRatio02','Ratio')
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

    print '==> Doubling time:', round(log(2)/evActive.parValues[1],1), 'days'

    return [ntuple,
            myCanvTotal,   myGraphTotal,
            myCanvActive,  myGraphActive, evActive, evActiveGraphN, statActive, now, willbe, myCanvActiveR0, evActiveGraphR0, myCanvActiveP, evActiveGraphP,
            myCanvDeaths,  myGraphDeaths,
            myCanvRatio01, myGraphRatio01,
            myCanvRatio02, myGraphRatio02]


def scanParameter(ntuple, N, minVal, maxVal):
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
    histo.GetXaxis().SetTitle('r')
    histo.GetYaxis().SetTitle('Entries')

    scatter = TH2D('Scatter', 'Scatter par-chi2', 100, minVal, maxVal, 100, 0, 10)
    scatter.GetXaxis().SetTitle('Transmission probability')
    scatter.GetYaxis().SetTitle('#chi^{2}/d.o.f.')

    scatterPar = TH2D('ScatterPar', 'Scatter par-r', 100, minVal, maxVal, 100, 0, 1)
    scatterPar.GetXaxis().SetTitle('Transmission probability')
    scatterPar.GetYaxis().SetTitle('r')

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
        val = minVal + (random() * (maxVal - minVal))
        transmissionProbability = val

        evolve = evolution([yValues[0], growthRate, recoveryRate, carryingCapacity], tStart, tStop, totalPopulation, symptomaticFraction, transmissionProbability, historyActive)
        status = evolve.runOptimization(xValues, yValues, erryValues, [2], -1)

        if status == 3:
            histo.     Fill(evolve.parValues[1])
            scatter.   Fill(val, evolve.chi2/evolve.dof)
            scatterPar.Fill(val, evolve.parValues[1])
            scatterC0. Fill(val, evolve.parValues[3])
            scatterCn. Fill(val, evolve.evolveActive(tStop, evolve.parValues)[3])

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


def analyzeWorld(country, province, tStart, tStop, totalPopulation, symptomaticFraction, transmissionProbability, recoveryRate):
    scaleError = 3
    growthRate = 0.13
    carryingCapacity = 4e5


    ###########################
    # Read data from database #
    ###########################
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


    ###############
    # Total cases #
    ###############
    myCanv01 = TCanvas('myCanv01','Total cases ' + country)
    myGraph01 = TGraphErrors()
    myGraph01.SetMarkerStyle(20)

    for k in sorted(total.keys()):
        myGraph01.SetPoint(myGraph01.GetN(), myGraph01.GetN(), total[k])
        myGraph01.SetPointError(myGraph01.GetN()-1, 0, sqrt(total[k]) * scaleError)

    myGraph01.Draw('APE1')
    myGraph01.GetHistogram().GetXaxis().SetTitle('Time (days)')
    myGraph01.GetHistogram().GetYaxis().SetTitle('Total cases')

    myCanv01.SetGrid()
    myCanv01.Modified()
    myCanv01.Update()


    ################
    # Active cases #
    ################
    myCanv02 = TCanvas('myCanv02','Active cases ' + country)
    myGraph02 = TGraphErrors()
    myGraph02.SetMarkerStyle(20)

    for k in sorted(active.keys()):
        myGraph02.SetPoint(myGraph02.GetN(), myGraph02.GetN(), active[k])
        myGraph02.SetPointError(myGraph02.GetN()-1, 0, sqrt(active[k]) * scaleError)
    myGraph02.Draw('APE1')
    myGraph02.GetHistogram().GetXaxis().SetTitle('Time (days)')
    myGraph02.GetHistogram().GetYaxis().SetTitle('Active cases')

    xValues    = [i                            for i in range(len(active.keys()))          if i >= tStart and i <= tStop]
    yValues    = [active[k]                    for i,k in enumerate(sorted(active.keys())) if i >= tStart and i <= tStop]
    erryValues = [sqrt(active[k]) * scaleError for i,k in enumerate(sorted(active.keys())) if i >= tStart and i <= tStop]

    historyActive = 0.
    for i,k in enumerate(sorted(active.keys())):
        if i < tStart:
            historyActive += active[k]
    print '==> History active cases:', historyActive

    evActive02 = evolution([yValues[0], growthRate, recoveryRate, carryingCapacity], tStart, tStop, totalPopulation, symptomaticFraction, transmissionProbability, historyActive)
    evActive02.runOptimization(xValues, yValues, erryValues, [2])
    evActiveGraph02N = evActive02.getGraphN()
    evActiveGraph02N.Draw('PL same')
    stat02 = evActive02.addStats()

    print '==> Active cases, history active cases * dt, p-infected, Carrying capacity', evActive02.evolveActive(tStop, evActive02.parValues), '@', tStop, 'day'

    myCanv02.SetGrid()
    myCanv02.Modified()
    myCanv02.Update()

    myCanv02R0 = TCanvas('myCanv02R0','R0 ' + country)

    evActiveGraph02R0 = evActive02.getGraphR0(evActiveGraph02N)
    evActiveGraph02R0.Draw('APL')
    evActiveGraph02R0.GetHistogram().GetXaxis().SetTitle('Time (days)')
    evActiveGraph02R0.GetHistogram().GetYaxis().SetTitle('R_{0}')

    myCanv02R0.SetGrid()
    myCanv02R0.Modified()
    myCanv02R0.Update()

    myCanv02P = TCanvas('myCanv02P','Probability infected ' + country)

    evActiveGraph02P = evActive02.getGraphPinfect()
    evActiveGraph02P.Draw('APL')
    evActiveGraph02P.GetHistogram().GetXaxis().SetTitle('Time (days)')
    evActiveGraph02P.GetHistogram().GetYaxis().SetTitle('Probability')

    myCanv02P.SetGrid()
    myCanv02P.Modified()
    myCanv02P.Update()


    ##########
    # Deaths #
    ##########
    myCanv03 = TCanvas('myCanv03','Deaths ' + country)
    myGraph03 = TGraphErrors()
    myGraph03.SetMarkerStyle(20)

    for k in sorted(deaths.keys()):
        myGraph03.SetPoint(myGraph03.GetN(), myGraph03.GetN(), deaths[k])
        myGraph03.SetPointError(myGraph03.GetN()-1, 0, sqrt(deaths[k]))

    myGraph03.Draw('APE1')
    myGraph03.GetHistogram().GetXaxis().SetTitle('Time (days)')
    myGraph03.GetHistogram().GetYaxis().SetTitle('Total deaths')

    myCanv03.SetGrid()
    myCanv03.Modified()
    myCanv03.Update()


    ########################
    # Ratio deaths / total #
    ########################
    myCanv04 = TCanvas('myCanv04','Ratio ' + country)
    myGraph04 = TGraphErrors()
    myGraph04.SetMarkerStyle(20)

    for k in sorted(deaths.keys()):
        myGraph04.SetPoint(myGraph04.GetN(), myGraph04.GetN(), (deaths[k]/total[k]) if total[k] != 0 else 0)
        myGraph04.SetPointError(myGraph04.GetN()-1, 0, (myGraph04.GetY()[myGraph04.GetN()-1] * sqrt(deaths[k]/(deaths[k]*deaths[k]) + total[k]/(total[k]*total[k]))) if total[k] != 0 and deaths[k] != 0 else 0 )

    myGraph04.Draw('APE1')
    myGraph04.GetHistogram().GetXaxis().SetTitle('Time (days)')
    myGraph04.GetHistogram().GetYaxis().SetTitle('Total deaths / Total cases')

    myCanv04.SetGrid()
    myCanv04.Modified()
    myCanv04.Update()


    ############################################
    # Ratio delta(deaths + recovered) / active #
    ############################################
    myCanv05 = TCanvas('myCanv05','Ratio ' + country)
    myGraph05 = TGraphErrors()
    myGraph05.SetMarkerStyle(20)

    sortedKeys = sorted(deaths.keys())
    for i,k in enumerate(sortedKeys[1:]):
        myGraph05.SetPoint(myGraph05.GetN(), myGraph05.GetN()+1, ((deaths[k] - deaths[sortedKeys[i]] + recovered[k] - recovered[sortedKeys[i]]) / active[k]) if active[k] != 0 else 0)
        numerator = deaths[k] + deaths[sortedKeys[i]] + recovered[k] + recovered[sortedKeys[i]]
        myGraph05.SetPointError(myGraph05.GetN()-1, 0, (myGraph05.GetY()[myGraph05.GetN()-1] * sqrt(numerator/pow(numerator,2) + active[k]/(active[k]*active[k]))) if active[k] != 0 and numerator != 0 else 0)

    myGraph05.Draw('AP')
    myGraph05.GetHistogram().GetXaxis().SetTitle('Time (days)')
    myGraph05.GetHistogram().GetYaxis().SetTitle('#Delta Recovered (alive + dead) / Active cases')

    myCanv05.SetGrid()
    myCanv05.Modified()
    myCanv05.Update()

    return [myCanv01, myGraph01,
            myCanv02, myGraph02, evActiveGraph02N, stat02,  myCanv02R0, evActiveGraph02R0, myCanv02P, evActiveGraph02P,
            myCanv03, myGraph03,
            myCanv04, myGraph04,
            myCanv05, myGraph05]


def runModel(totalPopulation, symptomaticFraction, transmissionProbability, recoveryRate):
    myCanvModels = TCanvas('myCanvModels','Models')
    myCanvModels.Divide(1,2)
    myGraphTmp1 = TGraph()
    myGraphTmp2 = TGraph()

    timeList = [9, 9+6, 9+6+11, 9+6+11+37, 9+6+11+37+60, 9+6+11+37+60+30, 9+6+11+37+60+30+60, 9+6+11+37+60+30+60+120]

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
    myGraphTmp2.GetYaxis().SetTitle('R_{0}')

    parList = [[ 2640, 0.318, recoveryRate,  51800],
               [ 8890, 0.239, recoveryRate, 211000],
               [48600, 0.172, recoveryRate, 397000],
               [0,     0.405, recoveryRate, 0],
               [0,     0.172, recoveryRate, 0],
               [0,     0.405, recoveryRate, 0],
               [0,     0.172, recoveryRate, 0]]

    evolve = evolution([223, 0.405, recoveryRate, 8870], 0, timeList[0], totalPopulation, symptomaticFraction, transmissionProbability)
    graphN = evolve.combineEvolutions(parList, timeList, totalPopulation, symptomaticFraction, transmissionProbability)
    # graphN = evolve.smearing(graphN)
    graphR0 = evolve.getGraphR0(graphN)
    graphN.SetLineColor(4)
    myCanvModels.cd(1)
    graphN.Draw('L same')
    myCanvModels.cd(2)
    graphR0.Draw('L same')
    """
    evolve1 = evolution([200, 0.17, recoveryRate, 8900], 0, 1000, totalPopulation, symptomaticFraction, transmissionProbability)
    graph1 = evolve1.getGraphN()
    graph1.Draw('L same')

    evolve2 = evolution([200, 0.17, 0.023, 8900], 0, 1000, 60e6, 0.3, 4.7e-3)
    graph2 = evolve2.getGraphN()
    graph2.SetLineColor(4)
    graph2.Draw('L same')

    evolve3 = evolution([200, 0.17, 0.023, 8900], 0, 1000, 60e6, 0.3, 4.7e-3)
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


def runToyMC(evolve, nEv, nToy):
    evolve.evolveActive(evolve.tStop, evolve.parValues, True)
    evolve.generateFunctionLookUpTable()

    nBins = int(round((evolve.tStop - evolve.tStart)/evolve.dt,1))

    xValues    = [evolve.tStart + i * evolve.dt for i in range(nBins)]
    erryValues = [0. for i in range(nBins)]

    histo01 = TH1D('Histo01', evolve.parNames[0], 100, -3, 3)
    histo01.GetXaxis().SetTitle('Pulls ' + evolve.parNames[0])
    histo01.GetYaxis().SetTitle('Entries')

    for i in range(nToy):
        print '\n==> Toy number:', i
        yValues = [0. for i in range(nBins)]

        for j in range(nEv):
            x = evolve.funLookUp.GetRandom()

            it = 0
            while it < nBins and xValues[it] < x:
                it += 1
            yValues[it-1] += 1

        for j in range(nBins):
            erryValues[j] = sqrt(yValues[j])

        evToy = evolution(evolve.parValues,
                          evolve.tStart,
                          evolve.tStop,
                          evolve.totalPopulation,
                          evolve.symptomaticFraction,
                          evolve.transmissionProbability,
                          evolve.historyActiveDt)
        evToy.runOptimization(xValues, yValues, erryValues, [2])

        histo01.Fill((evToy.parValues[0] - evolve.parValues[0]) / evToy.fitErr[0])

    canv01 = TCanvas('Canv01','Pulls ' + evolve.parNames[0])
    histo01.Draw()
    histo01.Fit('gaus')
    canv01.Modified()
    canv01.Update()

    return [canv01, histo01]


######################
# Start main program #
######################
SetStyle()

#graphModel = runModel(60e6, 0.3, 4.7e-3, 0.023)
graphItaly = analyzeItaly(27, 100, 60e6, 0.3, 4.7e-3, 0.023)
#graphItaly[3].cd()
#graphModel[3].Draw('same')
#scan = scanParameter(graphItaly[0], 100, 0.001, 0.02)
#graphWorld = analyzeWorld('China', 'Hubei', 22, 100, 60e6, 0.3, 4.7e-3, 0.031)

graphItaly[5].tStart = 0
runToyMC(graphItaly[5], 700000, 100)

raw_input('\nPress <ret> to end -> ')

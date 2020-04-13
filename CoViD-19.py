######################################
# Program to analyse CoViD-19 spread #
# Author: Mauro E. Dinardo           #
######################################


import csv
from datetime import datetime
from math     import sqrt, log, exp

from pyWget    import saveDataFromURL
from evolution import evolution

from ROOT import gROOT, gStyle, gApplication, TGaxis, TCanvas, TGraphErrors, TGraph, TLine


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




######################
# Start main program #
######################
SetStyle()


#############
# Constants #
#############
scaleError = 3
farFromMax = 0.95

deathFraction           = 0.11
totalPopulation         = 60e6
symptomaticFraction     = 0.3
transmissionProbability = 0.04

tStart   =   29 # [Day]
tStop    =  100 # [Day]
country  = 'China'
province = 'Hubei'


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

"""
##########
# Models #
##########
myCanvModels = TCanvas('myCanvModels','Models')
myCanvModels.Divide(1,2)
myGraphTmp1 = TGraph()
myGraphTmp2 = TGraph()

t0 = 8
timeList = [t0, t0+10, t0+10+9, t0+10+9+60, t0+10+9+60+20, t0+10+9+60+20+30, t0+10+9+60+20+30+20, 1000]
#timeList = [t0, t0+10, t0+10+9, t0+10+9+20, t0+10+9+20+20, t0+10+9+20+20+30, t0+10+9+20+20+30+20, 1000]

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

parList = [[000, 0.297, 0.023, 0.463],
           [000, 0.214, 0.023, 0.176],
           [000, 0.123, 0.023, 0.0078],
           [000, 0.452, 0.023, 0.589],
           [000, 0.123, 0.023, 0.0078],
           [000, 0.452, 0.023, 0.589],
           [000, 0.123, 0.023, 0.0078]]

evolve = evolution([200, 0.452, 0.023, 0.589], 0, timeList[0], deathFraction, totalPopulation, symptomaticFraction, transmissionProbability)
graphN, graphR0 = evolve.combineEvolutions(parList, timeList, deathFraction, totalPopulation, symptomaticFraction, transmissionProbability)
#graphN  = evolve.smearing(graphN)
#graphR0 = evolve.smearing(graphR0)
graphN.SetLineColor(4)
myCanvModels.cd(1)
graphN.Draw('L same')
myCanvModels.cd(2)
graphR0.Draw('L same')

#evolve1 = evolution([100, 0.123, 0.023, 0.0078], 0, 1000, 0.11, 60e6, 0.3, 0.04)
#graph1 = evolve1.getGraphN()
#graph1.Draw('L same')
#
#evolve2 = evolution([100, 0.123, 0.023, 0.0078], 0, 1000, 0.11, 60e6, 0.3, 0.04)
#graph2 = evolve2.getGraphN()
#graph2.SetLineColor(4)
#graph2.Draw('L same')
#
#evolve3 = evolution([100, 0.123, 0.023, 0.0078], 0, 1000, 0.11, 60e6, 0.3, 0.04)
#graph3 = evolve3.getGraphN()
#graph3.SetLineColor(1)
#graph3.Draw('L same')

myCanvModels.SetGrid()
myCanvModels.Modified()
myCanvModels.Update()
"""

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

#graphN.Draw('L same')

xValues    = [i                            for i in range(len(active.keys()))          if i >= tStart and i <= tStop]
yValues    = [active[k]                    for i,k in enumerate(sorted(active.keys())) if i >= tStart and i <= tStop]
erryValues = [sqrt(active[k]) * scaleError for i,k in enumerate(sorted(active.keys())) if i >= tStart and i <= tStop]
historyActive = 0.
for i,k in enumerate(sorted(active.keys())):
    if i < tStart:
        historyActive += active[k]
evActive = evolution([active[sorted(active.keys())[int(tStart)]], 0.13, 0.023, 5e5], tStart, tStop, deathFraction, totalPopulation, symptomaticFraction, transmissionProbability, 0., historyActive)
evActive.runOptimization(xValues, yValues, erryValues, [2])
evActiveGraph = evActive.getGraphN()
evActiveGraph.Draw('PL same')
statActive = evActive.addStats()

nowA = TLine(len(active)-1, 0, len(active)-1, evActive.fitFun.GetMaximum())
nowA.SetLineColor(4)
nowA.SetLineWidth(2)
nowA.Draw('same')

willbeA = TLine(evActive.fitFun.GetMaximumX(), 0, evActive.fitFun.GetMaximumX(), evActive.fitFun.GetMaximum())
willbeA.SetLineColor(6)
willbeA.SetLineWidth(2)
willbeA.Draw('same')

myCanvActive.SetGrid()
myCanvActive.Modified()
myCanvActive.Update()

myCanvActiveR0 = TCanvas('myCanvActiveR0','R0')

evActiveGraphR0 = evActive.getGraphR0()
evActiveGraphR0.Draw('APL')
evActiveGraphR0.GetHistogram().GetXaxis().SetTitle('Time (days)')
evActiveGraphR0.GetHistogram().GetYaxis().SetTitle('R_{0}')

myCanvActiveR0.SetGrid()
myCanvActiveR0.Modified()
myCanvActiveR0.Update()


#################
# CONTROL PLOTS #
#################

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


print 'Rth:', round(evActive.parValues[1] / evActive.parValues[2],1)
print 'Doubling time:', round(log(2)/evActive.parValues[1],1), 'days'


###################
# OTHER COUNTRIES #
###################

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
tStart =  22
tStop  = 100


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
evActive02 = evolution([active[sorted(active.keys())[int(tStart)]], 0.3, 0.031, 0.03], tStart, tStop, 0.04, totalPopulation, symptomaticFraction, transmissionProbability, 0., historyActive)
evActive02.runOptimization(xValues, yValues, erryValues, [2])

evActiveGraph02 = evActive02.getGraphN()
evActiveGraph02.Draw('PL same')
stat02 = evActive02.addStats()

myCanv02.SetGrid()
myCanv02.Modified()
myCanv02.Update()

myCanv02R0 = TCanvas('myCanv02R0','R0 ' + country)

evActiveGraph02R0 = evActive02.getGraphR0()
evActiveGraph02R0.Draw('APL')
evActiveGraph02R0.GetHistogram().GetXaxis().SetTitle('Time (days)')
evActiveGraph02R0.GetHistogram().GetYaxis().SetTitle('R_{0}')

myCanv02R0.SetGrid()
myCanv02R0.Modified()
myCanv02R0.Update()


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


raw_input('\nPress <ret> to end -> ')

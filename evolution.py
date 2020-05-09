#################################################
# Finite difference model for disease evolution #
# Author: Mauro E. Dinardo                      #
#################################################

from array   import array
from math    import sqrt, log, exp, factorial, pi
from decimal import *
import numpy as np

from ROOT import TMinuit, Long, Double, TString, TPaveText, TGraph, TF1

class evolution(object):
    def __init__(self, parValues, tStart, tStop, totalPopulation,
                 symptomaticFraction     = 0.3,
                 transmissionProbability = 4.7e-3,
                 historyActiveDt         = 0.):

        self.parNames  = ['Initial population', 'Carrying capacity', 'Recovery rate', 'Growth rate']
        self.parValues = parValues

        self.tStart = tStart
        self.tStop  = tStop
        self.doSmearing = False

        ####################
        # Model parameters #
        ####################
        self.totalPopulation         = totalPopulation
        self.symptomaticFraction     = symptomaticFraction
        self.transmissionProbability = transmissionProbability
        self.historyActiveDt         = historyActiveDt

        self.dt = 0.02 # [Day]

        self.nMeas, self.chi2, self.dof = 0., 0., 0.

        self.lookUpTable = np.zeros(0)
        self.bins = []

        self.funLookUp = TF1('EvolutionLookUp', self.evolveLookUpWrapper, self.tStart, self.tStop, len(self.parNames))
        self.setFitFun(self.evolveWrapper, self.tStart, self.tStop, self.parNames)

    def setFitFun(self, fun, tStart, tStop, parNames):
        self.fitFun = TF1('Evolution', fun, tStart, tStop, len(parNames))

    def evolveLookUpWrapper(self, t, par):
        return self.evolveLookUp(t[0])

    def evolveLookUp(self, t):
        return self.lookUpTable[np.digitize([t], self.bins)[0] - 1]

    def evolveWrapper(self, t, par):
        return self.evolve(t[0], par)[0]

    def evolve(self, t, par, doLookUp = False):
        T   = t - self.tStart
        N   = par[0] / self.symptomaticFraction
        CC  = par[1]
        CCn = 0.
        totalInfected = 0.
        historyActiveDt = self.historyActiveDt / self.symptomaticFraction
        if doLookUp == True:
            self.lookUpTable = np.zeros(int(round(T/self.dt,1)) + 1)
            self.bins = [self.tStart + i * self.dt for i in range(len(self.lookUpTable))]

        for n in range(int(round(T/self.dt,1))):
            if doLookUp == True:
                self.lookUpTable[n] = N * self.symptomaticFraction

            totalInfected = N + historyActiveDt * par[2]
            g = par[3] * (1. - totalInfected / CC)

            Nn = N
            N += self.dt * N * (g - par[2])

            CCn = CC
            CC += self.dt * par[3] / self.transmissionProbability * (N - Nn * (1. - self.dt * par[2])) * (1. - CC / self.totalPopulation)

            historyActiveDt += Nn * self.dt

        if doLookUp == True:
            self.lookUpTable[int(round(T/self.dt,1))] = N * self.symptomaticFraction

        return [N * self.symptomaticFraction, historyActiveDt * self.symptomaticFraction, (totalInfected / CCn) if CCn != 0. else 0., CC]

    def evolveGlobalWrapper(self, t, par):
        return self.evolveGlobal(self.evolutions, t[0], par)[0]

    def evolveGlobal(self, evolutions, t, par, doLookUp = False, doDefaults = False):
        self.parValues[0] = par[0]
        self.parValues[1] = par[1]
        self.parValues[2] = par[2]
        self.parValues[3] = par[3]
        T = 0

        if t <= self.tStop:
            T = t
        else:
            T = self.tStop

        [active, historyActiveDt, Pinf, CC] = self.evolve(T, self.parValues, doLookUp)

        if t <= self.tStop:
            return [active, historyActiveDt, Pinf, CC]

        for i,ev in enumerate(evolutions):
            ev.parValues[0] = ev.parValues[0] if doDefaults == True and ev.parValues[0] != 0. else active
            ev.parValues[1] = ev.parValues[1] if doDefaults == True and ev.parValues[1] != 0. else CC
            ev.parValues[2] = par[2]
            ev.parValues[3] = par[4+i]
            ev.historyActiveDt = historyActiveDt

            if t <= ev.tStop:
                T = t
            else:
                T = ev.tStop

            [active, historyActiveDt, Pinf, CC] = ev.evolve(T, ev.parValues, doLookUp)
            if doLookUp == True:
                self.lookUpTable = np.concatenate([self.lookUpTable[:-1], ev.lookUpTable])
                self.bins = self.bins[:-1] + ev.bins

            if t <= ev.tStop:
                return [active, historyActiveDt, Pinf, CC]

    def totalInfected(self, up2Time):
        return self.evolve(up2Time, self.parValues)[0] + self.totalRecovered(up2Time)

    def totalRecovered(self, up2Time):
        return (self.historyActiveDt + self.evolve(up2Time, self.parValues)[1]) * self.parValues[2]

    def costFunction(self, npar, grad, fval, par, iflag):
        self.nMeas, self.chi2, delta = 0., 0., 0.

        self.evolve(self.xValues[len(self.xValues)-1], par, True)
        if self.doSmearing == True:
            self.smearing()

        for i,erry in enumerate(self.erryValues):
            if erry != 0:
                delta = (self.yValues[i] - self.evolveLookUp(self.xValues[i])) / erry
                self.chi2 += delta * delta
                self.nMeas += 1

        fval[0] = self.chi2

    def costFunctionGlobal(self, npar, grad, fval, par, iflag):
        self.nMeas, self.chi2, delta = 0., 0., 0.

        self.evolveGlobal(self.evolutions, self.xValues[len(self.xValues)-1], par, True)
        if self.doSmearing == True:
            self.smearing()

        for i,erry in enumerate(self.erryValues):
            if erry != 0:
                delta = (self.yValues[i] - self.evolveLookUp(self.xValues[i])) / erry
                self.chi2 += delta * delta
                self.nMeas += 1

        fval[0] = self.chi2

    def prepareAndRunOptimizer(self, xValues, yValues, erryValues, fixParams, optFunction, parValues, parNames, printOutLevel = 0):
        self.xValues    = xValues
        self.yValues    = yValues
        self.erryValues = erryValues
        self.fixParams  = fixParams

        initialError = 0.1

        nBins   = len(self.xValues)
        gMinuit = TMinuit(nBins)
        gMinuit.SetFCN(optFunction)

        arglist = array('d', nBins*[0])
        ierflg  = Long(0)

        arglist[0] = 1 # 1 for chi2, 0.5 for likelihood
        gMinuit.mnexcm('SET ERR', arglist, 1, ierflg)

        # Set starting values and step sizes for parameters
        vStart = parValues
        vStep  = [initialError for i in range(len(parValues))]
        for i,name in enumerate(parNames):
            gMinuit.mnparm(i, name, vStart[i], vStep[i], 0, 0, ierflg)

        # Fix parameters (counting from 1)
        for i,fix in enumerate(self.fixParams):
            arglist[i] = fix + 1
        gMinuit.mnexcm('FIX', arglist, len(self.fixParams), ierflg)

        # Define printout level
        arglist[0] = printOutLevel
        # -1 = no output except from SHOW
        #  0 = minimum output (no starting values or intermediate results) default value, normal output
        #  1 = additional output giving intermediate results.
        #  2 = maximum output, showing progress of minimizations.
        gMinuit.mnexcm('SET PRI', arglist, 1, ierflg)

        # Now ready for minimization step
        arglist[0] = 5000 # Max calls
        arglist[1] =  0.1 # Tolerance
        gMinuit.mnexcm('MIGRAD', arglist, 2, ierflg)

        # Print results
        self.dof = self.nMeas - len(parValues) + len(self.fixParams)
        fmin, fedm, errdef  = Double(0.), Double(0.), Double(0.)
        npari, nparx, istat = Long(0), Long(0), Long(0)
        gMinuit.mnstat(fmin, fedm, errdef, npari, nparx, istat)
        print '\nFMIN:', round(fmin,2), '\tFEDM:', round(fedm,2), '\tERRDEF:', errdef, '\tNPARI:', npari, '\tNPARX:', nparx, '\tISTAT:', istat
        print 'chi-2:', round(self.chi2,2), '\td.o.f.:', self.dof, '\tchi-2/d.o.f.:', round(self.chi2/self.dof,2), '\n'

        return gMinuit, istat

    def runOptimization(self, xValues, yValues, erryValues, fixParams, doSmearing = False, printOutLevel = 0):
        self.doSmearing = doSmearing

        gMinuit, istat = self.prepareAndRunOptimizer(xValues, yValues, erryValues, fixParams, self.costFunction, self.parValues, self.parNames, printOutLevel)

        # Extract parameters and errors
        ierflg = Long(0)
        val, err, errLo, errHi = Double(0.), Double(0.), Double(0.), Double(0.)
        self.fitErr   = [0 for i in range(len(self.parValues))]
        self.fitErrLo = [0 for i in range(len(self.parValues))]
        self.fitErrHi = [0 for i in range(len(self.parValues))]
        for i in range(len(self.parValues)):
            gMinuit.mnpout(i, TString(''), val, err, errLo, errHi, ierflg)
            self.parValues[i] = float(val)
            self.fitErr[i]    = float(err)
            self.fitErrLo[i]  = float(errLo)
            self.fitErrHi[i]  = float(errHi)

        for i,value in enumerate(self.parValues):
            self.fitFun.SetParameter(i, value)

        return istat

    def runGlobalOptimization(self, evolutions, xValues, yValues, erryValues, fixParams, doSmearing = False, printOutLevel = 0):
        self.evolutions = evolutions
        self.doSmearing = doSmearing

        parNames = ['Initial population', 'Carrying capacity', 'Recovery rate', 'Growth rate-0']
        parNames.extend(['Growth rate-' + str(i+1) for i in range(len(evolutions))])
        parValues = [0 for i in range(len(parNames))]
        parValues[0] = self.parValues[0]
        parValues[1] = self.parValues[1]
        parValues[2] = self.parValues[2]
        parValues[3] = self.parValues[3]
        for i,ev in enumerate(self.evolutions):
            parValues[4+i] = ev.parValues[3]

        gMinuit, istat = self.prepareAndRunOptimizer(xValues, yValues, erryValues, fixParams, self.costFunctionGlobal, parValues, parNames, printOutLevel)

        # Extract parameters and errors
        ierflg = Long(0)
        val, err, errLo, errHi = Double(0.), Double(0.), Double(0.), Double(0.)
        self.fitErr   = [0 for i in range(len(parValues))]
        self.fitErrLo = [0 for i in range(len(parValues))]
        self.fitErrHi = [0 for i in range(len(parValues))]
        for i in range(len(parValues)):
            gMinuit.mnpout(i, TString(''), val, err, errLo, errHi, ierflg)
            parValues[i]     = float(val)
            self.fitErr[i]   = float(err)
            self.fitErrLo[i] = float(errLo)
            self.fitErrHi[i] = float(errHi)

        self.setFitFun(self.evolveGlobalWrapper, self.tStart, evolutions[-1].tStop, parNames)
        for i,value in enumerate(parValues):
            self.fitFun.SetParameter(i, value)
        return istat, parValues, parNames

    def addStats(self, parNames, parValues):
        pt = TPaveText(.15, .45, .55, .85, 'NDC')
        pt.SetBorderSize(1)
        pt.SetFillColor(0)
        pt.SetTextAlign(12)
        pt.SetTextFont(62)
        pt.AddText(' #chi^{2}/d.o.f.   ' + str(round(self.chi2,2)) + ' / ' + str(self.dof))

        for i,name in enumerate(parNames):
            if i not in self.fixParams:
                pt.AddText('')
                pt.AddText(' ' + name + '   ' + '{:.2e}'.format(parValues[i]) + ' #pm ' + '{:.2e}'.format(self.fitErr[i]))

        pt.Draw()

        return pt

    def getGraphN(self):
        graphN = TGraph()

        for x in self.bins:
            graphN.SetPoint(graphN.GetN(), x, self.evolveLookUp(x))

        graphN.SetLineColor(2)
        graphN.SetLineWidth(3)

        graphN.SetMarkerColor(2)
        graphN.SetMarkerSize(1.3)
        graphN.SetMarkerStyle(29)

        return graphN

    def getGraphR0(self, graphN):
        graphR0 = TGraph()

        for i in range(graphN.GetN() - 1):
            if Decimal(graphN.GetY()[i] * self.parValues[2]) != 0:
                graphR0.SetPoint(graphR0.GetN(), graphN.GetX()[i], Decimal((graphN.GetY()[i+1] - graphN.GetY()[i]) / self.dt) / Decimal(graphN.GetY()[i] * self.parValues[2]) + 1)
            else:
                graphR0.SetPoint(graphR0.GetN(), graphN.GetX()[i], 0)

        graphR0.SetLineColor(2)
        graphR0.SetLineWidth(3)

        graphR0.SetMarkerColor(2)
        graphR0.SetMarkerSize(1.3)
        graphR0.SetMarkerStyle(29)

        return graphR0

    def getGraphPinfect(self):
        graphP = TGraph()

        for i in range(1,int(round((self.tStop - self.tStart)/self.dt,1))):
            graphP.SetPoint(graphP.GetN(), self.tStart + i * self.dt, self.evolve(self.tStart + i * self.dt, self.parValues)[2])

        graphP.SetLineColor(2)
        graphP.SetLineWidth(3)

        graphP.SetMarkerColor(2)
        graphP.SetMarkerSize(1.3)
        graphP.SetMarkerStyle(29)

        return graphP

    def getGraphGlobalPinfect(self, evolutions, parValues):
        graphP = TGraph()

        for i in range(1,int(round((evolutions[-1].tStop - self.tStart)/self.dt,1))):
            graphP.SetPoint(graphP.GetN(), self.tStart + i * self.dt, self.evolveGlobal(evolutions, self.tStart + i * self.dt, parValues)[2])

        graphP.SetLineColor(2)
        graphP.SetLineWidth(3)

        graphP.SetMarkerColor(2)
        graphP.SetMarkerSize(1.3)
        graphP.SetMarkerStyle(29)

        return graphP

    def smearing(self, mean = 1.4, sigma = 0.3):
        nSigma = 100.
        nDaysBack = 3.

        g = [self.logNormal(i * self.dt, mean, sigma) for i in range(int(round(nSigma * sigma / self.dt,1)))]
        intro = np.ones(int(round(nDaysBack / self.dt,1))) * self.lookUpTable[0]
        self.lookUpTable = np.concatenate([intro, self.lookUpTable])
        self.lookUpTable = np.convolve(self.lookUpTable, g) * self.dt
        self.lookUpTable = self.lookUpTable[len(intro):]
        self.bins = [self.tStart + i * self.dt for i in range(len(self.lookUpTable))]

    def logNormal(self, x, mean, sigma):
        if x <= 0: return 0.
        arg = ((log(x) - mean) / sigma) if sigma != 0 else 0.
        return 1. / (x * sigma * sqrt(2.*pi)) * exp(-arg*arg / 2.)

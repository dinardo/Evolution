#################################################
# Finite difference model for disease evolution #
# Author: Mauro E. Dinardo                      #
#################################################

from array import array
from math  import sqrt, log, exp, pi

from ROOT  import TMinuit, Long, Double, TString, TPaveText, TGraph, TF1

class evolution(object):
    def __init__(self, parValues, tStart, tStop, totalPopulation,
                 symptomaticFraction     = 0.3,
                 transmissionProbability = 0.3,
                 historyActiveDt         = 0.):

        self.parNames  = ['Initial population', 'Growth rate', 'Recovery rate', 'Carrying capacity']
        self.parValues = parValues

        self.tStart = tStart
        self.tStop  = tStop

        ####################
        # Model parameters #
        ####################
        self.totalPopulation         = totalPopulation
        self.symptomaticFraction     = symptomaticFraction
        self.transmissionProbability = transmissionProbability
        self.historyActiveDt         = historyActiveDt

        self.dt = 0.2 # [Day]

        self.nMeas, self.chi2, self.dof = 0., 0., 0.

        self.fitFun = TF1('Evolution', self.evolveActiveWrapper, tStart, tStop, len(self.parNames))

    def evolveActiveWrapper(self, t, par):
        return self.evolveActive(t[0],par)[0]

    def evolveActive(self, t, par):
        T = t - self.tStart
        N = par[0] / self.symptomaticFraction
        CC = par[3]
        r = 0.
        CCn = 0.
        totalInfected = 0.
        historyActiveDt = self.historyActiveDt / self.symptomaticFraction

        for n in range(int(round(T/self.dt,1))):
            Nn = N
            totalInfected = N + historyActiveDt * par[2]
            r = par[1] * (1. - totalInfected / CC)
            N += self.dt * N * (r - par[2])
            CCn = CC
            CC += self.dt * par[1] / self.transmissionProbability * (N - Nn * (1. - par[2])) * (1. - CC / self.totalPopulation)
            historyActiveDt += Nn * self.dt

        return [N * self.symptomaticFraction, (totalInfected / CCn) if CCn != 0. else 0., CC]

    def sumHistoryActive(self, up2Time):
        return self.totalRecovered(up2Time) / self.parValues[2]

    def totalInfected(self, up2Time):
        return self.evolveActive(up2Time, self.parValues)[0] + self.totalRecovered(up2Time)

    def totalRecovered(self, up2Time):
        historyActiveDt = self.historyActiveDt

        for i in range(int(round((up2Time - self.tStart)/self.dt,1))):
            historyActiveDt += self.evolveActive(self.tStart + i*self.dt, self.parValues)[0] * self.dt

        return historyActiveDt * self.parValues[2]

    def myChi2(self, npar, grad, fval, par, iflag):
        self.nMeas, self.chi2, delta = 0., 0., 0.

        for i,erry in enumerate(self.erryValues):
            if erry != 0:
                delta = (self.yValues[i] - self.evolveActive(self.xValues[i], par)[0]) / erry
                self.chi2 += delta * delta
                self.nMeas += 1

        fval[0] = self.chi2

    def runOptimization(self, xValues, yValues, erryValues, fixParams):
        self.xValues    = xValues
        self.yValues    = yValues
        self.erryValues = erryValues
        self.fixParams  = fixParams

        initialError = 0.1

        nBins   = len(self.xValues)
        gMinuit = TMinuit(nBins)
        gMinuit.SetFCN(self.myChi2)

        arglist = array('d', nBins*[0])
        ierflg  = Long(0)

        arglist[0] = 1 # 1 for chi2, 0.5 for likelihood
        gMinuit.mnexcm('SET ERR', arglist, 1, ierflg)

        # Set starting values and step sizes for parameters
        vStart = self.parValues
        vStep  = [initialError for i in range(len(self.parNames))]
        for i,name in enumerate(self.parNames):
            gMinuit.mnparm(i, name, vStart[i], vStep[i], 0, 0, ierflg)

        # Fix parameters (counting from 1)
        for i,fix in enumerate(self.fixParams):
            arglist[i] = fix + 1
        gMinuit.mnexcm('FIX', arglist, len(self.fixParams), ierflg)

        # Now ready for minimization step
        arglist[0] = 1000 # Max calls
        arglist[1] =  1.0 # Tolerance
        gMinuit.mnexcm('MIGRAD', arglist, 2, ierflg)

        # Print results
        self.dof = self.nMeas - len(self.parNames) + len(self.fixParams)
        fmin, fedm, errdef  = Double(0.), Double(0.), Double(0.)
        npari, nparx, istat = Long(0), Long(0), Long(0)
        gMinuit.mnstat(fmin, fedm, errdef, npari, nparx, istat)
        print '\nFMIN:', round(fmin,2), '\tFEDM:', round(fedm,2), '\tERRDEF:', errdef, '\tNPARI:', npari, '\tNPARX:', nparx, '\tISTAT:', istat
        print 'chi-2:', round(self.chi2,2), '\td.o.f.:', self.dof, '\tchi-2/d.o.f.:', round(self.chi2/self.dof,2), '\n'

        # Extract parameters and errors
        val, err, errLo, errHi = Double(0.), Double(0.), Double(0.), Double(0.)
        self.fitErr   = [0 for i in range(len(self.parNames))]
        self.fitErrLo = [0 for i in range(len(self.parNames))]
        self.fitErrHi = [0 for i in range(len(self.parNames))]
        for i in range(len(self.parNames)):
            gMinuit.mnpout(i, TString(''), val, err, errLo, errHi, ierflg)
            self.parValues[i] = float(val)
            self.fitErr[i]    = float(err)
            self.fitErrLo[i]  = float(errLo)
            self.fitErrHi[i]  = float(errHi)

        for i,value in enumerate(self.parValues):
            self.fitFun.SetParameter(i, value)

    def addStats(self):
        pt = TPaveText(.15, .45, .55, .85, 'NDC')
        pt.SetBorderSize(1)
        pt.SetFillColor(0)
        pt.SetTextAlign(12)
        pt.SetTextFont(62)
        pt.AddText(' #chi^{2}/d.o.f.   ' + str(round(self.chi2,2)) + ' / ' + str(self.dof))

        for i,name in enumerate(self.parNames):
            if i not in self.fixParams:
                pt.AddText('')
                pt.AddText(' ' + name + '   ' + '{:.2e}'.format(self.parValues[i]) + ' #pm ' + '{:.2e}'.format(self.fitErr[i]))

        pt.Draw()

        return pt

    def getGraphN(self):
        graphN = TGraph()

        for i in range(int(round((self.tStop - self.tStart)/self.dt,1))):
            graphN.SetPoint(graphN.GetN(), self.tStart + i*self.dt, self.evolveActive(self.tStart + i*self.dt, self.parValues)[0])

        graphN.SetLineColor(2)
        graphN.SetLineWidth(3)

        graphN.SetMarkerColor(2)
        graphN.SetMarkerSize(1.3)
        graphN.SetMarkerStyle(29)

        return graphN

    def getGraphPinfect(self):
        graphP = TGraph()

        for i in range(1,int(round((self.tStop - self.tStart)/self.dt,1))):
            graphP.SetPoint(graphP.GetN(), self.tStart + i*self.dt, self.evolveActive(self.tStart + i*self.dt, self.parValues)[1])

        graphP.SetLineColor(2)
        graphP.SetLineWidth(3)

        graphP.SetMarkerColor(2)
        graphP.SetMarkerSize(1.3)
        graphP.SetMarkerStyle(29)

        return graphP

    def getGraphR0(self, graphN):
        graphR0 = TGraph()

        for i in range(graphN.GetN() - 1):
            graphR0.SetPoint(graphR0.GetN(), graphN.GetX()[i], (graphN.GetY()[i+1] - graphN.GetY()[i]) / self.dt / (graphN.GetY()[i]*self.parValues[2]) + 1)

        graphR0.SetLineColor(2)
        graphR0.SetLineWidth(3)

        graphR0.SetMarkerColor(2)
        graphR0.SetMarkerSize(1.3)
        graphR0.SetMarkerStyle(29)

        return graphR0

    def combineEvolutions(self, parList, timeList, totalPopulation, symptomaticFraction, transmissionProbability):
        myGraphN = self.getGraphN()

        historyActive = self.sumHistoryActive(timeList[0])
        [active, Pinf, CC] = self.evolveActive(timeList[0], self.parValues)

        for t,par in enumerate(parList):
            par[0] = active
            par[3] = CC if par[3] == 0 else par[3]
            evolve = evolution(par, timeList[t], timeList[t+1], totalPopulation, symptomaticFraction, transmissionProbability, historyActive)

            historyActive += evolve.sumHistoryActive(timeList[t+1])
            [active, Pinf, CC] = evolve.evolveActive(timeList[t+1], evolve.parValues)

            graphN = evolve.getGraphN()

            for i in range(graphN.GetN()):
                myGraphN.SetPoint(myGraphN.GetN(), graphN.GetX()[i], graphN.GetY()[i])

        return myGraphN

    def smearing(self, graph, sigma = 2.):
        mean    =   0.
        nSigma  = 100.
        myGraph = TGraph()

        normalize = 0.
        for i in range(int(round(graph.GetN() + nSigma*sigma/self.dt,1))):
            normalize += self.logNormal(i*self.dt, mean, sigma)
        normalize *= self.dt

        for i in range(int(round(graph.GetN() + nSigma*sigma/self.dt,1))):
            convolve = 0.
            for j in range(graph.GetN()):
                convolve += graph.GetY()[j] * self.logNormal(graph.GetX()[0] + i*self.dt - graph.GetX()[j], mean, sigma)

            myGraph.SetPoint(myGraph.GetN(), graph.GetX()[0] + i*self.dt, convolve * self.dt / normalize)

        myGraph.SetLineColor(2)
        myGraph.SetLineWidth(3)

        myGraph.SetMarkerColor(2)
        myGraph.SetMarkerSize(1.3)
        myGraph.SetMarkerStyle(29)

        return myGraph

    def logNormal(self, x, mean, sigma):
        if x <= 0: return 0.
        arg = ((log(x) - mean) / sigma) if sigma != 0 else 0.
        return 1. / (x * sigma * sqrt(2.*pi)) * exp(-arg*arg / 2.)

import ROOT
import sys
import math

# global counters
NAME_COUNT = 0

def cheb_fn(x, degree, kind):
  if kind == 1:
    if degree == 0:
        return 1
    elif degree == 1:
        return x
    else:
        t_minus_2 = 1
        t_minus_1 = x
        for n in range(2, degree + 1):
            current_term = 2 * x * t_minus_1 - t_minus_2
            t_minus_2 = t_minus_1
            t_minus_1 = current_term
        return current_term
  elif kind == 2:
    if degree == 0:
        return 1
    elif degree == 1:
        return 2 * x
    else:
        u_minus_2 = 1
        u_minus_1 = 2 * x
        for n in range(2, degree + 1):
            current_term = 2 * x * u_minus_1 - u_minus_2
            u_minus_2 = u_minus_1
            u_minus_1 = current_term
        return current_term
  else: raise ValueError('"kind" argument must be 1 or 2.')

def getname(prefix='obj'):
  '''
  helper to return unique names for ROOT objects or global variables
  '''
  global NAME_COUNT
  NAME_COUNT += 1
  return prefix+str(NAME_COUNT)

def TemplateToHistogram(func, bins, low, high, integral=False):
  '''
  Takes ROOT TF1 function template, and a binning

  Returns a ROOT TH1D histogram vesion
  '''
  name = getname()  
  hist = ROOT.TH1D(name, name, bins, low, high)
  for i in range(hist.GetNbinsX()):
    if not integral:
      val = func.Eval(hist.GetBinCenter(i+1))
      hist.SetBinContent(i+1, val) if not math.isnan(val) else hist.SetBinContent(i+1, 0)
    else:
      val = (func.Integral(hist.GetBinLowEdge(i+1), hist.GetBinLowEdge(i+1) + hist.GetBinWidth(i+1))) / hist.GetBinWidth(i+1)
      hist.SetBinContent(i+1, val) if not math.isnan(val) else hist.SetBinContent(i+1, 0)
  return hist

def HistogramToFunction(hist):
  '''
  Takes a ROOT TH1 histogram

  Returns a linearized version as a ROOT TF1
  '''
  def histfunc(x):
    return hist.GetBinContent(hist.FindBin(x[0])) 
  return histfunc

def MultiplyWithPolyToTF1(func, degree, range_low=0, range_high=50, cheb=0, parameters=None):
  '''
  Takes a python function

  Returns a TF1 object representing the input function times a polynomial, and the returns associated python function as well

  when cheb=0 (default) use regular polynomials (1, x^2, x^3, etc)
  when cheb=1 use Chebyshev polynomials of the first kind
  when cheb=2 use Chebyshev polynomials of the second kind
  '''
  if degree == 0 and cheb == 0:
    def func_after_mult(x, p):
      return func(x) * (p[0])
  if degree == 1 and cheb == 0:
    def func_after_mult(x, p):
      return func(x) * (p[0] + p[1]*x[0])
  if degree == 2 and cheb == 0:
    def func_after_mult(x, p):
      return func(x) * (p[0] + p[1]*x[0] + p[2]*(x[0]**2))
  if degree == 3 and cheb == 0:
    def func_after_mult(x, p):
      return func(x) * (p[0] + p[1]*x[0] + p[2]*(x[0]**2) + p[3]*(x[0]**3))
  if degree == 4 and cheb == 0:
    def func_after_mult(x, p):
      return func(x) * (p[0] + p[1]*x[0] + p[2]*(x[0]**2) + p[3]*(x[0]**3)
                        + p[4]*(x[0]**4))
  if degree == 5 and cheb == 0:
    def func_after_mult(x, p):
      return func(x) * (p[0] + p[1]*x[0] + p[2]*(x[0]**2) + p[3]*(x[0]**3)
                        + p[4]*(x[0]**4) + p[5]*x[0]**5)
  if degree == 6 and cheb == 0:
    def func_after_mult(x, p):
      return func(x) * (p[0] + p[1]*x[0] + p[2]*(x[0]**2) + p[3]*(x[0]**3)
                        + p[4]*(x[0]**4) + p[5]*x[0]**5 + p[6]*x[0]**6)
  if degree == 7 and cheb == 0:
    def func_after_mult(x, p):
      return func(x) * (p[0] + p[1]*x[0] + p[2]*(x[0]**2) + p[3]*(x[0]**3)
                        + p[4]*(x[0]**4) + p[5]*x[0]**5 + p[6]*x[0]**6 + p[7]*x[0]**7)
  if degree == 8 and cheb == 0:
    def func_after_mult(x, p):
      return func(x) * (p[0] + p[1]*x[0] + p[2]*(x[0]**2) + p[3]*(x[0]**3) 
                        + p[4]*(x[0]**4) + p[5]*x[0]**5 + p[6]*x[0]**6 + p[7]*x[0]**7
                        + p[8]*x[0]**8)

  if degree == 0 and cheb == 1:
    def func_after_mult(x, p):
      return func(x) * (p[0])
  if degree == 1 and cheb == 1:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 1))
  if degree == 2 and cheb == 1:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 1) + p[2]*cheb_fn(X, 2, 1))
  if degree == 3 and cheb == 1:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 1) + p[2]*cheb_fn(X, 2, 1)
                        + p[3]*cheb_fn(X, 3, 1))
  if degree == 4 and cheb == 1:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 1) + p[2]*chb_fn(X, 2, 1)
                        + p[3]*cheb_fn(X, 3, 1) + p[4]*cheb_fn(X, 4, 1))
  if degree == 5 and cheb == 1:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 1) + p[2]*cheb_fn(X, 2, 1)
                        + p[3]*cheb_fn(X, 3, 1) + p[4]*cheb_fn(X, 4, 1)
                        + p[5]*cheb_fn(X, 5, 1))
  if degree == 6 and cheb == 1:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 1) + p[2]*cheb_fn(X, 2, 1)
                        + p[3]*cheb_fn(X, 3, 1) + p[4]*cheb_fn(X, 4, 1)
                        + p[5]*cheb_fn(X, 5, 1) + p[6]*cheb_fn(X, 6, 1))
  if degree == 7 and cheb == 1:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 1) + p[2]*cheb_fn(X, 2, 1)
                        + p[3]*cheb_fn(X, 3, 1) + p[4]*cheb_fn(X, 4, 1)
                        + p[5]*cheb_fn(X, 5, 1) + p[6]*cheb_fn(X, 6, 1)
                        + p[7]*cheb_fn(X, 7, 1))
  if degree == 7 and cheb == 1:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 1) + p[2]*cheb_fn(X, 2, 1)
                        + p[3]*cheb_fn(X, 3, 1) + p[4]*cheb_fn(X, 4, 1)
                        + p[5]*cheb_fn(X, 5, 1) + p[6]*cheb_fn(X, 6, 1)
                        + p[7]*cheb_fn(X, 7, 1) + p[6]*cheb_fn(X, 8, 1))
  if degree == 8 and cheb == 1:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 1) + p[2]*cheb_fn(X, 2, 1)
                        + p[3]*cheb_fn(X, 3, 1) + p[4]*cheb_fn(X, 4, 1)
                        + p[5]*cheb_fn(X, 5, 1) + p[6]*cheb_fn(X, 6, 1)
                        + p[7]*cheb_fn(X, 7, 1) + p[6]*cheb_fn(X, 8, 1)
                        + p[8]*cheb_fn(X, 8, 1))

  if degree == 0 and cheb == 2:
    def func_after_mult(x, p):
      return func(x) * (p[0])
  if degree == 1 and cheb == 2:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 2))
  if degree == 2 and cheb == 2:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 2) + p[2]*cheb_fn(X, 2, 2))
  if degree == 3 and cheb == 2:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 2) + p[2]*cheb_fn(X, 2, 2)
                        + p[3]*cheb_fn(X, 3, 2))
  if degree == 4 and cheb == 2:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 2) + p[2]*chb_fn(X, 2, 2)
                        + p[3]*cheb_fn(X, 3, 2) + p[4]*cheb_fn(X, 4, 2))
  if degree == 5 and cheb == 2:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 2) + p[2]*cheb_fn(X, 2, 2)
                        + p[3]*cheb_fn(X, 3, 2) + p[4]*cheb_fn(X, 4, 2)
                        + p[5]*cheb_fn(X, 5, 2))
  if degree == 6 and cheb == 2:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 2) + p[2]*cheb_fn(X, 2, 2)
                        + p[3]*cheb_fn(X, 3, 2) + p[4]*cheb_fn(X, 4, 2)
                        + p[5]*cheb_fn(X, 5, 2) + p[6]*cheb_fn(X, 6, 2))
  if degree == 7 and cheb == 2:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 2) + p[2]*cheb_fn(X, 2, 2)
                        + p[3]*cheb_fn(X, 3, 2) + p[4]*cheb_fn(X, 4, 2)
                        + p[5]*cheb_fn(X, 5, 2) + p[6]*cheb_fn(X, 6, 2)
                        + p[7]*cheb_fn(X, 7, 2))
  if degree == 7 and cheb == 2:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 2) + p[2]*cheb_fn(X, 2, 2)
                        + p[3]*cheb_fn(X, 3, 2) + p[4]*cheb_fn(X, 4, 2)
                        + p[5]*cheb_fn(X, 5, 2) + p[6]*cheb_fn(X, 6, 2)
                        + p[7]*cheb_fn(X, 7, 2) + p[6]*cheb_fn(X, 8, 2))
  if degree == 8 and cheb == 2:
    def func_after_mult(x, p):
      X = x[0]
      return func(x) * (p[0] + p[1]*cheb_fn(X, 1, 2) + p[2]*cheb_fn(X, 2, 2)
                        + p[3]*cheb_fn(X, 3, 2) + p[4]*cheb_fn(X, 4, 2)
                        + p[5]*cheb_fn(X, 5, 2) + p[6]*cheb_fn(X, 6, 2)
                        + p[7]*cheb_fn(X, 7, 2) + p[6]*cheb_fn(X, 8, 2)
                        + p[8]*cheb_fn(X, 8, 2))

  globals()[getname('func')] = func_after_mult
  tf1 = ROOT.TF1(getname(), func_after_mult, range_low, range_high, degree+1)

  if degree>=0: tf1.SetParName(0, 'Constant') if cheb==0 else tf1.SetParName(0, 'Zero')
  if degree>=1: tf1.SetParName(1, 'Linear') if cheb==0 else tf1.SetParName(1, 'One')
  if degree>=2: tf1.SetParName(2, 'Quadratic') if cheb==0 else tf1.SetParName(2, 'Two')
  if degree>=3: tf1.SetParName(3, 'Cubic') if cheb==0 else tf1.SetParName(3, 'Three')
  if degree>=4: tf1.SetParName(4, 'Quartic') if cheb==0 else tf1.SetParName(4, 'Four')
  if degree>=5: tf1.SetParName(5, 'Quintic') if cheb==0 else tf1.SetParName(5, 'Five')
  if degree>=6: tf1.SetParName(6, 'Sextic') if cheb==0 else tf1.SetParName(6, 'Six')
  if degree>=7: tf1.SetParName(7, 'Septic') if cheb==0 else tf1.SetParName(7, 'Seven')
  if degree>=8: tf1.SetParName(8, 'Octic') if cheb==0 else tf1.SetParName(8, 'Eight')

  if not parameters:
    for i in range(degree+1): tf1.SetParameter(i, 0.1**i)
  else:
    tf1.SetParameters(*parameters)
  return tf1, func_after_mult

def fit_hist(hist, function, range_low, range_high, N=1, initial_guesses=None, integral=False):
  '''
  Takes a historgram, fits a function and returns a TF1 and the fit result

  function: can be set to 'exp' or 'landau' or 'full'
  N:        controls how many exponentials or landaus to use for 'exp' and 'landu'
            if function is 'full' N can be 11, 12, 13, 14, 22, 23, 23, 24
            where the tens digit is number of landaus, and ones digit is number of exps
  '''
  if function == 'landau' and N == 1:
    def python_func(x, p):
      norm = p[0]
      mpv = p[1]
      sigma = p[2]
      land = norm * ROOT.TMath.Landau(x[0], mpv, sigma)
      return land
    NPAR = 3
    if not initial_guesses: initial_guesses = [hist.GetEntries(), hist.GetMean(), 0.25]
    if not len(initial_guesses) == NPAR:
        raise AssertionError('Length of initial guesses list must be '+str(NPAR)+"!")
    tf1 = ROOT.TF1(getname('func'), python_func, range_low, range_high, NPAR)
    tf1.SetParNames("Constant", "MPV", "Sigma")
    tf1.SetParameters(*initial_guesses)

  if function == 'landau' and N == 2:
    def python_func(x, p):
      norm = p[0]
      mpv1 = p[1]
      sigma1 = p[2]
      bound = p[3]
      mpv2 = p[4]
      sigma2 = p[5]
      land1 = norm * ROOT.TMath.Landau(x[0], mpv1, sigma1)
      y1 = norm * ROOT.TMath.Landau(bound, mpv1, sigma1)
      y2 = ROOT.TMath.Landau(bound, mpv2, sigma2)
      land2 = (y1/y2) * ROOT.TMath.Landau(x[0], mpv2, sigma2)
      if x[0] < bound: return land1
      else: return land2
    NPAR = 6
    if not initial_guesses: initial_guesses = [
      hist.GetEntries(), hist.GetMean(), 0.25, (range_low+range_high)/2.0, hist.GetMean(), 0.25]
    if not len(initial_guesses) == NPAR:
        raise AssertionError('Length of initial guesses list must be '+str(NPAR)+"!")
    tf1 = ROOT.TF1(getname('func'), python_func, range_low, range_high, NPAR)
    tf1.SetParNames("Constant", "MPV1", "Sigma1", "bound", "MPV2", "Sigma2")
    tf1.SetParameters(*initial_guesses)
    tf1.SetParLimits(3, 0, range_high+1)
    tf1.SetParLimits(4, 0, hist.GetMean()*5)
    tf1.SetParLimits(5, 0, 1e5)

  if function == 'exp' and N == 1:
    def python_func(x, p):
      norm = p[0]
      C1 = p[1]
      exp1 = norm * ROOT.TMath.Exp(C1*x[0])
      return exp1
    NPAR = 2
    if not initial_guesses: initial_guesses = [hist.GetEntries(), -1]
    if not len(initial_guesses) == NPAR:
        raise AssertionError('Length of initial guesses list must be '+str(NPAR)+"!")
    tf1 = ROOT.TF1(getname('func'), python_func, range_low, range_high, NPAR)
    tf1.SetParNames("Constant", "C1")
    tf1.SetParameters(*initial_guesses)

  if function == 'exp' and N == 2:
    def python_func(x, p):
      norm = p[0]
      C1 = p[1]
      bound = p[2]
      C2 = p[3]
      y1 = norm * ROOT.TMath.Exp(C1*bound)
      y2 = ROOT.TMath.Exp(C2*bound)
      exp1 = norm * ROOT.TMath.Exp(C1*x[0])
      exp2 = ROOT.TMath.Exp(C2*x[0]) * (y1/y2)
      if x[0] < bound: return exp1
      else: return exp2
    NPAR = 4
    if not initial_guesses: initial_guesses = [
      hist.Integral(hist.FindBin(range_low), hist.FindBin(range_high)), -1, (range_low+range_high)/2.0, -1]
    if not len(initial_guesses) == NPAR:
        raise AssertionError('Length of initial guesses list must be '+str(NPAR)+"!")
    tf1 = ROOT.TF1(getname('func'), python_func, range_low, range_high, NPAR)
    tf1.SetParNames("Constant", "C1", "bound", "C2")
    tf1.SetParameters(*initial_guesses)
    tf1.SetParLimits(2, 0, hist.GetBinLowEdge(hist.GetNbinsX()))
    tf1.SetParLimits(3, -10, 0)

  if function == 'exp' and N == 3:
    def python_func(x, p):
      norm = p[0]
      C1 = p[1]
      bound1 = p[2]
      C2 = p[3]
      bound12 = p[4]
      C3 = p[5]
      bound2 = bound1 + bound12
      exp1 = norm * ROOT.TMath.Exp(C1*x[0])
      y1 = norm * ROOT.TMath.Exp(C1*bound1)
      y2 = ROOT.TMath.Exp(C2*bound1)
      exp2 = ROOT.TMath.Exp(C2*x[0]) * (y1/y2)
      y3 = ROOT.TMath.Exp(C2*bound2)
      y4 = ROOT.TMath.Exp(C3*bound2)
      exp3 = ROOT.TMath.Exp(C3*x[0]) * (y3/y4) * (y1/y2)
      if x[0] < bound1: return exp1
      elif x[0] < bound2: return exp2
      else: return exp3
    NPAR = 6
    if not initial_guesses: initial_guesses = [
      hist.Integral(hist.FindBin(range_low), hist.FindBin(range_high)), -1, (range_low+range_high)/2.0,
      -1, (range_low+range_high)/2.0, -1]
    if not len(initial_guesses) == NPAR:
        raise AssertionError('Length of initial guesses list must be '+str(NPAR)+"!")
    tf1 = ROOT.TF1(getname('func'), python_func, range_low, range_high, NPAR)
    tf1.SetParNames("Constant", "C1", "bound1", "C2", "bound12", "C3")
    tf1.SetParameters(*initial_guesses)
    tf1.SetParLimits(2, 0, hist.GetBinLowEdge(hist.GetNbinsX()))
    tf1.SetParLimits(3, -10, 0)
    tf1.SetParLimits(4, 0, hist.GetBinLowEdge(hist.GetNbinsX()))
    tf1.SetParLimits(5, -10, 0)

  if function == 'full' and N == 24:
    def python_func(x, p):
        norm1 = p[0]
        mpv1 = p[1]
        sigma1 = p[2]
        norm2 = p[3]
        mpv2 = p[4]
        sigma2 = p[5]
        C1 = p[6]
        C2 = p[7]
        C3 = p[8]
        C4 = p[9]
        bound1 = p[10]
        b12 = p[11]
        b23 = p[12]
        b34 = p[13]
        b45 = p[14]

        if bound1 < 0: bound1 = 0
        bound2 = bound1 + b12
        bound3 = bound2 + b23
        bound4 = bound3 + b34
        bound5 = bound4 + b45

        land1 = norm1 * ROOT.TMath.Landau(x[0], mpv1, sigma1)
      
        y11=norm1*ROOT.TMath.Landau(bound1, mpv1, sigma1)
        y12=norm2*ROOT.TMath.Landau(bound1, mpv2, sigma2)
        land2 = norm2*ROOT.TMath.Landau(x[0], mpv2, sigma2)*y11/y12
         
        y21=norm2*ROOT.TMath.Landau(bound2, mpv2, sigma2)*y11/y12
        y22=ROOT.TMath.Exp(C1*bound2)
        exp1=ROOT.TMath.Exp(C1*x[0])*y21/y22

        y31=ROOT.TMath.Exp(C1*bound3)*y21/y22
        y32=ROOT.TMath.Exp(C2*bound3)
        exp2=ROOT.TMath.Exp(C2*x[0])*y31/y32
         
        y41=ROOT.TMath.Exp(C2*bound4)*y31/y32
        y42=ROOT.TMath.Exp(C3*bound4)
        exp3=ROOT.TMath.Exp(C3*x[0])*y41/y42

        y51=ROOT.TMath.Exp(C3*bound5)*y41/y42
        y52=ROOT.TMath.Exp(C4*bound5)
        exp4=ROOT.TMath.Exp(C4*x[0])*y51/y52
        
        if x[0] < bound1: return land1
        elif x[0] < bound2: return land2
        elif x[0] < bound3: return exp1
        elif x[0] < bound4: return exp2
        elif x[0] < bound5: return exp3
        else: return exp4
    NPAR = 15
    if not initial_guesses:
        nEntries = hist.GetEntries()
        mean = hist.GetMean()
        initial_guesses = [
            nEntries, mean, 0.5, nEntries, mean, 0.5, -3, -1, -0.5, -0.25,
            mean, mean/2, mean/2, mean/2, mean/2]
    if not len(initial_guesses) == NPAR:
        raise AssertionError('Length of initial guesses list must be '+str(NPAR)+"!")
    tf1 = ROOT.TF1(getname('func'), python_func, range_low, range_high, NPAR)
    tf1.SetParNames("Constant1","MPV1","Sigma1","Constant2","MPV2","Sigma2","C1","C2","C3","C4","Boundary1")
    tf1.SetParName(11, "BoundDiff12")
    tf1.SetParName(12, "BoundDiff23")
    tf1.SetParName(13, "BoundDiff34")
    tf1.SetParName(14, "BoundDiff45")
    for i, guess in enumerate(initial_guesses): tf1.SetParameter(i, guess)
    tf1.SetParLimits(5, 0, 100)
    tf1.SetParLimits(6, -10, 0)
    tf1.SetParLimits(7, -10, 0)
    tf1.SetParLimits(8, -10, 0)
    tf1.SetParLimits(9, -10, 0)
    tf1.SetParLimits(10, 0, 25)
    tf1.SetParLimits(11, 0.2, 7)
    tf1.SetParLimits(12, 0.2, 7)
    tf1.SetParLimits(13, 0.2, 7)
    tf1.SetParLimits(14, 0.2, 7)

  if function == 'full' and N == 24:
    def python_func(x, p):
        norm1 = p[0]
        mpv1 = p[1]
        sigma1 = p[2]
        norm2 = p[3]
        mpv2 = p[4]
        sigma2 = p[5]
        C1 = p[6]
        C2 = p[7]
        C3 = p[8]
        C4 = p[9]
        bound1 = p[10]
        b12 = p[11]
        b23 = p[12]
        b34 = p[13]
        b45 = p[14]

        if bound1 < 0: bound1 = 0
        bound2 = bound1 + b12
        bound3 = bound2 + b23
        bound4 = bound3 + b34
        bound5 = bound4 + b45

        land1 = norm1 * ROOT.TMath.Landau(x[0], mpv1, sigma1)
      
        y11=norm1*ROOT.TMath.Landau(bound1, mpv1, sigma1)
        y12=norm2*ROOT.TMath.Landau(bound1, mpv2, sigma2)
        land2 = norm2*ROOT.TMath.Landau(x[0], mpv2, sigma2)*y11/y12
         
        y21=norm2*ROOT.TMath.Landau(bound2, mpv2, sigma2)*y11/y12
        y22=ROOT.TMath.Exp(C1*bound2)
        exp1=ROOT.TMath.Exp(C1*x[0])*y21/y22

        y31=ROOT.TMath.Exp(C1*bound3)*y21/y22
        y32=ROOT.TMath.Exp(C2*bound3)
        exp2=ROOT.TMath.Exp(C2*x[0])*y31/y32
         
        y41=ROOT.TMath.Exp(C2*bound4)*y31/y32
        y42=ROOT.TMath.Exp(C3*bound4)
        exp3=ROOT.TMath.Exp(C3*x[0])*y41/y42

        y51=ROOT.TMath.Exp(C3*bound5)*y41/y42
        y52=ROOT.TMath.Exp(C4*bound5)
        exp4=ROOT.TMath.Exp(C4*x[0])*y51/y52
        
        if x[0] < bound1: return land1
        elif x[0] < bound2: return land2
        elif x[0] < bound3: return exp1
        elif x[0] < bound4: return exp2
        elif x[0] < bound5: return exp3
        else: return exp4
    NPAR = 15
    if not initial_guesses:
        nEntries = hist.GetEntries()
        mean = hist.GetMean()
        initial_guesses = [
            nEntries, mean, 0.5, nEntries, mean, 0.5, -3, -1, -0.5, -0.25,
            mean, mean/2, mean/2, mean/2, mean/2]
    if not len(initial_guesses) == NPAR:
        raise AssertionError('Length of initial guesses list must be '+str(NPAR)+"!")
    tf1 = ROOT.TF1(getname('func'), python_func, range_low, range_high, NPAR)
    tf1.SetParNames("Constant1","MPV1","Sigma1","Constant2","MPV2","Sigma2","C1","C2","C3","C4","Boundary1")
    tf1.SetParName(11, "BoundDiff12")
    tf1.SetParName(12, "BoundDiff23")
    tf1.SetParName(13, "BoundDiff34")
    tf1.SetParName(14, "BoundDiff45")
    for i, guess in enumerate(initial_guesses): tf1.SetParameter(i, guess)
    tf1.SetParLimits(5, 0, 100)
    tf1.SetParLimits(6, -10, 0)
    tf1.SetParLimits(7, -10, 0)
    tf1.SetParLimits(8, -10, 0)
    tf1.SetParLimits(9, -10, 0)
    tf1.SetParLimits(10, 0, 25)
    tf1.SetParLimits(11, 0.2, 7)
    tf1.SetParLimits(12, 0.2, 7)
    tf1.SetParLimits(13, 0.2, 7)
    tf1.SetParLimits(14, 0.2, 7)

  globals()[getname('func')] = python_func
  fit_string = '0SL'
  if integral: fit_string += 'I'
  fit_result = hist.Fit(tf1, fit_string, "", range_low, range_high)
  return tf1, fit_result

if __name__ == '__main__':
  RANGE_LOW = 0 
  RANGE_HIGH = 20
  BINS = 20

  print('define a function')
  f1 = ROOT.TF1('f1', 'gaus(0)', RANGE_LOW, RANGE_HIGH)
  f1.SetParNames('first', 'second', 'third')
  f1.SetParameters(1.0, 1.0, 1.0)

  print('turn into histogram')
  bins, low, high = BINS, RANGE_LOW, RANGE_HIGH
  hist = TemplateToHistogram(f1, bins, low, high, integral=True)
  hist.Draw()
  raw_input()

  print('turn histo into python function')
  func = HistogramToFunction(hist)

  print('multiply with polynomial')
  func_with_poly = MultiplyWithPolyToTF1(func, 2)
  func_with_poly.Draw()
  raw_input()

  print('make second histogram')
  hist_target = ROOT.TH1D(getname(), 'target', bins, low, high)
  hist_target.FillRandom('f1', 100000)
  hist_target.Draw()
  raw_input() 

  print('fit second histogram with function')
  hist_target.Fit(func_with_poly, "L")
  hist_target.Draw()
  raw_input() 

  print('extract parameters of poly from fit')
  f = func_with_poly
  n = f.GetNpar()
  for i in range(n):
    print(f.GetParName(i), f.GetParameter(i))

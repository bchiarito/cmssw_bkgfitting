import ROOT
import sys
import math
from scipy.stats import f 

rand = ROOT.TRandom3()
BERN_UPPER_RANGE = 25
NAME_COUNT = 0 # global counter
pt_bins = [20,40,60,80,100,140,180,220,300,380]

# "scale" a control-region tight photon distribution to its corresponding signal distribution
def removeEntries(bkg_hist, sig_hist, seed=None, cutoff=100):
    if seed: rand.SetSeed(seed)

    if bkg_hist.Integral() == 0: return False
   
    p = float(sig_hist.Integral()) / bkg_hist.Integral() # probability of removing entry
    if p >= 1: return False
    
    for i in range(1, bkg_hist.GetNbinsX()+1):
        N = round(bkg_hist.GetBinContent(i))
        if N < cutoff:
            bkg_hist.SetBinContent(i, round(rand.Binomial(N, p)))
        else: 
            content = rand.Gaus(N*p, (N*p*(1-p))**0.5)
            if content < 0: bkg_hist.SetBinContent(i, 0) 
            else: bkg_hist.SetBinContent(i, round(content))

# check consecutive pull bins to ensure there are no significant bumps
def checkPullHist(pull_hist, nBins, sigma):
    bad_pull = False
    for i in range(1, pull_hist.GetNbinsX()+1):
        if i == pull_hist.GetNbinsX() - 6: break
        if pull_hist.GetBinContent(i) > sigma:
            bad_pull = True
            for j in range(i+1, i+nBins):
                if pull_hist.GetBinContent(j) > sigma: continue
                else: 
                    bad_pull = False
                    break
        if pull_hist.GetBinContent(i) < -sigma:
            bad_pull = True
            for j in range(i+1, i+nBins):
                if pull_hist.GetBinContent(j) < -sigma: continue
                else:
                    bad_pull = False
                    break
    return bad_pull

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

def bern_fn(x, n_degree, i_degree):
    x_adj = x / BERN_UPPER_RANGE
    if i_degree < 0 or i_degree > n_degree: return 0
    #return scipy.special.comb(n_degree, i_degree) * x_adj**i_degree * (1 - x_adj)**(n_degree - i_degree)
    return (math.factorial(n_degree))/(math.factorial(i_degree)*math.factorial(n_degree - i_degree)) * x_adj**i_degree * (1 - x_adj)**(n_degree - i_degree)

def getname(prefix='obj'):
  '''
  helper to return unique names for ROOT objects or global variables
  '''
  global NAME_COUNT
  NAME_COUNT += 1
  return prefix+str(NAME_COUNT)

def TemplateToHistogram(func, bins, low, high, integral=False):
    '''
    Convert a ROOT TF1 function template into a ROOT TH1D histogram.

    Args:
        func (ROOT.TF1): The ROOT TF1 function template to convert.
        bins (int): The number of bins for the histogram.
        low (float): The lower edge of the histogram.
        high (float): The upper edge of the histogram.
        integral (bool, optional): If True, compute the integral values instead
            of evaluating the function at bin centers.

    Returns:
        ROOT.TH1D: The converted histogram.

    Raises:
        None.

    Example:
        # Create a TF1 template
        template_func = ROOT.TF1("template", "x^2", 0, 10)

        # Convert the template to a histogram
        histogram = TemplateToHistogram(template_func, 100, 0, 10)

    Notes:
        - The function template is evaluated at each bin center and assigned to
          the corresponding bin in the histogram.
        - If the `integral` flag is True, the integral value within each bin is
          computed and assigned to the bin content instead.
        - If the function evaluation or integral computation results in NaN
          (Not a Number), the bin content is set to 0.
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

  Returns a linearized version as a python function
  '''
  def histfunc(x):
    return hist.GetBinContent(hist.FindBin(x[0])) 
  return histfunc

def MultiplyWithPolyToTF1(func, degree, range_low=0, range_high=50, poly=0, parameters=None):
    '''
    Takes a python function

    Returns a TF1 object representing the input function times a polynomial, and the returns associated python function as well

    when poly=0 (default) use regular polynomials (1, x^2, x^3, etc)
    when poly=1 use Chebyshev polynomials of the first kind
    when poly=2 use Chebyshev polynomials of the second kind
    '''
    if poly == 0:
        if degree == 0:
            def polynomial(x, p):
                return p[0]
        elif degree == 1:
            def polynomial(x, p):
                return p[0] + p[1]*x[0]
        elif degree == 2:
            def polynomial(x, p):
                return p[0] + p[1]*x[0] + p[2]*x[0]**2
        elif degree == 3:
            def polynomial(x, p):
                return p[0] + p[1]*x[0] + p[2]*x[0]**2 + p[3]*x[0]**3
        elif degree == 4:
            def polynomial(x, p):
                return p[0] + p[1]*x[0] + p[2]*x[0]**2 + p[3]*x[0]**3 \
                      + p[4]*x[0]**4
        elif degree == 5:
            def polynomial(x, p):
                return p[0] + p[1]*x[0] + p[2]*x[0]**2 + p[3]*x[0]**3 \
                        + p[4]*x[0]**4 + p[5]*x[0]**5
        elif degree == 6:
            def polynomial(x, p):
                return p[0] + p[1]*x[0] + p[2]*x[0]**2 + p[3]*x[0]**3 \
                        + p[4]*x[0]**4 + p[5]*x[0]**5 + p[6]*x[0]**6
        elif degree == 7:
            def polynomial(x, p):
                return p[0] + p[1]*x[0] + p[2]*x[0]**2 + p[3]*x[0]**3 \
                        + p[4]*x[0]**4 + p[5]*x[0]**5 + p[6]*x[0]**6 + p[7]*x[0]**7
        elif degree == 8:
            def polynomial(x, p):
                return p[0] + p[1]*x[0] + p[2]*x[0]**2 + p[3]*x[0]**3 \
                        + p[4]*x[0]**4 + p[5]*x[0]**5 + p[6]*x[0]**6 + p[7]*x[0]**7 \
                        + p[8]*x[0]**8
        else:
            raise ValueError('Got degree {} for poly type {}. Not Implemented!'.format(degree, poly))

    elif poly == 1:
        if degree == 0:
            def polynomial(x, p):
                return (p[0])
        elif degree == 1:
            def polynomial(x, p):
                X = x[0]
                return (p[0] + p[1]*cheb_fn(X, 1, 1))
        elif degree == 2:
            def polynomial(x, p):
                X = x[0]
                return (p[0] + p[1]*cheb_fn(X, 1, 1) + p[2]*cheb_fn(X, 2, 1))
        elif degree == 3:
            def polynomial(x, p):
                X = x[0]
                return (p[0] + p[1]*cheb_fn(X, 1, 1) + p[2]*cheb_fn(X, 2, 1)
                            + p[3]*cheb_fn(X, 3, 1))
        elif degree == 4:
            def polynomial(x, p):
                X = x[0]
                return (p[0] + p[1]*cheb_fn(X, 1, 1) + p[2]*cheb_fn(X, 2, 1)
                            + p[3]*cheb_fn(X, 3, 1) + p[4]*cheb_fn(X, 4, 1))
        elif degree == 5:
            def polynomial(x, p):
                X = x[0]
                return (p[0] + p[1]*cheb_fn(X, 1, 1) + p[2]*cheb_fn(X, 2, 1)
                            + p[3]*cheb_fn(X, 3, 1) + p[4]*cheb_fn(X, 4, 1)
                            + p[5]*cheb_fn(X, 5, 1))
        elif degree == 6:
            def polynomial(x, p):
                X = x[0]
                return (p[0] + p[1]*cheb_fn(X, 1, 1) + p[2]*cheb_fn(X, 2, 1)
                            + p[3]*cheb_fn(X, 3, 1) + p[4]*cheb_fn(X, 4, 1)
                            + p[5]*cheb_fn(X, 5, 1) + p[6]*cheb_fn(X, 6, 1))
        elif degree == 7:
            def polynomial(x, p):
                X = x[0]
                return (p[0] + p[1]*cheb_fn(X, 1, 1) + p[2]*cheb_fn(X, 2, 1)
                            + p[3]*cheb_fn(X, 3, 1) + p[4]*cheb_fn(X, 4, 1)
                            + p[5]*cheb_fn(X, 5, 1) + p[6]*cheb_fn(X, 6, 1)
                            + p[7]*cheb_fn(X, 7, 1))
        elif degree == 8:
            def polynomial(x, p):
                X = x[0]
                return (p[0] + p[1]*cheb_fn(X, 1, 1) + p[2]*cheb_fn(X, 2, 1)
                            + p[3]*cheb_fn(X, 3, 1) + p[4]*cheb_fn(X, 4, 1)
                            + p[5]*cheb_fn(X, 5, 1) + p[6]*cheb_fn(X, 6, 1)
                            + p[7]*cheb_fn(X, 7, 1) + p[8]*cheb_fn(X, 8, 1))
        else:
            raise ValueError('Got degree {} for poly type {}. Not Implemented!'.format(degree, poly))

    elif poly == 2:
        if degree == 0:
            def polynomial(x, p):
                return (p[0])
        elif degree == 1:
            def polynomial(x, p):
                X = x[0]
                return (p[0] + p[1]*cheb_fn(X, 1, 2))
        elif degree == 2:
            def polynomial(x, p):
                X = x[0]
                return (p[0] + p[1]*cheb_fn(X, 1, 2) + p[2]*cheb_fn(X, 2, 2))
        elif degree == 3:
            def polynomial(x, p):
                X = x[0]
                return (p[0] + p[1]*cheb_fn(X, 1, 2) + p[2]*cheb_fn(X, 2, 2)
                            + p[3]*cheb_fn(X, 3, 2))
        elif degree == 4:
            def polynomial(x, p):
                X = x[0]
                return (p[0] + p[1]*cheb_fn(X, 1, 2) + p[2]*cheb_fn(X, 2, 2)
                            + p[3]*cheb_fn(X, 3, 2) + p[4]*cheb_fn(X, 4, 2))
        elif degree == 5:
            def polynomial(x, p):
                X = x[0]
                return (p[0] + p[1]*cheb_fn(X, 1, 2) + p[2]*cheb_fn(X, 2, 2)
                            + p[3]*cheb_fn(X, 3, 2) + p[4]*cheb_fn(X, 4, 2)
                            + p[5]*cheb_fn(X, 5, 2))
        elif degree == 6:
            def polynomial(x, p):
                X = x[0]
                return (p[0] + p[1]*cheb_fn(X, 1, 2) + p[2]*cheb_fn(X, 2, 2)
                            + p[3]*cheb_fn(X, 3, 2) + p[4]*cheb_fn(X, 4, 2)
                            + p[5]*cheb_fn(X, 5, 2) + p[6]*cheb_fn(X, 6, 2))
        elif degree == 7:
            def polynomial(x, p):
                X = x[0]
                return (p[0] + p[1]*cheb_fn(X, 1, 2) + p[2]*cheb_fn(X, 2, 2)
                            + p[3]*cheb_fn(X, 3, 2) + p[4]*cheb_fn(X, 4, 2)
                            + p[5]*cheb_fn(X, 5, 2) + p[6]*cheb_fn(X, 6, 2)
                            + p[7]*cheb_fn(X, 7, 2))
        elif degree == 8:
            def polynomial(x, p):
                X = x[0]
                return (p[0] + p[1]*cheb_fn(X, 1, 2) + p[2]*cheb_fn(X, 2, 2)
                            + p[3]*cheb_fn(X, 3, 2) + p[4]*cheb_fn(X, 4, 2)
                            + p[5]*cheb_fn(X, 5, 2) + p[6]*cheb_fn(X, 6, 2)
                            + p[7]*cheb_fn(X, 7, 2) + p[8]*cheb_fn(X, 8, 2))
        else:
            raise ValueError('Got degree {} for poly type {}. Not Implemented!'.format(degree, poly))

    elif poly == 3:
        if degree == 0:
            def polynomial(x, p):
                X = x[0]
                return p[0]
        elif degree == 1:
            def polynomial(x, p):
                X = x[0]
                return p[0]*bern_fn(X, 1, 0) + p[1]*bern_fn(X, 1, 1)
        elif degree == 2:
            def polynomial(x, p):
                X = x[0]
                return p[0]*bern_fn(X, 2, 0) + p[1]*bern_fn(X, 2, 1) \
                       + p[2]*bern_fn(X, 2, 2)
        elif degree == 3:
            def polynomial(x, p):
                X = x[0]
                return p[0]*bern_fn(X, 3, 0) + p[1]*bern_fn(X, 3, 1) \
                       + p[2]*bern_fn(X, 3, 2) + p[3]*bern_fn(X, 3, 3)
        elif degree == 4:
            def polynomial(x, p):
                X = x[0]
                return p[0]*bern_fn(X, 4, 0) + p[1]*bern_fn(X, 4, 1) \
                       + p[2]*bern_fn(X, 4, 2) + p[3]*bern_fn(X, 4, 3) + p[4]*bern_fn(X, 4, 4)
        elif degree == 5:
            def polynomial(x, p):
                X = x[0]
                return p[0]*bern_fn(X, 5, 0) + p[1]*bern_fn(X, 5, 1) \
                       + p[2]*bern_fn(X, 5, 2) + p[3]*bern_fn(X, 5, 3) + p[4]*bern_fn(X, 5, 4) \
                       + p[5]*bern_fn(X, 5, 5)
        elif degree == 6:
            def polynomial(x, p):
                X = x[0]
                return p[0]*bern_fn(X, 6, 0) + p[1]*bern_fn(X, 6, 1) \
                       + p[2]*bern_fn(X, 6, 2) + p[3]*bern_fn(X, 6, 3) + p[4]*bern_fn(X, 6, 4) \
                       + p[5]*bern_fn(X, 6, 5) + p[6]*bern_fn(X, 6, 6)
        elif degree == 7:
            def polynomial(x, p):
                X = x[0]
                return p[0]*bern_fn(X, 7, 0) + p[1]*bern_fn(X, 7, 1) \
                       + p[2]*bern_fn(X, 7, 2) + p[3]*bern_fn(X, 7, 3) + p[4]*bern_fn(X, 7, 4) \
                       + p[5]*bern_fn(X, 7, 5) + p[6]*bern_fn(X, 7, 6) + p[7]*bern_fn(X, 7, 7)
        elif degree == 8:
            def polynomial(x, p):
                X = x[0]
                return p[0]*bern_fn(X, 8, 0) + p[1]*bern_fn(X, 8, 1) \
                       + p[2]*bern_fn(X, 8, 2) + p[3]*bern_fn(X, 8, 3) + p[4]*bern_fn(X, 8, 4) \
                       + p[5]*bern_fn(X, 8, 5) + p[6]*bern_fn(X, 8, 6) + p[7]*bern_fn(X, 8, 7) \
                       + p[8]*bern_fn(X, 8, 8)
        else:
            raise ValueError('Got degree {} for poly type {}. Not Implemented!'.format(degree, poly))

    else:
        raise ValueError('Got degree {} for poly type {}. Not Implemented!'.format(degree, poly))

    def func_after_mult(x, p):
        #return func(x) * polynomial(x, p)
        val = func(x) * polynomial(x, p)
        #if val < 10**-30: val = 10**-30
        return val

    globals()[getname('func')] = polynomial
    globals()[getname('func')] = func_after_mult

    if poly == 0 or poly == 1 or poly == 2:
        num_param = degree + 1
    if poly == 3 and degree == 0:
        num_param = 1
    if poly == 3 and degree > 0:
        num_param = degree + 1

    tf1 = ROOT.TF1(getname(), func_after_mult, range_low, range_high, num_param)

    # Set parameter names
    if num_param>=0: tf1.SetParName(0, 'Constant') if poly==0 else tf1.SetParName(0, 'Zero')
    if num_param>=1: tf1.SetParName(1, 'Linear') if poly==0 else tf1.SetParName(1, 'One')
    if num_param>=2: tf1.SetParName(2, 'Quadratic') if poly==0 else tf1.SetParName(2, 'Two')
    if num_param>=3: tf1.SetParName(3, 'Cubic') if poly==0 else tf1.SetParName(3, 'Three')
    if num_param>=4: tf1.SetParName(4, 'Quartic') if poly==0 else tf1.SetParName(4, 'Four')
    if num_param>=5: tf1.SetParName(5, 'Quintic') if poly==0 else tf1.SetParName(5, 'Five')
    if num_param>=6: tf1.SetParName(6, 'Sextic') if poly==0 else tf1.SetParName(6, 'Six')
    if num_param>=7: tf1.SetParName(7, 'Septic') if poly==0 else tf1.SetParName(7, 'Seven')
    if num_param>=8: tf1.SetParName(8, 'Octic') if poly==0 else tf1.SetParName(8, 'Eight')
    if num_param>=9: tf1.SetParName(9, 'Nonic') if poly==0 else tf1.SetParName(9, 'Nine')

    # Bernstein polynomials need positive coefficients
    if poly == 3:
        for i in range(num_param): tf1.SetParLimits(i, 0, 1e9)

    # For constant function require positive constant
    if degree == 0:
        tf1.SetParLimits(0, 0, 1e9)

    # Set parameter initial guesses
    if not parameters:
        if poly == 0 or poly == 1 or poly == 2:
            for i in range(num_param): tf1.SetParameter(i, 0.1**i)
        if poly == 3:
            for i in range(num_param): 
                tf1.SetParameter(i, 0.1)
                if num_param == 5:
                    test_params = [0.0415947, 0.00761743, 0.00217256, 0.00706772, 0.0125019, 0.101876]
                    tf1.SetParameter(i, test_params[i])
    else:
        tf1.SetParameters(*parameters)

    return tf1, func_after_mult, polynomial

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

  elif function == 'landau' and N == 2:
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
      if y2 == 0: land2 = (y1) * ROOT.TMath.Landau(x[0], mpv2, sigma2)
      else: land2 = (y1/y2) * ROOT.TMath.Landau(x[0], mpv2, sigma2)
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
    tf1.SetParLimits(3, range_low, range_high*2)
    tf1.SetParLimits(4, 0, hist.GetMean()*5)
    tf1.SetParLimits(5, 0, 1e5)

  elif function == 'exp' and N == 1:
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

  elif function == 'exp' and N == 2:
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

  elif function == 'exp' and N == 3:
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
    tf1.SetParLimits(2, 0, range_high)
    tf1.SetParLimits(3, -10, 0)
    tf1.SetParLimits(4, 0, range_high/2.0)
    tf1.SetParLimits(5, -10, 0)
  
  elif function == 'full' and N == 11: 
    def python_func(x, p):
        norm1 = p[0]
        mpv1 = p[1]
        sigma1 = p[2]
        C1 = p[3]
        bound1 = p[4]

        if bound1 < 0: bound1 = 0

        land1 = norm1 * ROOT.TMath.Landau(x[0], mpv1, sigma1)
      
        y11=norm1*ROOT.TMath.Landau(bound1, mpv1, sigma1)
        y12=ROOT.TMath.Exp(C1*bound1)
        if y12 == 0: exp1 = ROOT.TMath.Exp(C1*x[0])*y11
        else: exp1 = ROOT.TMath.Exp(C1*x[0])*y11/y12
         
        if x[0] < bound1: return land1
        else: return exp1

    NPAR = 5 
    if not initial_guesses:
        nEntries = hist.GetEntries()
        mean = hist.GetMean()
        initial_guesses = [
            nEntries, mean, 0.2, -3,
            mean]
    if not len(initial_guesses) == NPAR:
        raise AssertionError('Length of initial guesses list must be '+str(NPAR)+"!")
    tf1 = ROOT.TF1(getname('func'), python_func, range_low, range_high, NPAR)
    tf1.SetParNames("Constant1","MPV1","Sigma1","C1","Boundary1")
    for i, guess in enumerate(initial_guesses): tf1.SetParameter(i, guess)
    tf1.SetParLimits(3, -10, 0)
    tf1.SetParLimits(4, 0, 25)

  elif function == 'full' and N == 12: 
    def python_func(x, p):
        norm1 = p[0]
        mpv1 = p[1]
        sigma1 = p[2]
        C1 = p[3]
        C2 = p[4]
        bound1 = p[5]
        b12 = p[6]

        if bound1 < 0: bound1 = 0
        bound2 = bound1 + b12

        land1 = norm1 * ROOT.TMath.Landau(x[0], mpv1, sigma1)
      
        y11=norm1*ROOT.TMath.Landau(bound1, mpv1, sigma1)
        y12=ROOT.TMath.Exp(C1*bound1)
        if y12 == 0: exp1 = ROOT.TMath.Exp(C1*x[0])*y11
        else: exp1 = ROOT.TMath.Exp(C1*x[0])*y11/y12
         
        y21=ROOT.TMath.Exp(C1*bound2)*y11/y12
        y22=ROOT.TMath.Exp(C2*bound2)
        exp2=ROOT.TMath.Exp(C2*x[0])*y21/y22

        if x[0] < bound1: return land1
        elif x[0] < bound2: return exp1
        else: return exp2

    NPAR = 7
    if not initial_guesses:
        nEntries = hist.GetEntries()
        mean = hist.GetMean()
        initial_guesses = [
            nEntries, mean, 0.2, -3, -1,
            mean, mean/2]
    if not len(initial_guesses) == NPAR:
        raise AssertionError('Length of initial guesses list must be '+str(NPAR)+"!")
    tf1 = ROOT.TF1(getname('func'), python_func, range_low, range_high, NPAR)
    tf1.SetParNames("Constant1","MPV1","Sigma1","C1","C2","Boundary1","BoundDiff12")
    for i, guess in enumerate(initial_guesses): tf1.SetParameter(i, guess)
    tf1.SetParLimits(3, -10, 0)
    tf1.SetParLimits(4, -10, 0)
    tf1.SetParLimits(5, 0, 25)
    tf1.SetParLimits(6, 0.2, 7)
  
  elif function == 'full' and N == 13: 
    def python_func(x, p):
        norm1 = p[0]
        mpv1 = p[1]
        sigma1 = p[2]
        C1 = p[3]
        C2 = p[4]
        C3 = p[5]
        bound1 = p[6]
        b12 = p[7]
        b23 = p[8]

        if bound1 < 0: bound1 = 0
        bound2 = bound1 + b12
        bound3 = bound2 + b23

        land1 = norm1 * ROOT.TMath.Landau(x[0], mpv1, sigma1)
      
        y11=norm1*ROOT.TMath.Landau(bound1, mpv1, sigma1)
        y12=ROOT.TMath.Exp(C1*bound1)
        if y12 == 0: exp1 = ROOT.TMath.Exp(C1*x[0])*y11
        else: exp1 = ROOT.TMath.Exp(C1*x[0])*y11/y12
         
        y21=ROOT.TMath.Exp(C1*bound2)*y11/y12
        y22=ROOT.TMath.Exp(C2*bound2)
        exp2=ROOT.TMath.Exp(C2*x[0])*y21/y22

        y31=ROOT.TMath.Exp(C2*bound3)*y21/y22
        y32=ROOT.TMath.Exp(C3*bound3)
        exp3=ROOT.TMath.Exp(C3*x[0])*y31/y32
        
        if x[0] < bound1: return land1
        elif x[0] < bound2: return exp1
        elif x[0] < bound3: return exp2
        else: return exp3

    NPAR = 9
    if not initial_guesses:
        nEntries = hist.GetEntries()
        mean = hist.GetMean()
        initial_guesses = [
            nEntries, mean, 0.2, -3, -1, -0.5,
            mean, mean/2, mean/2]
    if not len(initial_guesses) == NPAR:
        raise AssertionError('Length of initial guesses list must be '+str(NPAR)+"!")
    tf1 = ROOT.TF1(getname('func'), python_func, range_low, range_high, NPAR)
    tf1.SetParNames("Constant1","MPV1","Sigma1","C1","C2","C3","Boundary1","BoundDiff12","BoundDiff23")
    for i, guess in enumerate(initial_guesses): tf1.SetParameter(i, guess)
    tf1.SetParLimits(3, -10, 0)
    tf1.SetParLimits(4, -10, 0)
    tf1.SetParLimits(5, -10, 0)
    tf1.SetParLimits(6, 0, 25)
    tf1.SetParLimits(7, 0.2, 7)
    tf1.SetParLimits(8, 0.2, 7)

  elif function == 'full' and N == 21:
    def python_func(x, p):
        norm1 = p[0]
        mpv1 = p[1]
        sigma1 = p[2]
        norm2 = p[3]
        mpv2 = p[4]
        sigma2 = p[5]
        C1 = p[6]
        bound1 = p[7]
        b12 = p[8]

        if bound1 < 0: bound1 = 0
        bound2 = bound1 + b12
        
        land1 = norm1 * ROOT.TMath.Landau(x[0], mpv1, sigma1)
      
        y11=norm1*ROOT.TMath.Landau(bound1, mpv1, sigma1)
        y12=norm2*ROOT.TMath.Landau(bound1, mpv2, sigma2)
        if y12 == 0: land2 = norm2*ROOT.TMath.Landau(x[0], mpv2, sigma2)*y11
        else: land2 = norm2*ROOT.TMath.Landau(x[0], mpv2, sigma2)*y11/y12
         
        if y12 == 0: y21=norm2*ROOT.TMath.Landau(bound2, mpv2, sigma2)*y11
        else: y21=norm2*ROOT.TMath.Landau(bound2, mpv2, sigma2)*y11/y12
        y22=ROOT.TMath.Exp(C1*bound2)
        exp1=ROOT.TMath.Exp(C1*x[0])*y21/y22

        if x[0] < bound1: return land1
        elif x[0] < bound2: return land2
        else: return exp1
    NPAR = 9
    if not initial_guesses:
        nEntries = hist.GetEntries()
        mean = hist.GetMean()
        initial_guesses = [
            nEntries, mean, 0.5, nEntries, mean, 0.5, -3, mean, mean/2]
    if not len(initial_guesses) == NPAR:
        raise AssertionError('Length of initial guesses list must be '+str(NPAR)+"!")
    tf1 = ROOT.TF1(getname('func'), python_func, range_low, range_high, NPAR)
    tf1.SetParNames("Constant1","MPV1","Sigma1","Constant2","MPV2","Sigma2","C1","Boundary1","BoundDiff12")
    for i, guess in enumerate(initial_guesses): tf1.SetParameter(i, guess)
    tf1.SetParLimits(6, -10, 0)
    tf1.SetParLimits(7, 0, 25)
    tf1.SetParLimits(8, 0.2, 7)
  
  elif function == 'full' and N == 24:
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
        if y12 == 0: land2 = norm2*ROOT.TMath.Landau(x[0], mpv2, sigma2)*y11
        else: land2 = norm2*ROOT.TMath.Landau(x[0], mpv2, sigma2)*y11/y12
         
        if y12 == 0: y21=norm2*ROOT.TMath.Landau(bound2, mpv2, sigma2)*y11
        else: y21=norm2*ROOT.TMath.Landau(bound2, mpv2, sigma2)*y11/y12
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

  else:
    raise ValueError('fit_hist(): Invalid type of fit: got function='+str(function)+', N='+str(N))

  globals()[getname('func')] = python_func
  fit_string = '0SL'
  if integral: fit_string += 'I'
  fit_result = hist.Fit(tf1, fit_string, "", range_low, range_high)
  return tf1, fit_result

def RSS(func, hist, bound=-1, error=0, integral=False, chi2=False, cutoff=None):
  '''
  helper for ftest()
  '''
  rss = 0
  by_bin = []
  for i in range(hist.GetNbinsX()):
    if hist.GetBinLowEdge(i+1) < bound: continue
    hist_val = hist.GetBinContent(i+1)
    if integral:
      fit_val = (func.Integral(hist.GetBinLowEdge(i+1), hist.GetBinLowEdge(i+1) + hist.GetBinWidth(i+1)))/hist.GetBinWidth(i+1)
    else:
      fit_val = func.Eval(hist.GetBinCenter(i+1))
    if cutoff and fit_val < cutoff: continue
    if error == 0:
        sr = (hist_val - fit_val)**2
    else:
        adj = abs(fit_val * error)
        sr = (max(abs(hist_val - fit_val) - adj, 0))**2
    if chi2:
        if hist_val > fit_val: hist_error = hist.GetBinErrorLow(i+1)
        else: hist_error = hist.GetBinErrorUp(i+1)
        if hist_val==0: hist_error = hist.GetBinErrorUp(i+1)
        sr_sigma = sr/(hist_error)**2
    else:
        sr_sigma = sr
    rss += sr_sigma
    by_bin.append((i, hist.GetBinLowEdge(i+1), hist_val, fit_val, sr_sigma, sr, hist_error))
  return rss, by_bin

def ExtractPolyFromTightFit(fitfunc, range_low=0, range_high=50, poly=0, debug=False):
    nparam = fitfunc.GetNpar()
    if poly == 0 or poly == 1 or poly == 2: poly_degree = nparam - 1
    elif poly == 3 and nparam == 1: poly_degree = 0
    elif poly == 3 and nparam > 1: poly_degree = nparam - 1
    _, _, poly_func = MultiplyWithPolyToTF1(lambda x : x[0], poly_degree, poly=poly)
    extracted_poly = ROOT.TF1(getname('func'), poly_func, range_low, range_high, nparam)
    for n in range(nparam):
        extracted_poly.SetParameter(n, fitfunc.GetParameter(n))
        if debug: print(n, fitfunc.GetParameter(n))
    return extracted_poly

def hist_to_numpy(hist):
    array = np.zeros(int(hist.GetNbinsX()))
    for j in range(hist.GetNbinsX()):
        array[j] = hist.GetBinContent(j+1)
    return array 

def count_nonzero_bins(hist, bound=-1):
  '''
  Helper for ftest()
  '''
  count = 0
  for i in range(hist.GetNbinsX()):
    if hist.GetBinLowEdge(i+1) < bound: continue
    if not hist.GetBinContent(i+1) == 0: count += 1
  return count

def lookup_ftest_target(dof1, dof2, sig):
    '''
    Helper for ftest()
    '''
    if 'scipy' in sys.modules:
        return f.ppf(1-sig, dof1, dof2)
    else:
        return 2

def ftest(hist, func2, fit2, func1, fit1, sig=0.1, integral=False):
    '''
    F-test comparing fit2 (more parameters) vs fit1 (less parameters)
    '''
    rss1, _ = RSS(func1, hist, integral=integral, chi2=True)
    rss2, _ = RSS(func2, hist, integral=integral, chi2=True)
    p1 = func1.GetNpar()
    p2 = func2.GetNpar()
    n = count_nonzero_bins(hist)
    F = ((rss1 - rss2)/(p2 - p1)) / (rss2/(n - p2))
    dof1, dof2 = p2 - p1, n - p2
    target = lookup_ftest_target(dof1, dof2, sig)
    decision = F > target
    return decision, F, target, rss1, rss2, dof1, dof2


def lookup_fit_guesses(control_region, eta_reg, pt_bin):
    # Default new method used; if there are less than ENTRIES_CUTOFF entries, the "old" method is used instead of the bulk fitting
    old_method = False
    landau_guess = None
    exp_guess = None
    guesses = None
    
    # Default guess for landau and exponential, all fits are adjusted around these
    nLandau = 1
    nExp = 3
    
    # Rarely any "bulk" region, i.e. few bins with more than ENTRIES_CUTOFF entries
    if control_region == "iso_sym" or control_region == "iso_asym": old_method = True

    # Iso-Sym Guesses
    if control_region == "iso_sym":
        if eta_reg == "barrel":
            if pt_bin == 20:
                nExp -= 1
                guesses = [1313, 0.7032, 0.1068, -5.47, -10, 1.014, 0.6314]
            if pt_bin == 40:
                nExp -= 2
                guesses = [1500, 0.9313, 0.1856, -3.766, 1.024]
            if pt_bin == 60:
                nExp -= 2
                guesses = [785.9, 1.117, 0.2567, -3.498, 1.184]
            if pt_bin == 80:
                nExp -= 1
                guesses = [1648, 1.205, 0.2847, -1.647, -3.407, 1.275, 0.321]
            if pt_bin == 100:
                nLandau += 1 
                nExp -= 2 
                #guesses = [2000, 0.5, 0.1, -0.5, -1.3, 1.1, 0.6]
                guesses = [2000, 0.5, 0.05, 1000, 1, 0.6, -1.7, 1.1, 0.6]
            if pt_bin == 140:
                nExp -= 1
                guesses = [1887, 1.314, 0.3242, -1.142, -2.368, 1.674, 0.2515]
            if pt_bin == 180:
                guesses = [2380, 1.361, 0.3352, -1.189, -2.023, -0.05, 1.775, 0.4189, 1.5]
            if pt_bin == 220:
                old_method = False
                landau_guess = [5458, 1.394, 0.3427]
                exp_guess = [5458, -1.012, 1.678, -1.784, 0.5775, -1]
            if pt_bin == 300:
                guesses = [4131, 1.451, 0.3552, -0.0001324, -1.054, -1.36, 1.369, 0.3256, 1.168]
            if pt_bin == 380:
                guesses = [1559, 1.54, 0.3837, -0.08148, -0.882, -1.123, 1.591, 0.2366, 0.337]

        elif eta_reg == "endcap":
            if pt_bin == 20:
                nExp -= 1
                guesses = [2198, 0.6681, 0.09261, -5.133, -10, 0.854, 0.721]
            if pt_bin == 40:
                nExp -= 2 
                guesses = [1000, 0.8865, 0.1719, -3.938, 1]
            if pt_bin == 60:
                nExp -= 1 
                guesses = [961.9, 1.039, 0.219, -0.06295, -3.632, 0.925, 0.2]
            if pt_bin == 80:
                nExp -= 1
                guesses = [1590, 1.129, 0.2523, -1.411, -3.112, 1.251, 0.2239]
            if pt_bin == 100:
                nExp -= 1
                guesses = [1480, 1.228, 0.2869, -1.143, -2.683, 1.387, 0.2375]
            if pt_bin == 140:
                nExp -= 1
                guesses = [1424, 1.296, 0.3049, -0.458, -2.143, 1.348, 0.4153]
            if pt_bin == 180:
                nExp -= 1
                guesses = [1559, 1.436, 0.3578, -0.8864, -2.058, 1.525, 0.6868]
            if pt_bin == 220:
                guesses = [3040, 1.432, 0.3386, -0.1704, -1.026, -1.618, 1.405, 0.4994, 0.6806]
            if pt_bin == 300:
                guesses = [1596, 1.523, 0.3731, -0.07532, -0.9641, -1.38, 1.405, 0.4691, 1.187]
            if pt_bin == 380:
                nExp -= 1
                guesses = [475.8, 1.65, 0.4243, -0.5718, -1.198, 1.825, 0.8539]
            if pt_bin == 460:
                guesses = [195.5, 1.389, 0.2903, -0.0001153, -0.7601, -0.9159, 1.305, 0.8198, 0.2123]

    # Iso-Asym Guesses
    if control_region == "iso_asym":
        if eta_reg == "barrel":
            if pt_bin == 20:
                nExp -= 1
                guesses = [246.5, 0.8312, 0.1213, -4.7, -10, 0.875, 0.4371]
            if pt_bin == 40:
                nExp -= 2
                guesses = [422.5, 0.9891, 0.1768, -4.297, 1.116]
            if pt_bin == 60:
                nExp -= 1
                guesses = [224.5, 1.061, 0.193, -2.913, -5.832, 1.141, 1.184]
            if pt_bin == 80:
                nExp -= 1
                guesses = [532.3, 1.19, 0.2404, -0.1541, -2.411, 1.075, 0.2]
            if pt_bin == 100:
                nExp -= 1
                guesses = [585.1, 1.328, 0.2896, -1.456, -2.425, 1.362, 0.2857]
            if pt_bin == 140:
                guesses = [772.1, 1.426, 0.3213, -0.09985, -2.09, -10, 1.379, 0.3031, 2.843]
            if pt_bin == 180:
                nExp -= 1
                guesses = [1236, 1.713, 0.4136, -1.181, -1.847, 1.674, 0.3422]
            if pt_bin == 220:
                old_method = False
                landau_guess = [3956, 1.6, 0.3]
                exp_guess = [2000, -0.7, 2.2, -1.487, 4, -0.8996]
            if pt_bin == 300:
                guesses = [1.706e+04, 1.529, 0.3766, -0.01528, -0.7213, -1.017, 4.5, 4, 2] 
            if pt_bin == 380:
                guesses = [1096, 2.116, 0.5494, -0.4469, -1.293, -0.8872, 2.428, 0.4568, 2.857]
            if pt_bin == 460:
                guesses = [793.4, 2.157, 0.5492, -0.01715, -0.9391, -0.5238, 2.047, 0.8956, 5.37]

        elif eta_reg == "endcap":
            if pt_bin == 20:
                nExp -= 1
                guesses = [246.5, 0.8312, 0.1213, -4.7, -10, 0.875, 0.4371]
            if pt_bin == 40:
                nExp -= 2
                guesses = [546.6, 0.8803, 0.1328, -4.339, 1.089]
            if pt_bin == 60:
                nExp -= 2
                guesses = [346.3, 1.128, 0.2239, -3.488, 0.9861]
            if pt_bin == 80:
                guesses = [787.1, 1.41, 0.3193, -0.9958, -2.985, -2.563, 1.049, 0.3013, 0.9245]
            if pt_bin == 100:
                nExp -= 1
                guesses = [678.8, 1.319, 0.2715, -2.841e-10, -2.461, 1.08, 0.2679]
            if pt_bin == 140:
                nExp -= 1
                guesses = [713.7, 1.365, 0.2823, -2.554e-11, -2.152, 1.299, 0.3368]
            if pt_bin == 180:
                guesses = [1041, 1.725, 0.4141, -1.863, -1.999, -0.9901, 1.807, 0.3725, 2.195]
            if pt_bin == 220:
                guesses = [2679, 2.01, 0.505, -1.64, -1.099, -1.108, 2.146, 2.514, 2.416]
            if pt_bin == 300:
                guesses = [1606, 1.895, 0.4265, -0.7, -1.483, -1, 2, 0.5, 2.5]
            if pt_bin == 380:
                guesses = [500, 2.301, 0.6017, -1.264, -1.592, -0.8502, 2.804, 0.971, 0.9496]

    # NonIso-Sym Guesses
    if control_region == "noniso_sym":
        if eta_reg == "barrel":
            if pt_bin == 20:
                landau_guess = [3.491e+06, 0.8838, 0.145] 
                exp_guess = [8.593e+05, -4.278, 1.595, -7.263, 0.3, -8]
            if pt_bin == 40:
                landau_guess = [1.705e+04, 1.33, 0.205]
                exp_guess = [3.636e+05, -2.809, 2.176, -4.8, 0.6, -1.807]
            if pt_bin == 60:
                nExp -= 1
                landau_guess = [2.212e+04, 0.6783, 0.09052]
                exp_guess = [4.467e+05, -2.789, 2.974, -5.139]
            if pt_bin == 80:
                landau_guess = [1.826e+04, 0.6239, 0.07162]
                exp_guess = [1.506e+04, -0.9609, 2.4, -2.24, 2.2, -0.4279]
            if pt_bin == 100:
                landau_guess = [1.705e+04, 0.6496, 0.0804] 
                exp_guess = [1.755e+04, -0.92, 3.125, -1.706, 2, -1]
            if pt_bin == 140:
                landau_guess = [8972, 0.8638, 0.1507]
                exp_guess = [9448, -1.041, 3.2, -1.423, 3, -3]
            if pt_bin == 180:
                landau_guess = [1.183e+04, 0.7977, 0.1293]
                exp_guess = [4.739e+04, -1.172, 3.5, -3.605, 5.5, -3.677]
            if pt_bin == 220:
                landau_guess = [8717, 0.8303, 0.1381]
                exp_guess = [9.701e+04, -1.193, 6.173, -0.6467, 6.452, -2.239]
            if pt_bin == 300:
                landau_guess = [1145, 1.149, 0.2453]
                exp_guess = [3.14e+04, -1.023, 6.75, -0.5922, 1.25, -1]
            if pt_bin == 380:
                old_method = True
                guesses = [4848, 1.63, 0.4089, -3.553e-13, -0.6753, -0.8161, 0.5, 1.404, 0.5423]

        elif eta_reg == "endcap":
            if pt_bin == 20:
                landau_guess = [1.705e+04, 1.108, 0.1949] 
                exp_guess = [4182, -0.9553, 1.4, -5.89, 0.5, -1]
            if pt_bin == 40:
                landau_guess = [1.123e+05, 0.875, 0.1524]
                exp_guess = [4.549e+04, -2.09, 1.915, -3.572, 0.6484, -6.328]
            if pt_bin == 60:
                landau_guess = [1.021e+04, 0.7237, 0.107]
                exp_guess = [5970, -0.9294, 1.702, -2.435, 1.262, -5.012]
            if pt_bin == 80:
                landau_guess = [1.826e+04, 0.6239, 0.07162]
                exp_guess = [1.506e+04, -0.9609, 2.4, -2.24, 2.2, -0.4279]
            if pt_bin == 100:
                landau_guess = [1.705e+04, 0.6496, 0.2]
                exp_guess = [1.755e+04, -0.92, 2.4, -1.706, 2.6, -1]
            if pt_bin == 140:
                landau_guess = [1.761e+04, 0.7322, 0.1067]
                exp_guess = [9448, -1.041, 3.2, -1.423, 3, -3]
            if pt_bin == 180:
                landau_guess = [1.183e+04, 0.7977, 0.1293]
                exp_guess = [4.739e+04, -1.172, 3.5, -3.605, 5.5, -3.677]
            if pt_bin == 220:
                landau_guess = [7937, 1.269, 0.2811]
                exp_guess = [4.457e+04, -1.4, 2.6, -1.9, 1.4, -1]
            if pt_bin == 300:
                old_method = True
                guesses = [1.706e+04, 1.529, 0.3766, -0.01528, -0.7213, -1.017, 0.5, 0.8, 0.5]  
            if pt_bin == 380:
                old_method = True
                guesses = [890.4, 1.723, 0.4403, -1.441, -0.9773, -0.2518, 2.5, 2.5, 5]

    # NonIso_Asym Guesses   
    if control_region == "noniso_asym":
        if eta_reg == "barrel":
            if pt_bin == 20:
                landau_guess = [9.069e+04, 0.9345, 0.1463]
                exp_guess = [4.834e+04, -3.378, 1.206, -8.486, 0.7373, -4.555]
            if pt_bin == 40:
                nExp -= 1
                landau_guess = [2.459e+04, 0.8719, 0.1261]
                exp_guess = [9328, -1.29, 1.599, -5.368]
            if pt_bin == 60:
                landau_guess = [1.021e+04, 0.7237, 0.107]
                exp_guess = [5970, -0.9294, 1.702, -2.435, 1.262, -5.012]
            if pt_bin == 80:
                landau_guess = [1.474e+04, 1.084, 0.1992]
                exp_guess = [7428, -0.8801, 2.325, -2.899, 1.044, -2.326]
            if pt_bin == 100:
                landau_guess = [1.705e+04, 0.6496, 0.0804] 
                exp_guess = [1.986e+04, -0.8552, 2.8, -2.188, 1.2, -1.68]
            if pt_bin == 140:
                landau_guess = [8211, 1.16, 0.2245] 
                exp_guess = [3.881e+04, -1.409, 3.225, -1.699, 2.5, -1.078]
            if pt_bin == 180:
                landau_guess = [1.183e+04, 0.7977, 0.1293]
                exp_guess = [4.739e+04, -1.172, 3.5, -3.605, 5.5, -3.677]
            if pt_bin == 220:
                landau_guess = [8442, 1.149, 0.2453]
                exp_guess = [3.14e+04, -1.023, 5, -0.5922, 3, -2]
            if pt_bin == 300:
                old_method = True
                nExp -= 1
                guesses = [3000, 2.157, 0.5678, -1.051, -0.58, 2, 6]
            if pt_bin == 380:
                old_method = True
                guesses = [1430, 2.107, 0.5122, -0.04034, -1.359, -0.5, 3, 3, 2]

        elif eta_reg == "endcap":
            if pt_bin == 20:
                nExp -= 1
                landau_guess = [5000, 2, 0.3]
                exp_guess = [4.834e+04, -3.378, 1.206, -7.486] 
            if pt_bin == 40:
                landau_guess = [1e+04, 0.9233, 0.1527]
                exp_guess = [3937, -0.4774, 0.9104, -2.799, 0.7293, -4.712]
            if pt_bin == 60:
                landau_guess = [1.021e+04, 0.7237, 0.107]
                exp_guess = [5970, -0.9294, 1.702, -2.435, 1.262, -5.012]
            if pt_bin == 80:
                landau_guess = [1.826e+04, 0.6239, 0.07162]
                exp_guess = [1.506e+04, -0.9609, 2.4, -2.24, 2.2, -0.4279]
            if pt_bin == 100:
                landau_guess = [1.705e+04, 0.6496, 0.2]
                exp_guess = [1.755e+04, -0.92, 2.4, -1.706, 2.6, -1]
            if pt_bin == 140:
                landau_guess = [6416, 1.397, 0.3019]
                exp_guess = [8019, -1.101, 1.875, -1.559, 3.2, -1.032]
            if pt_bin == 180:
                landau_guess = [7031, 1.004, 0.1965]
                exp_guess = [1.742e+04, -1.198, 3, -3.296, 3, -1.366]
            if pt_bin == 220:
                nExp -= 1
                old_method = True
                guesses = [4657, 2.046, 0.5207, -1.225, -0.7, 2.5, 3.5]
            if pt_bin == 300:
                old_method = True
                nExp -= 2 
                guesses = [1470, 1.687, 0.3483, -2, 2]
            if pt_bin == 380:
                old_method = True
                guesses = [1430, 2.107, 0.5122, -0.04034, -1.359, -0.5, 3, 3, 2]
    
    fit_init = {}
    fit_init['old_method'] = old_method
    fit_init['guesses'] = guesses
    fit_init['nLandau'] = nLandau
    fit_init['landau_guess'] = landau_guess
    fit_init['nExp'] = nExp
    fit_init['exp_guess'] = exp_guess
    return fit_init

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
    input()

    print('turn histo into python function')
    func = HistogramToFunction(hist)

    print('multiply with polynomial')
    func_with_poly, _, _ = MultiplyWithPolyToTF1(func, 2)
    func_with_poly.Draw()
    input()

    print('make second histogram')
    hist_target = ROOT.TH1D(getname(), 'target', bins, low, high)
    hist_target.FillRandom('f1', 100000)
    hist_target.Draw()
    input()

    print('fit second histogram with function')
    hist_target.Fit(func_with_poly, "L")
    hist_target.Draw()
    input()

    print('extract parameters of poly from fit')
    f = func_with_poly
    n = f.GetNpar()
    for i in range(n):
        print(f.GetParName(i), f.GetParameter(i))

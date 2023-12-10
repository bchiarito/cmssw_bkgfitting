import ROOT
import sys
import math

try:
    from scipy.stats import f 
except ImportError:
    pass

BERN_UPPER_RANGE = 25
NAME_COUNT = 0 # global counter

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
                #if i == 0: tf1.FixParameter(i, 0.150094)
                #if i == 1: tf1.FixParameter(i, 2.67687e-7)
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

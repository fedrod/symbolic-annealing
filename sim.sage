# https://omega0.xyz/sim.sage

def show_latex(expr):
    ''' %display latex
    doesn't work anymore.
    Use this function instead.
    USAGE:  show_latex(solve(x^2-2,x))'''
    from IPython.display import display, Latex
    display(Latex(r'$\displaystyle ' + latex(expr) + '$'))

ltx = show_latex

#
# Simplification by simulated annealing.
# Annealing by Richard Wagner.
# Simplification by Carlos Rodriguez.
# 
#
#
# Python module for simulated annealing - anneal.py - v1.0 - 2 Sep 2009
#
# Copyright (c) 2009, Richard J. Wagner <wagnerr@umich.edu>
#
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

"""
This module performs simulated annealing to find a state of a system that
minimizes its energy.

An example program demonstrates simulated annealing with a traveling
salesman problem to find the shortest route to visit the twenty largest
cities in the United States.

Notes:
    Matt Perry 6/24/12
        Changed to slicing lists instead of deepcopy-ing them.
        e.g. state = prevState[:] instead of state = deepcopy(prevState)
        Huge performance enhancement (~5-10x faster)
        Should be identical behavior if the items in the state list are immutable.
        (immutable objects include integers and strings so should be safe)
"""

# How to optimize a system with simulated annealing:
#
# 1) Define a format for describing the state of the system.
#
# 2) Define a function to calculate the energy of a state.
#
# 3) Define a function to make a random change to a state.
#
# 4) Choose a maximum temperature, minimum temperature, and number of steps.
#
# 5) Set the annealer to work with your state and functions.
#
# 6) Study the variation in energy with temperature and duration to find a
# productive annealing schedule.
#
# Or,
#
# 4) Run the automatic annealer which will attempt to choose reasonable values
# for maximum and minimum temperatures and then anneal for the allotted time.

import copy, math, sys, time
try:
    from numpy import random
except ImportError:
    import random

def round_figures(x, n):
    """Returns x rounded to n significant figures."""
    return round(x, int(n - math.ceil(math.log10(abs(x)))))

def time_string(seconds):
    """Returns time in seconds as a string formatted HHHH:MM:SS."""
    s = int(round(seconds))  # round to nearest second
    h, s = divmod(s, 3600)   # get hours and remainder
    m, s = divmod(s, 60)     # split remainder into minutes and seconds
    return '%4i:%02i:%02i' % (h, m, s)

class Annealer:
    """Performs simulated annealing by calling functions to calculate
    energy and make moves on a state.  The temperature schedule for
    annealing may be provided manually or estimated automatically.
    """
    def __init__(self, energy, move):
        self.energy = energy  # function to calculate energy of a state
        self.move = move      # function to make a random change to a state

    def anneal(self, state, Tmax, Tmin, steps, updates=0, on_accept=None):
        """Minimizes the energy of a system by simulated annealing.

        Keyword arguments:
        state -- an initial arrangement of the system
        Tmax -- maximum temperature (in units of energy)
        Tmin -- minimum temperature (must be greater than zero)
        steps -- the number of steps requested
        updates -- the number of updates to print during annealing

        Returns the best state and energy found."""

        step = 0
        start = time.time()

        def update(T, E, acceptance, improvement):
            """Prints the current temperature, energy, acceptance rate,
            improvement rate, elapsed time, and remaining time.

            The acceptance rate indicates the percentage of moves since the last
            update that were accepted by the Metropolis algorithm.  It includes
            moves that decreased the energy, moves that left the energy
            unchanged, and moves that increased the energy yet were reached by
            thermal excitation.

            The improvement rate indicates the percentage of moves since the
            last update that strictly decreased the energy.  At high
            temperatures it will include both moves that improved the overall
            state and moves that simply undid previously accepted moves that
            increased the energy by thermal excititation.  At low temperatures
            it will tend toward zero as the moves that can decrease the energy
            are exhausted and moves that would increase the energy are no longer
            thermally accessible."""

            elapsed = time.time() - start
            if step == 0:
                print(' Temperature        Energy    Accept   Improve     Elapsed   Remaining')
                print( '%12.2f  %12.2f                      %s            ' %                     (T, E, time_string(elapsed) ))
            else:
                remain = ( steps - step ) * ( elapsed / step )
                print( '%12.2f  %12.2f  %7.2f%%  %7.2f%%  %s  %s' % (T, E, 100.0*acceptance, 100.0*improvement,
                        time_string(elapsed), time_string(remain)))

        # Precompute factor for exponential cooling from Tmax to Tmin
        if Tmin <= 0.0:
            print( 'Exponential cooling requires a minimum temperature greater than zero.')
            sys.exit()
        Tfactor = -math.log( float(Tmax) / Tmin )

        # Note initial state
        T = Tmax
        E = self.energy(state)
        #prevState = copy.deepcopy(state)
        prevState = state[:]
        prevEnergy = E
        #bestState = copy.deepcopy(state)
        bestState = state[:]
        bestEnergy = E
        if on_accept is not None:
            on_accept(state, E)
        trials, accepts, improves = 0, 0, 0
        if updates > 0:
            updateWavelength = float(steps) / updates
            update(T, E, None, None)

        # Attempt moves to new states
        while step < steps:
            step += 1
            T = Tmax * math.exp( Tfactor * step / steps )
            self.move(state)
            E = self.energy(state)
            dE = E - prevEnergy
            trials += 1
            if dE > 0.0 and math.exp(-dE/T) < random.random():
                # Restore previous state
                #state = copy.deepcopy(prevState)
                state = prevState[:]
                E = prevEnergy
            else:
                # Accept new state and compare to best state
                accepts += 1
                if dE < 0.0:
                    improves += 1
                #prevState = copy.deepcopy(state)
                prevState = state[:]
                prevEnergy = E
                if on_accept is not None:
                    on_accept(state, E)
                if E < bestEnergy:
                    #bestState = copy.deepcopy(state)
                    bestState = state[:]
                    bestEnergy = E
            if updates > 1:
                if step // updateWavelength > (step-1) // updateWavelength:
                    update(T, E, float(accepts)/trials, float(improves)/trials)
                    trials, accepts, improves = 0, 0, 0

        # Return best state and energy
        return bestState, bestEnergy

    def auto(self, state, minutes, steps=2000):
        """Minimizes the energy of a system by simulated annealing with
        automatic selection of the temperature schedule.

        Keyword arguments:
        state -- an initial arrangement of the system
        minutes -- time to spend annealing (after exploring temperatures)
        steps -- number of steps to spend on each stage of exploration

        Returns the best state and energy found."""

        def run(state, T, steps):
            """Anneals a system at constant temperature and returns the state,
            energy, rate of acceptance, and rate of improvement."""
            E = self.energy(state)
            #prevState = copy.deepcopy(state)
            prevState = state[:]
            prevEnergy = E
            accepts, improves = 0, 0
            for step in range(steps):
                self.move(state)
                E = self.energy(state)
                dE = E - prevEnergy
                if dE > 0.0 and math.exp(-dE/T) < random.random():
                    #state = copy.deepcopy(prevState)
                    state = prevState[:]
                    E = prevEnergy
                else:
                    accepts += 1
                    if dE < 0.0:
                        improves += 1
                    #prevState = copy.deepcopy(state)
                    prevState = state[:]
                    prevEnergy = E
            return state, E, float(accepts)/steps, float(improves)/steps

        step = 0
        start = time.time()

        print( 'Attempting automatic simulated anneal...')

        # Find an initial guess for temperature
        T = 0.0
        E = self.energy(state)
        while T == 0.0:
            step += 1
            self.move(state)
            T = abs( self.energy(state) - E )

        print( 'Exploring temperature landscape:')
        print( ' Temperature        Energy    Accept   Improve     Elapsed')
        def update(T, E, acceptance, improvement):
            """Prints the current temperature, energy, acceptance rate,
            improvement rate, and elapsed time."""
            elapsed = time.time() - start
            print( '%12.2f  %12.2f  %7.2f%%  %7.2f%%  %s' % (T, E, 100.0*acceptance, 100.0*improvement, time_string(elapsed)))

        # Search for Tmax - a temperature that gives 98% acceptance
        state, E, acceptance, improvement = run(state, T, steps)
        step += steps
        while acceptance > 0.98:
            T = round_figures(T/1.5, 2)
            state, E, acceptance, improvement = run(state, T, steps)
            step += steps
            update(T, E, acceptance, improvement)
        while acceptance < 0.98:
            T = round_figures(T*1.5, 2)
            state, E, acceptance, improvement = run(state, T, steps)
            step += steps
            update(T, E, acceptance, improvement)
        Tmax = T

        # Search for Tmin - a temperature that gives 0% improvement
        while improvement > 0.0:
            T = round_figures(T/1.5, 2)
            state, E, acceptance, improvement = run(state, T, steps)
            step += steps
            update(T, E, acceptance, improvement)
        Tmin = T

        # Calculate anneal duration
        elapsed = time.time() - start
        duration = round_figures(int(60.0 * minutes * step / elapsed), 2)

        # MP: Don't perform anneal, just return params
        #return self.anneal(state, Tmax, Tmin, duration, 20)
        return {'tmax': Tmax, 'tmin': Tmin, 'steps': duration}


############################ routines for simplification follow this line


from mpmath import mp, pslq
import re


def clean_latex(ex, specials, multiplier):
    r"""
    Convert a SageMath expression to its LaTeX string and clean it by:
      - Replacing each special substring (or command) with '#' repeated 'multiplier' times.
        The specials list can contain items like r'\sqrt' or plain substrings like 'C' or 'a23'.
      - Removing all whitespace.
      - Replacing all other LaTeX commands with a single '#'.
      - Removing the LaTeX spacing command r'\,'.
      - Removing curly braces and commas.
    
    Parameters:
      ex         : SageMath expression.
      specials   : List of special strings to be replaced. For example:
                   [r'\sqrt'] or [r'\sqrt', 'C', 'a23'].
                   Defaults to [r'\sqrt'].
      multiplier : The number of '#' characters used for the special replacement.
    
    Returns:
      The cleaned LaTeX string.
    """
    if specials is None:
        specials = [r'\sqrt']
    
    sx = str(latex(ex))
    
    # Replace each special substring exactly as provided with '#' * multiplier.
    for special in specials:
        pattern = re.escape(special) + r'\b'
        sx = re.sub(pattern, '#' * multiplier, sx)
    
    # Remove all whitespace.
    sx = re.sub(r'\s+', '', sx)
    
    # Replace all remaining LaTeX commands with a single '#'
    sx = re.sub(r'\\[a-zA-Z]+', '#', sx)
    
    # Remove the spacing command r'\,'
    sx = re.sub(r'\\,', '', sx)
    
    # Remove curly braces and commas.
    sx = re.sub(r'[{},]+', '', sx)
    
    return sx

def lenl(ex, specials=None, multiplier=2):
    r"""
    Return the length of the cleaned LaTeX string representation of the expression,
    using the provided specials list and multiplier.
    """
    return len(clean_latex(ex, specials=specials, multiplier=multiplier))



def remove_redundancy(ex,tol=1e-15):
    try:
        if ex.operator() != (var('x')+1).operator():
            raise ValueError("The expression is not a sum\n")

        ijs = Combinations(range(len(ex)),2)
        exl = ex.operands()
        for i,j in ijs:
            if exl[i] == exl[j] or abs(mp.mpc(exl[i]-exl[j])) < tol:
                if len(exl[i]) < len(exl[j]): exl[i] = 2*exl[i]
                else: exl[i] = 2*exl[j]
                exl[j] = 0
            if exl[i] == -exl[j] or abs(mp.mpc(exl[i]+exl[j])) < tol:
                exl[i],exl[j] = 0,0
            if exl[i] == exl[j].conjugate() or abs(mp.mpc(exl[i]-exl[j].conjugate())) < tol:
                if len(real_part(exl[i])) < len(real_part(exl[j])): exl[i] = 2*real_part(exl[i])
                else: exl[i] = 2*real_part(exl[j])
                exl[j] = 0
        return sum(exl)
    except:
        return ex

rred = remove_redundancy

def decomplex(ex,tol=1e-15):
    try:
        if ex.operator() != (var('x')+1).operator():
            raise ValueError("The expression is not a sum\n")
        exl = ex.op[:]
        for i in xrange(len(exl)):
            if abs(mp.mpf(imag_part(exl[i]))) < tol: exl[i] = real_part(exl[i])

        E = [i for i in xrange(len(exl)) if not exl[i].n().is_real()]
        es = [exl[i] for i in E]
        S = sum(es)
        if abs(mp.mpf(imag_part(S))) < tol:
            if abs(mp.mpf(real_part(S))) < tol:
                for i in E: exl[i] = 0
            else:
                for i in E: exl[i] = real_part(exl[i])
        return(sum(exl))
    except:
        return ex

decom = decomplex

def iargmax(a):
    m = max(a)
    ii = [i for i,j in enumerate(a) if j == m]
    return ii[0]

def pslqc(ex,tol=1e-15):
    try:
        if ex.operator() != (var('x')+1).operator():
            raise ValueError("The expression is not a sum\n")
        exl = ex.op[:]
        n = len(exl)
        exr = map(real_part,exl)
        exc = map(imag_part,exl)
        lenths = map(len,exl)
        cj = pslq(exr)
        ex_new = ex
        try:
            if cj != [] and vector(cj).dot_product(vector(exc)) == 0:
                A = [i for i in xrange(n) if cj[i] != 0]
                lA = [lenths[i] for i in A]
                im = iargmax(lA)
                # remove im term
                ex_new = sum([(1-QQ(cj[j])/cj[im])*exl[j] for j in xrange(n)])
        except TypeError:
            return ex
        return ex_new
    except:
        return ex

###### NEW STUFF

# useful trig substitutions
w0 = SR.wild(0)
w1 = SR.wild(1)

S2PW = [cos(w1) + cos(w0) == 2*cos(1/2*w1 + 1/2*w0)*cos(-1/2*w1 + 1/2*w0),
 -cos(w1) + cos(w0) == -2*sin(1/2*w1 + 1/2*w0)*sin(-1/2*w1 + 1/2*w0),
 sin(w1) + sin(w0) == 2*cos(-1/2*w1 + 1/2*w0)*sin(1/2*w1 + 1/2*w0),
 -sin(w1) + sin(w0) == 2*cos(1/2*w1 + 1/2*w0)*sin(-1/2*w1 + 1/2*w0)]

P2SW = [cos(w1)*sin(w0) == 1/2*sin(w1 + w0) + 1/2*sin(-w1 + w0),
 cos(w0)*sin(w1) == 1/2*sin(w1 + w0) - 1/2*sin(-w1 + w0),
 cos(w1)*cos(w0) == 1/2*cos(w1 + w0) + 1/2*cos(-w1 + w0),
 sin(w1)*sin(w0) == -1/2*cos(w1 + w0) + 1/2*cos(-w1 + w0)]

SQSW = [cos(w0)^2 == 1/2*cos(2*w0) + 1/2,
 sin(w0)^2 == -1/2*cos(2*w0) + 1/2,
 tan(w0)^2 == -(cos(2*w0) - 1)/(cos(2*w0) + 1),
 cos(w0) + 1 == 2*cos(1/2*w0)^2,
 -cos(2*w0) + 1 == 2*sin(w0)^2]

S2E = {cos(w0): 1/2*exp(i*w0)+1/2*exp(-i*w0), sin(w1): 1/2/i*exp(i*w1)-1/2/i*exp(-i*w1)}


TL_available = ['TL_sim','TL_sumsim','TL_log','TL_simplify','TL_trig','TL_trig_subs','TL_pslqc']

TL_sim = ['_psf',"_el","_sl","_fs",'_Snj','_Sdj','_cb','_xp','_sr','_er','_es','_fa','_pf','_mrs','_pmrs']

TL_sumsim = ['_rr','_dc','_psf']

TL_pslqc = ['_ps'] # may take forever!

TL_log = ['_collect_common_factors','_combine','_distribute','_exp_simplify', '_expand', '_expand_log', '_expand_rational', '_expand_sum', '_expand_trig', '_factor', '_factorial_simplify','_gamma_normalize', '_log_expand','_log_simplify']

TL_simplify = ['_simplify', '_simplify_factorial', '_simplify_full', '_simplify_hypergeometric', '_simplify_log', '_simplify_radical_exp', '_simplify_rational', '_simplify_real', '_simplify_rectform', '_simplify_trig']

TL_trig =['_s2e','_trig_expand', '_trig_reduce', '_trig_simplify']
TL_trig_subs = ['_s2pw','_p2sw','_sqsw','_s2e']

TL_ALL = TL_sim+TL_sumsim+TL_log+TL_simplify+TL_trig
TL_default = TL_sim

def _s2pw(ex):
    return ex.subs(S2PW).full_simplify()
def _p2sw(ex):
    return ex.subs(P2SW).full_simplify()
def _sqsw(ex):
    return ex.subs(SQSW).full_simplify()

def s2e(ex):
    return Sim(ex.subs(S2E).canonicalize_radical().real_part())

def _s2e(ex):
    """
    transform sines and cosines to complex exponentials
    """
    try:
        return ex.subs(S2E).canonicalize_radical().real_part().full_simplify()
    except:
        return ex


def _simplify(ex):
    try:
        return ex.simplify()
    except:
        return ex


def _simplify_exp(ex):
    try:
        return ex.simplify_exp()
    except:
        return ex


def _simplify_factorial(ex):
    try:
        return ex.simplify_factorial()
    except:
        return ex


def _simplify_full(ex):
    try:
        return ex.simplify_full()
    except:
        return ex


def _simplify_hypergeometric(ex):
    try:
        return ex.simplify_hypergeometric()
    except:
        return ex


def _simplify_log(ex):
    try:
        return ex.simplify_log()
    except:
        return ex


def _simplify_radical_exp(ex):
    try:
        return ex.canonicalize_radical()
    except:
        return ex


def _simplify_rational(ex):
    try:
        return ex.simplify_rational()
    except:
        return ex


def _simplify_real(ex):
    try:
        return ex.simplify_real()
    except:
        return ex


def _simplify_rectform(ex):
    try:
        return ex.simplify_rectform()
    except:
        return ex


def _simplify_trig(ex):
    try:
        return ex.simplify_trig()
    except:
        return ex


def _trig_expand(ex):
    try:
        return ex.trig_expand()
    except:
        return ex


def _trig_reduce(ex):
    try:
        return ex.trig_reduce()
    except:
        return ex


def _trig_simplify(ex):
    try:
        return ex.trig_simplify()
    except:
        return ex

def _collect_common_factors(ex):
    try:
        return ex.collect_common_factors()
    except:
        return ex


def _combine(ex):
    try:
        return ex.combine()
    except:
        return ex


def _distribute(ex):
    try:
        return ex.distribute()
    except:
        return ex


def _exp_simplify(ex):
    try:
        return ex.canonicalize_radical()
    except:
        return ex


def _expand(ex):
    try:
        return ex.expand()
    except:
        return ex


def _expand_log(ex):
    try:
        return ex.expand_log()
    except:
        return ex


def _expand_rational(ex):
    try:
        return ex.expand_rational()
    except:
        return ex


def _expand_sum(ex):
    try:
        return ex.expand_sum()
    except:
        return ex


def _expand_trig(ex):
    try:
        return ex.expand_trig()
    except:
        return ex


def _factor(ex):
    try:
        return ex.factor()
    except:
        return ex


def _factorial_simplify(ex):
    try:
        return ex.factorial_simplify()
    except:
        return ex


def _gamma_normalize(ex):
    try:
        return ex.gamma_normalize()
    except:
        return ex


def _log_expand(ex):
    try:
        return ex.log_expand()
    except:
        return ex


def _log_simplify(ex):
    try:
        return ex.log_simplify()
    except:
        return ex


def _rr(ex):
    try:
        e = remove_redundancy(ex)
        return e
    except ValueError:
        return ex

def _dc(ex):
    try:
        e = decomplex(ex)
        return e
    except ValueError:
        return ex

def _ps(ex):
    try:
        e = pslqc(ex)
        return e
    except ValueError:
        return ex

def _el(ex): return ex.expand_log()
def _sl(ex): return ex.simplify_log()
def _fs(ex): return ex.full_simplify()

def _Snj(ex):
    ''' factor the numerator of the jth(random) subexpression'''
    try:
        exl = ex.op[:]
        j = random.randint(len(exl))
        so = ex.operator()
        exl[j] = exl[j].numerator().factor()/exl[j].denominator()
        return af2l(so,exl)
    except:
        return ex

def _Sdj(ex):
    ''' factor the denominator of the jth(random) subexpression'''
    try:
        exl = ex.op[:]
        j = random.randint(len(exl))
        so = ex.operator()
        exl[j] = exl[j].numerator()/exl[j].denominator().factor()
        return af2l(so,exl)
    except:
        return ex


def _cb(ex):
    try:
        return ex.combine()
    except:
        return ex

def _sr(ex):
    try:
        return ex.canonicalize_radical()
    except AttributeError:
        return ex

def _er(ex):
    try:
        return ex.expand_rational()
    except:
        return ex

def _es(ex):
    try:
        return ex.canonicalize_radical()
    except:
        return ex

def _fa(ex):
    try:
        return ex.factor()
    except:
        return ex

def _fand(e):
    try:
        return (e.numerator().factor()/e.denominator().factor()).simplify()
    except:
        return e

def _xp(ex):
    try:
        return ex.expand()
    except:
        return ex

# partial factorization
def _pf(e):
    try:
        n = len(e.op[:])
        k = random.randint(1,n,size=1)[0]
        p = random.choice(n,k,replace=False)
        notp = [i for i in range(n) if not i in p]
        q = (e.op[k] for k in p)
        q1 = (e.op[k] for k in notp)
        e0 = e.operator()(*q)
        return e.operator()(factor(e0),*q1)
    except:
        return e

# partial rootscontract

def _mrs(ex):
    try:
        e1 = ex.maxima_methods().rootscontract()
        return e1.operator()(factor(e1.op[0]),*e1.op[1:])
    except:
        return ex


def _pmrs(e):
    try:
        n = len(e.op[:])
        k = random.randint(0,n,size=1)[0]
        p = random.choice(n,k,replace=False)
        notp = [i for i in range(n) if not i in p]
        q = (e.op[k] for k in p)
        q1 = (e.op[k] for k in notp)
        e0 = e.operator()(*q)
        return e.operator()(_mrs(e0),*q1)
    except:
        return e

# partial sum simplification
def _psf(e):
    try:
        if e.operator() != (var('x')+1).operator():
            raise ValueError("The expression is not a sum\n")
        n = len(e.op[:])
        k = random.randint(0,n,size=1)[0]
        p = random.choice(n,k,replace=False)
        notp = [i for i in range(n) if not i in p]
        q = (e.op[k] for k in p)
        q1 = (e.op[k] for k in notp)
        return sum(q).factor() + sum(q1)
    except:
        return e

#### END NEW STUFF


def subsj(ex):
    try:
        exl = ex.op[:]
        j = random.randint(len(exl))
        return ex.subs({ex.op[j]:ex.op[j].combine()})
    except TypeError:
        return(ex)

def sdj(ex):
    try:
        exl = ex.op[:]
        j = random.randint(len(exl))
        return ex.subs({ex.op[j].denominator():ex.op[j].denominator().factor()})
    except:
        return ex

def snj(ex):
    try:
        exl = ex.op[:]
        j = random.randint(len(exl))
        return ex.subs({ex.op[j].numerator():ex.op[j].numerator().factor()})
    except:
        return ex

def af2l(f,L):
    """
    apply binary function f to list L
    """
    if len(L) == 2: return f(L[0],L[1])
    else:
        return f(L[0],af2l(f,L[1:]))


## empty TL
def _id(e): return e
TL_none = ['_id']

def inverse_dic(d): return {d[k]:k for k in d.keys()}

def discover_substitutions(ex):
    '''
    Automatically discover (variable + constant) binomial patterns in an expression
    and generate a substitution dictionary for use with Sim() and subSim().

    Walks the AST recursively, finding two-operand sums where one operand is a
    pure variable and the other is a constant. Generates short variable names
    (e.g. alpha - 1 -> a1, alpha + 1 -> ap1) and includes negated forms.

    Example:
        sage: dic1 = discover_substitutions(xp)
        sage: Sim(x1, subs_dic=dic1)
    '''
    subs_dic = {}
    seen_names = {}
    add_op = (var('x') + 1).operator()

    def find_binomials(e):
        if not hasattr(e, 'operator'): return
        op = e.operator()
        if op is None: return

        if op == add_op:
            ops = e.operands()
            vars_list = [o for o in ops if len(o.variables()) == 1 and o == o.variables()[0]]
            const_list = [o for o in ops if len(o.variables()) == 0]

            if len(vars_list) == 1 and len(const_list) == 1 and len(ops) == 2:
                v = vars_list[0]
                c = const_list[0]

                v_str = str(v)
                prefix = v_str[0] if v_str else 'v'

                if c < 0:
                    base_name = f"{prefix}{str(abs(c)).replace('/', '_')}"
                else:
                    base_name = f"{prefix}p{str(abs(c)).replace('/', '_')}"

                new_var_name = base_name
                counter = 1
                while new_var_name in seen_names and seen_names[new_var_name] != e:
                    new_var_name = f"{base_name}_{counter}"
                    counter += 1

                if e not in subs_dic:
                    from sage.symbolic.ring import SR
                    new_var = SR.symbol(new_var_name)

                    subs_dic[e] = new_var
                    subs_dic[-e] = -new_var
                    seen_names[new_var_name] = e

        for operand in e.operands():
            find_binomials(operand)

    find_binomials(ex)
    return subs_dic


def Sim(ex,Tmax=100,Tmin=1,steps=100,updates=0, energy = False,ops=(),subs_dic={},_TL=TL_default,specials=None,multiplier=2,on_accept=None):
    '''

    Simplification of the expression ex by Simulated Annealing. (see "Annealer").

    energy = length of the penalized LaTeX representation of ex.

    a list of 'specials' is used to incentivize or decentivize
    combining expressions (see help(clean_latex)).

    _TL is the list of available transformations for the Annealer. (see "TL_available").

    If either ops or subs_dic is non empty then subexpression simplification and
    substitution is applied by calling subSim which
    attempts to simplify e.op[ops] with 'Sim'.

    Example:

    sage: var('a b c d')
    
    sage: e0 = sqrt((a^3-b^3)/(a-b)+a*b)
    
    sage: Sim(e0)
    
    sage: ex = integral(1/((a*x+b)^2*(c*x+c)^2),x)
    
    sage: ex_sim = Sim(ex,energy= True)
    
    sage: show([ex_sim,lenl(ex,specials,multiplier)])
    
    sage: var('a a1 b x y')
    
    sage: t =  -2*a*b*x*y^2+3*b^2*x*y^2+(a^2*x-2*a*b*x+b^2*x)*y^2
    
    sage: Sim(t,ops=(2,0),subs_dic={(a-b):a1},updates=5)

    '''
    def ex_energy(state): return lenl(state[0],specials,multiplier)

    def ex_move(state):
        j = random.randint(len(_TL))
        state[0] = eval(_TL[j]+"(state[0])")

    ann = Annealer(ex_energy,ex_move)

    if ops != () or subs_dic != {}:
        return subSim(ex,Tmax,Tmin,steps,updates,energy,ops,subs_dic)
    state = [ex]
    state,e = ann.anneal(state,Tmax,Tmin,steps,updates,on_accept=on_accept)
    if energy: return state[0],e
    else:
        return state[0]

def subSim(e,Tmax=100,Tmin=1,steps=100,updates=0, energy = False,ops=(),subs_dic={}):
    '''
    subSim(e,ops=(),subs_dic={})
    subexpression simplification and substitution.
    It attempts to simplify e.op[ops] with 'Sim'.
    After simplification it applies the substitutions given by
    the subs_dic dictionary.
    It returns the original expression with the simplified subexpression and
    the substitutions.
    Example:

    sage: var('a a1 b x y')
    sage: t =  -2*a*b*x*y^2+3*b^2*x*y^2+(a^2*x-2*a*b*x+b^2*x)*y^2
    sage: subSim(t,ops=(2,0),subs_dic={(a-b):a1})
    '''

    if ops == (): return Sim(e,Tmax,Tmin,steps,updates,energy).subs(subs_dic)
    newpart = Sim(e.op[ops],Tmax,Tmin,steps,updates,energy)
    return e.subs({e.op[ops]:newpart}).subs(subs_dic)
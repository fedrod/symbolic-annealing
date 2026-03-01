# Sim: Symbolic Simplification by Annealing

A SageMath library that automatically simplifies symbolic mathematical expressions using **simulated annealing**. Instead of relying on a fixed sequence of algebraic rules, Sim treats simplification as an optimization problem — randomly applying transformations and keeping the ones that produce shorter, more elegant formulas.

## The Idea

Traditional computer algebra systems apply simplification rules in a predetermined order. This often misses the shortest form because the "best" sequence of transformations depends on the expression itself.

Sim reframes simplification as energy minimization:

- **State**: a symbolic expression
- **Energy**: the length of its LaTeX representation (optionally penalizing certain constructs like `\sqrt`)
- **Moves**: randomly chosen algebraic transformations — factor, expand, trig identities, log rules, radical simplification, and more
- **Cooling schedule**: the Metropolis criterion from simulated annealing accepts moves that shorten the expression, but also probabilistically accepts longer forms at high temperatures to escape local minima

This approach is inspired by [Kirkpatrick, Gelatt & Vecchi (1983)](https://doi.org/10.1126/science.220.4598.671), "Optimization by Simulated Annealing," which showed that the Metropolis algorithm could solve hard combinatorial optimization problems by analogy with statistical mechanics.

## Results

The paper (`Sim.pdf`) demonstrates the method on the **peak of a beta distribution's triangular envelope** — a formula involving parameters &alpha; and &beta; that starts as a 564-character LaTeX expression and is reduced to just 66 characters:

$$xp = \frac{1}{a_1^{-a_1} a_2^{a_2} b_1^{b_1} b_2^{-b_2} + 1}$$

where $a_1 = \alpha - 1$, $a_2 = \alpha - 2$, $b_1 = \beta - 1$, $b_2 = \beta - 2$.

This is an **8.5x reduction** that no single built-in SageMath simplification command can achieve.

## Quick Start

```python
# In a SageMath session:
load("sim.sage")

# Basic simplification
var('a b c d')
e0 = sqrt((a^3 - b^3)/(a - b) + a*b)
Sim(e0)

# With energy output (shows expression + its LaTeX length)
ex = integral(1/((a*x + b)^2 * (c*x + c)^2), x)
ex_sim, length = Sim(ex, energy=True)

# Subexpression simplification with substitution
var('a a1 b x y')
t = -2*a*b*x*y^2 + 3*b^2*x*y^2 + (a^2*x - 2*a*b*x + b^2*x)*y^2
Sim(t, ops=(2, 0), subs_dic={(a - b): a1})
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ex` | — | The symbolic expression to simplify |
| `Tmax` | 100 | Maximum annealing temperature |
| `Tmin` | 1 | Minimum annealing temperature |
| `steps` | 100 | Number of annealing steps |
| `updates` | 0 | Print progress every N updates |
| `energy` | False | If True, return `(expression, length)` tuple |
| `ops` | () | Tuple of operand indices for subexpression targeting |
| `subs_dic` | {} | Substitution dictionary applied after simplification |
| `_TL` | TL_default | Transformation list to use (see below) |
| `specials` | None | LaTeX constructs to penalize (e.g. `[r'\sqrt']`) |
| `multiplier` | 2 | Penalty weight for special constructs |

## Transformation Lists

The annealer draws random moves from configurable transformation lists:

| List | Contents |
|------|----------|
| `TL_sim` | Core simplifiers: full_simplify, factor, expand, combine, partial factorization, rootscontract, etc. |
| `TL_sumsim` | Sum-level: redundancy removal, decomplexification, partial sum factoring |
| `TL_simplify` | SageMath's built-in simplify variants (factorial, hypergeometric, log, radical, trig, ...) |
| `TL_log` | Logarithmic and exponential transforms |
| `TL_trig` | Trigonometric simplification and expansion |
| `TL_trig_subs` | Trig identity substitutions (sum-to-product, product-to-sum, half-angle) |
| `TL_pslqc` | PSLQ integer relation detection for eliminating redundant terms |

## Key References

- **Kirkpatrick, S., Gelatt, C. D., & Vecchi, M. P.** (1983). [Optimization by Simulated Annealing](https://doi.org/10.1126/science.220.4598.671). *Science*, 220(4598), 671–680.
- **Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., & Teller, E.** (1953). [Equation of State Calculations by Fast Computing Machines](https://doi.org/10.1063/1.1699114). *Journal of Chemical Physics*, 21(6), 1087–1092.
- **Wagner, R. J.** (2009). Python module for simulated annealing (`anneal.py`), included and adapted in this project.

## Requirements

- [SageMath](https://www.sagemath.org/) (includes Python, SymPy, Maxima, and mpmath)

## Author

**Carlos Rodriguez** — Department of Mathematics and Statistics, University at Albany (SUNY)

See `Sim.pdf` for the full paper with worked examples.

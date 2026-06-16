# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## About

A Python number theory library written as a learning exercise alongside *Elements of Number Theory* (Ireland & Rosen), *Number Fields* (Marcus), and *Lectures on Modular Forms* (Gunning).

## Commands

```bash
# Run all tests
python -m pytest

# Run a single test file
python -m pytest test/test_utilities.py

# Run a single test
python -m pytest test/test_utilities.py::UtilitiesTest::test_gcd

# Check code style
pycodestyle

# Build Sphinx docs
cd doc && make html
```

Modules import each other directly (e.g. `from utilities import isprime`), so scripts and the test suite must be run from the project root.

## Dependencies

`requirements.txt` lists `numpy` and `mpmath`. Several modules also use `sympy` (cyclotomic fields, Fibonacci) and `matplotlib` (modular forms, modular groups). `pynauty` is optional — `graph_algorithms.py` stubs it out gracefully if not installed.

## Code style

`pycodestyle` is enforced with these ignores (see `tox.ini`): E501 (line length), E741 (ambiguous variable names like `l`, `O`, `I`), and several continuation/whitespace rules. Single-character doubled variable names (`ii`, `kk`, etc.) are flagged by `test/test_variable_names.py`.

## Architecture

The codebase is a flat collection of Python modules. There is no package hierarchy — `__init__.py` is empty.

**Foundational utilities** (imported widely):
- `utilities.py` — gcd, bezout, isprime, factorize/factorize2, modpow, modinv, primitive_root, trial_division, Miller-Rabin witness test, and general helpers (prod, digits, timeit)
- `roots.py` — exact integer nth-root algorithms: `sqrtInt`, `cbrtInt`, `nrtInt`, `issq`
- `digits.py` — digit manipulation

**Primality and factoring**:
- `primality_tests.py` — Miller-Rabin and related probabilistic tests
- `factoring.py` — Fermat factoring, Pollard rho, quadratic sieve; imports from `utilities` and `mod2linalg`
- `prime_sieve.py` — segmented sieve yielding primes ≤ n in O(n) time / O(√n) space
- `lucy.py`, `sievecntsum.py` — Lucy Hedgehog algorithm for counting/summing primes up to n without enumerating them

**Multiplicative functions**:
- `multiplicative.py` — Euler's φ, Möbius μ, and summatory versions

**Algebraic number theory** (each implements its own element type with full arithmetic):
- `gaussian_integer.py` — Z[i]; `GaussianInteger` class
- `quadratic_extensions.py` — Q(√d) discriminants, Legendre symbol, Tonelli-Shanks
- `real_quadratic_fields.py`, `imaginary_complex_fields.py` — real/imaginary quadratic fields
- `cyclotomic_fields.py` — `CyclotomicInteger` class; uses sympy for polynomial arithmetic
- `finite_field.py` — `FiniteFieldElement`; irreducible polynomial machinery via numpy

**Analytic / modular**:
- `modular_forms.py` — Eisenstein series, modular form computations; uses matplotlib
- `modular_groups.py` — SL(2,Z) and related groups
- `gauss_sums.py`, `jacobi.py`, `jacobi_sums_over_arbitrary_ff.py` — character sums
- `fourier.py` — discrete Fourier analysis

**Other number theory**:
- `continued_fractions.py`, `stern_brocot.py` — continued fractions and the Stern-Brocot tree
- `partitions.py` — integer partitions
- `diophantine.py` — Diophantine equation solvers
- `faulhaber.py` — Faulhaber's formula (power sums)
- `fibonacci.py` — Fibonacci with optional modulus
- `almost_prime.py`, `biquadratic.py`, `desboves.py`, `gauss_last_entry.py` — specialised results from the source texts

**Combinatorics and data structures**:
- `combinatorics.py` — combinatorial helpers
- `polynomial.py` — `intpoly1d` (univariate integer polynomial)
- `graph_algorithms.py` — Dijkstra and graph utilities; optional pynauty for graph isomorphism
- `algorithmx.py` — Knuth's Dancing Links (Algorithm X) for exact cover
- `mod2linalg.py` — linear algebra over GF(2)
- `slv.py` — linear equation solver
- `cholesky.py` — Cholesky decomposition
- `heaps.py`, `fenwick.py`, `maxsum_segtree.py`, `numbertrees.py` — supporting data structures

**Testing helpers**:
- `testing.py` — legacy `assert_equal`/`assert_exception` helpers (predates pytest); some older code uses these directly
- `test/` — pytest test suite; one file per module

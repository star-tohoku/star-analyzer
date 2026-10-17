# candidate_expression.py -- the expression grammar of StCandidateFemtoMaker (CandidateExpression in
# StMaker/StCandidateFemtoMaker/CandidateBdtScorer.h), vectorised with numpy.
#
# A model's selections and input expressions are evaluated twice: here (training, apply, check) and in the Maker.
# Both sides must give the same numbers, so this file is a line-by-line port of the C++ parser:
#   + - * / ; < <= > >= == != ; && || ! ; log10 log exp sqrt abs max min ; numbers ; variable names.
#   Edge values as ROOT 5.34 TTreeFormula: log/log10 of x <= 0 = 0, sqrt(x) = sqrt(|x|), exp(x > 700) = exp(700),
#   exp(x < -700) = 0, x / 0 = 0. Comparisons and logic give 1.0 / 0.0.
#   && and || may not be mixed without parentheses; comparisons do not chain.
# Values are float64 (tree branches are float32 or int32, converted exactly, as in FillCandidateModelValues).
import re
import numpy as np

_NUMBER = re.compile(r'(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?')
_UNARY = {'log10': 1, 'log': 2, 'exp': 3, 'sqrt': 4, 'abs': 5}


class ExpressionError(ValueError):
    pass


class CandidateExpression:
    """Compile once (names = the variables the expression may use), evaluate on a dict of numpy arrays."""

    def __init__(self, text, names):
        self.text, self.names, self.code, self.pos = text, list(names), [], 0
        self._or()
        self._skip()
        if self.pos != len(self.text):
            raise ExpressionError(f"'{text}': unexpected '{self.text[self.pos:]}'")
        self.variables = sorted({v for k, v in self.code if k == 'var'})

    # --- parser (same structure and checks as the C++ class) ---
    def _skip(self):
        while self.pos < len(self.text) and self.text[self.pos].isspace():
            self.pos += 1

    def _accept(self, c):
        self._skip()
        if self.text.startswith(c, self.pos):
            self.pos += len(c)
            return True
        return False

    def _peek(self, c):
        self._skip()
        return self.text.startswith(c, self.pos)

    def _fail(self, msg):
        raise ExpressionError(f"'{self.text}': {msg}")

    def _or(self):
        n_and, n_or = self._and(), 0
        while self._accept('||'):
            n_and += self._and()
            self.code.append(('bin', '|'))
            n_or += 1
        if n_or > 0 and n_and > 0:
            self._fail('mixed && and || without parentheses')

    def _and(self):
        self._compare()
        n = 0
        while self._accept('&&'):
            self._compare()
            self.code.append(('bin', '&'))
            n += 1
        return n

    def _compare(self):
        self._sum()
        op = None
        for tok, o in (('<=', 'l'), ('>=', 'g'), ('==', '='), ('!=', 'n'), ('<', '<'), ('>', '>')):
            if self._accept(tok):
                op = o
                break
        if op is None:
            return
        self._sum()
        self.code.append(('bin', op))
        if self._peek('<') or self._peek('>') or self._peek('=') or self._peek('!='):
            self._fail('chained comparison (use &&)')

    def _sum(self):
        self._product()
        while True:
            if self._accept('+'):
                self._product(); self.code.append(('bin', '+'))
            elif self._peek('-'):
                self.pos += 1; self._product(); self.code.append(('bin', '-'))
            else:
                break

    def _product(self):
        self._unary()
        while True:
            if self._accept('*'):
                self._unary(); self.code.append(('bin', '*'))
            elif self._accept('/'):
                self._unary(); self.code.append(('bin', '/'))
            else:
                break

    def _unary(self):
        if self._accept('-'):
            self._unary(); self.code.append(('un', 0)); return
        if self._accept('+'):
            self._unary(); return
        self._skip()
        if self.text.startswith('!', self.pos) and not self.text.startswith('!=', self.pos):
            self.pos += 1
            self._unary(); self.code.append(('un', 6)); return
        self._primary()

    def _primary(self):
        self._skip()
        if self.pos >= len(self.text):
            self._fail('unexpected end of expression')
        if self._accept('('):
            self._or()
            if not self._accept(')'):
                self._fail("missing ')'")
            return
        c = self.text[self.pos]
        if c.isdigit() or c == '.':
            m = _NUMBER.match(self.text, self.pos)
            if not m:
                self._fail(f"bad number at '{self.text[self.pos:]}'")
            self.pos = m.end()
            self.code.append(('num', float(m.group(0))))
            return
        if not (c.isalpha() or c == '_'):
            self._fail(f"unexpected character '{c}'")
        m = re.compile(r'[A-Za-z0-9_]+').match(self.text, self.pos)
        ident = m.group(0)
        self.pos = m.end()
        if self._accept('('):
            if ident in _UNARY:
                self._or()
                if not self._accept(')'):
                    self._fail(f"missing ')' after {ident}")
                self.code.append(('un', _UNARY[ident]))
                return
            if ident in ('max', 'min'):
                self._or()
                if not self._accept(','):
                    self._fail(f'{ident} needs two arguments')
                self._or()
                if not self._accept(')'):
                    self._fail(f"missing ')' after {ident}")
                self.code.append(('bin', 'M' if ident == 'max' else 'm'))
                return
            self._fail(f"unknown function '{ident}'")
        if ident not in self.names:
            self._fail(f"unknown variable '{ident}'")
        self.code.append(('var', ident))

    # --- evaluation (same edge values as the C++ class) ---
    def evaluate(self, values):
        """values: dict name -> numpy array (any numeric dtype). Returns float64 array."""
        n = len(next(iter(values.values())))
        st = []
        with np.errstate(all='ignore'):
            for kind, arg in self.code:
                if kind == 'num':
                    st.append(np.full(n, arg, dtype=np.float64))
                elif kind == 'var':
                    st.append(np.asarray(values[arg], dtype=np.float64))
                elif kind == 'un':
                    a = st.pop()
                    if arg == 0: a = -a
                    elif arg == 1: a = np.where(a > 0., np.log10(np.where(a > 0., a, 1.)), 0.)
                    elif arg == 2: a = np.where(a > 0., np.log(np.where(a > 0., a, 1.)), 0.)
                    elif arg == 3: a = np.where(a > 700., np.exp(700.), np.where(a < -700., 0., np.exp(np.clip(a, -700., 700.))))
                    elif arg == 4: a = np.sqrt(np.abs(a))
                    elif arg == 5: a = np.abs(a)
                    else: a = np.where(a == 0., 1., 0.)
                    st.append(a)
                else:
                    b = st.pop(); a = st.pop()
                    if arg == '+': r = a + b
                    elif arg == '-': r = a - b
                    elif arg == '*': r = a * b
                    elif arg == '/': r = np.where(b == 0., 0., a / np.where(b == 0., 1., b))
                    elif arg == 'M': r = np.where(a > b, a, b)
                    elif arg == 'm': r = np.where(a < b, a, b)
                    elif arg == '<': r = (a < b).astype(np.float64)
                    elif arg == 'l': r = (a <= b).astype(np.float64)
                    elif arg == '>': r = (a > b).astype(np.float64)
                    elif arg == 'g': r = (a >= b).astype(np.float64)
                    elif arg == '=': r = (a == b).astype(np.float64)
                    elif arg == 'n': r = (a != b).astype(np.float64)
                    elif arg == '&': r = ((a != 0.) & (b != 0.)).astype(np.float64)
                    else: r = ((a != 0.) | (b != 0.)).astype(np.float64)
                    st.append(r)
        return st[-1] if st else np.zeros(n)

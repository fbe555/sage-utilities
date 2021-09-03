from IPython.display import Latex
from sage.misc.latex import LatexExpr
import operator
import traceback
from keyword import iskeyword
import numpy as np

latex.extra_preamble('')
latex.add_macro(r'\def\bs#1{{\boldsymbol #1}}')
latex.matrix_delimiters('[', ']')


def LMatrix(mat, name=None, hidden=False, wrapper='', **kwargs):
    if 'skip_zeroes' not in kwargs: kwargs['skip_zeroes'] = True
    if 'truncate' not in kwargs: kwargs['truncate'] = True
    
    s = wrapper + ' '
    if name is not None and not hidden:
        s += '{} = '.format(name)
    s += r'\begin{bmatrix}'
    for r in mat:
        for e in r:
            if parent(e) is RR:
                e_str = e.str(**kwargs)
            else:
                e_str = str(e)
            s += r' {} &'.format(e_str)
        s = s[:-2] + r' \\'
    s = s[:-2] + r'\\end{bmatrix} ' + wrapper
    if name is not None:
        globals()[name] = Matrix(mat)
    return s


def LVector(vec, name=None, *args, **kwargs):
    s = LMatrix([(e,) for e in vec], name, *args, **kwargs) 
    if name is not None:
        globals()[name] = vector(vec)
    return s


def LList(l, name=None, hidden=False, wrapper='$', **kwargs):
    if 'skip_zeroes' not in kwargs: kwargs['skip_zeroes'] = True
    if 'truncate' not in kwargs: kwargs['truncate'] = True
       
    s = wrapper + ' '
    if name is not None and not hidden:
        s += f'{name} = '
    l_str = str(l)
    s += f'\\left{l_str[0]}'
    for e in l:
        if parent(e) is RR:
            e_str = e.str(**kwargs)
        else:
            e_str = str(e)
        s += e_str + ', '
    s = s[:-2] + f'\right{l_str[-1]} '
    s += wrapper
    if name is not None:
        globals()[name] = l
    return s

    
def LScalarMul(mat, scalar, name=None):
    product = scalar*Matrix(mat)
    s = ''
    if name is not None:
        s += '{} = '.format(name)
    s += r'{} \cdot '.format(scalar)
    s += LMatrix(mat).data
    s += r' = \begin{bmatrix}'
    for r in mat:
        for e in r:
            s += r' {} \cdot {} &'.format(scalar, e)
        s = s[:-2] + r' \\'
    s = s[:-2] + r'\end{bmatrix} = ' + LMatrix(product).data 
    if name is not None:
        globals()[name] = product
    
    return Latex(s)


def LMatMul(mat1, mat2):
    mat1 = Matrix(mat1); mat2 = Matrix(mat2)
    s = LMatrix(mat1).data + r'\cdot' + LMatrix(mat2).data
    s += '=' + LMatrix([[''.join([r'{} \cdot {} + '.format(mat1[i,k], mat2[k,j]) 
                                  for k in range(mat1.ncols())])[:-3]
                        for j in range(mat2.ncols())]
                       for i in range(mat1.nrows())]).data
    s += '=' + LMatrix(mat1*mat2).data
    return Latex(s)


_operator_symbols = {'+': operator.add,
                     '-': operator.sub,
                     '/': operator.truediv,
                     '*': operator.mul}

def LElemWise(mat1, mat2, operator):
    s = LMatrix(mat1).data + str(operator) + LMatrix(mat2).data
    s += '= ' + LMatrix([['{} {} {}'.format(e1, operator, e2) 
                          for (e1, e2) in zip(r1, r2)] 
                         for (r1, r2) in zip(mat1, mat2)]).data
    s += '= ' + LMatrix(_operator_symbols[operator](Matrix(mat1), Matrix(mat2))).data
    return Latex(s)


def LCofactorDeterminant(mat, row=None, column=None, rhs_only=False):
    assert mat.is_square(), 'Matrix must be a square matrix for the determinant to be defined'
    
    if not rhs_only:
        s = r' \\left| {} \right|  = '.format(LMatrix(mat))
    else:
        s = ''
        
    if mat.nrows() == 1:
        return str(mat[0,0])
    elif mat.nrows() == 2:
        return r'{} \cdot {} - {} \cdot {}'.format(mat[0,0], mat[1,1], mat[1, 0], mat[0, 1])
    
    assert row is not None or column is not None, 'If matrix is larger than 2 x 2, either row or column must be specified'
    
    rows = range(mat.nrows()) if row is None else [row]*mat.nrows()
    columns = range(mat.ncols()) if column is None else [column]*mat.ncols()
    
    s += ' '.join([r'{} \cdot \\left| {} \right| + '
                   .format(mat[i, j], LMatrix(submatrix(mat, i, j)))
                   for i, j in zip(rows, columns) if mat[i, j] != 0])
    if mat.nrows() == 3:
        s = s[:-2] + '= '
        s += ' '.join([r'{} \cdot \\left( {} \right) + '
                       .format(mat[i, j], LCofactorDeterminant(submatrix(mat, i, j), 
                                                               rhs_only=True))
                       for i, j in zip(rows, columns) if mat[i,j] != 0])
        s = s[:-2] + '= ' + str(mat.determinant()) + '   '
    return s[:-3]


def submatrix(mat, row, col):
    return mat.delete_rows([row]).delete_columns([col])


def _is_fraction(value):
    print(value.str())
    print(RR(n(value, digits=3)).str(truncate=True, skip_zeroes=True))
    return value.str() != RR(n(value, digits=3)).str(truncate=True, skip_zeroes=True)

def _can_convert_to_ZZ(value):
    try:
        value.change_ring(ZZ)
        return True
    except AttributeError:
        try:
            ZZ(value)
            return True
        except TypeError:
            return False
    except TypeError:
        return False

def show_var(*names, approx=True, debug=False):
    if len(names) == 0:
        # First extract the line of source code where this was called
        source = traceback.extract_stack(limit=2)[0].line
        # Check if call was from a line starting with an assignment:
        if source.count('=') < 1:
            raise SyntaxError('The show_var() function can only be used on lines where a variable is assigned')
        # Then extract the variable name assigned to
        names = [name.strip() for name in source.split('=')[0].split(',')]
    for name in names:
        value = globals()[name]
        expr = LatexExpr(f'\text{{{name}}} = ') + latex(value)
        if approx:
            try:
                approx = n(value, digits=3)
                if value != approx:
                    expr += LatexExpr(f'\\sim') + latex(approx)
                elif parent(value) is QQ and not _can_convert_to_ZZ(value):
                    expr += LatexExpr(f'= {RR(approx).str(truncate=True, skip_zeroes=True)}')
                elif value.base_ring() is QQ and not _can_convert_to_ZZ(value):
                    expr += LatexExpr('=') + latex(approx)
            except Exception as e:
                if debug:
                    traceback.print_exc()
        show(expr)

        
        
def is_valid_var_name(name):
    return name.isidentifier() and not iskeyword(name)


class Table(object):
    def __init__(self, data,
                 row_label=None, row_categories=None,
                 column_label=None, column_categories=None):
        self.row_label = row_label
        self.row_categories = row_categories
        self.column_label = column_label
        self.column_categories = column_categories
        self._data = np.asarray(data, dtype=object)
    
    def __getitem__(self, key):
        return self._data[key]
    
    def __setitem__(self, key, item):
        self._data[key] = item
    
    def __call__(self, row_key, column_key):
        return self._data[self.row_categories.index(row_key)][self.column_categories.index(column_key)]
        
        
def draw_table(data,
               row_label=None, row_categories=None,
               column_label=None, column_categories=None,
               name=None, frame=False, **kwargs):
    if isinstance(data, Table):
        row_label = data.row_label
        row_categories = data.row_categories
        column_label = data.column_label
        column_categories = data.column_categories
        data = data[:]
    
    assert len(data) >= 1 and len(data[0]) >= 1
    assert name is None or is_valid_var_name(name)
    
    if 'skip_zeroes' not in kwargs: kwargs['skip_zeroes'] = True
    if 'truncate' not in kwargs: kwargs['truncate'] = True
        
    if name is not None:
        globals()[name] = Table(data,
                                row_label=row_label, 
                                row_categories=row_categories,
                                column_label=column_label, 
                                column_categories=column_categories)
    
    cat_or_label = any([row_label, column_label, 
                        row_categories, column_categories])
    
    s = r'\begin{array}{'
    if frame:
        s += '|'
    if cat_or_label:
        s += 'r|'
    s += 'c{}'.format('|' if frame else ' ') * len(data[0]) + '} '
    if frame:
        s += r'\hline '
    if row_categories:
        s += r'{\bf ' + row_label + r' \backslash ' + column_label + '} & '
        for l in column_categories:
            s += r'{\bf ' + str(l) + '}' + ' & '
        s = s[:-2] + r' \ \hline '
    for i, row in enumerate(data):
        if row_categories:
            s += r'{\bf' + str(row_categories[i]) + '} & '
        for e in row:
            if parent(e) is RR:
                s += e.str(**kwargs)
            else:
                s += str(e)
            s += ' & '
        s = s[:-2] + r'\ '
        if frame:
            s += r'\hline '
    s += r'\end{array}'
    return s
    
    
def multiply_piecewise(input_pw, input_regular):
    return piecewise([(rang, func*input_regular) for rang, func in input_pw.items()])

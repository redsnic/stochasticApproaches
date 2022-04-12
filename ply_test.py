from ply import yacc
from ply import lex

# LEXER

tokens = ('PLUS', 'MINUS', 'TIMES', 'DIV', 'NUMBER', 'OPEN', 'CLOSE')

precedence = (
    ('left', 'PLUS', 'MINUS'),
    ('left', 'TIMES', 'DIV'),
    ('right', 'UMINUS')
)

t_ignore = ' \t'

# basic tokens
t_PLUS = r'\+'
t_MINUS = r'-'
t_TIMES = r'\*'
t_DIV= r'/'
t_OPEN = r'\('
t_CLOSE = r'\)'

# tokens with action (in this case cast to int)
def t_NUMBER(t):
    r'[0-9]+'
    t.value = int(t.value)
    return t

# ignore newline and track line number
def t_ignore_newline(t):
    r'\n+'
    t.lexer.lineno += t.value.count('\n')

def t_error(t):
    print(f'Syntax Error {t.value[0]!r}')
    t.lexer.skip(1)

# PARSER

def p_expression(p):
    '''
    expression : expression PLUS expression
               | expression MINUS expression
               | expression DIV expression
               | expression TIMES expression
    '''
    p[0] = bin_op(p[2], p[1], p[3])

def p_expression_unop(p):
    '''
    expression : MINUS expression %prec UMINUS
    '''
    p[0] = -p[2]

def p_expression_brackets(p):
    '''
    expression : OPEN expression CLOSE
    '''
    p[0] = p[2]

def p_expression_value(p):
    '''
    expression : NUMBER
    '''
    p[0] = p[1]

def p_error(p):
    print(f'Syntax error at {p.value!r}')

# AUX

def bin_op(what, x1, x2):
    if what == '+':
        return x1+x2
    elif what == '-':
        return x1-x2
    elif what == '/':
        if x2 == 0:
            raise Exception(f'Division by zero detected!')
        return x1/x2
    elif what == '*':
        return x1*x2
    else:
        raise Exception(f'invalid operation {what}')


lexer = lex.lex()
parser = yacc.yacc()

while True:
    try:
        s = input('calc > ')
    except EOFError:
        break
    if not s: continue
    result = parser.parse(s)
    print(result)

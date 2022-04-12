from StochasticSimulations.ply import yacc
from StochasticSimulations.ply import lex
from StochasticSimulations.ReactionNetwork import ReactionNetwork

# LEXER

# symbol tables
param_lookup = {}
species_lookup = {}
# reaction network object
reaction_network = ReactionNetwork()

# keyword (as per documentation)
reserved = {
    'species' : 'SPECIESKW',
    'param' : 'PARAMKW',
    'emptyset' : 'EMPTYSET'
}

# tokens
tokens = [
    'PLUS', 'MINUS', 'TIMES', 'DIV', 'NUMBER', 'OPEN', 'CLOSE', 'EXP', 'FLOAT', # math expressions
    'PARAM', 'LEFTARROW', 'ARROWTAIL', 'RIGHTARROW', 'SPECIES', # reaction and param expressions
    'EQUAL', 'SEMICOLON', 'COLON', 'COMMA' # symbols 
] + list(reserved.values())

# mathematical precedence rules
precedence = (
    ('left', 'PLUS', 'MINUS'),
    ('left', 'TIMES', 'DIV'),
    ('left', 'EXP'),
    ('right', 'UMINUS')
)

# ignore whitespaces
t_ignore = ' \t'

# ignore comments (C-like)
def t_COMMENT(p):
    r'//.*'
    pass

def t_COMMENTBLOCK(p):
    r'/\*.*\*/'
    pass

# math tokens
t_PLUS = r'\+'
t_MINUS = r'-'
t_TIMES = r'\*'
t_DIV= r'/'
t_EXP= r'\^'
t_OPEN = r'\('
t_CLOSE = r'\)'
# single characters
t_EQUAL = r'='
t_SEMICOLON = r';'
t_COLON = r':'
t_COMMA = r','
# reaction arrows
t_LEFTARROW = r'<-'
t_ARROWTAIL = r'--'
t_RIGHTARROW = r'->'

# numbers

def t_FLOAT(t):
    r'[0-9]+\.[0-9]*'
    t.value = float(t.value)
    return t

def t_NUMBER(t):
    r'[0-9]+'
    t.value = int(t.value)
    return t

# species/param names

t_SPECIES = r'[A-Z][A-Za-z0-9_]*'

# manage both params and keywords (as they are lowercase)
def t_PARAM(t):
    r'[a-z][A-Za-z0-9_]*'
    t.type = reserved.get(t.value,'PARAM')
    return t

# ignore newline and track line number
def t_ignore_newline(t):
    r'\n+'
    t.lexer.lineno += t.value.count('\n')

def t_error(t):
    raise Exception(f'Syntax Error: {t.value[0]!r} at line {t.lexer.lineno}')

# PARSER

# --- core

# expression manages the initializzation of the reaction network
def p_expression(p):
    """
    expression : innerexpr
    """
    reaction_network.set_initial_state(unpack(species_lookup))
    reaction_network.initialize()
    p[0] = reaction_network

# innerexpr manages the core of the program
# there are three types of lines
# parameter definitions
# initialization definitions
# reaction definitions
def p_innerexpr(p):
    """
    innerexpr :  
            | paramdefinition SEMICOLON innerexpr
            | initialization SEMICOLON innerexpr
            | reaction SEMICOLON innerexpr
    """
    p[0] = None

# define a parameter and its value
# note that any othe parameter in the mathexpression must be defined before use
def p_paramdefinition(p):
    """
    paramdefinition : PARAMKW PARAM EQUAL mathexpression
    """
    try:
        param_lookup[p[2]]
        raise Exception(f'Error: multiple definitions of parameter {p[2]}')
    except:
        param_lookup[p[2]] = p[4]
    p[0] = None

# initialize a species to a certain value (default 0)
def p_initialization(p):
    """
    initialization : SPECIESKW SPECIES
                   | SPECIESKW SPECIES EQUAL intmathexpression
    """
    try:
        species_lookup[p[2]]
        raise Exception(f'Error: multiple inizializations of species {p[2]}')
    except:
        if len(p) == 3:
            param_lookup[p[2]] = 0
        else:
            species_lookup[p[2]] = p[4]
    p[0] = None

# define a reaction ( reactants -- (rate) -- reactants )
def p_reaction(p):
    """
    reaction : reactants ARROWTAIL mathexpression RIGHTARROW reactants
             | reactants LEFTARROW mathexpression ARROWTAIL reactants
             | reactants LEFTARROW mathexpression RIGHTARROW reactants
    """
    if p[2] == "--":
        reaction_network.add_reaction(unpack(p[1]), unpack(p[5]), p[3])
    elif p[4] == '--':
        reaction_network.add_reaction(unpack(p[5]), unpack(p[1]), p[3])
    else:
        reaction_network.add_reaction(unpack(p[1]), unpack(p[5]), p[3])
        reaction_network.add_reaction(unpack(p[5]), unpack(p[1]), p[3])
    p[0] = None

# reactants are comma separated lists of pairs SPECIES_NAME:SPECIES_QUANTITY
# if no reactants are involved, use the keyword emptyset
def p_reactants(p):
    """
    reactants : reactants COMMA SPECIES COLON intmathexpression
              | SPECIES COLON intmathexpression
              | EMPTYSET
    """
    current_reactant = {} if len(p) == 4 or len(p) == 2 else p[1]
    if len(p) != 2:
        try:
            current_reactant[p[1] if len(p) == 4 else p[3]] += (p[3] if len(p) == 4 else p[5])
        except:
            current_reactant[p[1] if len(p) == 4 else p[3]] = (p[3] if len(p) == 4 else p[5])
    p[0] = current_reactant


# --- math ---

# this section of the grammar implements basic mathematical operation

# this rule manages operations that must have integer results (type checking!)
def p_intmathexpression(p):
    '''
    intmathexpression : mathexpression
    '''
    if int(p[1]) - p[1] != 0:
        raise Exception('Error: non integer value detected')
    elif p[1] < 0:
        raise Exception('Error: negative integer detected')
    p[0] = int(p[1])

# binary operation mathematical expression
def p_mathexpression_binop(p):
    '''
    mathexpression : mathexpression PLUS mathexpression
               | mathexpression MINUS mathexpression
               | mathexpression DIV mathexpression
               | mathexpression TIMES mathexpression
               | mathexpression EXP mathexpression
    '''
    p[0] = bin_op(p[2], p[1], p[3])

# unary operation mathematical expression
def p_mathexpression_unop(p):
    '''
    mathexpression : MINUS mathexpression %prec UMINUS
    '''
    p[0] = -p[2]

# brackets
def p_mathexpression_brackets(p):
    '''
    mathexpression : OPEN mathexpression CLOSE
    '''
    p[0] = p[2]

# constant values
def p_mathexpression_value(p):
    '''
    mathexpression : numeric
    '''
    p[0] = p[1]

# parametric variables
def p_mathexpression_param(p):
    '''
    mathexpression : PARAM
    '''
    try:
        p[0] = param_lookup[p[1]]
    except KeyError:
        raise Exception(f'Error: variable {p[1]} is undefined')

# a numeric value can be both a float or an integer (NUMBER)
def p_numeric(p):
    '''
    numeric : NUMBER
            | FLOAT
    '''
    p[0] = p[1]

def p_error(p):
    raise Exception(f'Syntax error at {p.value!r}')


# --- Auxiliary functions ---

# compute a binary operation
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
    elif what == '^':
        return x1**x2
    else:
        raise Exception(f'invalid operation {what}')

# unroll a dictionary to a key value list
def unpack(dic):
    return [(k,v) for k,v in dic.items()]

### --- MAIN PROGRAM ---

def make_parser():
    lexer = lex.lex()
    return yacc.yacc()

if __name__ == '__main__':

    # create parser an lexer

    lexer = lex.lex()
    make_parser()

    # read input

    # f = open('./StochasticSimulations/example.rn')
    # text = ''.join(f.readlines())

    # uncomment to print the extracted tokens
    #lexer.input(text)
    #while True:
    #    tok = lexer.token()
    #    if not tok: 
    #        break      # No more input
    #    print(tok)

    # rn = parser.parse(text)
    # print(rn)




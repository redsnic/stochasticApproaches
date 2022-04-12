from StochasticSimulations.ReactionNetworkLanguage import make_parser
import traceback

def read_rn(parser, text):
    try:
        return parser.parse(text)
    except Exception:
        traceback.print_exc()
        raise IOError




import argparse
import sys
import logging
from typing import List, Sequence, Any, Dict, Optional

#Init logger
FMT = '%(asctime)s | [%(levelname)s] %(message)s'
DATEFMT = '%Y-%m-%d %H:%M:%S'
formatter = logging.Formatter(fmt=FMT, datefmt=DATEFMT)

h1 = logging.StreamHandler(sys.stdout)
h1.setLevel(logging.INFO)
h1.addFilter(lambda log: log.levelno == logging.INFO)
h1.setFormatter(formatter)

h2 = logging.StreamHandler(sys.stderr)
h2.setLevel(logging.WARNING)
h2.setFormatter(formatter)

logger = logging.getLogger(__name__)
logger.propagate = False
logger.setLevel(logging.DEBUG)
logger.handlers = [h1, h2]

class ParseKW(argparse.Action):
    '''
    A class to parse the dictionary-like input args.
    From https://sumit-ghosh.com/posts/parsing-dictionary-key-value-pairs-kwargs-argparse-python/

    Parameters:
    -----------
    parser : argparse.ArgumentParser
        The parser object. 
    namespace : argparse.Namespace
        An object, which holds attributes and returns it.
    values : str | Sequence[Any] | None
        A list of values which follow the argument.
    option_string : str | None
        The option string that is used to invoke the action.     
    '''
    def __call__(self,
                 parser: argparse.ArgumentParser,
                 namespace: argparse.Namespace,
                 values: str | Sequence[Any] | None,
                 option_string: str | None = None
                 ):

        setattr(namespace, self.dest, {})
        if values:
            for value in values:
                key, value = value.split('=')
                getattr(namespace, self.dest)[key] = value

def merge_args(d_args: Dict[str, Optional[str | int | bool | List[str] | Dict[str, int]]],
               config: Dict[str, Optional[str | int | bool | List[str] | Dict[str, int]]]) \
                       -> Dict[str, Optional[str | int | bool | List[str] | Dict[str, int]]]:
    '''
    A function to unite the parameters entered as the arguments
    from the console and the parameters loaded from a config file.

    Parameters:
    -----------
    d_args : Dict[str, Optional[str | int | bool | List[str] | Dict[str, int]]]
        input arguments represented as a dictionary.
    config : Dict[str, Optional[str | int | bool | List[str] | Dict[str, int]]]
        parameters loaded from a config file.
    '''
    for key in config.keys():
        if (not key in d_args.keys()) or (not d_args[key]):
            d_args[key] = config[key]
    return d_args.copy()

def check_sub_args(d_args: Dict[str, Optional[str]], 
                   required_sub_args: Dict[str, str]) -> None:
    '''
    A function to check the names of input 
    argument parameters.

    Parameters:
    -----------
    d_args : Dict[str, Optional[str]]
        Input arguments represented as a dictionary.
    required_sub_args : Dict[str, str]
        A dictionary of predetermined parameter names which
        should be used with arguments.
    '''
    for key in required_sub_args.keys():
        if d_args.get(key):
            for par in d_args[key].keys():
                if not par in required_sub_args[key]:
                    raise Exception(
                        f"The parameter's name {par} isn't appropriate, "\
                        f"should be one of {pars_appr}"
                        )


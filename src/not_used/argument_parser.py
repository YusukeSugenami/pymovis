"""Argument parser configuration for movis command-line interface.

Argument definitions are centralized in settings.py (PARSER_ARGS) for easier
maintenance and configuration management.
"""

import argparse
from settings import PARSER_ARGS

##############################################
# Argument parser
##############################################

parser = argparse.ArgumentParser(
    description="MO picture generator from Gaussian fchk or cube file."
)

# Apply all argument configurations from PARSER_ARGS
for arg_config in PARSER_ARGS:
    parser.add_argument(*arg_config["args"], **arg_config["kwargs"])


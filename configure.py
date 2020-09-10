import re


import argparse

parser = argparse.ArgumentParser()

cxx_choices = ['gcc','icc']


# --cxx=[name] argument
parser.add_argument(
    '--cxx',
     default='gcc',
     choices=cxx_choices,
     help='select C++ compiler and default set of flags')

args = vars(parser.parse_args())
makefile_input   = 'Makefile.in'
makefile_output  = 'Makefile'
makefile_options = dict()

makefile_options["LDLFLAGS"] = "-lgsl -lgslcblas"
if args["cxx"] == "gcc":
    makefile_options["LDLFLAGS"]+= " -lm"

# Read templates
with open(makefile_input, 'r') as current_file:
    makefile_template = current_file.read()

makefile_options["COMPILER_COMMAND"] = args["cxx"]

# Make substitutions
for key, val in makefile_options.items():
    makefile_template = re.sub(r'@{0}@'.format(key), val, makefile_template)

with open(makefile_output, 'w') as current_file:
    current_file.write(makefile_template)

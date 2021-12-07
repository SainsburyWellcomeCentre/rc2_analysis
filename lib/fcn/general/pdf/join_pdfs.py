import os
import sys
import argparse
from PyPDF2 import PdfFileMerger

parser = argparse.ArgumentParser(description='Join list of PDFs');
parser.add_argument('fnames', metavar='ifile', type=str, nargs='+',
                         help='list of input filenames')
parser.add_argument('-o', dest='output_fname', metavar='ofile', type=str, default='joined.pdf',
                         help='output filename (default="joined.pdf")')
parser.add_argument('-f', action='store_true')

args = parser.parse_args()

input_fnames = args.fnames
output_fname = args.output_fname

if args.f is True:
    if os.path.isfile(output_fname):
        os.remove(output_fname)

if os.path.isfile(output_fname):
    sys.exit('Output file already exists')

merger = PdfFileMerger()

for f in input_fnames:
    merger.append(f)

merger.write(output_fname)
merger.close()

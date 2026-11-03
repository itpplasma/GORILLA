#!/usr/bin/env python3
"""Run the native public-fixture driver through its loss-estimator assertions."""
import argparse
import ast
import os
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--applets', type=Path, required=True)
parser.add_argument('--core', type=Path, required=True)
parser.add_argument('--binary', type=Path, required=True)
parser.add_argument('--work', type=Path, required=True)
args = parser.parse_args()
if args.work.exists():
    parser.error('--work must name a new directory')
source = args.applets.resolve()/'TESTS/test_mono_energetic_transport_slow.py'
tree = ast.parse(source.read_text(), filename=str(source))
cutoffs = [index for index, node in enumerate(tree.body)
           if isinstance(node, ast.For) and isinstance(node.target, ast.Tuple)
           and [element.id for element in node.target.elts
                if isinstance(element, ast.Name)] == ['key', 'invalid', 'message']]
if len(cutoffs) != 1:
    raise RuntimeError('native rejection-probe boundary changed; inspect the driver')
# Rejection probes intentionally truncate the successful result files.
tree.body = tree.body[:cutoffs[0]]
os.environ.update(APPLETS_ROOT=str(args.applets.resolve()),
                  GORILLA_ROOT=str(args.core.resolve()),
                  GORILLA_APPLETS_BIN=str(args.binary.resolve()),
                  WORK_DIR=str(args.work.resolve()), TRANSPORT_COMPARE_BACKEND='cuda')
exec(compile(tree, str(source), 'exec'), {'__file__': str(source), '__name__': '__main__'})

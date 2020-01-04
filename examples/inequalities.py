from problemgen.container import *
import random
import argparse
import time

parser = argparse.ArgumentParser()
parser.add_argument("wn_number", help='Number of the worksheet')
args = parser.parse_args()

wn_number = args.wn_number
todaydate = time.strftime("%m-%d-%Y")

w = Worksheet('%s-%s.tex' % (todaydate, wn_number))
w.set_title('Worksheet %s for the week of %s' % (wn_number, todaydate))
w.set_message('Before you start this worksheet, outline (step by step) how you approach solving an inequality. Then, outline (step by step) how you approach solving an equation. Then, solve each equation for x. Show and check your work. Use interval notation for the solution to the inequalities.')
for _ in range(3):
    w.add_linear(num_lhs_terms = 3, num_rhs_terms=3, order_rhs=1, middle_sign='>')
    w.add_linear(num_lhs_terms = 3, num_rhs_terms=3, order_rhs=1, middle_sign='<')
    w.add_linear(num_lhs_terms = 3, num_rhs_terms=3, order_rhs=1, middle_sign='<=')
    w.add_linear(num_lhs_terms = 3, num_rhs_terms=3, order_rhs=1, middle_sign='>=')
    w.add_linear(num_lhs_terms = 3, num_rhs_terms=3, order_rhs=1, middle_sign='=')
w.shuffle()
w.make()

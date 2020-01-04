from problemgen import *
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
w.set_message('Solve each equation for x. Show and check your work.')
for _ in range(10):
    w.add_equation(num_lhs_terms=random.randint(1,4), num_rhs_terms=random.randint(1,4),
            order_lhs=1, order_rhs=1)
w.make()

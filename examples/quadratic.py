from problemgen import *

w = Worksheet('8-4-18.tex')
w.set_title('8-4-18')
w.set_message('Solve each quadratic equation by factoring (with grouping, if appropriate). Check all of your answers by plugging the solutions back into the original equation.')
for _ in range(100):
    w.add_quadratic(max_lowest_term=4, leading_coeff=True)
w.make()
w.show()

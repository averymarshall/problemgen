import random
import os
import math as m
import pdb
import sys
import linecache
import inflect

from sympy import *
from sympy.solvers.inequalities import solve_poly_inequalities
from sympy.solvers.solveset import linsolve

# From Apogentus on stackexchange
def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN (%s, LINE %d "%s"): %s' % (filename, lineno, line.strip(), exc_obj))

#TODO: add support for different bases
class Number:
    '''
    Object designed to store a number in several different forms.

    Member variables:
    num           -   Number stored as an float, e.g. 42
    word            -   Number stored as a string, e.g. forty-two
    expanded        -   Number stored as a list of tuples of powers of self.base.
                        428 would be stored as [(4, 2), (2, 1), (8, 0)]
                        (4 * 10^2 + 2 * 10^1 + 8 * 10^0)
    scientific      -   Number stored as scientific notation. Stored as a tuple,
                        e.g. (4.2, 1) to represent 42 = 4.2 * 10^1
    base            -   base to store the number in.
    '''

    def __init__(self, num):
        '''
        Arguments:

        num: number to use to initialize the Number
        base: base the number should be represented in (default 10)
        '''
        self.num = float(num)
        self.base = 10
        self.word = self.make_word()
        self.expanded = self.make_expanded()
        self.scientific = self.make_scientific()

    def make_word(self):
        '''
        Function to generate a new word given self.num.

        Returns a string representing self.word.
        '''
        e = inflect.engine()
        return e.number_to_words(self.num)

    def make_expanded(self):
        '''
        Function to generate a new expanded form given self.num.
        428 would be stored as [(4, 2), (2, 1), (8, 0)]
        (4 * 10^2 + 2 * 10^1 + 8 * 10^0)
        '''
        left_of_point = []
        right_of_point = []
        str_num = str(self.num)
        if '.' in str_num:
            # This float has components to the right and left of the
            # decimal point that need to be dealt with separately.
            str_left_of_point = str_num[:str_num.index('.')]
            str_right_of_point = str_num[str_num.index('.')+1:]
            for i, str_n in enumerate(str_left_of_point[::-1]):
                left_of_point.append((int(str_n), i))
            # Switching list to descending order (3, 2, 1)
            left_of_point[::-1]
            for i, str_n in enumerate(str_right_of_point):
                right_of_point.append((int(str_n), -i))
            return left_of_point + right_of_point
        else:
            # no decimal point
            for i, str_n in enumerate(str_num[::-1]):
                left_of_point.append((int(str_n), i))
            # Switching list to descending order (3, 2, 1)
            left_of_point[::-1]
            return left_of_point

    def make_scientific(self):
        '''
        Returns a number in scientific form from self.num.
        Returns as (float, exponent) to represent float * base^exponent.
        '''
        exp = 0
        num = self.num
        while num >= self.base:
            num = num / self.base
            exp += 1
        while num < 1:
            num = num * self.base
            exp -= 1
        return (num, exp)

    def set_num(self, num):
        self.num = num
        self.word = self.make_word()
        self.expanded = self.make_expanded()
        self.scientific = self.make_scientific()

class Term:
    '''
    Object designed to store several aspects of an individual term used to
    create an Expression.

    Member variables:

    term :          A SymPy expression of the term itself.
    latex_term :    A form of the term formatted as latex.
    str_term :      A form of the term formatted as a Python string.
    expand :        Determines if any operations between this term and
                    another should expand fully.
    '''

    def __init__(self, sympy_term, expand=True):
        '''
        Takes one Sympy expression as input to initalize the term.
        '''
        self.sympy_term = sympy_term
        self.latex_term = latex(self.sympy_term)
        self.str_term = str(self.sympy_term)

    # Overloaded operators

    def __add__(self, other):
        if not expand:
            return Term(self.sympy_term + other.sympy_term)
        return Term(expand(self.sympy_term + other.sympy_term))

    def __sub__(self, other):
        if not expand:
            return Term(self.sympy_term - other.sympy_term)
        return Term(expand(self.sympy_term - other.sympy_term))

    def __mul__(self, other):
        if not expand:
            return Term(self.sympy_term * other.sympy_term)
        return Term(expand(self.sympy_term * other.sympy_term))

    def __truediv__(self, other):
        if not expand:
            return Term(self.sympy_term / other.sympy_term)
        return Term(expand(self.sympy_term / other.sympy_term))

    def __pow__(self, other):
        if not expand:
            return Term(self.sympy_term ** other.sympy_term)
        return Term(expand(self.sympy_term ** other.sympy_term))

    def __str__(self):
        return self.str_term

class Expression:
    '''
    Object designed to combine a list of Terms into a single expression.

    Member variables:

    unreduced_terms :    A list of terms meant to combine in a manner
                        perserving their initial values, hence
                        disallowing the usage of SymPy to automatically
                        combine them. Consider the usage of UnevaluatedExpr
                        when creating these terms.
    reduced_terms :    A list of terms meant to combine into a single term,
                        hence SymPy can be used to combine them without
                        consequence.
    operations :        A list of the operations between the terms. The first
                        operation in this list corresponds to the first and
                        second terms in either of the terms lists, and so
                        forth. These operations are restricted to +-*/.
    '''

    def __init__(self, unreduced_terms, reduced_terms, operations):
        assert len(operations) == len(unreduced_terms) + 1
        assert len(unreduced_terms) == len(reduced_terms)
        self.unreduced_terms = unreduced_terms[:]
        self.reduced_terms = reduced_terms[:]
        self.operations = operations
        self.zero_clean()

    def zero_clean(self):
        '''
        Deleting any terms that are zero in unreduced form unless it is the only term.
        this is so that the expression doesn't show stuff like 7x + 0 + 1

        THIS FUNCTION IS SCREWED IF YOU ARE MULTIPLYING BY ZERO AND NOT JUST ADDING IT
        FIX BUG
        '''
        if len(self.unreduced_terms) <= 1:
            return
        for i in reversed(range(len(self.unreduced_terms))):
            if self.unreduced_terms[i].sympy_term == 0:
                del self.unreduced_terms[i]
                del self.reduced_terms[i]
                del self.operations[i]

    def copy(self):
        '''
        Returns a copy of the Expression.
        '''
        return Expression(self.unreduced_terms, self.reduced_terms, self.operations)

    def simplify(self):
        '''
        Reduces self.unreduced_terms and self.reduced_terms to lowest terms
        and reduces the operations list to ['', '']
        '''
        unreduced_term = self.combine_terms(self.unreduced_terms,
                self.operations)
        reduced_term = self.combine_terms(self.reduced_terms,
                self.operations)
        self.unreduced_terms = [unreduced_term]
        self.reduced_terms = [reduced_term]
        self.operations = ['', '']

    def get_sympy(self):
        '''
        Returns the sympy term that represents the reduced version of the
        reduced terms.
        '''
        e = self.copy()
        e.simplify()
        return e.reduced_terms[0].sympy_term

    def combine_terms(self, terms, ops):
        '''
        Takes a list of SymPy terms and a list of the operations between the terms
        and combines them in a manner that respects the order of operations.
        Is recursive.

        Arguments:
        terms - a list of Terms.
        ops - a list of operations (*, /, +, -, ^) that apply to the terms. (The first
        op corresponds to the operation between the first and second terms, and so
        forth.) This list should be one element longer than the list of terms.
        Parentheticals are denoted as *(, )*, or )*( for any of the five operations
        listed above.

        Returns:
        A single Term
        '''
        # Copying lists to avoid strange bugs
        # hashtaghonestcomments hashtagthistookmeadaytofigureout
        terms = terms[:]
        ops = ops[:]

        assert len(ops) == len(terms) + 1

        #########################
        # Base case, a single term
        if len(terms) == 1:
            return terms[0]

        #########################
        # Enclosure marks case
        # Resolving parentheticals
        # list of all locations of '(' character
        left_most_locs = [i for i in range(len(ops)) if '(' in ops[i]]
        # list of all locations of ')' character
        right_most_locs = [i for i in range(len(ops)) if ')' in ops[i]]
        # Making sure all parentheticals are closed
        assert len(left_most_locs) == len(right_most_locs)
        if len(left_most_locs) != 0:
            assert left_most_locs[0] < right_most_locs[0]
            # Finding the parentheses closest to the last (
            closest_dist = 10000
            closest_loc = -1
            for r in right_most_locs:
                distance = r - left_most_locs[-1]
                if distance > 0 and distance < closest_dist:
                    closest_dist = distance
                    closest_loc = r
            # Now that the enclosed operation has been denoted, it can be
            # combined.
            new_ops = [''] + ops[left_most_locs[-1]+1:closest_loc] + ['']
            new_term = self.combine_terms(terms[left_most_locs[-1]:closest_loc], new_ops)
            # Removing the enclosed operation and replacing it with the resulting term
            del ops[left_most_locs[-1] + 1: closest_loc]
            del terms[left_most_locs[-1] : closest_loc]
            # Inserting the new term
            terms.insert(left_most_locs[-1], new_term)
            # Removing enclosure marks
            ops[left_most_locs[-1]] = ops[left_most_locs[-1]].replace('(', '')
            ops[left_most_locs[-1]+1] = ops[left_most_locs[-1]+1].replace(')', '')
            return self.combine_terms(terms, ops)

        #########################
        # Exponent case
        # Finding any instances of ^ symbol
        if '^' in ops:
            caret_loc = ops.index('^')
            # deleting the caret
            del ops[caret_loc]
            # terms[caret_loc - 1] is being raised to terms[caret_loc].
            terms[caret_loc] = terms[caret_loc - 1] ** terms[caret_loc]
            del terms[caret_loc - 1]
            # Recursion is now ready.
            return self.combine_terms(terms, ops)

        #########################
        # Multiplication/Division case
        # Finding all multiplication and division first:
        # Mul/Division is done left to right, so the first occurence of either must
        # be found
        mul_found = False
        div_found = False
        if '*' in ops:
            mul_loc = ops.index('*')
            mul_found = True
        if '/' in ops:
            div_loc = ops.index('/')
            div_found = True
        # If the mul/div was found, the operation is performed, and the two input
        # lists are modified accordingly.
        if (mul_found and div_found and mul_loc < div_loc) or \
                (mul_found and not div_found):
            # mul is first
            # deleting the multiplication
            del ops[mul_loc]
            # terms[mul_loc] and terms[mul_loc - 1] are the terms being mulitplied.
            # terms[mul_loc] is replaced with the product of these two terms.
            terms[mul_loc] = terms[mul_loc - 1] * terms[mul_loc]
            # terms[mul_loc-1] can now be removed.
            del terms[mul_loc-1]
            # The recursion is now ready.
            return self.combine_terms(terms, ops)
        if (mul_found and div_found and div_loc < mul_loc) or \
                (div_found and not mul_found):
            # div is first
            # deleting the division
            del ops[div_loc]
            # terms[div_loc] and terms[div_loc - 1] are the two terms being divided.
            # terms[div_loc] is replaced with the quotient of these two terms.
            terms[div_loc] = terms[div_loc-1] / terms[div_loc]
            # terms[div_loc-1] can now be removed.
            del terms[div_loc - 1]
            # The recursion is now ready.
            return self.combine_terms(terms, ops)

        #########################
        # Addition/Subtraction case
        # Now the addition and subtraction terms are ready to be found.
        # Add/Sub is done left to right, so the first occurence of either must
        # be found
        add_found = False
        sub_found = False
        if '+' in ops:
            add_loc = ops.index('+')
            add_found = True
        if '-' in ops:
            sub_loc = ops.index('-')
            sub_found = True
        # If the add/sub was found, the operation is performed, and the two input
        # lists are modified accordingly.
        if (add_found and sub_found and add_loc < sub_loc) or \
                (add_found and not sub_found):
            # add is first
            # deleting the addition
            del ops[add_loc]
            # terms[add_loc] and terms[add_loc - 1] are the two terms being summed.
            # terms[add_loc] is replaced with the product of these two terms.
            terms[add_loc] = terms[add_loc-1] + terms[add_loc]
            # terms[add_loc-1] can now be removed.
            del terms[add_loc-1]
            # The recursion is now ready.
            return self.combine_terms(terms, ops)
        if (sub_found and add_found and sub_loc < add_loc) or \
                (sub_found and not add_found):
            # sub is first
            # deleting the subtraction
            del ops[sub_loc]
            # terms[sub_loc] and terms[sub_loc - 1] are the two terms being divided.
            # terms[sub_loc + 1] is replaced with the quotient of these two terms.
            terms[sub_loc] = terms[sub_loc-1] - terms[sub_loc]
            # terms[sub_loc] can now be removed.
            del terms[sub_loc-1]
            # The recursion is now ready.
            return self.combine_terms(terms, ops)

    # Overloaded operators
    def __str__(self):
        return str(Problem(self))

    def __add__(self, other):
        unreduced_terms = self.unreduced_terms[:]
        reduced_terms = self.reduced_terms[:]
        unreduced_terms.extend(other.unreduced_terms)
        reduced_terms.extend(other.reduced_terms)
        # adding parentheticals to preserve order of operations
        ops = ['('] + self.operations[1:-1] + [')+('] + other.operations[1:-1] + \
                [')']
        return Expression(unreduced_terms, reduced_terms, ops)

    def __sub__(self, other):
        unreduced_terms = self.unreduced_terms[:]
        reduced_terms = self.reduced_terms[:]
        unreduced_terms.extend(other.unreduced_terms)
        reduced_terms.extend(other.reduced_terms)
        # adding parentheticals to preserve order of operations
        ops = ['('] + self.operations[1:-1] + [')-('] + other.operations[1:-1] + \
                [')']
        return Expression(unreduced_terms, reduced_terms, ops)

    def __mul__(self, other):
        unreduced_terms = self.unreduced_terms[:]
        reduced_terms = self.reduced_terms[:]
        unreduced_terms.extend(other.unreduced_terms)
        reduced_terms.extend(other.reduced_terms)
        # adding parentheticals to preserve order of operations
        ops = ['('] + self.operations[1:-1] + [')*('] + other.operations[1:-1] + \
                [')']
        return Expression(unreduced_terms, reduced_terms, ops)

    def __truediv__(self, other):
        unreduced_terms = self.unreduced_terms[:]
        reduced_terms = self.reduced_terms[:]
        unreduced_terms.extend(other.unreduced_terms)
        reduced_terms.extend(other.reduced_terms)
        # adding parentheticals to preserve order of operations
        ops = ['('] + self.operations[1:-1] + [')/('] + other.operations[1:-1] + \
                [')']
        return Expression(unreduced_terms, reduced_terms, ops)

    def __pow__(self, other):
        unreduced_terms = self.unreduced_terms[:]
        reduced_terms = self.reduced_terms[:]
        unreduced_terms.extend(other.unreduced_terms)
        reduced_terms.extend(other.reduced_terms)
        # adding parentheticals to preserve order of operations
        ops = ['('] + self.operations[1:-1] + [')^('] + other.operations[1:-1] + \
                [')']
        return Expression(unreduced_terms, reduced_terms, ops)

class Equation:
    '''
    Object designed to take two Expressions and store them as part of an equation.
    Member variables:

    lhs             -   The Expression stored on the left hand side of the equation.
    rhs             -   The Expression stored on the right hand side of the equation.
    middle_sign     -   The sign equating the lhs and rhs. Typically '=', but supports
                        inequalities.
    variable        -   The variable used in the equation.
    '''

    def __init__(self, lhs, rhs, variable='x', middle_sign='='):
        assert middle_sign == '=' or middle_sign == '>' or middle_sign == '<' or \
                middle_sign == '>=' or middle_sign == '<='
        assert isinstance(lhs, Expression)
        assert isinstance(rhs, Expression)
        self.lhs = lhs
        self.rhs = rhs
        self.middle_sign = middle_sign
        self.variable = Symbol(variable)

class System:
    '''
    Object designed to store multiple Equations to represent a system of equations.
    Member variables:

    equations       -   A list of Equations defining the system of equations.
    variables       -   A list of all the variables contained within the system.
    '''

    # TODO: add support for systems of inequalities
    def __init__(self, equations):
        assert isinstance(equations, list) or isinstance(equations, Equation)
        if isinstance(equations, list):
            self.equations = equations
        elif isinstance(equations, Equation):
            self.equations = [equations]
        self.variables = [e.variable for e in equations]

class Problem:
    '''
    Object designed to take an Expression, Equation, or System and store it
    and its solution in a print ready format (typically in a typesetting
    language like LaTeX or Markdown).

    Member variables:

    latex_question  :   The question formatted as a LaTeX string.
    latex_solution  :   The solution formatted as a LaTeX string.
    str_question    :   The question formatted as a Python string.
    str_solution    :   The solution formatted as a Python string.
    '''

    def __init__(self, data):
        '''
        Initializes the Problem.

        Arguments:

        data :      The Expression, Equation, System, or tuple intended to be
                    stored in the problem object.
                    If data is a tuple, it should be in the form
                    (latex_question, latex_solution, str_question, str_solution)
        '''
        self.latex_question = ''
        self.latex_solution = ''
        self.str_question = ''
        self.str_solution = ''
        if isinstance(data, Expression):
            self.str_question = self.str_question_from_expression(data)
            self.str_solution = self.str_solution_from_expression(data)
            self.latex_question = self.latex_question_from_expression(data)
            self.latex_solution = self.latex_solution_from_expression(data)
        elif isinstance(data, Equation):
            self.str_question = self.str_question_from_equation(data)
            self.str_solution = self.str_solution_from_equation(data)
            self.latex_question = self.latex_question_from_equation(data)
            self.latex_solution = self.latex_solution_from_equation(data)
        elif isinstance(data, System):
            self.str_question = self.str_question_from_system(data)
            self.str_solution = self.str_solution_from_system(data)
            self.latex_question = self.latex_question_from_system(data)
            self.latex_solution = self.latex_solution_from_system(data)
        elif isinstance(data, tuple):
            self.latex_question = data[0]
            self.latex_solution = data[1]
            self.str_question = data[2]
            self.str_solution = data[3]

    def __str__(self):
        '''
        String representation of the Problem.
        '''
        return "Question: " + self.str_question + "\n" + "Solution: " + \
                self.str_solution

    def str_question_from_system(self, s):
        str_question = ''
        for e in s.equations:
            str_question += self.str_question_from_equation(e) + '\n'
        return str_question

    # TODO: make the equals signs in the equations be aligned. Use
    # \begin{align*} ?
    def latex_question_from_system(self, s):
        latex_question = ''
        for e in s.equations:
            latex_question += self.latex_question_from_equation(e) + ' \\\\ '
        return latex_question

    def str_solution_from_system(self, s):
        str_solution = ''
        # Getting list of equations such that they are equal to 0
        eqs = [(e.lhs - e.rhs).get_sympy() for e in s.equations]
        solutions = linsolve(eqs, tuple(s.variables))
        if len(solutions) == 0:
            return 'No solution'
        # Iterating through solution
        for solution in solutions:
            for i, sol in enumerate(solution):
                str_solution += str(s.variables[i]) +  ' = ' + str(sol) + ', '
            str_solution = str_solution[:-1] # deleting last extra comma
            str_solution += '\n'
        return str_solution

    def latex_solution_from_system(self, s):
        latex_solution = ''
        # Getting list of equations such that they are equal to 0
        eqs = [(e.lhs - e.rhs).get_sympy() for e in s.equations]
        solutions = linsolve(eqs, tuple(s.variables))
        if len(solutions) == 0:
            return '\\text{No solution}'
        # Iterating through solution
        for solution in solutions:
            for i, sol in enumerate(solution):
                latex_solution += latex(s.variables[i]) +  ' = ' + latex(sol) + ','
            latex_solution = latex_solution[:-1] # deleting last extra comma
            latex_solution += '\\\\'
        return latex_solution

    def str_question_from_equation(self, e):
        return self.str_question_from_expression(e.lhs) + \
                e.middle_sign + \
                self.str_question_from_expression(e.rhs)

    def latex_question_from_equation(self, e):
        return self.latex_question_from_expression(e.lhs) + \
                self.convert_op_to_latex(e.middle_sign) + \
                self.latex_question_from_expression(e.rhs)

    def str_solution_from_equation(self, e):
        lhs_term = e.lhs.combine_terms(e.lhs.reduced_terms, e.lhs.operations).sympy_term
        rhs_term = e.rhs.combine_terms(e.rhs.reduced_terms, e.rhs.operations).sympy_term
        if e.middle_sign == '=':
            solutions = solve(lhs_term - rhs_term, e.variable)
        else:
            # This equation defines an inequality
            solutions = solve_poly_inequality(Poly(lhs_term - rhs_term,
                e.variable, domain='ZZ'), e.middle_sign)

        if len(solutions) == 0:
            # no solutions
            return 'No solution'

        str_solution = ''
        if e.middle_sign == '=':
            for s in solutions:
                str_solution += str(e.variable) + e.middle_sign + ' '  + str(s) + ', '
            # Deleting last comma and space
            str_solution = str_solution[:-2]
        else:
            # Printing the interval notation appropriately
            # u220a is the element of symbol, u222a is the union symbol
            str_solution = str(e.variable) + ' \u220a '
            for s in solutions:
               str_solution += str(s) + ' \u222a '
            str_solution = str_solution[:-3] # deleting last union symbol

        return str_solution

    def latex_solution_from_equation(self, e):
        lhs_term = e.lhs.combine_terms(e.lhs.reduced_terms, e.lhs.operations).sympy_term
        rhs_term = e.rhs.combine_terms(e.rhs.reduced_terms, e.rhs.operations).sympy_term
        solutions = solve(lhs_term - rhs_term, e.variable)
        if e.middle_sign == '=':
            solutions = solve(lhs_term - rhs_term, e.variable)
        else:
            # This equation defines an inequality
            solutions = solve_poly_inequality(Poly(lhs_term - rhs_term,
                e.variable, domain='ZZ'), e.middle_sign)

        if len(solutions) == 0:
            # no solutions
            return '\\text{No solution}'

        latex_solution = ''
        if e.middle_sign == '=':
            for s in solutions:
                latex_solution += latex(e.variable) + e.middle_sign + ' '  + latex(s) + ', '
            # Deleting last comma and space
            latex_solution = latex_solution[:-2]
        else:
            # Printing the interval notation appropriately
            # u220a is the element of symbol, u222a is the union symbol
            latex_solution = latex(e.variable) + ' \\in '
            for s in solutions:
               latex_solution += latex(s) + ' \\cup '
            latex_solution = latex_solution[:-6] # deleting last union symbol
        if latex_solution == '':
            return '\\text{ No solution }'
        return latex_solution

    def latex_question_from_expression(self, e):
        '''
        Creates a latex string representing the question in question_terms
        and returns it.
        '''
        # Adding each operation
        latex_question = ''
        for i in range(len(e.unreduced_terms)):
            latex_question += \
                    self.convert_op_to_latex(e.operations[i]) + \
                    e.unreduced_terms[i].latex_term
        # Adding final operation
        latex_question += self.convert_op_to_latex(e.operations[-1])

        # Culling repeated signs from question
        latex_question = latex_question.replace('+-', '-')
        latex_question = latex_question.replace('-+', '-')
        latex_question = latex_question.replace('--', '+')
        latex_question = latex_question.replace('++', '+')

        return latex_question

    def latex_solution_from_expression(self, e):
        '''
        Creates a latex string representing the solution in an Expression
        and returns it.
        '''
        return e.combine_terms(e.reduced_terms, e.operations).latex_term

    def str_question_from_expression(self, e):
        '''
        Creates a Python string representing the question in question_terms
        and returns it.
        '''
        # Adding each operation
        str_question = ''
        for i in range(len(e.unreduced_terms)):
            str_question += \
                    e.operations[i] + \
                    e.unreduced_terms[i].str_term
        # Adding final operation
        str_question += e.operations[-1]
        # Culling any doubled signs
        str_question = str_question.replace('+-', '-')
        str_question = str_question.replace('-+', '-')
        str_question = str_question.replace('--', '+')
        str_question = str_question.replace('++', '+')
        # Culling addition or subtraction of 0
        # add this in, use regular expressions

        return str_question

    def str_solution_from_expression(self, e):
        '''
        Creates a Python string representing the solution in an Expression
        and returns it.
        '''
        return e.combine_terms(e.reduced_terms, e.operations).str_term

    def convert_op_to_latex(self, op):
        '''
        Converts the input string (which should be one of '+-*/') into its
        latex equivalent. Returns the conversion or the input string (if it
        didn't satisfy the input parameters).
        '''
        if op == '*':
            return '\\cdot '
        elif op == '/':
            return '\\div '
        elif op == '>=':
            return '\\geq '
        elif op == '<=':
            return '\\leq '
        elif op == '!=':
            return '\\neq '
        else:
            return op


class Generator:
    '''
    This object is designed to create Expressions representing different
    kinds of mathematical problems, and compile them into a list to be
    printed (typically as a latex document).

    Member variables:

    problem_list -      a list of Problems designed to maintain
                        every problem generated through one of the member
                        functions.
    worksheet_fn -      the name of the pdf document for the generated
                        worksheet. Is initially blank.

    '''
    def __init__(self):
        self.problem_list = []
        self.worksheet_fn = ''

    def gen_factorable_expression(self, order=2, factor_order=1, max_lowest_term=10,
            symbols='x', leading_coeff=False, mixed_var=False, len_factor=2):
        '''
        Generates an expression to the given order that is the
        product of a series of factors (default binominal) terms.

        Arguments:
        order       -   The order of the resulting polynomial.
        symbols     -   a string containing the variables that
                        should be used.
        leading_coeff   -   Set to true to have a leading coeff
                            for the binomial terms that is
                            greater than 0.
        factor_order    -   order of the binominal terms
        len_factor      -   length of individual factors, default 2.

        Returns an expression where the expanded form is in the
        reduced terms, and the factored expression is stored in the
        unreduced terms.
        '''
        # Generating factors
        if leading_coeff:
            factors = [self.gen_algebraic_expression(num_terms=len_factor, order=factor_order, symbols=symbols, mixed_var=mixed_var, max_lowest_term=max_lowest_term) for i in range(order)]
        else:
            factors = [self.gen_algebraic_expression(\
                    num_terms=len_factor, order=factor_order, coeff=[1, random.randint(0, max_lowest_term)], symbols=symbols, mixed_var=mixed_var, max_lowest_term=max_lowest_term) \
                    for i in range(order)]
        # Multiplying factors
        for f in factors:
            f.simplify()
        expr = factors[0]
        for i in range(1, len(factors)):
            expr = expr * factors[i]
        expr.simplify()
        # Changing the reduced term to reflect the initial factors
        expanded_term = expr.reduced_terms[0]
        expr.reduced_terms[0] = Term(factor(expanded_term.sympy_term))

        return expr

    def gen_equation(self, num_lhs_terms=2, num_rhs_terms=1, types='i',
            symbols='x', order_lhs=1, order_rhs=0, lhs_coeff=[],
            rhs_coeff=[], variable='x',
            mixed_var=False, max_lowest_term=10, middle_sign='=',
            max_multiple=1, same_base_root=True):
        '''
        Generates an Equation involving denoted variables to the order
        specified.

        num_lhs_terms   -   Number of terms on the left hand side of the
                            equation. Default 2.
        num_rhs_terms   -   Number of terms on the right hand side of the
                            equation. Default 1.
        types           -   Types of coefficients. 'i' -> Integers,
                            'f' -> fractions, 'r' -> square roots.
        order_lhs       -   Order of the left hand side expression.
                            Default 1
        order_rhs       -   Order of the right hand side expression.
                            Default 0.
        lhs_coeff       -   A list of the coefficients for the lhs. Should
                            be of the same length as the number of terms
                            on the left hand side, and will override any
                            automatic number generation. Default [], allowing
                            automatic number generation.
        rhs_coeff       -   A list of the coefficients for the rhs. Should
                            be of the same length as the number of terms
                            on the left hand side, and will override any
                            automatic number generation. Default [], allowing
                            automatic number generation.
        mixed_var       -   boolean determining if variables should be mixed
                            or not (xy vs. x^2 + y^2). Default False.
        max_lowest_term -   the maximum coefficient obtained through automatic
                            coefficient generation. Defautl 10.
        max_multiple    -   Maximum multiple to mulitply coefficients by. This
                            will increase the threshold for automatic numbers
                            beyond max_lowest_term. Default 1.
        same_base_root  -   Boolean determining if radical expressions should have
                            the same base root so they can reduce to a single term.
                            Default True.
        symbols         -   Variables in the equation. E.g. 'xy' will include x and y
                            terms. Default 'x'.
        variable        -   Variable to solve for.
        '''
        lhs_expr = self.gen_algebraic_expression(num_terms=num_lhs_terms,
                types=types, symbols=symbols, order=order_lhs, coeff=lhs_coeff,
                mixed_var=mixed_var, max_lowest_term=max_lowest_term,
                max_multiple=max_multiple, same_base_root=same_base_root)
        rhs_expr = self.gen_algebraic_expression(num_terms=num_rhs_terms,
                types=types, symbols=symbols, order=order_rhs, coeff=rhs_coeff,
                mixed_var=mixed_var, max_lowest_term=max_lowest_term,
                max_multiple=max_multiple, same_base_root=same_base_root)
        equation = Equation(lhs_expr, rhs_expr, middle_sign=middle_sign,
                variable=variable)
        return equation

    def gen_algebraic_expression(self, num_terms=2, types='i',
            symbols='x', order=1, mixed_var=False, coeff=[],
            max_lowest_term=10, max_multiple=1, same_base_root=True):
        '''
        Generates an Expression involving denoted variables to the order
        specified.

        Arguments:

        num_terms       -   the number of terms in the Expression.
                            Default 2.
        types           -   The types of coefficients allowed to appear.
                            This is expresssed as a string, containing only
                            the characters 'i' (positive integers), 'r'
                            (square roots), and/or 'f' (fractions).
                            Default 'i'.
        symbols         -   A string listing the variables to be used in the
                            expression. Must be single character variables.
                            Default 'x'.
        order           -   The order of the expression to be generated.
                            Default 1.
        mixed_var       -   Denotes if the expression should mix different
                            variables together (within the bounds of the order)
                            or not.
        coeff           -   A list containing coefficients to be used in creating
                            the expression. They should be in order from highest
                            order to lowest order
        max_lowest_term -   the largest number that can appear in the reduced
                            expression. (Of course, larger numbers may appear due
                            to the operations being done.)
        max_multiple    -   the maximum multiplier used in the creation of fractions
                            and radicals.
        same_base_root  -   bool determining if all radicals in the expression
                            should reduce to the same base root.
        '''
        assert order >= 0
        assert num_terms >= 1
        assert len(coeff) == num_terms or len(coeff) == 0
        # other arguments are checked when fed to gen_numerical_expression.

        # Initialzing sympy variables
        variables = []
        for i in range(len(symbols)):
            variables.append(Term(Symbol(symbols[i])))

        # Creating algebraic terms
        algebraic_terms = []

        # Making the terms
        for o in range(order+1)[::-1]:
            # Creating the terms with the highest order first
            term = Term(Rational(1, 1)) # the multiplicative identity
            if mixed_var:
                for i in range(0, o):
                    term *= random.choice(variables)
            else:
                variable = random.choice(variables)
                for i in range(0, o):
                    term *= variable
            algebraic_terms.append(term)

        # Adjusting algebraic terms to fit the size specified by num_terms
        if len(algebraic_terms) < num_terms:
            # Adding more terms
            for i in range(num_terms - len(algebraic_terms)):
                term = Term(Rational(1, 1))
                if random.choice([True, False]) and order >= 1:
                    random_order = random.randint(1, order)
                else:
                    random_order = 0
                if mixed_var:
                    for j in range(random_order):
                        term *= random.choice(variables)
                else:
                    variable = random.choice(variables)
                    for j in range(random_order):
                        term *= variable
                algebraic_terms.append(term)

        while len(algebraic_terms) > num_terms:
            # Chopping terms off from the end.
            # Use the coefficients argument and set some coefficients to 0
            # to make some other terms disappear.
            del algebraic_terms[-1]

        assert len(algebraic_terms) == num_terms

        if len(coeff) is 0:
            # Generating numeric expression
            expression = self.gen_numerical_expression(num_terms,
                    types=types, max_lowest_term=max_lowest_term,
                    max_multiple=max_multiple, same_base_root=same_base_root)
            # Combining algebraic terms into expression
            for i in range(num_terms):
                expression.unreduced_terms[i] *= algebraic_terms[i]
                expression.reduced_terms[i] *= algebraic_terms[i]
        else:
            # Generating numeric expression where every constant is 1
            expression = self.gen_numerical_expression(num_terms, op='+-',
                    types='i', max_lowest_term=1, max_multiple=1)
            # Combining algebraic terms and coeff into expression
            for i in range(num_terms):
                expression.unreduced_terms[i] *= algebraic_terms[i] * \
                        Term(coeff[i])
                expression.reduced_terms[i] *= algebraic_terms[i] * \
                        Term(coeff[i])
            # Removing zero terms
            expression.zero_clean()

        return expression

    def gen_num_conv(self, q_type='num', s_type='word', types='i',
            lower_num_bound=1, upper_num_bound=1e9):
        '''
        Generates a number conversion problem, i.e. taking a number written
        as a decimal, in word form, in expanded form, or in scientific form
        and rewriting it in one of the aforementioned forms.

        Arguments:
        q_type      -   Type of number to be converted from. Should be 'num'
                        (for a decimal print out), 'word' (for a phrase
                        representing the number), 'expand' (for expanded form),
                        'sci' (for scientific form).
        s_type      -   Type of number to be converted to. Should be 'num'
                        (for a decimal print out), 'word' (for a phrase
                        representing the number), 'expand' (for expanded form),
                        'sci' (for scientific form).
        types       -   Determines if the number should be an integer, decimal,
                        etc. should be 'i' for integers, 'd' for decimals.
        lower_num_bound     -   Lower bound of number generated
        upper_num_bound     -   Upper bound of number generated

        Returns a Problem of the given specifications.
        '''
        assert type(q_type) == str
        assert type(s_type) == str
        assert type(types) == str
        assert type(lower_num_bound) == int or type(lower_num_bound) == float
        assert type(upper_num_bound) == int or type(upper_num_bound) == float
        if types == 'i':
            num = random.randint(int(lower_num_bound), int(upper_num_bound))
        elif types == 'd':
            num = random.uniform(lower_num_bound, upper_num_bound)
        number = Number(num)

        # Setting solution based upon argument
        if s_type == 'num':
            str_solution = str(number.num)
            latex_solution = str(number.num)
            conv_to = 'number'
        elif s_type == 'word':
            str_solution = number.word
            latex_solution = '\\text{%s}' % number.word
            conv_to = 'phrase'
        elif s_type == 'expand':
            str_solution = ''
            latex_solution = ''
            for pair in number.expanded:
                str_solution += '%d * %d**%d + ' % (pair[0], number.base, pair[1])
                latex_solution += '%d \\cdot %d^{%d} + ' % (pair[0],
                        number.base, pair[1])
            # Removing extra ' + '
            str_solution = str_solution[:-3]
            latex_solution = latex_solution[:-3]
            conv_to = 'expanded form'
        elif s_type == 'sci':
            str_solution = '%d * %d**%d' % (number.scientific[0], number.base,
                    number.scientific[1])
            latex_solution = '%d \\cdot %d^{%d}' % (number.scientific[0], number.base,
                    number.scientific[1])
            conv_to = 'scientific form'
        # Setting question based upon argument
        if q_type == 'num':
            str_question = str(number.num)
            latex_question = str(number.num)
            str_question += '->(' + conv_to + ')'
            latex_question += '\\rightarrow\\text{(' + conv_to + ')}'
        elif q_type == 'word':
            str_question = number.word
            latex_question = ('\\text{%s}' % number.word)
            str_question += '->(' + conv_to + ')'
            latex_question += '\\rightarrow\\text{(' + conv_to + ')}'
        elif q_type == 'expand':
            str_question = ''
            latex_question = ''
            for pair in number.expanded:
                str_question += '%d * %d**%d + ' % (pair[0], number.base, pair[1])
                latex_question += '%d \\cdot %d^{%d} + ' % (pair[0],
                        number.base, pair[1])
            # Removing extra ' + '
            str_question = str_question[:-3]
            latex_question = latex_question[:-3]
            str_question += '->(' + conv_to + ')'
            latex_question += '\\rightarrow\\text{(' + conv_to + ')}'
        elif q_type == 'sci':
            str_question = '%d * %d**%d' % (number.scientific[0], number.base,
                    number.scientific[1])
            latex_question = '%d \\cdot %d^{%d}' % (number.scientific[0], number.base,
                    number.scientific[1])
            str_question += '->(' + conv_to + ')'
            latex_question += '\\rightarrow\\text{(' + conv_to + ')}'
        return Problem((latex_question, latex_solution, str_question,
            str_solution))

    def gen_numerical_expression(self, num_terms=2, op='+-', types='i',
            max_lowest_term=10, max_multiple=1, same_base_root=True):
        '''
        Generates an expression involving a chosen quantity of numbers and their
        solution.

        * Note : all radicals involved in the expression

        Arguments:
        num_terms   -   the number of terms in the Expression. Default 2.
        op          -   a string expressing which operations should be used.
                        String should only contain the characters +, -, *, ^, or
                        /. For example, the string '+/' would indicate the
                        problem should include addition and/or division.
                        Default '+-'.
        types       -   The types of numbers that should be allowed to appear.
                        This is expressed as a string, containing only the
                        characters 'i' (positive integers), 'r' (square roots),
                        'f' (fractions).
                        Default 'i'.
        max_lowest_term
                    -   the largest number that can appear in the reduced
                        expression. (Of course, larger numbers may appear due
                        to the operations being done.)
        max_multiple
                    -   the maximum multiplier used in the creation of fractions
                        and radicals.
        same_base_root
                    -   bool determining if all radicals in the expression
                        should reduce to the same base root.

        Returns an Expression.
        '''
        assert num_terms >= 1
        assert '+' in op or '-' in op or '*' in op or '/' in op
        assert max_lowest_term >= 1
        assert max_multiple >= 1
        assert 'i' in types or 'r' in types or 'f' in types

        # Generating list of perfect squares for use as multipliers
        perfect_squares = [x**2 for x in range(1,
            int(m.sqrt(max_multiple)) + 1)]

        # Generating terms in the expression
        terms = []
        reduced_terms = [] # these are terms allowed to be reduced by sympy
        # the root appearing in the reduced expression
        operations = ['']
        if same_base_root:
            base_root = random.randint(1, max_lowest_term)
        for n in range(num_terms):
            # Generating some constants for later use
            # 0 and 1 are random integers, 2 is a random perfect square,
            # 3 is a random multiplier
            constants = [random.randint(1, max_lowest_term),
                    random.randint(1, max_lowest_term),
                    random.choice(perfect_squares),
                    random.randint(1, max_multiple)]
            # There's a 1 in 4 chance of generating a negative number
            if not random.randint(0, 3):
                constants[0] *= -1
            if not random.randint(0, 3):
                constants[1] *= -1
            # Making sure constants[0] is smaller than constants[1]
            if constants[0] > constants[1]:
                # Swapping the two
                temp = constants[0]
                constants[0] = constants[1]
                constants[1] = temp
            # Choosing a type of number to generate
            type_of_num = random.choice(types)

            if type_of_num is 'i':
                # Generating integer
                terms.append(Term(constants[0]))
                reduced_terms.append(Term(constants[0]))
            elif type_of_num is 'r':
                # UnevaluatedExpr has to be used here so the term isn't fully
                # reduced
                terms.append(Term(constants[0] * \
                        sqrt(UnevaluatedExpr(constants[2]*base_root))))
                reduced_terms.append(Term(constants[0] * \
                        sqrt(constants[2] * base_root)))
            elif type_of_num is 'f':
                # Mul is used in this way so the fraction doesn't reduce
                terms.append(Term(Mul(constants[0] * constants[3],
                    Rational(1, constants[1] * constants[3]), evaluate=False)))
                reduced_terms.append(Term(Rational(constants[0]*constants[3],
                    constants[1]*constants[3])))
            operations.append(random.choice(op))
        # Deleting the last operation in the operations list, since it's
        # unnecessary
        del operations[-1]
        # Replacing it with a null character
        operations.append('')

        # Returning expression
        return Expression(terms, reduced_terms, operations)

    def gen_linear(self, max_lowest_term=10, max_multiple=1, types='i',
            num_lhs_terms=2, num_rhs_terms=1, lhs_coeff=[], rhs_coeff=[],
            middle_sign='=', mixed_var=False, symbols='x',
            same_base_root=True, order_lhs=1, order_rhs=0):
        '''
        Adds a linear equation to the problem list.

        num_lhs_terms   -   Number of terms on the left hand side of the
                            equation. Default 2.
        num_rhs_terms   -   Number of terms on the right hand side of the
                            equation. Default 1.
        types           -   Types of coefficients. 'i' -> Integers,
                            'f' -> fractions, 'r' -> square roots.
        symbols         -   lists the variables that should be used as a string.
                            Default 'x'.
        lhs_coeff       -   A list of the coefficients for the lhs. Should
                            be of the same length as the number of terms
                            on the left hand side, and will override any
                            automatic number generation. Default [], allowing
                            automatic number generation.
        rhs_coeff       -   A list of the coefficients for the rhs. Should
                            be of the same length as the number of terms
                            on the left hand side, and will override any
                            automatic number generation. Default [], allowing
                            automatic number generation.
        mixed_var       -   boolean determining if variables should be mixed
                            or not (xy vs. x^2 + y^2). Default False.
        max_lowest_term -   the maximum coefficient obtained through automatic
                            coefficient generation. Defautl 10.
        max_multiple    -   Maximum multiple to mulitply coefficients by. This
                            will increase the threshold for automatic numbers
                            beyond max_lowest_term. Default 1.
        same_base_root  -   Boolean determining if radical expressions should have
                            the same base root so they can reduce to a single term.
                            Default True.

        Returns Equation.
        '''
        assert order_lhs == 0 or order_lhs == 1
        assert order_rhs == 0 or order_rhs == 1
        equation = self.gen_equation(num_lhs_terms=num_lhs_terms,
                num_rhs_terms=num_rhs_terms, types=types,
                symbols=symbols, order_lhs=order_lhs, order_rhs=order_rhs,
                lhs_coeff=lhs_coeff, rhs_coeff=rhs_coeff,
                mixed_var=mixed_var, max_lowest_term=max_lowest_term,
                middle_sign=middle_sign, max_multiple=max_multiple,
                same_base_root=same_base_root)
        return equation

    def gen_quadratic(self, max_lowest_term=10, factorable=True,
            solvable=True, leading_coeff=False, middle_sign='='):
        '''
        Generates a quadratic equation that can be factorable, unfactorable, or unsolvable.

        Arguments:

        max_lowest_term     -   The lowest coefficient that can appear in the problem.
        factorable          -   bool determining if the equation should be factorable.
        solvable            -   bool determining if the equation should be solvable.
        leading_coeff       -   bool determining if a in ax^2 + bx + c should be 1 or not.
                                False means a = 1.
        middle_sign         -   Determines if this is an equation or an inequality.
                                Appropriate arguments are '=', '>=', '<=', '<', '>', '!='.

        Returns an Equation.
        '''
        if not factorable and solvable:
            # This requires b^2 - 4ac to not be a perfect square
            # Generating list of perfect squares
            perfect_squares = [i**2 for i in range(max_lowest_term)]
            # Generating list of non perfect squares
            not_perfect = []
            for i in range(max_lowest_term**2):
                if i not in perfect_squares:
                    not_perfect.append(i)
            # this is b^2 - 4ac
            discrim = random.choice(not_perfect)
            # this is b^2
            b2 = random.choice(perfect_squares)
            # this is -4ac
            product = (discrim - b2)
            if product % 4 != 0:
                # discrim is not a perfect square, so multiplying it
                # by 4 can't make it a perfect square
                discrim *= 4
                # b2 is a perfect square, so by multiplying by 4 it will
                # stay a perfect square
                b2 *= 4
                product = (discrim - b2)
            # Choosing two random divisors of -4ac to be a and c
            ac_choices = divisors(product)
            # Eliminating all divisors above max
            for i in reversed(range(len(ac_choices))):
                if ac_choices[i] > max_lowest_term:
                    del ac_choices[i]
            a = random.choice(ac_choices)
            if len(ac_choices) > 1:
                del ac_choices[ac_choices.index(a)]
            # Rememeber product is -4ac, not just ac
            c = int(random.choice(ac_choices) / (-4))
            # divisors returns all positive numbers regardless of the
            # input. If product (-4ac) is negative, a and c must have same
            # signs.
            if product < 0 and ((a < 0 and c > 0) or (a > 0 and c < 0)):
                c *= -1
            # setting b
            b = int(m.sqrt(b2))
            # b can be positive or negative and a and c can either be
            # both positive or both negative.
            if random.randint(0, 1):
                c *= -1
                a *= -1
            if random.randint(0, 1):
                b *= -1
            equation = self.gen_equation(num_lhs_terms=3, order_lhs=2,
                    lhs_coeff=[a, b, c], rhs_coeff=[0],
                    max_lowest_term=max_lowest_term, middle_sign=middle_sign)
        elif not solvable:
            # This requires b^2 - 4ac < 0. We do a procedure similiar to
            # the one above.
            # Generating list of perfect squares
            perfect_squares = [i**2 for i in range(max_lowest_term)]
            # this is b^2
            b2 = random.choice(perfect_squares)
            # This is b^2 - 4ac
            discrim = random.randint(-max_lowest_term**2, -1)
            # this is -4ac
            product = (discrim - b2)
            # Choosing two random divisors of ac to be a and c
            ac_choices = divisors(product)
            # Eliminating all divisors above max
            for i in reversed(range(len(ac_choices))):
                if ac_choices[i] > max_lowest_term:
                    del ac_choices[i]
            a = random.choice(ac_choices)
            if len(ac_choices) > 1:
                del ac_choices[ac_choices.index(a)]
            # Rememeber product is -4ac, not just ac
            c = random.choice(ac_choices) / (-4)
            # divisors returns all positive numbers regardless of the
            # input. If product (-4ac) is negative, a and c must have same
            # signs.
            if product < 0 and ((a < 0 and c > 0) or (a > 0 and c < 0)):
                c *= -1
            # setting b
            b = m.sqrt(b2)
            # b can be positive or negative and a and c can either be
            # both positive or both negative.
            if random.randint(0, 1):
                c *= -1
                a *= -1
            if random.randint(0, 1):
                b *= -1
            equation = self.gen_equation(num_lhs_terms=3, order_lhs=2,
                    lhs_coeff=[a, b, c], rhs_coeff=[0],
                    max_lowest_term=max_lowest_term, middle_sign=middle_sign)
        else:
            # Creating factorable quadratic
            rhs = Expression([Term(0)], [Term(0)], ['',''])
            lhs = self.gen_factorable_expression(order=2,
                    max_lowest_term=max_lowest_term,
                    leading_coeff=leading_coeff)
            equation = Equation(lhs, rhs, middle_sign=middle_sign)

        return equation

    def gen_frac_to_dec(self, max_lowest_term=10, max_multiple=1):
        '''
        Generates a Problem for converting fractions to decimals.
        Arguments:
        max_lowest_term
                    -   the largest number that can appear in the reduced
                        expression. (Of course, larger numbers may appear due
                        to the operations being done.)
        max_multiple
                    -   the maximum multiplier used in the creation of fractions
                        and radicals.
        Returns Problem.
        '''
        # Generating Expression
        expression = self.gen_numerical_expression(num_terms=1, types='f',
                max_lowest_term=max_lowest_term, max_multiple=max_multiple)
        # turning expression into problem
        problem = Problem(expression)
        # Getting decimal expression of the single fraction created
        decimal = float(expression.reduced_terms[0].sympy_term)
        # Formatting into Latex
        solution = '%.5f' % decimal
        # Substituting latex expression into problem
        problem.latex_solution = solution
        return problem

    def gen_dec_to_frac(self, max_lowest_term=10, max_multiple=1):
        '''
        Adds a Problem for converting decimals to fractions.
        Arguments:
        max_lowest_term
                    -   the largest number that can appear in the reduced
                        expression. (Of course, larger numbers may appear due
                        to the operations being done.)
        max_multiple
                    -   the maximum multiplier used in the creation of fractions
                        and radicals.
        Returns Problem.
        '''
        # Generating Expression
        expression = self.gen_numerical_expression(num_terms=1, types='f',
                max_lowest_term=max_lowest_term, max_multiple=max_multiple)
        # turning expression into problem
        problem = Problem(expression)
        # Getting decimal expression of the single fraction created
        decimal = float(expression.reduced_terms[0].sympy_term)
        # Formatting into Latex
        question = '%.5f' % decimal
        # Substituting latex expression into problem
        problem.latex_question = question

        return problem

    # TODO: decrease the rates of problems generated without solutions
    def gen_system(self, num_equations=2, num_lhs_terms=2, num_rhs_terms=1,
            types='i', symbols='xy', order_lhs=1, order_rhs=0, lhs_coeff=[],
            rhs_coeff=[], mixed_var=False, max_lowest_term=10, middle_sign='=',
            max_multiple=1, same_base_root=True):
        '''
        Generates a System of Equations involving denoted variables to the
        order specified.

        Arguments:

        num_equations   -   Number of Equations in the System. Default 2.
        num_lhs_terms   -   Number of terms on the left hand side of the
                            equation. Default 2.
        num_rhs_terms   -   Number of terms on the right hand side of the
                            equation. Default 1.
        types           -   Types of coefficients. 'i' -> Integers,
                            'f' -> fractions, 'r' -> square roots.
        order_lhs       -   Order of the left hand side expression.
                            Default 1
        order_rhs       -   Order of the right hand side expression.
                            Default 0.
        lhs_coeff       -   A list of the coefficients for the lhs. Should
                            be of the same length as the number of terms
                            on the left hand side, and will override any
                            automatic number generation. Default [], allowing
                            automatic number generation.
        rhs_coeff       -   A list of the coefficients for the rhs. Should
                            be of the same length as the number of terms
                            on the left hand side, and will override any
                            automatic number generation. Default [], allowing
                            automatic number generation.
        mixed_var       -   boolean determining if variables should be mixed
                            or not (xy vs. x^2 + y^2). Default False.
        max_lowest_term -   the maximum coefficient obtained through automatic
                            coefficient generation. Defautl 10.
        max_multiple    -   Maximum multiple to mulitply coefficients by. This
                            will increase the threshold for automatic numbers
                            beyond max_lowest_term. Default 1.
        same_base_root  -   Boolean determining if radical expressions should have
                            the same base root so they can reduce to a single term.
                            Default True.
        symbols         -   Variables in the equation. E.g. 'xy' will include x and y
                            terms. Default 'x'.
        '''
        equations = []
        for i in range(num_equations):
            if i > len(symbols):
                # There's more equations than variables in this equation,
                # so the latter equations are set up with the first variable
                variable = symbols[0]
            else:
                variable = symbols[i]
            equations.append(self.gen_equation(num_lhs_terms=num_lhs_terms,
                num_rhs_terms=num_rhs_terms, types=types, symbols=symbols,
                order_lhs=order_lhs, order_rhs=order_rhs, lhs_coeff=lhs_coeff,
                rhs_coeff=rhs_coeff, mixed_var=mixed_var,
                max_lowest_term=max_lowest_term, middle_sign=middle_sign,
                max_multiple=max_multiple, same_base_root=same_base_root,
                variable=variable))
        return System(equations)

################################### Error classes
class Error(Exception):
    '''Base class for exceptions in this module.'''
    pass

class GeneratorError(Error):
    '''Exception raised for errors with generation of Problems.

    Member variables:
    type_prob        -   type of problem being generated
    message     -   explanation of error
    '''

    def __init__(self, type_prob, message):
        self.type_prob = type_prob
        self.message = message

# Add problem container with all of the add methods and a problem list
# have worksheet be a child class
class ProblemContainer:
    '''
    Class designed to add and maintain a list of problems.

    Member variables:
    gen         -   Generator used to create all of the problems.
    problems    -   List of problems contained within the container.

    Constants:
    NUM_ATTEMPTS -  number of times an attempt at generating a non-duplicate
                    problem can be made before raising an Exception
    '''

    def __init__(self):
        self.gen = Generator()
        self.problems = []
        self.NUM_ATTEMPTS = 200

    def __str__(self):
        problems_str = ''
        for p in self.problems:
            problems_str += str(self.problems.index(p)) + '. ' + str(p) + '\n'
        return problems_str

    def clear_problems(self):
        '''
        Resets all problems.
        '''
        self.problem_list = []

    def shuffle(self):
        '''
        Shuffles the order of the problems.
        '''
        random.shuffle(self.problems)

    def add_problem(self, p):
        '''
        Adds a problem to problem list, not allowing duplicates.
        Returns a boolean (True if problem was added successfully,
        false otherwise)
        '''
        # Checking for duplicate
        for problem in self.problems:
            if str(problem) == str(p):
                return False
        # Adding problem
        self.problems.append(p)
        return True

    def add_algebraic_expression(self, num_terms=2, types='i',
            symbols='x', order=1, mixed_var=False, coeff=[],
            max_lowest_term=10, max_multiple=1, same_base_root=True):
        '''
        Adds a generic algebraic expression (involving denoted variables to the order
        specified) as a problem (where the goal is to simplify).

        Arguments:

        num_terms       -   the number of terms in the Expression.
                            Default 2.
        types           -   The types of coefficients allowed to appear.
                            This is expresssed as a string, containing only
                            the characters 'i' (positive integers), 'r'
                            (square roots), and/or 'f' (fractions).
                            Default 'i'.
        symbols         -   A string listing the variables to be used in the
                            expression. Must be single character variables.
                            Default 'x'.
        order           -   The order of the expression to be generated.
                            Default 1.
        mixed_var       -   Denotes if the expression should mix different
                            variables together (within the bounds of the order)
                            or not.
        coeff           -   A list containing coefficients to be used in creating
                            the expression. They should be in order from highest
                            order to lowest order
        max_lowest_term -   the largest number that can appear in the reduced
                            expression. (Of course, larger numbers may appear due
                            to the operations being done.)
        max_multiple    -   the maximum multiplier used in the creation of fractions
                            and radicals.
        same_base_root  -   bool determining if all radicals in the expression
                            should reduce to the same base root.
        '''

        try:
            for i in range(self.NUM_ATTEMPTS):
                # Generating expression
                expr = self.gen.gen_algebraic_expression(num_terms=num_terms, types=types,
                        symbols=symbols, order=order, mixed_var=mixed_var, coeff=coeff,
                        max_lowest_term=max_lowest_term, max_multiple=max_multiple, same_base_root=same_base_root)
                # Attempting to add it
                if self.add_problem(Problem(expr)):
                    return
                # Problem was a dupe, looping back
            # All of the problems generated were dupes, there likely aren't many unique
            # problems for the parameters given
            raise GeneratorError('algebraic', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            PrintException()

    def add_num_conv(self, q_type='num', s_type='word', types='i',
            lower_num_bound=1, upper_num_bound=1e9):
        '''
        Adds a number conversion problem, i.e. taking a number written
        as a decimal, in word form, in expanded form, or in scientific form
        and rewriting it in one of the aforementioned forms.

        Arguments:
        q_type      -   Type of number to be converted from. Should be 'num'
                        (for a decimal print out), 'word' (for a phrase
                        representing the number), 'expand' (for expanded form),
                        'sci' (for scientific form).
        s_type      -   Type of number to be converted to. Should be 'num'
                        (for a decimal print out), 'word' (for a phrase
                        representing the number), 'expand' (for expanded form),
                        'sci' (for scientific form).
        types       -   Determines if the number should be an integer, decimal,
                        etc. should be 'i' for integers, 'd' for decimals.
        lower_num_bound     -   Lower bound of number generated
        upper_num_bound     -   Upper bound of number generated

        '''
        try:
            for i in range(self.NUM_ATTEMPTS):
                # Generating expression
                prob = self.gen.gen_num_conv(q_type=q_type, s_type=s_type,
                        types=types, lower_num_bound=lower_num_bound,
                        upper_num_bound=upper_num_bound)
                # Attempting to add it
                if self.add_problem(prob):
                    return
                # Problem was a dupe, looping back
            # All of the problems generated were dupes, there likely aren't many unique
            # problems for the parameters given
            raise GeneratorError('num_conv', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            PrintException()
    def add_dec_to_frac(self, max_lowest_term=10, max_multiple=1):
        '''
        Adds a Problem for converting decimals to fractions.
        Arguments:
        max_lowest_term
                    -   the largest number that can appear in the reduced
                        expression. (Of course, larger numbers may appear due
                        to the operations being done.)
        max_multiple
                    -   the maximum multiplier used in the creation of fractions
                        and radicals.
        '''

        try:
            for i in range(self.NUM_ATTEMPTS):
                # Generating expression
                prob = self.gen.gen_dec_to_frac(max_lowest_term=max_lowest_term,
                        max_multiple=max_multiple)
                # Attempting to add it
                if self.add_problem(prob):
                    return
                # Problem was a dupe, looping back
            # All of the problems generated were dupes, there likely aren't many unique
            # problems for the parameters given
            raise GeneratorError('dec_to_frac', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            PrintException()

    # Fix bug for repeating decimals being truncated
    def add_frac_to_dec(self, max_lowest_term=10, max_multiple=1):
        '''
        Adds a Problem for converting fractions to decimals.

        Arguments:
        max_lowest_term
                    -   the largest number that can appear in the reduced
                        expression. (Of course, larger numbers may appear due
                        to the operations being done.)
        max_multiple
                    -   the maximum multiplier used in the creation of fractions
                        and radicals.
        '''
        try:
            for i in range(self.NUM_ATTEMPTS):
                # Generating expression
                prob = self.gen.gen_frac_to_dec(max_lowest_term=max_lowest_term,
                        max_multiple=max_multiple)
                # Attempting to add it
                if self.add_problem(prob):
                    return
                # Problem was a dupe, looping back
            # All of the problems generated were dupes, there likely aren't many unique
            # problems for the parameters given
            raise GeneratorError('frac_to_dec', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            PrintException()

    def add_equation(self, num_lhs_terms=2, num_rhs_terms=1, types='i',
            symbols='x', order_lhs=1, order_rhs=0, lhs_coeff=[],
            rhs_coeff=[],
            mixed_var=False, max_lowest_term=10, middle_sign='=',
            max_multiple=1, same_base_root=True):
        '''
        Adds an equation involving denoted variables to the order
        specified (as a Problem).

        num_lhs_terms   -   Number of terms on the left hand side of the
                            equation. Default 2.
        num_rhs_terms   -   Number of terms on the right hand side of the
                            equation. Default 1.
        types           -   Types of coefficients. 'i' -> Integers,
                            'f' -> fractions, 'r' -> square roots.
        order_lhs       -   Order of the left hand side expression.
                            Default 1
        order_rhs       -   Order of the right hand side expression.
                            Default 0.
        lhs_coeff       -   A list of the coefficients for the lhs. Should
                            be of the same length as the number of terms
                            on the left hand side, and will override any
                            automatic number generation. Default [], allowing
                            automatic number generation.
        rhs_coeff       -   A list of the coefficients for the rhs. Should
                            be of the same length as the number of terms
                            on the left hand side, and will override any
                            automatic number generation. Default [], allowing
                            automatic number generation.
        mixed_var       -   boolean determining if variables should be mixed
                            or not (xy vs. x^2 + y^2). Default False.
        max_lowest_term -   the maximum coefficient obtained through automatic
                            coefficient generation. Defautl 10.
        max_multiple    -   Maximum multiple to mulitply coefficients by. This
                            will increase the threshold for automatic numbers
                            beyond max_lowest_term. Default 1.
        same_base_root  -   Boolean determining if radical expressions should have
                            the same base root so they can reduce to a single term.
                            Default True.
        '''

        try:
            for i in range(self.NUM_ATTEMPTS):
                # Generating expression
                eq = self.gen.gen_equation(num_lhs_terms=num_lhs_terms,
                        num_rhs_terms=num_rhs_terms, types=types,
                        symbols=symbols, order_lhs=order_lhs, order_rhs=order_rhs,
                        lhs_coeff=lhs_coeff, rhs_coeff=rhs_coeff,
                        mixed_var=mixed_var, max_lowest_term=max_lowest_term,
                        middle_sign=middle_sign, max_multiple=max_multiple,
                        same_base_root=same_base_root)
                # Attempting to add it
                if self.add_problem(Problem(eq)):
                    return
                # Problem was a dupe, looping back
            # All of the problems generated were dupes, there likely aren't many unique
            # problems for the parameters given
            raise GeneratorError('generic_equation', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            PrintException()


    def add_factorable_expression(self, order=2, max_lowest_term=10, factor_order=1,
            symbols='x', leading_coeff=False, mixed_var=False, len_factor=2):
        '''
        Adds a problem to the given order requiring the factoring of
        a polynominal.

        Arguments:
        order       -   The order of the resulting polynomial.
        factor_order-   The order of the binominal factors of the resulting polynominal.
        symbols     -   a string containing the variables that
                        should be used.
        leading_coeff   -   Set to true to have a leading coeff
                            for the binomial terms that is
                            greater than 0.

        Returns a tuple of expressions. The first expression is still
        factored, the second expression is expanded.
        '''

        try:
            for i in range(self.NUM_ATTEMPTS):
                # Generating expression
                expr = self.gen.gen_factorable_expression(factor_order=factor_order, order=order, leading_coeff=leading_coeff,
                        max_lowest_term=max_lowest_term, symbols=symbols, mixed_var=mixed_var, len_factor=len_factor)
                prob = Problem(expr)
                # Attempting to add it
                if self.add_problem(prob):
                    return
                # Problem was a dupe, looping back
            # All of the problems generated were dupes, there likely aren't many unique
            # problems for the parameters given
            raise GeneratorError('factorable', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            PrintException()

    def add_expandable_expression(self, order=2, max_lowest_term=10, factor_order=1,
            symbols='x', leading_coeff=False, mixed_var=False, len_factor=2):
        '''
        Adds a problem to the given order requiring the expansion of
        a polynominal.

        Arguments:
        order       -   The order of the resulting polynomial.
        factor_order -  The order of the binominal terms.
        symbols     -   a string containing the variables that
                        should be used.
        leading_coeff   -   Set to true to have a leading coeff
                            for the binomial terms that is
                            greater than 0.

        Returns a tuple of expressions. The first expression is still
        factored, the second expression is expanded.
        '''

        try:
            for i in range(self.NUM_ATTEMPTS):
                # Generating expression
                expr = self.gen.gen_factorable_expression(len_factor=len_factor, factor_order=factor_order, order=order, leading_coeff=leading_coeff,
                        max_lowest_term=max_lowest_term, symbols=symbols, mixed_var = mixed_var)
                # swap the reduced and unreduced terms
                temp = expr.unreduced_terms
                expr.unreduced_terms = expr.reduced_terms
                expr.reduced_terms = temp
                # Setting up problem
                prob = Problem(expr)
                # Attempting to add it
                if self.add_problem(prob):
                    return
                # Problem was a dupe, looping back
            # All of the problems generated were dupes, there likely aren't many unique
            # problems for the parameters given
            raise GeneratorError('expandable', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            PrintException()

    def add_linear(self, max_lowest_term=10, max_multiple=1, types='i',
            num_lhs_terms=2, num_rhs_terms=1, lhs_coeff=[], rhs_coeff=[],
            middle_sign='=', mixed_var=False, symbols='x',
            same_base_root=True, order_lhs=1, order_rhs=0):
        '''
        Adds a linear equation to the problem list.

        num_lhs_terms   -   Number of terms on the left hand side of the
                            equation. Default 2.
        num_rhs_terms   -   Number of terms on the right hand side of the
                            equation. Default 1.
        types           -   Types of coefficients. 'i' -> Integers,
                            'f' -> fractions, 'r' -> square roots.
        symbols         -   lists the variables that should be used as a string.
                            Default 'x'.
        lhs_coeff       -   A list of the coefficients for the lhs. Should
                            be of the same length as the number of terms
                            on the left hand side, and will override any
                            automatic number generation. Default [], allowing
                            automatic number generation.
        rhs_coeff       -   A list of the coefficients for the rhs. Should
                            be of the same length as the number of terms
                            on the left hand side, and will override any
                            automatic number generation. Default [], allowing
                            automatic number generation.
        mixed_var       -   boolean determining if variables should be mixed
                            or not (xy vs. x^2 + y^2). Default False.
        max_lowest_term -   the maximum coefficient obtained through automatic
                            coefficient generation. Defautl 10.
        max_multiple    -   Maximum multiple to mulitply coefficients by. This
                            will increase the threshold for automatic numbers
                            beyond max_lowest_term. Default 1.
        same_base_root  -   Boolean determining if radical expressions should have
                            the same base root so they can reduce to a single term.
                            Default True.

        Returns Equation.
        '''

        try:
            for i in range(self.NUM_ATTEMPTS):
                # Generating expression
                eq = self.gen.gen_linear(max_lowest_term=max_lowest_term,
                        max_multiple=max_multiple, types=types,
                        num_lhs_terms=num_lhs_terms, num_rhs_terms=num_rhs_terms,
                        lhs_coeff=lhs_coeff, rhs_coeff=rhs_coeff,
                        middle_sign=middle_sign, mixed_var=mixed_var,
                        symbols=symbols, same_base_root=same_base_root,
                        order_lhs=order_lhs, order_rhs=order_rhs)
                # Attempting to add it
                if self.add_problem(Problem(eq)):
                    return
                # Problem was a dupe, looping back
            # All of the problems generated were dupes, there likely aren't many unique
            # problems for the parameters given
            raise GeneratorError('linear', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            PrintException()

    def add_numerical_expression(self, num_terms=2, op='+-', types='i',
            max_lowest_term=10, max_multiple=1, same_base_root=True):
        '''
        Adds an expression involving a chosen quantity of numbers and their
        solution.

        * Note : all radicals involved in the expression

        Arguments:
        num_terms   -   the number of terms in the Expression. Default 2.
        op          -   a string expressing which operations should be used.
                        String should only contain the characters +, -, *, ^, or
                        /. For example, the string '+/' would indicate the
                        problem should include addition and/or division.
                        Default '+-'.
        types       -   The types of numbers that should be allowed to appear.
                        This is expressed as a string, containing only the
                        characters 'i' (positive integers), 'r' (square roots),
                        'f' (fractions).
                        Default 'i'.
        max_lowest_term
                    -   the largest number that can appear in the reduced
                        expression. (Of course, larger numbers may appear due
                        to the operations being done.)
        max_multiple
                    -   the maximum multiplier used in the creation of fractions
                        and radicals.
        same_base_root
                    -   bool determining if all radicals in the expression
                        should reduce to the same base root.

        Returns an Expression.
        '''

        try:
            for i in range(self.NUM_ATTEMPTS):
                # Generating expression
                expr = self.gen.gen_numerical_expression(num_terms=num_terms,
                        op=op, types=types, max_lowest_term=max_lowest_term,
                        max_multiple=max_mulitple, same_base_root=same_base_root)
                # Attempting to add it
                if self.add_problem(Problem(expr)):
                    return
                # Problem was a dupe, looping back
            # All of the problems generated were dupes, there likely aren't many unique
            # problems for the parameters given
            raise GeneratorError('numerical', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            PrintException()

    # coeffs are generated like max_lowest_term^2, not like max_lowest_term
    # TODO: bug
    def add_quadratic(self, max_lowest_term=4, factorable=True,
            solvable=True, leading_coeff=False, middle_sign='='):
        '''
        Adds a quadratic equation that can be factorable, unfactorable, or unsolvable.

        Arguments:

        max_lowest_term     -   The lowest coefficient that can appear in the problem.
        factorable          -   bool determining if the equation should be factorable.
        solvable            -   bool determining if the equation should be solvable.
        leading_coeff       -   bool determining if a in ax^2 + bx + c should be 1 or not.
                                False means a = 1.
        middle_sign         -   Determines if this is an equation or an inequality.
                                Appropriate arguments are '=', '>=', '<=', '<', '>', '!='.

        Returns an Equation.
        '''

        try:
            for i in range(self.NUM_ATTEMPTS):
                # Generating expression
                eq = self.gen.gen_quadratic(max_lowest_term=max_lowest_term,
                        factorable=factorable, solvable=solvable,
                        leading_coeff=leading_coeff, middle_sign=middle_sign)
                # Attempting to add it
                if self.add_problem(Problem(eq)):
                    return
                # Problem was a dupe, looping back
            # All of the problems generated were dupes, there likely aren't many unique
            # problems for the parameters given
            raise GeneratorError('quadratic', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            PrintException()

    def add_system(self, num_equations=2, num_lhs_terms=2, num_rhs_terms=1,
            types='i', symbols='xy', order_lhs=1, order_rhs=0, lhs_coeff=[],
            rhs_coeff=[], mixed_var=False, max_lowest_term=10, middle_sign='=',
            max_multiple=1, same_base_root=True):
        '''
        Generates a System of Equations involving denoted variables to the
        order specified.

        Arguments:

        num_equations   -   Number of Equations in the System. Default 2.
        num_lhs_terms   -   Number of terms on the left hand side of the
                            equation. Default 2.
        num_rhs_terms   -   Number of terms on the right hand side of the
                            equation. Default 1.
        types           -   Types of coefficients. 'i' -> Integers,
                            'f' -> fractions, 'r' -> square roots.
        order_lhs       -   Order of the left hand side expression.
                            Default 1
        order_rhs       -   Order of the right hand side expression.
                            Default 0.
        lhs_coeff       -   A list of the coefficients for the lhs. Should
                            be of the same length as the number of terms
                            on the left hand side, and will override any
                            automatic number generation. Default [], allowing
                            automatic number generation.
        rhs_coeff       -   A list of the coefficients for the rhs. Should
                            be of the same length as the number of terms
                            on the left hand side, and will override any
                            automatic number generation. Default [], allowing
                            automatic number generation.
        mixed_var       -   boolean determining if variables should be mixed
                            or not (xy vs. x^2 + y^2). Default False.
        max_lowest_term -   the maximum coefficient obtained through automatic
                            coefficient generation. Default 10.
        max_multiple    -   Maximum multiple to mulitply coefficients by. This
                            will increase the threshold for automatic numbers
                            beyond max_lowest_term. Default 1.
        same_base_root  -   Boolean determining if radical expressions should have
                            the same base root so they can reduce to a single term.
                            Default True.
        symbols         -   Variables in the equation. E.g. 'xy' will include x and y
                            terms. Default 'x'.
        '''

        try:
            for i in range(self.NUM_ATTEMPTS):
                # Generating expression
                syst = self.gen.gen_system(num_equations=num_equations,
                        num_lhs_terms=num_lhs_terms, num_rhs_terms=num_rhs_terms,
                        types=types, symbols=symbols, order_lhs=order_lhs,
                        order_rhs=order_rhs, lhs_coeff=lhs_coeff, rhs_coeff=rhs_coeff,
                        mixed_var=mixed_var, max_lowest_term=max_lowest_term,
                        middle_sign=middle_sign, max_multiple=max_multiple,
                        same_base_root=same_base_root)
                # Attempting to add it
                if self.add_problem(Problem(syst)):
                    return
                # Problem was a dupe, looping back
                # All of the problems generated were dupes, likely aren't many unique
                # problems for the parameters given
                raise GeneratorError('system', 'Unable to generate additional ' +
                        'unique problems after trying ' + str(self.NUM_ATTEMPTS) +
                        ' times.' + 'Input parameters may be too restrictive.')
        except GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            PrintException()

class Worksheet(ProblemContainer):
    '''
    Class for creating a worksheet.

    Member variables:

    worksheet_fn: Filename of the worksheet file. Should end in .tex.
    output_fn: Filename of output pdf.
    title:  Title of the worksheet. should be a string.
    message: Message to be displayed at the beginning of the worksheet.
    author: Author of the worksheet. Should be a string.

    Member variables from ProblemContainer (parent class)

    gen: Generator used to create problems.
    problems: list of Problems.
    '''

    def __init__(self, worksheet_fn):
        assert type(worksheet_fn) == str
        ProblemContainer.__init__(self)
        self.worksheet_fn = worksheet_fn
        self.title = ''
        self.author = ''
        self.message = ''
        self.output_fn = ''

    def set_title(self, title):
        assert type(title) == str
        self.title = title

    def set_author(self, author):
        assert type(author) == str
        self.author = author

    def set_message(self, message):
        assert type(message) == str
        self.message = message

    # TODO: improve line break errors for fields in \text{}
    def make(self, num_cols=2):
        '''
        Takes a predefined latex template, an author, a title, and a list of
        problems with their solutions and generates a latex worksheet.

        Returns nothing.
        '''
        assert num_cols >= 1
        # Opening template
        if num_cols == 1:
            template = TEMPLATE1COL
        else:
            template = TEMPLATE

        # Formatting problems to fit into a latex enumerate environment

        question_str = ''
        solution_str = ''
        for p in self.problems:
            q = self.make_line_breaks(p.latex_question, 200)
            s = self.make_line_breaks(p.latex_solution, 200)
            question_str += '\\item $ %s $\n \\vspace{10mm}\n' % q
            solution_str += '\\item $ %s $\n \\vspace{10mm}\n' % s

        # Creating author and title strings
        title_str = '\\chead{\\textbf{\\LARGE %s }}\n' % self.title
        author_str = '\\rhead{%s}\n' % self.author

        # Subbing strings into template to make the worksheet. The template
        # should have 5 separate %s characters marking the locations of each
        # of these substitutions in order, and a %d character for the number of columns
        worksheet = template % (num_cols, title_str, author_str, self.message, question_str,
                solution_str)

        # Saving worksheet
        worksheet_file = open(self.worksheet_fn, 'w')
        worksheet_file.write(worksheet)
        worksheet_file.close()

        # Compiling worksheet
        os.system("pdflatex %s" % self.worksheet_fn)

        # Cleaning files and organizing
        os.system("rm *.aux *.log")
        worksheet_dir = 'worksheets'
        tex_dir = 'tex'
        if not os.path.exists(worksheet_dir):
            os.makedirs(worksheet_dir)
        if not os.path.exists(tex_dir):
            os.makedirs(tex_dir)
        os.system("mv " + self.worksheet_fn + ' ' +  tex_dir + '/' + self.worksheet_fn)
        output_pdf = self.worksheet_fn.replace('.tex', '.pdf')
        os.system("mv " + output_pdf + ' ' +  worksheet_dir + '/' + output_pdf)

        # Saving worksheet name
        self.output_fn = worksheet_dir + '/' + output_pdf

    def make_line_breaks(self, string, num_char=70):
        '''
        Makes line breaks every num_char characters, returns formatted string.
        '''
        current_line_len = len(string)
        break_loc = num_char - 1
        while current_line_len > num_char:
            string = string[:break_loc] + '\\\\' + string[break_loc:]
            break_loc += num_char
            current_line_len -= num_char
        return string

    def show(self):
        '''
        Uses default pdf viewer to display a filename.

        Arguments:

        filename :      filename of the file to display.

        Returns nothing.
        '''
        os.system('xdg-open "%s"' % self.output_fn)
# Latex template
TEMPLATE = '''
\\documentclass[11pt]{article}

\\usepackage{graphicx}
\\usepackage{seqsplit}
\\usepackage{caption}
\\usepackage{subcaption}
\\usepackage{float}
\\usepackage{fancyhdr}
\\usepackage[utf8]{inputenc}
\\usepackage[english]{babel}
\\usepackage{lastpage}
\\usepackage{multicol}
\\usepackage{blindtext}
\\usepackage{amstext}


\\pagestyle{fancy}
\\fancyhf{}

\\rfoot{Page \\thepage}
\\begin{document}

\\newcommand{\\numcols}{%d}
\\lhead{}
%s
%s
%s
\\begin{multicols}{\\numcols{}}[\section*{Problems}]
\\begin{enumerate}
%s
\\end{enumerate}
\\end{multicols}
\\newpage
\\begin{multicols}{\\numcols{}}[\section*{Answers}]
\\begin{enumerate}
%s
\\end{enumerate}
\\end{multicols}

\\end{document}
'''

TEMPLATE1COL = '''
\\documentclass[11pt]{article}

\\usepackage{graphicx}
\\usepackage{caption}
\\usepackage{seqsplit}
\\usepackage{subcaption}
\\usepackage{float}
\\usepackage{fancyhdr}
\\usepackage[utf8]{inputenc}
\\usepackage[english]{babel}
\\usepackage{lastpage}
\\usepackage{multicol}
\\usepackage{blindtext}
\\usepackage{amstext}


\\pagestyle{fancy}
\\fancyhf{}

\\rfoot{Page \\thepage}
\\begin{document}

\\newcommand{\\numcols}{%d}
\\lhead{}
%s
%s
%s
\section*{Problems}
\\begin{enumerate}
%s
\\end{enumerate}
\\newpage
\section*{Answers}
\\begin{enumerate}
%s
\\end{enumerate}

\\end{document}
'''

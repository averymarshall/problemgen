import problemgen.backend as backend
import random
import os
import subprocess

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
        self.gen = backend.Generator()
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
                if self.add_problem(backend.Problem(expr)):
                    return
                # Problem was a dupe, looping back
            # All of the problems generated were dupes, there likely aren't many unique
            # problems for the parameters given
            raise backend.GeneratorError('algebraic', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except backend.GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            backend.PrintException()

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
            raise backend.GeneratorError('num_conv', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except backend.GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            backend.PrintException()
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
            raise backend.GeneratorError('dec_to_frac', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except backend.GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            backend.PrintException()

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
            raise backend.GeneratorError('frac_to_dec', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except backend.GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            backend.PrintException()

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
                if self.add_problem(backend.Problem(eq)):
                    return
                # Problem was a dupe, looping back
            # All of the problems generated were dupes, there likely aren't many unique
            # problems for the parameters given
            raise backend.GeneratorError('generic_equation', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except backend.GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            backend.PrintException()


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
                prob = backend.Problem(expr)
                # Attempting to add it
                if self.add_problem(prob):
                    return
                # Problem was a dupe, looping back
            # All of the problems generated were dupes, there likely aren't many unique
            # problems for the parameters given
            raise backend.GeneratorError('factorable', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except backend.GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            backend.PrintException()

# TODO: Sometimes this generates monomials, strange behavior
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
                prob = backend.Problem(expr)
                # Attempting to add it
                if self.add_problem(prob):
                    return
                # Problem was a dupe, looping back
            # All of the problems generated were dupes, there likely aren't many unique
            # problems for the parameters given
            raise backend.GeneratorError('expandable', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except backend.GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            backend.PrintException()

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
                if self.add_problem(backend.Problem(eq)):
                    return
                # Problem was a dupe, looping back
            # All of the problems generated were dupes, there likely aren't many unique
            # problems for the parameters given
            raise backend.GeneratorError('linear', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except backend.GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            backend.PrintException()

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
                if self.add_problem(backend.Problem(expr)):
                    return
                # Problem was a dupe, looping back
            # All of the problems generated were dupes, there likely aren't many unique
            # problems for the parameters given
            raise backend.GeneratorError('numerical', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except backend.GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            backend.PrintException()

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
                if self.add_problem(backend.Problem(eq)):
                    return
                # Problem was a dupe, looping back
            # All of the problems generated were dupes, there likely aren't many unique
            # problems for the parameters given
            raise backend.GeneratorError('quadratic', 'Unable to generate additional ' +
                    'unique problems after trying ' + str(self.NUM_ATTEMPTS) + ' times.' +
                    'Your input parameters may be too restrictive.')
        except backend.GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            backend.PrintException()

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
                if self.add_problem(backend.Problem(syst)):
                    return
                # Problem was a dupe, looping back
                # All of the problems generated were dupes, likely aren't many unique
                # problems for the parameters given
                raise backend.GeneratorError('system', 'Unable to generate additional ' +
                        'unique problems after trying ' + str(self.NUM_ATTEMPTS) +
                        ' times.' + 'Input parameters may be too restrictive.')
        except backend.GeneratorError as e:
            print('GeneratorError: %s' % e.message)
        except:
            backend.PrintException()

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
        self.output_fn_no_answers = ''

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
    # TODO: Use PyLatex to avoid external dependencies
    def make(self, num_cols=2, separate_answers=True):
        '''
        Takes a predefined latex template, an author, a title, and a list of
        problems with their solutions and generates a latex worksheet.

        num_cols: number of columns in the worksheet.
        separate_answers: boolean describing if answers should be generated in
        a separate worksheet.

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

        filename = self.worksheet_fn[:-4] + '--with-answers' + self.worksheet_fn[-4:]
        # Saving worksheet
        worksheet_file = open(filename, 'w')
        worksheet_file.write(worksheet)
        worksheet_file.close()

        # Compiling worksheet
        try:
            subprocess.check_output("pdflatex %s" % filename, stderr=subprocess.STDOUT)
            # os.system("pdflatex %s" % filename) 
            # Cleaning files and organizing
            os.system("rm *.aux *.log")
            worksheet_dir = 'worksheets'
            tex_dir = 'tex'
            if not os.path.exists(worksheet_dir):
                os.makedirs(worksheet_dir)
            if not os.path.exists(tex_dir):
                os.makedirs(tex_dir)
            os.system("mv " + filename + ' ' +  tex_dir + '/' + filename)
            output_pdf = filename.replace('.tex', '.pdf')
            os.system("mv " + output_pdf + ' ' +  worksheet_dir + '/' + output_pdf)
        except (OSError, IOError, subprocess.CalledProcessError) as e:
            print e
            

        # Saving worksheet name
        self.output_fn = worksheet_dir + '/' + output_pdf

        if separate_answers:
            # Generating another sheet with no Answers

            worksheet = TEMPLATE_NO_ANSWERS % (num_cols, title_str, author_str, self.message, question_str)

            filename = self.worksheet_fn[:-4] + '--without-answers' + self.worksheet_fn[-4:]
            # Saving worksheet
            worksheet_file = open(filename, 'w')
            worksheet_file.write(worksheet)
            worksheet_file.close()

            try:
                subprocess.check_output("pdflatex %s" % filename, stderr=subprocess.STDOUT)
                # os.system("pdflatex %s" % filename) 
                # Cleaning files and organizing
                os.system("rm *.aux *.log")
                worksheet_dir = 'worksheets'
                tex_dir = 'tex'
                if not os.path.exists(worksheet_dir):
                    os.makedirs(worksheet_dir)
                if not os.path.exists(tex_dir):
                    os.makedirs(tex_dir)
                os.system("mv " + filename + ' ' +  tex_dir + '/' + filename)
                output_pdf = filename.replace('.tex', '.pdf')
                os.system("mv " + output_pdf + ' ' +  worksheet_dir + '/' + output_pdf)
            except (OSError, IOError, subprocess.CalledProcessError) as e:
                print e

            self.output_fn_no_answers = worksheet_dir + '/' + output_pdf

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

TEMPLATE_NO_ANSWERS = '''
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

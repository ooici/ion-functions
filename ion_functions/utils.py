#!/usr/bin/env python

"""
@package ion_functions.qc_functions
@file ion_functions/qc_functions.py
@author Christopher Mueller
@brief Module containing helper functions, ported from matlab, as described in DPS documents
"""

import numpy as np

ALL_KINDS = ('i', 'u', 'f', 'c', 'S', 'a', 'U')  # Does not include 'V' which is raw (void) or O which is object
NUMERIC_KINDS = ('i', 'u', 'f', 'c')
REAL_KINDS = ('i', 'u', 'f', 'S', 'a', 'U')  # All kinds but complex


def isnumeric(dat):
    """
    isnumeric - Determine whether input is numeric
    Syntax
    tf = isnumeric(A)
    Description
    tf = isnumeric(A) returns logical 1 (true) if A is a numeric array and logical 0 (false)
    otherwise. For example, sparse arrays and double-precision arrays are numeric, while strings,
    cell arrays, and structure arrays and logicals are not.
    Examples
    Given the following cell array,
    C{1,1} = pi;                 % double
    C{1,2} = 'John Doe';         % char array
    C{1,3} = 2 + 4i;             % complex double
    C{1,4} = ispc;               % logical
    C{1,5} = magic(3)            % double array
    C =
       [3.1416] 'John Doe' [2.0000+ 4.0000i] [1][3x3 double]
    isnumeric shows that all but C{1,2} and C{1,4} are numeric arrays.
    for k = 1:5
    x(k) = isnumeric(C{1,k});
    end
    x
    x =
         1     0     1     0     1

    """

    return np.array([np.atleast_1d(d).dtype.kind in NUMERIC_KINDS for d in np.nditer(np.atleast_1d(dat))]).astype('int8')


def isreal(dat):
    """
    isreal - Check if input is real array
    Syntax
    TF = isreal(A)
    Description
    TF = isreal(A) returns logical 1 (true) if A does not have an imaginary part. It returns logical
    0 (false) otherwise. If A has a stored imaginary part of value 0, isreal(A) returns logical 0
    (false).
    Note   For logical and char data classes, isreal always returns true. For numeric data
    types, if A does not have an imaginary part isreal returns true; if A does have an imaginary
    part isreal returns false. For cell, struct, function_handle, and object data types,
    isreal always returns false.
    ~isreal(x) returns true for arrays that have at least one element with an imaginary
    component. The value of that component can be 0.
    Tips
    If A is real, complex(A) returns a complex number whose imaginary component is 0, and
    isreal(complex(A)) returns false. In contrast, the addition A + 0i returns the real value A,
    and isreal(A + 0i) returns true.
    If B is real and A = complex(B), then A is a complex matrix and isreal(A) returns false,
    while A(m:n) returns a real matrix and isreal(A(m:n)) returns true.
    Because MATLAB software supports complex arithmetic, certain of its functions can introduce
    significant imaginary components during the course of calculations that appear to be limited to
    real numbers. Thus, you should use isreal with discretion.
    Example 1
    If a computation results in a zero-value imaginary component, isreal returns true.
    x=3+4i;
    y=5-4i;
    isreal(x+y)
    ans =
         1
    Example 2
    These examples use isreal to detect the presence or absence of imaginary numbers in an
    array. Let
    x = magic(3);
    y = complex(x);
    isreal(x) returns true because no element of x has an imaginary component.
    isreal(x)
    ans =
         1
    isreal(y) returns false, because every element of x has an imaginary component, even
    though the value of the imaginary components is 0.
    isreal(y)
    ans =
         0
    This expression detects strictly real arrays, i.e., elements with 0-valued imaginary components
    are treated as real.
    ~any(imag(y(:)))
    ans =
         1
    Example 3
    Given the following cell array,
    C{1} = pi;                 % double
    C{2} = 'John Doe';         % char array
    C{3} = 2 + 4i;             % complex double
    C{4} = ispc;               % logical
    C{5} = magic(3);           % double array
    C{6} = complex(5,0)        % complex double
    C =
      [3.1416]  'John Doe'  [2.0000+ 4.0000i]  [1]  [3x3 double]  [5]
    isreal shows that all but C{1,3} and C{1,6} are real arrays.
    for k = 1:6
    x(k) = isreal(C{k});
    end
    x
    x =
         1     1     0     1     1    0
    """

    return np.array([np.atleast_1d(d).dtype.kind in REAL_KINDS for d in np.nditer(np.atleast_1d(dat))]).astype('int8')


def isscalar(dat):
    """
    isscalar - Determine whether input is scalar
    Syntax
    isscalar(A)
    Description
    isscalar(A) returns logical 1 (true) if size(A) returns [1 1], and logical 0 (false) otherwise.
    Examples
    Test matrix A and one element of the matrix:
    A = rand(5);
    isscalar(A)
    ans =
         0
    isscalar(A(3,2))
    ans =
         1


    """
    return np.atleast_1d(dat).size == 1


def isvector(dat):
    """
    isvector - Determine whether input is vector
    Syntax
    isvector(A)
    Description
    isvector(A) returns logical 1 (true) if size(A) returns [1 n] or [n 1] with a nonnegative
    integer value n, and logical 0 (false) otherwise.
    Examples
    Test matrix A and its row and column vectors:
    A = rand(5);
    isvector(A)
    ans =
         0
    isvector(A(3, :))
    ans =
         1
    isvector(A(:, 2))
    ans =
         1


    """
    return np.atleast_1d(dat).size > 1


def ismatrix(dat):
    """
    ismatrix - test if input array is formatted as a matrix
    
    Syntax
    
        flag = ismatrix(dat)
        
    Description
    
        ismatrix(dat) returns logical 1 (true) if np.atleast_1d(dat).shape
        returns (m, n) with nonnegative integer values m and n, and logical 0
        (false) otherwise.

    Examples
        
        dat = np.array([0, 1]);
        ismatrix(dat)
        0
        dat = np.array([[0,1]])
        ismatrix(dat)
        1
        
    """
    return len(np.atleast_1d(dat).shape) == 2


def isempty(dat):
    """
    isempty - Test if array is empty
    Syntax
    tf = isempty(A)
    Description
    tf = isempty(A) returns logical true (1) if A is an empty array and logical false (0) otherwise.
    An empty array has at least one dimension of size zero, for example, 0-by-0 or 0-by-5.
    Examples
    B = rand(2,2,2);
    B(:,:,:) = [];
    isempty(B)
    ans =
         1
    """

    return np.atleast_1d(dat).size == 0


def islogical(inflags):
    """
    islogical - test if input array is a boolean array (all values are either a
    "0" or "1")
    
    Syntax
    
        flag = islogical(inflags)
        
    Description
    
        ismatrix(inflags) returns logical 1 (true) if all values in the inflag
        array are either a 0 or a 1 with the type encoding set to an int8

    Examples
        
        inflags = np.array([0, 1]).astype('np.int8')
        islogical(inflags)
        True
        inflags = np.array([[0, 2]]).astype('np.int8')
        islogical(inflags)
        False
        inflags = np.array([[0, 1]]).astype('np.float')
        islogical(inflags)
        False
        
    """
    inflags = np.atleast_1d(inflags)
    flag = (all(np.in1d(inflags.flatten(), [0,1])) and
            all(isinstance(n, np.int8) for n in inflags.flatten()))
    return flag

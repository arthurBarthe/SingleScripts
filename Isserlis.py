# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 10:59:32 2018

@author: Arthur Guillaumin
This script derives the partitions of a set into pairs. It can be used to
derive the Isserlis formula.


Example: we obtain that (s() autocovariance of the stationary Gaussian process
X_t):
E{X_{t}^2 X_{t+t_1}^2 X_{t+t_2}^2 X_{t+t_3}^2}
=
s(0)^4  
+  2*s(0)^2s(t_3-t_1)^2  
+  2*s(0)^2s(t_2-t_1)^2  
+  2*s(0)^2s(t_3-t_2)^2  
+  2*s(0)^2s(t_1)^2  
+  2*s(0)^2s(t_3)^2  
+  2*s(0)^2s(t_2)^2  
+  4*s(t_2)^2s(t_3-t_1)^2  
+  4*s(t_2-t_1)^2s(t_3)^2  
+  4*s(t_1)^2s(t_3-t_2)^2  
+  8*s(0)s(t_1)s(t_3-t_1)s(t_3)  
+  8*s(0)s(t_2-t_1)s(t_3-t_2)s(t_3-t_1)  
+  8*s(0)s(t_1)s(t_2-t_1)s(t_2)  
+  8*s(0)s(t_2)s(t_3-t_2)s(t_3)  
+  16*s(t_1)s(t_2)s(t_3-t_2)s(t_3-t_1)  
+  16*s(t_2-t_1)s(t_2)s(t_3-t_1)s(t_3)  
+  16*s(t_1)s(t_2-t_1)s(t_3-t_2)s(t_3) 
"""
from itertools import permutations
from collections import OrderedDict

def isserlis(set_):
    """This function takes a list of elements which represent a list of
    variables. It then returns a list of lists, where the ordering of the 
    inner list is to be interpreted as a way to partition into pairs.
    For instance if set_ = [1,1,2,3], the [1,2,1,3] is to be interpreted as 
    the partition {{1,2}, {1,3}}."""
    n_elements = len(set_)
    #We shall go through the permutations an keep only a fraction of them.
    #This is far from the fastest way but works well for not too large sets
    #of variables.
    def test(n_els, t):
        i = t.index(0)
        for j in range(0,n_els,2):
            temp = t.index(j)
            if temp < i or temp > t.index(j+1):
                return False
            else:
                i = temp
        return True
    for p in permutations(range(n_elements)):
        if test(n_elements, p):
            yield([set_[p.index(i)] for i in range(n_elements)])

def isserlis_stationary(set_):
    """This function derives Isserlis for a stationary time series. This
    generator yields elements of the final sum, where each element is a list
    to be interpreted as follows: 
    (1,0,1) -> s(0)s(2)
    (1,0,2) -> s(0)s(2)^2
    (0,1,0,10) -> s(1)s(3)^10"""
    n_elements = len(set_)
    for el in isserlis(set_):
        l = [abs(el[j+1]-el[j]) for j in range(0,n_elements, 2)]
        yield tuple([l.count(i) for i in range(max(l)+1)])

def isserlis_stationary_analytical(set_):
    """This function allows to derive expectations such as 
    E{X_{t}^2 X_{t+t_1}^2 X_{t+t_2}^2 X_{t+t_3}^2} (see example lower).
    For this it converts each unique string of the list set_ into an integer,
    in a way that ensures that differences between two integers are uniquely
    related to those integers. For instance, [0, 1, 3, 7, 15, 31, 63, ...]"""
    #First we get unique strings
    set_elements = set(set_)
    #Then we build a unique correspondance.
    correspondances = {el: 2**k-1 for k,el in enumerate(set_elements)}
    set_to_numbers = [correspondances[el] for el in set_]
    return set_to_numbers, isserlis_stationary(set_to_numbers)

def counts(l):
    """This function takes a list as an argument and returns a count dictionary
    of the elements in the list."""
    l = tuple(l)
    return {j:l.count(j) for j in l}

def beautiful_print(set_, dic, option="Screen"):
    if option == ".tex":
        times_ = '\\ '
        endline_ = '\\\\'
        newline_ = '&'
        nonumber_ = '\\nonumber'
        s = ['', '\\tau_1', '\\tau_2', '\\tau_3', '\\tau_4', '\\tau_5']
        s2 = ['', '-\\tau_1', '-\\tau_2', '-\\tau_3', '-\\tau_4', '-\\tau_5']
        cov = 's_X'
    else:
        times_ = '*'
        endline_ = ''
        newline_ = ''
        s = ['', 't_1', 't_2', 't_3', 't_4', 't_5']
        s2 = ['', '-t_1', '-t_2', '-t_3', 't_4', 't_5']
        cov = 's'
        nonumber_ = ''
    set_ = list(set(set_))
    differences = {}
    for k, i in enumerate(set_):
        for l,j in enumerate(set_[k+1:]):
            differences[abs(j-i)] = cov + '(' + s[k+l+1] + s2[k] + ')'
    differences[0] = cov +'(0)'
    dic = OrderedDict(sorted(dic.items(), key=lambda t: t[1]))
    for j,k in enumerate(dic.keys()):
        txt = ''
        if dic[k] != 1:
            txt += str(dic[k]) + times_
        for _,l in enumerate(k):
            if l != 0:
                if l==1:
                    txt += differences[_]
                else:
                    txt += differences[_] + '^' + str(l)
        if j==0:
            print(newline_, txt, endline_,  nonumber_)
        else:
            print('+', newline_, txt, endline_, nonumber_)

#set_ = [0, 0, 1, 1, 3, 3, 7,7]
#Here we will get E{X_{t}^2 X_{t+t_1}^2 X_{t+t_2}^2 X_{t+t_3}^2}
set_analytic = ['t', 't', 't+t_1', 't+t_1', 't+t_2', 't+t_2',
                't+t_3', 't+t_3']#, 't+t_2', 't+t_3']
#print(list(isserlis(set_)), '\n')
#print(list(isserlis_stationary(set_)), '\n')
#iss = counts(isserlis_stationary(set_))

#beautiful_print(set_, iss)

set_to_nbs, iss = isserlis_stationary_analytical(set_analytic)
iss = counts(iss)
beautiful_print(set_to_nbs, iss)


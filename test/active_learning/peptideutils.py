# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 23:24:20 2019

@author: Alexander van Teijlingen
"""

peptideutils_letters1 = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
peptideutils_letters3 = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HSE', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
    
def translate1to3(string):
    global peptideutils_letters1
    global peptideutils_letters3
    code = list(string)
    new_string = ""
    for letter in code:
        index = peptideutils_letters1.index(letter)
        new_string = new_string + peptideutils_letters3[index] + "-"
    new_string = new_string[:-1]
    return new_string

def translate3to1(string):
    global peptideutils_letters1
    global peptideutils_letters3
    code = string.split("-")
    new_string = ""
    for AA in code:
        if AA == "HIS":
            AA = "HSE"
        index = peptideutils_letters3.index(AA)
        new_string = new_string + peptideutils_letters1[index]
    return new_string


#!/usr/bin/python3
s = open('rosalind_dna.txt','r').read()
print(s.count("A"), s.count("C"), s.count("G"), s.count("T"))

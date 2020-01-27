# Rewrite from Liang's original code in Python 2.5

# 1, just use 2to3 *.py to convert python2 code to python3 code
1, print xxx  -> print(xxx)
2, not equal to, <> (py2.5) -> != (py3)

# 2, auto reformat, fix indentation issues
autopep8 --in-place --aggressive --aggressive bpch2_rw_v2.py


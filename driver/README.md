# Adapted from Liang's code in Python 2.5, circa 2007

### 1, use 2to3 xxx.py to convert python2 code to python3 code
Solves issues like these:

1, print xxx  -> print(xxx)

2, not equal to, <> (py2.5) -> != (py3)

### 2, auto reformat, fix indentation issues

autopep8 --in-place --aggressive --aggressive xxx.py

### 3, replace outdated modules with new ones

E.g., basemap -> cartopy

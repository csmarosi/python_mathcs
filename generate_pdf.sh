#!/bin/sh

rm hw_py.tex
./continuum_mechanics_hw2.py > hw_py.tex
rm hw.tex
cat latex_header.tex hw_py.tex > hw.tex
pdflatex -halt-on-error hw.tex


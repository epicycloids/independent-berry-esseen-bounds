"""Interval bounds for the variance range in the Gaussian mixture comparison."""
from fractions import Fraction as Q

import hashlib

import json

from pathlib import Path

from flint import arb,arb_series,ctx

import flint

def number(q):
    q=Q(q)
    return arb(q.numerator)/q.denominator


def interval(lo,hi):
    m=(lo+hi)/2;r=(hi-lo)/2
    return arb(f"{m.numerator}/{m.denominator}",f"{r.numerator}/{r.denominator}")


def function(x,kind):
    return 1/(x*x.sqrt()) if kind=="power" else -x.log()


def derivative(x,kind):
    return -arb(3)/(2*x*x*x.sqrt()) if kind=="power" else -1/x


def chord(l,d,kind):
    a=1-l;b=1+d-l;r=d-l
    return (r*function(a,kind)+l*function(b,kind))/d-function(1+0*l,kind)


def chord_derivative(l,d,kind):
    a=1-l;b=1+d-l;r=d-l
    return (function(b,kind)-function(a,kind)
            -r*derivative(a,kind)-l*derivative(b,kind))/d


def upper(lo,hi,d,kind):
    mid=(lo+hi)/2
    g=chord_derivative(interval(lo,hi),number(d),kind)
    magnitude=max(abs(g.lower()),abs(g.upper()))
    value=chord(number(mid),number(d),kind)+number((hi-lo)/2)*magnitude
    return value.upper()



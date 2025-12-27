# TSSpy package init
"""
TSSpy: Python CLI for TSS Analysis

A comprehensive command-line tool for transcription start site (TSS) data analysis.
Inspired by TSSr (R/Bioconductor) but implemented in pure Python.
"""

__version__ = "0.3.0"
__author__ = "Johnny Chen"

from TSSpy.main import app, main

__all__ = [
    'app',
    'main',
    '__version__',
]

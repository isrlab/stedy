# Contributing guidelines

## How to contribute

This project accepts contributions from anyone interested in improving it. To do so, follow these steps:

1. [Fork it](https://github.com/uqLab/stedy/fork)
2. Improve it, respecting the <a href="#codstand">coding standards</a> and
using proper [commit messages](https://chris.beams.io/posts/git-commit/)
3. Submit a [pull request](https://help.github.com/articles/creating-a-pull-request)

Your work will then be reviewed as soon as possible (suggestions about some
changes, improvements or alternatives may be given).

<a name="codstand" />

## Coding standards

* Code must work in [MATLAB] (>= 2017a).
* MATLAB-style documentation and comments.
* Encoding: UTF-8.
* Newlines: Unix style, i.e. LF or sprintf('\n').

## Issues and support

Problems with this software can be reported in the [issues section](https://github.com/uqLab/stedy/issues).
Support for *stedy* is provided on best effort basis by emailing the author at vaishnavtv@tamu.edu.

## Test Case
2 test cases ([*TestPendulum.m*] and [*TestTBar.m*]) has been provided in the [Tests] folder of the repository for the user to test their version of the code.
On running the files, the results of the three checks - Bar Length violations, total energy violations, and motion errors - will be displayed on the command window.

Please note that the main contributions of the software comprise an efficient Lagrangian formulation of the dynamics and a novel constraint correction method. Improvements to plotting functions are not the focus of these test cases.

#### Note

This document is partially adapted from the [TempoSimple](https://github.com/gnugat-legacy/tempo-simple/blob/master/CONTRIBUTING.md)
contributing guidelines.

[MATLAB]: http://www.mathworks.com/products/matlab/
[Tests]: https://github.com/uqLab/stedy/tree/master/Tests
[*TestPendulum.m*]: https://github.com/uqLab/stedy/blob/master/Tests/TestPendulum.m
[*TestTBar.m*]: https://github.com/uqLab/stedy/blob/master/Tests/TestTBar.m

MateMD
======

Molecular Dynamics code in python with c bindings.

This program is inspired on
[CoffeMD](https://github.com/pabloalcain/coffeemd), work done during
the 2014 ICTP [Workshop on Advanced Techniques for Scientific
Programming and Management of Open Source Software
Packages](http://cdsagenda5.ictp.it/full_display.php?ida=a13190),
which is in turned based on the original c code from Axel Kohlmeyer,
[ljmd-c](https://github.com/akohlmey/ljmd-c) and refactored to add
python bindings. It is distributed under the terms of the GNU General
Public License.

In this new work, we'll explicitly make everything from Python and
just then rewrite in c the hard math. The code is object oriented in
its design, meant to be easily integrated to Python scripts.

Introduction
------------

As with CoffeeMD, the main goal of this program is to *build*
molecular dynamics simulation with the flexibility that Python gives,
leaving the hard math to c. Unlike CoffeMD, however, this code will be
object oriented from Python.
earth/src/tests/README
----------------------

The tests in this directory are much more comprehensive than
tests\test.earth.R but also take much longer to execute.

For regression tests on earth.c with the Microsoft compiler: run test.earthc.bat
For regression tests on earth.c with gcc: run test.earthc.bat
For regression tests of earth with a variety of models: run test.earth.full.bat
For regression tests of earth with factors and glms: run test.earth.glm.bat
For regression tests of earth with a big model: run test.earth.big.bat
For regression tests on plotmo: run test.plotmo.bat
For some old timing test results: see timing-tests.txt
For timing tests (as shown on the earth web page): run earth.times.bat

None of these do much testing of appropriate behaviour for bad arguments.

See earth-data-structures.jpg for a description of earth's 
internal data structures.

Stephen Milborrow Apr 2007 Petaluma

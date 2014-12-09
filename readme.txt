Spektro - a small tool to assist in analysis of mass-spectrometry run results

See http://code.google.com/p/ouspg/wiki/Spektro for more details.

GLOBAL INSTALL:
  $ tar -zxvf spektro-0.1.tar.gz
  $ cd spektro-0.1
  $ make
  $ sudo make install

I JUST WANT TO SOLVE ONE THING REAL QUICK WITHOUT INSTALL:
  $ cd spektro-0.1
  $ cc -O2 -o spektro spektro.c
  # have a look at examples from the above mentioned web page
  # and run ./spektro <your data>


tetemoko
========

Tetemoko is a dynamic rupture code based on the Chombo AMR framework. It is
authored by Jeremy Edward Kozdon <jkozdon@stanford.edu>.

Chombo 3.1 is a prerequisite for to use this package and can be obtained from
    commons.lbl.gov/display/chombo/
Tetemoko is currently configure to run only with version 3.1.

First make sure that you can build Chombo which will require setting up your
environment according the Chombo installation instructions. Once this is done
(i.e. you can compile the library) go the Chombo-3.1/lib/src and run the command

```
  cd ${CHOMBO_PATH}/lib/src
  patch -p0 < ${TETEMOKO_PATH}/tetemoko_chombo_patchfile
```
which will make the necessary modifications to the Chombo sources files to all
you to use tetemoko.

Modify the makefile 

```
  ${TETEMOKO_PATH}/execLinElast/GNUmakefile
```

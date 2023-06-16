Small test cases for TsunAWI
============================

Usage
-----
* First clone or download and build the [TsunAWI source code](https://gitlab.awi.de/tsunawi/tsunawi). From the root 
  directory tsuawi/, clone or download this verification path. If you put it into another location, fine, but 
  remember to adjust the exepath, see below.
  
* Check ./test_tsunawi.sh if
  ```
   exe=ompTsuna.x  
   exepath=$(pwd)/../src  
   runcommand=""  
  ``` 
  fit to your system, adjust if necessary. An example for a runcommand on AWI's cluster with Slurm is included.

  If you want to employ OpenMP, set OMP_NUM_THREADS or adjust the runcommand according to your system.

* Start the first run to produce first checksums with
  ``` 
  $ ./test_tsunawi -r 
  ```

* Work on the TsunAWI code, and check if the checksums are reproducable after your changes:
  ```
  $ ./test_tsunawi -c      # no output, only in case of changed checksums
  $ ./test_tsunawi -cv     # a bit verbose
  $ ./test_tsunawi -cV2    # with TsunAWI output - helpful if the code breaks
  ```
* You can also run each test individually, calling TsunAWI directly, e.g.,
  ```
   $ cd bm2_okushiri  
   $ ../../src/ompTsuna.x
  ```

Test cases
----------

All test cases are kept short with a very reduced number of timesteps. Also the 
setups itself are chosen with regard to small mesh size.
IO is (mostly) switched off, only checksums are compared, netcdf or XML support
is not needed. 
However, each directory contains a namelist that can be adjusted for more
realistic setups: increase the number of timesteps, switch on IO.

* **bm1_sloping_beach**
  Synthetic benchmark for the runup on a plain beach
  https://nctr.pmel.noaa.gov/benchmark/Analytical/index.html

* **bm2_okushiri**
  Channel experiment simulating the 1983 tsunami that heavily struck the Monai Valley on the Japanese island of Okushiri.
  https://nctr.pmel.noaa.gov/benchmark/Laboratory/Laboratory_MonaiValley/index.html
* **global_cosinebell**
  A very coarse global mesh (far too coarse for realistic simulations!), 
  initialized with a cosine bell at the approximate location of the Lissabon 1755 earthquake and tsunami,
  though with increased magnitude of 9.0.
* **global_cosinebell_restart**
  As above, but the simulation is split half way and restarted once, to test this functionality.

Adding test cases
-----------------

Basically, just add a directory with a namelist.tsunami, with write_final_checksums=.true. switched on.
The script test_tsunawi.sh will automatically enter this directory and run TsunAWI. 

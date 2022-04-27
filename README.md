Unispring
=================
Unispring implementation with uniform distribution.

Dependencies
--------
python 3.5+ with :
- scipy
- numpy
- python-osc

How to launch
----------
Launch the osc server in your terminal:

	python /*gitDirectory*/Auto-unispring/Python/osc-server-unispring.py
	
Launch the Max patch : catart_unispring.maxpat

How to use
----------
1) Import sounds (file by file or whole folder)
2) Create the distribution
3) Change the region using the drawing area point by point
4) Close the region and fit the distribution to the region

Known bug
----------
udpsend limit : if a buffer has more than ~900 rows (depending on system buffersize), udpsend start dropping packets.

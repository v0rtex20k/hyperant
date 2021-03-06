## Author: Victor Arsenescu
## Date: 10/18/21
## Project Title: Boston Medical Center (BMC)
## Research Programmer Technical Test – "Ant on a Cube"

## To run from barebones repo:
### `$ python3 ui_hypercube.py` // for GUI version
### `$ ./scripts/run_hypercube <n_dimensions> <start node idx> <end node idx> <n_threads>`
### // headless - e.x `./scripts/run_hypercube 3 0 7 50`
### Remember, nodes are zero indexed!!

----------------------------->>> Introduction <<--------------------------------

First and foremost, I would like to thank BMC for giving me this problem. I've
done a lot of coding interviews, and most of them just slap down some algorithm
implementation question that you're supposed to summon from rote memory or throw
a bunch of boring HackerRank questions at you that force you to stay in very
tight margins of complexity, time, memory allocation, function arguments, etc.

I think above all, I admire freedom and creativity, and even if I don't get the
job, I really respect the fact that BMC let me have free reign on this and
choose how I wanted to implement it. Guiding without commanding is a nuanced
skill, and I think you've done a great job of it. Thank you for taking the time
to read all this and run everything :).

Anyway, enough waffling. Here's how you start the code:


**To run this code, all you need to do is install it with pip, then run
`$ hyperant-gui` or `$ hyperant-gui` in terminal from *anywhere* on your local
machine (symlinked to the same UI version in /bin, so they're both global).
That's it! After that, the GUI should open and you can enter your preferences.**

As a quick note, the vertex indices start at *ZERO*, so for example we would run
the original problem (ant traveling to opposite corner on a regular 3D cube ) by
entering 3, 0, 7 for the dimensions, entry index, and exit index. The number of
threads is up to you - I usually pick 50, but you can try out different settings.
Then just hit "Compute", or hit "End" if you want to try new settings.


For dimensions > 6, it takes a while to actually draw the hypercube, so for a
while you might just see a bunch of dots (all black except for one blue and one
red) appearing on the screen - this is expected.

------------------------------>>> Overview <<-----------------------------------

I debated which language to use for a while, and ultimately decided that C++ had
the most robust and generalizeable parallelization and linear algebra tools, so
I went with C++ for the backend (I'll describe what that actually entails in a
second). I actually had a frontend piece in Java for a bit, but it was very slow
and often quite ugly, even after I tried to parallelize it. Also, I found out
that turtles (or each language's equivalent) do a ton of extra stuff under the
hood, making operations involving them very poorly parallelizable. I really did
not like the available libraries in Java and C/C++, even though they were
purportedly meant to mirror Python's turtle library. I therefore decided to use
Python for the front end and make a simple GUI that will display the result and
render the surface (more on that in a sec).

In a sentence:
	I built a tool which generates a k-dimensional hypercube (where k is set by
	the user through the UI, or as a cmdline argument) and computes the mean
	time it would take to travel between any two vertices, where the specific
	vertices are again chosen by the user.

This is a twist on what was originally asked for (ant on 3-dimensional cube). I
am not 100% sure why I chose this, but here we are. It seemed more fun and
challenging. To clarify, my code works for k-dimensional hypercubes up to like
k == 20 ish. The exact number depends on your computer and how much space in
memory terminal is allowed to reserve during runtime. After k = 8 though it
takes a really long time to actually plot the hypercubes, no matter what
language I tried. I think this is just part-and-parcel with the complexity of
hypercubes themselves. I will now describe the back and front end components so
you can get a better idea of what I did.

------------------------------->>> Backend <<-----------------------------------

This part was coded entirely in C++. It is extremely parallelized and supports
multithreading. I tried to use as few external libaries as possible (besides
libomp of course, otherwise it wouldn't support multithreading) but one library
which my code does rely on for the very last step is Armadillo, which is
basically C++'s version of numpy for Python (sorry, numpy is still the GOAT). As
soon as I saw this problem, I thought "I can model this as a random walk on an
undirected k-regular graph", and then from there it was obvious that I could
represent the walk as an Absorbing Markov Chain where the destination vertex
was constructed to be the only absorbing node each time. I then implemented all
this as follows:

	* First, I had to solve the labeling problem. Given an abitrary dimension
	  number k, how do you label each node of your k-dimensional hypercube such
	  that each label is unique and all the edges make sense? I was originally
	  planning on doing something with the Hamming distance between the
	  bit-strings of each node in the range [0,2^k) but I figured out this
	  super quick and dirty way to get all the neighbors of a node. For node n,
	  all one needs to compute is [n^(1<<s)], s in [0,k). This is important,
	  because now I didn't have to make a whole separate function with a doubly
	  nested for loop. Instead, I could compute this atomically and parallelize
	  the original nested for loop directly! (which is what I did)

	* I then allocated and seeded my arrays with some small random values, as
	  is done with most modern linalg libraries. This left me with a matrix Q
	  (among others) that represented the submatrix of my original transition
	  matrix wihtout the rows and columns with the same index as the destination
	  node. A bit of randomness is important because if we happen to produce
	  extremely ill-conditioned matrices for a particular choice of settings,
	  there is at least some hope that within a few iterations, the matrix will
	  be a bit less ugly and a solver can at least approximate the value.

	* Next, I needed to get the inverse of Q so I could solve for the steady
	  state vector of my Markov Chain. Obviously, it's very inefficient to just
	  invert Q directly, so I decided to take the more efficient route and compute
	  the LU decomposition of A first. I did this **by hand**, using the famous
	  Dolittle Algorithm. This allowed me to parallelize a few inner loops more
	  directly. I checked the validity of these neighbors by seeing if LU indeed
	  equaled Q (pretty close in most cases).

	* The last part of my code needed to invert L and Q (which I sanitized
	  and fored to be perfectly triangular) relies on Armadillo's solve function
	  which, supposedly, is implemented as a thin wrapper over the core DGETRF,
	  SGETR, and LAPACK libraries, which honestly probably makes it way faster
	  than even a well parallelized solution.

	* With all this, I had enough to compute the steady state vector, from which
	 I then took the average of the i^th entry and the entries of its k-nearest
	 neighbors to try and avoiding any rounding errors or outliers that may have
	 appeared. I also checked my result for the steady state vector against the
	 "true" one generated by a naive Python helper script I wrote (not included)

The script `scripts/run-hypercube.sh` acts as a makefile for my backend code. Do
note that you will need clang installed to compile the code. Furthermore, to
run it on a Linux machine (or at least Linux OS in a VM or container), you will
need brew in order to install armadillo and openmp (installs are already
embedded in the build script).

I should also note I did a thorough analysis of my code with valgrind and there
are no memory leaks or bugs.

------------------------------>>> Frontend <<-----------------------------------

As I mentioned, I used Python for the frontend simple UI because it was fast and
effective. I didn't realize the sorry state that C++ and Java are in when it
comes to simple graphics rendering! Arcane functions and very counter-intuitive.
Basically, I made a simple UI that framed the problem as a Multi-dimensional
Bus or Train or some form of metropolitan ant transport that could take you from
point A to point B in the k-dimensional hypercube and returned the expected time
to do so. I then passed args to my C++ through Python internals and collected
the output so I could display it. My UI now takes the user's preferences and
actually draws the k-dimensional hypercube (after k=8 it gets very slow).
I should also note that I took care not to import any non-standard modules
(everything I imported is part of Python's out-of-the-box Lib).

As before, file "scripts/run-hypercube.sh" acts as a makefile for my backend C++
code, which is freshly compiled each time I run it. To run the python file, do
`$ hyperant-ui` or `$ hyperant-gui` and wait for the GUI to pop up. Note that
sometimes it takes a minute or two for the code to start drawing as you travel
to higher dimensions. To start the animation, just fill in the fields in the 
UI popup (ndim >=1, nthreads >=1, etc) and hit "Compute". If you want to change
certain parameter settings, hit "End" and repeat these steps to launch again.
Technically if you spam "Compute" with new parameters it will actually work,
but it's a bit wonky, so the recommended behavior is just to terminate the
script with the "End" button (avoid Control-C) and start again.


----------------------------->>> Conlusion <<-----------------------------------

That's it! Thanks for reading and for considering me for the position - I really
hope I get the privilige of working with you all :)


Twine upload instructions:

// Increment the version number in setup.py

// (inside hyperant)
`$ python3.7 setup.py sdist bdist_wheel`
`$ twine check dist/*` // check PASSED

// (inside \~/.pypirc )
// [pypi]
// username = __token__
// password = API_KEY from [here][https://pypi.org/manage/account/token/]

`$ twine upload --repository pypi --skip-existing dist/*  --verbose`


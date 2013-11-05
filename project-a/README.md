# Computational Physics
=======================

This repo contains all the code I porduced during the Computational Physics course at Imperial College London. Each project has its own dedicated directory. The data directory under each of these cannot be relied upon to stay constant as I'm constantly re-running the code with different parameters.

### Before You Start

To accompany this course I built a handy data viewer web-app which will monitor a file on you hard disk (typicallly the data outputted from the programs I've developed here) and periodically refresh a view of it in a web browser. It has some bespoke features such as a pendulum visualiser built in also. You can take a look at it [here](http://github.com/jeshuamaxey/data-viewer).

Instructions to get it set up can be found in the README of its repo, but feel free to direct any questions about it to me on [twitter](http://twitter.com/jeshuamaxey).

## Project-A

Project A involves evaluating different finite difference methods for the porblems of a single pendulum and double pendulum system. The Euler method, Leapfrog method (a special case of 2nd order Runga-Kutta method) and 4th order Runga-Kutta method are all used for modelling the single pendulum. Only 4th order Runga-Kutta is used for the double pedulum.


### Aside

I will endeavour to include the instructional PDF files I used while conducting the course in the appropriate directories in a future commit. They help in understanding the problem being solved, but you will require a firm understanding in numerical methods, differential calculus and finite difference approximations to solve the problems yourself

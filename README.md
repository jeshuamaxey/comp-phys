# Computational Physics

This repo contains all the code I produced during the Computational Physics course at Imperial College London. Each project has its own dedicated directory.

If you run the code yourself, be sure to add a `data` and `data/dp` directory in the `project-a` directory before trying to generate the output files. The aliases file include bash aliases that I've found I keep repeating during development. They're intended as a copy-paste utility to speed things up a little.

### Before You Start

To accompany this course I built a handy data viewer web-app which will monitor a file on you hard disk (typicallly the data outputted from the programs I've developed here) and periodically refresh a view of it in a web browser. It has some bespoke features such as a pendulum visualiser built in also. You can take a look at it [here](http://github.com/jeshuamaxey/data-viewer).

Instructions to get it set up can be found in the README of its repo, but feel free to direct any questions about it to me on [twitter](http://twitter.com/jeshuamaxey).

## Project-A

Project A involves evaluating different finite difference methods for the problems of a single pendulum and double pendulum system. The Euler method, Leapfrog method (a special case of 2nd order Runga-Kutta method) and 4th order Runga-Kutta method are all used for modelling the single pendulum. Only 4th order Runga-Kutta is used for the double pedulum.


### Aside

I will endeavour to include the instructional PDF files I used while conducting the course in the appropriate directories in a future commit. They help in understanding the problem being solved, but you will require a firm understanding in numerical methods, differential calculus and finite difference approximations to solve the problems yourself.
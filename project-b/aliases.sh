#!/bin/bash
alias go="echo compiling project-b.cpp && g++ project-b.cpp  -lgsl -o bin/project-b && echo compilation complete. Running project-b now... && ./bin/project-b && node update-file-list.js"

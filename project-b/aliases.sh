#!/bin/bash
alias go="echo compiling project-b.cpp && g++ project-b.cpp  -lgsl -o bin/project-b && echo compilation complete. Running project-b now... && ./bin/project-b && node update-file-list.js && say 'all done'"
echo "alias go all set :)"

# To run this script so that the aliases apply to this shell run:
# $ source ./aliases.sh


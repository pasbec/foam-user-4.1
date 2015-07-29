#!/bin/bash

###############################################################################

remote='dotfiles'
prefix='pascal-3.1/etc/bashrc.d'
branch='subtree/dot/bashrc.d'

if [ -e "$prefix" ]; then

    git subtree pull -P "$prefix" "$remote" "$branch" --squash

else

    echo "WARNING: Subtree for prefix '$prefix' soes not exist. Skipping."
    exit 1

fi

###############################################################################

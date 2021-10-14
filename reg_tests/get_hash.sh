#!/bin/bash

#set -x

commit_string=$(git log -1 --oneline)

commit_num=$(echo $commit_string | cut -c1-7)

echo ${commit_num}

export commit_num

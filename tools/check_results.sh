#!/bin/bash
LOGS=$(find $1 -name "*.out" -type f | grep ours)
# LOGS=$(find $1 -name "*.err" -type f | grep ours)
for f in $LOGS
do
    echo "$f:"
    errors=$(cat $f | grep -e "error" | grep -v -e "collisions_unsolved")
    warnings=$(cat $f | grep -e "warning")
	if [ -n "$errors" ]; then
        # echo -e "\033[0;31mFailing\033[0m: $(echo "$errors" | head -1)"
        echo -e "\033[0;31mErrors\033[0m: $errors"
    elif [ -n "$warnings" ]; then
        echo -e "\033[1;33mWarnings\033[0m: $warnings"
    else
        echo -e "\033[0;32mPassing\033[0m"
    fi
    echo
done

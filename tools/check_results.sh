#!/bin/bash
# Script to check the logs of our simulation to see if the simulation had
# errors, warnings, or passed.

LOGS=$(find $1 -name "*.out" -type f | grep ours)
# LOGS=$(find $1 -name "*.err" -type f | grep ours)

for f in $LOGS
do
    finished=$(cat "$f" | grep postprocess)
    if [ -z "$finished" ]; then
        continue
    fi
    echo "$f:"
    errors=$(cat $f | grep -e "error" | grep -v -e "collisions_unsolved")
    warnings=$(cat $f | grep -e "warning" | grep -v -e "collisions_unsolved")
	if [ -n "$errors" ]; then
        # echo -e "\033[0;31mFailing\033[0m: $(echo "$errors" | head -1)"
        first_error=$(echo "$errors" | head -n 1)
        num_errors=$(echo "$errors" | wc -l)
        printf "\033[0;31mErrors\033[0m: %s\n+%d more\n" "$first_error" "$num_errors"
    elif [ -n "$warnings" ]; then
        first_warning=$(echo "$warnings" | head -n 1)
        num_warnings=$(echo "$warnings" | wc -l)
        printf "\033[1;33mWarnings\033[0m: %s\n+%d more\n" "$first_warning" "$num_warnings"
    else
        printf "\033[0;32mPassing\033[0m\n"
    fi
    echo
done

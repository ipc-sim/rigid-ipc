#!/bin/bash

for i in {1..5}
do bash +x tools/run_eb_test.sh . tools/generate_falling_diamond_fixture.py "--e_b 1e-8 --c 1e-$i" data/results.local/ "eb_1e-8_c_1e-$i"
done

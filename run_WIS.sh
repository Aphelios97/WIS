#!/bin/bash

g++ -std=c++11 -O2 WIS.cpp -o WIS

(
echo $1
echo $2
echo $3
) | ./WIS

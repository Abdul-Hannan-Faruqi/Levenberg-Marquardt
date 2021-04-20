#! /bin/sh

a=1
while [ "$a" -le 100000 ]    # this is loop1
do
#   b=$((2 * $a))
#   while [ "$b" -le 20 ]  # this is loop2
#   do
      ./LM 1 $a 3 19
#      b=$(( $b + 1))
#   done
   a=$((10 * $a))
done

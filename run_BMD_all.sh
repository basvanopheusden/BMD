#!/bin/bash
for i in {1..5}

   do
    ./bmd x$i.txt  integral_table.txt params$i.txt changepoints$i.txt >output$i.txt
   done

exit
#!/bin/bash

for ex in $@
do
    echo ${ex}
    ldd ${ex} | grep "=> /" | awk '{print $3}'
done | sort -u | tar -ch -T - -f gathered.tar

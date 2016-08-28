#!/bin/csh

if ($#argv != 2) then
    echo "Usage: $0 <New sample_sheet, Old sample_sheet>"
        exit 0
endif
set nrows = `cat $1 | wc -l`
set nrows_keep = `expr $nrows - 1`
tail -$nrows_keep $1 > ss.tmp
cat $2 ss.tmp > ss.tmp.2
mv ss.tmp.2 $2
rm -f ss.tmp

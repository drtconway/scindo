# Stabby - an implementation of efficient interval stabbing

The stabby class used here is an implementation of the data structure
described in the paper

> Interval Stabbing Problems in Small Integer Ranges
> - Jens M. Schmidt, 2010.

In this paper, the author presents a clever data structure that
supports efficient queries over a set of (closed) intervals in a
dense integer domain.  The last bit is important, because it means
we can use a O(1) lookup in an array to find things, and during
construction, we can pre-compute the starting points for all queries.

The downside, is that in the context of genomic annotations, we
want to query over genomic positions which are sparse, with respect
to annoations.  To address this, we introduce a pair of rank/select
data structures to allow efficient mapping between the sparse and
dense domains.

Consider the following set of intervals from the annotation of the
gene TRAPPC4 on Chromosome 11:

```
    |-----------------------------------------------------------| [119018763,119025454]
    |-----------------------------------------------------------| [119018766,119025454]
    |-----------------------------------------------|             [119018766,119024134]
    |-------------------------------------------|                 [119018766,119023672]
    |------------------------------------------|                  [119018766,119023665]
    |------------------------------------------|                  [119018766,119023642]
    |----------------------------------------|                    [119018766,119023421]
    ||                                                            [119018766,119018970]
       ||                                                         [119019143,119019317]
       ||                                                         [119019143,119019210]
       |                                                          [119019143,119019188]
                ||                                                [119020150,119020253]
                ||                                                [119020205,119020253]
                              |-|                                 [119021760,119021886]
                               ||                                 [119021874,119021886]
                                            |-------|             [119023321,119024134]
                                            |---|                 [119023321,119023672]
                                            |--|                  [119023321,119023665]
                                            |--|                  [119023321,119023642]
                                            ||                    [119023321,119023421]
                                                          |-----| [119024857,119025454]
   -+---------+---------+---------+---------+---------+---------+-
    + 119018763         + 119020993         + 119023223         + 119025454
              + 119019878         + 119022108         + 119024338
```

If we take the positions, and condense them into ranks, the set of
intervals becomes:

```
    |-------------------|     [119018763,119025454]/[0,20]
     |------------------|     [119018766,119025454]/[1,20]
     |----------------|       [119018766,119024134]/[1,18]
     |---------------|        [119018766,119023672]/[1,17]
     |--------------|         [119018766,119023665]/[1,16]
     |-------------|          [119018766,119023642]/[1,15]
     |------------|           [119018766,119023421]/[1,14]
     ||                       [119018766,119018970]/[1,2]
       |--|                   [119019143,119019317]/[3,6]
       |-|                    [119019143,119019210]/[3,5]
       ||                     [119019143,119019188]/[3,4]
           |-|                [119020150,119020253]/[7,9]
            ||                [119020205,119020253]/[8,9]
              |-|             [119021760,119021886]/[10,12]
               ||             [119021874,119021886]/[11,12]
                 |----|       [119023321,119024134]/[13,18]
                 |---|        [119023321,119023672]/[13,17]
                 |--|         [119023321,119023665]/[13,16]
                 |-|          [119023321,119023642]/[13,15]
                 ||           [119023321,119023421]/[13,14]
                       ||     [119024857,119025454]/[19,20]
   -+---------+---------+-----
    + 0       + 10      + 20
```

This is now in a suitable form for the the algorithm discussed in
the Schmidt paper. We can construct a sparse array using the
implentation of Okanohara and Sadakane's *sdarray* in the SDSL
library by Simon Gog et al. This contains the positions of the
coordinates of interest, and then we can use `rank`/`select` to map
between the two domains. The beauty of this approach is that both
operations are _O(1)_.

The problem is that it doesn't readily support all the queries we
want, since we might want to query for annotations at the genomic
position chr11:119020256, which lies between two exons (i.e. in an
intron). Since, it doesn't coincide with with any of the positions
used to build the dense domain, and lies between positions 9 and
10, we can't efficiently query it directly. One solution would be
to use the dense domain, and stab both positions 9 and 10, and take
the intersection of the results, but this doubles the amount of
work.

Instead, we use an augmented dense domain. When we construct the
list of positions from the intervals, we note any positions that
are not contiguous in the sparse (genomic) domain, and add an element
to the dense domain to stand for all positions between the two.  In
practice, this approximately doubles the size of the domain, since
most positions we want to mention have gaps between them.

Taking the example above, this yields the set of intervals:

```
     |---------------------------------------|     [119018763,119025454]/[1,41]
       |-------------------------------------|     [119018766,119025454]/[3,41]
       |---------------------------------|         [119018766,119024134]/[3,37]
       |-------------------------------|           [119018766,119023672]/[3,35]
       |-----------------------------|             [119018766,119023665]/[3,33]
       |---------------------------|               [119018766,119023642]/[3,31]
       |-------------------------|                 [119018766,119023421]/[3,29]
       |-|                                         [119018766,119018970]/[3,5]
           |-----|                                 [119019143,119019317]/[7,13]
           |---|                                   [119019143,119019210]/[7,11]
           |-|                                     [119019143,119019188]/[7,9]
                   |---|                           [119020150,119020253]/[15,19]
                     |-|                           [119020205,119020253]/[17,19]
                         |---|                     [119021760,119021886]/[21,25]
                           |-|                     [119021874,119021886]/[23,25]
                               |---------|         [119023321,119024134]/[27,37]
                               |-------|           [119023321,119023672]/[27,35]
                               |-----|             [119023321,119023665]/[27,33]
                               |---|               [119023321,119023642]/[27,31]
                               |-|                 [119023321,119023421]/[27,29]
                                           |-|     [119024857,119025454]/[39,41]
   -+---------+---------+---------+---------+------
    + 0       + 10      + 20      + 30      + 40
```

Now, our query can map to dense position 19, which correctly stabs
between the relevant exons.

We implement this by generating the list of positions of interest
from the intervals as before, and we construct the same sparse
array, but in addition, we construct a list of the positions of the
coordinates in the augmented dense domain. We can now construct a
dense array with `rank`/`select` support. We can take ranks from
the original sparse domain, and use `select` to determine the
position in the augmented domain, and similarly, when we get an
interval back from the the data structure, we can use `rank` on the
end-points to get their ranks in the sparse domain, and we can then
use `select` to recover the original genomic position. Again, we
can perform these operations in _O(1)_.

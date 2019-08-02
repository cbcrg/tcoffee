#############
Release Notes
#############
.. Warning:: This log of recent modifications is not as thorough and accurate as it should be.


-9.86 New data structure for the primary library that results in highly improved running times for mcoffee and significantly decreased memory usage.


-5.80 Novel assembly algorithm (linked_pair_wise) and the primary library is now made of probcons style pairwise alignments (proba_pair)


-4.30 and upward: the FAQ has moved into a new tutorial document


-4.30 and upward: -in has will be deprecated and replaced by the flags: -profile,-method,-aln,-seq,-pdb


-4.02: -mode=dna is still available but not any more needed or supported. Use type=protein or dna if you need to force things


-3.28: corrected a bug that prevents short sequences from being correctly aligned


-Use of @ as a separator when specifying methods parameters


-The most notable modifications have to do with the structure of the input. From version 2.20, all files must be tagged to indicate their nature (A: alignment, S: Sequence, L: Library...). We are becoming stricter, but that's for your own good...


Another important modification has to do with the flag -matrix: it now controls the matrix being used for the computation

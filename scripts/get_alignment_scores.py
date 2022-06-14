import sys

for line in sys.stdin:
    A = line.rstrip().split()
    region = A[0].split(':') # contig:start-end
    contig = region[0]
    start, end = region[1].split('-')
    score = A[13].split(':')[2] # AS:i:score
    print('\t'.join([contig, start, end, score]))

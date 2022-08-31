import sys

# TODO add param for choosing score type
# either AS (alignment score) or ND (edit distance)
# TODO normalize scores
for line in sys.stdin:
    line = line.rstrip().split()
    region = line[0].split(':') # contig:start-end
    contig = region[0]
    start, end = region[1].split('-')

    # find the alignment score amongst the optional fields
    for field in line[11:]:
        field = field.split(':')
        if field[0] == 'AS':
            score = field[2]
            print('\t'.join([contig, start, end, score]))
            break

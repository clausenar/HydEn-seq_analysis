import sys

unhit=open(sys.argv[1])
pair=open(sys.argv[2])
out=sys.argv[3]

i=0
p=0

def read_id(header_line):
    # Compare only the base read ID, not the full header line: R1 and R2
    # headers carry different mate-number annotations (e.g. "... 1:N:0:..."
    # vs "... 2:N:0:...") even for the same underlying read pair.
    return header_line.split()[0]

with open(out, 'w') as the_file:

    for line in unhit:

        pair_line = pair.readline()
        while read_id(line) != read_id(pair_line):
            next(pair)
            next(pair)
            next(pair)
            pair_line = pair.readline()

        the_file.write(line)
        the_file.write(next(pair))
        the_file.write(next(pair))
        the_file.write(next(pair))

        next(unhit)
        next(unhit)
        next(unhit)


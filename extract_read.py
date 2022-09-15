import sys

unhit=open(sys.argv[1])
pair=open(sys.argv[2])
out=sys.argv[3]

i=0
p=0

with open(out, 'w') as the_file:

    for line in unhit:

        while line!=pair.readline():
            next(pair)
            next(pair)
            next(pair)
            
        the_file.write(line)
        the_file.write(next(pair))
        the_file.write(next(pair))
        the_file.write(next(pair))

        next(unhit)
        next(unhit)
        next(unhit)
   

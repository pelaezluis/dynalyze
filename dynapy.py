from dynalyze import Distance
from sys import argv

input_file = argv[1]
frames = int(argv[2])

a = Distance(input_file, frames)
a.graph()


from dynalyze import Dies, Distance, Energy
from sys import argv

input_file = argv[1]

distances = Distance(input_file)
energy = Energy(input_file) 
dies = Dies(distances.distances_dataframe(), energy.convert_structure(), input_file)
dies.save_dies()


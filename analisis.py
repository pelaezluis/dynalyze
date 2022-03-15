from dynalyze import Dies, Distance, Energy
from sys import argv

input_file = argv[1]

distances = Distance(input_file)
energy = Energy(input_file) 
dies = Dies(distances.distances_dataframe(), energy.convert_structure(), input_file)
#df = dies.join_energy_distances()
# print(dies.ordered_by_mean(df))
# print(dies.ordered_by_median(df))

#dies.energy_boxplot(dies.ordered_by_median(df), 'F0T')
#dies.save_dies()

# Colores
# "#2A7FFF"
# "#FF0000"
# "#00FF00"

dies.density_graph('F0T A', 'c3a', 'FAD603', 'red')


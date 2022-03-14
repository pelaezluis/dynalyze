# Importando librerías necesarias

import pandas as pd
import matplotlib.pyplot as plt
# import numpy as np
# import seaborn


class Distance:

    """
        Esta clase manejará los archivos de distancias
        Creará promedios cada cierto tiempo
        Y sus respectivas graficas
    """
    
    def __init__(self, input_file, frames):
        self.input_file = input_file
        self.output_file = 'avg_' + input_file.split('.')[0]
        self.frames = frames
        
    
    # Leer datos de vmd
    def load_data(self):
        """
            Load data from vmd txt file
        """

        # Verify if the file exists
        file = open(self.input_file, 'r')
        return [round(float(distance.split()[1]), 3) for distance in file]
    
    

    def mean(self):
        """
            Calculate the mean every n frames
        """
        data = self.load_data()
        counter = 0
        avg = []
        intervals = []
        frames = (len(data) // self.frames) + 2
        interval = self.frames / 100
        print(interval)
        for frame in range(frames):
            if counter == 0:
                avg.append(data[0])
                intervals.append(frame * interval)
                counter += 1
            elif counter > 0:
                data_slice = data[counter: counter + self.frames]
                counter += self.frames
                avg.append(round(sum(data_slice) / len(data_slice), 3))
                intervals.append(frame * interval)
        print(intervals)
        return (intervals, avg)


    def distances(self):
        """
            Create a dataframe with all the distance parameters
        """
        pass


    def graph(self):
        interval, avg = self.mean()
        plt.plot(interval, avg)
        plt.show()


class Energy:
    """
        Clase para manejo de las energias de interacción
        provenientes del multipdb
    """

    def __init__(self, input_file):
        self.input_file = input_file


    def load_data(self):
        """
            Load data from multipdb
        """
        data = pd.read_csv(self.input_file)
        #print(data.head())
        return data


    def extract_residues(self, dataframe, residue):
        """
            Extract the list of aminoacid residues in the multipdb
            and its data
        """
        res = dataframe[dataframe['fragment'] == residue]
        return res


    def residues(self):
        """
            Extract the list of ligand, residues and create a list of the 
            residues
        """
        dataframe = self.load_data()
        residue = list(dataframe[dataframe['frame'] == 0].fragment)
        residue = [res for res in residue]
        residue.sort()
        list_residue = []
        for res in residue:
            list_residue.append(self.extract_residues(dataframe, res))
        ligand = list(dataframe['ligand'])
        
        return ligand, residue, list_residue


    def convert_structure(self):
        """
            Convert the input file to a easy structure to read
        """
        ligand_, residues_, list_residue = self.residues()
        ligand = [lig.split('_')[0] for lig in ligand_]
        mode = [mode.split('_')[1] for mode in ligand_]
        residues = [res.split('_')[0] + res.split('_')[2] for res in residues_]
        rows = len(list(list_residue[0]['ie']))
        df = {'FNQ': ligand[:rows], 'Mode': mode[:rows]}
        for res in range(len(residues)):
            residue = list(map( lambda x: round(x, 3), list_residue[res]['ie']))
            df[residues[res]] = residue
        df = pd.DataFrame(df)
        return df
        

    def save(self):
        """
            Save the new dataframe in a csv file
        """
        e = self.input_file.split("_")
        output = f'{e[0]}_{e[1]}_ie.csv'
        dataframe = self.convert_structure()
        # verify if file already exists
        dataframe.to_csv(output, index=False)


class Dies:
    """
        Analyze energy and distance from Energy and Distance class
    """

    def __init__(self):
        pass


    def join_energy_distances(self):
        pass


    def energy_graph(self):
        pass


    def energy_boxplot(self):
        pass


    def energy_distance_graph(self):
        pass




#################################################
#                     Pruebas                   #
#################################################

# file = 'test.txt'
# test = Distance(file, 100)
# test.graph()

file = 'F0T_A_multi_ie_matrix.csv'
test = Energy(file)
print(test.convert_structure())
#test.save()

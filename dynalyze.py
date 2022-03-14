# Importando librerías necesarias
from os.path import isfile
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
    
    def __init__(self, input_file):
        self.input_file = input_file
        

    def load_input_data(self):
        """
            read input file with data and configuration
        """
        files = []
        if isfile(self.input_file):
            print('\n>>> Input file found...')
            print('>>> Reading parameters...\n')
            data = open(self.input_file, 'r')

            for parameter in data:
                line = parameter.split()
                if len(line) == 0:
                    continue
                elif line[0] == 'file:':
                    if isfile(line[1]):
                        print(f'* File {line[1]} found!')
                        files.append(line[1])
                    else:
                        print(f'File {line[1]} does not exist!')
                        exit()
                elif line[0] == 'interval:':
                    interval = int(line[1])
            return (files, interval)

        else:
            print('File input.dyn does not exist...')
            exit()


    def load_data(self, vmd_file):
        """
            Load data from vmd txt file and create a dataframe
        """
        file = open(vmd_file, 'r')
        return [round(float(distance.split()[1]), 3) for distance in file]
    

    def distances_dataframe(self):
        df = {}
        files = self.load_input_data()[0]
        print('\n>>> Creating distances dataframe...')
        for file in files:
            column = file.split('.')[0].split('_')[1]
            distances = self.load_data(file)
            df[column] = distances

        df = pd.DataFrame(df)
        return df


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


    def load_input_data(self):
        if isfile(self.input_file):
            data = open(self.input_file, 'r')
            for line in data:
                csv_file = line.split()
                if len(csv_file) == 0:
                    continue
                elif csv_file[0] == 'ie_file:':
                    print(f'* File {csv_file[1]} found!\n')
                    return csv_file[1]

        else:
            print('File input.dyn does not exist...')
            exit()


    def load_data(self):
        """
            Load data from multipdb
        """
        input_file = self.load_input_data()
        data = pd.read_csv(input_file)
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

    def __init__(self, ie_df, d_df, input_file):
        self.ie_df = ie_df
        self.d_df = d_df
        self.input_file = input_file


    def load_input_data(self):
        if isfile(self.input_file):
            data = open(self.input_file, 'r')
            for line in data:
                output_file = line.split()
                if len(output_file) == 0:
                    continue
                elif output_file[0] == 'output:':
                    print(f'* File {output_file[1]} found!\n')
                    return output_file[1]
            
        else:
            print('File input.dyn does not exist...')
            exit()


    def join_energy_distances(self):
        print('>>> Joining energy and distance dataframes...\n')
        return pd.concat([self.d_df, self.ie_df], axis=1)

    
    def save_dies(self):
        df = self.join_energy_distances()
        output = self.load_input_data()
        if isfile(output):
            option = input(f'*** {output} already exists. Overwrite it [y/n] ***')
            if option.lower() == 's':
                df.to_csv(output, index=False)
            elif option.lower() == 'n':
                pass
            else:
                print('*** Invalid option. Closing the script...')
        else:
            df.to_csv(output, index=False)


    def energy_graph(self):
        pass


    def energy_boxplot(self):
        pass


    def energy_distance_graph(self):
        pass

    



#################################################
#                     Pruebas                   #
#################################################

# file = 'input.dyn'
# test = Distance(file)
# test.distances_dataframe()

# test = Energy(file)
# test.load_input_data()
# print(test.convert_structure())
# test.save()

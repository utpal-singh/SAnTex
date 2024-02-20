# sage/ebsd/ebsd.py
from .ctf_parser import Ctf
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
class EBSD:

    def __init__(self, filename) -> None:
        self._filename = filename
        self.ctf = Ctf(self._filename)
    
    def get_euler_angles(self, phase):
        """
        phase should be integer number, not a phase name. This is mapped in another function
        """
        data, header_data = self.ctf.get_data()
        phase_df = data[data['Phase'] == phase]
        euler_angles = phase_df[['Euler1', 'Euler2', 'Euler3']]
        euler_angles = euler_angles.reset_index(drop = True)
        # print(euler_angles.iloc[0]["Euler1"])
        # print(euler_angles.loc[3450]["Euler1"])
        return euler_angles
    
    def plot(self):
        data, header_data = self.ctf.get_data()
        # df_new = data[data['Phase'] == 1]
        df_new = data[data['Phase'].isin([0, 1, 2])]
        plt.scatter(df_new['X'], df_new['Y'], c=df_new['Phase'], cmap='viridis')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('EBSD Map with Color Based on Phase')
        # plt.colorbar(label='Phase')
        plt.show()

    def phases_names(self):
        phases_info =  self.ctf.phases_info()[6]
        column_keep = ["phase"]
        phases_ = phases_info[column_keep].copy()
        return phases_
    
    def phases_info(self):
        phases_info_ =  self.ctf.phases_info()[6]
        return phases_info_
    
    def get_euler_angles_by_phase(self, phase):
        """
        This will be written further so that the euler angles can be extracted using phase names
        """
        pass

    def get_phases_numbers(self):
        return self.ctf.get_phases_numbers()
    
    def phases(self):
        df_phases = self.ctf.get_data()[0]
        phases_list = list(self.phases_names()["phase"])
        phases_list.insert(0, 'NaN')
        # Convert list to dictionary with keys starting from 0
        phase_names = {i: item for i, item in enumerate(phases_list)}

        # Map phase numbers to their names
        df_phases['Phase'] = df_phases['Phase'].map(phase_names)

        # Group by Phase and calculate the count
        phase_counts = df_phases['Phase'].value_counts().sort_index()

        # Calculate the total count
        total_count = phase_counts.sum()

        # Calculate the percentage of each phase
        phase_percentages = (phase_counts / total_count) * 100

        # Sort by percentage in descending order
        phase_counts_sorted = phase_counts.sort_values(ascending=False)
        phase_percentages_sorted = phase_percentages.sort_values(ascending=False)

        result_df = pd.concat([phase_counts_sorted, phase_percentages_sorted], axis=1)
        result_df.columns = ['Count', 'Percentage']

        return result_df
    
    def get_ebsd_data(self):
        data, header_data = self.ctf.get_data()
        return data, header_data


if __name__ == "__main__":
    ctfobj = EBSD("ebsd.ctf")
    print(ctfobj.get_euler_angles(3))
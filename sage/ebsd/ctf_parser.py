import pandas as pd
import csv

class Ctf:
    def __init__(self, filename) -> None:
        self._filename = filename

    def phases_info(self):

        lookup_string = ["Euler angles refer to Sample Coordinate system (CS0)!", "Phases"]

        with open(self._filename, 'r') as f:
            reader = csv.reader(f, delimiter='\t')

            for i, line in enumerate(reader):
                line_str = '\t'.join(line)

                if lookup_string[0] in line_str:

                    fields = line_str.split('\t')

                    self.mag = fields[fields.index('Mag') + 1]
                    self.coverage = fields[fields.index('Coverage') + 1]
                    self.device = fields[fields.index('Device') + 1]
                    self.kv = fields[fields.index('KV') + 1]
                    self.tilt_angle = fields[fields.index('TiltAngle') + 1]

                    line_number = i + 1


                if lookup_string[1] in line_str:
                    fields = line_str.split('\t')
                    self.nphases = int(fields[fields.index('Phases') + 1])
                    phases_start_line = i + 1


        phases_columns = ["cell_length_a", "cell_length_b", "cell_length_c", "cell_angle_alpha", "cell_angle_beta", "cell_angle_gamma", "phase", "symmetry", "nan1", "nan2", "nan3", "author"]
        self.phases = pd.read_csv(self._filename, delimiter='[;\t]', skiprows=phases_start_line, nrows=self.nphases, names=phases_columns, index_col=None, engine="python")
        self.phases.index = self.phases.index + 1

        return self.nphases, self.mag, self.coverage, self.device, self.kv, self.tilt_angle, self.phases
            

    def get_data(self):
        with open(self._filename, 'r') as file:
            lines = file.readlines()

    # print(lines[:15])
        header_start = lines.index('JobMode\tGrid\n')
        # print(header_start)
        data_start = lines.index('Phase\tX\tY\tBands\tError\tEuler1\tEuler2\tEuler3\tMAD\tBC\tBS\n')
        # print(data_start)

        header_data = pd.read_csv(self._filename, delimiter='\t', skiprows=header_start+1, nrows=7, names=["info", "value"])
        # print(header_data)

        data = pd.read_csv(self._filename, delimiter='\t', skiprows=data_start)
        # print(data)
        # print(data)
        # header_info = pd.read_csv(filename, sep="\t", nrows=data_start)
        # print(data)
        return data, header_data

    def get_phases_numbers(self):
        return self.get_data

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    ctf = Ctf("ebsd.ctf")
    # ctf.get_data()
    ctf.phases_info()
    data, header_data = ctf.get_data()
    # df_new = data[data['Phase'] == 1]
    df_new = data[data['Phase'].isin([0, 1, 2])]
    # define phase colors
    # colors = [plt.cm.tab10(0), plt.cm.tab10(1), plt.cm.tab10(2)]
    # # plot phase map
    # plt.scatter(df_new.X, df_new.Y, c=df_new.Phase, cmap=colors, s=1)
    # plt.title('Phase Map')
    # plt.xlabel('X')
    # plt.ylabel('Y')
    # plt.show()


    import matplotlib.pyplot as plt

    # Assuming df_new is your DataFrame
    plt.scatter(df_new['X'], df_new['Y'], c=df_new['Phase'], cmap='viridis')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('2D Map with Color Based on Phase')
    # plt.colorbar(label='Phase')
    plt.savefig('2d_map_2.png')
    plt.show()

    # import matplotlib.pyplot as plt

    # # Assuming df_new is your DataFrame
    # df_new["Phase"] = df_new["Phase"].map({0: "phase0", 1: "phase1", 2: "phase2"})

    # # Create a scatter plot for each phase
    # for phase in df_new["Phase"].unique():
    #     phase_data = df_new[df_new['Phase'] == phase]
    #     plt.scatter(phase_data['X'], phase_data['Y'], label=phase)

    # plt.xlabel('X')
    # plt.ylabel('Y')
    # plt.title('2D Map with Color Based on Phase')
    # plt.legend()
    # plt.savefig('2d_map.png')
    # plt.show()




    # print(data.loc[0, "Phase"])
    # print(len(data))
    # for i in range(len(data)):
    #     if data.loc[i, "Phase"] == 1:
    #         print(data.iloc[i, "Euler1"])

# print(data)


# class Ctf:
#     def __init__(self, filepath) -> None:
#         self.filepath = filepath
    
#     def get_database(self):
#         df = pd.read_csv
#         print(df)


# if __name__ == "__main__":
#     ctf = Ctf("myFile.ctf")
#     ctf.get_database()
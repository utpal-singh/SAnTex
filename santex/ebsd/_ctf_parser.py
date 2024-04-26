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

        header_start = lines.index('JobMode\tGrid\n')
        data_start = lines.index('Phase\tX\tY\tBands\tError\tEuler1\tEuler2\tEuler3\tMAD\tBC\tBS\n')
        header_data = pd.read_csv(self._filename, delimiter='\t', skiprows=header_start+1, nrows=7, names=["info", "value"])
        data = pd.read_csv(self._filename, delimiter='\t', skiprows=data_start)
        return data, header_data

    def get_phases_numbers(self):
        return self.get_data

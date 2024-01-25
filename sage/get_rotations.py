from ctf_parser import Ctf

class EulerAngles:
    def __init__(self, filename) -> None:
        self._filename = filename
    
    def get_euler_angles(self, phase):
        """
        phase should be integer number, not a phase name. This is mapped in another function
        """
        ctf = Ctf(self._filename)
        data, header_data = ctf.get_data()
        phase_df = data[data['Phase'] == phase]
        euler_angles = phase_df[['Euler1', 'Euler2', 'Euler3']]
        euler_angles = euler_angles.reset_index(drop = True)
        # print(euler_angles.iloc[0]["Euler1"])
        # print(euler_angles.loc[3450]["Euler1"])
        return euler_angles
    
    def get_euler_angles_by_phase(self, phase):
        """
        This will be written further so that the euler angles can be extracted using phase names
        """
        pass

    

if __name__ == "__main__":
    ctfobj = EulerAngles("ebsd.ctf")
    print(ctfobj.get_euler_angles(3))
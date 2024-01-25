from get_rotations import EulerAngles


ctfobj = EulerAngles("ebsd.ctf")
print(ctfobj.get_euler_angles(3))


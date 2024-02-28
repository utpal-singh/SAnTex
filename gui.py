# # sage/gui.py
# import tkinter as tk
# from sage.material.material import Material
# from sage.ebsd.ebsd import EBSD
# from sage.tensor.tensor import Tensor
# from sage.anisotropy.anisotropy import Anisotropy

# class SageGUI:
#     def __init__(self, master):
#         self.master = master
#         master.title("Sage Library GUI")

#         self.material_instance = Material(name="Example Material", properties={"density": 2.5})
#         self.ebsd_instance = EBSD(data="Sample EBSD Data")
#         self.tensor_instance = Tensor(values=[1, 2, 3, 4])
#         self.anisotropy_instance = Anisotropy(material=self.material_instance, tensor=self.tensor_instance)

#         # GUI components
#         self.label = tk.Label(master, text="Sage Library GUI")
#         self.label.pack()

#         self.text_display = tk.Text(master, height=10, width=40)
#         self.text_display.pack()

#         self.button_display_properties = tk.Button(master, text="Display Material Properties", command=self.display_material_properties)
#         self.button_display_properties.pack()

#         self.button_analyze_ebsd = tk.Button(master, text="Analyze EBSD Data", command=self.analyze_ebsd_data)
#         self.button_analyze_ebsd.pack()

#         self.button_calculate_magnitude = tk.Button(master, text="Calculate Tensor Magnitude", command=self.calculate_tensor_magnitude)
#         self.button_calculate_magnitude.pack()

#         self.button_calculate_effect = tk.Button(master, text="Calculate Anisotropy Effect", command=self.calculate_anisotropy_effect)
#         self.button_calculate_effect.pack()

#     def display_material_properties(self):
#         properties_text = f"Material: {self.material_instance.name}\nProperties: {self.material_instance.properties}"
#         self.text_display.insert(tk.END, properties_text + "\n\n")

#     def analyze_ebsd_data(self):
#         analysis_text = f"Analyzing EBSD Data: {self.ebsd_instance.data}"
#         self.text_display.insert(tk.END, analysis_text + "\n\n")

#     def calculate_tensor_magnitude(self):
#         magnitude = self.tensor_instance.calculate_magnitude()
#         magnitude_text = f"Tensor Magnitude: {magnitude}"
#         self.text_display.insert(tk.END, magnitude_text + "\n\n")

#     def calculate_anisotropy_effect(self):
#         effect_text = f"Calculating Anisotropy Effect for {self.material_instance.name}"
#         self.text_display.insert(tk.END, effect_text + "\n\n")

# if __name__ == "__main__":
#     root = tk.Tk()
#     app = SageGUI(root)
#     root.mainloop()

import tkinter as tk
from sage import Material

# class MaterialGUI:
    # def __init__(self):
    #     self.material_instance = Material()

    #     self.root = tk.Tk()
    #     self.root.title("Material Properties GUI")

    #     self.phase_label = tk.Label(self.root, text="Enter Material Phase:")
    #     self.phase_entry = tk.Entry(self.root)

    #     self.get_properties_button = tk.Button(self.root, text="Get Material Properties", command=self.get_material_properties)
    #     self.get_matrix_button = tk.Button(self.root, text="Get Voigt Matrix", command=self.get_voigt_matrix)

    #     self.result_label = tk.Label(self.root, text="Result will be displayed here.")

    #     self.phase_label.pack()
    #     self.phase_entry.pack()
    #     self.get_properties_button.pack()
    #     self.get_matrix_button.pack()
    #     self.result_label.pack()

    #     self.root.mainloop()

    # def get_material_properties(self):
    #     phase = self.phase_entry.get()
    #     properties = self.material_instance.get_properties_by_phase(phase)
    #     self.display_result(properties)

    # def get_voigt_matrix(self):
    #     phase = self.phase_entry.get()
    #     voigt_matrix = self.material_instance.get_voigt_matrix(phase)
    #     self.display_result(voigt_matrix)

    # def display_result(self, result):
    #     if result is not None:
    #         self.result_label.config(text=str(result))
    #     else:
    #         self.result_label.config(text="Phase not found in the materials database.")

# Run the GUI

import json
import tabulate
import tkinter as tk
from tkinter import ttk

class MaterialGUI:
    def __init__(self, master):
        self.master = master
        self.master.title("Material Information GUI")

        self.material_instance = Material()

        self.create_widgets()

    def create_widgets(self):
        tab_control = ttk.Notebook(self.master)

        phases_tab = ttk.Frame(tab_control)
        primary_phases_tab = ttk.Frame(tab_control)
        crystal_systems_tab = ttk.Frame(tab_control)

        tab_control.add(phases_tab, text="Available Phases")
        tab_control.add(primary_phases_tab, text="Available Primary Phases")
        tab_control.add(crystal_systems_tab, text="Available Crystal Systems")

        self.create_table(phases_tab, self.material_instance.availablePhases())
        self.create_table(primary_phases_tab, self.material_instance.availablePrimaryPhases())
        self.create_table(crystal_systems_tab, self.material_instance.availableCrystalSystems())

        tab_control.pack(expand=1, fill="both")

    def create_table(self, tab, data):
        table_label = ttk.Label(tab, text=data, font=("Arial", 12), anchor="center", padding=(10, 10))
        table_label.pack(expand=1, fill="both")

def main():
    root = tk.Tk()
    app = MaterialGUI(root)
    root.geometry("800x600")
    root.mainloop()

if __name__ == "__main__":
    main()

MaterialGUI()

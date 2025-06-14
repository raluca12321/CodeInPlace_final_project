import customtkinter as ctk
from rdkit import RDLogger, Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, AllChem, DataStructs
import seaborn as sns
import matplotlib.pyplot as plt

#Disable RDKit warnings
RDLogger.DisableLog('rdApp.*') #Useful when the user inputs an invalid string in order to not have errors from rdkit and too much text in the terminal

#Part 1 - get the SMILES from the user and convert them to mol object
class MoleculesAnalyzer(ctk.CTk):
    def __init__(self):
        super().__init__()

        ctk.set_appearance_mode("Light")
        ctk.set_default_color_theme("green")
        self.geometry("600x500")
        self.title("Molecules Analyzer")

        self.molecules_list = [] #Create an empty list to store the inputs from the user

        self.init_widgets()

    def init_widgets(self):
        self.label = ctk.CTkLabel(self, text="Enter at least 2 SMILES strings one at a time.") 
        self.label.pack(pady=10)

        self.entry = ctk.CTkEntry(self, width=300, placeholder_text="e.g. CCO")
        self.entry.pack(pady=10)

        self.add_button = ctk.CTkButton(self, text="Add SMILES", command=self.add_smiles)
        self.add_button.pack(pady=10)

        self.smiles_list_label = ctk.CTkLabel(self, text="", justify="left", wraplength=500)
        self.smiles_list_label.pack(pady=10)

        self.status_label = ctk.CTkLabel(self, text="", text_color="gray")
        self.status_label.pack(pady=5)

        self.finish_button = ctk.CTkButton(self, text="Finish", command=self.finish_input) #Finish Button to appear after 2 valid SMILES strings
        self.finish_button.pack(pady=10)

    def finish_input(self):
        if len(self.molecules_list) < 2:
            self.status_label.configure(text="⚠️ Please enter at least 2 SMILES before finishing.")
        
        else:
            smiles_display = "\n".join([s for s, _ in self.molecules_list])
            self.smiles_list_label.configure(text=f"Final list of SMILES:\n{smiles_display}")
            self.status_label.configure(text=" Input finished. You can proceed!")


    def add_smiles(self):
        smiles = self.entry.get().strip()
        mol = Chem.MolFromSmiles(smiles)

        if mol: #Condition if there's a valid SMILES string 
            self.molecules_list.append((smiles, mol))
            self.entry.delete(0, 'end')
            self.status_label.configure(text="SMILES added")

            smiles_display = "\n".join([s for s, _ in self.molecules_list])
            self.smiles_list_label.configure(text=f"SMILES entered:\n{smiles_display}")

        else: #Condition if the SMILES string is not valid
            self.status_label.configure(text=" Invalid SMILES string. Try again.")


if __name__ == "__main__":
    app = MoleculesAnalyzer()
    app.mainloop()

import customtkinter as ctk
from rdkit import RDLogger, Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, AllChem, DataStructs
import seaborn as sns
import matplotlib.pyplot as plt
from PIL import Image, ImageTk
import io

#Disable RDKit warnings
RDLogger.DisableLog('rdApp.*') 

class MoleculesAnalyzer(ctk.CTk):
    def __init__(self):
        super().__init__()

        ctk.set_appearance_mode("Light")
        ctk.set_default_color_theme("blue")
        self.geometry("600x400")
        self.title("Molecules Analyzer")

        self.molecules_list = [] #Create an empty list to store the inputs from the user
        self.default_smiles = ["CCO", "c1ccccc1", "CC(=O)O"]

        self.show_input_screen() # it ensures the first screen the user sees is the one for inputting SMILES.

  #Define a function to clear the screen
    def clear_screen(self):
        for widget in self.winfo_children():
            widget.pack_forget()


    def show_input_screen(self):
        self.clear_screen()
        
        self.label = ctk.CTkLabel(self, text = "Enter SMILES (or use defaults):")  #Label to guide the user
        self.label.pack(pady=10)

        self.entry = ctk.CTkEntry(self, width=300, placeholder_text="e.g. CCO") #Entry box for use to input SMILES
        self.entry.pack(pady=10)

        self.add_button = ctk.CTkButton(self, text="Add SMILES", command=self.add_smiles) #Button to add SMILES
        self.add_button.pack(pady=5)

        self.default_button = ctk.CTkButton(self, text="Use Default SMILES", command=self.use_default_smiles) #Button to use default list of SMILES
        self.default_button.pack(pady=5)

        self.finish_button = ctk.CTkButton(self, text="Finish", command=self.finish_input) #Button when user finishes the input of SMILES
        self.finish_button.pack(pady=5)

        self.smiles_list_label = ctk.CTkLabel(self, text="", justify="left") #Displays the list of SMILES strings the user has entered (or the default list if used)
        self.smiles_list_label.pack(pady=5)

        self.status_label = ctk.CTkLabel(self, text="", text_color="gray") #Provides user feedback/messages, like whether a SMILES was added successfully or if it was invalid.
        self.status_label.pack(pady=5)


    def finish_input(self): #For finish button
        if len(self.molecules_list) < 2:
            self.status_label.configure(text="⚠️ Add at least 2 SMILES before finishing.")
        else:
            self.show_function_menu()  

    def show_function_menu(self):
        self.clear_screen()

        self.menu_label = ctk.CTkLabel(self, text = "Select what you want do to:")
        self.menu_label.pack(pady = 10)

        self.prop_button = ctk.CTkButton(self, text="1. Show Properties", command=self.show_properties)
        self.prop_button.pack(pady=5)

        self.structures_button = ctk.CTkButton(self, text="2. Show 2D Structures", command=self.show_structures)
        self.structures_button.pack(pady=5)

        self.similarity_button = ctk.CTkButton(self, text="3. Show Similarity Heatmap", command=self.show_similarity)
        self.similarity_button.pack(pady=5)

        self.back_button = ctk.CTkButton(self, text="Back to Input", command=self.show_input_screen)
        self.back_button.pack(pady=10)


    def show_similarity(self): #PART 3 - heatmap similarity
        self.clear_screen()

    #Compute fingerprints
        fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) for smi, mol in self.molecules_list]

    #Compute similarity matrix - Tanimoto similarity
        similarity_matrix = []         #2D list of similarity scores
        for fp1 in fingerprints:
            similarities = [DataStructs.TanimotoSimilarity(fp1, fp2) for fp2 in fingerprints]
            similarity_matrix.append(similarities)
        
        smiles_labels = [smi for smi, _ in self.molecules_list]

    #Create heatmap figure
        fig, ax = plt.subplots(figsize=(5, 4))
        sns.heatmap(similarity_matrix, annot=True, fmt=".2f", cmap="coolwarm", xticklabels=smiles_labels, yticklabels=smiles_labels, ax=ax)
        plt.title("Molecular Similarity Heatmap")
        plt.tight_layout()

    #Save figure to bytes buffer
        with io.BytesIO() as buffer:
            fig.savefig(buffer, format="PNG")
            buffer.seek(0)
            pil_image = Image.open(buffer).copy()

        plt.close(fig)  # Close the figure to free memory

    #Convert to PhotoImage for Tkinter
        tk_image = ImageTk.PhotoImage(pil_image)

    #Display in CTkLabel (with image)
        image_label = ctk.CTkLabel(self, image=tk_image, text="")
        image_label.image = tk_image  # keep reference!
        image_label.pack(pady=10)

    #Back button
        self.back_button = ctk.CTkButton(self, text="Back", command=self.show_function_menu)
        self.back_button.pack(pady=10)


    def show_structures(self): #PART 2 - show the 2D chemical structures of the molecules
        self.clear_screen()

    #Extract only mol objects for drawing
        mols = [mol for _, mol in self.molecules_list]
        legends = [smi for smi, _ in self.molecules_list]
        
    #Create a single grid image of all molecules
        img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(200,200), legends=[smi for smi, _ in self.molecules_list], useSVG=False)

    #Convert to Tkinter-compatible image
        with io.BytesIO() as buffer:
            img.save(buffer, format="PNG")
            buffer.seek(0)
            pil_image = Image.open(buffer)
            tk_image = ImageTk.PhotoImage(pil_image)
    
    # Display the image
        image_label = ctk.CTkLabel(self, image=tk_image, text="")
        image_label.image = tk_image  # Prevent garbage collection
        image_label.pack(pady=10)

        self.back_button = ctk.CTkButton(self, text="Back", command=self.show_function_menu)
        self.back_button.pack(pady=10)
    

    def show_properties(self):
        self.clear_screen()

        output = ""
        for smi, mol in self.molecules_list:
            MW = Descriptors.MolWt(mol) #molecular weight
            LogP = Descriptors.MolLogP(mol) #octanol - water partition coefficient
            HBA = Descriptors.NOCount(mol)  #hydrogen bond acceptors
            HBD = Descriptors.NHOHCount(mol) #hydrogen bond donors
            aromatic = rdMolDescriptors.CalcNumAromaticRings(mol) #number of aromatic rings


            output += f"SMILES: {smi}\n"
            output += f" - MW: {MW:.2f}, LogP: {LogP:.2f}, HBA: {HBA}, HBD: {HBD}, Aromatic Rings: {aromatic}\n\n"

        self.prop_label = ctk.CTkLabel(self, text=output, justify="left", wraplength=600)
        self.prop_label.pack(pady=10)

        self.back_button = ctk.CTkButton(self, text="Back", command=self.show_function_menu)
        self.back_button.pack(pady=10)



    def use_default_smiles(self): #function to make the button "use dafault smiles list" work
        self.molecules_list = [(smi, Chem.MolFromSmiles(smi)) for smi in self.default_smiles]
        smiles_display = "\n".join(self.default_smiles)
        self.smiles_list_label.configure(text=f"Using default SMILES:\n{smiles_display}")
        self.status_label.configure(text="Default SMILES loaded.")    


    def add_smiles(self): #function to make the button "add SMILES" work
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
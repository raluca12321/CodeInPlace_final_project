#Importing Libraries
from rdkit import RDLogger, Chem # first one is to disable RDKit warnings
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, AllChem, DataStructs
import seaborn as sns
import matplotlib.pyplot as plt
import io
from PIL import Image

#Disable RDKit warnings
RDLogger.DisableLog('rdApp.*') #Useful when the user inputs an invalid string in order to not have errors from rdkit and too much text in the terminal

#Explanation of what this tool does
print("This tool analyzes molecules based on their SMILES notation (Simplified Molecular Input Line Entry System). ")
print("It's a notation which translates a molecule's 3D structure into text strings.")
print("For the molecules you want to analyze, search on the internet their SMILES)
print("")

molecules_list = [] #Create an empty list to store the inputs from the user
print("Enter at least 2 SMILES strings one at a time. Type 'done' to finish.")

#Part 1 - get the SMILES from the user and convert them to mol object
while True:
    if len(molecules_list) >= 2:
        SMILES = input("Enter a SMILES string (or 'done' to finish): ").strip()
        if SMILES.lower() == 'done': #.lower converts for example DOne to done
            break  # Exit the loop
    else:
        SMILES = input("Enter a SMILES string: ").strip()
    
    mol = Chem.MolFromSmiles(SMILES)
    while mol is None: #If the SMILES string is not valid - Keep prompting until a valid SMILES string is given
        SMILES = input("Please enter a valid SMILES string: ")
        mol = Chem.MolFromSmiles(SMILES)

    molecules_list.append((SMILES, mol))    #Store both the SMILES string and the mol object
print("")



# Part 2 - Calculate and display descriptors for each molecule
print("1. You will get some properties of these molecules: ")
for smi, mol in molecules_list:
    MW = Descriptors.MolWt(mol) #molecular weight
    LogP = Descriptors.MolLogP(mol) #octanol - water partition coefficient
    HBA = Descriptors.NOCount(mol)  #hydrogen bond acceptors
    HBD = Descriptors.NHOHCount(mol) #hydrogen bond donors
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol) #number of aromatic rings
    print("")
    print(f"SMILES: {smi} - properties: MW= {MW:.2f}, LogP= {LogP:.4f}, HBA= {HBA}, HBD= {HBD}, Aromatic Rings=  {num_aromatic_rings}") #Print the properties
    print("")

#Part 3 - Display image of each molecule if the user wants it
images_of_molecules = input("2. Do you want to get images the molecular structures of the previous SMILES strings? - yes/no: ")
if images_of_molecules.lower() == "yes":
    mols = [mol for _, mol in molecules_list]
    img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(200,200), legends=[smi for smi, _ in molecules_list])
    img.show()
print("")

#Part 4 - Molecules similarity
#Convert the molecule objects from the molecules_list to Morgan fingerprints
print("3. You will get the similarity matrix of these molecules and also visualize the matrix using a color-coded heatmap ")
fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) for smi, mol in molecules_list]

#Compute similarity matrix - Tanimoto similarity
similarity_matrix = []         #2D list of similarity scores
for fp1 in fingerprints:
    similarities = [DataStructs.TanimotoSimilarity(fp1, fp2) for fp2 in fingerprints]
    similarity_matrix.append(similarities)
    
#Adjust to 4 decimals and print the similarity matrix
for (smi, _), row in zip(molecules_list, similarity_matrix):
    formatted_row = [f"{val:.4f}" for val in row]
    print(f"{smi}: {formatted_row}")
print("")

# Create the heatmap visualization for similarity
smiles_labels = [smi for smi, mol in molecules_list]
sns.heatmap(similarity_matrix, annot=True, fmt=".2f", cmap="coolwarm",
            xticklabels=smiles_labels, yticklabels=smiles_labels)
plt.title("Molecular Similarity Heatmap")
plt.show()

#Import libraries
import sys
import logging
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from databaseManager import DatabaseManager
import base64
from PIL import Image

#Setting up the logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)
fh_log = logging.FileHandler('2743102_log.txt')
fh_log.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(fh_log)

#Function to convert the image into base64
def image_to_base64(image_path):
    try:
        with open(image_path, "rb") as img_file:
            # Read the image file
            img_data = img_file.read()
            # Encode the image data as base64
            base64_encoded = base64.b64encode(img_data).decode('utf-8')
            logger.info(f"Image {image_path} successfully converted.")
            return base64_encoded
    except Exception as e:
        logger.error(f"Error in converting image to Base64 : {e}")
        

#Initialise the parameter list for insert statement
params= []
#Initialise the database path
db_file = '/Users/harshitasrivastava/Desktop/Chemical Structures /Chemical Databases/protein.db'

try:
    #Create the database     
    db  = DatabaseManager(db_file)
    fd = open('2743102_Protein_Database.sql', 'r')
    db.createdb(fd.read())

    #Parse the .sdf file and create parameter list for Proteins table from it.
    suppl = Chem.SDMolSupplier('Molecules20.sdf')
except Exception  as e:
    logger.error(f"Error in reading sdf file : {e} ")
    raise SystemExit(1)
else:
    logger.info(f"Successfully created DB at {db_file}")


try:
    #Iterate through the file
    for mol in suppl:
        if mol is None: continue
        #Set a parameter list for each variable
        param = []

        #Column 1: molname
        name = mol.GetProp('_Name')
        param.append((name))

        #Column 2: structure .Get the 2-D structure of the coumpounds.
        img = AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol)
        img.save('structure.png')
        path = 'structure.png'
        base64_string = image_to_base64(path)
        param.append(base64_string)
        
            
        #Column 3: SMILES .Print the smiles coordinates
        smiles = Chem.MolToSmiles(mol)
        #print(f'SMILES String {smiles}')
        param.append(smiles)

        #Column 4: formula of the compund
        formula = mol.GetProp('Formula')
        param.append(formula)
        
        #Column 5: molecular wt of the compund
        mol_wt = Descriptors.MolWt(mol)
        param.append(round(mol_wt,2))

        #cOLUMN 6: LogP of the compound
        logP = Descriptors.MolLogP(mol)
        #print(f'LogP = {logP}')
        param.append(round(logP,2))

        #Column 7: LogD of the compound
        logD = mol.GetProp('LogD')
        param.append(float(logD))
        
        #Column 8: Number of atoms
        atoms = mol.GetNumAtoms()
        param.append(atoms)

        #Column 9: Number of rotatable bonds
        rotatable_bond = Descriptors.NumRotatableBonds(mol)
        param.append(rotatable_bond)

        #Column 10: number of rings
        ring_Count = Descriptors.RingCount(mol)
        param.append(ring_Count)

        #Column 11: number of H donors
        h_donors = Descriptors.NumHDonors(mol)
        param.append(h_donors)

        #Column 12: number of H acceptors
        h_acceptors = Descriptors.NumHAcceptors(mol)
        param.append(h_acceptors)
        
        #Column 13: number of valence electrons 
        valence_electrons = Descriptors.NumValenceElectrons(mol)
        param.append(valence_electrons)

        #Column 14: total  polar surface area (TPSA)
        tpsa = Descriptors.TPSA(mol)
        param.append(round(tpsa,2))

        #Column 15: Number of fused aromatic rings
        # Get the smallest set of symmetrical rings 
        sssr = Chem.GetSymmSSSR(mol)
        # Set Fused Aromatic Ring (FAR) count to 0 
        FAR_count = 0
        # Get the length of both rings  
        for x in range(len(sssr)):
                for j in range(x + 1, len(sssr)):
                    # Check whether all the bonds are aromatic 
                    if all(mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in sssr[x]) and all(mol.GetAtomWithIdx(atoms).GetIsAromatic() for atoms in sssr[j]):
                        # Check whether the rings share more than one atom
                        shared_atoms = len(set(sssr[x]) & set(sssr[j])) # both sets share common atom 
                        if shared_atoms > 1 : 
                            FAR_count += 1
        param.append(FAR_count)

        #Check if the molecular name is present
        if not param[0] == '':
            #Append to params list
            params.append(param)
        else:
            #If molecular name is not present, log the parameter list
            logger.info(f'Following parameters were not inserted into the database.\n{param}\n')
except Exception as e:
    logger.error(f'Error occured while during parameter value generation:{e}')
    raise SystemExit(1) 
else:
    logger.info(f"Parameters generated successfully.")

try:
    #Insert columns into database
    sql =  "INSERT INTO PROTEIN_DATABASE VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);"
    db  = DatabaseManager(db_file)
    db.insertdb(sql,params)
except Exception  as e:
    logger.error(f'Error occurred while inserting data into the database:\n{e}')
    raise SystemExit(1)
else:
    logger.info(f'Data inserted into {db_file}.\n')
#==============================================================================================================================#








  








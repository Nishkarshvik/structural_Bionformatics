
                                                    #################################################################################################
                                                    #                                                                                               #                                                                                                           #                 Python script for Prediction of ADME Properties using RDKIT                   #                                                                                                           #                                                                                               #                                                                                                           #                          DEVELOPED BY: Niskarsh Vikram Singh                                  #                                                                                                           #                                JRF at CARI, Kolkata                                           #                                                                                                           #                                                                                               #                                                                                                           #                    Prediction of Physiochemical Properties using RdKit                        #                                                                                                           #                                                                                               #                                                                                                           #                                                                                               #                                                                                                           #                                                                                               #                                                                                                           #                                                                                               #                                                                                                           #                                                                                               #                                                                                                           #                                                                                               #                                                                                                           #################################################################################################                                                                                                           
#!/bin/bash

# script usage : python3 script_name smiles.txt output.csv
# smiles.txt format expected : Phytochemical Name, Phytochemical ID, SMILES (Seperated by TAB as delimiter)
# Dependencies : Rdkit, python3, Pandas

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen
import pandas as pd
import sys

def calculate_admet_properties(smiles):
    """
    Calculate ADMET-related physicochemical properties for a given SMILES string using RDKit.

    Args:
        smiles (str): SMILES representation of the compound.

    Returns:
        dict: A dictionary containing the calculated properties.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string provided: {smiles}")

    # Calculate properties
    properties = {
        "Phytochemical Name": "",  # Placeholder for name, will be added later
        "Phytochemical ID": "",    # Placeholder for ID, will be added later
        "SMILES": smiles,
        "Formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
        "MW": round(Descriptors.MolWt(mol), 2),
        "Heavy atoms": round(Descriptors.HeavyAtomCount(mol), 2),
        "Aromatic heavy atoms": sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic()),
        "Fraction Csp3": round(Descriptors.FractionCSP3(mol), 2),
        "Rotatable bonds": Lipinski.NumRotatableBonds(mol),
        "H-bond acceptors": Lipinski.NumHAcceptors(mol),
        "H-bond donors": Lipinski.NumHDonors(mol),
        "MR": round(Crippen.MolMR(mol), 2),  # Molar Refractivity
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "LogP": round(Crippen.MolLogP(mol), 2),
        "Lipinski Violations": sum([
            Descriptors.MolWt(mol) > 500,
            Lipinski.NumHAcceptors(mol) > 10,
            Lipinski.NumHDonors(mol) > 5,
            Crippen.MolLogP(mol) > 5
        ]),
        "Ghose violations": sum([
            Descriptors.MolWt(mol) < 160 or Descriptors.MolWt(mol) > 480,
            Crippen.MolLogP(mol) < -0.4 or Crippen.MolLogP(mol) > 5.6,
            Descriptors.HeavyAtomCount(mol) < 20 or Descriptors.HeavyAtomCount(mol) > 70
        ]),
        "Veber violations": sum([
            Lipinski.NumRotatableBonds(mol) > 10,
            Descriptors.TPSA(mol) > 140
        ]),
        "Egan violations": sum([
            Crippen.MolLogP(mol) < -1 or Crippen.MolLogP(mol) > 5.2,
            Descriptors.TPSA(mol) > 131.6
        ]),
        "Muegge violations": calculate_muegge_violations(mol),  # Muegge Violations
        "Bioavailability Score": round(calculate_bioavailability_score(mol), 2),
        "HIA": calculate_hia(mol),
        "BBB": calculate_bbb(mol)
    }

    return properties

def calculate_muegge_violations(mol):
    """
    Calculate Muegge rule violations based on the following criteria:
    1. MW > 500
    2. LogP > 5
    3. Number of rotatable bonds > 10
    4. Number of hydrogen bond donors > 5
    5. Number of hydrogen bond acceptors > 10
    6. TPSA > 140

    Args:
        mol: RDKit molecule object.

    Returns:
        int: The number of violations.
    """
    violations = 0
    if Descriptors.MolWt(mol) > 500:
        violations += 1
    if Crippen.MolLogP(mol) > 5:
        violations += 1
    if Lipinski.NumRotatableBonds(mol) > 10:
        violations += 1
    if Lipinski.NumHDonors(mol) > 5:
        violations += 1
    if Lipinski.NumHAcceptors(mol) > 10:
        violations += 1
    if Descriptors.TPSA(mol) > 140:
        violations += 1
    return violations

def calculate_bioavailability_score(mol):
    """
    Calculate the bioavailability score based on certain physicochemical properties.
    
    Bioavailability score is derived from Lipinski’s rule and other relevant descriptors.
    
    Args:
        mol: RDKit molecule object.
    
    Returns:
        float: Bioavailability score, ranging from 0 (low) to 1 (high).
    """
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    rot_bonds = Lipinski.NumRotatableBonds(mol)

    # Bioavailability score logic
    score = 0
    if mw < 500:
        score += 0.2  # MW < 500 is favorable
    if logp >= -1 and logp <= 5:
        score += 0.2  # LogP within range is favorable
    if tpsa <= 140:
        score += 0.2  # Low TPSA is favorable
    if rot_bonds <= 10:
        score += 0.2  # Fewer rotatable bonds is favorable
    return score

def calculate_hia(mol):
    """
    Estimate the human intestinal absorption (HIA) of a compound.
    Placeholder function to demonstrate the concept.
    
    Args:
        mol: RDKit molecule object.
    
    Returns:
        str: "High" or "Low" absorption based on simple rules.
    """
    logp = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)

    if logp >= -0.5 and logp <= 5 and tpsa <= 140:
        return "High"
    else:
        return "Low"

def calculate_bbb(mol):
    """
    Estimate the blood-brain barrier (BBB) permeability of a compound.
    Placeholder function to demonstrate the concept.
    
    Args:
        mol: RDKit molecule object.
    
    Returns:
        str: "Yes" or "No" indicating BBB permeability.
    """
    logp = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)

    if logp >= 1.5 and logp <= 3 and tpsa <= 60:
        return "Yes"
    else:
        return "No"

def main(input_file, output_file):
    print("Reading input file...")
    with open(input_file, "r") as file:
        lines = [line.strip().split() for line in file if line.strip()]

    print(f"Found {len(lines)} lines in the input file.")

    results = []
    print("Calculating properties for each compound...")
    for line in lines:
        if len(line) < 3:
            print(f"Skipping invalid line: {line}. Expected format: PhytochemicalName PhytochemicalID SMILES")
            continue

        phytochemical_name, phytochemical_id, smiles = line[0], line[1], line[2]
        try:
            properties = calculate_admet_properties(smiles)
            # Update with the Phytochemical Name and ID
            properties["Phytochemical Name"] = phytochemical_name
            properties["Phytochemical ID"] = phytochemical_id
            results.append(properties)
        except ValueError as e:
            print(f"Skipping invalid SMILES: {smiles}. Error: {e}")

    # Reorder the DataFrame to ensure Phytochemical Name, ID, and SMILES come first
    df = pd.DataFrame(results)
    cols = ["Phytochemical Name", "Phytochemical ID", "SMILES"] + [col for col in df.columns if col not in ["Phytochemical Name", "Phytochemical ID", "SMILES"]]
    df = df[cols]

    print("Saving results to CSV file...")
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 ADME_rdkit.py <input_file> <output_csv_file>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        main(input_file, output_file)

print ("Prediction of the physiochemical properties of all the compounds has been completed.")

exit


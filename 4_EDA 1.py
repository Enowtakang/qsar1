import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski


path = "C:/Users/HP/PycharmProjects/MachineLearningEnow/Comp_Drug_Discovery/labeled_curated_data.csv"
df = pd.read_csv(path)
# print(df.head())

"""
We'll compute molecular descriptors, then 
perform Exploratory Data Analysis (EDA)
"""

"""
Clean the canonical_smiles column 
"""
df_no_smiles = df.drop(columns='canonical_smiles')

smiles = []

for i in df.canonical_smiles.tolist():
    cpd = str(i).split('.')
    cpd_longest = max(cpd, key=len)
    smiles.append(cpd_longest)

smiles = pd.Series(smiles, name='canonical_smiles')

df_clean_smiles = pd.concat(
    [df_no_smiles, smiles], axis=1)

# print(df_clean_smiles.head())

"""
Calculate Lipinski descriptors

Christopher Lipinski, a scientist at Pfizer, 
came up with a set of rule-of-thumb for 
evaluating the drug-likeness of compounds. 

Such drug-likeness is based on the: 
Absorption, Distribution, 
Metabolism and Excretion (ADME) 
that is also known as the pharmacokinetic 
profile. 

Lipinski analyzed all orally active 
FDA-approved drugs in the formulation 
of what is to be known as the Rule-of-Five or 
Lipinski's Rule.

The Lipinski's Rule stated the following:

    Molecular weight < 500 Dalton
    Octanol-water partition coefficient (LogP) < 5
    Hydrogen bond donors < 5
    Hydrogen bond acceptors < 10
"""


# Inspired by: https://codeocean.com/explore/capsules?query=tag:data-curation

def lipinski(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)

        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])

        if (i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i + 1

    columnNames = ["MW", "LogP", "NumHDonors",
                   "NumHAcceptors"]
    descriptors = pd.DataFrame(data=baseData,
                               columns=columnNames)

    return descriptors


df_lipinski = lipinski(
    df_clean_smiles.canonical_smiles)

# print(df_lipinski.head())

"""
Combine both:
            - df 
            - df_lipinski
side-by-side
"""
df_combined = pd.concat(
    [df, df_lipinski], axis=1)

# print(df_combined)

"""
Convert IC50 to pIC50

To allow IC50 data to be more uniformly 
distributed, we will convert IC50 to the 
negative logarithmic scale which is essentially 
-log10(IC50).

This custom function pIC50() will accept a 
DataFrame as input and will:

    - Take the IC50 values from the standard_value column 
        and converts it from nM to M by multiplying 
        the value by 10.
   -  Take the molar value and apply -log10.
   -  Delete the standard_value column and create 
        a new pIC50 column.

# https://github.com/chaninlab/estrogen-receptor-alpha-qsar/blob/master/02_ER_alpha_RO5.ipynb

Point to note: 
Values greater than 100,000,000 
will be fixed at 100,000,000 otherwise the 
negative logarithmic value will become negative.
"""


def pIC50(input):
    pIC50 = []

    for i in input['standard_value_norm']:
        molar = i * (10 ** -9)  # Converts nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', 1)

    return x


def norm_value(input):
    norm = []

    for i in input['standard_value']:
        if i > 100000000:
            i = 100000000
        norm.append(i)

    input['standard_value_norm'] = norm
    x = input.drop('standard_value', 1)

    return x


# We will first apply the norm_value()
# function so that the values in the
# standard_value column is normalized.
df_norm = norm_value(df_combined)

# print(df_norm.head())
# print(df_norm.standard_value_norm.describe())


df_final = pIC50(df_norm)

# print(df_final.head())

# print(df_final.pIC50.describe())

"""
Let's write this to CSV file.
"""
# df_final.to_csv('pIC50_data.csv')

"""
Removing the 'intermediate' bioactivity class:
    Here, we will be removing the intermediate 
        class from our data set.
"""
df_2_class = df_final[
    df_final['class'] != 'intermediate']

# print(df_2_class)

"""

"""
df_2_class.to_csv('pIC50_no_intermediate.csv')

import pandas as pd

path = "C:/Users/HP/PycharmProjects/MachineLearningEnow/Comp_Drug_Discovery/bioactivity_data_1.csv"
data = pd.read_csv(path)
# print(data.head())

"""
Handling missing data:
    If any compounds has 'missing value' for the 
    standard_value and canonical_smiles column 
    then drop it.
"""
df_2 = data[data.standard_value.notna()]
df_2 = df_2[data.canonical_smiles.notna()]
# print(df_2.head())

"""
Handling Duplicates
"""
df_2_nr = df_2.drop_duplicates(['canonical_smiles'])
# print(df_2_nr.head())

"""
Combine the 3 columns: 
        - molecule_chembl_id, 
        - canonical_smiles,
        - standard_value 
 into a DataFrame
"""
selection = ['molecule_chembl_id',
             'canonical_smiles',
             'standard_value']
df_3 = df_2_nr[selection]

# print(df_3.head())

"""
Save dataframe to CSV file
"""
df_3.to_csv('bioactivity_data_2.csv', index=False)

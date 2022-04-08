import pandas as pd
from chembl_webresource_client.new_client import new_client

"""
Target search for coronavirus
"""
target = new_client.target
target_query = target.search('acetylcholinesterase')
targets = pd.DataFrame.from_dict(target_query)

# print(targets.head())

"""
   - Select and retrieve bioactivity data for 
     Human Acetylcholinesterase (first entry)

   - We will assign the fifth entry 
     (which corresponds to the target protein, 
     Human Acetylcholinesterase) to the 
     selected_target variable
"""
selected_target = targets.target_chembl_id[0]

# print(selected_target) # CHEMBL 220

"""
Here, we will retrieve 'only' bioactivity 
data for Human Acetylcholinesterase 
(CHEMBL220) that are reported as pChEMBL values.
"""
activity = new_client.activity
res = activity.filter(
    target_chembl_id=selected_target).filter(
    standard_type="IC50")

df = pd.DataFrame.from_dict(res)

# print(df.head(3))

"""
Finally we will save the resulting bioactivity 
data to a CSV file bioactivity_data.csv.
"""
df.to_csv('bioactivity_data_1.csv', index=False)

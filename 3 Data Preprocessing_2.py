import pandas as pd

path = "C:/Users/HP/PycharmProjects/MachineLearningEnow/Comp_Drug_Discovery/bioactivity_data_2.csv"
df_4 = pd.read_csv(path)
# print(df_4.head())

"""
Labeling compounds as either being: 
active, inactive or intermediate

The bioactivity data is in the IC50 unit. 
Compounds having values of less than 1000 nM 
will be considered to be active while 
those greater than 10,000 nM will be considered 
to be inactive. As for those values in between 
1,000 and 10,000 nM will be referred to as 
intermediate.
"""
bioactivity_threshold = []
for i in df_4.standard_value:
    if float(i) >= 10000:
        bioactivity_threshold.append("inactive")
    elif float(i) <= 1000:
        bioactivity_threshold.append("active")
    else:
        bioactivity_threshold.append("intermediate")

bioactivity_class = pd.Series(
    bioactivity_threshold, name='class')
df_5 = pd.concat(
    [df_4, bioactivity_class], axis=1)

# print(df_5.head())

"""
Save data to CSV file
"""
df_5.to_csv('labeled_curated_data.csv', index=False)

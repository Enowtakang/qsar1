import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import VarianceThreshold


"""
Load the dataset
"""
path = "https://github.com/dataprofessor/data/raw/master/acetylcholinesterase_06_bioactivity_data_3class_pIC50_pubchem_fp.csv"

df = pd.read_csv(path)
# print(df.head())

"""
Define Input features

The dataset contains 307 input features
and 1 output variable (pIC50).
Remember to drop the 'Name' variable
"""
drop_variable = 'pIC50'
X = df.drop(drop_variable, axis=1)

# convert all values in all columns of the new X
# dataframe to type 'float'
# X = X.astype(float)

# print(X.head())

"""
Output features
"""
y = df.pIC50
# print(y.head())

"""
Examine data dimensions
"""
def data_dimensions():
    print(X.shape)
    print(y.shape)

# data_dimensions()

"""
Remove low variance features

Only 12 features would be left!
# """
selection = VarianceThreshold(
    threshold=(.8 * (1 - .8)))
X = selection.fit_transform(X)
# print(X.shape)    # 12 features left

"""
Data split (80/20 ratio)
"""
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2)

def train_test_data_shapes():
    print(f"X_train.shape, y_train.shape = {X_train.shape}, {y_train.shape}")
    print(f"X_test.shape, y_test.shape = {X_test.shape}, {y_test.shape}")

# train_test_data_shapes()

"""
Building a Regression Model using 
Random Forest regressor algorithm
"""
np.random.seed(100)
model = RandomForestRegressor(n_estimators=100)
model.fit(X_train, y_train)
r2 = model.score(X_test, y_test)
# print(r2)

"""
Get all predicted values
"""
Y_pred = model.predict(X_test)

"""
Scatter Plot of Experimental vs Predicted 
pIC50 Values
"""
def scatter_last():
    sns.set(color_codes=True)
    sns.set_style("white")

    ax = sns.regplot(y_test, Y_pred,
                     scatter_kws={'alpha':0.4})
    ax.set_xlabel('Experimental pIC50',
                  fontsize='large', fontweight='bold')
    ax.set_ylabel('Predicted pIC50',
                  fontsize='large', fontweight='bold')
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 12)
    ax.figure.set_size_inches(5, 5)
    plt.show()

# scatter_last()

"""
It may be necessary to run this
script in 'google collaborator',
since PyCharm summarizes results tables
in such a way that they may not be able
to be copied and pasted in a word document.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import VarianceThreshold
from lazypredict.Supervised import LazyRegressor


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
Define and build the lazy classifier
"""

clf = LazyRegressor(verbose=0, ignore_warnings=True, custom_metric=None)
models_train, predictions_train = clf.fit(X_train, X_train, y_train, y_train)
models_test, predictions_test = clf.fit(X_train, X_test, y_train, y_test)

# print('Predictions_train > ')
# print('# Performance table of the training set (80% subset) > ')
# print(predictions_train)
# # try code below in collaborator?
# print('')
# print('Predictions_test > ')
# print('# Performance table of the test set (20% subset) > ')
# print(predictions_test)



"""
Model Performance Visualization
1. R-squared
"""


def view_r2_performance():
    # Bar plot of R-squared values
    plt.figure(figsize=(5, 10))
    sns.set_theme(style="whitegrid")
    ax = sns.barplot(y=predictions_train.index, x="R-Squared", data=predictions_train)
    ax.set(xlim=(0, 1))
    plt.show()


# view_r2_performance()

"""
Model Performance Visualization
2. RMSE
"""


def view_rmse_performance():
    # Bar plot of RMSE values
    plt.figure(figsize=(5, 10))
    sns.set_theme(style="whitegrid")
    ax = sns.barplot(y=predictions_train.index, x="RMSE", data=predictions_train)
    ax.set(xlim=(0, 10))
    plt.show()


# view_rmse_performance()

"""
Model Performance Visualization
3. Calculation Time
"""


def view_calculation_time():
    # Bar plot of calculation time
    plt.figure(figsize=(5, 10))
    sns.set_theme(style="whitegrid")
    ax = sns.barplot(y=predictions_train.index, x="Time Taken", data=predictions_train)
    ax.set(xlim=(0, 10))
    plt.show()


# view_calculation_time()

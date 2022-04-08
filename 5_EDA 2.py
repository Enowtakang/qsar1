"""
Exploratory Data Analysis
(Chemical Space Analysis)
(Chemical Space Analysis)
via Lipinski descriptors
"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


path = "C:/Users/HP/PycharmProjects/MachineLearningEnow/Comp_Drug_Discovery/pIC50_no_intermediate.csv"
df_2_class = pd.read_csv(path)
# print(df_2_class.head())

"""
Frequency plot of the 2 bioactivity classes
"""
def fig_1():
    sns.countplot(
        x='class', data=df_2_class,
        edgecolor='black')

    plt.xlabel('Bioactivity class',
               fontsize=14, fontweight='bold')
    plt.ylabel('Frequency',
               fontsize=14, fontweight='bold')

    # plt.savefig('plot_bioactivity_class.pdf')

    plt.show()

# fig_1()

"""
Scatter plot of MW versus LogP

It can be seen that the 2 bioactivity 
classes are spanning similar chemical 
spaces as evident by the scatter plot of MW vs LogP.
"""
def fig_2():
    plt.figure(figsize=(5.5, 5.5))

    sns.scatterplot(x='MW', y='LogP',
                    data=df_2_class,
                    hue='class', size='pIC50',
                    edgecolor='black', alpha=0.7)

    plt.xlabel('MW', fontsize=14,
               fontweight='bold')
    plt.ylabel('LogP', fontsize=14,
               fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc=2, borderaxespad=0)
    # plt.savefig('plot_MW_vs_LogP.pdf')

    plt.show()

# fig_2()

"""
Box plots:

pIC50 value
"""
def fig_3():
    plt.figure(figsize=(5.5, 5.5))

    sns.boxplot(x='class', y='pIC50',
                data=df_2_class)

    plt.xlabel('Bioactivity class',
               fontsize=14, fontweight='bold')
    plt.ylabel('pIC50 value',
               fontsize=14, fontweight='bold')

    # plt.savefig('plot_ic50.pdf')

    plt.show()

# fig_3()

"""
Statistical analysis | Mann-Whitney U Test
"""
def mannwhitney(descriptor, verbose=False):
    # https://machinelearningmastery.com/
    # nonparametric-statistical-significance-
    # tests-in-python/
    from numpy.random import seed
    from numpy.random import randn
    from scipy.stats import mannwhitneyu

    # seed the random number generator
    seed(1)

    # actives and inactives
    selection = [descriptor, 'class']
    df = df_2_class[selection]
    active = df[df['class'] == 'active']
    active = active[descriptor]

    selection = [descriptor, 'class']
    df = df_2_class[selection]
    inactive = df[df['class'] == 'inactive']
    inactive = inactive[descriptor]

    # compare samples
    stat, p = mannwhitneyu(active, inactive)
    # print('Statistics=%.3f, p=%.3f' % (stat, p))

    # interpret
    alpha = 0.05
    if p > alpha:
        interpretation = 'Same distribution (fail to reject H0)'
    else:
        interpretation = 'Different distribution (reject H0)'

    results = pd.DataFrame({'Descriptor': descriptor,
                            'Statistics': stat,
                            'p': p,
                            'alpha': alpha,
                            'Interpretation': interpretation}, index=[0])
    filename = 'p_mannwhitneyu_' + descriptor + '.csv'
    results.to_csv(filename)

    return results


# mannwhitney('pIC50')

"""
Molecular Weight
"""
def fig_4():
    plt.figure(
        figsize=(5.5, 5.5))
    sns.boxplot(x='class', y='MW',
                data=df_2_class)
    plt.xlabel('Bioactivity class',
               fontsize=14, fontweight='bold')
    plt.ylabel('MW', fontsize=14,
               fontweight='bold')
    # plt.savefig('plot_MW.pdf')

    plt.show()

# fig_4()

"""
Mann-Whitney Molecular Weight
"""
# mannwhitney('MW')

"""
LogP
"""
def fig_5():
    plt.figure(
        figsize=(5.5, 5.5))
    sns.boxplot(x='class',
                y='LogP', data=df_2_class)
    plt.xlabel('Bioactivity class',
               fontsize=14, fontweight='bold')
    plt.ylabel('LogP', fontsize=14,
               fontweight='bold')
    # plt.savefig('plot_LogP.pdf')

    plt.show()

# fig_5()

"""
Statistical analysis | Mann-Whitney U Test
"""
# mannwhitney('LogP')

"""
NumHDonors
"""
def fig_6():
    plt.figure(figsize=(5.5, 5.5))
    sns.boxplot(x='class', y='NumHDonors',
                data=df_2_class)
    plt.xlabel('Bioactivity class', fontsize=14,
               fontweight='bold')
    plt.ylabel('NumHDonors', fontsize=14,
               fontweight='bold')
    # plt.savefig('plot_NumHDonors.pdf')
    plt.show()

# fig_6()

"""
Statistical analysis | Mann-Whitney U Test
"""
# mannwhitney('NumHDonors')

"""
NumHAcceptors
"""
def fig_7():
    plt.figure(figsize=(5.5, 5.5))
    sns.boxplot(x='class', y='NumHAcceptors',
                data=df_2_class)
    plt.xlabel('Bioactivity class',
               fontsize=14, fontweight='bold')
    plt.ylabel('NumHAcceptors',
               fontsize=14, fontweight='bold')
    # plt.savefig('plot_NumHAcceptors.pdf')

    plt.show()

# fig_7()

"""
Statistical analysis | Mann-Whitney U Test
"""
# mannwhitney('NumHAcceptors')

"""
Interpretation of Statistical Results

Box Plots

pIC50 values:

Taking a look at pIC50 values, the actives 
and inactives displayed statistically 
significant difference, which is to be 
expected since threshold values 
(IC50 < 1,000 nM = Actives while IC50 > 10,000 
nM = Inactives, corresponding to 
pIC50 > 6 = Actives and pIC50 < 5 = Inactives) 
were used to define actives and inactives.


Lipinski's descriptors:

All of the 4 Lipinski's descriptors exhibited 
statistically significant difference between 
the actives and inactives.

"""

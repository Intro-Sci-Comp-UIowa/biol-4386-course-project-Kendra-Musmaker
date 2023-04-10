Making a bar plot of frequency bt genotype with error bars

import pandas as pd
import math
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# specify the file path
file_path = '/data_for_comp_sci.csv'

# read the CSV file using pandas read_csv() function
df = pd.read_csv(file_path)

# print the first 5 rows of the dataframe
#print(df)


# change the data type of the 'my_column' column from int to float
df = df.astype({'survivors': 'float'})
df = df.astype({'cells_plated': 'float'})

df['Individual_Frequency_true'] = df.survivors/df.cells_plated
df['Individual_Frequency'] = df.survivors/df.cells_plated * (10**10)

# print the first 5 rows of the dataframe
print(df.head())

# 'msh2Δtlc1Δ', 'tlc1Δ', 'mph1Δtlc1Δ', 'mlh1Δtlc1Δ'

genotypes_set = (set(df['genotype'].values))


counter = 0
for genotype in genotypes_set:
  # sum the values in 'my_column' where 'my_string_column' contains 'foo'
  total_survivors = df.loc[df['genotype'].str.contains(genotype), 'survivors'].sum()
  #print(total_survivors)

  total_cells = df.loc[df['genotype'].str.contains(genotype), 'cells_plated'].sum()
  #print(total_cells)

  total_frequency = total_survivors/total_cells
  #print(total_frequency)

  critical_value = 1.96 * math.sqrt( (total_frequency*(1-total_frequency)) / total_cells ) 

  lower_ci = total_frequency - critical_value
  upper_ci = total_frequency + critical_value

  adjusted_total_frequency = total_frequency * (10**10)
  adj_lower_ci = lower_ci * (10**10)
  adj_upper_ci = upper_ci * (10**10)
  error_bar = adj_upper_ci - adj_lower_ci

  if counter == 0:
    new_row = ({'Genotype': genotype,
                'total_survivors': total_survivors,
                'total_cells': total_cells,
                'frequency_true': total_frequency,
                'Frequency (10^10)': adjusted_total_frequency,
                'error': error_bar,
                'Lower_CI_true': lower_ci,
                'Upper_CI_true': upper_ci,
                'Lower_CI': adj_lower_ci,
                'Upper_CI': adj_upper_ci})
    df_summary = pd.DataFrame(data=new_row, index=[0])
    counter += 1
  else:
    new_row =  pd.Series({'Genotype': genotype,
                          'total_survivors': total_survivors,
                          'total_cells': total_cells,
                          'frequency_true': total_frequency,
                          'Frequency (10^10)': adjusted_total_frequency,
                          'error': error_bar,
                          'Lower_CI_true': lower_ci,
                          'Upper_CI_true': upper_ci,
                          'Lower_CI': adj_lower_ci,
                          'Upper_CI': adj_upper_ci})
    df_summary = pd.concat([df_summary, new_row.to_frame().T], ignore_index=True)
#print(df_summary)

def plot_errorbars(arg, **kws):
    np.random.seed(sum(map(ord, "error_bars")))
    x = np.random.normal(0, 1, 100)
    f, axs = plt.subplots(2, figsize=(7, 2), sharex=True, layout="tight")
    sns.pointplot(x=x, errorbar=arg, **kws, capsize=.3, ax=axs[0])
    sns.stripplot(x=x, jitter=.3, ax=axs[1])

sns.set(style="whitegrid")
ax = sns.barplot(data=df_summary, x='Genotype', y='Frequency (10^10)', errcolor='k', errwidth=10, errorbar=('ci', 95))
x_coords = [p.get_x() + 0.5 * p.get_width() for p in ax.patches]
y_coords = [p.get_height() for p in ax.patches]
#plot_errorbars(lambda x: (x.min(), x.max()))
ax.errorbar(x=x_coords, y=y_coords, yerr=df_summary["error"], fmt="none", c="k")
ax.set(ylim=(0,8000000))

sns.scatterplot(data=df, x='genotype', y='Individual_Frequency', size=4, color='w', edgecolor='k', linewidth=1, legend=False)

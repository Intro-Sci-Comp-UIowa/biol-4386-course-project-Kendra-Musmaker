Text for the report
In this figure, the strain AM3692(TLC1/tlc1::BSD) and its derivatives are used. AM3692 is heterozygous for a deletion of the telomerase component TLC1. Tetrad dissection, in which haploid spores from a single mother cell are separated, is used to isolate haploid tlc1Δ spores with/without the deletion of a second gene of interest. These haploid tlc1Δ colonies are then started in 1 ml YEPD (yeast extract-peptone-dextrose) overnight, after which they are passaged every 24 hours in 4 ml YEPD until senescence is reached. Senescence is reached after 2-3 passages depending on the growth pattern of the strain, with cells being plated after passage 1 or 2, respectively. Cultures with RAD52, RAD51 and SRS2 deficient strains require plating after the first passage, all others can be plated after passage 2. The number of cell/ml passaged is determined by the growth pattern of the strain, with cultures that divide poorly having higher numbers of cells passaged in order to allow enough cells to be plated to detect the formation of survivors. In all cases, cells passaged is kept low enough that only one survivor is expected to form per culture. This allows for the researcher to be confident that survivors are independent from one another. The total number of survivors for a genotype is monitored and then divided by the total adjusted cells plated to determine a frequency for ALT.

Rather than recreating this figure, I am interested in the effect of deleting mismatch repair genes such as MSH2 and MLH1. In the main branch of the github repository biol-4386-course-project-Kendra-Musmaker/Data at main · Intro-Sci-Comp-UIowa/biol-4386-course-project-Kendra-Musmaker (github.com)is a file data_for_comp_sci.Each row in this file is an individual culture. Column 1A is strain number from which the culture was derived, column B is the genotype of the culture. Column C is the survivors counted for that colony. Column D I the adjusted cells plated for that colony, where adjusted cells plated is calculated as the cells plated multiplied by 2(remaining population doublings). This file contains cultures from strain 6169 (TLC1/tlc1::BSD MSH2/msh2::kan), 6205 TLC1/tlc1::BSD MLH1/mlh1::kan) and AM6212 (TLC1/tlc1::BSD MPH1/mph1::kan). These cultures were obtained from experiments run between 2020 and early 2023.

To plot overall frequency for each genotype a table df_summary was created to have total survivors (total_survivors), total cells plated (total_cells), frequency, lower and upper CI and error for each genotype using the following code. Total_frequency is calculated as total survivors/total cells plated, while lower_CI and upper_CI are calculated as 𝑓±1.96√((𝑓(1−𝑓)𝑛)), respectively, with f=total_frequency and n=total cells. Adjusted _total frequencies, adj_lower_ci and adj_upper_ci are the previous values multiplied by 10^10 and are used for graphing purposes.

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

new_row = pd.Series({'Genotype': genotype,

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

Adjusted_total_frequency was then plotted using seaborn with genotypes along the x-axis and frequency (x10^10) as the y axis. Error bars were set to equal ‘error’ from df_summary, with ‘error’ being adj_upper_ci – adj_lower_ci. The following code was used for this step.

sns.set(style="whitegrid")

ax = sns.barplot(data=df_summary, x='Genotype', y='Frequency (10^10)', errcolor='k', errwidth=10, errorbar=('ci', 95))

x_coords = [p.get_x() + 0.5 * p.get_width() for p in ax.patches]

y_coords = [p.get_height() for p in ax.patches]

#plot_errorbars(lambda x: (x.min(), x.max()))

ax.errorbar(x=x_coords, y=y_coords, yerr=df_summary["error"], fmt="none", c="k")

ax.set(ylim=(0,8000000))

The next steps for this data is to overlay the individual frequencies on top of the bars for overall frequency. To do this individual frequencies will need to be multiplied by 10^10. The graph will also be reordered so that tlc1Δ is first on the graph. 

Making a bar plot of frequency bt genotype with error bars code

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

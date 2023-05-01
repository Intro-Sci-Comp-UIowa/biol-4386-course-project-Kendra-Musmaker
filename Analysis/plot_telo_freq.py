import pandas as pd
import math
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

###PARAMS###
DF_PATH = '../Data/data_for_comp_sci.csv'
OUTPUT_F_NAME = 'expt_1_summary.csv'
OUTPUT_FIG_NAME = 'plot_telomer.png'
SHOW_FIG = False
###END_PARAMS###

def df_analysis(df_path):

  # read the CSV file using pandas read_csv() function
  df = pd.read_csv(df_path)

  # print the first 5 rows of the dataframe
  #print(df)

  # change the data type of the 'my_column' column from int to float
  df = df.astype({'survivors': 'float'})
  df = df.astype({'cells_plated': 'float'})

  df['Individual_Frequency_true'] = df.survivors/df.cells_plated
  df['Individual_Frequency'] = df.survivors/df.cells_plated * (10**10)

  # print the first 5 rows of the dataframe
  print(f"Read data from {df_path}. Here is a peek:")
  print(df.head())

  # 'msh2Δtlc1Δ', 'tlc1Δ', 'mph1Δtlc1Δ', 'mlh1Δtlc1Δ'

  genotypes_set = (set(df['genotype'].values))

  counter = 0
  for genotype in genotypes_set:
    print(f"Analyzing {genotype}")
    # sum the values in 'my_column' where 'my_string_column' contains 'foo'
    total_survivors = df.loc[df['genotype'].str.contains(genotype), 'survivors'].sum()
    #print(total_survivors)

    total_cells = df.loc[df['genotype'].str.contains(genotype), 'cells_plated'].sum()
    #print(total_cells)

    total_frequency = total_survivors/total_cells
    #print(total_frequency)
    #determining 95% confidence interval
    #critical width is half the 95% CI
    critical_value = 1.96 * math.sqrt( (total_frequency*(1-total_frequency)) / total_cells ) 
    #defining lower and upper bounds
    lower_ci = total_frequency - critical_value
    upper_ci = total_frequency + critical_value
    #multipplying by 10^10 for graphing purposes. Setting error to be the width of the interval
    adjusted_total_frequency = total_frequency * (10**10)
    adj_lower_ci = lower_ci * (10**10)
    adj_upper_ci = upper_ci * (10**10)
    error_bar = adj_upper_ci - adj_lower_ci

    #creating a datframe of summary statistics for the genotypes
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

  return(df_summary, df)

def plot_fig(df_summary, df, output_f_name, fig_f_name="plot_telomer.png",show=False):
  #creating a barplot of total frequency for each genotype with error bars being our 95% CI
  sns.set(style="whitegrid")
  ax = sns.barplot(data=df_summary, x='Genotype', y='Frequency (10^10)', errcolor='k', errwidth=10, errorbar=('ci', 95))
  x_coords = [p.get_x() + 0.50 * p.get_width() for p in ax.patches]
  y_coords = [p.get_height() - 6000 for p in ax.patches]
  ax.errorbar(x=x_coords, y=y_coords, yerr=df_summary["error"], fmt="none", c="k")
  ax.set(ylim=(0,2000000))
  #plotting individual frequencies on top of the bar graph of overall frequencies
  sns.scatterplot(data=df, x='genotype', y='Individual_Frequency', size=200, color='w', edgecolor='k', linewidth=1, legend=False)

  plt.savefig(fig_f_name, dpi=300)
  print(f"Saved figure as {fig_f_name}")

  if show:
    plt.show()

  df_summary.to_csv(output_f_name, sep='\t')

if __name__ == "__main__":

  # specify the file path
  df_path = DF_PATH
  output_f_name = OUTPUT_F_NAME
  df_summary, df = df_analysis(df_path)
  plot_fig(df_summary, df, output_f_name, fig_f_name=OUTPUT_FIG_NAME, show=SHOW_FIG)

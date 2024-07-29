#!/usr/bin/env python
# coding: utf-8

# # Pymaceuticals Inc.
# ---
# 
# ### Analysis
# 
# - For mouse l509, the overall tumor volume decreased during its treatment with Capomulin. But the decrease was not linear. It initially increased
# until a steep dropoff at around 20 days. The volume continued to on and off decrease. But with dropoffs being greater than increases.
# - Mice taking Ketapril had the largest tumors on average, with a volume of 55.23(mm3). While mice taking Capomulin had the smallest, with an average    volume of 40.68(mm3).
# - The overall mouse dataset was nearly equal in terms of genders tested. With males making up 51% percent, and Females making up the other 49% of
# the 249 total mice.
#  

# In[21]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
import os

# Study data files
mouse_metadata_path = "./data/Mouse_metadata.csv"
study_results_path = "./data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single DataFrame
completed_study_df = pd.merge(study_results, mouse_metadata, on = "Mouse ID", how = "left")
#completed_study_df
# Display the data table for preview


# In[23]:


# Checking the number of mice.
num_mice = len(completed_study_df["Mouse ID"].unique())
print(f"Number of mice: {num_mice}")


# In[25]:


# Our data should be uniquely identified by Mouse ID and Timepoint
# Get the duplicate mice by ID number that shows up for Mouse ID and Timepoint.
dup_mouse_id_num = completed_study_df[completed_study_df.duplicated(subset = ["Mouse ID", 
                                                           "Timepoint"])]["Mouse ID"].unique()
dup_mouse_id_num


# In[11]:


# Optional: Get all the data for the duplicate mouse ID.
dup_mouse_data = completed_study_df[completed_study_df['Mouse ID'] == 'g989']
dup_mouse_data


# In[13]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_completed_study_df = completed_study_df[completed_study_df['Mouse ID'].isin(dup_mouse_id_num) == False]
## Removed all duplicated values
clean_completed_study_df


# In[15]:


# Checking the number of mice in the clean DataFrame.
len(clean_completed_study_df['Mouse ID'].unique())


# ## Summary Statistics

# In[17]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen:
# mean, median, variance, standard deviation, and SEM of the tumor volume.
# Assemble the resulting series into a single summary DataFrame.
mean_tumor_vol_per_drug = clean_completed_study_df.groupby('Drug Regimen').agg({'Tumor Volume (mm3)': 'mean'})
mean_tumor_vol_per_drug


# In[19]:


# A more advanced method to generate a summary statistics table of mean, median, variance, standard deviation,
# and SEM of the tumor volume for each regimen (only one method is required in the solution)

# Using the aggregation method, produce the same summary statistics in a single line
tumor_stats = clean_completed_study_df.groupby('Drug Regimen').agg({'Tumor Volume (mm3)':
                                                                    ['mean', 'median', 'var',
                                                                     'std', 'sem']})
summary_table = tumor_stats.round(2)
summary_table


# ## Bar and Pie Charts

# In[45]:


# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using Pandas.
dr_counts = clean_completed_study_df['Drug Regimen'].value_counts()
#dr_counts
dr_counts.plot(kind='bar')
plt.xlabel('Drug Regimen')
plt.xticks(rotation=45)
plt.ylabel('Number of Mice Tested')
plt.show


# In[43]:


# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using pyplot.
## Value is reset
dr_counts = clean_completed_study_df['Drug Regimen'].value_counts()
plt.bar(dr_counts.index.values, dr_counts.values)
plt.xlabel('Drug Regimen')
plt.xticks(rotation=45)
plt.ylabel('Number of Mice Tested')
plt.show


# In[53]:


# Generate a pie chart, using Pandas, showing the distribution of unique female versus male mice used in the study

# Get the unique mice with their gender
dr_counts = clean_completed_study_df.Sex.value_counts()
#dr_counts

# Make the pie chart
dr_counts.plot(kind='pie', autopct='%1.1f%%')


# In[61]:


# Generate a pie chart, using pyplot, showing the distribution of unique female versus male mice used in the study

# Get the unique mice with their gender
dr_counts = clean_completed_study_df.Sex.value_counts()

# Make the pie chart
plt.pie(dr_counts.values, labels=dr_counts.index.values, autopct='%1.1f%%')
plt.ylabel('Sex')
plt.show


# ## Quartiles, Outliers and Boxplots

# In[73]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:
# Capomulin, Ramicane, Infubinol, and Ceftamin

# Start by getting the last (greatest) timepoint for each mouse
largest_tumor = clean_completed_study_df.groupby(['Mouse ID'])['Timepoint'].max()
#largest_tumor
largest_tumor = largest_tumor.reset_index()
# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
merged_data = largest_tumor.merge(clean_completed_study_df, on=['Mouse ID', 'Timepoint'], how='left')
#merged_data


# In[91]:


# Put treatments into a list for for loop (and later for plot labels)
treatments = ['Capomulin', 'Ramicane', 'Infubinol', 'Ceftamin']

# Create empty list to fill with tumor vol data (for plotting)
tumor_vol_data = []

# Calculate the IQR and quantitatively determine if there are any potential outliers.
for drug in treatments:

    # Locate the rows which contain mice on each drug and get the tumor volumes
    final_tumor_vol = merged_data.loc[merged_data['Drug Regimen'] == drug, 'Tumor Volume (mm3)']
    #final_tumor_vol
    # add subset
    tumor_vol_data.append(final_tumor_vol)

    # Determine outliers using upper and lower bounds
    quartiles = final_tumor_vol.quantile([.25, .5, .75])
    lower_quartile = quartiles[.25]
    upper_quartile = quartiles[.75]
    iqr = upper_quartile-lower_quartile
    lower_bound = lower_quartile-(1.5*iqr)
    upper_bound = upper_quartile+(1.5*iqr)
    outliers = final_tumor_vol.loc[(final_tumor_vol < lower_bound) | (final_tumor_vol > upper_bound)]
    print(f"{drug}'s potential outliers {outliers}")


# In[99]:


# Generate a box plot that shows the distribution of the tumor volume for each treatment group.
tumor_distr = dict(markerfacecolor='red', markersize=10)
plt.boxplot(tumor_vol_data, labels = treatments, flierprops = tumor_distr)
plt.ylabel("Final Tumor Volume (mm3)")
plt.show


# ## Line and Scatter Plots

# In[129]:


# Generate a line plot of tumor volume vs. time point for a single mouse treated with Capomulin
capomulin_data = clean_completed_study_df[clean_completed_study_df['Drug Regimen'] == "Capomulin"]
#capomulin
mouse_data = capomulin_data[capomulin_data['Mouse ID'] == "l509"]

plt.plot(mouse_data['Timepoint'], mouse_data['Tumor Volume (mm3)'])
plt.xlabel("Timepoint (days)")
plt.ylabel("Tumor Volume (mm3)")
plt.title("Capomulin Treatment of Mouse l509")


# In[161]:


# Generate a scatter plot of mouse weight vs. the average observed tumor volume for the entire Capomulin regimen
capomulin_table = clean_completed_study_df[clean_completed_study_df['Drug Regimen'] == "Capomulin"]
#capomulin_table
#capomulin_table.groupby(["Mouse ID"])
capomulin_av = capomulin_table.groupby('Mouse ID').agg({'Tumor Volume (mm3)': 'mean',
                                               'Timepoint': 'mean',
                                              'Metastatic Sites': 'mean',
                                              'Age_months': 'mean',
                                              'Weight (g)': 'mean'})
#capomulin_av
plt.scatter(capomulin_av['Weight (g)'], capomulin_av['Tumor Volume (mm3)'])
plt.xlabel('Weight (g)')
plt.ylabel('Average Tumor Volume (mm3)')



# ## Correlation and Regression

# In[183]:


# Calculate the correlation coefficient and a linear regression model
# for mouse weight and average observed tumor volume for the entire Capomulin regimen
correlation = st.pearsonr(capomulin_av['Weight (g)'], capomulin_av['Tumor Volume (mm3)'])

print(f"The correlation between mouse weight and the average tumor volume is {round(correlation[0],2)}")
chart = st.linregress(capomulin_av['Weight (g)'], capomulin_av['Tumor Volume (mm3)'])
#chart
slope = chart[0]
b = chart[1]
y_values = capomulin_av['Weight (g)'] * slope + b
plt.scatter(capomulin_av['Weight (g)'], capomulin_av['Tumor Volume (mm3)'])
plt.plot(capomulin_av['Weight (g)'], y_values, color="red")
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.show


# In[ ]:





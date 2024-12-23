"""
## About


STAGE 1: Unique taxanomy analysis.
STEP  2: identification and tabulation


This code is to find the core, sub core microbes and cross check with Chennai to find which microbes are prevalent and were are they found
"""
import pandas as pd
import os
import copy
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.ticker as ticker
import warnings
import json

warnings.filterwarnings("ignore")
'''
## Lets identify the species of interest to Chennai
FIND = os.path.join("findings")
CHENNAI = os.path.join("processed_data", "chennai_taxa_num")
IMG_OUTPUT = os.path.join("images", "core_microbe")
OTHER_CITY_TAXA = os.path.join(
    "dataset", "MetaSUB_core", "full", "processed_taxa_num", "merged_taxa"
)

species_df = pd.read_csv(os.path.join(CHENNAI, "filtered_species.csv"))
print("Loaded the Chennai species dataframe")
## prevalance
chennai_core_df, chennai_sub_core_df, chennai_peripheral_df = prevalance_extractor(
    species_df, bounds=[0.97, 0.7]
)

chennai_core_df.to_csv(os.path.join(FIND, "chennai_core.csv"), index=False)
chennai_sub_core_df.to_csv(os.path.join(FIND, "chennai_sub_core.csv"), index=False)
chennai_peripheral_df.to_csv(os.path.join(FIND, "chennai_peripheral.csv"), index=False)


print("Calculated core sub-core and peripheral microbes of Chennai samples")
## other cities
other_city_species = pd.read_csv(
    os.path.join(OTHER_CITY_TAXA, "merged_species_filtered.csv")
)
print("Loaded the Other city species dataframe")
other_city_species = other_city_species.drop(
    columns=[
        "uuid",
        "core_project",
        "project",
        "surface_material",
        "continent",
        "surface_ontology_fine",
        "control_type",
    ]
)

other_city = list(other_city_species.city.value_counts().index.values)
"""
other_city = [
    "singapore",
    "porto",
    "london",
    "hong_kong",
    "fairbanks",
    "new_york_city",
]
"""
print("Total No of Other Cities", len(other_city))

city_taxa_dict = {}
city_taxa_dict["chennai"] = {
    "core": set(list(chennai_core_df["Species"].values)),
    "sub_core": set(list(chennai_sub_core_df["Species"].values)),
    "peripheral": set(list(chennai_peripheral_df["Species"].values)),
}
for city in other_city:
    temp = other_city_species.query("city==@city")
    temp = temp.drop(
        columns=[
            "city",
        ]
    )
    core_df, sub_core_df, peripheral_df = prevalance_extractor(temp, bounds=[0.97, 0.7])

    city_taxa_dict[city] = {
        "core": set(list(core_df["Species"].values)),
        "sub_core": set(list(sub_core_df["Species"].values)),
        "peripheral": set(list(peripheral_df["Species"].values)),
    }
print("Calculated core sub-core and peripheral microbes of Other city samples")
'''
"""
### Q1

if we combine all the cities core microbes, does chennai has anything unique in its core microbes?
"""

DATASET = os.path.join(
    "dataset",
)
CHENNAI_DATASET = os.path.join(DATASET, "chennai_v1_data")
#THRESHOLD = 0.01
FOLDER = '0001'
FINDINGS = os.path.join(
    "findings",
    "chennai_special_core_4",
    FOLDER
)

def filter_all_zero(df1, col_not_to_include=["Samples"]):
    df = copy.copy(df1)
    col_list = [i for i in df.columns if i not in col_not_to_include]
    temp = df[col_list].sum().to_frame()#.reset_index()
    temp.rename(columns={0: "non_zeros"}, inplace=True)
    print("before ", df.shape)
    print(temp.head())
    df = df.drop(
        columns=list(temp.query("non_zeros<=0.000001").index.values)
    )  ## The name of the species which are not having any read mapped. These columns can be dropped

    print("after ", df.shape)
    #print('Removed',list(temp.query("non_zeros<=0.000001").index.values))
    return df

city_taxa_dict = json.load(open(os.path.join(FINDINGS,'all_city_core.json')))
other_city = [i for i in city_taxa_dict.keys() if i != 'Chennai']

print("############################################################################")
print("Working on Q1")
## join all city cores

level = "core"
all_cities = other_city + ["Chennai"]
unique_taxas_dict = {}
for city_of_int in all_cities:
    all_other_city_core_set = set()
    for city in other_city:
        if city != city_of_int:
            all_other_city_core_set = all_other_city_core_set.union(
                city_taxa_dict[city][level]
            )

    c_difference = list(
        set(city_taxa_dict[city_of_int][level]).difference(all_other_city_core_set)
    )
    unique_taxas_dict[city_of_int] = c_difference

    print(
        f" Total {city_of_int.capitalize()} {level} ",
        len(city_taxa_dict[city_of_int]["core"]),
        "Where",
        c_difference.__len__(),
        "are unique to the city",
    )
temp = pd.DataFrame.from_dict({"Chennai_spl": unique_taxas_dict["Chennai"]})
temp.to_csv(os.path.join(FINDINGS, "chennai_spl_core_only_core_considered.csv"), index=False)
"""
### Q2

if we combine all the cities core microbes and sub core, does chennai has anything unique in its core microbes?
"""
print("############################################################################")
print("Working on Q2")

## join all city cores

level = "core"
unique_taxas_dict = {}
for city_of_int in all_cities:
    all_other_city_core_subcore_set = set()
    for city in other_city:
        if city != city_of_int:
            all_other_city_core_subcore_set = all_other_city_core_subcore_set.union(
                city_taxa_dict[city]["core"]
            )
            all_other_city_core_subcore_set = all_other_city_core_subcore_set.union(
                city_taxa_dict[city]["sub_core"]
            )

    c_difference = list(
        set(city_taxa_dict[city_of_int][level]).difference(all_other_city_core_subcore_set)
    )
    unique_taxas_dict[city_of_int] = c_difference
    print(
        f" Total {city_of_int.capitalize()} {level} ",
        len(city_taxa_dict[city_of_int]["core"]),
        "Where",
        c_difference.__len__(),
        "are unique to the city",
    )
temp = pd.DataFrame.from_dict({"Chennai_spl": unique_taxas_dict["Chennai"]})
temp.to_csv(os.path.join(FINDINGS, "chennai_spl_core_considered_core_subcore.csv"), index=False)


####################### plot the special Chennai core ###############################
species_df = pd.read_csv(os.path.join(CHENNAI_DATASET, 'bracken_species_frac_processed.csv'))
species_df.set_index('Samples', inplace= True)
s_list = species_df.columns

## Threshold the cell 
#def abd_threshold(row):
#    for column in s_list:
#        if row[column] < THRESHOLD:
#            row[column] = 0
#    return row

#species_df = species_df.apply(abd_threshold, axis=1)
#species_df = filter_all_zero(species_df)

unique_df = pd.read_csv(os.path.join(FINDINGS, "chennai_spl_core_considered_core_subcore.csv"))
print(
    "No of unique core species in Chennai samples", len(unique_df["Chennai_spl"].values)
)

unique_species_df = species_df[list(unique_df["Chennai_spl"].values)]
print('================== Special core basic statistics ==================')
print('================== Minimum ==================')
print(unique_species_df.min())
print('================== Maximum ==================')
print(unique_species_df.max())
print('================== Mean ==================')
print(unique_species_df.mean())
print('================== Median ==================')
print(unique_species_df.median())
#print('================== Minimum ==================')
#print(unique_species_df.min())
num_unique_core_species = len(unique_df["Chennai_spl"].values)

if num_unique_core_species % 2 == 0:
    row = int(num_unique_core_species / 2)
    col = 2
else:
    row = int((num_unique_core_species + 1) / 2)
    col = 2


# Assuming your dataframe is named 'df' and the three columns are 'column1', 'column2', and 'column3'
fig, axes = plt.subplots(
    row, col, figsize=(8, 16)
)  # Create subplot with 3 rows and 1 column
# Increase the font size
plt.rcParams.update({"font.size": 18})
x_min = 0
x_max = 0.045
y_min = 0
y_max = 50000
# Loop through each column and create a KDE plot in each subplot
for i, column in enumerate(unique_df["Chennai_spl"].values[:row]):
    ax = axes[i][0]  # Select the current subplot
    sns.kdeplot(
        data=unique_species_df, x=column, ax=ax
    )  # Create KDE plot for the current column
    ax.set_xlim(x_min, x_max)  # Set the x-axis limits
    ax.set_xlabel(column, fontsize=14)   # Set the x-axis label as the column name
    ax.set_ylabel("")
    ax.set_yticks([])  # Remove the y-axis ticks
    ax.set_yticklabels([])  # Remove the y-axis tick labels
    ax.xaxis.set_major_locator(ticker.FixedLocator([x_min, (x_min + x_max) / 2, x_max]))

# Loop through each column and create a KDE plot in each subplot
for i, column in enumerate(unique_df["Chennai_spl"].values[row:]):
    ax = axes[i][1]  # Select the current subplot
    sns.kdeplot(
        data=unique_species_df, x=column, ax=ax
    )  # Create KDE plot for the current column
    ax.set_xlim(x_min, x_max)  # Set the x-axis limits
    ax.set_xlabel(column, fontsize=14)  # Set the x-axis label as the column name
    ax.set_ylabel("")
    ax.set_yticks([])  # Remove the y-axis ticks
    ax.set_yticklabels([])  # Remove the y-axis tick labels
    # Set custom locator for x-axis ticks
    ax.xaxis.set_major_locator(ticker.FixedLocator([x_min, (x_min + x_max) / 2, x_max]))

plt.tight_layout()  # Adjust spacing between subplots
plt.savefig(os.path.join(FINDINGS, "Chennai_unique_core_microbes.png"))
plt.savefig(os.path.join(FINDINGS, "Chennai_unique_core_microbes.svg"), format='svg')
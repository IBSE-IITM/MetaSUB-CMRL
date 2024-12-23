"""
About
: Used Currently
This code is about cleaning the bracken species relative abundance table of Chennai and other cities to find Chennai specific core taxa.
The filtering step removes the non-core samples from MetaSUB 
This produces a json called all_city_core.json whcih has the core sub-core and peri taxa names for the respective prevalance and rel abd filters within the city
"""
import pandas as pd
import os
import copy
import json
from tqdm import tqdm

THRESHOLD = 0.0001 # 0.0001
FOLDER = '0001'
bounds = [0.97, 0.8]
peri = 0.15
MIN_SAMPLES = 12
CITIES_TO_IGNORE = []
DATASET = os.path.join(
    "dataset",
)
FINDINGS = os.path.join(
    "findings",
    "chennai_special_core_4",
    FOLDER
)
os.makedirs(FINDINGS, exist_ok = True)
CHENNAI_DATASET = os.path.join(DATASET, "chennai_v1_data")
CORE_TAXA = os.path.join(
    DATASET,
    "MetaSUB_core",
    "full_w_redown",
)

print('Cleaning for threshold',THRESHOLD)

def prevalance_extractor(df, bounds=[], col_not_to_include = []):
    taxa_df = copy.copy(df)
    taxa_list = [i for i in taxa_df.columns if i not in col_not_to_include]

    for i in taxa_list:
        taxa_df[i] = taxa_df[i] > 0
        taxa_df[i] = taxa_df[i] * 1

    v = taxa_df[taxa_list].sum() / taxa_df.shape[0]
    v = v.to_frame()
    v = v.reset_index()
    v.rename(columns={"index": "Species", 0: "prevalance"}, inplace=True)
    #v["Species"] = v["Species"].str.replace("prev_", "")
    v.sort_values('prevalance', ascending = False, inplace = True)
    core_df = v.query("prevalance > @bounds[0]")
    sub_core_df = v.query("prevalance > @bounds[1] and prevalance < @bounds[0]")
    peripheral_df = v.query("prevalance < @peri and prevalance > 0.00001")
    v = v.rename(columns={'prevalance':'Prevalance'})
    return v, core_df, sub_core_df, peripheral_df

def merger_function(taxa_df, metadata_df):
    taxa_temp = copy.copy(taxa_df)
    taxa_meta = copy.copy(metadata_df)
    taxa_temp = taxa_temp.merge(
        taxa_meta, how="inner", left_on="Samples", right_on="uuid"
    )
    return taxa_temp

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

"""
def core_species(df1, col_not_to_include=["Samples"], core_thresh = 0.97):
    df = copy.copy(df1)
    num_of_samples = df.shape[0]
    col_list = [i for i in df.columns if i not in col_not_to_include]
    ## For prevalance check if > 0 and set to 1 # converted to presence and absence
    for i in col_list:
        df[i] = df[i] > 0 # 
        df[i] = df[i] * 1
    temp = df[col_list].sum().to_frame()#.reset_index()
    temp = temp/num_of_samples
    temp.rename(columns={0: "Prevalance"}, inplace=True)
    #print("before ", temp.shape)
    temp.sort_values('Prevalance', ascending = False, inplace = True)
    temp = temp.query("Prevalance >= @core_thresh")#.index.values
    #print("after ", temp.shape)
    #print('Removed',list(temp.query("non_zeros<=0.000001").index.values))
    return temp

"""
####################### S1 Chennai core species ###############################

chennai_species_df = pd.read_csv(os.path.join(CHENNAI_DATASET, 'bracken_species_non_human_processed.csv'))
chennai_species_df.set_index('Samples', inplace= True)

s_list = chennai_species_df.columns

####################### S2 Abundance Filtering ###############################
## Threshold the cell 
def abd_threshold(row):
    for column in s_list:
        if row[column] < THRESHOLD:
            row[column] = 0
    return row

chennai_species_df = chennai_species_df.apply(abd_threshold, axis=1)
chennai_species_df = filter_all_zero(chennai_species_df)
################################## S2 ########################################

####################### S3 Prevalance Filtering ###############################
#chennai_core_df = core_species(chennai_species_df, col_not_to_include=["Samples"], core_thresh = CORE_PRE)
full_prevalance_df, chennai_core_df, chennai_sub_core_df, chennai_peripheral_df = prevalance_extractor(chennai_species_df, bounds=bounds, col_not_to_include = ['Samples'])
print('Prevalance filtered')
print(chennai_core_df.head())
print(chennai_core_df.tail())
print("subcore")
print(chennai_sub_core_df.head())
print(chennai_sub_core_df.tail())
print("peri")
print(chennai_peripheral_df.head())
print(chennai_peripheral_df.tail())

full_prevalance_df.to_csv(os.path.join(FINDINGS,'species_prevalance.csv'),index=False)
#################################### S3 ######################################

#################################### S1 ######################################

## Dictionary with city as key and core species as values

core_species_dict = {'Chennai': { 
                                    'core':list(chennai_core_df["Species"].values),
                                    'sub_core':list(chennai_sub_core_df["Species"].values),
                                    'peripheral':list(chennai_peripheral_df["Species"].values) 
                                  } 
                                  }

####################### Cont. for all cities ###############################

####################### S1 Filter the Metadata ###############################
METADATA = os.path.join(DATASET,"MetaSUB_core","complete_metadata.csv",)
metadata_df = pd.read_csv(METADATA)
# surface_material has air
# control_type has positive, negative and other control samples
metadata_mini_df = metadata_df[
    [
        "uuid",
        "core_project",
        "project",
        "city",
        "surface_material",
        "continent",
        "surface_ontology_fine",
        "control_type",
    ]
]

## lets remove all the ctrl samples from the metadata
metadata_mini_df["control_type"] = metadata_mini_df["control_type"].fillna(
    "not_control"
)
print("Actual", metadata_mini_df.shape)
metadata_mini_df = metadata_mini_df.query("core_project=='core'")
print("only core", metadata_mini_df.shape)
metadata_mini_df = metadata_mini_df.query("control_type=='not_control'")
print("only not control samples", metadata_mini_df.shape)
metadata_mini_df = metadata_mini_df.query("surface_material!='air'")
print("No Air samples", metadata_mini_df.shape)
metadata_mini_df = metadata_mini_df.query("surface_ontology_fine!='biological'")
print("No biological samples", metadata_mini_df.shape)

#################################### S1 ######################################

####################### S1 City core species ###############################
metasub_species = pd.read_csv(os.path.join(CORE_TAXA, 'bracken_species_non_human_processed.csv'))
s_list = [i for i in metasub_species.columns if i != 'Samples']
print('Core apply threshold')
metasub_species = metasub_species.apply(abd_threshold, axis=1)
print('Core filter all zero')
metasub_species = filter_all_zero(metasub_species)
print('Core merge metadata')
metasub_species = merger_function(metasub_species, metadata_mini_df)
print('########## do we have non core project sampels ###############')
print(metasub_species['core_project'].value_counts())
print(metasub_species['city'].value_counts())
## Filter the cities based on number of samples
vc = metasub_species['city'].value_counts()
filtered_vc = vc[vc >= MIN_SAMPLES] 

## Pritty printing Just for printing purpose
temp_ms_df = copy.copy(metasub_species)
temp_ms_df['city'] = temp_ms_df['city'].apply(lambda x: x.replace('_',' ').title())
temp_vc = metasub_species['city'].value_counts()
temp_filtered_vc = temp_vc[temp_vc >= MIN_SAMPLES] 
print(temp_filtered_vc)
print('Total',temp_filtered_vc.sum())

## Iterate through each cities
for city in tqdm(filtered_vc.index):
    
    city_species_df = metasub_species.query('city == @city')
    print('For the city',city)
    print(city_species_df['core_project'].value_counts())
    col_not_to_include=["Samples","uuid",
                        "core_project",
                        "project",
                        "city",
                        "surface_material",
                        "continent",
                        "surface_ontology_fine",
                        "control_type",]
    _, _core_df, _sub_core_df, _peripheral_df = prevalance_extractor(city_species_df, bounds=bounds, col_not_to_include = col_not_to_include)
    core_species_dict[city] = { 
                                    'core':list(_core_df["Species"].values),
                                    'sub_core':list(_sub_core_df["Species"].values),
                                    'peripheral':list(_peripheral_df["Species"].values) 
                                  } 
                                  

with open(os.path.join(FINDINGS,'all_city_core.json'), 'w') as fp:
    json.dump(core_species_dict, fp)

#################################### S1 ######################################

"""
Which are the cities that have species at core level
"""

cities_with_core = []
for i in core_species_dict.keys():
    if len(core_species_dict[i]['core']) > 0:
        cities_with_core.append(i)

pd.DataFrame.from_dict({f'Cities with core microbes at rel abd {THRESHOLD}': cities_with_core}).to_csv(os.path.join(FINDINGS,'cities_with_core_microbes.csv'), index = False)

#!/usr/bin/env python
# coding: utf-8

# # Supplementary Notebook 3: Stratifying dataset into "low" and "high" abundant bacteria
# ## Paper: Novel Approach for Microbiome Analysis Using Bacterial Replication Rates and Causal Inference to Determine Resistome Potential
# ### Vitalii Stebliankin, Musfiqur Sazal, Camilo Valdes, Kalai Mathee, and GiriNarasimhan
# 
# #### Dataset: Gibson et al. (BioProject ID: PRJNA301903)
# 
# We define "high" abundant bacteria as 25% most abundant species, and "low" abundant as 25% least abundant taxa.

# In[1]:


import pandas #as pd

class DivideHighLowPlugin:
 def input(self, inputfile):
  self.PTR_file = inputfile#"analysis-out/1-FilteringPTR/PTR_species_filtered_metadata_major.csv"

 def run(self):
     pass

 def output(self, outputfile):
  out_dir = outputfile
  ptr_df = pandas.read_csv(self.PTR_file, index_col=0)
  def get_species_list(ptr_df):
    columns = ptr_df.columns
    species_list=[]
    for col in columns:
        if "PTR" in col:
            species = col.replace("#PTR", "")
            species_list.append(species)
    print("Total of {} species in the dataset.".format(len(species_list)))
    return species_list
  species_list = get_species_list(ptr_df)
  import seaborn as sns
  from matplotlib import pyplot as plt

  all_abundance = []
  for species in species_list:
    tmp_df =  ptr_df[ptr_df[species+"#PTR"].notnull()]
    #tmp_df = ptr_df
    tmp_df = tmp_df[tmp_df[species+"#abundance"]>0]
    all_abundance += list(tmp_df[species+"#abundance"])
  print("Total number of values: {}".format(len(all_abundance)))

  sns.distplot(all_abundance)


  all_abundance_sorted = sorted(all_abundance)
  low_cutoff = all_abundance_sorted[int(len(all_abundance)*0.25)]
  high_cutoff = all_abundance_sorted[int(len(all_abundance)*0.75)]

  # low_cutoff = 0.001
  # high_cutoff = 0.02

  print("Cutoff for low relative abundance: {}".format(low_cutoff))
  print("Cutoff for high relative abundance: {}".format(high_cutoff))
  print("")
  print("Distribution plot:")


  # In[2]:


  def get_species_list(ptr_df):
    columns = ptr_df.columns
    species_list=[]
    for col in columns:
        if "PTR" in col:
            species = col.replace("#PTR", "")
            species_list.append(species)
    print("Total of {} species in the dataset.".format(len(species_list)))
    return species_list
  species_list = get_species_list(ptr_df)
  import seaborn as sns
  from matplotlib import pyplot as plt

  all_abundance = []
  for species in species_list:
    tmp_df =  ptr_df[ptr_df[species+"#PTR"].notnull()]
    #tmp_df = ptr_df
    tmp_df = tmp_df[tmp_df[species+"#abundance"]>0]
    all_abundance += list(tmp_df[species+"#abundance"])
  print("Total number of values: {}".format(len(all_abundance)))

  sns.distplot(all_abundance)


  all_abundance_sorted = sorted(all_abundance)
  low_cutoff = all_abundance_sorted[int(len(all_abundance)*0.25)]
  high_cutoff = all_abundance_sorted[int(len(all_abundance)*0.75)]

  print("Cutoff for low relative abundance: {}".format(low_cutoff))
  print("Cutoff for high relative abundance: {}".format(high_cutoff))
  print("")
  print("Distribution plot:")

  from scipy.stats import ttest_ind
  import matplotlib.style as style
  import matplotlib
  import pandas as pd

  def high_low_prepare(PTR_file, out_all):
    ptr_df = pd.read_csv(PTR_file, index_col=0)
    def get_antibiotics(x):
        try:
            if "TC" in x:
                return "Ticarcillin-Clavulanate"
            elif "Mero" in x or "Amp" in x:
                return 'Ampicillin/Meropenem'
        except TypeError:
            return "Other"
        
    
    ptr_df["Antibiotics_group"] = ptr_df["Antibiotic_Treatment_unfiltered"].apply(lambda x: get_antibiotics(x))
    def get_species_list(ptr_df):
        columns = ptr_df.columns
        species_list=[]
        for col in columns:
            if "PTR" in col:
                species = col.replace("#PTR", "")
                species_list.append(species)
        print("Total of {} species in the dataset.".format(len(species_list)))
        return species_list
    species_list = get_species_list(ptr_df)
    def get_abundance_group(x):
        if x<low_cutoff:
            return "low"
        elif x>high_cutoff:
            return "high"
        else:
            return "medium"

    def get_before_after(x):
        x=str(x)
        if "efore" in x:
            return "Before"
        elif "fter" in x:
            return "After"
        else:
            return "Before"

    group_dict = {"PTR":[], "abundance":[], "Abundance Level":[], "Cohort":[], "Treatment":[], "TreatmentType":[], "species":[], "Antibiotics_group":[]}
    antibiotics = ["Gentamicin","Ampicillin","Meropenem","Vancomycin","Ticarcillin-Clavulanate",
                  "r_Gentamicin","r_Ampicillin","r_Meropenem","r_Vancomycin","r_Ticarcillin-Clavulanate"]
    for ant in antibiotics:
        group_dict[ant]=[]

    
    def get_groups(ptr_df):
        for species in species_list:
            tmp_df = ptr_df
            #tmp_df = ptr_df[ptr_df[species+".PTR"].notnull()]
            tmp_df["group"] = tmp_df[species+"#abundance"].apply(lambda x: get_abundance_group(x))
            tmp_df["Treatment"] = tmp_df["Antibiotic_Treatment_unfiltered"].apply(lambda x: get_before_after(x))
            
            group_dict["PTR"]+=list(tmp_df[species+"#PTR"])
            group_dict["abundance"]+=list(tmp_df[species+"#abundance"])
            group_dict["Abundance Level"] +=list(tmp_df["group"])
            group_dict["Cohort"]+= list(tmp_df["Cohort"])
            group_dict["Treatment"]+= list(tmp_df["Treatment"])
            group_dict["TreatmentType"]+=list(tmp_df["Antibiotic_Treatment_unfiltered"])
            group_dict["species"]+=[species for x in range(0,len(tmp_df))]
            group_dict["Antibiotics_group"]+=list(tmp_df["Antibiotics_group"])
            for ant in antibiotics:
                group_dict[ant]+=list(tmp_df[ant])

        group_df = pd.DataFrame(group_dict)
        return group_df

    # For all samples:
    group_df = get_groups(ptr_df)
    #group_df = group_df.fillna('noRecentTreatment')
    group_df.to_csv(out_all, index=False)
    return group_df

  out_all = out_dir+"/high-low.csv"#"./analysis-out/3-Divide_High-Low/high-low.csv"
  high_low_prepare(self.PTR_file, out_all)


  # In[3]:


  group_df = pd.read_csv(out_all)
  species_list = ["Klebsiella pneumoniae",
  "Escherichia coli",
  "Enterococcus faecalis",
  "Staphylococcus epidermidis",
  "Enterobacter cloacae",
  "Enterobacter hormaechei",
  "Klebsiella oxytoca",
  "Enterococcus faecium",
  "Bifidobacterium longum",
  "Klebsiella quasipneumoniae",
  "Klebsiella variicola",
  "Klebsiella michiganensis",
  "Klebsiella aerogenes",
  "Citrobacter freundii",
  "Veillonella parvula",
  "Clostridium perfringens",
  "Enterobacter roggenkampii",
  "Enterobacter asburiae",
  "Enterobacter sp. CRENT-193",
  "Klebsiella sp. M5al",
  "Enterobacter sp. DKU_NT_01",
  "Enterobacter sp. HK169",
  "Citrobacter"]
  group_df = group_df[group_df['species'].isin(species_list)]


  # In[4]:


  avg_df = group_df[['species', 'PTR', 'abundance']].groupby('species').mean()
  avg_df = avg_df.sort_values(by='abundance', ascending=False)


  # In[5]:


  avg_df


  # In[6]:


  import seaborn as sns
  sns.set_theme()
  sns.set(rc={'figure.figsize':(7,11.7)})

  g = sns.barplot(x = 'PTR', y = 'species', hue = 'Cohort', data = group_df,
            palette = 'hls', order=avg_df.index)
  g.set(xlim=(1, 3))
  g


  # In[7]:


  import seaborn as sns
  sns.set_theme()
  sns.set(rc={'figure.figsize':(7,11.7)})

  g = sns.barplot(x = 'abundance', y = 'species', hue = 'Cohort', data = group_df,
            palette = 'hls', order=avg_df.index)
  #g.set(xlim=(1, 3))
  g


  # 

  # In[8]:


  # Check signigicance Abundance
  from scipy import stats

  for sp in avg_df.index:
    tmp_df = group_df[group_df['species']==sp]
    control = tmp_df[tmp_df['Cohort']=='Control']['abundance']
    antibiotics = tmp_df[tmp_df['Cohort']=='Antibiotic']['abundance']
    
    U1, p = stats.mannwhitneyu(control, antibiotics)
    if p<0.05:
        print(sp, p)


  # In[9]:


  # Check signigicance PTR
  from scipy import stats

  for sp in avg_df.index:
    tmp_df = group_df[group_df['species']==sp]
    tmp_df = tmp_df[tmp_df['PTR'].notnull()]
    control = tmp_df[tmp_df['Cohort']=='Control']['PTR']
    antibiotics = tmp_df[tmp_df['Cohort']=='Antibiotic']['PTR']
    
    U1, p = stats.mannwhitneyu(control, antibiotics)
    if p<0.05:
        print(sp, p)
    


  # In[10]:


  # Check signigicance PTR-TIM
  from scipy import stats

  unique_significant_species = set()

  for sp in avg_df.index:
    tmp_df = group_df[group_df['species']==sp]
    tmp_df = tmp_df[tmp_df['PTR'].notnull()]
    control = tmp_df[tmp_df['Cohort']=='Control']['PTR']
    antibiotics = tmp_df[tmp_df['Antibiotics_group']=='Ticarcillin-Clavulanate']
    antibiotics = antibiotics[antibiotics['Treatment']=='After']['PTR']

    U1, p = stats.mannwhitneyu(control, antibiotics)
    if p<0.05 and p!=0 and sp!='Klebsiella variicola':
        print(sp, p)
        unique_significant_species.add(sp)
    


  # In[11]:


  # Check signigicance PTR-AMP
  from scipy import stats

  for sp in avg_df.index:
    tmp_df = group_df[group_df['species']==sp]
    tmp_df = tmp_df[tmp_df['PTR'].notnull()]

    control = tmp_df[tmp_df['Cohort']=='Control']['PTR']
    antibiotics = tmp_df[tmp_df['Antibiotics_group']=='Ampicillin/Meropenem']
    antibiotics = antibiotics[antibiotics['Treatment']=='After']['PTR']

    U1, p = stats.mannwhitneyu(control, antibiotics)
    if p<0.05 and p!=0:
        unique_significant_species.add(sp)
        print(sp, p)
    


  # In[12]:


  # Check signigicance Abundance-TIM
  from scipy import stats

  for sp in avg_df.index:
    tmp_df = group_df[group_df['species']==sp]
    tmp_df = tmp_df[tmp_df['PTR'].notnull()]
    control = tmp_df[tmp_df['Cohort']=='Control']['abundance']
    antibiotics = tmp_df[tmp_df['Antibiotics_group']=='Ticarcillin-Clavulanate']
    antibiotics = antibiotics[antibiotics['Treatment']=='After']['abundance']

    U1, p = stats.mannwhitneyu(control, antibiotics)
    if p<0.05 and p!=0:
        unique_significant_species.add(sp)
        print(sp, p)
    


  # In[13]:


  # Check signigicance Abundance-AMP
  from scipy import stats

  for sp in avg_df.index:
    tmp_df = group_df[group_df['species']==sp]
    tmp_df = tmp_df[tmp_df['PTR'].notnull()]
    control = tmp_df[tmp_df['Cohort']=='Control']['abundance']
    antibiotics = tmp_df[tmp_df['Antibiotics_group']=='Ampicillin/Meropenem']
    antibiotics = antibiotics[antibiotics['Treatment']=='After']['abundance']

    U1, p = stats.mannwhitneyu(control, antibiotics)
    if p<0.05 and p!=0:
        unique_significant_species.add(sp)
        print(sp, p)


  # In[14]:


  print(unique_significant_species)
  print(len(unique_significant_species))


  # In[15]:


  import seaborn as sns
  from matplotlib import pyplot as plt

  group_df_tmp = group_df
  species_list = ['Klebsiella pneumoniae','Enterobacter cloacae', 'Enterobacter hormaechei', 'Klebsiella oxytoca','Klebsiella quasipneumoniae', 'Klebsiella aerogenes', 'Citrobacter freundii']

  group_df_tmp['Antibiotics_group'] = group_df_tmp.apply(lambda row: 'control' if row['Cohort']=='Control' else row['Antibiotics_group'], axis=1)
  group_df_tmp = group_df_tmp[group_df_tmp['Antibiotics_group']!='Other']
  group_df_tmp = group_df_tmp[group_df_tmp['species'].isin(unique_significant_species)]

  sns.set(rc={'figure.figsize':(3, 8)})
  sns.set_theme(style="whitegrid")

  g = sns.boxplot(x = 'PTR', y = 'species', hue = 'Antibiotics_group', data = group_df_tmp,
            palette = 'hls', order=species_list)
  g.set(xlim=(1, 3.5))
  plt.legend([],[], frameon=False)

  plt.savefig(out_dir+'/ptr_amp_tim.png', dpi=1000)

  g


  # In[16]:


  g = sns.boxplot(x = 'abundance', y = 'species', hue = 'Antibiotics_group', data = group_df_tmp,
            palette = 'hls', order=species_list)
  g.set(xlim=(0, 0.5))
  plt.legend([],[], frameon=False)

  plt.savefig(out_dir+'/abundance_amp_tim.png', dpi=1000)

  g


  # In[ ]:





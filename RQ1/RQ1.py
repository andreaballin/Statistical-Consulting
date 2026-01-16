import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# DATA 
df_2024 = pd.read_excel('bioblitz_export-WUR2024.xlsx')
df_2025 = pd.read_excel('bioblitz_export-WUR2025.xlsx')
df_2024['year'] = 2024
df_2025['year'] = 2025
plants_2024 = df_2024[(df_2024['species_group_name'] == 'Plants')|(df_2024['species_group_name'].str.lower() == 'mosses and lichens')]
plants_2025 = df_2025[(df_2025['species_group_name'] == 'Plants')|(df_2025['species_group_name'].str.lower() == 'mosses and lichens')]
richness_2024 = plants_2024['species_name_scientific'].nunique()
richness_2025 = plants_2025['species_name_scientific'].nunique()
obs_freq_2024 = plants_2024['species_name_scientific'].value_counts()
obs_freq_2025 = plants_2025['species_name_scientific'].value_counts()
print(f"  2024 most frequent species: {obs_freq_2024.index[0]}")
print(f"  2025 Most frequent species: {obs_freq_2025.index[0]}")


# Jaccard Index
un_spec_2024 = set(plants_2024['species_name_scientific'].unique())
un_spec_2025 = set(plants_2025['species_name_scientific'].unique())


shared_species = un_spec_2024.intersection(un_spec_2025)
lost_species = un_spec_2024.difference(un_spec_2025)
new_species = un_spec_2025.difference(un_spec_2024)

jaccard_similarity = len(shared_species) / len(un_spec_2024 | un_spec_2025)

print("\nJaccard Index:")
print(f"Jaccard Similarity: {jaccard_similarity:.4f}")

# top species and top species change

total_observations_df = pd.DataFrame({'2024': obs_freq_2024,'2025': obs_freq_2025}).fillna(0)
total_observations_df['Change over Time'] = total_observations_df['2025'] - total_observations_df['2024']

print("Year 2025 Top 10 Most Observed:")
print(obs_freq_2024.head(10))
print("Year 2025 Top 10 Most Observed:")
print(obs_freq_2025.head(10))


top_increases = total_observations_df.nlargest(10, 'Change over Time')
# print(top_increases[['2024', '2025', 'Change over Time']])

top_decreases = total_observations_df.nsmallest(10, 'Change over Time')
# print(top_decreases[['2024', '2025', 'Change over Time']])

# Rarefaction Analysis

def rarefaction(data, sample_size):
    results = []

    for each in sample_size:
        equal_efforts = []

        for i in range(100):
            size = min(each, len(data))
            sample = data.sample(n=size, replace=False)
            equal_efforts.append(sample['species_name_scientific'].nunique())

        results.append({
            'Sample Size': each,
            'Mean Richness': np.mean(equal_efforts),
            'Std Richness': np.std(equal_efforts),
            'CI Low': np.percentile(equal_efforts, 2.5),
            'CI High': np.percentile(equal_efforts, 97.5)
        })

    return pd.DataFrame(results)

#Minimum sample size
sample_sizes = list(range(50, min(len(plants_2024), len(plants_2025)) +1, 50))


es_2024 = rarefaction(plants_2024, sample_sizes)
es_2025 = rarefaction(plants_2025, sample_sizes)

es_2024_eq = es_2024[es_2024['Sample Size'] <= min(len(plants_2024), len(plants_2025))].iloc[-1]
es_2025_eq = es_2025[es_2025['Sample Size'] <= min(len(plants_2024), len(plants_2025))].iloc[-1]

#Permutation Test
def permutation_test(df1, df2, n_permutations=10000):
 
    observed_diff = len(np.unique(df2)) - len(np.unique(df1))

    # Combine both groups
    combined = np.concatenate((df1, df2))
    n1 = len(df1)

    perm_diffs = []
    for _ in range(n_permutations):
        np.random.shuffle(combined)
        perm1 = combined[:n1]
        perm2 = combined[n1:]

        diff = len(np.unique(perm2)) - len(np.unique(perm1))
        perm_diffs.append(diff)

    perm_diffs = np.array(perm_diffs)
    p_value = np.sum(np.abs(perm_diffs) >= np.abs(observed_diff)) / n_permutations

    return observed_diff, p_value

richness_diff, p_value= permutation_test(plants_2024['species_name_scientific'].to_numpy(), plants_2025['species_name_scientific'].to_numpy(), n_permutations=10000)

print(f"p-value: {p_value:.4f}")
if p_value < 0.05:
    print("Result: Significant (p < 0.05), we can reject null hypothesis")
else:
    print("Result: Not Significant (p ≥ 0.05), we can't reject null hypothesis")


# visualz

#Fig1: Showing species Richness
fig, ax = plt.subplots()
ax.set_title('Wageningen University Plant Species Richness')
ax.set_ylabel('Number of Unique Species')


labels = ax.bar(['2024', '2025'], [richness_2024, richness_2025])
for label in labels:
    height = label.get_height()
    ax.text(label.get_x() + label.get_width()/2., height, f'{int(height)}',ha='center', va='bottom')

#Fig2: Species Turnover
fig, ax = plt.subplots()
ax.set_ylabel('Number of Species')
ax.set_title('Species Turnover 2024-2025')
categories = ['Shared Species', 'Species present in 2024 only', 'Species New in 2025']
values = [len(shared_species), len(lost_species), len(new_species)]


labels = ax.bar(categories, values)
for label, val in zip(labels, values):
    height = label.get_height()
    pct = val / len(un_spec_2024 | un_spec_2025) * 100
    ax.text(label.get_x() + label.get_width()/2., height, f'{int(val)}\n({pct:.1f}%)',
            ha='center', va='bottom')


#Fig3-a: Number of observations per species 2024
fig, ax = plt.subplots()
bins_2024 = [1, 2, 5, 10, 20, 50, max( obs_freq_2024.max(), 50) + 1]
labels_2024 = ['1', '2-4', '5-9', '10-19', '20-49', '50+']
plt.xticks(range(len(labels_2024)), labels_2024)
plt.xlabel('Observations per Species')
plt.ylabel('Number of Species')
plt.title('2024 Observation Frequency')

cat_2024 = pd.cut(obs_freq_2024, bins=bins_2024, labels=labels_2024, right=False)
countss= cat_2024.value_counts().sort_index()

labels = plt.bar(range(len(countss)), countss.values)
for i, label in enumerate(labels):
    plt.text(label.get_x() + label.get_width()/2, label.get_height(),str(int(label.get_height())), ha='center', va='bottom')

#Fig3-b: Number of observations per species 2025
fig, ax = plt.subplots()
bins_2025 = [1, 2, 5, 10, 20, 50, max( obs_freq_2025.max(), 50) + 1]
labels_2025 = ['1', '2-4', '5-9', '10-19', '20-49', '50+']
plt.xticks(range(len(labels_2025)), labels_2025)
plt.xlabel('Observations per Species')
plt.ylabel('Number of Species')
plt.title('2025 Observation Frequency')

cat_2025 = pd.cut(obs_freq_2025, bins=bins_2025, labels=labels_2025, right=False)
countsss= cat_2025.value_counts().sort_index()

labels = plt.bar(range(len(countsss)), countsss.values)
for i, label in enumerate(labels):
    plt.text(label.get_x() + label.get_width()/2, label.get_height(),str(int(label.get_height())), ha='center', va='bottom')


#Fig4: Rarefaction Analysis
fig, ax = plt.subplots()
ax.set_title('Species Accumulation Curves with CI')
ax.set_xlabel('Number of Observation Records')
ax.set_ylabel('Expected Species Richness')


ax.plot(es_2024['Sample Size'], es_2024['Mean Richness'],label='2024')
ax.fill_between(es_2024['Sample Size'], es_2024['CI Low'], es_2024['CI High'],alpha=0.4)

ax.plot(es_2025['Sample Size'], es_2025['Mean Richness'],label='2025')
ax.fill_between(es_2025['Sample Size'], es_2025['CI Low'],es_2025['CI High'], alpha=0.4)

ax.legend()


#Fi5: Top Observation frequency changes in species

# Increases
plt.figure()
plt.xlabel('Change in Observation Frequency')
plt.title('Top 10 Increases in Species Observations')

top_inc = np.arange(len(top_increases))
label_inc = plt.barh(top_inc, top_increases['Change over Time'])
plt.yticks(top_inc, [name for name in top_increases.index])


for label, value in zip(label_inc, top_increases['Change over Time']):
    plt.text(value, label.get_y() + label.get_height() / 2,f'+{int(value)}', va='center')

# Decreases
plt.figure()
plt.xlabel('Change in Observation Frequency')
plt.title('Top 10 Decreases in Species Observations')

top_dec = np.arange(len(top_decreases))
label_dec = plt.barh(top_dec, top_decreases['Change over Time'])
plt.yticks(top_dec, [name for name in top_decreases.index])


for label, value in zip(label_dec, top_decreases['Change over Time']):
    plt.text(value, label.get_y() + label.get_height() / 2,f'{int(value)}', va='center')







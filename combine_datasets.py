'''
datasetcombiner.py
purpose: combine data from different CDR3 datasets into one simple table.
author: Yuta Nagano
ver: 1.0.1
'''


from itertools import combinations
import pandas as pd
from polyleven import levenshtein as leven


# Main code
if __name__ == '__main__':
    # Load the datasets
    vdjdb = pd.read_csv('datasets/vdjdb.tsv',sep=r'\t')
    iedb = pd.read_csv('datasets/iedb.csv')


    # Clean the iedb dataset
    # 1) Remove all rows with odd/modified antigen epitopes
    iedb = iedb[~iedb['Description'].str.contains(r'\+')]

    # 2) Get a unified column for beta chain CDR3 data prioritising any
    # 'curated' sequences, and also add in any missing cysteine and
    # phenylalanine residues.
    cdr3s = iedb['Chain 2 CDR3 Curated'].where(
        ~iedb['Chain 2 CDR3 Curated'].isnull(),
        iedb['Chain 2 CDR3 Calculated']
    )

    def sandwich(cdr3: str) -> str:
        if not cdr3.startswith('C'): cdr3 = 'C'+cdr3+'F'
        return cdr3

    iedb['CDR3'] = cdr3s.map(sandwich)

    # 3) Drop any rows where either there is no CDR3 information or epitope
    # information detected
    iedb = iedb.dropna(subset=['CDR3', 'Description'])


    # Combine the two tables (keep columns for the CDR3, epitope, and add a new
    # column to keep track of which database a particular entry came from)
    vdjdb = vdjdb[['Epitope','CDR3']]
    vdjdb['Dataset'] = 'vdjdb'

    iedb = iedb[['Description','CDR3']]
    iedb = iedb.rename(columns={'Description': 'Epitope'})
    iedb['Dataset'] = 'iedb'

    combined = pd.concat([vdjdb, iedb])


    # Filter out duplicate rows
    combined = combined.drop_duplicates(subset=['CDR3','Epitope'])


    # Filter out promiscuous T cells that appear to be responding to very
    # dissimilar epitopes
    grouped_cdr3 = combined.groupby('CDR3')
    num_epitopes = grouped_cdr3.count()['Epitope']
    promiscuous = num_epitopes[num_epitopes > 1].index.tolist()
    true_promiscuous = []

    for cdr3 in promiscuous:
        epitopes = grouped_cdr3.get_group(cdr3)['Epitope']
        for e1, e2 in combinations(epitopes, 2):
            if leven(e1,e2,1) > 1:
                true_promiscuous.append(cdr3)
                break

    combined = combined[~combined['CDR3'].isin(true_promiscuous)]


    # Remove any epitope groups with size less than 10
    num_cdr3s = combined.groupby('Epitope').count()['CDR3']
    tiny = num_cdr3s[num_cdr3s < 10].index.tolist()

    combined = combined[~combined['Epitope'].isin(tiny)]


    # Sort by epitope and re-index the table
    combined = combined.sort_values(by=['Epitope','CDR3'])
    combined = combined.reset_index(drop=True)


    # Save the cleaned, combined table
    combined.to_csv('datasets/vdjdb_iedb_combined.csv',index=False)

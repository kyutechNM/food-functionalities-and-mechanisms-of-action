import pandas as pd
from gmpy2 import comb
from itertools import chain
import argparse
import numpy as np
# from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection
from functools import lru_cache
from itertools import combinations_with_replacement
from itertools import chain


# import piper

class Pathway():
    def __init__(self,
            pathway_file='your_file1',
            excluded_pathway_file='excluded_pathway_file',
            FDR_each_compound=True,
            FDR_all=False,
            detail=False,
            pathway_id_col_name="pathway_columname",
            protein_id_col_name="protein_columname",
            ):

        self.FDR_each_compound = FDR_each_compound
        self.FDR_all = FDR_all
        self.detail = detail

        # Reading the data required for enrichment analysis.
        self.pathway_use = pd.read_table(pathway_file)

        if excluded_pathway_file:
            self.excluded_pathway = set(pd.read_table(excluded_pathway_file, header=None)[0])
            self.pathway_use = self.pathway_use[~self.pathway_use[pathway_id_col_name].isin(self.excluded_pathway)]

        self.pathway2protein = self.pathway_use.groupby(pathway_id_col_name)[protein_id_col_name].agg(list).to_dict()
        self.all_genes = set(self.pathway_use[protein_id_col_name])

    @lru_cache(maxsize=None)
    def p_value(self, l, r, k, z):
# Functions for which values are obtained using Fisher's exact probability test.
# l: Number of all genes of interest
# r: Number of genes altered by addition of compounds etc.
# k: Number of genes included in the pathway
# z: Number of genes in the product set of all target genes and the set of genes that have been changed by addition of compounds etc.
        p = sum(map(lambda i: comb(k, i)*comb(l-k, r-i), range(z, min(r, k)+1)))/comb(l, r)
        return np.float64(p)

#     def p_value_scipy(self, l, r, k, z):
#         return fisher_exact([[z, r-z],[k-z, l-k+z-r]])[1]
    
    # pathwayIDに対してp値を求める関数
    def p_value_pathway(self, pathway, selected_genes):
        l = len(self.all_genes)
        r = len(selected_genes&self.all_genes)
        k = len(set(self.pathway2protein[pathway])&self.all_genes)
        z = len(selected_genes&set(self.pathway2protein[pathway])&self.all_genes)
        return {
            "ID":pathway,
            "p-value":self.p_value(l, r, k, z),
            }

    def odds_ratio(self, A, B, C, D):
        try:
            return (A*D)/(B*C)
        except:
            return np.nan

    def p_value_pathway_detail(self, pathway, selected_genes):
        l = len(self.all_genes)
        r = len(selected_genes&self.all_genes)
        k = len(set(self.pathway2protein[pathway])&self.all_genes)
        z = len(selected_genes&set(self.pathway2protein[pathway])&self.all_genes)
        A = z
        B = r-z
        C = k-z
        D = l-k+z-r
        return {
            "ID":pathway,
            "p-value":self.p_value(l, r, k, z),
            "A(number_of_selected_genes_in_payhway)":A,
            "B(number_of_selected_genes_not_in_pathway)":B,
            "C(number_of_non-selected_genes_in_pathway)":C,
            "D(number_of_non-selected_genes_not_in_pathway)": D,
            "A+C": k,
            "A/(A+C)(filling_rate)": z/k,
            "(A*D)/(B*C)(odds_ratio)":self.odds_ratio(A, B, C, D),
            }

# Function to do FDR correction (np.array -> np.array)
    def FDR_statsmodels(self, p):
        return fdrcorrection(p, is_sorted=False)[1]

# ====================================================================
# Function to calculate FDR-corrected p-values in pathway enrichment analysis.
# (set -> pd.DataFrame)
    def enrichiment_analysis(self, selected_genes):

        p_value_pathway = self.p_value_pathway_detail if self.detail else self.p_value_pathway

        # Calculate p-values for all pathways
        df = pd.DataFrame(map(lambda x: p_value_pathway(x, selected_genes), self.pathway2protein))

        # FDR correction
        if self.FDR_each_compound:
            df.insert(1,'p-value (FDR-corrected in each_compound)', self.FDR_statsmodels(df['p-value'].values))
        return df

    def enrichiment_analysis_with_name(self, selected_genes, name, col_name):
        df = self.enrichiment_analysis(selected_genes)
        df.insert(0, col_name, name)
        return df


    def enrichiment_analysis_all(self, selected_genes_dict, col_name):
        df = pd.concat(map(lambda x: self.enrichiment_analysis_with_name(selected_genes_dict[x], x, col_name), selected_genes_dict))
        if self.FDR_all:
            df.insert(2,'p-value (FDR-corrected in all p-values)', self.FDR_statsmodels(df['p-value'].values))

        return df  

    def enrichiment_analysis_all_comb(self, selected_genes_dict,  col_name):
        comb_key = combinations_with_replacement(selected_genes_dict, 2)
        df = pd.concat(map(lambda x: self.enrichiment_analysis_with_name(
            selected_genes_dict[x[0]] | selected_genes_dict[x[1]], x[0]+ '&' + x[1], col_name), comb_key))
        name_df = df[col_name].str.split("&", expand=True).rename(columns={0:col_name + "1", 1:col_name + "2"})
        return pd.concat([name_df, df.drop(columns=[col_name])], axis=1)

    def enrichiment_analysis_all_df(self, df, threshold, compound="compound", protein="protein", score="score"):
        selected_genes_dict = dict.fromkeys(df[compound].drop_duplicates(), set())
        selected_genes_dict.update(df[df[score] >= threshold].groupby(compound)[protein].agg(set).to_dict())
        return self.enrichiment_analysis_all(selected_genes_dict, compound)
    
    def enrichiment_analysis_all_tsv(self, tsv_file_name, threshold, compound="compound", protein="protein", score="score", sep="\t"):
        df = pd.read_csv(tsv_file_name, sep=sep)
        return self.enrichiment_analysis_all_df(df, threshold, compound=compound, protein=protein, score=score)

    def __call__(self, selected_genes_dict, col_name="compound"):
        return self.enrichiment_analysis_all(selected_genes_dict, col_name)

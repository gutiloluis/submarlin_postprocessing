#%%
import goatools.base
from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
from goatools.anno.gaf_reader import GafReader
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.semantic import semantic_similarity
from goatools.semantic import TermCounts, get_info_content
import submarlin_postprocessing.clustering_viz as clustering_viz
import submarlin_postprocessing.filepaths as filepaths
import numpy as np
import pandas as pd
#%%
def search_go(ns2assoc, obodag, inv_gene_to_id, go_term):
    namespace_abbv = {"biological_process":"BP","molecular_function":"MF","cellular_component":"CC"}
    
    print("Searching for " + str(obodag[go_term].name))
    namespace = namespace_abbv[obodag[go_term].namespace]
    child_goterms = list(obodag[go_term].get_all_children())
    gene_list = [inv_gene_to_id[key] for key,val in ns2assoc[namespace].items() if go_term in val]
    for child_goterm in child_goterms:    
        gene_list += [inv_gene_to_id[key] for key,val in ns2assoc[namespace].items() if child_goterm in val]
    gene_list = sorted(list(set(gene_list)))
    return gene_list

def highlight_go(an_df, ns2assoc, obodag, inv_gene_to_id, go_term, go_mode = "all", color="red", fontsize=18):
    
    gene_list = search_go(ns2assoc, obodag, inv_gene_to_id, go_term)
    if go_mode == "all":
        fig = highlight_regex(an_df,gene_list,column_label="Gene",color_all_mode=color,highlight="all",match_mode="exact",fontsize=fontsize)
    elif go_mode == "each":
        fig = highlight_regex(an_df,gene_list,column_label="Gene",highlight="each",match_mode="exact",fontsize=fontsize)
        
    plt.title(go_term)
    return fig

def get_enriched_GO_terms(background_gene_list,gene_list,obodag,objanno,ns2assoc,pval=0.05,GO_type="BP"):
    
    gene_to_id = {assoc.DB_Symbol:assoc.DB_ID for assoc in objanno.associations}
    synonym_dict = {synonym:assoc.DB_ID for assoc in objanno.associations for synonym in assoc.DB_Synonym}
    gene_to_id.update(synonym_dict)

    #background gene set

    all_genes_uniprot = [gene_to_id[item] for item in background_gene_list if item in gene_to_id.keys()]
    selected_genes_uniprot = [gene_to_id[item] for item in gene_list if item in gene_to_id.keys()]
    
    print(len(all_genes_uniprot))
    print(len(selected_genes_uniprot))
    
    goeaobj = GOEnrichmentStudy(
            all_genes_uniprot, # List of mouse protein-coding genes
            ns2assoc[GO_type], # geneid/GO associations
            obodag, # Ontologies
            propagate_counts = True,
            alpha = pval, # default significance cut-off
            methods=["fdr_bh"]); # defult multipletest correction method

    goea_results_all = goeaobj.run_study(selected_genes_uniprot, prt=None)
    goea_quiet_sig = [r for r in goea_results_all if r.p_fdr_bh < pval]
    goea_quiet_enriched = [r for r in goea_quiet_sig if r.enrichment == "e"]
    return goea_quiet_enriched

def pick_exemplar(go1,go2,termcounts,obodag,info_thr,pval_factor=2.):
    
    info_1_low = get_info_content(go1.GO, termcounts)<info_thr
    info_2_low = get_info_content(go2.GO, termcounts)<info_thr
    if info_1_low and not info_2_low:
        return go2
    elif info_2_low and not info_1_low:
        return go1
    elif info_2_low and info_1_low:
        return go1
    
    pval_ratio = go1.p_fdr_bh/go2.p_fdr_bh
    
    if pval_ratio > pval_factor:
        return go2
    elif pval_ratio < (1./pval_factor):
        return go1
    
    go1_parents = list(obodag[go1.GO].get_all_parents())
    go2_parents = list(obodag[go2.GO].get_all_parents())
    
    if go2.GO in go1_parents:
        return go2
    
    elif go1.GO in go2_parents:
        return go1
    
    return go1

def get_filtered_go_terms(obodag,objanno,goea_list,sim_thr = 0.05,info_thr = 1.,GO_type="BP"):
    
    termcounts = TermCounts(obodag, objanno.get_ns2assc()[GO_type])

    go_term_list = [item.GO for item in goea_list]
    sim_arr = np.zeros((len(go_term_list),len(go_term_list)))
    for i in range(len(go_term_list)):
        for j in range(len(go_term_list)):
            sim_arr[i,j] = semantic_similarity(go_term_list[i], go_term_list[j], obodag)
    np.fill_diagonal(sim_arr, 0.)
    
    working_group_idx = 0
    grouped_terms = {}
    group_exemplars = {}
    go_term_indices = list(range(len(go_term_list)))

    while len(go_term_indices)>0:
        i = go_term_indices[0]
        most_sim_arg = np.argmax(sim_arr[i])
        sim_score = sim_arr[i,most_sim_arg]
        if sim_score > sim_thr:
            if len(grouped_terms)>0:
                in_other_group_keys = [key for key,val in grouped_terms.items() if most_sim_arg in val]
                if len(in_other_group_keys)==1:
                    other_group_idx = in_other_group_keys[0]
                    grouped_terms[other_group_idx] = grouped_terms[other_group_idx] + [i]
                    group_exemplars[other_group_idx] = pick_exemplar(group_exemplars[other_group_idx],goea_list[i],termcounts,obodag,info_thr)
                else:        
                    grouped_terms[working_group_idx] = [i,most_sim_arg]
                    group_exemplars[working_group_idx] = pick_exemplar(goea_list[i],goea_list[most_sim_arg],termcounts,obodag,info_thr)
                    working_group_idx += 1
                    go_term_indices.remove(most_sim_arg)
            else:
                grouped_terms[working_group_idx] = [i,most_sim_arg]
                group_exemplars[working_group_idx] = pick_exemplar(goea_list[i],goea_list[most_sim_arg],termcounts,obodag,info_thr)
                working_group_idx += 1
                go_term_indices.remove(most_sim_arg)
        go_term_indices.remove(i)
    
    group_exemplars = list(group_exemplars.values())
    
    return group_exemplars

def get_GO_assign_dict(selected_goea,cluster_genes_uniprot):
    all_study_items = copy.copy(cluster_genes_uniprot)
    depth_list = sorted(set([item.depth for item in selected_goea]))[::-1]
    assign_dict = {}
    for depth in depth_list:
        go_terms_at_level = [item for item in selected_goea if item.depth == depth]
        for go_term in go_terms_at_level:
            study_item_list = list(go_term.study_items)
            for study_item in study_item_list:
                if study_item in all_study_items:
                    assign_dict[study_item] = go_term.name
                    all_study_items.remove(study_item)

    for remaining_item in all_study_items:
        assign_dict[remaining_item] = "Unassigned"
        
    return assign_dict

def _safe_div(a, b):
    return float(a) / float(b) if b else 0.0

COLUMN_SPECS = {
    "GO":                   lambda r: r.GO,
    "GO Term":              lambda r: r.goterm.name,
    "N group in":           lambda r: r.ratio_in_study[0],
    "N background in":      lambda r: r.ratio_in_pop[0],
    "N group total":        lambda r: r.ratio_in_study[1],
    "N background total":   lambda r: r.ratio_in_pop[1],
    "Percent of GO in group": lambda r: _safe_div(r.ratio_in_study[0], r.ratio_in_pop[0]),
    "Percent of group in GO": lambda r: _safe_div(r.ratio_in_study[0], r.ratio_in_study[1]),
    "FDR":                  lambda r: r.p_fdr_bh,
}

def get_go_enrichment_df(goea_results, column_specs=None):
    """
    Build a GO enrichment summary DataFrame using only your chosen columns.
    """
    specs = column_specs or COLUMN_SPECS

    rows = []
    for r in goea_results:
        rows.append({col: func(r) for col, func in specs.items()})

    df = pd.DataFrame(rows)
    return df

def find_genes_with_go_term_in_cluster(
    genes_in_go_term,
    genes_in_cluster
):
    '''
    Find which genes associated with a GO term are present in a given cluster.
    Returns a pandas DataFrame with a boolean column 'in_cluster' indicating presence.
    '''
    return (
        pd.DataFrame({
            'Gene': genes_in_go_term
        })
        .assign(in_cluster = lambda df_: df_['Gene'].isin(genes_in_cluster))
    )

def _concatenate_go_dfs(dfs:dict) -> pd.DataFrame:
        '''
        Concatenate GO enrichment DataFrames from multiple clusters into a single DataFrame.
        '''
        return (
            pd.concat(dfs.values(), keys=dfs.keys(), names=['Cluster ID', 'Row ID'])
            .reset_index()
            .assign(neglog10FDR = lambda df_: -np.log10(df_['FDR']))
        )

def get_all_genes_in_genome() -> list:
    '''
    Get the list of all genes in the genome from the sgrna timeseries file.
    Currently using sgrna_timeseries_filename from filepaths as reference table.
    '''
    sgrna_timeseries_filename = filepaths.sgRNA_timeseries_filenames['merged_all']
    return (
        pd.read_pickle(sgrna_timeseries_filename)
        ['Gene']
        .unique()
        .tolist()
    )

def get_all_genes_in_clustering_df(clustering_df:pd.DataFrame) -> list:
    '''
    Get the list of all genes in the clustering DataFrame.
    '''
    COLS_TO_KEEP = ['Gene', 'locus_tag', 'L0', 'L1', 'L2', 'L3']
    return(
        clustering_df
        .loc[:, COLS_TO_KEEP]
        .groupby('Gene')
        .first()
        .index
        .to_list()
    )

class GOEnrichmentAnalysis:
    def __init__(self):
        obo_fname = download_go_basic_obo()
        self.obodag = GODag(obo_fname)
        gaf_url = 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/6.B_subtilis_168.goa'
        gaf_fname = './6.B_subtilis_168.goa'

        goatools.base.http_get(gaf_url, gaf_fname)
        self.objanno = GafReader(gaf_fname)
        self.ns2assoc = self.objanno.get_ns2assc()

        self.gene_to_id = {assoc.DB_Symbol:assoc.DB_ID for assoc in self.objanno.associations}
        self.inv_gene_to_id = {assoc.DB_ID:assoc.DB_Symbol for assoc in self.objanno.associations}
        synonym_dict = {synonym:assoc.DB_ID for assoc in self.objanno.associations for synonym in assoc.DB_Synonym}
        self.gene_to_id.update(synonym_dict)
    
        self.namespace_abbv = {"biological_process":"BP","molecular_function":"MF","cellular_component":"CC"}

    def get_enriched_GO_terms(self, background_gene_list, gene_list, pval=0.05, GO_type="BP"):
        return get_enriched_GO_terms(background_gene_list, gene_list, self.obodag, self.objanno, self.ns2assoc, pval, GO_type)

    def get_filtered_go_terms(self, goea_list, sim_thr = 0.05, info_thr = 1., GO_type="BP"):
        return get_filtered_go_terms(self.obodag, self.objanno, goea_list, sim_thr, info_thr, GO_type)

    def search_go(self, go_term):
        return search_go(self.ns2assoc, self.obodag, self.inv_gene_to_id, go_term)

    def get_gene_list_for_cluster(self, clustering_df, cluster_id) -> list:
        '''
        Get the list of genes in a given cluster from the clustering DataFrame.
        '''
        COLS_TO_KEEP = ['Gene', 'locus_tag', 'L0', 'L1', 'L2', 'L3']
        subset = (
            clustering_df
            .loc[lambda df_: df_['L3'] == cluster_id, COLS_TO_KEEP]
            .groupby('Gene')
            .first()
            .index
            .to_list()
        )
        return subset

    def run_go_enrichment_analysis_and_filtering_single_group(
        self,
        background_gene_list,
        gene_list,
        pval=0.05,
        GO_type="BP"
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        '''
        Run GO enrichment analysis and filtering for a single gene list.
        '''
        goea_quiet_enriched = self.get_enriched_GO_terms(
            background_gene_list = background_gene_list,
            gene_list = gene_list,
            pval = pval,
            GO_type = GO_type
        )
        quiet_df = get_go_enrichment_df(goea_quiet_enriched)

        group_exemplars = self.get_filtered_go_terms(
            goea_quiet_enriched,
            sim_thr = 0.05,
            info_thr = 1.,
            GO_type = GO_type
        )
        exemplar_df = get_go_enrichment_df(group_exemplars)

        return quiet_df, exemplar_df

    def run_go_enrichment_analysis_and_filtering_single_cluster(
        self,
        clustering_df,
        background_gene_list,
        cluster_id,
        pval=0.05,
        GO_type="BP"
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        '''
        Run GO enrichment analysis and filtering for a single cluster.
        '''
        subset = self.get_gene_list_for_cluster(clustering_df, cluster_id)

        return self.run_go_enrichment_analysis_and_filtering_single_group(
            background_gene_list=background_gene_list,
            gene_list=subset,
            pval=pval,
            GO_type=GO_type
        )

    def run_go_enrichment_analysis_and_filtering_multiple_clusters(
        self,
        clustering_df,
        background_gene_list,
        clusters_to_include,
        pval=0.05,
        GO_type="BP"
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        '''
        Run GO enrichment analysis and filtering over multiple clusters.
        '''
        dfs_go = {}
        dfs_go_exemplar = {}

        for cluster_id in clusters_to_include:
            df_go, df_go_exemplar = self.run_go_enrichment_analysis_and_filtering_single_cluster(
                clustering_df=clustering_df,
                background_gene_list=background_gene_list,
                cluster_id=cluster_id,
                pval=pval,
                GO_type=GO_type
            )
            dfs_go[cluster_id] = df_go
            dfs_go_exemplar[cluster_id] = df_go_exemplar

        df_go = _concatenate_go_dfs(dfs_go)
        df_go_exemplar = _concatenate_go_dfs(dfs_go_exemplar)

        return df_go, df_go_exemplar

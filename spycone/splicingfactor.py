from collections import defaultdict
import numpy as np
from Bio import motifs
from Bio.Seq import Seq
import pandas as pd
import os
import pickle
import requests
import seaborn as sns
import matplotlib.pyplot as plt
from functools import reduce
from scipy.stats import pearsonr
from scipy.stats import mannwhitneyu, fisher_exact, kruskal
from Bio.SeqUtils import GC
from joblib import Parallel, delayed
import gc


from ._util_stat.multipletesting import pvalue_correction

dir_path = os.path.dirname(os.path.realpath(__file__))

class SF_coexpression():
    """
    Return the coexpression between Splicing factor and isoforms
    
    Parameters 
    ------------
    dataset : input dataset object

    padj_method : Multiple testing method. Default=Bonferonni

    corr_cutoff :  Correlation coefficient cutoff. Dafault=0.7

    padj_cutoff : Adjusted p-value cutoff. Default=0.05

    method : Method to calculate the correlation value.
    
    Return
    -------
    Create an instance for co-expression analysis of splicing factors and transcript abundance.

    """
    def __init__(self, dataset, padj_method="bonf", corr_cutoff=0.7, padj_cutoff = 0.05, method="pearson"):
        self.dataset=dataset
        self.padj_method = padj_method
        self.corr_cutoff = corr_cutoff
        self.method= method
        

    def _get_splicing_factors(self):
        sfs=[]
        with open(os.path.join(dir_path, "data/network/splicingfactors.txt"), "r") as f:
            tmp=f.readlines()
            for ss in tmp:
                sfs.append(ss.replace("\n", ""))

        sfs = pd.DataFrame(sfs, columns=['GeneSymbol'])
        mapper = dict(zip(self.dataset.symbs, self.dataset.gene_id))
        sfs['entrezid'] = list(map(lambda x: mapper[x] if x in self.dataset.symbs else None, sfs['GeneSymbol']))
        sfs = sfs[sfs['GeneSymbol'].isin(self.dataset.symbs)]
        print(str(np.sum(sfs['GeneSymbol'].isin(self.dataset.symbs)))+" of the splicing factor are in the dataset")
    
        return sfs

    def _get_matrices(self, sfs, timelag=1):
        ##calculate isoform abundances
        #
        isoabun = np.empty(shape=(self.dataset.shape[0], 0, self.dataset.shape[2]))
        transcript_idx = []
        geneid = []
        for gene, info in self.dataset.isoobj.normdf.items():
            if len(info['indices'])>1:
                transcript_idx.append(info['indices'])
                isoabun = np.concatenate((isoabun, info['normarr']), axis=1)
                for _ in range(info['normarr'].shape[1]):
                    geneid.append(gene)
        transcript_idx = reduce(lambda x,y: x+y, transcript_idx)
        transcript = self.dataset.transcript_id[transcript_idx]
        isoabun = np.concatenate(isoabun, axis=1) ##isoform abundance of isoforms shape = (isoforms, time points*replicates)

        #correlate splicing factor with isoform abundance
        ###splicing factor expression shape = (splicing factor, time points*replicates)
        sfs_expression = np.empty(shape=(self.dataset.shape[0], 0, self.dataset.shape[2]))
        sf_id = []
        for e,sf in enumerate(sfs['entrezid']):
            tmp = self.dataset.isoobj.normdf[sf]
            sfexp = np.expand_dims(np.sum(tmp['array'], axis=1),axis=1)
            sfs_expression=np.concatenate([sfs_expression, sfexp], axis=1)
            sf_id.append(sfs['GeneSymbol'].tolist()[e])

        sfs_expression = np.concatenate(sfs_expression, axis=1)

        ##filter for only IS genes
        isgenes_idx = np.where(np.isin(geneid, self.dataset.isoobj.is_result['gene']))[0]
        final_isoabun = isoabun[isgenes_idx,:]
        final_transcript = transcript[isgenes_idx]
        final_geneid = [geneid[x] for x in isgenes_idx]


        ##filter splicing factors variances
        highvar = np.where(np.var(sfs_expression, axis=1)>1)[0]
        final_sfsexp =sfs_expression[highvar,:] 
        final_sfid = [sf_id[x] for x in highvar]

        ##
        final_mat = np.concatenate([final_isoabun, final_sfsexp])
        mat_index = np.append(final_transcript, final_sfid)

        assert final_mat.shape[0] == mat_index.shape[0]
        return final_mat, mat_index, final_transcript, final_geneid, final_sfid

    def _corrcoef_loop(self, matrix):
        rows, cols = matrix.shape[0], matrix.shape[1]
        r = np.ones(shape=(rows, rows))
        p = np.ones(shape=(rows, rows))
        for i in range(rows):
            for j in range(i+1, rows):
                if self.method=="pearson":
                    r_, p_ = pearsonr(matrix[i], matrix[j])
                if self.method=="softdtw":
                    pass
                r[i, j] = r[j, i] = r_
                p[i, j] = p[j, i] = p_
        return r, p


    def _calculate_coexpression(self, final_mat, final_transcript, final_geneid, final_sfid):
        ##correlation with p values #######tmp read the results below
        sf_iso_mat, sf_iso_pval = self._corrcoef_loop(final_mat)
        
        ####multiple correction
        mapper = dict(zip(self.dataset.gene_id, self.dataset.symbs))

        sf_iso_07 = sf_iso_mat[final_transcript.shape[0]+1:, 0:final_transcript.shape[0]]
        sf_iso_dict = {"sf":[], "iso":[], "iso_gene":[], "iso_genesymb":[], "corr":[], "pval":[]}
        for sf in range(sf_iso_07.shape[0]):
            for iso in range(sf_iso_07.shape[0]):
                sf_iso_dict['sf'].append(final_sfid[sf])
                sf_iso_dict['iso'].append(final_transcript[iso])
                sf_iso_dict['iso_gene'].append(final_geneid[iso])
                sf_iso_dict['iso_genesymb'].append(mapper[final_geneid[iso]])
                sf_iso_dict['corr'].append(sf_iso_07[sf,iso])
                sf_iso_dict['pval'].append(sf_iso_pval[sf][iso])


        ##correct for multiple testing
        sf_iso_df = pd.DataFrame(sf_iso_dict)
        sf_iso_df = sf_iso_df.sort_values("pval")
        mp = pvalue_correction(sf_iso_df['pval'], method=self.padj_method)
        sf_iso_df['padj']=mp.corrected_pvals
        
        ##cleaning data, 
        sf_iso_df['cluster'] = list(map(self._map_to_clusters, sf_iso_df['iso_gene']))
        sf_iso_df = sf_iso_df[(np.abs(sf_iso_df['corr'])>self.corr_cutoff) & (sf_iso_df['padj']<0.05)].sort_values("iso_gene")
        sf_iso_df['bin_corr'] = [1 if x > 0 else -1 for x in sf_iso_df['corr']]

        sf_group = sf_iso_df.groupby(["bin_corr","iso_gene","iso_genesymb", "cluster"]).aggregate(lambda x: tuple(x))
        if sf_group.shape[0] == 0:
            return "No significant correlated SF is found."

        everycluster = pd.DataFrame(sf_group.groupby(["cluster", "bin_corr"]).aggregate(lambda x: tuple(x))['sf'].apply(lambda x: np.unique(reduce(lambda a,y: a+y, x)))).reset_index()
        #union_sf : the sf that are only found in one direction

        

        ####get list of SF that is positively and negatively correlated to at least one of the isoform in the cluster
        
        everycluster['diff_sf']=everycluster.groupby("cluster").aggregate(lambda x: list(x))['sf'].apply(lambda x: np.repeat(set(x[0]).symmetric_difference(set(x[1])),2) if len(x)>1 else [x[0]]).explode('sf')

        unique_sf = []
        union_sf =[]
        for x,s in everycluster.iterrows():
            unique_sf.append(set(s['sf']).intersection(set(s['diff_sf'])))
            union_sf.append(set(s['sf']).difference(set(s['diff_sf'])))

        everycluster['unique_sf']=unique_sf
        everycluster['union_sf']=union_sf
        return sf_iso_df, everycluster

    def _map_to_clusters(self, gene):
        #sf_iso_clu = []
        #for gene in sf_iso_df['iso_gene']:
        match=False
        for u,v in self.dataset.clusterobj.genelist_clusters.items():
            if gene in set(v):
                match=True
                break
        if match:
            return u
        #sf_iso_df['cluster'] = sf_iso_clu


    def coexpression_with_sf(self):
        sfs = self._get_splicing_factors()
        final_mat, mat_index, list_transcript, list_geneid, list_sfid = self._get_matrices(sfs)
        sf_iso_df, cluster_sf = self._calculate_coexpression(final_mat, list_transcript, list_geneid, list_sfid)

        return sf_iso_df, cluster_sf

def create_motifsobj():
        ########this only read the consensus
        # motifs_df = pd.read_csv(os.path.join(dir, "sfanalysis/motif_info.csv"), sep="\t")
        # print(motifs_df.head())
        # motifs_df['motifs']=motifs_df['motifs'].apply(eval)
        # instances_pre = motifs_df[motifs_df['motifname']==self.list_SF[0]]['motifs'].to_list()

        # motif_obj = []
        # for instance in instances_pre:
        #     instance = list(map(lambda x: x.replace("N",""), instance))
        #     print(instance)
        #     instances = list(map(Seq, instance))
        #     m = motifs.create(instances)

        #     motif_obj.append(m)
        # return motif_obj
        motif_obj = defaultdict(list)
        for mat in os.listdir(os.path.join(dir_path, "data/eCLIP_PWM")):
            if ".txt" in mat :
                with open(os.path.join(dir_path, "data/eCLIP_PWM",mat)) as handle:
                    record = motifs.parse(handle, "TRANSFAC")
                motif_obj[mat.split(".")[1]].append(record)

        with open(os.path.join(dir_path, "data/network/splicingfactors.txt"),"w") as f:
            for u in motif_obj.keys():
                f.write(u+"\n")
            f.close()

        # with open(os.path.join(dir,"spycone_pkg/spycone/data/motif_obj.pkl"), "wb") as f:
        #     pickle.dump(motif_obj, f, pickle.HIGHEST_PROTOCOL)
        
        # f = open(os.path.join(dir,"spycone_pkg/spycone/data/motif_obj.pkl"), "rb")
        # motif_obj = pickle.load(f)
        motif_pssm = defaultdict(list)
        for sf, motif_o in motif_obj.items():
            for e,bpmotif in enumerate(motif_o):
                if len(bpmotif)==1:          
                    
                    pwm=bpmotif[0].counts.normalize()
                    # background = {"A":(1-gc_ratio)/2,"C":gc_ratio/2,"G":gc_ratio/2,"T":(1-gc_ratio)/2}
                    # pssm=pwm.log_odds(background)
                    motif_pssm[sf].append([pwm])

        with open(os.path.join(dir,"spycone_pkg/spycone/data/motif_pwm.pkl"), "wb") as f:
            pickle.dump(motif_pssm, f, pickle.HIGHEST_PROTOCOL)

        return motif_obj, motif_pssm


class SF_motifsearch():
    '''
    Return the PSSM score of target SF and exons binding.

    Parameters
    ------------
    list_SF : list
        List of splicing factors / RBPs. E.g. the SF that is co-expressed to the isoform abundance. 

    list_genes : list
        List of genes you want to check for SF binding sites. E.g. the cluster of genes that the input SF is co-expressed. 

    gtf_df : str or dataframe
        gtf file path / dataframe

    gc_ratio : GC content that makes up the background for PSSM score calculation. Default=0.6.

    flanking : Flanking region size 

    Return
    -------
    Create an instance for SF motif enrichment analysis.

    '''
    def __init__(self, list_SF, list_genes, dataset, gtf, gc_ratio=0.6, flanking=400):
        self.list_SF = list_SF
        self.list_genes = list_genes
        self.dataset = dataset
        self.gtf = gtf
        self.flanking = flanking
        #self.motif_pssm = pickle.load(open(os.path.join(dir_path,"data/motif_pssm.pkl"), "rb"))
        self.motif_obj = pickle.load(open(os.path.join(dir_path,"data/motif_obj.pkl"), "rb"))
        self.motif_pwm = pickle.load(open(os.path.join(dir_path,"data/motif_pwm.pkl"), "rb"))
        self.motif_thres = pickle.load(open(os.path.join(dir_path,"data/motif_thres.pkl"), "rb"))

    ###perform motifs search
    def get_seq(self, exid, site="5p", bg=False):
        target = self.gtf[self.gtf['exon_id']==exid]
        if target.shape[0]!=0:
            targetpos=target[['seqname','start','end','exon_id']]
            
            targetpos = targetpos.set_index(['exon_id'])
            targetpos = targetpos[['seqname', 'start','end']].astype(str)

            # if bg:
            #     exon_len = int(targetpos['end'].values[0])-int(targetpos['start'].values[0])+1
            #     random_start = np.random.randint(1,45000000, size=1)
            #     url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr{targetpos['seqname'].values[0]};start={str(int(random_start)-self.flanking)};end={str(int(random_start)+exon_len+self.flanking)}"
            #     response = requests.get(url)
            # else:
            if site=="5p":
                url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr{targetpos['seqname'].values[0]};start={str(int(targetpos['start'].values[0])-self.flanking)};end={str(int(targetpos['start'].values[0])+self.flanking)}"
            
            elif site=="3p":
                url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr{targetpos['seqname'].values[0]};start={str(int(targetpos['end'].values[0])-self.flanking)};end={str(int(targetpos['end'].values[0])+self.flanking)}"

            response = requests.get(url)
            
            if 'dna' in response.json().keys():
                return response.json()['dna']
            else:
                return None

    def _exon_search(self, exons, sf, site):
        sf_results=defaultdict(list)
        motifs_pwm = self.motif_pwm[sf]
        for exon in exons:
            exon_seq = self.get_seq(exid=exon, site=site)
            if exon_seq is None:
                continue
            else:
                for e,pwm in enumerate(motifs_pwm):  
                    if len(pwm)==1:          
                        test_seq = Seq(exon_seq)
                        submotstr = f'{sf}_{e}'
                        gc_ratio = round(GC(test_seq)/100,2)
      
                        background = {"A":(1-gc_ratio)/2,"C":gc_ratio/2,"G":gc_ratio/2,"T":(1-gc_ratio)/2}
                        pssm=pwm[0].log_odds(background)

                        # distribution = pssm[0].distribution(background=background, precision=10**4)
                        # threshold = distribution.threshold_patser()
                        submotifs = None
                        submotstr = None
                        # for pos, sc in pssm[0].search(test_seq, self.threshold):
                        #     dire = np.sign(pos)

                        # #relpos_start = np.abs(pos) - self.flanking
                        # #relpos_end = np.abs(pos) - len(exon_seq) + self.flanking
                        #     submotifs = f"{exon};{dire};{pos};{sc}"
                        #     submotstr = f'{sf}_{e}'
                        
                        submotstr = f'{sf}_{e}'
                        allscores = pssm.calculate(test_seq)
                        # if self.motif_thres[submotstr] is not None:
                        #     allscores = pssm[0].search(test_seq, threshold=self.motif_thres[submotstr])
                            
                        # for pos, sc in enumerate(allscores):
                        #     submotifs = f"{exon};1;{pos};{sc}"
                        #     sf_results[submotstr].append(submotifs)
                        #### this only save the scores not the position  
                        this_threshold = self.motif_thres[submotstr] 
                        sf_results[submotstr] = [sc if not np.isnan(sc) else next if sc > this_threshold else next for sc in allscores]
                    else:
                        continue
            gc.collect()
        return sf_results

    def _parallel_search(self, sf, exon_dicts, site, bg=False):
        sf_results = []
        motif_str = []
        #for each SF search in every exons
        #for gene, exons in self.target_exons_dict.items():
        # with concurrent.futures.ProcessPoolExecutor() as executor:
        #     for gene, each_gene_res in zip(exon_dicts.values(), executor.map(self._exon_search, exon_dicts.values(), np.repeat(sf, len(exon_dicts.values())))):
        #         sf_results.append(each_gene_res)
      
        sf_results = Parallel()(delayed(self._exon_search)(exons=v, sf=sf, site=site) for v in exon_dicts.values())
        # for v in exon_dicts.values():
        #     res = self._exon_search(v, sf, site=site)
        #     sf_results.append(res)
  
        sf_per_results = defaultdict(list)
        for eachgene in sf_results:
            for eachmot, eachres in eachgene.items():
                if eachmot is not None:
                    for x in eachres:
                        sf_per_results[eachmot].append(x)
        return sf_per_results
            
    def search_motifs(self, exons="loss", site="5p"):
        """
        Calculate PSSM scores on all exons in the input list of genes. 

        Parameters
        -----------
        

        """
        #get exons from IS results
        target_exons = self.dataset.isoobj.is_result[self.dataset.isoobj.is_result['gene_symb'].isin(self.list_genes)]
        target_exons_dict = defaultdict(list)
        n_exons =0
        if exons=="loss":
            for x in target_exons.iterrows():
                target_exons_dict[x[1]['gene_symb']]=x[1]['loss_exons']
                n_exons+=len(x[1]['loss_exons'])
        else:
            for x in target_exons.iterrows():
                target_exons_dict[x[1]['gene_symb']]=x[1]['gain_exons']
                n_exons+=len(x[1]['gain_exons'])

        self.target_exons_dict = target_exons_dict

        
        sf_results = defaultdict(list)
        motif_str =defaultdict(list)
        #get motif objects from list SF
        # motif_obj, motif_pssm = self.create_motifsobj()
        # self.motif_obj = motif_obj
        # self.motif_pssm = motif_pssm
        # parallelize motif searching 
        ## return a list of result in the following format:
         ### exonID;direction of the strand;matched position of the seq;score from biopython;relative pos to the start site (including flanking region);relative pos to the end + flanking region; length of the exon; 
        for sf in self.list_SF:
            each_sf_res = self._parallel_search(sf, target_exons_dict, site=site)
            sf_results[sf].append(each_sf_res)
        # iterate the results to create a dictionary where keys are SF and values are the results for each exon
        ### note that each SF might have multiple motifs, in this case the motifs are denoted as SF_1, SF_2 from the same SF.
        self.motif_result = sf_results
        gc.collect()
        return sf_results, n_exons

    def search_background_motifs(self, site="5p"):
        allexons = self.gtf["exon_id"].to_numpy()
        #get exons from IS results
        if "self.target_exons_dict" not in locals():
            ValueError("Please run search_motifs first.")

        bg_exons_dict=defaultdict(list)
        n_exons=0
        target_exons = self.dataset.isoobj.is_result[self.dataset.isoobj.is_result['gene_symb'].isin(self.list_genes)]
        for x in target_exons.iterrows():
             bg_exons_dict[x[1]['gene_symb']]=x[1]['constant_exons']
             n_exons+=len(x[1]['constant_exons'])

        self.bg_exons_dict =  bg_exons_dict

        print("get exons bg dict")
        #get motif objects from list SF
        sf_results = defaultdict(list)
        motif_str =defaultdict(list)
        #get motif objects from list SF
        # motif_obj, motif_pssm = self.create_motifsobj()
        # self.motif_obj = motif_obj
        # self.motif_pssm = motif_pssm
        # parallelize motif searching 
        ## return a list of result in the following format:
         ### exonID;direction of the strand;matched position of the seq;score from biopython;relative pos to the start site (including flanking region);relative pos to the end + flanking region; length of the exon; 
        
       # with concurrent.futures.ProcessPoolExecutor() as executor:
        #    for each_sf_res, sf in zip(self.list_SF, executor.map(self._parallel_search, self.list_SF)):
        for sf in self.list_SF:
            each_sf_res = self._parallel_search(sf, bg_exons_dict, site=site)
            sf_results[sf].append(each_sf_res)
        # iterate the results to create a dictionary where keys are SF and values are the results for each exon
        ### note that each SF might have multiple motifs, in this case the motifs are denoted as SF_1, SF_2 from the same SF.

        self.motif_background = sf_results
        gc.collect()
        return sf_results, n_exons
    

    def _df_forvis(self, lossres, gainres, bg=None):
        '''
        MannWhitney U test or Fisher exact test

        '''
        def _checkscore(x):
            if np.isfinite(x):
                return x
            else:
                return 0

        vis_sf = []
        ldf = {#'allpos':[],
        'logscores':[],
        #'allscores':[],
        'sf_motif':[],
        'threshold':[]}
        for sf,v in lossres[0].items():
            for eachmot, exons in v[0].items():
                for ex in exons:
                    # signs = np.sign(float(ex.split(";")[2]))
                    # signss = 1 if signs == 0 else signs
                    # signsss = "+ strand" if signss==1 else "- strand"
                    #ldf['allpos'].append(np.abs(float(ex.split(";")[2]))-self.flanking)
                    thisscore=float(ex)
                    ldf['logscores'].append(thisscore)
                    # df['alldir'].append(signsss)
                    ldf['sf_motif'].append(eachmot)
                    ldf['threshold'].append(self.motif_thres[eachmot])

        ldf = pd.DataFrame(ldf)
        ldf['rawtype'] = "lost exons"
        ldf['type'] = f"LE ({lossres[1]})"
        vis_sf.append(ldf)

        gdf = {#'allpos':[],
        'logscores':[],
        #'allscores':[],
        'sf_motif':[],
        'threshold':[]}
        #countexons=[]
        for sf,v in gainres[0].items():
            for eachmot, exons in v[0].items():
                for ex in exons:
                    #countexons.append(ex)
                    # signs = np.sign(float(ex.split(";")[2]))
                    # signss = 1 if signs == 0 else signs
                    # signsss = "+ strand" if signss==1 else "- strand"
                    #gdf['allpos'].append(np.abs(float(ex.split(";")[2]))-self.flanking)
                    thisscore=float(ex)
                    gdf['logscores'].append(thisscore)
                    #gdf['allscores'].append(10**thisscore)
                    # df['alldir'].append(signsss)
                    gdf['sf_motif'].append(eachmot)
                    gdf['threshold'].append(self.motif_thres[eachmot])

        gdf = pd.DataFrame(gdf)
        gdf['rawtype'] = "gained exons"
        gdf['type'] = f"GE ({gainres[1]})"
        
        vis_sf.append(gdf)

        if bg:
            bgdf = {#'allpos':[],
            'logscores':[],
        #'allscores':[],
        'sf_motif':[],
        'threshold':[]}
            for sf,v in bg[0].items():
                for eachmot, exons in v[0].items():
                    for ex in exons:
                        #countexons.append(ex.split(";")[0])
                        # signs = np.sign(float(ex.split(";")[2]))
                        # signss = 1 if signs == 0 else signs
                        # signsss = "+ strand" if signss==1 else "- strand"
                        #bgdf['allpos'].append(np.abs(float(ex.split(";")[2]))-self.flanking)
                        thisscore=float(ex)
                        bgdf['logscores'].append(thisscore)
                        #bgdf['allscores'].append(10**thisscore)
                        # bgdf['alldir'].append(signsss)
                        bgdf['sf_motif'].append(eachmot)
                        bgdf['threshold'].append(self.motif_thres[eachmot])

            bgdf = pd.DataFrame(bgdf)
            bgdf['rawtype'] = "background"
            bgdf['type'] = f"UE ({bg[1]})"
            vis_sf.append(bgdf)
        
        alldf = pd.concat(vis_sf)
        alldf = alldf[alldf['threshold'].notnull()]
        alldf=alldf[(np.isfinite(alldf['logscores']))].reset_index()
        #alldf['titles'] = alldf['sf_motif']+", "+alldf['alldir']
        return alldf

    def _get_pvals(self, alldf, method="mannwhitney"):
        alldf = alldf[alldf['logscores']>alldf['threshold']]
        ### cal pvalue
        if method=="mannwhitney":
            pvals =defaultdict(list)
            #for stra in ["- strand", "+ strand"]:
                #ddf = df[df['alldir']==stra].reset_index()
                #bgddf = bgdf[bgdf['alldir']==stra].reset_index()
            for sf in alldf['sf_motif'].unique():
                idx=alldf[(alldf['sf_motif']==sf) & (alldf['rawtype']=="gained exons")].index
                lidx=alldf[(alldf['sf_motif']==sf) & (alldf['rawtype']=="lost exons")].index
                bgidx=alldf[(alldf['sf_motif']==sf) & (alldf['rawtype']=="background")].index
                gexons = [alldf['logscores'][x] for x in idx]
                lexons = [alldf['logscores'][x] for x in lidx]
                bgs = [alldf['logscores'][x] for x in bgidx]
                
                try:
                    pvals[sf].append(mannwhitneyu(lexons, bgs, alternative="greater")[1]) #[lost, gain]
                except ValueError:
                    pvals[sf].append(1)

                try:
                    pvals[sf].append(mannwhitneyu(gexons, bgs, alternative="greater")[1]) #[lost, gain]
                except ValueError:
                    pvals[sf].append(1)

                try:
                    pvals[sf].append(mannwhitneyu(lexons, gexons, alternative="greater")[1]) #[lost, gain]
                except ValueError:
                    pvals[sf].append(1)


        return pvals

    def vis_motif_density(self, res=None, style="boxplot", bg=None, method="mannwhitney", pval_cutoff=None, merge_motif=False, alldf=None, pvals=None, palette=None, sf_list=None):
        if alldf is None and pvals is None:
            if res is None:
                raise("Please input result from search_motif() as res")
            alldf, pvals = self._df_pval_forvis(res, bg=bg, method=method)
        

        filtered_alldf = alldf[alldf['logscores']>alldf['threshold']]
        
        def print_pvalues(g, i, j, hdiff, pvid):
            for x,ax in enumerate(g.axes_dict.values()):
                text=pvals[ax.get_title()][pvid]

                px = []
                ph=[]
                phmin = []
                for p in alldf['type'].unique():
                    tmp = filtered_alldf[(filtered_alldf['type']==p)&(filtered_alldf['sf_motif']==ax.get_title())]['logscores']
                    if not np.isnan(tmp.max()):
                        ph.append(tmp.max()) ##max of each group
                        phmin.append(filtered_alldf[(filtered_alldf['type']==p)&(filtered_alldf['sf_motif']==ax.get_title())]['logscores'].min())
                    else:
                        ph.append(0)
                        phmin.append(0)

                pxx = (i+j)/2
                phh = max(ph[i],ph[j])
                
                props = {'connectionstyle':'bar','arrowstyle':'-',
                            'shrinkA':20,'shrinkB':20,'linewidth':2}
                ax.plot([i,i,j,j], [phh+hdiff, phh+hdiff+0.2, phh+hdiff+0.2, phh+hdiff], lw=1.5, c="k")
                ax.text(s=f"p={str('{:.2e}'.format(text))}", x=pxx-0.3, y=phh+hdiff+0.3)
                ax.set_ylim((np.min(phmin)-1, phh+5))

        if sf_list:
            filtered_alldf = filtered_alldf[filtered_alldf['sf_motif'].isin(sf_list)]
            
        if pval_cutoff:
            sigs_sf = []
            for u,v in pvals.items():
                checksig = False
                for vv in v:
                    if vv<=pval_cutoff:
                        checksig=True

                if checksig == True:
                    sigs_sf.append(u)

            
            sigs_df = filtered_alldf[filtered_alldf['sf_motif'].isin(sigs_sf)]

            palette1 = defaultdict(str)
            for e,x in enumerate(sigs_df['type'].unique()):
                palette1[x] = list(palette.values())[e]
                
            if style=="boxplot":
                plt.figure(figsize=[8,9])
                g= sns.catplot(data=sigs_df, x="type",  y="logscores",  col="sf_motif", col_wrap=5, kind="box",linewidth=3, palette=palette1, sharey=False)
                g.set_titles("{col_name}")
                print_pvalues(g, 0, 1, 0.5, 2) #loss vs gain
                print_pvalues(g, 1, 2, 1, 1) #gain vs bg
                print_pvalues(g, 0, 2, 2, 0) #loss vs bg
                # for x,ax in enumerate(g.axes_dict.values()):
                #     thistitle = ax.get_title()
                #     pvvv = pvals[thistitle][3]
                #     ax.set_title(thistitle+"\n"+f"kruskal p-value: {str('{:.2e}'.format(pvvv))}", y=0.95)
                g.set_axis_labels(x_var="Type", y_var="log(Score)")
                g.set_xticklabels()

        else:
            if style=="boxplot":
                palette1 = defaultdict(str)
                for e,x in enumerate(filtered_alldf['type'].unique()):
                    palette1[x] = list(palette.values())[e]

                plt.figure(figsize=[8,9])
                g= sns.catplot(data=filtered_alldf, x="type",  y="logscores",  col="sf_motif", col_wrap=5, kind="box",linewidth=3, palette=palette1, sharey=False)
                g.set_titles("{col_name}")
                print_pvalues(g, 0, 1, 0.5, 2) #loss vs gain
                print_pvalues(g, 1, 2, 1, 1) #gain vs bg
                print_pvalues(g, 0, 2, 1.5, 0) #loss vs bg
                # for x,ax in enumerate(g.axes_dict.values()):
                #     thistitle = ax.get_title()
                #     pvvv = pvals[thistitle][3]
                #     ax.set_title(thistitle+"\n"+f"kruskal p-value: {str('{:.2e}'.format(pvvv))}", y=0.95)
                g.set_axis_labels(x_var="Type", y_var="log(Score)")
                g.set_xticklabels()
            
        
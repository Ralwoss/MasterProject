labels = "preparations/10_0_10000/labels"
interaction_matrices = "preparations/10_0_10000/InteractionMatrices"
model = "model/5_0_10kb_batch_normalized.h5"
hic_matrix = "data/GM12878-MboI-allreps-filtered.10kb.cool"
TAD_domains = "data/GM12878-MboI-allreps-filtered-TAD-domains.bed"
binsize = 10000
window_size = 25
overlap_size = 10
save_suffix = "center"
save_preparation_id = str(window_size)+'_'+ str(overlap_size)+'_'+ str(binsize)+'_'+str(save_suffix)
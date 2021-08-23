
hic_matrix = "data/GM12878-MboI-allreps-filtered.10kb.cool"
TAD_domains = "data/GM12878-MboI-allreps-filtered-TAD-domains.bed"
binsize = 10000
window_size = 15
overlap_size = 5
detection_range = 2
save_suffix = "center"
balanceData = 0 # 0: nothing, 1: class weights, 2: oversampling, 3: undersampling
save_preparation_id = str(window_size)+'_'+ str(overlap_size)+'_'+ str(binsize)+'_'+str(save_suffix)
interaction_matrices_pos = "preparations/"+ save_preparation_id+"/InteractionMatricesPos"
interaction_matrices_neg = "preparations/"+ save_preparation_id+"/InteractionMatricesNeg"
interaction_matrices_pos_npz = 'preparations/' + save_preparation_id + '/InteractionMatricesPosNp.npz'
interaction_matrices_neg_npz = 'preparations/' + save_preparation_id + '/InteractionMatricesNegNp.npz'
windows_bed = 'preparations/' + save_preparation_id + '/windows.bed'
results_dir = 'results/'
balancing = {0:"", 1:"class_weights", 2:"oversampling", 3:"undersampling"}
model = "model/" + save_preparation_id + "_"+ balancing[balanceData] +"_batch_normalized.h5"
TRAINCHORMS = ['chr5', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
VALCHROMS = ['chr1', 'chr3', 'chr7']
TESTCHROMS = ['chr2', 'chr4', 'chr6']
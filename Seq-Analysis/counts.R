dirname <- "/home/RNA_Seq/scRNA_seq_data/seq_data/Control/10x_data/outs/filtered_feature_bc_matrix/"
counts_matrix_filename = paste0(dirname,"/home/RNA_Seq/scRNA_seq_data/seq_data/Control/10x_data/outs/filtered_feature_bc_matrix/")
counts <- Read10X(data.dir = "/home/RNA_Seq/scRNA_seq_data/seq_data/Control/10x_data/outs/filtered_feature_bc_matrix")

rule seurat:
	input:
		metaQuery = "resources/data/{dataset}/query_meta.rda",
		expQuery = "resources/data/{dataset}/query_exp.rda",
		metaTrain = "resources/data/{dataset}/train_meta.rda",
		expTrain = "resources/data/{dataset}/train_exp.rda"
	output:
		transfdata_umap = "results/{dataset}/seurat/transfdata_umap.png",
		refAndQuery_umap = "results/{dataset}/seurat/refAndQuery_umap.png",
		query_violin = "results/{dataset}/seurat/query_violin.png",
		query_heatmap = "results/{dataset}/seurat/query_heatmap.png"
	params:
		Cell_type_colname = config["Cell_type_colname"],
		Cell_ID_colname = config["Cell_ID_colname"],
		Cell_description_colname = config["Cell_description_colname"],
		Seurat_features = config["Suerat_features"],
		fig_width = config["fig_width"],
		fig_height = config["fig_height"]
	script:
		"scripts/seurat.R"

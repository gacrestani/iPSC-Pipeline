rule fusca:
	input:
		metaQuery = "resources/data/{dataset}/query_meta.rda",
		expQuery = "resources/data/{dataset}/query_exp.rda",
		metaTrain = "resources/data/{dataset}/train_meta.rda",
		expTrain = "resources/data/{dataset}/train_exp.rda"
	output:
		umap_queryAndReference = "results/{dataset}/fusca/umap_queryAndReference.png"
	params:
		Cell_type_colname = config["Cell_type_colname"],
		Cell_ID_colname = config["Cell_ID_colname"],
		Cell_description_colname = config["Cell_description_colname"],
        fig_width = config["fig_width"],
		fig_height = config["fig_height"]
	script:
		"scripts/fusca.R"

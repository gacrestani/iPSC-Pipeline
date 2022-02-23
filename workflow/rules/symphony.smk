rule symphony:
	input:
		metaQuery = "resources/data/{dataset}/query_meta.rda",
		expQuery = "resources/data/{dataset}/query_exp.rda",
		metaTrain = "resources/data/{dataset}/train_meta.rda",
		expTrain = "resources/data/{dataset}/train_exp.rda"
	output:
		reference_umap = "results/{dataset}/symphony/reference_umap.png",
		umap_mixed_queryAndReference = "results/{dataset}/symphony/umap_mixed_queryAndReference.png",
		umap_query = "results/{dataset}/symphony/umap_query.png"
	params:
		Cell_type_colname = config["Cell_type_colname"],
		Cell_ID_colname = config["Cell_ID_colname"],
		Cell_description_colname = config["Cell_description_colname"],
		fig_width = config["fig_width"],
		fig_height = config["fig_height"]
	script:
		"scripts/symphony.R"

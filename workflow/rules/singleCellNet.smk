rule singleCellNet:
	input:
		metaQuery = "resources/data/{dataset}/query_meta.rda",
 		expQuery = "resources/data/{dataset}/query_exp.rda",
 		metaTrain = "resources/data/{dataset}/train_meta.rda",
		expTrain = "resources/data/{dataset}/train_exp.rda"
	output:
		assessment_curves = "results/{dataset}/singleCellNet/assessment-curves.png",
		assessment_metrics = "results/{dataset}/singleCellNet/assessment_metrics.png",
		test_data_heatmap = "results/{dataset}/singleCellNet/test-data_heatmap.png",
		attr_plot = "results/{dataset}/singleCellNet/atrr-plot.png",
		average_top_genepair = "results/{dataset}/singleCellNet/average-top-genepair.png",
		query_classification_heatmap = "results/{dataset}/singleCellNet/query_classification-heatmap.png",
		query_violin_plot = "results/{dataset}/singleCellNet/query-violin-plot.png",
		query_subcluster_violin_plot = "results/{dataset}/singleCellNet/query-subcluster-violin-plot.png",
		query_umap = "results/{dataset}/singleCellNet/query_umap.png",
		score_heatmap = "results/{dataset}/singleCellNet/score_heatmap.png"
	params:
		Cell_type_colname = config["Cell_type_colname"],
		Cell_ID_colname = config["Cell_ID_colname"],
		Cell_description_colname = config["Cell_description_colname"],
		Cell_type_to_analyze = config["Cell_type_to_analyze"],
		fig_width = config["fig_width"],
		fig_height = config["fig_height"]
	script:
		"scripts/singleCellNet.R"

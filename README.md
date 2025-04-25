# Camila_Cuadrado_Final_Project

This project is focused on parsing, analyzing, and visualizing genetic variant data from VCF files using the vcfR, dplyr, ggplot2, and other supporting R packages. The goal is to extract functional annotations (using the ANN field from snpEff), summarize variant effects, and visualize their distribution by type and gene impact.

📁 Project Structure
	•	Sample_1.vcf, Sample_2.vcf, …: Input VCF files containing annotated variant data.
	•	R Markdown / .R scripts: Scripts for data extraction, transformation, statistical analysis, and visualization.
	•	README.md: Overview of the pipeline, tools used, and how to run the analysis.

📦 Required R Packages
install.packages(c("vcfR", "dplyr", "ggplot2", "stringr", "ggpubr"))
🔁 Workflow Overview

1. Load and Parse VCF
vcf_data <- read.vcfR("Sample_1.vcf", verbose = TRUE)
annotations <- extract.info(vcf_data, "ANN")
2. Extract and Format Annotations
	•	Split complex ANN field.
	•	Parse each annotation into separate fields.
	•	Convert the list into a data frame with named columns.

3. Count Variant Effects
	•	Summarize the number of mutations per Effect type.
	•	Join effect counts back to the full annotation table.

4. Visualizations
	•	Bar plots of variant effects and gene impacts using ggplot2.
	•	Faceted views per impact type.
	•	Annotated statistical results (e.g., Chi-square test).

5. Statistical Analysis
	•	Chi-square test to assess distribution of variant effects.
	•	T-tests or ANOVA (if applicable) for group comparisons.

6. Summarize Gene-Level Information
	•	Group variant data by gene.
	•	Create gene type categories (e.g., M, F, P).
	•	Visualize total variants per gene or category.

📊 Example Plot
ggplot(effect_counts_df, aes(x = Effect, y = Count, fill = Effect)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("Variant Effect") +
  ylab("Count") +
  ggtitle("Distribution of Variant Effects")
  📌 Notes
	•	The script processes only the first annotation for each variant (ANN[1]) for simplicity.
	•	Designed for snpEff-compatible VCF files.
	•	You can easily extend this pipeline for additional VCFs by wrapping the script into a loop or function.

🧪 Contact

If you need help customizing this pipeline for your data or adding new statistical analyses, feel free to reach out!

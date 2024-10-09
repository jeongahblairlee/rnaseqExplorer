# <img src="https://raw.githubusercontent.com/jeongahblairlee/rnaseqExplorer/refs/heads/main/notebook/logo_pic.png" alt="logo" style="width:400px; height: auto;">



Welcome to the **rnaseqExplorer**. This software is designed to help users easily upload and analyze RNA sequencing data. The application provides a user-friendly interface that allows you to visualize your data through various plots and explore the results in an intuitive way.

**rnaseqExplorer** provides functionalities for data analysis, including filtering, normalization, and understanding gene expression through visualizations like:

- 🎻 Violin plots
- 🗺️ Heatmaps
- 📈 PCA
- 🔍 Gene annotation

With its dynamic, interactive design, **rnaseqExplorer** is an essential tool for RNA-seq dataset analysis. It empowers bench biologists to conduct exploratory data analysis with ease, while delivering in-depth insights for seasoned data analysts. 💡

## Getting Started

To use this application, clone the repository and follow the setup instructions in the documentation.

### In R

You can install the package using the following command:

```r
devtools::install_github("jeongahblairlee/rnaseqExplorer")
library(rnaseqExplorer)
launchApp()
```

### In web interface

You can also access the application via the web interface:

**https://jeongahblairlee.shinyapps.io/rnaseqExplorer/**

> **Warning:** In case the dataset is large, I recommend using R due to limitations from the web server.




## Workflow

In the **Workflow**, you can see the overall workflow of the RNAseq Explorer application:

![Workflow](https://raw.githubusercontent.com/jeongahblairlee/rnaseqExplorer/refs/heads/main/notebook/workflow.png)

## How to Use

### Step 1: Upload Your Data

To get started, upload your RNA sequencing data using the provided functionality:

![Upload Data](https://raw.githubusercontent.com/jeongahblairlee/rnaseqExplorer/refs/heads/main/notebook/function1.png)

### Step 2: Explore the Data

Once your data is uploaded, you can explore it through various visualization options:

#### 2-1. Violin Plot

Generate violin plots to visualize the distribution of your data:

![Violin Plot](https://raw.githubusercontent.com/jeongahblairlee/rnaseqExplorer/refs/heads/main/notebook/function2.png)

#### 2-2. Heatmap

Create heatmaps for a detailed view of the gene expression levels:

![Heatmap](https://raw.githubusercontent.com/jeongahblairlee/rnaseqExplorer/refs/heads/main/notebook/function3.png)

#### 2-3. PCA

Perform Principal Component Analysis (PCA) to understand the variance in your data:

![PCA](https://raw.githubusercontent.com/jeongahblairlee/rnaseqExplorer/refs/heads/main/notebook/function4.png)


## Contributions

We welcome contributions to enhance the functionality and usability of RNAseq Explorer. Please feel free to submit a pull request or open an issue to discuss improvements.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For any inquiries or issues, please contact [Jeongah Lee](jeongahblair@gmail.com). 

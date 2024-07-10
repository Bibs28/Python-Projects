# Python-Projects

## Introduction
## Welcome to my Python portfolio where I would have hopefully been able to demostrate and apply my knowledge of Python
- **If theres any feedback you can give to improve my current understanding and make my projects better do not hesistate to email me: nabibmalik@gmail.com

## Project 1 : Expense Tracker  [View Expense Tracker on GitHub](https://github.com/Bibs28/Python-Projects/blob/5cfefe79be4fcfe4e042996bb112b62e5621839d/expense_tracker.py)

## Features

- **Add Expenses:** Users can add their expenses to the tracker by entering the date, description, and amount spent. The expenses are then stored and used for calculations.

- **View Expenses:** The application provides an option to view all the expenses listed with their respective dates and amounts.

- **View Disposable Income:** Users can check their disposable income, which is calculated by subtracting the total expenses from the entered income. If the disposable income is negative, the user will receive a helpful message suggesting that they slow down their spending.

- **Detailed Breakdown:** There is an additional option to view a detailed breakdown of all the listed expenses along with the disposable income.



## Project 2:  Slot Machine Game  [View Slot Machine on GitHub](https://github.com/Bibs28/Python-Projects/blob/5cfefe79be4fcfe4e042996bb112b62e5621839d/Slot_machine)

Welcome to the Slot Machine Game! This is a simple text-based casino-style game developed in Python. The game allows you to try your luck by spinning the slot machine and matching symbols to win virtual coins.

## Instructions

1. Enter your desired initial balance when prompted to set up your custom starting point.
2. Place bets on each spin by entering the amount you want to wager. You can enter 0 to quit the game.
3. The slot machine consists of three spinning reels, each displaying one of six symbols: Cherry, Bell, Lemon, Orange, Star, and Skull.
4. If you get three matching symbols, you win based on the symbol's payout ratio. For example, getting three Cherries pays out 3 times your bet.
5. If you get two matching symbols, you receive your bet back.
6. If there are no matching symbols, you lose your bet.
7. The game continues until you decide to quit or until your balance runs out.

## Features

- **Custom Initial Balance:** You can set your preferred initial balance to start the game.
  
- **Interactive Gameplay:** The game prompts you to place bets and displays the result of each spin.
  
- **Random Symbol Generation:** The slot machine generates random symbols for each spin.
  
- **Payout Calculations:** The game calculates payouts based on the matched symbols and payout ratios.



## Project 3: Stock Ticker Dashboard: Stock Ticker Dashboard  [View Stock Ticker on GitHub](https://github.com/Bibs28/Python-Projects/blob/bdedf3747ec1cc7f76c6a4902c0efe520543c799/Stock%20Ticker%20dashboard.ipynb)

The Stock Ticker Dashboard is a web application built using Dash that allows users to visualize stock price data for multiple stocks over a selected date range. This dashboard provides an interactive and user-friendly interface for exploring stock market trends.

## Prequisites

1. Before using the Stock Ticker Dashboard, ensure that you have the following prerequisites installed on your system:
- [Jupyter Notebook](https://jupyter.org/install)
  
2. Download lates version of Nasdaq stock screener as a csv and make sure it is saved within the same folder as your stock ticker

## Features

- **Select Stocks:** Users can select one or more stocks from a predefined list to visualize their price data.

- **Date Range:** Users can choose a date range, specifying a start and end date, to view historical stock prices within that period.

- **Interactive Graphs:** The application displays stock price data as interactive graphs, making it easy to analyze trends and patterns.

- **Submit and Update:** By clicking the "Submit" button, users can update the graph to reflect their stock and date range selections.

- **Clear Visualization:** Users can clear the graph to start fresh and select new stocks and date ranges.

## Project 4: Heart Disease Prediction Using Machine Learning [View Heart Disease Prediction on GitHub](https://github.com/Bibs28/Python-Projects/blob/c696807c7eec29580bb0e0f7a9c88411c93fedf2/ML_Assignment.ipynb)

Welcome to the Heart Disease Prediction project! This Python-based project applies machine learning techniques to predict heart disease using a well-structured dataset. The aim is to leverage data analytics for early detection and diagnosis, enhancing healthcare outcomes through technological innovation.

### Features

- **Data Preprocessing:** Rigorous data cleaning, encoding categorical data, and scaling numerical features to prepare the dataset for effective modelling.
- **Exploratory Data Analysis (EDA):** Utilizing Matplotlib and Seaborn for visual analysis, creating histograms, box plots, and correlation matrix heatmaps to explore data relationships and distributions.
- **Model Building and Evaluation:** Training models such as K-Nearest Neighbors, Random Forest, and Logistic Regression. Models are evaluated based on accuracy, precision, recall, F1-score, and Matthews correlation coefficient to ensure robust predictions.
- **Principal Component Analysis (PCA):** Dimensionality reduction to visualize the most influential features and understand the data structure in two-dimensional space.

### Ethical Considerations

- Ensuring privacy and data security, this project uses anonymized data, emphasizing ethical AI practices and responsible data handling.

### Computational Details

- **Environment:** Python 3.7.6 and Google Colab for high-performance computing.
- **Libraries:** Scikit-learn for machine learning, Pandas for data manipulation, and visualization tools Matplotlib and Seaborn.

### Results and Future Directions

The project highlights the potential of machine learning in preventive health diagnostics, indicating areas for improvement such as class imbalance resolution and model optimization. Future enhancements may include integrating more sophisticated machine learning models and expanding data sources for better generalizability and accuracy.

### Conclusion

Machine learning's capacity to analyze complex patterns makes it invaluable for medical diagnostics. This project not only demonstrates the practical application of machine learning in predicting heart disease but also sets the stage for future advancements in medical technology.

## Project 5: Missing Biological Data Challenge [View Missing Data Challenge on GitHub]([BIO732P Missing Data Code.R](https://github.com/Bibs28/Python-Projects/blob/caf8698895f8e8e1a25c4b48a0a5118b4fb29191/BIO732P%20Missing%20Data%20Code.R))

Welcome to the Missing Biological Data Challenge project! This R-based project focuses on addressing missing values in biological datasets, specifically gene expression data in model organisms like mice. The goal is to apply various imputation methods to enhance data quality and reliability for accurate differential gene expression analysis.


### Data Preprocessing
- **Data Loading and Cleaning:** Loading datasets (`countdata5`, `countdata10`, and `countdata30`) and transposing them so that genes become columns and samples (mice) become rows. Columns containing only zero values are identified and removed to ensure meaningful analysis.
- **Visualizing Missing Data:** Utilizing the `vis_miss` function from the `visdat` package to create heatmap-like visualizations, highlighting the extent and pattern of missing values.

### Imputation Methods
1. **Mean Imputation:** Replacing each missing value with the mean of the non-missing values in the same column.
2. **Hot-Deck Imputation:** Replacing missing values with observed values from similar records.
3. **k-Nearest Neighbors (kNN) Imputation:** Using the `impute.knn` function to find the k-nearest neighbors of an observation with missing data and impute the missing values with the average of these neighbors.
4. **Random Forest Imputation:** Utilizing the `missForest` package for iterative imputation based on random forests, dividing the dataset into manageable subsets to handle computational complexity.

### Model Evaluation
- **Error Metrics:** Evaluating the accuracy of imputation methods using Mean Absolute Error (MAE) and Root Mean Squared Error (RMSE) to compare imputed values to actual values.
- **Evaluation Process:** Splitting datasets into training (70%) and testing (30%) sets, introducing artificial missing values, applying imputation methods, and calculating MAE and RMSE.

### Differential Expression Analysis
- **Data Preparation:** Loading metadata and assigning groups based on geographical latitude, ensuring distinct groups for analysis.
- **DESeq2 Analysis:** Performing differential expression analysis using DESeq2, creating a `DESeqDataSet` object and analyzing gene expression differences.
- **Identification of Significant Genes:** Identifying significant genes based on adjusted p-values (padj < 0.05) and visualizing them with volcano plots using the `EnhancedVolcano` package.

### Visualization and Comparison
- **Venn Diagrams:** Creating Venn diagrams to visualize the overlap of upregulated and downregulated genes between different imputation methods and datasets using `ggvenn`.
- **Bar Charts:** Comparing MAE and RMSE of different imputation methods across datasets using bar charts for performance evaluation.

### Ethical Considerations
This project emphasizes responsible data handling, ensuring privacy and data security by using anonymized datasets, and adhering to ethical AI practices.

### Computational Details
- **Environment:** R 4.0.2
- **Libraries:** Includes `naniar`, `ggplot2`, `dplyr`, `caret`, `VIM`, `tibble`, `impute`, `visdat`, `missForest`, `BiocManager`, `gprofiler2`, `EnhancedVolcano`, `ggvenn`, `RColorBrewer`, and `gridExtra`.

### Results and Future Directions
The project demonstrates effective handling of missing data, crucial for accurate downstream analyses. Key findings include:
- **Mean Imputation:** Exhibited high MAE values (up to 700.71) and RMSE values (up to 5494.11), indicating poor performance.
- **Random Forest Imputation:** Showed the lowest MAE (as low as 108.31) and RMSE (as low as 1282.36), highlighting superior accuracy.
- **Volcano Plots:** Significant genes identified and visualized, with Random Forest imputation preserving a robust set of differentially expressed genes.

Future research directions include exploring more advanced imputation techniques, increasing sample sizes, improving computational resources, and developing efficient imputation algorithms.

### Conclusion
Handling missing data is vital for the integrity of biological analyses. This project showcases multiple imputation techniques, providing a robust framework for improving data quality and reliability in biological research. By addressing these challenges, researchers can enhance the accuracy and validity of their findings, contributing to a deeper understanding of complex biological systems and the genetic basis of various conditions.

For more details, refer to the scripts and visualizations provided in the repository.

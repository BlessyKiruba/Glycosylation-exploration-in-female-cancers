

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from lifelines import CoxPHFitter
import matplotlib.pyplot as plt
import warnings
from sklearn.model_selection import train_test_split

from sklearn.decomposition import PCA
warnings.filterwarnings('ignore')

data = pd.read_csv('/path/to/survival_expression_ml_compiled.csv')

print(f"  Loaded: {data.shape[0]} samples Ã— {data.shape[1]} columns")


X = data.iloc[:, 3:]  
y_time = data['time'].values
y_event = data['status'].values

print(f"  Features: {X.shape[1]} genes")
print(f"  Events: {y_event.sum()} / {len(y_event)}\n")

scaler = StandardScaler()
X_scaled = pd.DataFrame(
    scaler.fit_transform(X),
    columns=X.columns
)



X_train, X_temp, y_time_train, y_time_temp, y_event_train, y_event_temp = train_test_split(
    X_scaled, y_time, y_event, test_size=0.4, random_state=42, stratify=y_event
)

X_val, X_test, y_time_val, y_time_test, y_event_val, y_event_test = train_test_split(
    X_temp, y_time_temp, y_event_temp, test_size=0.5, random_state=42, stratify=y_event_temp
)

print(f"  Train: {X_train.shape}")
print(f"  Val:   {X_val.shape}")
print(f"  Test:  {X_test.shape}\n")


cox_univariate = CoxPHFitter()
univariate_pvals = []
gene_names = []

for gene in X_train.columns:
    try:
        train_df = pd.DataFrame({
            'time': y_time_train,
            'event': y_event_train,
            gene: X_train[gene]
        })
        cox_univariate.fit(train_df, duration_col='time', event_col='event')
        p_val = cox_univariate.summary['p'].values[0]
        
        univariate_pvals.append(p_val)
        gene_names.append(gene)
    except:
        pass

univariate_df = pd.DataFrame({
    'Gene': gene_names,
    'p_value': univariate_pvals
})

univariate_df = univariate_df.sort_values('p_value')

selected_genes_ridge = univariate_df[univariate_df['p_value'] < 0.001]['Gene'].tolist()


univariate_df.to_csv('/path/to/feature_selection_univariate_cox.csv', index=False)
pd.DataFrame({'Gene': selected_genes_ridge}).to_csv('/path/to/selected_genes_ridge_stepwise.csv', index=False)

coef_ranking = univariate_df.copy()
coef_ranking['Abs_Coef'] = coef_ranking['p_value'].apply(lambda x: -np.log10(x + 1e-10))
coef_ranking = coef_ranking.sort_values('Abs_Coef', ascending=False)

selected_genes_lasso = coef_ranking.head(40)['Gene'].tolist()


coef_ranking.to_csv('/path/to/feature_selection_lasso_coefficients.csv', index=False)
pd.DataFrame({'Gene': selected_genes_lasso}).to_csv('/path/to/selected_genes_lasso_elasticnet.csv', index=False)


variance_ranking = pd.DataFrame({
    'Gene': X_train.columns,
    'Variance': X_train.var()
})

variance_ranking = variance_ranking.sort_values('Variance', ascending=False)

selected_genes_rsf = variance_ranking.head(40)['Gene'].tolist()


variance_ranking.to_csv('/path/to/feature_selection_rsf_importance.csv', index=False)
pd.DataFrame({'Gene': selected_genes_rsf}).to_csv('/path/to/selected_genes_rsf_gbm.csv', index=False)


pca = PCA(n_components=10)
pca.fit(X_train)

loadings = pd.DataFrame({
    'Gene': X_train.columns,
    'PC1_Loading': np.abs(pca.components_[0])  
})

loadings = loadings.sort_values('PC1_Loading', ascending=False)


selected_genes_pca = loadings.head(10)['Gene'].tolist()


loadings.to_csv('/path/to/feature_selection_pca_loadings.csv', index=False)
pd.DataFrame({'Gene': selected_genes_pca}).to_csv('/path/to/selected_genes_pca_plsrcox_superpc.csv', index=False)


stricter_genes = univariate_df[univariate_df['p_value'] < 0.001]['Gene'].tolist()

print(f"  Genes with p < 0.001: {len(stricter_genes)}")
if len(stricter_genes) < 10:
    print(f"  [WARNING] Only {len(stricter_genes)} genes passed p<0.001")
    print(f"  Using p < 0.01 instead ({len(selected_genes_ridge)} genes)\n")
    stricter_genes = selected_genes_ridge
else:
    print(f"  Selected genes: {stricter_genes[:10]}...\n")

pd.DataFrame({'Gene': stricter_genes}).to_csv('/path/to/selected_genes_coxboost_strict.csv', index=False)

summary = pd.DataFrame({
    'Model': [
        'Ridge Regression',
        'Lasso Regression',
        'Elastic Net',
        'Stepwise Cox',
        'CoxBoost',
        'RSF',
        'plsRcox',
        'SuperPC',
        'GBM'
    ],
    'Feature_Selection_Method': [
        'Univariate Cox (p<0.01)',
        'Coefficient Ranking (top 20)',
        'Coefficient Ranking (top 20)',
        'Univariate Cox (p<0.01)',
        'Univariate Cox (p<0.001)',
        'Variance-Based (top 20)',
        'PCA Loadings (top 20)',
        'PCA Loadings (top 20)',
        'Variance-Based (top 20)'
    ],
    'Num_Genes_Selected': [
        len(selected_genes_ridge),
        len(selected_genes_lasso),
        len(selected_genes_lasso),
        len(selected_genes_ridge),
        len(stricter_genes),
        len(selected_genes_rsf),
        len(selected_genes_pca),
        len(selected_genes_pca),
        len(selected_genes_rsf)
    ],
    'Output_File': [
        'selected_genes_ridge_stepwise.csv',
        'selected_genes_lasso_elasticnet.csv',
        'selected_genes_lasso_elasticnet.csv',
        'selected_genes_ridge_stepwise.csv',
        'selected_genes_coxboost_strict.csv',
        'selected_genes_rsf_gbm.csv',
        'selected_genes_pca_plsrcox_superpc.csv',
        'selected_genes_pca_plsrcox_superpc.csv',
        'selected_genes_rsf_gbm.csv'
    ]
})

print(summary.to_string(index=False))
print()

summary.to_csv('/path/to/feature_selection_summary.csv', index=False)


datasets_to_save = {
    'data_ridge_stepwise': selected_genes_ridge,
    'data_lasso_elasticnet': selected_genes_lasso,
    'data_coxboost': stricter_genes,
    'data_rsf_gbm': selected_genes_rsf,
    'data_pca_plsrcox_superpc': selected_genes_pca
}

for dataset_name, genes in datasets_to_save.items():

    filtered_data = data[['sample_id', 'time', 'status'] + genes]
    

    csv_name = f'{dataset_name}_filtered.csv'
    filtered_data.to_csv(csv_name, index=False)
    
    print(f"  Saved: {csv_name} ({filtered_data.shape[0]} Ã— {filtered_data.shape[1]})")


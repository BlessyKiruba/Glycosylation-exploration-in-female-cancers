
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve
from lifelines.utils import concordance_index
from lifelines import KaplanMeierFitter
import warnings
warnings.filterwarnings('ignore')

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import torch.nn.functional as F

import optuna
from optuna.visualization import plot_optimization_history, plot_param_importances

import os
import json
from datetime import datetime
import argparse

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(RANDOM_SEED)


def load_and_preprocess_data(file_path):

    df = pd.read_csv(file_path)
    
    print(f"Dataset shape: {df.shape}")
    print(f"\nSurvival status distribution:\n{df['status'].value_counts()}")
    print(f"Time range: {df['time'].min()} - {df['time'].max()} days")
    

    sample_ids = df['sample_id'].values
    time = df['time'].values
    status = df['status'].values
    

    feature_cols = [col for col in df.columns if col not in ['sample_id', 'time', 'status']]
    X = df[feature_cols].values
    
    print(f"\nNumber of genes (features): {len(feature_cols)}")
    print(f"Number of samples: {len(X)}")
    
    return X, time, status, feature_cols, sample_ids


def create_time_bins(time, status, num_bins=10):

    event_times = time[status == 1]
    bins = np.percentile(event_times, np.linspace(0, 100, num_bins + 1))
    bins = np.unique(bins)
    bins[0] = 0
    bins[-1] = time.max() + 1
    return bins



class SurvivalDataset(Dataset):
    def __init__(self, X, time, status, bins=None):
        self.X = torch.FloatTensor(X)
        self.time = torch.FloatTensor(time)
        self.status = torch.LongTensor(status)
        self.bins = bins
        
        if bins is not None:
            self.time_bin = torch.LongTensor(np.digitize(time, bins) - 1)
            self.time_bin = torch.clamp(self.time_bin, 0, len(bins) - 2)
    
    def __len__(self):
        return len(self.X)
    
    def __getitem__(self, idx):
        if self.bins is not None:
            return self.X[idx], self.time[idx], self.status[idx], self.time_bin[idx]
        return self.X[idx], self.time[idx], self.status[idx]



class DeepSurv(nn.Module):

    def __init__(self, input_dim, hidden_dims=[128, 64, 32], dropout=0.3, 
                 use_batch_norm=True, activation='relu'):
        super(DeepSurv, self).__init__()
        
        layers = []
        prev_dim = input_dim
        
        if activation == 'relu':
            act_fn = nn.ReLU
        elif activation == 'elu':
            act_fn = nn.ELU
        elif activation == 'selu':
            act_fn = nn.SELU
        else:
            act_fn = nn.ReLU
        
        for hidden_dim in hidden_dims:
            layers.append(nn.Linear(prev_dim, hidden_dim))
            
            if use_batch_norm:
                layers.append(nn.BatchNorm1d(hidden_dim))
            
            layers.append(act_fn())
            layers.append(nn.Dropout(dropout))
            prev_dim = hidden_dim
        
        layers.append(nn.Linear(prev_dim, 1))
        
        self.network = nn.Sequential(*layers)
        
        self._initialize_weights()
    
    def _initialize_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.kaiming_normal_(m.weight)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
    
    def forward(self, x):
        return self.network(x)


def cox_loss(risk_scores, time, status):

    idx = torch.argsort(time, descending=True)
    risk_scores = risk_scores[idx]
    status = status[idx]
    
    hazard_ratio = torch.exp(risk_scores)
    log_risk = torch.log(torch.cumsum(hazard_ratio, dim=0) + 1e-7)
    
    uncensored_likelihood = risk_scores - log_risk
    loss = -torch.sum(uncensored_likelihood * status) / (torch.sum(status) + 1e-7)
    
    return loss



class DeepHit(nn.Module):

    def __init__(self, input_dim, num_bins, hidden_dims=[128, 64], dropout=0.3,
                 use_batch_norm=True, activation='relu'):
        super(DeepHit, self).__init__()
        
        self.num_bins = num_bins
        

        if activation == 'relu':
            act_fn = nn.ReLU
        elif activation == 'elu':
            act_fn = nn.ELU
        elif activation == 'selu':
            act_fn = nn.SELU
        else:
            act_fn = nn.ReLU
        
        shared_layers = []
        prev_dim = input_dim
        
        for hidden_dim in hidden_dims:
            shared_layers.append(nn.Linear(prev_dim, hidden_dim))
            
            if use_batch_norm:
                shared_layers.append(nn.BatchNorm1d(hidden_dim))
            
            shared_layers.append(act_fn())
            shared_layers.append(nn.Dropout(dropout))
            prev_dim = hidden_dim
        
        self.shared = nn.Sequential(*shared_layers)
        self.output = nn.Linear(prev_dim, num_bins)
        

        self._initialize_weights()
    
    def _initialize_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.kaiming_normal_(m.weight)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
    
    def forward(self, x):
        shared_repr = self.shared(x)
        pmf = torch.softmax(self.output(shared_repr), dim=1)
        return pmf


def deephit_loss(pmf, time_bin, status, alpha=0.5, beta=0.1):

    batch_size = pmf.shape[0]
    
    event_pmf = pmf[range(batch_size), time_bin]
    likelihood_loss = -torch.mean(status.float() * torch.log(event_pmf + 1e-7))
    
    cif = torch.cumsum(pmf, dim=1)
    
    ranking_loss = 0
    num_pairs = 0
    
    for i in range(batch_size):
        if status[i] == 1:
            for j in range(batch_size):
                if time_bin[j] > time_bin[i]:
                    eta = cif[j, time_bin[i]] - cif[i, time_bin[i]]
                    ranking_loss += torch.exp(eta)
                    num_pairs += 1
    
    if num_pairs > 0:
        ranking_loss = ranking_loss / num_pairs
    

    smooth_loss = torch.mean(torch.abs(pmf[:, 1:] - pmf[:, :-1]))
    

    total_loss = likelihood_loss + alpha * ranking_loss + beta * smooth_loss
    
    return total_loss




def train_deepsurv(model, train_loader, val_loader, num_epochs=100, lr=0.001, 
                   device='cpu', patience=20, use_multi_gpu=False):

    if use_multi_gpu and torch.cuda.is_available() and torch.cuda.device_count() > 1:
        print(f"Using {torch.cuda.device_count()} GPUs!")
        model = nn.DataParallel(model)
    
    model = model.to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', 
                                                      patience=10, factor=0.5, verbose=True)
    
    train_losses = []
    val_losses = []
    val_cindices = []
    best_val_cindex = 0
    best_model_state = None
    patience_counter = 0
    
    for epoch in range(num_epochs):

        model.train()
        train_loss = 0
        for X_batch, time_batch, status_batch in train_loader:
            X_batch = X_batch.to(device)
            time_batch = time_batch.to(device)
            status_batch = status_batch.to(device)
            
            optimizer.zero_grad()
            risk_scores = model(X_batch).squeeze()
            loss = cox_loss(risk_scores, time_batch, status_batch)
            loss.backward()
            

            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            
            optimizer.step()
            train_loss += loss.item()
        
        train_loss /= len(train_loader)
        train_losses.append(train_loss)
        
 
        model.eval()
        val_loss = 0
        all_risks = []
        all_times = []
        all_status = []
        
        with torch.no_grad():
            for X_batch, time_batch, status_batch in val_loader:
                X_batch = X_batch.to(device)
                time_batch = time_batch.to(device)
                status_batch = status_batch.to(device)
                
                risk_scores = model(X_batch).squeeze()
                loss = cox_loss(risk_scores, time_batch, status_batch)
                val_loss += loss.item()
                
                all_risks.extend(risk_scores.cpu().numpy())
                all_times.extend(time_batch.cpu().numpy())
                all_status.extend(status_batch.cpu().numpy())
        
        val_loss /= len(val_loader)
        val_losses.append(val_loss)
        
        cindex = concordance_index(all_times, -np.array(all_risks), all_status)
        val_cindices.append(cindex)
        
        scheduler.step(val_loss)
        

        if cindex > best_val_cindex:
            best_val_cindex = cindex
            best_model_state = model.state_dict().copy()
            patience_counter = 0
        else:
            patience_counter += 1
        
        if patience_counter >= patience:
            print(f"Early stopping at epoch {epoch+1}")
            break
        
        if (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch+1}/{num_epochs} - Train Loss: {train_loss:.4f}, "
                  f"Val Loss: {val_loss:.4f}, Val C-index: {cindex:.4f}")
    

    model.load_state_dict(best_model_state)
    
    return model, train_losses, val_losses, val_cindices


def train_deephit(model, train_loader, val_loader, num_epochs=100, lr=0.001,
                  alpha=0.5, beta=0.1, device='cpu', patience=20, use_multi_gpu=False):

    if use_multi_gpu and torch.cuda.is_available() and torch.cuda.device_count() > 1:
        print(f"Using {torch.cuda.device_count()} GPUs!")
        model = nn.DataParallel(model)
    
    model = model.to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min',
                                                      patience=10, factor=0.5, verbose=True)
    
    train_losses = []
    val_losses = []
    val_cindices = []
    best_val_cindex = 0
    best_model_state = None
    patience_counter = 0
    
    for epoch in range(num_epochs):

        model.train()
        train_loss = 0
        for X_batch, time_batch, status_batch, time_bin_batch in train_loader:
            X_batch = X_batch.to(device)
            time_batch = time_batch.to(device)
            status_batch = status_batch.to(device)
            time_bin_batch = time_bin_batch.to(device)
            
            optimizer.zero_grad()
            pmf = model(X_batch)
            loss = deephit_loss(pmf, time_bin_batch, status_batch, alpha=alpha, beta=beta)
            loss.backward()
            

            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            
            optimizer.step()
            train_loss += loss.item()
        
        train_loss /= len(train_loader)
        train_losses.append(train_loss)
        
        model.eval()
        val_loss = 0
        all_risks = []
        all_times = []
        all_status = []
        
        with torch.no_grad():
            for X_batch, time_batch, status_batch, time_bin_batch in val_loader:
                X_batch = X_batch.to(device)
                time_batch = time_batch.to(device)
                status_batch = status_batch.to(device)
                time_bin_batch = time_bin_batch.to(device)
                
                pmf = model(X_batch)
                loss = deephit_loss(pmf, time_bin_batch, status_batch, alpha=alpha, beta=beta)
                val_loss += loss.item()
                
                cif = torch.cumsum(pmf, dim=1)
                risk_scores = torch.sum(cif, dim=1)
                
                all_risks.extend(risk_scores.cpu().numpy())
                all_times.extend(time_batch.cpu().numpy())
                all_status.extend(status_batch.cpu().numpy())
        
        val_loss /= len(val_loader)
        val_losses.append(val_loss)
        

        cindex = concordance_index(all_times, all_risks, all_status)
        val_cindices.append(cindex)
        
        scheduler.step(val_loss)
        

        if cindex > best_val_cindex:
            best_val_cindex = cindex
            best_model_state = model.state_dict().copy()
            patience_counter = 0
        else:
            patience_counter += 1
        
        if patience_counter >= patience:
            print(f"Early stopping at epoch {epoch+1}")
            break
        
        if (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch+1}/{num_epochs} - Train Loss: {train_loss:.4f}, "
                  f"Val Loss: {val_loss:.4f}, Val C-index: {cindex:.4f}")
    

    model.load_state_dict(best_model_state)
    
    return model, train_losses, val_losses, val_cindices




def objective_deepsurv(trial, X_train, time_train, status_train, X_val, time_val, status_val, device):

    n_layers = trial.suggest_int('n_layers', 2, 4)
    hidden_dims = [trial.suggest_int(f'n_units_l{i}', 32, 256, log=True) for i in range(n_layers)]
    dropout = trial.suggest_float('dropout', 0.1, 0.5)
    lr = trial.suggest_float('lr', 1e-4, 1e-2, log=True)
    batch_size = trial.suggest_categorical('batch_size', [32, 64, 128])
    activation = trial.suggest_categorical('activation', ['relu', 'elu', 'selu'])
    

    train_dataset = SurvivalDataset(X_train, time_train, status_train)
    val_dataset = SurvivalDataset(X_val, time_val, status_val)
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
    

    model = DeepSurv(
        input_dim=X_train.shape[1],
        hidden_dims=hidden_dims,
        dropout=dropout,
        activation=activation
    )
    

    model, _, _, val_cindices = train_deepsurv(
        model, train_loader, val_loader,
        num_epochs=50, 
        lr=lr,
        device=device,
        patience=10
    )
    

    return max(val_cindices)


def objective_deephit(trial, X_train, time_train, status_train, X_val, time_val, 
                     status_val, time_bins, device):

    n_layers = trial.suggest_int('n_layers', 2, 4)
    hidden_dims = [trial.suggest_int(f'n_units_l{i}', 32, 256, log=True) for i in range(n_layers)]
    dropout = trial.suggest_float('dropout', 0.1, 0.5)
    lr = trial.suggest_float('lr', 1e-4, 1e-2, log=True)
    batch_size = trial.suggest_categorical('batch_size', [32, 64, 128])
    alpha = trial.suggest_float('alpha', 0.1, 1.0)
    beta = trial.suggest_float('beta', 0.0, 0.3)
    activation = trial.suggest_categorical('activation', ['relu', 'elu', 'selu'])
    

    train_dataset = SurvivalDataset(X_train, time_train, status_train, bins=time_bins)
    val_dataset = SurvivalDataset(X_val, time_val, status_val, bins=time_bins)
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)

    model = DeepHit(
        input_dim=X_train.shape[1],
        num_bins=len(time_bins)-1,
        hidden_dims=hidden_dims,
        dropout=dropout,
        activation=activation
    )
    

    model, _, _, val_cindices = train_deephit(
        model, train_loader, val_loader,
        num_epochs=50,  
        lr=lr,
        alpha=alpha,
        beta=beta,
        device=device,
        patience=10
    )
    

    return max(val_cindices)


def tune_hyperparameters(model_name, X_train, time_train, status_train, 
                        X_val, time_val, status_val, time_bins=None,
                        n_trials=50, device='cpu', save_dir=None):


    study = optuna.create_study(direction='maximize', study_name=f'{model_name}_tuning')
    

    if model_name == 'DeepSurv':
        objective = lambda trial: objective_deepsurv(
            trial, X_train, time_train, status_train, X_val, time_val, status_val, device
        )
    else:  
        objective = lambda trial: objective_deephit(
            trial, X_train, time_train, status_train, X_val, time_val, status_val, time_bins, device
        )
    

    study.optimize(objective, n_trials=n_trials, show_progress_bar=True)
    

    print(f"\nBest trial:")
    print(f"  C-index: {study.best_trial.value:.4f}")
    print(f"  Params:")
    for key, value in study.best_trial.params.items():
        print(f"    {key}: {value}")
    

    if save_dir:

        with open(f'{save_dir}/{model_name}_best_params.json', 'w') as f:
            json.dump(study.best_trial.params, f, indent=4)

        try:
            fig1 = plot_optimization_history(study)
            fig1.write_image(f'{save_dir}/{model_name}_optimization_history.png')
            
            fig2 = plot_param_importances(study)
            fig2.write_image(f'{save_dir}/{model_name}_param_importances.png')
        except:
            print("Could not save Optuna visualizations (requires kaleido package)")
    
    return study.best_trial.params



def evaluate_model(model, test_loader, is_deephit=False, device='cpu'):

    model.eval()
    all_risks = []
    all_times = []
    all_status = []
    
    with torch.no_grad():
        if is_deephit:
            for X_batch, time_batch, status_batch, time_bin_batch in test_loader:
                X_batch = X_batch.to(device)
                pmf = model(X_batch)
                cif = torch.cumsum(pmf, dim=1)
                risk_scores = torch.sum(cif, dim=1)
                
                all_risks.extend(risk_scores.cpu().numpy())
                all_times.extend(time_batch.cpu().numpy())
                all_status.extend(status_batch.cpu().numpy())
        else:
            for X_batch, time_batch, status_batch in test_loader:
                X_batch = X_batch.to(device)
                risk_scores = model(X_batch).squeeze()
                
                all_risks.extend(-risk_scores.cpu().numpy())
                all_times.extend(time_batch.cpu().numpy())
                all_status.extend(status_batch.cpu().numpy())
    
    all_risks = np.array(all_risks)
    all_times = np.array(all_times)
    all_status = np.array(all_status)
    
    cindex = concordance_index(all_times, all_risks, all_status)
    

    median_time = np.median(all_times[all_status == 1])
    risk_at_median = (all_times > median_time).astype(int)
    
    if len(np.unique(risk_at_median)) > 1:
        auc = roc_auc_score(risk_at_median, all_risks)
    else:
        auc = None
    
    return cindex, auc, all_risks, all_times, all_status


def calculate_feature_importance(model, X, feature_names, is_deephit=False, device='cpu'):

    model.eval()
    X_tensor = torch.FloatTensor(X).to(device)
    
    with torch.no_grad():
        if is_deephit:
            pmf = model(X_tensor)
            cif = torch.cumsum(pmf, dim=1)
            baseline_pred = torch.sum(cif, dim=1).cpu().numpy()
        else:
            baseline_pred = model(X_tensor).squeeze().cpu().numpy()
    
    importances = []
    
    for i, feature in enumerate(feature_names):
        X_permuted = X.copy()
        np.random.shuffle(X_permuted[:, i])
        X_permuted_tensor = torch.FloatTensor(X_permuted).to(device)
        
        with torch.no_grad():
            if is_deephit:
                pmf = model(X_permuted_tensor)
                cif = torch.cumsum(pmf, dim=1)
                permuted_pred = torch.sum(cif, dim=1).cpu().numpy()
            else:
                permuted_pred = model(X_permuted_tensor).squeeze().cpu().numpy()
        
        importance = np.abs(baseline_pred - permuted_pred).mean()
        importances.append(importance)
    
    importances = np.array(importances)
    importances = importances / (importances.sum() + 1e-7)
    
    importance_df = pd.DataFrame({
        'Gene': feature_names,
        'Importance': importances
    }).sort_values('Importance', ascending=False)
    
    return importance_df



def plot_training_curves(train_losses, val_losses, val_cindices, model_name, save_dir):

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    axes[0].plot(train_losses, label='Training Loss', linewidth=2)
    axes[0].plot(val_losses, label='Validation Loss', linewidth=2)
    axes[0].set_xlabel('Epoch', fontsize=12)
    axes[0].set_ylabel('Loss', fontsize=12)
    axes[0].set_title(f'{model_name} - Training Curves', fontsize=14)
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    axes[1].plot(val_cindices, label='Validation C-index', linewidth=2, color='green')
    axes[1].axhline(y=0.5, color='r', linestyle='--', label='Random', alpha=0.5)
    axes[1].set_xlabel('Epoch', fontsize=12)
    axes[1].set_ylabel('C-index', fontsize=12)
    axes[1].set_title(f'{model_name} - C-index', fontsize=14)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{save_dir}/{model_name}_training_curves.png', dpi=300, bbox_inches='tight')
    plt.close()


def plot_risk_stratification(risks, times, status, model_name, save_dir, n_groups=3):

    risk_quantiles = np.percentile(risks, np.linspace(0, 100, n_groups + 1))
    risk_groups = np.digitize(risks, risk_quantiles[1:-1])
    
    plt.figure(figsize=(10, 6))
    kmf = KaplanMeierFitter()
    
    colors = ['green', 'orange', 'red']
    labels = ['Low Risk', 'Medium Risk', 'High Risk']
    
    for i in range(n_groups):
        mask = risk_groups == i
        if mask.sum() > 0:
            kmf.fit(times[mask], status[mask], label=labels[i])
            kmf.plot_survival_function(ci_show=True, color=colors[i], linewidth=2)
    
    plt.xlabel('Time (days)', fontsize=12)
    plt.ylabel('Survival Probability', fontsize=12)
    plt.title(f'{model_name} - Risk Stratification (Kaplan-Meier)', fontsize=14)
    plt.legend(loc='best')
    plt.grid(True, alpha=0.3)
    
    plt.savefig(f'{save_dir}/{model_name}_risk_stratification.png', dpi=300, bbox_inches='tight')
    plt.close()


def plot_feature_importance(importance_df, model_name, save_dir, top_n=30):

    top_features = importance_df.head(top_n)
    
    plt.figure(figsize=(10, 8))
    plt.barh(range(len(top_features)), top_features['Importance'].values, color='steelblue')
    plt.yticks(range(len(top_features)), top_features['Gene'].values)
    plt.xlabel('Importance Score', fontsize=12)
    plt.ylabel('Gene', fontsize=12)
    plt.title(f'{model_name} - Top {top_n} Prognostic Genes', fontsize=14)
    plt.gca().invert_yaxis()
    plt.grid(True, alpha=0.3, axis='x')
    plt.tight_layout()
    
    plt.savefig(f'{save_dir}/{model_name}_feature_importance.png', dpi=300, bbox_inches='tight')
    plt.close()



def main():
    parser = argparse.ArgumentParser(description='Pan-Cancer Survival Analysis')
    parser.add_argument('--tune', action='store_true', help='Enable hyperparameter tuning')
    parser.add_argument('--n_trials', type=int, default=30, help='Number of tuning trials')
    parser.add_argument('--multi_gpu', action='store_true', help='Use multiple GPUs if available')
    parser.add_argument('--epochs', type=int, default=150, help='Maximum training epochs')
    args = parser.parse_args()
    

    DATA_FILE = '/path/to/survival_expression_ml_compiled.csv'
    OUTPUT_DIR = '/path/to/outputs'
    RESULTS_DIR = os.path.join(OUTPUT_DIR, 'survival_results_tuned')
    os.makedirs(RESULTS_DIR, exist_ok=True)
    

    if torch.cuda.is_available():
        device = torch.device('cuda')
        print(f"Using device: CUDA GPU")
        print(f"GPU: {torch.cuda.get_device_name(0)}")
        print(f"Number of GPUs: {torch.cuda.device_count()}")
    elif torch.backends.mps.is_available():
        device = torch.device('mps')
        print(f"Using device: MPS (Apple Silicon)")
        print(f"Metal GPU acceleration enabled")
    else:
        device = torch.device('cpu')
        print(f"Using device: CPU")
    
    print(f"Device: {device}")
    

    X, time, status, feature_names, sample_ids = load_and_preprocess_data(DATA_FILE)
    

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    

    TEST_SIZE = 0.2
    VAL_SIZE = 0.2
    
    X_temp, X_test, time_temp, time_test, status_temp, status_test = train_test_split(
        X_scaled, time, status, test_size=TEST_SIZE, random_state=RANDOM_SEED, stratify=status
    )
    
    X_train, X_val, time_train, time_val, status_train, status_val = train_test_split(
        X_temp, time_temp, status_temp, test_size=VAL_SIZE/(1-TEST_SIZE), 
        random_state=RANDOM_SEED, stratify=status_temp
    )
    
    print(f"\nData split:")
    print(f"Train: {len(X_train)} samples")
    print(f"Validation: {len(X_val)} samples")
    print(f"Test: {len(X_test)} samples")
    

    NUM_TIME_BINS = 15  
    time_bins = create_time_bins(time_train, status_train, num_bins=NUM_TIME_BINS)

    if args.tune:
        ds_best_params = tune_hyperparameters(
            'DeepSurv', X_train, time_train, status_train,
            X_val, time_val, status_val,
            n_trials=args.n_trials,
            device=device,
            save_dir=RESULTS_DIR
        )
        

        dh_best_params = tune_hyperparameters(
            'DeepHit', X_train, time_train, status_train,
            X_val, time_val, status_val,
            time_bins=time_bins,
            n_trials=args.n_trials,
            device=device,
            save_dir=RESULTS_DIR
        )
    else:

        ds_best_params = {
            'n_layers': 3,
            'n_units_l0': 256,
            'n_units_l1': 128,
            'n_units_l2': 64,
            'dropout': 0.4,
            'lr': 0.0005,
            'batch_size': 64,
            'activation': 'relu'
        }
        
        dh_best_params = {
            'n_layers': 3,
            'n_units_l0': 256,
            'n_units_l1': 128,
            'n_units_l2': 64,
            'dropout': 0.4,
            'lr': 0.001,
            'batch_size': 64,
            'alpha': 0.3,
            'beta': 0.1,
            'activation': 'relu'
        }


    n_layers = ds_best_params['n_layers']
    hidden_dims = [ds_best_params[f'n_units_l{i}'] for i in range(n_layers)]
    

    train_dataset_ds = SurvivalDataset(X_train, time_train, status_train)
    val_dataset_ds = SurvivalDataset(X_val, time_val, status_val)
    test_dataset_ds = SurvivalDataset(X_test, time_test, status_test)
    
    train_loader_ds = DataLoader(train_dataset_ds, batch_size=ds_best_params['batch_size'], shuffle=True)
    val_loader_ds = DataLoader(val_dataset_ds, batch_size=ds_best_params['batch_size'], shuffle=False)
    test_loader_ds = DataLoader(test_dataset_ds, batch_size=ds_best_params['batch_size'], shuffle=False)
    

    deepsurv_model = DeepSurv(
        input_dim=X_train.shape[1],
        hidden_dims=hidden_dims,
        dropout=ds_best_params['dropout'],
        activation=ds_best_params['activation']
    )
    

    deepsurv_model, ds_train_losses, ds_val_losses, ds_val_cindices = train_deepsurv(
        deepsurv_model, train_loader_ds, val_loader_ds,
        num_epochs=args.epochs,
        lr=ds_best_params['lr'],
        device=device,
        patience=30,
        use_multi_gpu=args.multi_gpu
    )

    ds_cindex, ds_auc, ds_risks, ds_times, ds_status = evaluate_model(
        deepsurv_model, test_loader_ds, is_deephit=False, device=device
    )
    
    print(f"\nDeepSurv Test Results:")
    print(f"C-index: {ds_cindex:.4f}")
    if ds_auc is not None:
        print(f"AUC: {ds_auc:.4f}")

    print("\nCalculating DeepSurv feature importance...")
    ds_importance = calculate_feature_importance(
        deepsurv_model, X_test, feature_names, is_deephit=False, device=device
    )
    

    n_layers = dh_best_params['n_layers']
    hidden_dims = [dh_best_params[f'n_units_l{i}'] for i in range(n_layers)]
    

    train_dataset_dh = SurvivalDataset(X_train, time_train, status_train, bins=time_bins)
    val_dataset_dh = SurvivalDataset(X_val, time_val, status_val, bins=time_bins)
    test_dataset_dh = SurvivalDataset(X_test, time_test, status_test, bins=time_bins)
    
    train_loader_dh = DataLoader(train_dataset_dh, batch_size=dh_best_params['batch_size'], shuffle=True)
    val_loader_dh = DataLoader(val_dataset_dh, batch_size=dh_best_params['batch_size'], shuffle=False)
    test_loader_dh = DataLoader(test_dataset_dh, batch_size=dh_best_params['batch_size'], shuffle=False)
    

    deephit_model = DeepHit(
        input_dim=X_train.shape[1],
        num_bins=len(time_bins)-1,
        hidden_dims=hidden_dims,
        dropout=dh_best_params['dropout'],
        activation=dh_best_params['activation']
    )
    

    deephit_model, dh_train_losses, dh_val_losses, dh_val_cindices = train_deephit(
        deephit_model, train_loader_dh, val_loader_dh,
        num_epochs=args.epochs,
        lr=dh_best_params['lr'],
        alpha=dh_best_params['alpha'],
        beta=dh_best_params['beta'],
        device=device,
        patience=30,
        use_multi_gpu=args.multi_gpu
    )

    dh_cindex, dh_auc, dh_risks, dh_times, dh_status = evaluate_model(
        deephit_model, test_loader_dh, is_deephit=True, device=device
    )
    
    print(f"\nDeepHit Test Results:")
    print(f"C-index: {dh_cindex:.4f}")
    if dh_auc is not None:
        print(f"AUC: {dh_auc:.4f}")
    

    print("\nCalculating DeepHit feature importance...")
    dh_importance = calculate_feature_importance(
        deephit_model, X_test, feature_names, is_deephit=True, device=device
    )
    

    torch.save(deepsurv_model.state_dict(), f'{RESULTS_DIR}/deepsurv_model.pth')
    torch.save(deephit_model.state_dict(), f'{RESULTS_DIR}/deephit_model.pth')
    

    ds_importance.to_csv(f'{RESULTS_DIR}/deepsurv_feature_importance.csv', index=False)
    dh_importance.to_csv(f'{RESULTS_DIR}/deephit_feature_importance.csv', index=False)
    

    with open(f'{RESULTS_DIR}/deepsurv_hyperparams.json', 'w') as f:
        json.dump(ds_best_params, f, indent=4)
    with open(f'{RESULTS_DIR}/deephit_hyperparams.json', 'w') as f:
        json.dump(dh_best_params, f, indent=4)

    results = {
        'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'hyperparameter_tuning': args.tune,
        'multi_gpu': args.multi_gpu,
        'dataset': {
            'total_samples': len(X),
            'n_features': len(feature_names),
            'train_samples': len(X_train),
            'val_samples': len(X_val),
            'test_samples': len(X_test),
            'event_rate': float(status.mean())
        },
        'deepsurv': {
            'test_cindex': float(ds_cindex),
            'test_auc': float(ds_auc) if ds_auc is not None else None,
            'best_val_cindex': float(max(ds_val_cindices)),
            'hyperparameters': ds_best_params,
            'top_10_genes': ds_importance.head(10)['Gene'].tolist()
        },
        'deephit': {
            'test_cindex': float(dh_cindex),
            'test_auc': float(dh_auc) if dh_auc is not None else None,
            'best_val_cindex': float(max(dh_val_cindices)),
            'hyperparameters': dh_best_params,
            'top_10_genes': dh_importance.head(10)['Gene'].tolist()
        }
    }
    
    with open(f'{RESULTS_DIR}/results_summary.json', 'w') as f:
        json.dump(results, f, indent=4)
    

    plot_training_curves(ds_train_losses, ds_val_losses, ds_val_cindices, 'DeepSurv', RESULTS_DIR)
    plot_training_curves(dh_train_losses, dh_val_losses, dh_val_cindices, 'DeepHit', RESULTS_DIR)
    plot_risk_stratification(ds_risks, ds_times, ds_status, 'DeepSurv', RESULTS_DIR)
    plot_risk_stratification(dh_risks, dh_times, dh_status, 'DeepHit', RESULTS_DIR)
    plot_feature_importance(ds_importance, 'DeepSurv', RESULTS_DIR, top_n=30)
    plot_feature_importance(dh_importance, 'DeepHit', RESULTS_DIR, top_n=30)
    

    print(f"\nDeepSurv (Optimized):")
    print(f"  Test C-index: {ds_cindex:.4f}")
    print(f"  Architecture: {hidden_dims}")
    print(f"  Top 5 prognostic genes:")
    for i, row in ds_importance.head(5).iterrows():
        print(f"    {row['Gene']}: {row['Importance']:.4f}")
    
    print(f"\nDeepHit (Optimized):")
    print(f"  Test C-index: {dh_cindex:.4f}")
    print(f"  Architecture: {hidden_dims}")
    print(f"  Top 5 prognostic genes:")
    for i, row in dh_importance.head(5).iterrows():
        print(f"    {row['Gene']}: {row['Importance']:.4f}")
    
    print(f"\nAll results saved to: {RESULTS_DIR}")
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
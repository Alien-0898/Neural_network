import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler,MinMaxScaler
import numpy as np
from sklearn.metrics import mean_squared_error
import math
# ==== Preprocessing ====
def preprocess(X, y):
    x_scaler = StandardScaler()
    y_scaler = StandardScaler()
    X_scaled = x_scaler.fit_transform(X)
    y_scaled = y_scaler.fit_transform(y)
    return X_scaled, y_scaled, x_scaler, y_scaler

# ==== Network ====
class PESNet(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 80),
            nn.Tanh(),
            nn.Dropout(0.1),  # Optional: better generalization
            nn.Linear(80, 80),
            nn.Tanh(),
            nn.Linear(80, 20),
            nn.Tanh(),
            nn.Linear(20, 1)
        )

    def forward(self, x):
        return self.net(x)

# ==== Training ====
def train_nn(X, y, max_epochs=1000, patience=25):
    print(X.shape)
    X, y, x_scaler, y_scaler = preprocess(X, y)

    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.1, random_state=42)
    X_train_t = torch.tensor(X_train, dtype=torch.float32)
    y_train_t = torch.tensor(y_train, dtype=torch.float32)
    X_val_t = torch.tensor(X_val, dtype=torch.float32)
    y_val_t = torch.tensor(y_val, dtype=torch.float32)

    model = PESNet(X.shape[1])
    criterion = nn.SmoothL1Loss()
    optimizer = optim.AdamW(model.parameters(), lr=1e-2, weight_decay=1e-4)

    best_model = None
    best_loss = float('inf')
    patience_counter = 0

    for epoch in range(max_epochs):
        model.train()
        optimizer.zero_grad()
        output_train = model(X_train_t)
        loss = criterion(output_train, y_train_t)
        loss.backward()
        optimizer.step()

        model.eval()
        with torch.no_grad():
            output_val = model(X_val_t)
            val_loss = criterion(output_val, y_val_t)

            # Inverse transform to original scale for RMSD calculation
            train_pred = y_scaler.inverse_transform(output_train.numpy())
            val_pred = y_scaler.inverse_transform(output_val.numpy())
            y_train_true = y_scaler.inverse_transform(y_train_t.numpy())
            y_val_true = y_scaler.inverse_transform(y_val_t.numpy())

            train_rmsd = np.sqrt(mean_squared_error(y_train_true, train_pred))
            val_rmsd = np.sqrt(mean_squared_error(y_val_true, val_pred))

        print(f"Epoch {epoch:03d}: "
              f"Train Loss = {loss.item():.6f}, Train RMSD = {train_rmsd:.4f} | "
              f"Val Loss = {val_loss.item():.6f}, Val RMSD = {val_rmsd:.4f}")

        if val_loss < best_loss:
            best_loss = val_loss
            best_model = PESNet(X.shape[1])
            best_model.load_state_dict(model.state_dict())
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= patience:
                print("Early stopping.")
                break

    if best_model is None:
        best_model = model

    return best_model, x_scaler, y_scaler

# ==== LBFGS Fine-Tuning (Optional) ====
def finetune_lbfgs(model, X, y, x_scaler, y_scaler):
    X_scaled = torch.tensor(x_scaler.transform(X), dtype=torch.float32)
    y_scaled = torch.tensor(y_scaler.transform(y), dtype=torch.float32)

    criterion = nn.SmoothL1Loss()
    optimizer = optim.LBFGS(model.parameters(), lr=0.8, max_iter=2000, max_eval=None, tolerance_grad=1e-8, tolerance_change=1e-12, history_size=100)

    def closure():
        optimizer.zero_grad()
        output = model(X_scaled)
        loss = criterion(output, y_scaled)
        loss.backward()
        return loss

    print("Starting LBFGS fine-tuning...")
    model.train()
    optimizer.step(closure)
    print("LBFGS fine-tuning complete.")
    return model

# ==== Evaluation ====
def evaluate(model, X, y, x_scaler, y_scaler):
    X_scaled = torch.tensor(x_scaler.transform(X), dtype=torch.float32)
    with torch.no_grad():
        pred_scaled = model(X_scaled)
        pred = y_scaler.inverse_transform(pred_scaled.numpy())
        rmsd = np.sqrt(np.mean((pred - y) ** 2))
        print(f"RMSD = {rmsd:.5f}")
    return rmsd




def compute_G_matrix(X):
    assert X.shape[1] == 15, "Input must have 15 columns (x1 to x15)"
    x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15 = X.T

    G = np.empty((X.shape[0], 37))

    G[:, 0]  = x15
    G[:, 1]  = x14
    G[:, 2]  = x10
    G[:, 3]  = x11 + x12 + x13
    G[:, 4]  = x7 + x8 + x9
    G[:, 5]  = x4 + x5 + x6
    G[:, 6]  = x1 + x2 + x3
    G[:, 7]  = (x11**2 + x12**2 + x13**2)**0.5
    G[:, 8]  = (x7*x11 + x8*x12 + x9*x13)**0.5
    G[:, 9]  = (x4*x11 + x5*x12 + x6*x13)**(1/2)
    G[:,10]  = (x1*x11 + x2*x11 + x1*x12 + x3*x12 + x2*x13 + x3*x13)**(1/2)
    G[:,11]  = (x7**2 + x8**2 + x9**2)**(1/2)
    G[:,12]  = (x4*x7 + x5*x8 + x6*x9)**(1/2)
    G[:,13]  = (x1*x7 + x2*x7 + x1*x8 + x3*x8 + x2*x9 + x3*x9)**(1/2)
    G[:,14]  = (x4**2 + x5**2 + x6**2)**(1/2)
    G[:,15]  = (x1*x4 + x2*x4 + x1*x5 + x3*x5 + x2*x6 + x3*x6)**(1/2)
    G[:,16]  = (x1**2 + x2**2 + x3**2)**(1/2)
    G[:,17]  = (x11**3 + x12**3 + x13**3)**(1/3)
    G[:,18]  = (x7*x11**2 + x8*x12**2 + x9*x13**2)**(1/3)
    G[:,19]  = (x4*x11**2 + x5*x12**2 + x6*x13**2)**(1/3)
    G[:,20]  = (x1*x11**2 + x2*x11**2 + x1*x12**2 + x3*x12**2 + x2*x13**2 + x3*x13**2)**(1/3)
    G[:,21]  = (x7**2*x11 + x8**2*x12 + x9**2*x13)**(1/3)
    G[:,22]  = (x4*x7*x11 + x5*x8*x12 + x6*x9*x13)**(1/3)
    G[:,23]  = (x1*x7*x11 + x2*x7*x11 + x1*x8*x12 + x3*x8*x12 + x2*x9*x13 + x3*x9*x13)**(1/3)
    G[:,24]  = (x4**2*x11 + x5**2*x12 + x6**2*x13)**(1/3)
    G[:,25]  = (x1*x4*x11 + x2*x4*x11 + x1*x5*x12 + x3*x5*x12 + x2*x6*x13 + x3*x6*x13)**(1/3)
    G[:,26]  = (x1**2*x11 + x2**2*x11 + x1**2*x12 + x3**2*x12 + x2**2*x13 + x3**2*x13)**(1/3)
    G[:,27]  = (x7**3 + x8**3 + x9**3)**(1/3)
    G[:,28]  = (x4*x7**2 + x5*x8**2 + x6*x9**2)**(1/3)
    G[:,29]  = (x1*x7**2 + x2*x7**2 + x1*x8**2 + x3*x8**2 + x2*x9**2 + x3*x9**2)**(1/3)
    G[:,30]  = (x4**2*x7 + x5**2*x8 + x6**2*x9)**(1/3)
    G[:,31]  = (x1*x4*x7 + x2*x4*x7 + x1*x5*x8 + x3*x5*x8 + x2*x6*x9 + x3*x6*x9)**(1/3)
    G[:,32]  = (x1**2*x7 + x2**2*x7 + x1**2*x8 + x3**2*x8 + x2**2*x9 + x3**2*x9)**(1/3)
    G[:,33]  = (x4**3 + x5**3 + x6**3)**(1/3)
    G[:,34]  = (x1*x4**2 + x2*x4**2 + x1*x5**2 + x3*x5**2 + x2*x6**2 + x3*x6**2)**(1/3)
    G[:,35]  = (x1**2*x4 + x2**2*x4 + x1**2*x5 + x3**2*x5 + x2**2*x6 + x3**2*x6)**(1/3)
    G[:,36]  = (x1**3 + x2**3 + x3**3)**(1/3)

    return G


def cosine_cutoff_np(R, Rc):
    """Vectorized cosine cutoff function for NumPy arrays"""
    fc = 0.5 * (np.cos(np.pi * R / Rc) + 1.0)
    fc[R >= Rc] = 0.0
    return fc

# # Example input: shape (300, 15)
# # Replace this with your actual data
# X = np.random.uniform(0.5, 7.0, size=(300, 15))  # sample distances

# Rc = 6.0  # cutoff radius
# fc = cosine_cutoff_np(X, Rc)  # shape (300, 15)

# # Now compute exp(-X) * fc(X)
# result = np.exp(-X) * fc  # shape (300, 15)

# print(result.shape)       # (300, 15)


hatree_to_cm=219474.6
# X = np.array([...])  # Shape: (60000, D)
# y = np.array([...])  # Shape: (60000, 1)
data=np.loadtxt('/home/Tapish/neural_network/2.0.1/pes_data_active_learn')

# X=np.exp(-data[:,:-1])
X=np.copy(data[:,:-1])
fc = cosine_cutoff_np(X, 10)
X = np.exp(-0.5*X) * fc
# X=np.copy(1/data[:,:-1])
# X=X*np.exp(-0.5*X)
y=np.copy(hatree_to_cm*data[:,15])
y=y.reshape(-1, 1)
X=compute_G_matrix(X)
model, x_scaler, y_scaler = train_nn(X, y)
evaluate(model, X, y, x_scaler, y_scaler)
model = finetune_lbfgs(model, X, y, x_scaler, y_scaler)
evaluate(model, X, y, x_scaler, y_scaler)
test_data=np.loadtxt('/home/Tapish/neural_network/deleted_points_faiss_exp_015.txt')
# X=np.copy(1/test_data[:,:-1])
y=np.copy(hatree_to_cm*test_data[:,15])
y=y.reshape(-1, 1)
X=np.copy(test_data[:,:-1])
# X=X*np.exp(-0.5*X)
fc = cosine_cutoff_np(X, 10)
X = np.exp(-0.5*X) * fc
X=compute_G_matrix(X)
print('test--------')
evaluate(model, X, y, x_scaler, y_scaler)

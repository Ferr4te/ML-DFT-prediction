from sklearn.svm import LinearSVR, SVR
from sklearn.linear_model import SGDRegressor, Ridge
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import cross_validate,  GridSearchCV, train_test_split, ShuffleSplit
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit import Chem
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np
import tqdm 

def get_data_fingerprint(column, bits, fingerprint):
    excel_file_path = "dataset_with_redox_HOMO_LUMO.xlsx"
    database = pd.read_excel(excel_file_path)
    print(database.head())
    database.columns = database.columns.str.strip()  # remove empty row

    condition = database["redox_potential_V"].notna() & np.isfinite(database["redox_potential_V"])
    database_sel = database.loc[condition].copy()

    redox = database_sel["redox_potential_V"]               # Series
    gap = database_sel["HOMO_LUMO_gap_Eh"].to_numpy()       # 1D np.array (n_samples,)
    smiles = database_sel["SMILES"]

    # Plot data distribution graph
    plt.figure(figsize=(8, 4))
    plt.hist(redox, bins=50)
    plt.xlabel("Redox potential (V)")
    plt.ylabel("Count")
    plt.title("Redox potential distribution")
    plt.savefig("redox_distribution.png", dpi=300)
    plt.close()

    plt.figure(figsize=(8, 4))
    plt.hist(gap, bins=50)
    plt.xlabel("HOMO–LUMO gap (Eh)")
    plt.ylabel("Count")
    plt.title("HOMO–LUMO gap distribution")
    plt.savefig("gap_distribution.png", dpi=300)
    plt.close()

    # Scatter plot: Redox potential vs HOMO–LUMO gap
    plt.figure(figsize=(6, 5))
    plt.scatter(gap, redox, alpha=0.6, edgecolors='w', linewidth=0.5)
    plt.xlabel("HOMO–LUMO gap (Eh)")
    plt.ylabel("Redox potential (V)")
    plt.title("Redox potential vs HOMO–LUMO gap")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig("redox_vs_gap.png", dpi=300)
    plt.close()

    # Fingerprint Generation
    mols = [Chem.MolFromSmiles(sm) for sm in smiles]

    # Remove any molecules that RDKit failed to parse (resulting in None)
    valid_mols_indices = [i for i, mol in enumerate(mols) if mol is not None]
    mols_sel = [mols[i] for i in valid_mols_indices]

    # Re-filter redox and gap if any molecules failed to parse
    if len(valid_mols_indices) < len(smiles):
        redox = redox.iloc[valid_mols_indices]
        gap = gap[valid_mols_indices]                     # align gap with mols_sel

    # gap as column vector after filtering
    gap = gap.reshape(-1, 1)

    if fingerprint == "morgan":
        fp = [AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=int(bits)) for mol in tqdm.tqdm(mols_sel)]
    else:
        fp = [rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=int(bits)) for mol in tqdm.tqdm(mols_sel)]
    
    # Convert fingerprints to a NumPy array (rows = molecules, cols = bits)
    X_fingerprint = np.array([np.array(fp_i) for fp_i in fp], dtype=int)
    Y = redox.to_numpy()
    # Combine fingerprint + HOMO–LUMO gap into one feature matrix

    if column == "with":
        X = np.hstack((X_fingerprint, gap))  # shape: (n_samples, n_bits + 1)
    else:
        X = X_fingerprint


    print("\n")
    print(f"X (Features, Fingerprints{' + HOMO-LUMO gap' if column == 'with' else ''}) shape: {X.shape}")
    print(f"Y (Target, Redox Potential) shape: {Y.shape}")
    print("\n")
    return X, Y


def learning_model(X, Y, input_regressor, column, bits, fingerprint):
    lsvr = LinearSVR(C=0.1, max_iter=5000, random_state=42, verbose=0)  # Linear SVR
    sgd = SGDRegressor(penalty='elasticnet', alpha=0.0001, eta0=0.05, max_iter=500, random_state=42, verbose=0)  # SGD
    rr = Ridge(alpha=5, solver='auto', random_state=42)  # Ridge Regression
    rfr = RandomForestRegressor(max_depth=None, random_state=42, verbose=0, n_jobs=-1)  # Random Forest
    svr = SVR(kernel='rbf', C=100, epsilon=0.5, gamma="auto", max_iter=500)  # SVR

    if input_regressor == "lsvr":
        regressor = lsvr
    elif input_regressor == "sgd":
        regressor = sgd
    elif input_regressor == "rr":
        regressor = rr
    elif input_regressor == "rfr":
        regressor = rfr
    elif input_regressor == "svr":
        regressor = svr
    else:
        raise ValueError("Invalid regression model selected")

    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.1, random_state=42)
    regressor.fit(X_train, Y_train)
    train_r2 = regressor.score(X_train, Y_train)
    test_r2 = regressor.score(X_test, Y_test)
    print('Training R2 =', train_r2)
    print(' Testing R2 =', test_r2)
    Y_pred = regressor.predict(X_test)

    mae = mean_absolute_error(Y_test, Y_pred)
    rmse = np.sqrt(mean_squared_error(Y_test, Y_pred))

    print(f"\n")
    print(f"Mean Absolute Error (MAE): {mae:.4f} V")
    print(f"Root Mean Squared Error (RMSE): {rmse:.4f} V")
    print(f"\n")

    plt.figure(figsize=(8, 8))
    plt.scatter(Y_test, Y_pred, c='dodgerblue', alpha=0.6, edgecolors='w', linewidth=0.5)
    min_val = min(Y_test.min(), Y_pred.min())
    max_val = max(Y_test.max(), Y_pred.max())
    plt.axline((min(Y_test), min(Y_test)), slope=1.0, color="black", linewidth=0.5, linestyle='--')
    plt.axline((min(Y_test), min(Y_test)+2), slope=1.0, color="red",linewidth=0.5, linestyle='-.')
    plt.axline((min(Y_test), min(Y_test)-2), slope=1.0, color="red",linewidth=0.5, linestyle='-.')
    plt.xlabel('True Redox Potential (V)')
    plt.ylabel('Predicted Redox Potential (V)')
    if column == "with":
        plt.title(f'{type(regressor).__name__} (trained with HOMO-LUMO gap)')
    else:
        plt.title(f'{type(regressor).__name__} ')
    plt.text(
        min_val,
        max_val * 0.95,
        f'Train $R^2$: {train_r2:.3f}\nTest $R^2$: {test_r2:.3f}\nMAE: {mae:.3f}',
        verticalalignment='top',
        fontsize=12,
        bbox=dict(boxstyle="round,pad=0.5", fc="white", alpha=0.7)
    )
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.axis('equal')  # Ensures the plot is square
    filename = f"{type(regressor).__name__}_{bits}_{fingerprint}_{column}.png".replace(" ", "")
    plt.savefig(filename, dpi=300)
    plt.close()
    return mae, regressor.score(X_train, Y_train), regressor.score(X_test, Y_test)


def cross_valditation(X, Y):
    cv = ShuffleSplit(n_splits=5, test_size=0.1, random_state=42)

    parameters_rfr = {'n_estimators': [100, 150],
                      'max_depth':[5,None],
                      'max_features':[None,]}

    params = parameters_rfr
    estimator = RandomForestRegressor()

    cv_grid = GridSearchCV(estimator, params, cv=cv, verbose=3) # Get the best parameters
    cv_grid.fit(X, Y)
    cv_scores = cross_validate(cv_grid.best_estimator_, X, Y, cv=cv, scoring=('r2', 'neg_mean_squared_error','neg_mean_absolute_error'), return_train_score=True) # Get the score

    print(f"best hyperparameter = {cv_grid.best_params_}\nmean cross-validated testing R2 = {cv_grid.best_score_}\nTraning R2: {cv_scores['train_r2']}\nTesting R2: {cv_scores['test_r2']}")
    print(f"Training MSE: {-cv_scores['train_neg_mean_squared_error']}\nTesting MSE: {-cv_scores['test_neg_mean_squared_error']}")
    print(f"Training MAE: {-cv_scores['train_neg_mean_absolute_error']}\nTesting MAE: {-cv_scores['train_neg_mean_absolute_error']}")
    return cv_grid.best_params_, cv_grid.best_score_, cv_scores


def main():
    with open("result.txt", "w") as file:
        column = ["with"] #train with HOMO-LUMO gap? (with/without)
        bits = [2048]
        fingerprint = ["atom"]
        input_regressor = ["rfr"]
        for c in column:
            for b in bits:
                for f in fingerprint:
                    for inp in input_regressor:
                        X, Y = get_data_fingerprint(c, b, f)
                        mae, traning, testing = learning_model(X, Y, inp, c, b, f)
                        best_params, best_score ,cv_scores= cross_valditation(X, Y)
                        file.write(f"{c}_{b}_{f}_{inp}.txt\n")
                        file.write(f"Traning R2 = {traning}\nTesting R2 = {testing}\nMAE = {mae}\n")
                        file.write("Cross validate:\n")
                        file.write(f"best hyperparameter = {best_params}\nmean cross-validated testing R2 = {best_score}\nTraning R2: {cv_scores['train_r2']}\nTesting R2: {cv_scores['test_r2']}\n")
                        file.write(f"Training MSE: {-cv_scores['train_neg_mean_squared_error']}\nTesting MSE: {-cv_scores['test_neg_mean_squared_error']}\n")
                        file.write(f"Training MAE: {-cv_scores['train_neg_mean_absolute_error']}\nTesting MAE: {-cv_scores['train_neg_mean_absolute_error']}\n")
                        file.write("\n")
main()
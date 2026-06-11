import numpy as np
import os, glob, re, pickle

def compute_yeo7_metrics(fc, yeo_labels):
    # Step 1 – restrict to cortical parcels (Yeo 1-7)
    cortical_mask = (yeo_labels >= 1) & (yeo_labels <= 7)
    fc_cx  = np.array(fc[np.ix_(cortical_mask, cortical_mask)], dtype=float, copy=True)
    yeo_cx = yeo_labels[cortical_mask].astype(int)

    # Step 2 – Fisher z-transform, remove self-connections
    fc_cx = np.clip(fc_cx, -0.999999, 0.999999)
    fc_z  = np.arctanh(fc_cx)
    np.fill_diagonal(fc_z, np.nan)

    # Step 3 – per-network segregation / integration
    K = 7
    segregation = np.full(K, np.nan)
    integration = np.full(K, np.nan)

    for k in range(1, K + 1):
        mask_k = (yeo_cx == k)
        mask_other = (yeo_cx != k)

        if mask_k.sum() >= 2:
            within = fc_z[np.ix_(mask_k, mask_k)].copy()
            np.fill_diagonal(within, np.nan)
            segregation[k - 1] = np.nanmean(within)

        if mask_k.sum() > 0 and mask_other.sum() > 0:
            integration[k - 1] = np.nanmean(fc_z[np.ix_(mask_k, mask_other)])

    segregation  = np.tanh(segregation)
    integration  = np.tanh(integration)
    normalized_segregation = (segregation - integration) / (segregation + integration)

    # Step 4 – positive-weight matrix for nodal graph metrics
    fc_pos = np.where(fc_z > 0, fc_z, 0.0)
    fc_pos = np.nan_to_num(fc_pos, nan=0.0)
    np.fill_diagonal(fc_pos, 0.0)

    N = fc_pos.shape[0]

    # Step 5 – participation coefficient (positive weights only)
    pc_parcel = np.full(N, np.nan)

    for i in range(N):
        k_i = np.sum(fc_pos[i, :])
        if k_i <= 0:
            continue

        sum_sq = 0.0
        for s in range(1, K + 1):
            k_is = np.sum(fc_pos[i, yeo_cx == s])
            sum_sq += (k_is / k_i) ** 2

        pc_parcel[i] = 1.0 - sum_sq

    pc_parcel = np.clip(pc_parcel, 0.0, 1.0)

    # Step 6 – within-module degree z-score (WMD)
    wmd_parcel = np.full(N, np.nan)

    for s in range(1, K + 1):
        mask_s = (yeo_cx == s)
        idx_s = np.where(mask_s)[0]

        if idx_s.size < 2:
            continue

        # within-module strength for each node in module s
        kappa_s = np.sum(fc_pos[np.ix_(idx_s, idx_s)], axis=1)

        mu_s = np.mean(kappa_s)
        sigma_s = np.std(kappa_s)

        if sigma_s <= 1e-12:
            wmd_parcel[idx_s] = 0.0
        else:
            wmd_parcel[idx_s] = (kappa_s - mu_s) / sigma_s

    return segregation, integration, normalized_segregation, pc_parcel, wmd_parcel, yeo_cx


def compute_all_metrics(all_data):
    for subj in all_data.values():
        yeo = subj["yeo_labels"]
        for sess_runs in subj["sessions"].values():
            for rd in sess_runs.values():
                seg, intg, normseg, pc_parcel, wmd_parcel, yeo_cx = compute_yeo7_metrics(rd["FC"], yeo)
                rd["segregation"] = seg
                rd["integration"] = intg
                rd["normalized_segregation"] = normseg
                rd["pc_parcel"] = pc_parcel
                rd["wmd_parcel"] = wmd_parcel
                rd["yeo_labels_cortical"] = yeo_cx


def load_subject(data_root, subject_id):
    """
    Load one subject's data.

    Returns
    -------
    sessions   : dict  {session_id -> {run_key -> {"condition": str, "FC": ndarray}}}
    yeo_labels : (N,) int array  Yeo-7 labels for the r200 parcellation
    """
    subj_dir = os.path.join(data_root, subject_id, "NF")

    # Yeo labels
    utils_path = os.path.join(subj_dir, f"{subject_id}_utils.pkl")
    with open(utils_path, "rb") as f:
        utils = pickle.load(f)
    yeo_labels = np.array(utils["yeo7"]["r200"]["yeoROIs"]).astype(int).ravel()

    # Session FC files
    fc_files = sorted(glob.glob(os.path.join(subj_dir, f"{subject_id}_V*_FC_Schaefer.pkl")))
    print(f"  {subject_id}: {len(fc_files)} session file(s)")

    sessions = {}
    for fc_path in fc_files:
        session_id = os.path.basename(fc_path).split("_")[1]   # V01, V02, ...
        with open(fc_path, "rb") as f:
            data = pickle.load(f)

        run_map = {}
        for key, fc in data["FC_FD"].items():
            if not key.endswith("_200"):
                continue
            m = re.match(r"V\d+_Run(\d+)(?:_(Transfer|NoFeedback))?_200", key)
            if m is None:
                continue
            run_key   = f"Run{m.group(1)}"
            condition = m.group(2) if m.group(2) else "Feedback"
            run_map[run_key] = {"condition": condition, "FC": fc}

        sessions[session_id] = run_map

    return sessions, yeo_labels
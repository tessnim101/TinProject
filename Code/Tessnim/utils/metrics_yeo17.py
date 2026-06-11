import numpy as np
import pandas as pd
import os, glob, re, pickle

# ── Yeo-17 sub-network order (matches yeo_17_to_7.ipynb) ─────────────────────
# Integer IDs are 1-indexed (1 = VisCent, ..., 17 = TempPar)
YEO17_ORDER = [
    "VisCent", "VisPeri", "SomMotA", "SomMotB",
    "DorsAttnA", "DorsAttnB", "SalVentAttnA", "SalVentAttnB",
    "LimbicB", "LimbicA", "ContA", "ContB", "ContC",
    "DefaultA", "DefaultB", "DefaultC", "TempPar",
]
YEO17_NAME_TO_ID = {name: idx + 1 for idx, name in enumerate(YEO17_ORDER)}
YEO17_ID_TO_NAME = {idx + 1: name for idx, name in enumerate(YEO17_ORDER)}


# ── Lookup table loader ───────────────────────────────────────────────────────
def load_yeo17_labels(yeo7_labels, lookup_csv_path):
    """
    Convert a Yeo-7 label array to a Yeo-17 integer label array
    using the matched_labels_exact.csv lookup table.

    The lookup table has 200 rows — one per cortical Schaefer parcel —
    in the same Yeo-7 parcel order as yeo7_labels.
    Subcortical parcels (yeo7_labels == 0 or > 7) are assigned label 0.

    Parameters
    ----------
    yeo7_labels     : (N,) int array  — Yeo-7 labels (1-7 cortical, 0/8+ subcortical)
    lookup_csv_path : str             — path to matched_labels_exact.csv

    Returns
    -------
    yeo17_labels : (N,) int array  — Yeo-17 labels (1-17 cortical, 0 subcortical)
    """
    df_matched = pd.read_csv(lookup_csv_path)

    # Extract sub-network name from full label string:
    # "17Networks_LH_VisCent_ExStr_1" → "VisCent"  (index 2 after split on "_")
    net17_per_parcel = df_matched["label_df2"].apply(
        lambda x: x.split("_")[2]
    ).values                                              # length 200, Yeo-7 order

    id17_per_parcel = np.array(
        [YEO17_NAME_TO_ID[n] for n in net17_per_parcel], dtype=int
    )

    # Identify cortical parcels in yeo7_labels
    cortical_mask    = (yeo7_labels >= 1) & (yeo7_labels <= 7)
    cortical_indices = np.where(cortical_mask)[0]

    if len(cortical_indices) != len(id17_per_parcel):
        raise ValueError(
            f"Cortical parcel count mismatch: "
            f"yeo7_labels has {len(cortical_indices)} cortical entries "
            f"but matched_labels_exact.csv has {len(id17_per_parcel)} rows. "
            f"Check that the lookup table corresponds to the same parcellation."
        )

    # Build output array — 0 for subcortical, 1-17 for cortical
    yeo17_labels = np.zeros_like(yeo7_labels, dtype=int)
    yeo17_labels[cortical_indices] = id17_per_parcel

    return yeo17_labels


# ── Core metric computation ───────────────────────────────────────────────────
def compute_yeo17_metrics(fc, yeo17_labels):
    """
    Compute segregation, integration, participation coefficient, and
    within-module degree z-score using Yeo-17 sub-networks as modules.

    Everything is identical to compute_yeo7_metrics except:
      - K = 17  (17 modules instead of 7)
      - yeo17_labels (integer 1-17) drives all module assignments

    Parameters
    ----------
    fc           : (N, N) float array  — functional connectivity matrix
    yeo17_labels : (N,) int array      — Yeo-17 labels (1-17 cortical, 0 subcortical)

    Returns
    -------
    segregation            : (17,) float — mean within-network FC per sub-network
    integration            : (17,) float — mean between-network FC per sub-network
    normalized_segregation : (17,) float — (seg - int) / (seg + int)
    pc_parcel              : (N_cortical,) float — participation coefficient from positive edges (0-1)
    pc_parcel_neg          : (N_cortical,) float — participation coefficient from negative edges (0-1)
    wmd_parcel             : (N_cortical,) float — within-module degree z-score
    yeo_cx                 : (N_cortical,) int   — Yeo-17 labels for cortical parcels
    net_names              : list of str          — sub-network name for each of the 17 entries
    """
    K = 17

    # Step 1 — restrict to cortical parcels
    cortical_mask = (yeo17_labels >= 1) & (yeo17_labels <= K)
    fc_cx  = np.array(fc[np.ix_(cortical_mask, cortical_mask)], dtype=float, copy=True)
    yeo_cx = yeo17_labels[cortical_mask].astype(int)

    # Step 2 — Fisher z-transform, remove self-connections
    fc_cx = np.clip(fc_cx, -0.999999, 0.999999)
    fc_z  = np.arctanh(fc_cx)
    np.fill_diagonal(fc_z, np.nan)

    # Step 3 — per-sub-network segregation and integration
    segregation = np.full(K, np.nan)
    integration = np.full(K, np.nan)

    for k in range(1, K + 1):
        mask_k     = (yeo_cx == k)
        mask_other = (yeo_cx != k)

        if mask_k.sum() >= 2:
            within = fc_z[np.ix_(mask_k, mask_k)].copy()
            np.fill_diagonal(within, np.nan)
            segregation[k - 1] = np.nanmean(within)

        if mask_k.sum() > 0 and mask_other.sum() > 0:
            integration[k - 1] = np.nanmean(fc_z[np.ix_(mask_k, mask_other)])

    # Back-transform from Fisher z
    segregation            = np.tanh(segregation)
    integration            = np.tanh(integration)
    normalized_segregation = (segregation - integration) / (segregation + integration)

    # Step 4 — positive- and negative-weight matrices for nodal graph metrics
    fc_pos = np.where(fc_z > 0, fc_z, 0.0)
    fc_pos = np.nan_to_num(fc_pos, nan=0.0)
    np.fill_diagonal(fc_pos, 0.0)

    fc_neg = np.where(fc_z < 0, -fc_z, 0.0)   # absolute values of negative edges
    fc_neg = np.nan_to_num(fc_neg, nan=0.0)
    np.fill_diagonal(fc_neg, 0.0)

    N = fc_pos.shape[0]

    # Step 5 — participation coefficient using 17 modules
    # PC_i = 1 - sum_s( k_is / k_i )^2
    # where k_is = strength to module s, k_i = total strength
    pc_parcel = np.full(N, np.nan)
    pc_parcel_neg = np.full(N, np.nan)

    for i in range(N):
        k_i = np.sum(fc_pos[i, :])
        if k_i > 0:
            sum_sq = 0.0
            for s in range(1, K + 1):
                k_is    = np.sum(fc_pos[i, yeo_cx == s])
                sum_sq += (k_is / k_i) ** 2
            pc_parcel[i] = 1.0 - sum_sq

        k_i_neg = np.sum(fc_neg[i, :])
        if k_i_neg > 0:
            sum_sq_neg = 0.0
            for s in range(1, K + 1):
                k_is_neg    = np.sum(fc_neg[i, yeo_cx == s])
                sum_sq_neg += (k_is_neg / k_i_neg) ** 2
            pc_parcel_neg[i] = 1.0 - sum_sq_neg

    # Clip to [0, 1] — any value outside this range is a numerical error
    pc_parcel     = np.clip(pc_parcel,     0.0, 1.0)
    pc_parcel_neg = np.clip(pc_parcel_neg, 0.0, 1.0)

    # Step 6 — within-module degree z-score using 17 modules
    wmd_parcel = np.full(N, np.nan)

    for s in range(1, K + 1):
        mask_s = (yeo_cx == s)
        idx_s  = np.where(mask_s)[0]

        if idx_s.size < 2:
            continue

        kappa_s = np.sum(fc_pos[np.ix_(idx_s, idx_s)], axis=1)
        mu_s    = np.mean(kappa_s)
        sigma_s = np.std(kappa_s)

        if sigma_s <= 1e-12:
            wmd_parcel[idx_s] = 0.0
        else:
            wmd_parcel[idx_s] = (kappa_s - mu_s) / sigma_s

    # Sub-network name for each of the 17 output slots
    net_names = [YEO17_ID_TO_NAME[k] for k in range(1, K + 1)]

    return segregation, integration, normalized_segregation, pc_parcel, pc_parcel_neg, wmd_parcel, yeo_cx, net_names


# ── Apply to all subjects ─────────────────────────────────────────────────────
def compute_all_metrics_yeo17(all_data, lookup_csv_path):
    """
    Same interface as compute_all_metrics() in metrics.py,
    but uses Yeo-17 sub-networks as modules.

    Adds to each run dict:
      segregation17, integration17, normalized_segregation17,
      pc_parcel17, wmd_parcel17, yeo_labels_cortical17, net_names17
    """
    for subj in all_data.values():
        yeo7  = subj["yeo_labels"]                          # original Yeo-7 labels
        yeo17 = load_yeo17_labels(yeo7, lookup_csv_path)   # converted to Yeo-17
        subj["yeo17_labels"] = yeo17                        # store on subject

        for sess_runs in subj["sessions"].values():
            for rd in sess_runs.values():
                (seg, intg, normseg,
                 pc_parcel, pc_parcel_neg, wmd_parcel,
                 yeo_cx, net_names) = compute_yeo17_metrics(rd["FC"], yeo17)

                rd["segregation17"]            = seg
                rd["integration17"]            = intg
                rd["normalized_segregation17"] = normseg
                rd["pc_parcel17"]              = pc_parcel
                rd["pc_parcel_neg17"]          = pc_parcel_neg
                rd["wmd_parcel17"]             = wmd_parcel
                rd["yeo_labels_cortical17"]    = yeo_cx
                rd["net_names17"]              = net_names


# ── Subject loader (identical to metrics.py) ──────────────────────────────────
def load_subject(data_root, subject_id):
    """
    Identical to load_subject() in metrics.py.
    Loads FC matrices and Yeo-7 labels — Yeo-17 conversion happens
    in compute_all_metrics_yeo17() using the lookup table.
    """
    subj_dir = os.path.join(data_root, subject_id, "NF")

    utils_path = os.path.join(subj_dir, f"{subject_id}_utils.pkl")
    with open(utils_path, "rb") as f:
        utils = pickle.load(f)
    yeo_labels = np.array(utils["yeo7"]["r200"]["yeoROIs"]).astype(int).ravel()

    fc_files = sorted(glob.glob(os.path.join(subj_dir, f"{subject_id}_V*_FC_Schaefer.pkl")))
    print(f"  {subject_id}: {len(fc_files)} session file(s)")

    sessions = {}
    for fc_path in fc_files:
        session_id = os.path.basename(fc_path).split("_")[1]
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

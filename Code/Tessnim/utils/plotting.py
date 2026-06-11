
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
# -------------------------------------------------------------------------
# Helper: collect per-subject rows (list of dicts, one per session×run)
# -------------------------------------------------------------------------
from matplotlib import axes


def collect_rows(all_data, condition_filter=None):
    """
    Build a list of row-dicts from all subjects.
    Each row: subject, session, run, condition, segregation (7,), integration (7,),
              normalized_segregation (7,).

    condition_filter : str or None  – keep only this condition if given
    """
    all_rows = {}
    for sid, subj in all_data.items():
        rows = []
        for sess_id, sess_runs in subj["sessions"].items():
            for run_key in sorted(sess_runs.keys(), key=lambda x: int(x.replace("Run", ""))):
                rd = sess_runs[run_key]
                if condition_filter and rd["condition"] != condition_filter:
                    continue
                rows.append({
                    "subject":               sid,
                    "session":               sess_id,
                    "run":                   run_key,
                    "label":                 f"{sess_id}_{run_key}",
                    "condition":             rd["condition"],
                    "segregation":           rd["segregation"],
                    "integration":           rd["integration"],
                    "normalized_segregation": rd["normalized_segregation"],
                })
        all_rows[sid] = rows
    return all_rows


# -------------------------------------------------------------------------
# Plot 1: all runs, x = Session_Run
# -------------------------------------------------------------------------
def plot_runs_all_subjects(
    all_data,
    yeo_names,
    condition_filter=None,
    average_only=False,
    shade_conditions=False,
    show_subject_counts=False,
    save_path=None,  # ← NEW ARGUMENT: path where to save (e.g. "plots/seg_int_norm.png")
):
    """
    3×1 stacked plot (Segregation / Integration / Norm. Segregation).
    X-axis: every Session_Run label (union across subjects, in sorted order).
    Each subject is a semi-transparent line; group average is bold black.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    subject_rows = collect_rows(all_data, condition_filter)

    # Build the union of all Session_Run labels in sorted order
    all_labels_set = set()
    for rows in subject_rows.values():
        for r in rows:
            all_labels_set.add(r["label"])
    all_labels = sorted(all_labels_set)
    label_to_x = {lbl: i for i, lbl in enumerate(all_labels)}
    x_all = np.arange(len(all_labels))

    # Map each label to its condition (take from any subject that has it)
    label_to_condition = {}
    for rows in subject_rows.values():
        for r in rows:
            label_to_condition[r["label"]] = r["condition"]

    shade_colors = {
        "Transfer": ("orange", 0.18),
        "NoFeedback": ("green", 0.18),
    }

    # Accumulate values per label for the group average
    group_acc = {
        lbl: {"segregation": [], "integration": [], "normalized_segregation": []}
        for lbl in all_labels
    }
    label_subjects = {lbl: set() for lbl in all_labels}

    for sid, rows in subject_rows.items():
        for r in rows:
            label_subjects[r["label"]].add(sid)

    label_counts = np.array([len(label_subjects[lbl]) for lbl in all_labels])

    for rows in subject_rows.values():
        for r in rows:
            lbl = r["label"]
            group_acc[lbl]["segregation"].append(r["segregation"])
            group_acc[lbl]["integration"].append(r["integration"])
            group_acc[lbl]["normalized_segregation"].append(r["normalized_segregation"])

    # Group average arrays  (n_labels × 7)
    seg_avg = np.vstack([np.nanmean(group_acc[l]["segregation"], axis=0) for l in all_labels])
    int_avg = np.vstack([np.nanmean(group_acc[l]["integration"], axis=0) for l in all_labels])
    normseg_avg = np.vstack(
        [np.nanmean(group_acc[l]["normalized_segregation"], axis=0) for l in all_labels]
    )

    cond_str = f" [{condition_filter}]" if condition_filter else ""
    fig, axes = plt.subplots(3, 1, figsize=(18, 14), sharex=True)

    if shade_conditions:
        for lbl in all_labels:
            cond = label_to_condition.get(lbl)
            if cond in shade_colors:
                color, alpha = shade_colors[cond]
                xi = label_to_x[lbl]
                for ax in axes:
                    ax.axvspan(xi - 0.5, xi + 0.5, color=color, alpha=alpha, zorder=0)
        # Add shading legend entries
        from matplotlib.patches import Patch

        shade_handles = [
            Patch(color=c, alpha=a, label=cond)
            for cond, (c, a) in shade_colors.items()
        ]
        axes[1].legend(handles=shade_handles, loc="center left", bbox_to_anchor=(1, 0.5))

    if not average_only:
        for sid, rows in subject_rows.items():
            if not rows:
                continue
            xs = [label_to_x[r["label"]] for r in rows]
            segs = np.vstack([r["segregation"] for r in rows])
            ints = np.vstack([r["integration"] for r in rows])
            nrms = np.vstack([r["normalized_segregation"] for r in rows])

            for i in range(7):
                axes[0].plot(xs, segs[:, i], marker="o", alpha=0.35, linewidth=1)
                axes[1].plot(xs, ints[:, i], marker="o", alpha=0.35, linewidth=1)
                axes[2].plot(xs, nrms[:, i], marker="o", alpha=0.35, linewidth=1)

    # Group average — bold, labeled
    for i in range(7):
        axes[0].plot(x_all, seg_avg[:, i], marker="o", linewidth=2, label=yeo_names[i])
        axes[1].plot(x_all, int_avg[:, i], marker="o", linewidth=2, label=yeo_names[i])
        axes[2].plot(x_all, normseg_avg[:, i], marker="o", linewidth=2, label=yeo_names[i])

    axes[0].set_ylabel("Segregation")
    axes[0].set_title(f"Within-network connectivity (Segregation){cond_str}")

    axes[1].set_ylabel("Integration")
    axes[1].set_title(f"Between-network connectivity (Integration){cond_str}")

    axes[2].set_ylabel("Normalized segregation")
    axes[2].set_title(f"Normalized segregation across runs{cond_str}")
    axes[2].axhline(0, linestyle="--", linewidth=1)

    axes[2].set_xticks(x_all)
    axes[2].set_xticklabels(all_labels, rotation=90)
    axes[2].set_xlabel("Session / Run")

    if show_subject_counts:
        ax_count = axes[2].twinx()
        ax_count.plot(x_all, label_counts, color="black", linestyle="--", marker="s")
        ax_count.set_ylabel("Subjects")

    axes[0].legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.tight_layout()

    # Save the figure if save_path is given
    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Plot saved to {save_path}")

    plt.show()


# -------------------------------------------------------------------------
# Plot 2: session averages, x = session
# -------------------------------------------------------------------------
def plot_session_averages_all_subjects(all_data, yeo_names, condition_filter=None, average_only=False):
    """
    3×1 stacked plot of session-averaged metrics.
    X-axis: session labels (union across subjects).
    Each subject is a semi-transparent line; group average is bold black.
    """
    subject_rows = collect_rows(all_data, condition_filter)

    # Union of sessions
    all_sessions_set = set()
    for rows in subject_rows.values():
        for r in rows:
            all_sessions_set.add(r["session"])
    all_sessions = sorted(all_sessions_set)
    x_sessions = np.arange(len(all_sessions))

    # Per-subject session averages
    def session_avg_for_subject(rows):
        """Returns (n_sessions × 7) arrays for seg / int / normseg."""
        session_groups = defaultdict(list)
        for r in rows:
            session_groups[r["session"]].append(r)

        seg_sv, int_sv, normseg_sv = [], [], []
        for sess in all_sessions:
            sess_rows = session_groups.get(sess, [])
            if sess_rows:
                seg_sv.append(np.nanmean(np.vstack([r["segregation"]            for r in sess_rows]), axis=0))
                int_sv.append(np.nanmean(np.vstack([r["integration"]            for r in sess_rows]), axis=0))
                normseg_sv.append(np.nanmean(np.vstack([r["normalized_segregation"] for r in sess_rows]), axis=0))
            else:
                seg_sv.append(np.full(7, np.nan))
                int_sv.append(np.full(7, np.nan))
                normseg_sv.append(np.full(7, np.nan))

        return (np.vstack(seg_sv), np.vstack(int_sv), np.vstack(normseg_sv))

    per_subject = {sid: session_avg_for_subject(rows)
                   for sid, rows in subject_rows.items() if rows}

    # Group average across subjects
    seg_group     = np.nanmean(np.stack([v[0] for v in per_subject.values()]), axis=0)
    int_group     = np.nanmean(np.stack([v[1] for v in per_subject.values()]), axis=0)
    normseg_group = np.nanmean(np.stack([v[2] for v in per_subject.values()]), axis=0)

    cond_str = f" [{condition_filter}]" if condition_filter else ""
    fig, axes = plt.subplots(3, 1, figsize=(18, 14), sharex=True)

    # Individual subjects (faded)
    if not average_only:
        for sid, (seg_sv, int_sv, normseg_sv) in per_subject.items():
            for i in range(7):
                axes[0].plot(x_sessions, seg_sv[:, i],     marker="o", alpha=0.3, linewidth=1)
                axes[1].plot(x_sessions, int_sv[:, i],     marker="o", alpha=0.3, linewidth=1)
                axes[2].plot(x_sessions, normseg_sv[:, i], marker="o", alpha=0.3, linewidth=1)

    # Group average — bold, labeled
    for i in range(7):
        axes[0].plot(x_sessions, seg_group[:, i],     marker="o", linewidth=2, label=yeo_names[i])
        axes[1].plot(x_sessions, int_group[:, i],     marker="o", linewidth=2, label=yeo_names[i])
        axes[2].plot(x_sessions, normseg_group[:, i], marker="o", linewidth=2, label=yeo_names[i])

    axes[0].set_ylabel("Segregation")
    axes[0].set_title(f"Within-network connectivity (Segregation) — Session averages (all runs){cond_str}")

    axes[1].set_ylabel("Integration")
    axes[1].set_title(f"Between-network connectivity (Integration) — Session averages (all runs){cond_str}")

    axes[2].set_ylabel("Normalized segregation")
    axes[2].set_title(f"Normalized segregation across sessions — Session averages (all runs){cond_str}")
    axes[2].axhline(0, linestyle="--", linewidth=1)

    axes[2].set_xticks(x_sessions)
    axes[2].set_xticklabels(all_sessions)
    axes[2].set_xlabel("Session")

    axes[0].legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.show()


print("Plotting functions defined.")
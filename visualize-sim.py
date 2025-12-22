import numpy as np
from pathlib import Path
import re
import matplotlib.pyplot as plt
from matplotlib import colors
from pathlib import Path

cmap_bv = colors.LinearSegmentedColormap.from_list('bv',["black", "red"],256)
cmap_bv._init()

alphas = np.linspace(0.0, 1.0, cmap_bv.N+3)
cmap_bv._lut[:,-1] = alphas

tissues = ["Breast", "Bones", "Lungs", "Liver"]

def _find_cell_array(data, cell_type):
    key = cell_type.lower()[:3]
    if key in data:
        return data[key]
    for k in data.files:
        if k.lower().startswith(cell_type.lower().split()[0][:3]):
            return data[k]
    for k in data.files:
        arr = data[k]
        if isinstance(arr, np.ndarray) and arr.ndim == 2:
            return arr
    return None

def make_2x2_plot(iteration, cell_type, base_dir=Path("./arrays")):
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    axs = axs.flatten()
    for ax, tissue in zip(axs, tissues):
        p = base_dir / f"iter-{tissue.lower()}-{iteration}.npz"
        if not p.exists():
            ax.set_title(f"{tissue} - missing file")
            ax.axis("off")
            continue
        try:
            data = np.load(p)
            grid = _find_cell_array(data, cell_type)
            if grid is None:
                ax.set_title(f"{tissue} - no matching cell array")
                ax.axis("off")
                continue
            ax.imshow(grid, cmap="inferno", vmin=0, vmax=4)
            if "bv" in data:
                ax.imshow(data["bv"], cmap=cmap_bv, vmin=0, vmax=1)
            ax.axis("off")
            ax.set_title(f"{tissue} - {cell_type}")
        except Exception:
            ax.set_title(f"{tissue} - error loading")
            ax.axis("off")
    plt.suptitle(f"Iteration {iteration} - {cell_type}", fontsize=16)
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.show()


make_2x2_plot(100, "Mesenchymal Cells")
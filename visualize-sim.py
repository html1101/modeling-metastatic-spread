import numpy as np
from pathlib import Path
import re

import matplotlib.pyplot as plt

TISSUE = "Bones"
ITERATION = 13250
CELL_TYPE = "Epithelial Cells"

# Load and display .npz files
data_file = F"./arrays/iter-{TISSUE.lower()}-{ITERATION}.npz"

data = np.load(data_file)
grid = data[CELL_TYPE.lower()[:3]]  # Adjust key name if needed

plt.figure(figsize=(8, 6))
plt.imshow(grid, cmap='inferno', vmin=0, vmax=4)
plt.axis("off")
plt.colorbar()
plt.title(F"{TISSUE} Tissue - {CELL_TYPE} (Iteration {ITERATION})")
plt.show()

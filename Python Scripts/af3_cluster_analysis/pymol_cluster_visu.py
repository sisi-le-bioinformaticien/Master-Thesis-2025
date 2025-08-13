import numpy as np

model_name = "mut"
model_count = 10
ref_model = f"{model_name}1"

# Define a list of distinguishable colors (PyMOL built-ins or custom RGBs)
colors = [
    "red", "green", "blue", "yellow", "cyan", "magenta", 
    "orange", "purple", "teal", "slate"
]

# Load and align all models
for i in range(1, model_count + 1):
    model = f"{model_name}{i}"
    cmd.load(f"./{model}.cif", model)

for i in range(2, model_count + 1):
    cmd.align(f"{model_name}{i}", ref_model)

# Style 
for i in range(1, model_count + 1):
    model = f"{model_name}{i}"
    color = colors[(i - 1) % len(colors)]  

    # Show protein chains as cartoon with transparency
    cmd.show_as("cartoon", f"{model} and not chain Z")
    cmd.set("cartoon_transparency", 0.2, model)
    cmd.color(color, f"{model} and not chain Z")

    # Show drug (Z chain) as lines and keep opaque
    cmd.show_as("lines", f"{model} and chain Z")
    cmd.color("white", f"{model} and chain Z")

# Clean up view
cmd.orient()

---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Appendix Z. Orders of Magnitude

## Volume

```{code-cell} ipython3
import matplotlib.pyplot as plt

# Define a few course-relevant examples (volume in m³) and labels.
# These examples span from atomic to macroscopic scales.
# (You may adjust the values or add more points as desired.)
examples = [
    (1e-30, 'Atom'),
    (1e-27, 'Molecule'),
    (1e-23, 'Protein'),
    (1e-15, 'Cell'),
    (2.24e-2, 'Molar Gas'),
    (1e-1, 'Human')
]

# Create the figure and axis.
fig, ax = plt.subplots(figsize=(10, 2))

# Set the x-axis to logarithmic scale.
ax.set_xscale('log')

# Hide the y-axis (we only need a horizontal line).
ax.get_yaxis().set_visible(False)

# Set x-axis limits to cover the data range.
ax.set_xlim(1e-30, 1e0)

# Draw a horizontal axis line at y = 0.
ax.axhline(0, color='black', lw=1)

# Plot each example as a vertical dash with a text label.
# We alternate the label placement above and below the axis to reduce crowding.
for i, (x, label) in enumerate(examples):
    # Draw a vertical dash at x.
    # ymin and ymax are in axis-relative coordinates (0 to 1)
    ax.axvline(x=x, ymin=0.45, ymax=0.55, color='black')
    
    # Alternate label positions: even index above, odd below.
    if i % 2 == 0:
        y_text = 0.6  # above the axis
        va = 'bottom'
    else:
        y_text = 0.4  # below the axis
        va = 'top'
        
    # Place the text with a rotation so that the "bottom" of the text faces right.
    # A rotation of -90 degrees rotates the text 90° clockwise.
    ax.text(x, y_text, label, rotation=-90, ha='center', va=va, fontsize=10)

# Optionally, remove other spines for a cleaner look.
for spine in ['top', 'right', 'left']:
    ax.spines[spine].set_visible(False)

# Label the x-axis.
ax.set_xlabel("Volume (m³, log scale)")

plt.tight_layout()
plt.show()
```

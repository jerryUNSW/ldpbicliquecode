import matplotlib.pyplot as plt
import subprocess

# Coordinates for nodes
coords = {
    'u_1': (0, 1),
    'u_2': (2, 1),
    'v_1': (-1, 0),
    'v_2': (1, 0),
    'v_3': (3, 0),
}

motifs = [
    ([], r"$i=0$"),
    ([("u_1", "v_1")], r"$i=1$"),
    ([("u_1", "v_1"), ("u_1", "v_2")], r"$i=2$"),
    ([("u_1", "v_1"), ("u_2", "v_2")], r"$i=2$"),
    ([("u_1", "v_1"), ("u_2", "v_1")], r"$i=2$"),
    ([("u_1", "v_1"), ("u_1", "v_2"), ("u_1", "v_3")], r"$i=3$"),
    ([("u_1", "v_1"), ("u_1", "v_2"), ("u_2", "v_3")], r"$i=3$"),
    ([("u_1", "v_1"), ("u_2", "v_1"), ("u_1", "v_2")], r"$i=3$"),
    ([("u_1", "v_1"), ("u_2", "v_1"), ("u_2", "v_2")], r"$i=3$"),
    ([("u_1", "v_1"), ("u_1", "v_2"), ("u_2", "v_2"), ("u_2", "v_3")], r"$i=4$"),
    ([("u_1", "v_1"), ("u_1", "v_2"), ("u_2", "v_1"), ("u_2", "v_2")], r"$i=4$"),
    ([("u_1", "v_1"), ("u_1", "v_2"), ("u_2", "v_1"), ("u_2", "v_3")], r"$i=4$"),
    ([("u_1", "v_1"), ("u_1", "v_3"), ("u_2", "v_1"), ("u_2", "v_2")], r"$i=4$"),
    ([("u_1", "v_1"), ("u_1", "v_2"), ("u_1", "v_3"), ("u_2", "v_1"), ("u_2", "v_2")], r"$i=5$"),
    ([("u_1", "v_1"), ("u_1", "v_2"), ("u_2", "v_1"), ("u_2", "v_2"), ("u_2", "v_3")], r"$i=5$"),
    ([("u_1", "v_1"), ("u_1", "v_2"), ("u_1", "v_3"), ("u_2", "v_1"), ("u_2", "v_2"), ("u_2", "v_3")], r"$i=6$"),
]

# Set up an 8×2 grid of subplots
fig, axes = plt.subplots(2, 8, figsize=(20, 6))  # wider figure
axes = axes.flatten()

# Draw each motif
for ax, (edges, label) in zip(axes, motifs):
    for u, v in edges:
        x = [coords[u][0], coords[v][0]]
        y = [coords[u][1], coords[v][1]]
        ax.plot(x, y, linewidth=2)
    
    for node, (x, y) in coords.items():
        # Set color based on node type
        if node.startswith('u_'):
            facecolor = 'lightgrey'
        else:  # v_ nodes
            facecolor = 'white'
        
        ax.scatter(x, y, s=600, edgecolors='black', facecolors=facecolor, linewidths=2, zorder=3)
        ax.text(x, y, f"${node}$", fontsize=18, ha='center', va='center')
    
    ax.set_title(label, fontsize=18)
    ax.axis('off')
    ax.set_xlim(-2, 3.5)
    ax.set_ylim(-0.5, 1.5)

# Hide unused axes
for i in range(len(motifs), len(axes)):
    axes[i].axis('off')

# Tighter layout
plt.subplots_adjust(
    left=0.02, right=0.98, 
    top=0.9, bottom=0.1, 
    wspace=0.2, hspace=0.4
)

# Save and close
plt.savefig("motif_diagrams.pdf")
plt.close()

subprocess.run(["bash", "transfer.sh"])
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
import matplotlib.patches as mpatches
import subprocess

import matplotlib.patheffects as pe


# Create figure and axis
fig, ax = plt.subplots(1, 1, figsize=(12, 8))

# Define positions for nodes
# Upper row (u nodes)
U_level = 4.5
V_level = 4
bottom_level = 3.5

# U_level = 5
# V_level = 4
# bottom_level = 3.5

node_size = 0.15
node_font = 25
equation_font = 25

u_positions = {
    'u2': (2, U_level),
    'u3': (4, U_level),
    'u4': (6, U_level),
}
# Lower row (v nodes)
v_positions = {
    'v1': (3, V_level),
    'v2': (5, V_level)
}



# Bottom node (u1)
u1_position = (4, bottom_level)

# Draw dotted blue lines from u nodes to v nodes
for u_label, u_pos in u_positions.items():
    for v_label, v_pos in v_positions.items():
        ax.plot([u_pos[0], v_pos[0]], [u_pos[1], v_pos[1]],
                color='#00B0F0', linestyle=':', linewidth=3, zorder=1)

# Draw solid black lines from v nodes to u1
for v_label, v_pos in v_positions.items():
    ax.plot([v_pos[0], u1_position[0]], [v_pos[1], u1_position[1]],
            color='black', linewidth=3, zorder=1)

# Draw upper nodes (white circles)
for label, pos in u_positions.items():
    circle = Circle(pos, node_size, color='grey', ec='black', linewidth=3, zorder=3)
    ax.add_patch(circle)
    # if label == 'up':
    #     display_label = '$u_p$'
    # else:
    # Extract the number from the label (e.g., 'u2' -> '2')
    idx = label[1:]
    display_label = f'$u_{idx}$'

    ax.text(pos[0], pos[1]+node_size*2, display_label, 
        ha='center', va='center', fontsize=node_font, 
        zorder=4)
    
    idx = label[1:]
    # annotation = rf'$\phi(v_1, u_{idx})\phi(v_2, u_{idx})+$'
    print("idx = ", idx)

    if idx == '4':
        annotation = rf'$\phi(v_1, u_{idx})\phi(v_2, u_{idx}) = f_1$'
    else:
        annotation = rf'$\phi(v_1, u_{idx})\phi(v_2, u_{idx})+$'
    
    ax.text(pos[0], pos[1] + 0.7, annotation, ha='center', va='center', fontsize=equation_font, zorder=4
        # path_effects=[pe.withStroke(linewidth=0.1, foreground='black')]
        )

# Draw middle nodes (gray circles)
for label, pos in v_positions.items():
    circle = Circle(pos, node_size, color='white', ec='black', linewidth=3, zorder=3)
    ax.add_patch(circle)
    # Extract the number from the label (e.g., 'v1' -> '1')
    idx = label[1:]
    display_label = f'$v_{idx}$'
    ax.text(pos[0], pos[1]-node_size*2, display_label, ha='center', va='center', fontsize=node_font, zorder=4)

# Draw bottom node (white circle)
circle = Circle(u1_position, node_size, color='grey', ec='black', linewidth=3, zorder=3)
ax.add_patch(circle)
ax.text(u1_position[0], u1_position[1]-node_size*2, '$u_1$', ha='center', va='center', fontsize=node_font, zorder=4)

# Add mathematical expressions
# Top equation
# eq_text = r'$\phi(v_1, u_2)\phi(v_2, u_2) + \phi(v_1, u_3)\phi(v_2, u_3) + \phi(v_1, u_4)\phi(v_2, u_4) + \phi(v_1, u_p)\phi(v_2, u_p) = f_1$'
# ax.text(5, 8.5, eq_text, ha='center', va='center', fontsize=14)

# Bottom equation with arrow
ax.text(6, bottom_level, 
    r'$\tilde{f}_1= f_1 + Lap(\frac{\gamma^{p-1}}{\varepsilon_2})$', 
        ha='left', va='center', fontsize=equation_font, 
        # path_effects=[pe.withStroke(linewidth=0.1, foreground='black')]
        )

ax.annotate(
    '', 
    xy=(7, 4), 
    xytext=(7, 4.5),
    arrowprops=dict(
        arrowstyle='-|>',     # sleek arrowhead
        color='black',
        lw=2.5,
        mutation_scale=20,      # <<< Bigger arrowhead
        # connectionstyle="arc3,rad=-0.2"  # slight curve
    ),
    zorder=2
)


# ax.text(7, U_level +0.7, r'$Lap(\frac{\gamma^{p-1}}{\varepsilon_2})$', 
#     ha='left', va='center', 
#         fontsize=equation_font,
#         fontweight='bold',  # optional
#         )




# # Set axis properties
# ax.set_xlim(0, 12)
ax.set_xlim(1, 7)
ax.set_ylim(bottom_level-0.5, U_level+0.5)
ax.set_aspect('equal')
ax.axis('off')
# Set tighter axis limits
# ax.set_xlim(1.5, 9.5)  # Adjusted to fit content
# ax.set_ylim(2.5, 6.0)  # Adjusted to fit content
# ax.set_aspect('equal')
# ax.axis('off')

# Add a border around the plot
# border = mpatches.Rectangle((0.2, 0.2), 11.6, 9.6, linewidth=1,
#                            edgecolor='black', facecolor='none')
# ax.add_patch(border)

plt.tight_layout()
plt.savefig("example.pdf", bbox_inches='tight')
plt.close()

subprocess.run(["bash", "transfer.sh"])
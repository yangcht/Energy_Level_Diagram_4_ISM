#### MUST USE BEFORE (almost) ANY CODE!
import matplotlib.pyplot as plt
import numpy as np
plt.style.use(['science'])
####

from matplotlib.path import Path
from astropy.io import ascii 
from astropy.table import join
from matplotlib.ticker import MultipleLocator, NullLocator

from adjustText import adjust_text


### Reading the Table
### This is suitable for any table downloaded from Cologne Molecular Data Base
col_ends = (13, 20, 29, 31, 41, 43, 54, 57, 59, 61, 63, 69, 71, 73, 75)
names = ['freq', 'freq_error', 'inten', 'col4', 'Energy_lower', 'col6', 'col7', 'J_u', 'Ka_u', 'Kc_u', 'F_u', 'J_l', 'Ka_l', 'Kc_l', 'F_l']

# Read the ASCII table using Astro Table
table = ascii.read('h2o.cat', format='fixed_width_no_header', col_ends=col_ends, names=names)

# adding the upper Energy column of to the table, 
# looping over all rows and find the J_l, Ka_l and Kc_l 
# combination maching J_u, Ka_u and Kc_u

table['Energy_upper'] = table['Energy_lower']

for idx, row in enumerate(table):
    mask = (table["J_l"] == row['J_u']) & (table["Ka_l"] == row['Ka_u']) & (table["Kc_l"] == row['Kc_u'])
    if sum(mask)>=1:
        table['Energy_upper'][idx] = table[mask]['Energy_lower'][0]
    else:
        table['Energy_upper'][idx] = np.nan
        
k_B = 0.695035 # in unit of cm-1/K

J_up = range(7)

# defining the segemant line for plotting the Energy levels

path = Path([(-2, 0.5), (2, 0.5)])

fig, ax = plt.subplots(ncols=2, figsize=(10,6))

# For plotting the energy levels for each, ortho and para
def plot_energy_levels(ax, J_values, E_values, Ka_values, Kc_values, offset):
    ax.plot(J_values * np.ones(len(E_values)), E_values, linestyle='', marker=path, markersize=17, markeredgewidth=2)
    for i, (x, y) in enumerate(zip(J_values * np.ones(len(E_values)), E_values)):
        ax.annotate(f'{J_values[i]}$_{Ka_values[i]}$$_{Kc_values[i]}$', xy=(x, y), xytext=(x + offset, y - 5), fontsize=18)


# Function for filtering for Ortho or Para tables --> derive E_K, J, Ka and Kc
def extract_variables(filter_table, condition):
    filtered_table = filter_table[condition]
    E_K = filtered_table['Energy_lower'] / k_B
    J = filtered_table['J_l']
    Ka = filtered_table['Ka_l']
    Kc = filtered_table['Kc_l']
    return E_K, J, Ka, Kc


# Plotting the diagram
for i, J in enumerate(J_up):
    # Plot the data for the current J_up value
    mask = (table["J_l"] == J) & (table["J_l"] <= 7)
    filter_table = table[mask]
    filter_table['ka+kc'] = filter_table['Ka_l'] + filter_table['Kc_l']

    ortho_mask = filter_table['ka+kc'] % 2 == 1
    para_mask  = filter_table['ka+kc'] % 2 == 0
    
    o_E_K, o_J, o_Ka, o_Kc = extract_variables(filter_table, ortho_mask)
    p_E_K, p_J, p_Ka, p_Kc = extract_variables(filter_table, para_mask)
    
    plot_energy_levels(ax[0], o_J, o_E_K, o_Ka, o_Kc, 0.75)
    plot_energy_levels(ax[1], p_J, p_E_K, p_Ka, p_Kc, 0.24) 
    
    
# Put the J Ka Kc annotation
ax[0].annotate(r'${\rm{H_2O}}\;\,J_{K{\rm a}\,K{\rm c}}$', xy=(max(J_up)*1.1, 28), fontsize=30)

for axis in ax:
    axis.set_ylim(-30, max(o_E_K)*1.1)
    axis.set_xlim(-0.5, max(J_up)*1.2)

# Set the shared X axis
ax[0].set_ylabel(r'$E/k_B {\rm(K)}$')

# Only show one sets of y tick values
ax[1].set_yticklabels([])

# Invert the left figure horizontally
ax[0].invert_xaxis()

# Set the x axis major tick locator to show only every other tick
major_tick_locator = MultipleLocator(1)
minor_tick_locator = NullLocator()
for axis in ax:
    axis.xaxis.set_major_locator(major_tick_locator)
    axis.xaxis.set_minor_locator(minor_tick_locator)

# Put annotation of the ortho and para 
ax[0].annotate(r'ortho', xy=(2.1, max(o_E_K)), fontsize=34)
ax[1].annotate(r'para',  xy=(0.2, max(o_E_K)), fontsize=34)

# Plot the shared X title
fig.suptitle(r'Angular momentum quantum number $J$', x=0.5, y=0.04, fontsize=22, ha='center', va='center')

# remove the space between two subfigures
plt.subplots_adjust(wspace=0)

##### plotting arrows
### Lines to plot, read from the file
col_ends = (0, 1, 2, 3, 4, 5, 6, 7)
names = ['por', 'J_u', 'Ka_u', 'Kc_u', 'rm', 'J_l', 'Ka_l', 'Kc_l']
transitions_to_plot = ascii.read('transitions.txt', format='fixed_width_no_header', col_ends=col_ends, names=names)
transitions_to_plot.remove_column('rm')

#### extract the table for plotting the transtions
transition_table = join(table, transitions_to_plot, keys=['J_u', 'Ka_u', 'Kc_u', 'J_l', 'Ka_l', 'Kc_l'])

#for row in transition_table, where only the transitions are selected based on the "transition.txt"
texts_ortho = []
texts_para = []
for ind, row in enumerate(transition_table):
    x_up = row['J_u']
    x_lo = row['J_l']
    y_up = row['Energy_upper'] / k_B
    y_lo = row['Energy_lower'] / k_B
    freq = row['freq'] / 1000.0
    
    ax[int(row['por']=='p')].annotate('', xy=(x_lo, y_lo + 6), xytext=(x_up, y_up + 6), 
                    arrowprops=dict(facecolor='#9c9c9c', width=0.4, headwidth=6, edgecolor='black', lw=0.5), zorder=2)
    
    text = ax[int(row['por']=='p')].annotate(f'{freq:.1f}', xy=((x_up+x_lo)/2 - 1.1*int(row['por']=='p'), (y_up+y_lo)/2), color='#888888', fontsize=14)
    
    if row['por'] == 'o':
        texts_ortho.append(text)
    else:
        texts_para.append(text)

# Optimize the positions of text annotations to avoid overlapping
adjust_text(texts_ortho, ax=ax[0], force_text=(2, 0))
adjust_text(texts_para, ax=ax[1], force_text=(2, 0))

# Show the plot
plt.savefig('H2O_diag.png')

plt.show()
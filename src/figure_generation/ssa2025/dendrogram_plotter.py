
# # Load precalculated linkage matrix
# linkmat = np.load(LMD)


# # Create dendrogram rendering elements
# dend = dendrogram(linkmat,labels=df_cat[df_cat.CAT3].index.values,
#                   distance_sort=False, count_sort='ascending',
#                   color_threshold=1 - CT, above_threshold_color=at_color,
#                   no_plot=True, get_leaves=True)


### PLOTTING SECTION ###

### DENDROGRAMS ###
## DO DENDROGRAM UP HERE TO GET LEAF COLOR ASSIGNMENTS FOR GROUPS

# fig = plt.figure(figsize=(16, 8))
# ax = fig.add_subplot(111)
# for ic, dc, cc in zip(dend['icoord'], dend['dcoord'], dend['color_list']):
#     _ic = np.array(ic) - np.min(np.min(dend['icoord']))
#     _ic /= (np.max(np.max(dend['icoord'])) - np.min(np.min(dend['icoord'])))
#     _ic *= len(dend['leaves']) - 1
#     if cc == at_color:
#         ax.plot(_ic, dc, lw=0.5, color=cc)
#     else:
#         ax.plot(_ic, dc, lw=0.75, color=cc)

# ax.set_ylabel('Average Coherence Level')
# ax.set_yticks(np.arange(0,0.9, 0.1), [f'{_e:.1f}' for _e in np.arange(1.0,0.1, -0.1)])
# ax.set_ylim([0,ax.get_ylim()[1]])

# # Format X-axis
# ticks = []
# labels = []
# for _gn, _gc in df_cat.group.value_counts().items():
#     if _gc > 1:
#         _df = df_cat[df_cat.group==_gn]
#         _lx = _df.leafpos.mean()
#         _ll = int(_df.tidy_group.mean())
#         ticks.append(_lx)
#         labels.append(f'{_ll:d}')
# ax.set_xticks(ticks, labels=labels, fontsize=4)

# ax.set_xlim([-1, len(df_cat[df_cat.CAT3])])

# breakpoint()

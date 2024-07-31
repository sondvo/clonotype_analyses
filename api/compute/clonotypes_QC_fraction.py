import numpy as np
import pandas as pd


def _merge_metadata_fields(VDJ_10X, keys=['Condition']):
	if len(keys) == 1 and keys[0] in VDJ_10X.columns:
		return VDJ_10X, keys[0]

	meta_key = '<<>>'.join(keys)
	VDJ_10X[meta_key] = [
		'<<>>'.join(i)
		for i in VDJ_10X[keys].values.astype('str')
	]
	return VDJ_10X, meta_key


def _count_TRA_TRB(chain_names, chain_counts):
	if 'TRA' not in chain_names:
		n_TRA = 0
	else:
		n_TRA = chain_counts[chain_names=='TRA'][0]

	if 'TRB' not in chain_names:
		n_TRB = 0
	else:
		n_TRB = chain_counts[chain_names=='TRB'][0]

	return n_TRA, n_TRB


def _init_clonotypes_types(all_cells, meta_key):
	n_cells = len(all_cells)
	fraction_clo_dct = {
		'barcode': all_cells,
		'1 TRA - 1 TRB': np.zeros(n_cells, dtype=np.bool_),
		'0 TRA - 1 TRB': np.zeros(n_cells, dtype=np.bool_),
		'1 TRA - 0 TRB': np.zeros(n_cells, dtype=np.bool_),
		'2 TRA - 1 TRB': np.zeros(n_cells, dtype=np.bool_),
		'1 TRA - 2 TRB': np.zeros(n_cells, dtype=np.bool_),
		'2 TRA - 2 TRB': np.zeros(n_cells, dtype=np.bool_),
		'others': np.zeros(n_cells, dtype=np.bool_),
		meta_key: np.array([''] * n_cells),
	}

	return fraction_clo_dct


## TODO: optimize, remove for loop!
def create_clonotype_fraction_df(VDJ_10X, keys=['Condition']):
	VDJ_10X, meta_key = _merge_metadata_fields(VDJ_10X, keys=keys)
	all_cells = VDJ_10X['barcode'].unique()
	n_cells = len(all_cells)

	fraction_clo_dct = _init_clonotypes_types(all_cells, meta_key)
	for i in range(n_cells):
		bc = all_cells[i]
		tmp_df = VDJ_10X.iloc[
			VDJ_10X['barcode'].values == bc,
			:
		]

		all_chains = tmp_df['chain'].values
		chain_names, chain_counts = np.unique(all_chains, return_counts=True)
		n_TRA, n_TRB = _count_TRA_TRB(chain_names, chain_counts)

		tmp_key = '{} TRA - {} TRB'.format(n_TRA, n_TRB)
		if tmp_key not in fraction_clo_dct:
			fraction_clo_dct['others'][i] = True
		else:
			fraction_clo_dct[tmp_key][i] = True

		fraction_clo_dct[meta_key][i] = tmp_df[meta_key].values[0]

	fraction_clo_df = pd.DataFrame(fraction_clo_dct)
	fraction_clo_df = fraction_clo_df.set_index('barcode')

	return fraction_clo_df, meta_key


def grouping_ratio_clonotype_types(fraction_clo_df, meta_key):
	visualize_fraction_clo_df = fraction_clo_df.groupby(meta_key).sum()
	visualize_fraction_clo_df = visualize_fraction_clo_df.astype(np.float32)

	for i in range(len(visualize_fraction_clo_df)):
		visualize_fraction_clo_df.iloc[i, :] /= np.sum(visualize_fraction_clo_df.iloc[i, :].values)

	visualize_fraction_clo_df = visualize_fraction_clo_df.sort_index(ascending=False)
	return visualize_fraction_clo_df

def plotly_ratio_clonotype_types(fraction_clo_df, meta_key):
	fraction_clo_df = grouping_ratio_clonotype_types(fraction_clo_df, meta_key)
	reformated_dct = {
		'Groups': [],
		'Clonotypes types': [],
		'Ratio': []
	}
	for i in fraction_clo_df.index:
		for j in fraction_clo_df.columns:
			reformated_dct['Groups'].append(i)
			reformated_dct['Clonotypes types'].append(j)
			reformated_dct['Ratio'].append(
				fraction_clo_df.loc[i, j]
			)
	reformated_df = pd.DataFrame(reformated_dct)
	reformated_df = reformated_df.set_index('Groups')
	return reformated_df


def visualize_ratio_clonotype_types(fraction_clo_df, meta_key):
	visualize_fraction_clo_df = grouping_ratio_clonotype_types(fraction_clo_df, meta_key)

	ax = visualize_fraction_clo_df.plot(
		kind='bar',
		stacked=True,
		title='Ratio of clonotypes types',
		xlabel=''
	)
	ax.legend(
		loc='center left',
		bbox_to_anchor=(1.04, 0.5),
		fancybox=True,
		shadow=True,
		ncol=1
	)
	plt.show()
	return
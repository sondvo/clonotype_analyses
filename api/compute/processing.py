import numpy as np

from .. import common
from ..common.constants import VDJ_10X_COLUMNS


def _preprocessing_clonotypes(df):
	df['raw_clonotype_id'] = df['raw_clonotype_id'].fillna('N/A')
	df = df.iloc[
		df['raw_clonotype_id'].values != 'N/A',
		:
	]
	return df


def reformat_clonotypes(vdj_path, preprocessing=True):
	chosen_columns = [i.value for i in VDJ_10X_COLUMNS]

	vdj_df = common.read_csv(vdj_path)
	vdj_df.columns = [i.lower() for i in vdj_df.columns]

	try:
		vdj_df = vdj_df[chosen_columns]
	except Exception as e:
		raise (e)

	if preprocessing:
		return _preprocessing_clonotypes(vdj_df)

	return vdj_df


def matching_barcodes(vdj_df, meta_df):
	print ('WARNING: only keep intersect barcodes')
	vdj_original_bc = vdj_df[VDJ_10X_COLUMNS.BARCODES.value].values
	bc_idx = common.matching_barcodes_idx(
		vdj_original_bc,
		meta_df.index.values
	)
	matched_bc = np.unique(vdj_original_bc[bc_idx])
	vdj_df = vdj_df.iloc[
		bc_idx,
		:
	]
	meta_df = meta_df.loc[
		matched_bc,
		:
	]
	meta_df[VDJ_10X_COLUMNS.BARCODES.value] = meta_df.index.values
	return vdj_df, meta_df


def write_hdf5_vdj(vdj_df, output_h5_path):
	with common.H5AtomicWriter(output_h5_path) as f:
		for i in vdj_df.columns:
			try:
				f.create_dataset(i, data=vdj_df[i].values)
			except:
				f.create_dataset(i, data=vdj_df[i].values.astype('S'))
	return output_h5_path
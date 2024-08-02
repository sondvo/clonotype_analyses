import h5py
import numpy as np
import pandas as pd

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

	df = df.set_index(VDJ_10X_COLUMNS.BARCODES.value)

	if preprocessing:
		return _preprocessing_clonotypes(vdj_df)

	return vdj_df


def matching_barcodes(vdj_df, meta_df):
	print ('WARNING: only keep intersect barcodes')
	bc_idx = common.matching_barcodes_idx(
		vdj_df.index.values,
		meta_df.index.values
	)
	matched_bc = np.unique(vdj_df.index.values[bc_idx])
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


def h5_to_pandas(h5_path, columns):
	meta_dct = {}
	with h5py.File(h5_path) as f:
		bc = f[VDJ_10X_COLUMNS.BARCODES.value][:].astype('str')
		for i in columns:
			if i == VDJ_10X_COLUMNS.BARCODES.value:
				continue
			arr = f[i][:]
			try:
				arr = arr.astype(np.float32)
			except:
				arr = arr.astype('str')
			meta_dct[i] = arr
	df = pd.DataFrame(meta_dct, index=bc)
	return df


def merge_vdj_and_clinical_meta(vdj_df, clinical_df):
	if not np.all(
		vdj_df.index.unique() == clinical_df.index.unique()
	):
		raise Exception('Barcodes of 2 dataframes must matched')

	vdj_df[clinical_df.columns.values] = clinical_df.loc[
		vdj_df.index.values,
		clinical_df.columns.values
	]

	return vdj_df
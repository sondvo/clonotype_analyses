import os
import json
import h5py
import uuid
import numpy as np
import pandas as pd

from .constants import VDJ_10X_COLUMNS


def mkdir(path):
	if not os.path.isdir(path):
		os.system(f'mkdir -p "{path}"')


def isfile(path):
	 os.path.isfile(path)


def isdir(path):
	return os.path.isdir(path)


def list_dir(path):
	return os.listdir(path)


def list_dir_fullpath(path):
	return [
		os.path.join(path, x)
		for x in os.listdir(path)
		if os.path.isdir(os.path.join(path, x))
	]


def write_json(o, path):
	with open(path, "w", encoding="utf-8") as f:
		json.dump(o, f)


def read_json(path):
	if not os.path.isfile(path):
		return {}

	with open(path, "r", encoding="utf-8") as f:
		data = json.load(f)

	return data


def read_csv(path, **kwargs):
	df = pd.read_csv(filepath_or_buffer=path, sep="\t", **kwargs)

	if "index_col" in kwargs:
		if len(df.columns) == 0:
			return pd.read_csv(filepath_or_buffer=path, sep=",", **kwargs)
	else:
		if len(df.columns) < 2:
			return pd.read_csv(filepath_or_buffer=path, sep=",", **kwargs)
	return df


def write_csv(path, content, **kwargs):
	df = pd.DataFrame(content)
	df.to_csv(path_or_buff=path, **kwargs)
	return True


def join_path(*args):
	if len(args) == 0:
		return ""
	return os.path.join(*args)


class H5AtomicWriter(h5py.File):
	def __init__(self, path):
		self._path = path
		self._temp_path = f"{path}.TEMP{uuid.uuid4().hex}"
		super().__init__(self._temp_path, mode="w")

	def close(self):
		super().close()
		os.rename(self._temp_path, self._path)

	def __exit__(self, exc_type, exc_value, traceback):
		if exc_type is None and exc_value is None and traceback is None:
			self.close()
		else:
			super().close()
			os.remove(self._temp_path)


# TODO: hashtable -> Optimize memory
def matching_barcodes_idx(arr1, arr2):
	tmp_dct = {i: j for j, i in enumerate(arr2)}
	chosen_idx_arr1 = np.array([True if i in tmp_dct else False for i in arr1])
	return chosen_idx_arr1


def get_file_name(file_path):
	return file_path.split('/')[-1]


def extract_h5(h5_path, key):
	with h5py.File(h5_path) as f:
		if key not in h5_path:
			raise Exception('cannot find {} in h5'.format(key))
		res = f[key][:]
	return res


def concat_df(df_list, **kwargs):
	final_df = None
	for df in df_list:
		if final_df is None:
			final_df = df
		else:
			final_df = pd.concat([final_df, df], **kwargs)
	return final_df
import pandas as pd
from typing import List

from . import common
from .compute import processing


class ClonotypePreprocessing(object):
	def __init__(
		self,
		vdj_path: List[str],
		ouput_h5_path: str,
		preprocessing: bool = True,
		add_prefix: bool = True,
		prefixes: List[str] = None,
	):
		self.__vdj_path = vdj_path
		self.__ouput_h5_path = ouput_h5_path
		self.__preprocessing = preprocessing
		self.__add_prefix = add_prefix
		self.__prefixes = prefixes


	def _write_h5_vdj(
			self,
			vdj_df
		):
		processing.write_hdf5_vdj(vdj_df, self.__ouput_h5_path)
		return self.__ouput_h5_path


	def processing_vdj(self):
		DF = None
		for i in self.__vdj_path:
			df_path = self.__vdj_path[i]
			df = processing.reformat_clonotypes(df_path, self.__preprocessing)

			if self.__add_prefix:
				if self.__prefixes is None:
					df_file_name = common.get_file_name(df_path)
					prefix = df_file_name.replace('_contig_annotations.csv', '')
				else:
					prefix = self.__prefixes[i]
				df.index = prefix + '_' + df.index.values

			if DF is None:
				DF = df
			else:
				DF = pd.concat([DF, df], axis=0)

		return DF


	def merge_with_clinical_meta(
			self,
			reformated_vdj_df,
			input_clincal_path: str,
			output_clincal_path: str,
		):
		vdj_df, clinical_df = processing.merge_with_clinical_meta(
			reformated_vdj_df,
			input_clincal_path,
			output_clincal_path
		)
		common.write_csv(clinical_df, output_clincal_path)
		return vdj_df, clinical_df


	def ingest_data(
			self,
			input_clincal_path: str = None,
			output_clincal_path: str = None,
		):
		vdj_df = self.processing_vdj()
		if input_clincal_path is not None:
			vdj_df, _ = self.merge_with_clinical_meta(
				vdj_df,
				input_clincal_path,
				output_clincal_path
			)
		h5_path = self._write_h5_vdj(vdj_df)
		return {
			'10X_VDJ': h5_path,
			'clinical_meta': output_clincal_path
		}


class ClonotypeToolkits(object):
	def __init__(
		self,
		vdj_h5_path: str,
		clinical_meta_path: str = None,
	):
		"""
		Parameters
		----------
		vdj_path : List[str]
			The path to ORIGINAL 10X VDJ files (could be multiple files for multiple samples)
		clinical_meta_path : str
			The path to associate metadata file (eg. include 'condition', 'cell types', ...)\n
			First column must be barcodes\n
			First row must be field names
		"""
		self.__vdj_h5_path = vdj_h5_path
		self.__clinical_meta = None
		if clinical_meta_path is not None:
			self.__clinical_meta = common.read_csv(clinical_meta_path, index_col=0)

import pandas as pd
from typing import List
from typing import Dict

from .compute import processing
from .compute import clonotypes_tracing
from .compute import clonotypes_expansion
from .compute import clonotypes_QC_fraction
from .compute import clonotypes_diversity_rate

from . import common
from .common.constants import VDJ_10X_COLUMNS


class ClonotypePreprocessing(object):
	def __init__(
		self,
		batch_info: List[Dict],
		ouput_h5_path: str,
		preprocessing: bool = True,
	):
		"""
		Create object for preprocessing clonotypes data

		Parameters
		----------
		batch_info : List[dict]
			Allow merging multiple VDJ samples\n
			For each sample, barcodes of vdj_df and clinical_meta_df MUST have similar items \n
			Example:
				[{
					'sample_name': 'sample_1',\n
					'vdj_path': 'root/fol/GSM123_vdj_t_filtered_contig_annotations.csv',\n
					'clinical_meta_path': 'root/fol/GSM123_clinical_meta.tsv',
				}, {...}]
		"""
		self.__batch_info = batch_info
		self.__ouput_h5_path = ouput_h5_path
		self.__preprocessing = preprocessing


	def processing_vdj(self):
		"""
		Only keep necessary columns in VDJ file\n
		Match barcodes between VDJ and clinical metadata
		"""
		final_vdj_df = None
		final_clinical_df = None
		for info in self.__batch_info:
			prefix = info.get('sample_name', '')
			vdj_path = info['vdj_path']
			clinical_path = info['clinical_meta_path']

			vdj_df = processing.reformat_clonotypes(vdj_path, self.__preprocessing)
			clinical_df = common.read_csv(clinical_path, index_col=0)

			if len(prefix):
				vdj_df[
					VDJ_10X_COLUMNS.BARCODES.value
				] = prefix + '_' + vdj_df[VDJ_10X_COLUMNS.BARCODES.value].values
				vdj_df['sample_name'] = prefix

				clinical_df.index = prefix + '_' + clinical_df.index.values
				clinical_df['sample_name'] = prefix

			final_vdj_df = common.concat_df([final_vdj_df, vdj_df], axis=0)
			final_clinical_df = common.concat_df([final_clinical_df, clinical_df], axis=0)
		return final_vdj_df, final_clinical_df


	def ingest_data(
			self,
			output_clincal_path: str = None,
		):
		"""
		Save VDJ as h5 format\n
		Keep clinical metadata as pandas dataframe
		"""
		vdj_df, clinical_df = self.processing_vdj()
		vdj_df, clinical_df = processing.matching_barcodes(
			vdj_df, clinical_df
		)
		processing.write_hdf5_vdj(vdj_df, self.__ouput_h5_path)
		common.write_csv(clinical_df, output_clincal_path)
		return {
			'10X_VDJ': self.__ouput_h5_path,
			'clinical_meta': output_clincal_path
		}


class ClonotypeToolkits(object):
	def __init__(
		self,
		vdj_h5_path: str,
		clinical_meta_path: str,
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

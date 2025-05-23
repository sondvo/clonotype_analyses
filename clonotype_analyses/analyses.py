import pandas as pd
from typing import List
from typing import Dict

from .compute import processing
from .compute import clonotypes_tracing
from .compute import clonotypes_expansion
from .compute import clonotypes_QC_fraction
from .compute import clonotypes_diversity_rate

from . import common
from .common.constants import SAMPLE_NAME
from .common.constants import VDJ_10X_COLUMNS


class ClonotypePreprocessing(object):
	def __init__(
		self,
		batch_info: List[Dict],
	):
		"""
		Create object for preprocessing clonotypes data

		Parameters
		----------
		batch_info : List[dict]
			Allow merging multiple VDJ samples\n
			For each sample: barcodes of vdj_df and clinical_meta_df MUST have similar items \n
			Example:
				[{
					'sample_name': 'sample_1',\n
					'vdj_path': 'root/fol/GSM123_vdj_t_filtered_contig_annotations.csv',\n
					'clinical_meta_path': 'root/fol/GSM123_clinical_meta.tsv',
				}, {...}]
		"""
		print ('NOTE: barcode (first column) of vdj and clinical_meta MUST have similar items')
		self.__batch_info = batch_info
		self.__sample_specific_columns = [
			VDJ_10X_COLUMNS.RAW_CLONOTYPE_ID.value,
			VDJ_10X_COLUMNS.RAW_CONSENSUS_ID.value
		]


	def _processing_vdj(self, preprocessing=True):
		"""
		Only keep necessary columns in VDJ file
		"""
		final_vdj_df = None
		final_clinical_df = None
		for info in self.__batch_info:
			prefix = info.get('sample_name', '')
			vdj_path = info['vdj_path']
			clinical_path = info['clinical_meta_path']

			vdj_df = processing.reformat_clonotypes(
				vdj_path, preprocessing
			)
			clinical_df = common.read_csv(
				clinical_path, index_col=0
			)

			if len(prefix):
				for col in self.__sample_specific_columns:
					vdj_df[col] = prefix + '_' + vdj_df[col].values
				vdj_df.index = prefix + '_' + vdj_df.index.values
				vdj_df[SAMPLE_NAME] = prefix

				clinical_df.index = prefix + '_' + clinical_df.index.values
				clinical_df[SAMPLE_NAME] = prefix

			final_vdj_df = common.concat_df(
				[final_vdj_df, vdj_df],
				axis=0
			)
			final_clinical_df = common.concat_df(
				[final_clinical_df, clinical_df],
				axis=0
			)
		return final_vdj_df, final_clinical_df


	# For clonotypes, each function require some spefic columns, not all of them.
	# Eg: shannon entropy: [barcodes, raw_clonotype_id]
	# -> save each column as key in hdf5 -> optimize data on memory
	##
	# Why not setting barcodes as the outer-most layer?
	# Eg: f['ATGCAT...'].keys() = ['TRA', 'TRB', 'cdr3', ...] ?
	# -> very few functions require picking some specific barcodes -> not optimize for analyses.
	def ingest_data(
			self,
			ouput_vdj_h5_path: str,
			ouput_clinical_h5_path: str,
			preprocessing: bool = True,
			intersect_barcodes: bool = False,
		):
		"""
		Save VDJ as h5 format\n
		Keep clinical metadata as pandas dataframe

		Parameters
		----------
		ouput_vdj_h5_path : str
		ouput_clinical_h5_path : str
		preprocessing : boolean
		intersect_barcodes : boolean

		Returns
		----------
		Saving location: Dict
			In which:
				'10X_VDJ': location of VDJ h5\n
				'clinical_meta': location of clinical meta h5
		"""
		vdj_df, clinical_df = self._processing_vdj(preprocessing)

		if intersect_barcodes:
			vdj_df, clinical_df = processing.matching_barcodes(
				vdj_df, clinical_df
			)

		processing.store_df_as_h5(vdj_df, ouput_vdj_h5_path)
		processing.store_df_as_h5(clinical_df, ouput_clinical_h5_path)
		return {
			'10X_VDJ': ouput_vdj_h5_path,
			'clinical_meta': ouput_clinical_h5_path
		}


class ClonotypeToolkits(object):
	def __init__(
		self,
		vdj_h5_path: str,
		clinical_h5_path: str,
	):
		"""
		Load necessary files for analyses\n

		Parameters
		----------
		vdj_h5_path : str
			vdj_h5_path output from ClonotypePreprocessing.ingest_data()
		clinical_h5_path : str
			clinical_h5_path output from ClonotypePreprocessing.ingest_data()
		"""
		print ("NOTE: vdj_h5_path, clinical_h5_path are the outputs from ClonotypePreprocessing.ingest_data()")
		self.__vdj_h5_path = vdj_h5_path
		self.__clinical_h5_path = clinical_h5_path


	def _prepare_clonotypes_QC_fraction(self, meta_keys):
		if not len(meta_keys):
			raise Exception('Need at least one meta column to compare')

		vdj_df = processing.h5_to_pandas(
			self.__vdj_h5_path, [VDJ_10X_COLUMNS.CHAIN.value]
		)
		clinical_df = processing.h5_to_pandas(
			self.__clinical_h5_path, meta_keys
		)
		vdj_df, clinical_df = processing.matching_barcodes(
			vdj_df, clinical_df
		)
		merged_df = processing.merge_vdj_and_clinical_meta(
			vdj_df, clinical_df
		)
		return clonotypes_QC_fraction.create_clonotype_fraction_df(
			merged_df, meta_keys
		)


	def matplotlib_clonotypes_QC_fraction(self, meta_keys=[]):
		"""
		Visualize fraction of paired / ambiguous clonotypes separated by metadata fields.

		Parameters
		----------
		meta_keys : List[str]
			List of metadata keys for grouping \n
			If more than one keys provided, the values will be combined together \n
			Eg: ['condition', 'patient']

		Returns
		----------
		Matplotlib image
		"""
		fraction_clo_df, merged_meta_key = self._prepare_clonotypes_QC_fraction(meta_keys)
		return clonotypes_QC_fraction.visualize_ratio_clonotype_types(
			fraction_clo_df, merged_meta_key
		)


	def plotly_clonotypes_QC_fraction(self, meta_keys=[]):
		"""
		Produce dataframe for plotly: fraction of paired / ambiguous clonotypes separated by metadata fields.

		Parameters
		----------
		meta_keys : List[str]
			List of metadata keys for grouping \n
			If more than one keys provided, the values will be combined together \n
			Eg: ['condition', 'patient']

		Returns
		----------
		Pandas dataframe using for plotly
		"""
		fraction_clo_df, merged_meta_key = self._prepare_clonotypes_QC_fraction(meta_keys)
		return clonotypes_QC_fraction.plotly_ratio_clonotype_types(
			fraction_clo_df, merged_meta_key
		)




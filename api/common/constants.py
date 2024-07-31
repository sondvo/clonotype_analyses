from enum import Enum


class VDJ_10X_COLUMNS(Enum):
	BARCODES='barcodes'
	IS_CELL='is_cell'
	HIGH_CONFIDENCE='high_confidence'
	LENGTH='length'
	CHAIN='chain'
	V_GENE='v_gene'
	D_GENE='d_gene'
	J_GENE='j_gene'
	C_GENE='c_gene'
	FULL_LENGTH='full_length'
	PRODUCTIVE='productive'
	CDR3='cdr3'
	UMIS='umis'
	RAW_CLONOTYPE_ID='raw_clonotype_id'
	RAW_CONSENSUS_ID='raw_consensus_id'
	EXACT_SUBCLONOTYPE_ID='exact_subclonotype_id'


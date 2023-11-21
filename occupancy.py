import gzip, shutil, subprocess
import os, csv, time
from functools import wraps
from cell_types import CELL_TYPES

import matplotlib.pyplot as plt


# Column            	Example             	Description
# ------            	-------             	-----------
# 1.  seqname       	chr10               	Chromosome
# 2.  start         	97320044            	Start position (0-based)
# 3.  end           	97320056            	End position
# 4.  identifier    	10.754379.4         	Unique identifier (DHS_chom#.DHS_position%.fp_position%; DHS_chrom#.DHS_position% = DHS index identifier)
# 5.  mean_signal   	55.865317           	Mean footprint -log(1-posterior) across biosamples ("confidence score")
# 6.  num_samples   	9                   	Number of unique biosamples contributing to this index footprint (posterior >= 0.99)
# 7.  num_fps       	9                   	Number of unique footprints across all samples that contributed to this consensus footprint 
# 8.  width         	12                  	Width of consensus footprint (column 3-column 2)
# 9.  summit        	97320049            	Estimated footprint summit position
# 10. core_start    	97320042            	Start position of core-region containing 95% of per-biosample summits
# 11. core_end      	97320053            	End position of core-region containing 95% of per-biosample summits
# 12. motif_clusters	CTCF;KLF/SP/2;ZNF563	Overlapping motif clusters (delimited with semi-colon)

CONSENSUS_LEGEND = [
	'seqname', 'start', 'end', 'identifier', 'mean_signal',
	'num_samples', 'num_fps', 'width', 'summit',
	'core_start', 'core_end', 'motif_clusters'
]

# For writing occupancies we only need a subset of that data
OCCUPANCIES_LEGEND = [
	'seqname', 'start', 'end', 'identifier', 'occupancy'
]


# Per-nucleotide data
#  	Column	Example	Description
# 1	contig	chr1	Chromosome
# 2	start	9823494	Start position (0-based)
# 3	stop	9823495	End position (start+1)
# 4	obs	24	Observed DNase I cleavages (both strands combined)
# 5	exp	56	Expected DNase I cleavages
# 6	lnp	2.1	–log p-value lower-tail negative binomial
# 7	winlp	14.5	–log p-value windowed test (Stouffer’s Z)
# 8	fdr	0.0014	Empircal false-discovery rate

# NOTE: the official "observed" and "expected" columns  in the data from Vierstra 
# are swapped relative to what they should be. This is fixed in the nucleotide 
# legend below
NUCLEOTIDE_LEGEND = [
	'contig', 'start', 'stop', 'exp',
	'obs', 'lnp', 'winlp', 'fdr'
]

# Helper for logging the time each function takes
def timed(func):
	@wraps(func)
	def timed_wrapper(*args, **kwargs):
		print(f'Starting {func.__name__}')
		start_time = time.perf_counter()
		result = func(*args, **kwargs)
		end_time = time.perf_counter()
		total_time = end_time - start_time
		print(f'Function {func.__name__} took {total_time:.4f} seconds')
		return result
	return timed_wrapper

def gzip_and_clean(out_path):
	with open(out_path, 'rb') as f_in:
		gz_path = f'{out_path}.gz'
		with gzip.open(gz_path, 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
	os.remove(out_path)

@timed
def write_consensus_footprints_for_cell_type(data_path, consensus_data_path, cell_type):
	cell_type_index = CELL_TYPES.index(cell_type)
	matrix_path = os.path.join(data_path, 'consensus.index', 'consensus_index_matrix_binary_hg38.txt.gz')
	cell_type_dir = os.path.join(data_path, 'per.dataset', cell_type)
	if not os.path.exists(cell_type_dir):
		os.makedirs(cell_type_dir)
	out_path = os.path.join(cell_type_dir, 'only_consensus.bed')
	with gzip.open(matrix_path, 'rb') as matrix, gzip.open(consensus_data_path) as consensus_data:
		with open(out_path, 'w', newline='', encoding='utf-8') as tsv_file:
			tsv_writer = csv.writer(tsv_file, delimiter='\t')
			# We can lean on the fact that both the matrix and the consensus data are sorted and reference
			# the same set of consensus identifiers in the same order
			for i, (matrix_row, consensus_row) in enumerate(zip(matrix, consensus_data)):
				vals = matrix_row.strip().split()
				# Each row consists of the footprint identifier followed by a 0 or 1 for every cell type
				matrix_identifier = vals[0]
				present_in_cell_type = bool(int(vals[cell_type_index + 1]))
				if present_in_cell_type:
					footprint = {col: val.decode() for (col, val) in zip(CONSENSUS_LEGEND, consensus_row.strip().split()) if col in OCCUPANCIES_LEGEND}
					tsv_writer.writerow(footprint.values())
	gzip_and_clean(out_path)

@timed
def get_nucleotides(data_path, cell_type):
	consensus_data_path = os.path.join(data_path, 'per.dataset', cell_type, 'only_consensus.bed.gz')
	base_url = "https://resources.altius.org/~jvierstra/projects/footprinting.2020"
	retrieve_data_url = f'{base_url}/per.dataset/{cell_type}/interval.all.bedgraph.gz'
	output_data_file = os.path.join(data_path, 'per.dataset', cell_type, 'interval.only_consensus_loci.gz')
	if not os.path.exists(output_data_file):
		command = f'tabix --targets {consensus_data_path} {retrieve_data_url} | bgzip > {output_data_file}'
		os.system(command)
		command = f'tabix -0 {output_data_file}'
		os.system(command)
	else:
		print(f'Nucleotide data already exists for {cell_type}')


@timed
def write_occupancies(data_path, cell_type):
	# This assumes that both the consensus data and the nucleotide data is sorted
	# And furthermore that the nucleotide data consists of a set of loci that is 
	# a subset of the loci covered by consensus footprints
	nucleotide_data_path = os.path.join(data_path, 'per.dataset', cell_type, 'interval.only_consensus_loci.gz')
	consensus_data_path = os.path.join(data_path, 'per.dataset', cell_type, 'only_consensus.bed.gz')
	out_path = os.path.join(data_path, 'per.dataset', cell_type, 'occupancy.bed')
	with gzip.open(consensus_data_path) as consensus_data, gzip.open(nucleotide_data_path) as nucleotide_data:
		with open(out_path, 'w', newline='', encoding='utf-8') as tsv_file:
			tsv_writer = csv.writer(tsv_file, delimiter='\t')
			next_nuc = {col: val for (col, val) in zip(NUCLEOTIDE_LEGEND, next(nucleotide_data).strip().split())}
			for i, line in enumerate(consensus_data):
				if i % 200_000 == 0:
					print(f'Processed {i} footprints')
				footprint = {col: val for (col, val) in zip(CONSENSUS_LEGEND, line.strip().split())}
				start = int(footprint['start'])
				end = int(footprint['end'])
				total_expected = 0
				total_observed = 0
				for pos in range(start, end):
					if pos == int(next_nuc['start']):
						total_expected += float(next_nuc['exp'])
						total_observed += float(next_nuc['obs'])
						try:
							next_nuc = {col: val for (col, val) in zip(NUCLEOTIDE_LEGEND, next(nucleotide_data).strip().split())}
						except StopIteration:
							continue
				if total_observed > 0 and total_expected > 0:
					occupancy = 1 - (total_observed / total_expected)
					# Many occupancies are <0 or >1. Round those off.
					occupancy = max(occupancy, 0)
					occupancy = min(occupancy, 1)
					footprint['occupancy'] = str(occupancy).encode('utf8')
					footprint = {key: val.decode() for (key, val) in footprint.items() if key in OCCUPANCIES_LEGEND}
					tsv_writer.writerow(footprint.values())
	# gzip file and delete unzipped version
	gzip_and_clean(out_path)
	os.remove(consensus_data_path)

def occupancy_path(data_path, cell_type):
	return os.path.join(data_path, 'per.dataset', cell_type, 'occupancy.bed.gz')

def occupancy_already_written(data_path, cell_type):
	return os.path.exists(occupancy_path(data_path, cell_type))

def get_occupancies(data_path, cell_type):
	occupancies = []
	path = occupancy_path(data_path, cell_type)
	with gzip.open(path, 'rb') as f:
		for line in f:
			data = {col: val for (col, val) in zip(OCCUPANCIES_LEGEND, line.strip().split())}
			occupancies.append(float(data['occupancy']))
	return occupancies

def plot_occupancies(data_path, cell_types):
	for cell_type in cell_types:
		occupancies = get_occupancies(data_path, cell_type)
		plt.hist(occupancies, bins=50, histtype="step", label=cell_type, density=True)
	plt.legend()
	plt.xlabel("occupancy")
	plt.ylabel("probability density")
	plt.title(f'Probability distribution of footprint occupancy by cell type')
	path = os.path.join(data_path, 'all_occupancies_histogram.png')
	plt.savefig(path)
	plt.close()

import occupancy
import time, os, argparse
from cell_types import CELL_TYPES


if __name__ == "__main__":
	print("Starting script")
	start_time = time.time()
	
	arg_parser = argparse.ArgumentParser()
	arg_parser.add_argument("-A", "--all", action='store_true', help="Check to rebuild all occupancy files")
	arg_parser.add_argument("-n", "--num", help="Maximum number of cell types to populate")
	arg_parser.add_argument("-p", "--path", help="Base path for writing data")

	args = arg_parser.parse_args()
	import_all = bool(args.all)
	num_cell_types = args.num
	if num_cell_types is not None:
		num_cell_types = int(num_cell_types)
	data_path = args.path
	if data_path is None:
		data_path = os.path.join('data')

	# We assume that all files under data_path are organized according to the same
	# folder schema as in the Vierstra online repository
	consensus_data_path = os.path.join(data_path, 'consensus.index', 'consensus_footprints_and_motifs_hg38.bed.gz')

	cell_types_imported = 0
	for cell_type in CELL_TYPES:
		if num_cell_types and cell_types_imported >= num_cell_types:
			break
		if import_all or not occupancy.occupancy_already_written(data_path, cell_type):
			print(f'Populating {cell_type}')
			occupancy.write_consensus_footprints_for_cell_type(data_path, consensus_data_path, cell_type)
			occupancy.get_nucleotides(data_path, cell_type)
			occupancy.write_occupancies(data_path, cell_type)
			cell_types_imported += 1

	end_time = time.time()
	print(f'Completed in {end_time - start_time} s')

import occupancy
import time, os, argparse
from cell_types import CELL_TYPES


if __name__ == "__main__":
	print("Starting script")
	start_time = time.time()

	arg_parser = argparse.ArgumentParser()
	arg_parser.add_argument("-p", "--path", help="Base path for writing data")
	args = arg_parser.parse_args()
	data_path = args.path
	if data_path is None:
		data_path = os.path.join('data')

	cell_types = []
	for cell_dir in os.listdir(os.path.join(data_path, 'per.dataset')):
		if cell_dir in CELL_TYPES:
			cell_types.append(cell_dir)
	
	occupancy.plot_occupancies(data_path, cell_types)
	end_time = time.time()
	print(f'Completed in {end_time - start_time} s')
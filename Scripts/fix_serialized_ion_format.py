import ast
import csv
import glob
import json
import os
import sys

serialized_columns = ["Oxonium_ions", "bare_b_ions", "bare_y_ions", "b_ions_with_HexNAc", \
					"y_ions_with_HexNAc", "b_ion_coverage", "y_ion_coverage", "Stub_ions"]

def main(infile, outfile = None):
	if outfile is None:
		outfile = os.path.splitext(infile)[0] + ".fix.csv"
	f = open(infile,'rb')
	reader = csv.DictReader(f)
	data = [r for r in reader]
	for row in data:
		for s in serialized_columns:
			row[s] = json.dumps(ast.literal_eval(row[s]))

	f.close()
	f = open(outfile, 'wb')
	writer = csv.DictWriter(f, reader.fieldnames)
	writer.writeheader()
	writer.writerows(data)
	return outfile

if __name__ == '__main__':
	for f in sys.argv[1:]:
		for g in glob.glob(f):
			main(g)


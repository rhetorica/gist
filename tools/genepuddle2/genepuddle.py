# inputs: dir of ffns, tsv associating each ffn with a quantity, template PAR, output project path (opp)
# outputs: .par files in <opp>/pars, .gtf files in <opp>/pars, and .fasta files in <opp>/genes
#          run script suitable for feeding to GNU parallel

from __future__ import print_function
import sys, os

if 5 == len(sys.argv):
	input_dir = sys.argv[1].strip()
	manifest = sys.argv[2].strip()
	template_PAR = sys.argv[3].strip()
	opp = sys.argv[4].strip()
else:
	print("Incorrect parameter format. RTFS.", file=sys.stderr)
	sys.exit(1)

#create opp
#create opp/genes
#create opp/pars

gene_dir = opp + '/genes'
par_dir = opp + '/pars'

if not os.path.exists(gene_dir):
	os.makedirs(gene_dir)

if not os.path.exists(par_dir):
	os.makedirs(par_dir)

with open(template_PAR, 'r') as tempfi:
	template = ''.join(tempfi.readlines())

with open(opp + '/execute.sh', 'w') as exefo:
	with open(manifest, 'r') as manifi:
		for r in manifi:
			line = r.strip().split('\t')
			genome_name = line[0]
			reads = int(line[1])
	
			read_name = ""
			read_sequence = ""
			fasta_output = ""
			i = 0

			with open(par_dir + '/' + genome_name + '.par', 'w') as pfo:
				print(template.replace('number', str(reads)).replace('genome', genome_name), file=pfo)

			print('flux-simulator -x -l -s -p ' + par_dir + '/' + genome_name + '.par', file=exefo);

			with open(input_dir + '/' + genome_name, 'r') as genfi:
				with open(par_dir + '/' + genome_name + '.gtf', 'w') as gtfo:
					for gr in genfi:
						gline = gr.strip()
						if gline[0] == ">":
							i = i + 1
							if read_name != '' and read_sequence != "":
								print(genome_name + "\t" \
									  + "genepuddle\t" \
								      + "exon\t" \
								      + str(len(fasta_output) + 1) + "\t" \
								      + str(len(fasta_output) + len(read_sequence)) + "\t" \
								      + ".\t+\t0\t" \
								      + 'gene_name \"' + read_name + '\"; '\
									  + 'gene_id \"gene_' + str(i) + '\"; ' \
									  + 'transcript_id \"gene_' + str(i) + '\"', file=gtfo)
								fasta_output = fasta_output + read_sequence
	
							read_name = gline[1:].replace(" ", "_")
							read_sequence = ""
						else:
							read_sequence = read_sequence + gline
	
					i = i + 1

					print("condensed_" + genome_name + "\t" \
						  + "genepuddle\t" \
					      + "exon\t" \
					      + str(len(fasta_output) + 1) + "\t" \
					      + str(len(fasta_output) + len(read_sequence)) + "\t" \
					      + ".\t+\t0\t" \
					      + 'gene_name \"' + read_name + '\"; '\
						  + 'gene_id \"gene_' + str(i) + '\"; ' \
						  + 'transcript_id \"gene_' + str(i) + '\"', file=gtfo)
					fasta_output = fasta_output + read_sequence
	
			with open(gene_dir + '/' + genome_name + '.fa', 'w') as gfo:
				fasta_output = ">condensed_" + genome_name + "\n" + fasta_output
				print(fasta_output, file=gfo)
	
			print(genome_name)

# todo: generate shell script

print("Done.")

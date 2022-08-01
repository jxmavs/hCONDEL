hybrid_genome_file_name=$1
hybrid_genome_size_file_name=$2

faSize ${hybrid_genome_file_name} -detailed > ${hybrid_genome_size_file_name}

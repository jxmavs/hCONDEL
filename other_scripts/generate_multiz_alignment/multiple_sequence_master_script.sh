
#http://genomewiki.ucsc.edu/index.php/Hg38_100-way_conservation_lastz_parameters
chain_param_dir=/cluster_path/ape_project/deletions_project/msa/chainParamFiles 
lastz_param_dir=/cluster_path/ape_project/deletions_project/msa/lastzParamFiles

#look at this for parameters used http://genomewiki.ucsc.edu/index.php/Hg38_100-way_conservation_lastz_parameters
 
#set the target species/assembly name
#genomes_of_interest.txt is a file that contains the query genomes of interest for the msa
# multiple_sequence_alignment.sh is in /cluster_path/ape_project/deletions_project/msa
TARGET=chimp
TARGETASSEMBNAME=panTro4
TDIR=/cluster_path/ape_project/deletions_project/reference_genome_files/${TARGET}

TARGETDIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way

#first download the target data/sort the files if not already done
wget http://hgdownload.cse.ucsc.edu/goldenPath/${refGenome}/bigZips/${TARGETASSEMBNAME}.2bit -O /cluster_path/ape_project/deletions_project/reference_genome_files/${species}/${TARGETASSEMBNAME.2bit
wget http://hgdownload.cse.ucsc.edu/goldenPath/${refGenome}/bigZips/${TARGETASSEMBNAME}.chrom.sizes -O /cluster_path/ape_project/deletions_project/reference_genome_files/${species}/${TARGETASSEMBNAME}.chrom.sizes
	
sort -k1,1 /cluster_path/ape_project/deletions_project/reference_genome_files/${TARGET}/${TARGETASSEMBNAME}.chrom.sizes > /cluster_path/ape_project/deletions_project/reference_genome_files/${TARGET}/${TARGETASSEMBNAME}.chrom.sizes.sorted
mv /cluster_path/ape_project/deletions_project/reference_genome_files/${TARGET}/${TARGETASSEMBNAME}.chrom.sizes.sorted /cluster_path/ape_project/deletions_project/reference_genome_files/${TARGET}/${TARGETASSEMBNAME}.chrom.sizes
	
mkdir -p ${TARGETDIR}

#this file stores the species used 
cat <<_EOF_ > ${TARGETDIR}/genomes_of_interest.txt
macaque	rheMac8	${lastz_param_dir}/primateLastzParamFile.txt 	${chain_param_dir}/primateChainParamFile.txt 
bonobo	panPan1	${lastz_param_dir}/primateLastzParamFile.txt 	${chain_param_dir}/primateChainParamFile.txt 
gorilla	gorGor4	${lastz_param_dir}/primateLastzParamFile.txt 	${chain_param_dir}/primateChainParamFile.txt 
orangutan	ponAbe2	${lastz_param_dir}/primateLastzParamFile.txt 	${chain_param_dir}/primateChainParamFile.txt 
mouse	mm10	${lastz_param_dir}/defaultLastzParamFile.txt	${chain_param_dir}/chainMediumParams.txt
dog	canFam3	${lastz_param_dir}/defaultLastzParamFile.txt	${chain_param_dir}/chainMediumParams.txt
cow	bosTau8	${lastz_param_dir}/defaultLastzParamFile2.txt	${chain_param_dir}/chainMediumParams.txt 
opposum	monDom5	${lastz_param_dir}/HoxD55LastzParamFile.txt	${chain_param_dir}/chainFarParams.txt
chicken	galGal4	${lastz_param_dir}/HoxD55LastzParamFile.txt	${chain_param_dir}/chainFarParams.txt
platypus	ornAna2	${lastz_param_dir}/HoxD55LastzParamFile.txt	${chain_param_dir}/chainFarParams.txt
_EOF_

#download necessary files first 2bit and chrom.sizes
while read -a array
do
	species=${array[0]}
	refGenome=${array[1]}

	mkdir -p /cluster_path/ape_project/deletions_project/reference_genome_files/${species}
	wget http://hgdownload.cse.ucsc.edu/goldenPath/${refGenome}/bigZips/${refGenome}.2bit -O /cluster_path/ape_project/deletions_project/reference_genome_files/${species}/${refGenome}.2bit
	wget http://hgdownload.cse.ucsc.edu/goldenPath/${refGenome}/bigZips/${refGenome}.chrom.sizes -O /cluster_path/ape_project/deletions_project/reference_genome_files/${species}/${refGenome}.chrom.sizes
	
done< /cluster_path/ape_project/deletions_project/msa/chimp_10_way/genomes_of_interest.txt

#first loop generates lastz script and chain script
while read -a array
do
	species=${array[0]}
	refGenome=${array[1]}
	lastzParamFile=${array[2]}
	chainParamFile=${array[3]}
	
	#sort the chromosome file
	sort -k1,1 /cluster_path/ape_project/deletions_project/reference_genome_files/${species}/${refGenome}.chrom.sizes > /cluster_path/ape_project/deletions_project/reference_genome_files/${species}/${refGenome}.chrom.sizes.sorted
	mv /cluster_path/ape_project/deletions_project/reference_genome_files/${species}/${refGenome}.chrom.sizes.sorted /cluster_path/ape_project/deletions_project/reference_genome_files/${species}/${refGenome}.chrom.sizes
	
	QDIR=/cluster_path/ape_project/deletions_project/reference_genome_files/${species}
	WRKDIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}/work
	mkdir -p ${WRKDIR}
	
	#obtain the respective lastz/chain parameter files
	cp ${lastzParamFile} ${WRKDIR}/lastzParamFile.txt
	cp ${chainParamFile} ${WRKDIR}/chainParamFile.txt
	
	qsub RunLastzChain_revised.sh ${TARGETASSEMBNAME} ${refGenome} ${TDIR} ${QDIR} ${WRKDIR} ${WRKDIR}/lastzParamFile.txt ${WRKDIR}/chainParamFile.txt
	
done< /cluster_path/ape_project/deletions_project/msa/chimp_10_way/genomes_of_interest.txt

#check to see if we have all the right number of jobs for chainJobList
#refGenome=rheMac8
#wc -l /cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}/work/target.list
#wc -l /cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}/work/chainJobsList 
#${WRKDIR}/${QNAME}PartList
#second loop runs lastz
refGenome=rheMac8
species=macaque
QDIR=/cluster_path/ape_project/deletions_project/reference_genome_files/${species}
WRKDIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}/work
mkdir -p ${WRKDIR}

chain_param_dir=/cluster_path/ape_project/deletions_project/msa/chainParamFiles 
lastz_param_dir=/cluster_path/ape_project/deletions_project/msa/lastzParamFiles

cp ${lastz_param_dir}/primateLastzParamFile.txt ${WRKDIR}/lastzParamFile.txt
cp ${chain_param_dir}/primateChainParamFile.txt ${WRKDIR}/chainParamFile.txt
qsub RunLastzChain_revised.sh ${TARGETASSEMBNAME} ${refGenome} ${TDIR} ${QDIR} ${WRKDIR} ${WRKDIR}/lastzParamFile.txt ${WRKDIR}/chainParamFile.txt

while read -a array
do
	refGenome=bosTau8
	WRKDIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}/work
	cd ${WRKDIR}
	qsub ${WRKDIR}/runAllLastz.sh ${WRKDIR}/jobList
done

#now run post_lastz_check.sh on each genome to ensure that there were no errors
#then, run post_lastz_check_exists.sh to make sure that none of the runs have error
while read -a array
do
	refGenome=${array[1]} 
	#qsub post_pairalign_check_exists.sh ${refGenome} 

	#refGenome=rheMac8
	WRKDIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}/work
	#outFile=${WRKDIR}/check_file.txt
	#wc -l ${outFile}
	
	#qsub post_pairalign_check_runs.sh ${refGenome} 
	#head ${WRKDIR}/run_check.txt
done< /cluster_path/ape_project/deletions_project/msa/chimp_10_way/genomes_of_interest.txt

#for the runs that didn't finish, rerun them again, with a bit more memory for them to finish
while read -a array
do
	refGenome=${array[1]}
	WRKDIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}/work
	cd ${WRKDIR}
	if [[ -e ${WRKDIR}/check_file.txt ]];
	then


		cat ${WRKDIR}/check_file.txt | awk -F"[/.]" '{print $(NF-5)"."$(NF-4)"-"$(NF-3)"."$(NF-2)}' > ${WRKDIR}/idlist.txt
		awk '{print $4"-"$5"\t"$0}' ${WRKDIR}/jobList > ${WRKDIR}/jobList.id
		join -11 -21 -t $'\t' <(cat ${WRKDIR}/jobList.id | sort) <(cat ${WRKDIR}/idlist.txt | sort) -o '1.2 1.3 1.4 1.5' > ${WRKDIR}/jobList.v2
		
		del_dir=/cluster_path/ape_project/deletions_project
		totalMem=2G
		totalNumTasks=$(wc -l ${WRKDIR}/jobList.v2 | awk '{print $1}')	
		sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l m_mem_free=totalMem/#$ -l m_mem_free=${totalMem}/g" ${del_dir}/runTasks.sh > ${WRKDIR}/runAllLastz.2.sh

		qsub ${WRKDIR}/runAllLastz.2.sh ${WRKDIR}/jobList.v2
	fi
	
done< /cluster_path/ape_project/deletions_project/msa/chimp_10_way/genomes_of_interest.txt


refGenome=rheMac8
WRKDIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}/work
outFile=${WRKDIR}/check_file.txt
wc -l ${outFile}
/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_bosTau8/work/psl/panTro4.2bit:chr2B:0-20010000/panTro4.2bit:chr2B:0-20010000.bosTau8.2bit:chrX:20000000-40000000.psl.gz
/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_bosTau8/work/psl/panTro4.2bit\:chr2B\:0-20010000/                                  
/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_bosTau8/work/psl/panTro4.2bit\:chr2B\:0-20010000/panTro4.2bit\:chr2B\:0-20010000.bosTau8.2bit\:chrX\:20000000-40000000.psl.gz
grep chr1:0-20000000 ${WRKDIR}/jobList | grep chr2B:0-20010000 | head 
speciescheck=macaque
cd /cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}/work


#now after lastz finishes running, run a task array to generate chains for all target sequences 
while read -a array
do
	species=${array[0]}
	refGenome=${array[1]}
	
	TDIR=/cluster_path/ape_project/deletions_project/reference_genome_files/${TARGET}
	QDIR=/cluster_path/ape_project/deletions_project/reference_genome_files/${species}
	TARGETF=${TDIR}/${TARGETASSEMBNAME}.2bit
	QUERYF=${QDIR}/${refGenome}.2bit
	
	WRKDIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}/work
	#sed 's/$(awk -v runNum=${runNum} \'NR==runNum{print}\' ${tasklist})/eval $(awk -v runNum=${runNum} \'NR==runNum{print}\' ${tasklist})/g' ${WRKDIR}/runAllChains.sh

	del_dir=/cluster_path/ape_project/deletions_project
	totalMem=2G
	totalNumTasks=$(wc -l ${WRKDIR}/chainJobsList | awk '{print $1}')	
	#submit to long because you need all jobs to finish, can't tell if job dies in short queue due to time limit or it just finishes...
	sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l m_mem_free=totalMem/#$ -l m_mem_free=${totalMem}/g" ${del_dir}/runTasksEval.sh > ${WRKDIR}/runAllChains.sh
	
	rm -f ${WRKDIR}/chainJobsList
	rm ${WRKDIR}/chain/*
	for partFile in $(ls ${WRKDIR}/${TNAME}PartList)
	do

		echo "gunzip -c ${WRKDIR}/psl/${partFile}/${partFile}.*.psl.gz | axtChain -psl -verbose=0 ${chainParams} stdin ${TARGETF} ${QUERYF} stdout | chainAntiRepeat ${TARGETF} ${QUERYF} stdin ${WRKDIR}/chain/${partFile}.chain && touch ${WRKDIR}/chain/${partFile}.chain.finished" >> ${WRKDIR}/chainJobsList
	done 
	cd ${WRKDIR}
	
	qsub ${WRKDIR}/runAllChains.sh ${WRKDIR}/chainJobsList
	
done< /cluster_path/ape_project/deletions_project/msa/chimp_10_way/genomes_of_interest.txt


####check again for error messages, check_chains_test.sh#####
#gap locations aren't aligned
while read -a array
do
	refGenome=${array[1]}
	WRKDIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}/work
	#cat $id > ${del_dir}/check_log.txt 
	#check to see if chain file exists
	#qsub check_chains.sh ${refGenome} ${TARGETASSEMBNAME}
	wc -l ${WRKDIR}/chain_check.txt
	
	#look at output chain run
	#rm -f ${WRKDIR}/all_chains_log.txt
	for fileID in runAllChains.sh.*
	do
		#echo $fileID >>${WRKDIR}/all_chains_log.txt
		#cat $fileID >>${WRKDIR}/all_chains_log.txt
		paste <(echo $fileID) <(awk '{print $1}' $fileID) >> ${WRKDIR}/all_chains_log.txt
	done
done< /cluster_path/ape_project/deletions_project/msa/chimp_10_way/genomes_of_interest.txt

while read -a array
do
	refGenome=${array[1]}
	WRKDIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}/work
	#cat $id > ${del_dir}/check_log.txt 
	#check to see if chain file exists
	#qsub check_chains.sh ${refGenome} ${TARGETASSEMBNAME}
	wc -l ${WRKDIR}/chain_check.txt
	#echo $refGenome
	#awk '$2!=""{print}' ${WRKDIR}/all_chains_log.txt
done< /cluster_path/ape_project/deletions_project/msa/chimp_10_way/genomes_of_interest.txt

refGenome=monDom5
WRKDIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}/work
head ${WRKDIR}/chain_check.txt
grep panTro4.2bit:chr19:40000000-60010000.chain ${WRKDIR}/chainJobsList 
awk 'NR==96{print}' ${WRKDIR}/chainJobsList 

while read -a array
do
	refGenome=${array[1]}
	WRKDIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}/work
	cat ${WRKDIR}/runAllChains.sh.* >> ${WRKDIR}/all_chains_log.txt
	
	sort ${WRKDIR}/all_chains_log.txt | uniq | head
done< /cluster_path/ape_project/deletions_project/msa/chimp_10_way/genomes_of_interest.txt

rm -f all_chains_log.txt
for fileID in runAllChains.sh.*
do
	echo $fileID >>all_chains_log.txt
	cat $fileID >>all_chains_log.txt
done
cat runAllChains.sh.* >> all_chains_log.txt

awk '$1=="stdin"{print}' all_chains_log.txt
sort all_chains_log.txt | uniq | head

#######

#now combine all chain files into a larger file
while read -a array
do
	species=${array[0]}
	refGenome=${array[1]}
	WRKDIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}/work
	OUTDIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}
	qsub chainMerge.sh ${TARGETASSEMBNAME} ${refGenome} ${WRKDIR} ${OUTDIR}
	
done< /cluster_path/ape_project/deletions_project/msa/chimp_10_way/genomes_of_interest.txt


refGenome=ponAbe2
WRKDIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}/work
OUTDIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}
qsub chainMerge.sh ${TARGETASSEMBNAME} ${refGenome} ${WRKDIR} ${OUTDIR}


#after chaining, convert to Net format, runAllNets.sh generates synNet files also for MAF later
while read -a array
do
	species=${array[0]}
	refGenome=${array[1]}
	QDIR=/cluster_path/ape_project/deletions_project/reference_genome_files/${species}
	PAIR_DIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}
	qsub runAllNets.sh ${TARGETASSEMBNAME} ${refGenome} ${TDIR} ${QDIR} ${PAIR_DIR}
done< /cluster_path/ape_project/deletions_project/msa/chimp_10_way/genomes_of_interest.txt

cat runAllNets.sh.o* > look.txt

species=orangutan
refGenome=ponAbe2
QDIR=/cluster_path/ape_project/deletions_project/reference_genome_files/${species}
PAIR_DIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}
#PAIR_DIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/panTro4_hg38
qsub runAllNets.sh ${TARGETASSEMBNAME} ${refGenome} ${TDIR} ${QDIR} ${PAIR_DIR}
#after netting, prepare for MAF run

#now run MAF
#run the commands in runMultiz.sh
#http://redmine.soe.ucsc.edu/forum/index.php?t=msg&goto=12389&S=aa617af7c6ed6eaa5c46f4e8da33cf8f


#this is to delete stuff
while read -a array
do
	#species=${array[0]}
	refGenome=${array[1]}
	WRKDIR=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/query_${refGenome}/work
	rm ${WRKDIR}/runAllLastz.sh.*
done< /cluster_path/ape_project/deletions_project/msa/chimp_10_way/genomes_of_interest.txt

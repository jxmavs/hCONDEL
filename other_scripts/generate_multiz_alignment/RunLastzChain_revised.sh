#!/bin/sh
#$ -cwd 
#$ -l m_mem_free=2G
#$ -q short
#$ -V
#$ -j y
source /broad/software/scripts/useuse
reuse -q GCC-4.9
reuse -q Python-2.7
LD_LIBRARY_PATH=/broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.3.0/lib64/R/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/curl_7.47.1/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/pcre_8.38/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/xz_5.2.2/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/bzip2_1.0.6/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_4.9.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/perl_5.8.9/lib:/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/oracle_instantclient-10.2.0.4.0/instantclient_10_2:/broad/uge/8.4.0/lib/lx-amd64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/zlib_1.2.6/lib64:/broad/software/free/Linux/redhat_6_x86_64/pkgs/hdf5_1.8.9/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/graphviz_2.28.0/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/db_4.7.25/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/tcltk8.5.9/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/jre/lib/amd64/server:/broad/software/free/Linux/redhat_6_x86_64/pkgs/jdk1.8.0_92/lib:/cluster_path/bin/jellyfish-1.1.11/lib:/cluster_path/bin/openmpi/lib:/cluster_path/bin/libncurses5/include:/cluster_path/bin/packages/lib:/broad/software/free/Linux/redhat_6_x86_64/pkgs/libiconv_1.13.1/lib

# This script generates the psl formats from running lastz, then performs chaining
# this script is analagous to the original RunLastzChain.sh file
# except this is modified to on the Broad cluster

#  WRKDIR is where you want this to work
# TNAME is the ref. directory of the target sequence
# QNAME is the ref. directory of the query sequence

TNAME=$1
QNAME=$2
TDIR=$3
QDIR=$4
WRKDIR=$5
lastzParamFile=$6
chainParamFile=$7

#TNAME=panTro4
#QNAME=hg38
#TDIR=/cluster_path/ape_project/deletions_project/reference_genome_files/chimp
#QDIR=/cluster_path/ape_project/deletions_project/reference_genome_files/human
#WRKDIR=/cluster_path/ape_project/deletions_project/msa/work
#lastzParamFile=/cluster_path/ape_project/deletions_project/msa/lastzParamFile.txt
#chainParamFile=/cluster_path/ape_project/deletions_project/msa/chainParamFile.txt


chainParams=""
while read -a array
do
	key=${array[0]}
	value=${array[1]}
	chainParams=$(echo "${chainParams} ${key}=${value}")

done < ${chainParamFile}

lastzParams=""
while read -a array
do
	key=${array[0]}
	value=${array[1]}
	lastzParams=$(echo "${lastzParams} ${key}=${value}")

done < ${lastzParamFile}



TARGET=${TDIR}/${TNAME}.2bit
QUERY=${QDIR}/${QNAME}.2bit

mkdir -p ${WRKDIR}/${QNAME}PartList
mkdir -p ${WRKDIR}/${TNAME}PartList

#divide up sequences into smaller sequences of 10000000 bp, note that the 5th parameter is the limit on number of sequences lumped together, set it to 0

partitionSequence.pl 20000000 10000 ${TARGET} ${TDIR}/${TNAME}.chrom.sizes -xdir ${WRKDIR}/xdir.sh -rawDir ${WRKDIR}/psl 2000 -lstDir ${WRKDIR}/${TNAME}PartList > ${WRKDIR}/${TNAME}.part.list


#create directory for each partition file, this will store the pairwise outputs
bash ${WRKDIR}/xdir.sh


#partition the query sequences
partitionSequence.pl 20000000 0 ${QUERY} ${QDIR}/${QNAME}.chrom.sizes 5000 -lstDir ${WRKDIR}/${QNAME}PartList > ${WRKDIR}/${QNAME}.part.list

#see number of lines in each file
#how long does it take for a target run of 5000 files?
for f in ${WRKDIR}/${QNAME}PartList/*.lst
do
	wc -l ${f}
done

for f in ${WRKDIR}/${TNAME}PartList/*.lst
do
	wc -l ${f}	
done

#grep part015.lst jobList


#put all files into a single list for both query and target 
#################################################
#this part extracts out the chromosome coordinates from each target file for reading later

rm -f ${WRKDIR}/target.list

grep -v PartList ${WRKDIR}/${TNAME}.part.list > ${WRKDIR}/target_single_file.list

while read -a array
do

	read_path="${array[0]}"

	fa_file_end=$(echo ${read_path} | awk -F'/' '{print $NF}')
	#get chr position from file end from read path
	position=$(echo ${fa_file_end} | awk -F":" '{print $(NF-1)":"$NF}')
	
	echo "${position}" > "${WRKDIR}/${TNAME}PartList/${fa_file_end}"
	echo "${WRKDIR}/${TNAME}PartList/${fa_file_end}" >> ${WRKDIR}/target.list
	
done <${WRKDIR}/target_single_file.list

rm ${WRKDIR}/target_single_file.list

for F in ${WRKDIR}/${TNAME}PartList/*.lst
do

    
    awk -F'[/:]' '{print $(NF-1)":"$NF}' ${F} > ${F}.pos
    mv ${F}.pos ${F}
    echo "${F}" >> ${WRKDIR}/target.list
done 



#################################################
#

rm -f ${WRKDIR}/query.list

grep -v PartList ${WRKDIR}/${QNAME}.part.list > ${WRKDIR}/query_single_file.list

while read -a array
do

	read_path="${array[0]}"
	#echo ${read_path}
	fa_file_end=$(echo ${read_path} | awk -F'/' '{print $NF}')
	position=$(echo ${fa_file_end} | awk -F":" '{print $(NF-1)":"$NF}')
	
	echo "${position}" > "${WRKDIR}/${QNAME}PartList/${fa_file_end}"
	echo "${WRKDIR}/${QNAME}PartList/${fa_file_end}" >> ${WRKDIR}/query.list
done <${WRKDIR}/query_single_file.list

rm ${WRKDIR}/query_single_file.list

for F in ${WRKDIR}/${QNAME}PartList/*.lst
do

    awk -F'[/:]' '{print $(NF-1)":"$NF}' ${F} > ${F}.pos
    mv ${F}.pos ${F}
    echo "${F}" >> ${WRKDIR}/query.list
done 

##################

echo "#LOOP" > ${WRKDIR}/template
echo "${WRKDIR}/runLastz \$(path1) \$(path2) \$(file1) \$(file2)" >> ${WRKDIR}/template
echo "#ENDLOOP" >> ${WRKDIR}/template

#T is file path
#FT is name at end of file path
#same concept with Q

#now you may have multiple query sequences in each file, but only one target sequence

cat <<_EOF_ > ${WRKDIR}/runLastz

#T contains the direct path to the target sequences
#Q contains the fasta coordinates, used for -seqList

T=\$1
Q=\$2
FT=\$3
FQ=\$4


tmpDir=/cluster_path_temp/msa/${TNAME}_\${FT}_${QNAME}_\${FQ}
mkdir -p \${tmpDir} ${WRKDIR}/psl ${WRKDIR}/psl/\${FT}
	
#only one target sequence run at a time
twoBitToFa -seqList=\${T} ${TDIR}/${TNAME}.2bit \${tmpDir}/${TNAME}_\${FT}.fa

#get a bunch of fasta files to one file
twoBitToFa -seqList=\${Q} ${QDIR}/${QNAME}.2bit \${tmpDir}/${QNAME}_\${FQ}.fa

#get chromosome headers from file
awk -F'[:-]' '{print \$2"\t"\$1"\t"\$3-\$2"\t"\$1}' \${T} | sort -k2,2 > \${tmpDir}/target_list.txt
awk -F'[:-]' '{print \$2"\t"\$1"\t"\$3-\$2"\t"\$1}' \${Q} | sort -k2,2 > \${tmpDir}/query_list.txt

#intersect chromosome headers with list of chromosomes from sizes
#make sure that the chromosome files are already sorted
join -12 -21 -t $'\t' <( cat \${tmpDir}/target_list.txt ) <( cat ${TDIR}/${TNAME}.chrom.sizes ) -o '1.1 1.2 1.3 1.4 2.2' > \${tmpDir}/target.lift
join -12 -21 -t $'\t' <( cat \${tmpDir}/query_list.txt ) <( cat ${QDIR}/${QNAME}.chrom.sizes ) -o '1.1 1.2 1.3 1.4 2.2' > \${tmpDir}/query.lift

#run lastz
lastz --format=axt \${tmpDir}/${TNAME}_\${FT}.fa[multiple]  \${tmpDir}/${QNAME}_\${FQ}.fa ${lastzParams} > \${tmpDir}/\${FT}.\${FQ}.pre.axt
	

#liftUp the target and query coordinates to the right ones, it changes
liftUp -type=.axt \${tmpDir}/\${FT}.\${FQ}.pre.targetlift.axt  \${tmpDir}/target.lift error \${tmpDir}/\${FT}.\${FQ}.pre.axt 
liftUp -nohead -axtQ -type=.axt \${tmpDir}/\${FT}.\${FQ}.axt  \${tmpDir}/query.lift error \${tmpDir}/\${FT}.\${FQ}.pre.targetlift.axt 

#convert axt to psl	
axtToPsl \${tmpDir}/\${FT}.\${FQ}.axt ${TDIR}/${TNAME}.chrom.sizes ${QDIR}/${QNAME}.chrom.sizes ${WRKDIR}/psl/\${FT}/\${FT}.\${FQ}.psl 
 
gzip -f ${WRKDIR}/psl/\${FT}/\${FT}.\${FQ}.psl  

rm -fr \${tmpDir} 
	

_EOF_

chmod +x ${WRKDIR}/runLastz
#gensub2 creates the jobList to run
echo "ready to run lastz cluster job:"
gensub2 ${WRKDIR}/target.list ${WRKDIR}/query.list ${WRKDIR}/template ${WRKDIR}/jobList

#time ./runLastz  /cluster_path/ape_project/deletions_project/msa/work/panTro4PartList/part007.lst /cluster_path/ape_project/deletions_project/msa/work/hg38PartList/part006.lst part007.lst part006.lst
#echo "para make jobList"
del_dir=/cluster_path/ape_project/deletions_project
totalMem=1G
totalNumTasks=$(wc -l ${WRKDIR}/jobList | awk '{print $1}')	
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l m_mem_free=totalMem/#$ -l m_mem_free=${totalMem}/g" ${del_dir}/runTasks.sh > ${WRKDIR}/runAllLastz.sh

#run the folllowing below later
#cd ${WRKDIR}
#qsub ${WRKDIR}/runAllLastz.sh ${WRKDIR}/jobList
echo "total number of jobs"
wc -l ${WRKDIR}/jobList

mkdir -p ${WRKDIR}/chain

#decompress all psl files and chain all the files
#each directories has multiple targets
rm -f ${WRKDIR}/chainJobsList

for partFile in $(ls ${WRKDIR}/${TNAME}PartList)
do
	rm -f ${WRKDIR}/chain/${partFile}.chain.finished
	echo "gunzip -c ${WRKDIR}/psl/${partFile}/${partFile}.*.psl.gz | axtChain -psl -verbose=0 ${chainParams} stdin ${TARGET} ${QUERY} stdout | chainAntiRepeat ${TARGET} ${QUERY} stdin ${WRKDIR}/chain/${partFile}.chain && touch ${WRKDIR}/chain/${partFile}.chain.finished" >> ${WRKDIR}/chainJobsList
done 

#run post_pairalign_check_exists.sh to ensure that lastz files exist
#then run reRunLastz.sh to finish up the remaining files

del_dir=/cluster_path/ape_project/deletions_project
totalMem=1G
totalNumTasks=$(wc -l ${WRKDIR}/chainJobsList | awk '{print $1}')	
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l m_mem_free=totalMem/#$ -l m_mem_free=${totalMem}/g" ${del_dir}/runTasksEval.sh > ${WRKDIR}/runAllChains.sh


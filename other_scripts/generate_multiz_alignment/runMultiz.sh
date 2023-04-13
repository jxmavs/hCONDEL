#https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/makeDb/doc/hg38/multiz100way.txt
#https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/makeDb/doc/galVar1/multiz5way.txt
#https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/makeDb/doc/panTro3.txt
#12-way multiz

msa_dir=/cluster_path/ape_project/deletions_project/msa/chimp_10_way
#multiz_dir=/cluster_path/bin/multiz-tba.012109
del_dir=/cluster_path/ape_project/deletions_project
targetRefGenome=panTro4


#mkdir ${msa_dir}/multiz12way
#cd ${msa_dir}/multiz12way

wget https://raw.githubusercontent.com/ucscGenomeBrowser/kent/master/src/hg/utils/phyloTrees/191way.nh -O ${msa_dir}/191way.nh

awk '{print $2}' /cluster_path/ape_project/deletions_project/msa/chimp_10_way/genomes_of_interest.txt > ${msa_dir}/other_genomes.txt

#edit it a bit

#${msa_dir}/other_genomes.txt - the names in the file corresponds to the tree file ${msa_dir}/191way.nh, use this
echo ${targetRefGenome} > ${msa_dir}/all_genomes.txt
awk '{print $1}' ${msa_dir}/other_genomes.txt >> ${msa_dir}/all_genomes.txt

#	All distances remain as specified in the 46way.nh
tree_doctor --prune-all-but $(awk -v nlines=$(wc -l ${msa_dir}/all_genomes.txt | awk '{print $1}')  'NR!=nlines{printf "%s,", $1 } END{printf "%s", $1}' ${msa_dir}/all_genomes.txt) ${msa_dir}/191way.nh | sed -e "s/rheMac3/rheMac8/; s/gorGor3/gorGor4/; s/ponAbe3/ponAbe2/; s/ornAna5/ornAna2/ " > ${msa_dir}/revisedtree.nh

#	what that looks like:
cat ${msa_dir}/revisedtree.nh


mv ${msa_dir}/revisedtree.nh ${msa_dir}/${targetRefGenome}.revisedtree.nh


all_dists ${msa_dir}/${targetRefGenome}.revisedtree.nh > ${msa_dir}/revisedtree.distances.txt
 #	Use this output to create the table below
grep -i ${targetRefGenome} ${msa_dir}/revisedtree.distances.txt | sort -k3,3n


# None of this concern for distances matters in building the first step, the
# maf files.
# create species list and stripped down tree for autoMZ
#strip down tree is just a tree with names
sed 's/[a-z][a-z]*_//g; s/:[0-9\.][0-9\.]*//g; s/;//; /^ *$/d' ${msa_dir}/${targetRefGenome}.revisedtree.nh > ${msa_dir}/tmp.nh

#gets rid of commas
cat ${msa_dir}/tmp.nh | sed 's/ //g; s/,/ /g' > ${msa_dir}/tree.nh

#species.list is just a list of species, with a space between them
sed 's/[()]//g; s/,/ /g' ${msa_dir}/tree.nh > ${msa_dir}/species.list
    
    # split the maf files into a set of hashed named files
    # this hash named split keeps the same chr/contig names in the same
    # named hash file. !!! IMPORTANT
########### stopped here!!!!!!! 8/4/16 ##########
rm -r ${msa_dir}/mafSplit
mkdir -p ${msa_dir}/mafSplit


#for species in $(sed -e "s/${targetRefGenome} //" ${msa_dir}/species.list)
#do
rm -f ${msa_dir}/mafSplitJobList

while read -a array
do
	species=${array[0]}
	refGenome=${array[1]}
    echo "${species}"
    rm -fr ${msa_dir}/mafSplit/${refGenome}
    mkdir -p ${msa_dir}/mafSplit/${refGenome}
    #cd $D
	echo "mafSplit -byTarget -useHashedName=10 /dev/null ${msa_dir}/mafSplit/${refGenome}/ ${msa_dir}/query_${refGenome}/${targetRefGenome}.${refGenome}.synNet.maf.gz" >> ${msa_dir}/mafSplitJobList
    #cd ..
    #example of axt conversion to MAF
done< /cluster_path/ape_project/deletions_project/msa/chimp_10_way/genomes_of_interest.txt

del_dir=/cluster_path/ape_project/deletions_project
totalMem=5G
totalNumTasks=$(wc -l ${msa_dir}/mafSplitJobList | awk '{print $1}')	
#submit to long because you need all jobs to finish, can't tell if job dies in short queue due to time limit or it just finishes...
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l m_mem_free=totalMem/#$ -l m_mem_free=${totalMem}/g" ${del_dir}/runTasks.sh > ${msa_dir}/runAllMAFSplit.sh

qsub ${msa_dir}/runAllMAFSplit.sh ${msa_dir}/mafSplitJobList
cat ${msa_dir}/runAllMAFSplit.sh.*> ${msa_dir}/look.txt
grep -v Splitting ${msa_dir}/look.txt | head 


# construct a list of all possible maf file names.
# they do not all exist in each of the species directories
find ${msa_dir}/mafSplit -type f | wc -l
# 2792
find ${msa_dir}/mafSplit -type f | grep ".maf$" | xargs -L 1 basename | sort -u > ${msa_dir}/mafSplit/maf.list

wc -l ${msa_dir}/mafSplit/maf.list
# 484 maf.list
rm -fr ${msa_dir}/splitRun
mkdir -p ${msa_dir}/splitRun

mkdir -p ${msa_dir}/splitRun/maf ${msa_dir}/splitRun/run
#cd run
#bash ${del_dir}/createMz.sh
cat > ${msa_dir}/splitRun/run/autoMultiz.sh <<EOF
#!/bin/sh 
db=${targetRefGenome}
c=\$1
result=\$2

tmp=/cluster_path_temp/\$db/multiz.\$c
pairs=${msa_dir}/mafSplit
rm -fr \$tmp
mkdir -p \$tmp
cp ${msa_dir}/tree.nh ${msa_dir}/species.list \$tmp
pushd \$tmp > /dev/null

for s in \$(sed -e "s/\$db //" ${msa_dir}/species.list); 
do
	
    in=\$pairs/\$s/\$c.maf
    out=\$db.\$s.sing.maf
    #echo \$in
    #decompress maf file if necessary
    if [[ -e \$in.gz ]]; then
        gunzip -c \$in.gz > \$out
		if [[ ! -s \$out ]]; then
	    	echo "##maf version=1 scoring=autoMZ" > \$out
		fi

	#if already decmompressed, then create a symbolic link to real file
    elif [[ -e \$in ]]; then
        ln -s \$in \$out
    else
        echo "##maf version=1 scoring=autoMZ" > \$out
    fi
done

tree=\$(cat ${msa_dir}/tree.nh)

roast + T=\$tmp E=\$db "\$tree" \$db.*.sing.maf \$c.maf  > /dev/null

popd > /dev/null

cp -p \$tmp/\$c.maf \$result
rm -fr \$tmp

EOF

# << happy emacs
chmod +x ${msa_dir}/splitRun/run/autoMultiz.sh

cat  <<EOF > ${msa_dir}/splitRun/run/template
#LOOP
${msa_dir}/splitRun/run/autoMultiz.sh \$(root1) ${msa_dir}/splitRun/maf/\$(root1).maf
#ENDLOOP
EOF


#rm -fr /cluster_path_temp/\$db

cd ${msa_dir}/splitRun/run
#gensub 2 creates joblist, each line has a job to ru 
gensub2 ${msa_dir}/mafSplit/maf.list single ${msa_dir}/splitRun/run/template ${msa_dir}/splitRun/run/jobList
#tac jobList > smallJobsFirst

#need more time, try submitting to long queue?
totalMem=5G
totalNumTasks=$(wc -l ${msa_dir}/splitRun/run/jobList | awk '{print $1}')	
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l m_mem_free=totalMem/#$ -l m_mem_free=${totalMem}/g; s/#$ -q short/#$ -q long/g" ${del_dir}/runTasks.sh > ${msa_dir}/splitRun/run/runAllMAF.sh
qsub ${msa_dir}/splitRun/run/runAllMAF.sh ${msa_dir}/splitRun/run/jobList

#checked number of files is 493
ls -l ${msa_dir}/splitRun/maf | wc -l
cat ${msa_dir}/splitRun/run/runAllMAF.sh.* > ${msa_dir}/splitRun/run/look.txt
head ${msa_dir}/splitRun/run/look.txt

    


# assemble into a single maf file, this just removes the # tags
    
head -1 ${msa_dir}/splitRun/maf/000.maf >${msa_dir}/target_${targetRefGenome}_multiz.maf
for F in ${msa_dir}/splitRun/maf/*.maf
do
    egrep -v "^#" ${F}  >> ${msa_dir}/target_${targetRefGenome}_multiz.maf
done

tail -1 ${msa_dir}/splitRun/maf/000.maf >> ${msa_dir}/target_${targetRefGenome}_multiz.maf


##############phastcons################

#get the ensembl gene list 
#ensgene can also be downloaded here 
#wget http://hgdownload.cse.ucsc.edu/goldenpath/panTro4/database/ensGene.txt.gz

#msa_dir=/cluster_path/ape_project/deletions_project/msa/multiz12way 
msa_dir=/cluster_path/ape_project/deletions_project/msa/chimp_10_way

dir4d=${msa_dir}/4d
tempdir4d=/cluster_path_temp/msa
rm -fr ${tempdir4d}
mkdir -p ${tempdir4d}
rm -fr ${dir4d}
mkdir -p ${dir4d}
mkdir -p ${dir4d}/run
mkdir -p ${dir4d}/mfa
mkdir -p ${dir4d}/mafSplit

#parse the ensGene file first
#hgsql -N -e "select * from refGene" panTro3 | cut -f 2-20 | egrep -E -v "chrM|random|chrUn" > ${dir4d}/ensGene.gp
targetRefGenome=panTro4
hgsql ${targetRefGenome} -Ne "select * from ensGene WHERE cdsEnd > cdsStart;" | cut -f 2-20 | egrep -E -v "chrM|random|chrUn" > ${dir4d}/ensGene.gp
#	verify chromosome selection, should just be the ordinary chroms:
cut -f2 ${dir4d}/ensGene.gp | sort | uniq -c
wc -l *.gp

#genePredSingleCover removes gene isoforms - is that what it means by getting a single gene coverage?
#https://groups.google.com/a/soe.ucsc.edu/forum/#!msg/genome/WF-Qt1d9HV0/nXrgRmjlRzAJ
genePredSingleCover ${dir4d}/ensGene.gp stdout | sort > ${dir4d}/ensGeneNR.gp



#cd /hive/data/genomes/panTro3/bed/multiz12way/4d
#mkdir mafSplit
#cd mafSplit
#time mafSplit -byTarget -useFullSequenceName /dev/null . ../../anno/target_${targetRefGenome}_multiz.maf

time mafSplit -byTarget -useFullSequenceName /dev/null ${dir4d}/mafSplit/ ${msa_dir}/target_${targetRefGenome}_multiz.maf


cat <<_EOF_ > ${dir4d}/run/4d.sh
#!/bin/bash 
c=\$1
infile=${dir4d}/mafSplit/\$2
outfile=\$3
#cd /scratch/tmp
# 'clean' maf, perl -wpe is pretty much same as sed
#perl -wpe "s/^s ([^.]+)\.\S+/s \$c/" \$infile > ${tempdir4d}/\$c.maf  
cat \$infile > ${tempdir4d}/\$c.maf  

#gets only genes for specific chromosome, and changes the header name up a bit 
awk -v C=\$c '\$2 == C {print}' ${dir4d}/ensGeneNR.gp | sed -e "s/\t\$c\t/\t${targetRefGenome}.\$c\t/" > ${tempdir4d}/\$c.gp

#get number of lines
NL=\$(wc -l ${tempdir4d}/\$c.gp| awk '{print \$1}')

if [[ "\$NL" != "0" ]]; then
   echo ${tempdir4d}/\$c.maf
   echo ${tempdir4d}/\$c.gp 
   #--features reads file as a gff file, -i input is maf, --4d is extract 4-fold degenerate sites
   #ss file format is each line is a "column" of alignment containing 3 columns, which represent a codon
   msa_view --4d --features ${tempdir4d}/\$c.gp -i MAF ${tempdir4d}/\$c.maf -o SS > ${tempdir4d}/\$c.ss
   
   #this just obtains the 4d sites
   msa_view -i SS --tuple-size 1 ${tempdir4d}/\$c.ss > \$outfile
else
    echo "" > ${dir4d}/run/\$outfile
fi
#rm -f ${tempdir4d}/\$c.gp ${tempdir4d}/\$c.maf ${tempdir4d}/\$c.ss
_EOF_

##### testing  ######
c=chr22
perl -wpe "s/^s ([^.]+)\.\S+/s $c/" /cluster_path/ape_project/deletions_project/msa/multiz12way/4d/mafSplit/chr22.maf > testloo.txt 
# << happy emacs
chmod +x ${dir4d}/run/4d.sh
ls -1S ${dir4d}/mafSplit/*.maf | \
egrep -E -v "chrM|random|chrUn" | sed -e "s#.*multiz12way/4d/mafSplit/##" > ${dir4d}/run/maf.list
head ${dir4d}/mafSplit/chrY.maf
head /cluster_path_temp/msa/chrY.gp
head /cluster_path_temp/msa/chrY.ss
#################


chmod +x ${dir4d}/run/4d.sh

ls -1S ${dir4d}/mafSplit/*.maf | egrep -E -v "chrM|random|chrUn" | sed -e "s#${dir4d}/mafSplit/##" > ${dir4d}/run/maf.list

cat <<_EOF_ > ${dir4d}/run/template
#LOOP
${dir4d}/run/4d.sh \$(root1) \$(path1) ${dir4d}/mfa/\$(root1).mfa
#ENDLOOP
_EOF_

gensub2 ${dir4d}/run/maf.list single ${dir4d}/run/template stdout | tac > ${dir4d}/run/jobList


del_dir=/cluster_path/ape_project/deletions_project
totalMem=7G
totalNumTasks=$(wc -l ${dir4d}/run/jobList | awk '{print $1}')	
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l m_mem_free=totalMem/#$ -l m_mem_free=${totalMem}/g; s/#$ -q short/#$ -q long/g" ${del_dir}/runTasks.sh > ${dir4d}/run/runAll4d.sh
cd ${dir4d}/run
qsub ${dir4d}/run/runAll4d.sh ${dir4d}/run/jobList

cat ${dir4d}/run/runAll4d.sh.* > ${dir4d}/run/check_runs.txt 
sort ${dir4d}/run/check_runs.txt | uniq | head


#want comma-less species.list
msa_view --aggregate "$(cat ${msa_dir}/species.list)" ${dir4d}/mfa/*.mfa | sed s/"> "/">"/  > ${dir4d}/4d.all.mfa

#run phyloFit to obtain tree
#${msa_dir}/tmp.nh is just the tree without commas
phyloFit --EM --precision MED --msa-format FASTA --subst-mod REV -o ${dir4d}/all --tree ${msa_dir}/tmp.nh ${dir4d}/4d.all.mfa

#mv ${dir4d}/phyloFit.mod ${dir4d}/all.mod


####now make ss file ######

consSSDir=${msa_dir}/cons/SS
rm -fr ${consSSDir}
mkdir -p ${consSSDir}
#cd ${consSSDir}

#this script creates sufficient statistics from alignment, by splitting alignment
#each maf file will have its own set of windows

    cat <<_EOF_ > ${consSSDir}/mkSS.sh
#!/bin/bash 
c=\$1
MAF=${dir4d}/mafSplit/\$c.maf
WINDOWS=${consSSDir}/\$c
WC=\$(cat \$MAF | wc -l)
#NL is just the number of meta info headers, if a file is empty then it will just consist of only these lines
NL=\$(grep "^#" \$MAF | wc -l)



rm -fr \$WINDOWS
mkdir -p \$WINDOWS

#if file is not empty then split, NL should always be 1
#partitions maf into windows of 10000000, converts it to ss format
if [[ \$WC != \$NL ]]; then
msa_split \$MAF -i MAF -o SS -r \$WINDOWS/\$c -w 10000000,0 -I 1000 -B 5000
fi

_EOF_


chmod +x ${consSSDir}/mkSS.sh

cat <<_EOF_ > ${consSSDir}/template
#LOOP
${consSSDir}/mkSS.sh \$(root1)
#ENDLOOP
_EOF_


#	do the easy ones first to see some immediate results
ls -1S -r ${dir4d}/mafSplit | sed -e "s/.maf//" > ${consSSDir}/maf.list

gensub2 ${consSSDir}/maf.list single ${consSSDir}/template ${consSSDir}/jobList

cd ${consSSDir}
del_dir=/cluster_path/ape_project/deletions_project
totalMem=2G
totalNumTasks=$(wc -l ${consSSDir}/jobList | awk '{print $1}')	
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l m_mem_free=totalMem/#$ -l m_mem_free=${totalMem}/g; s/#$ -q short/#$ -q long/g" ${del_dir}/runTasks.sh > ${consSSDir}/mkAllSS.sh
qsub ${consSSDir}/mkAllSS.sh ${consSSDir}/jobList

#look for error messages
cat ${consSSDir}/mkAllSS.sh.* > ${consSSDir}/check_runs.txt 

#######checking the warning#######

grep -v Writing ${consSSDir}/check_runs.txt | grep -v Reading | grep -v Creating | grep -v Done | uniq

rm -f ${consSSDir}/all_warnings.txt
for fileIter in ${consSSDir}/mkAllSS.sh.*
do
	paste <(echo $fileIter) <(grep WARNING $fileIter) >> ${consSSDir}/all_warnings.txt
done

grep WARNING ${consSSDir}/all_warnings.txt

awk 'NR==3732{print}' ${consSSDir}/jobList
awk 'NR==3755{print}' ${consSSDir}/jobList
awk '($3>50000000 && $3<60000000 && $2~"panTro4"){print }'  ${dir4d}/mafSplit/chr9.maf | head 

####### run phastcons next########


#phastconsdir=/cluster_path/ape_project/deletions_project/msa/multiz12way/cons/run.phast
phastconsdir=/cluster_path/ape_project/deletions_project/msa/chimp_10_way/cons/run.phast
rm -fr ${phastconsdir}
mkdir -p ${phastconsdir}/pp

mkdir -p ${phastconsdir}/bed

    #	there are going to be only one phastCons run using
    #	this same script.  It triggers off of the current working directory
    #	$cwd:t which is the "grp" in this script.  It is:
    #	all 

    cat <<_EOF_ > ${phastconsdir}/doPhast.sh
#!/bin/bash 

c=\$1
f=\$2

len=\$3
cov=\$4
rho=\$5
tree_path=\$6

consSS=${consSSDir}

tmp=/cluster_path_temp/msa/phastrun/$f
mkdir -p \$tmp

#change to just be all.mod tree from useGrp
#in ss format, tuple IDX order lists order of alignment, so can use sufficient file format
cd \${tmp}
phastCons \${consSS}/\$c/\$f.ss \${tree_path} \
    --rho \$rho --expected-length \$len --target-coverage \$cov --quiet \
    --seqname \$c --idpref \$c --most-conserved ${phastconsdir}/bed/\$f.bed --score > ${phastconsdir}/pp/\$f.pp
    

sleep 4

_EOF_

chmod +x ${phastconsdir}/doPhast.sh

tree_path=${dir4d}/all.mod

cat <<_EOF_ > ${phastconsdir}/template
#LOOP
${phastconsdir}/doPhast.sh \$(root1) \$(file1) 45 0.3 0.3 ${tree_path} 
#ENDLOOP
_EOF_

#### stopped here #####

find ${consSSDir} -type f | grep ".ss$" | sed -e 's/.ss$//' > ${phastconsdir}/ss.list
wc -l ${phastconsdir}/ss.list

gensub2 ${phastconsdir}/ss.list single ${phastconsdir}/template ${phastconsdir}/jobList

cd ${phastconsdir}
del_dir=/cluster_path/ape_project/deletions_project
totalMem=5G
totalNumTasks=$(wc -l ${phastconsdir}/jobList | awk '{print $1}')	
sed "s/#$ -t 1-totalNumTasks/#$ -t 1-${totalNumTasks}/g; s/#$ -l m_mem_free=totalMem/#$ -l m_mem_free=${totalMem}/g; s/#$ -q short/#$ -q long/g" ${del_dir}/runTasks.sh > ${phastconsdir}/doPhastAll.sh
qsub ${phastconsdir}/doPhastAll.sh ${phastconsdir}/jobList

#check no errors
cat ${phastconsdir}/doPhastAll.sh.* > ${phastconsdir}/check_runs.txt 
sort ${phastconsdir}/check_runs.txt | uniq 
ls -l  ${phastconsdir}/bed/* | wc -l
wc -l ${phastconsdir}/jobList

#find ${phastconsdir}/bed -type f | grep ".bed$" | xargs cat | sort -k1,1 -k2,2n | awk '{printf "%s\t%d\t%d\tlod=%d\t%s\n", $1, $2, $3, $5, $5;}' > ${phastconsdir}/tmpMostConserved.bed

#lodToBedScore ${phastconsdir}/tmpMostConserved.bed > ${phastconsdir}/mostConserved.bed

find ${phastconsdir}/bed -type f | grep ".bed$" | xargs cat | sort -k1,1 -k2,2n > ${phastconsdir}/mostConserved.bed

#number of conserved elements
wc -l ${phastconsdir}/mostConserved.bed

head ${phastconsdir}/mostConserved.bed

#maximum length is 999 bp
#awk '{print $3-$2}' ${phastconsdir}/mostConserved.bed | sort -k1,1n | tail
#awk '{print $3-$2}' ${phastconsdir}/mostConserved.bed | awk '$1<20{print}' | wc -l
#awk '{print $3-$2}' ${phastconsdir}/mostConserved.bed | awk '{sum+=$1} END {print sum}'
#125070838
#125668334
	
	



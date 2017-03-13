# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 11:23:11 2016

@author: lenovo
"""

# version v 0.1 2016-09-01
# revieve torrent raw data
# working flow contains, samtools index, tmap, samtools sort, samtools depth, retmapping, BP filtration
# by default we spare one cpu for each torrent sequencing data
# use a bi-fork method for boost depth-add step
# use ordered dict to avoid for loop
# considerring soft-clipping head, duplication and deletion stt and stp filtration standard should be modified
# add soft-clipping head or tail tag for R script analysis

################
# changing log #
# if softclipping < 20:
#	filtra
# if rawreads mapping FLAG != raw-reads cut mapping FLAG:
		# filtration
################
# apply hierachy calssification for hotspot
# 0906 corrected
################
# use new mapping method
# change R plotting method to '/home/ljzhang/project/dmd_breakpoint_caller_v1.05.R'
################
# traceback.print_exc()
################
# 0908, changing retmap method by removing -g 3 option to maximize member for classification
################
# 0914 changing retmapping reference to dmd region
# be careful for remapping positions
# samtools faidx ucsc.hg19.fasta chrX:31132344-33234673 > ucsc.hg19.dmd.5kb.fasta
################
# there maybe small inversion for Ion001 sample, allow soft-clipping on both side to avoid this problem
# tmap mapall -f %s -r %s -s %s stage1 map4
################


#==============================================================================
# 20161123 modification:
# 	use apply_async instead of map
# 	use int cpuNbrs
#	try writing error info to a error log
# 20170313 modification:
# 	change error log to stdout
#==============================================================================


from __future__ import division
import shutil
import os,glob,sys,re
from multiprocessing import Pool
import subprocess as sp
from optparse import OptionParser
#from collections import OrderedDict
#from math import sqrt
from collections import Counter
import traceback

parser=OptionParser()

parser.add_option(
	'-D',
	'--raw-dir',
	dest='rawdir',
	help='dir of raw sequencing data, default: /home/ljzhang/project/dmd/0707_bam_test',
	default='/home/ljzhang/project/dmd/0707_bam_test'
	)
	
parser.add_option(
	'-T',
	'--tmp-dir',
	dest='workdir',
	help='working directory name for storing tempory files'
	)
	
parser.add_option(
	'-F',
	'--fas-dir',
	dest='fasdir',
	help='directory of fasta reference files, default: /home/ljzhang/data/hg19',
	default='/home/ljzhang/data/hg19'
	)
	
parser.add_option(
	'-P',
	'--cpu-nbr',
	dest='cpus',
	help='maximum cpus called, default set to raw bam file numbers'
	)
	
parser.add_option(
	'-E',
	'--exon-table',
	dest='exons',
	help='dmd.exons.txt, default: /home/ljzhang/data/dmdExons/dmd.exons.txt',
	default='/home/ljzhang/data/dmdExons/dmd.exons.txt'
	)
	
parser.add_option(
	'-R',
	'--r-script',
	dest='rscrpt',
	help='candidate r filtration script, default: /home/ljzhang/project/dmd_breakpoint_caller_v1.05.R',
	default='/home/ljzhang/project/dmd_breakpoint_caller_v1.05.R'
	)
	
parser.add_option(
	'-M',
	'--processing-mode',
	dest='mode',
	help='input file mode, 1 for single file mode, 0 for multiple file mode, default set to 0',
	default='0'
	)
	
parser.add_option(
	'-I',
	'--input-bam',
	dest='input',
	help='input bam file under single file mode'
	)
	
parser.add_option(
	'-Z',
	'--reprocessing-procedure',
	dest='rep',
	help='reprocessing soft-clipping procedure,1: ture 0: false, default set to true',
	default='1'
	)

(options,args)=parser.parse_args()

if not options.rawdir or not options.workdir:
	parser.print_help()
	sys.exit(1)

if options.mode=='1' and not options.input:
	parser.print_help()
	sys.exit(1)

rawdir=options.rawdir
workdir=rawdir+'/'+options.workdir
fadir=options.fasdir
cpus=(int(options.cpus) if options.cpus else 1)

if not os.path.exists(workdir):
	os.mkdir(workdir)

chrxfa=fadir+'/ucsc.hg19.dmd.5kb.fasta'

def main():
	global rawdir,workdir,fadir,chrxfa,cpus
	os.chdir(rawdir)
	files=glob.glob('*rawlib.basecaller.bam')
	
	if options.mode=='0':
		print('main starts')
#		cpus=(int(options.cpus) if options.cpus else len(files))
		pool=Pool(cpus)
		for subfile in files:
			pool.apply_async(mainExc,args=(subfile,rawdir,workdir,fadir,chrxfa))
		pool.close()
		pool.join()
	else:
		## mainExc(filename=files[1])
		mainExc(options.input,rawdir,workdir,fadir,chrxfa)
	return
	
def mainExc(filename,rawdirMe,workdirMe,fadirMe,chrxfaMe):
	'''
	prepare for sorted.chrX.sam and DMD depth file
	global rawdirMe,workdirMe,fadirMe
	'''
	print('mainExc STARTS')
	
	os.chdir(rawdirMe)
	rawPat=re.compile(r'(.*?rawlib.basecaller)')
	rawname=rawPat.findall(filename)[0]
	
	print(filename+' raw samtools index starts')
	
	if not os.path.exists(filename+'.bai'):
		cmd='samtools index %s'	%	(filename)
		try:
			sp.call(cmd,shell=True)
		except Exception as error:
			print error
			pass
		
	print(filename+' tmap mapall for raw bam starts')
		
	if not os.path.exists(workdirMe+'/mapped.'+filename):
		cmd='tmap mapall -f %s -r %s -s %s stage1 map4'	%	(fadirMe+'/ucsc.hg19.fasta',filename,workdirMe+'/mapped.'+filename)
		try:
			sp.call(cmd,shell=True)
		except Exception as error:
			print error
			pass
		
	print('changed to workdir')
	os.chdir(workdirMe)
	
	if not os.path.exists(workdirMe+'/sortTmp'):
		os.mkdir(workdirMe+'/sortTmp')
	
	print('samtools sort for raw bam')
	if not os.path.exists(rawname+'.sorted.bam'):
		cmd='samtools sort -T ./sortTmp/%s -o %s %s'	%	(rawname,rawname+'.sorted.bam','mapped.'+filename)
		try:
			sp.call(cmd,shell=True)
		except Exception as error:
			print error
			pass
		
	print('samtools index for sorted bam')
	
	cmd='samtools index %s'	%	(rawname+'.sorted.bam')
	try:
		sp.call(cmd,shell=True)
	except Exception as error:
		print error
		pass
		
	print('samtools cut chrX info')
	if not os.path.exists(rawname+'.chrX.sorted.sam'):
		cmd='samtools view %s chrX > %s'	%	(rawname+'.sorted.bam',rawname+'.chrX.sorted.sam')
		try:
			sp.call(cmd,shell=True)
		except Exception as error:
			print error
			pass
	
	### getSoftclipStats(samfile=rawname+'.chrX.sorted.sam')
	
	## processing starts here, no judgement added
	if options.rep=='1':
		try:
			getSoftclipStats(samfileGss=rawname+'.chrX.sorted.sam',workdirGss=workdirMe)
		except Exception:
			traceback.print_exc()
			pass
#		getSoftclipStats(samfileGss=rawname+'.chrX.sorted.sam',workdirGss=workdirMe)
		
		if not os.path.exists(workdirMe+'/reTmapdir'):
			os.mkdir(workdirMe+'/reTmapdir')
		
		try:
			reTmap(stsfile=rawname+'.chrX.sorted.sam.sft.statistics',chrxfaRt=chrxfaMe,workdirRt=workdirMe)
		except Exception:
			traceback.print_exc()
			pass
#		reTmap(stsfile=rawname+'.chrX.sorted.sam.sft.statistics',chrxfaRt=chrxfaMe,workdirRt=workdirMe)
	
	if not os.path.exists(workdirMe+'/reTmapdir'):
		os.mkdir(workdirMe+'/reTmapdir')
	
	print('changed to workdir/reTmapdir')
	os.chdir(workdirMe+'/reTmapdir')
	
	# changing log starts here
	try:
		hierachyMain(filename=rawname+'.chrX.sorted.sam.sft.statistics.retmapped')
	except Exception:
		traceback.print_exc()
		pass
#	hierachyMain(filename=rawname+'.chrX.sorted.sam.sft.statistics.retmapped')
		
	try:
		potentFiltra(reshapedFile=rawname+'.chrX.sorted.sam.sft.statistics.retmapped.posReshaped')
	except Exception:
		traceback.print_exc()
		pass
#	potentFiltra(reshapedFile=rawname+'.chrX.sorted.sam.sft.statistics.retmapped.posReshaped')
		
	try:
		classifier(filenameClas=rawname+'.chrX.sorted.sam.sft.statistics.retmapped.posReshaped.potentialDel')
	except Exception:
		traceback.print_exc()
		pass
#	classifier(filenameClas=rawname+'.chrX.sorted.sam.sft.statistics.retmapped.posReshaped.potentialDel')
		
	try:
		classifier(filenameClas=rawname+'.chrX.sorted.sam.sft.statistics.retmapped.posReshaped.potentialDup')
	except Exception:
		traceback.print_exc()
		pass
#	classifier(filenameClas=rawname+'.chrX.sorted.sam.sft.statistics.retmapped.posReshaped.potentialDup')
		
	try:
		getMedian(filenameGm=rawname+'.chrX.sorted.sam.sft.statistics.retmapped.posReshaped.potentialDel.mergedInfo')
	except Exception:
		traceback.print_exc()
		pass
#	getMedian(filenameGm=rawname+'.chrX.sorted.sam.sft.statistics.retmapped.posReshaped.potentialDel.mergedInfo')
	
	try:
		getMedian(filenameGm=rawname+'.chrX.sorted.sam.sft.statistics.retmapped.posReshaped.potentialDup.mergedInfo')
	except Exception:
		traceback.print_exc()
		pass
#	getMedian(filenameGm=rawname+'.chrX.sorted.sam.sft.statistics.retmapped.posReshaped.potentialDup.mergedInfo')
		
	try:
		mapRegion(filenameMapR=rawname+'.chrX.sorted.sam.sft.statistics.retmapped.posReshaped.potentialDel.mergedInfo.InfoSum')
	except Exception:
		traceback.print_exc()
		pass
#	mapRegion(filenameMapR=rawname+'.chrX.sorted.sam.sft.statistics.retmapped.posReshaped.potentialDel.mergedInfo.InfoSum')
		
	try:
		mapRegion(filenameMapR=rawname+'.chrX.sorted.sam.sft.statistics.retmapped.posReshaped.potentialDup.mergedInfo.InfoSum')
	except Exception:
		traceback.print_exc()
		pass
#	mapRegion(filenameMapR=rawname+'.chrX.sorted.sam.sft.statistics.retmapped.posReshaped.potentialDup.mergedInfo.InfoSum')
	
	## axrange=[31137344-500, 33229673+500]
	
	print('samtools dmd depth')
	# if not os.path.exists(rawname+'.sorted.bam.dmddepth'):
	cmd='samtools depth %s -r chrX:%d-%d > %s'	%	('../'+rawname+'.sorted.bam',31137344-8000, 33229673+8000, rawname+'.sorted.bam.dmddepth')
	try:
		sp.call(cmd,shell=True)
	except Exception as error:
		traceback.print_exc()
		pass
	
	## disable R filtration for now
	rpat=re.compile(r'.*?(dmd_breakpoint.*?\.R)')
	
	try:
		localRscpt=rpat.findall(options.rscrpt)[-1]
	except Exception:
		traceback.print_exc()
		pass
	
	# if not os.path.exists(localRscpt):
		# shutil.copy(options.rscrpt, localRscpt)
	
	shutil.copy(options.rscrpt, localRscpt)
	
	dupfilenameforR=rawname+'.chrX.sorted.sam.sft.statistics.retmapped.posReshaped.potentialDup.mergedInfo.InfoSum.exonAnno'
	
	try:
		sp.call('R --vanilla --slave < %s %s %s'	%	(localRscpt,dupfilenameforR,'DUP'),shell=True)
	except Exception:
		traceback.print_exc()
		pass
		
	delfilenameforR=rawname+'.chrX.sorted.sam.sft.statistics.retmapped.posReshaped.potentialDel.mergedInfo.InfoSum.exonAnno'
	
	try:
		sp.call('R --vanilla --slave < %s %s %s'	%	(localRscpt,delfilenameforR,'DEL'),shell=True)
	except Exception:
		traceback.print_exc()
		pass
	
	return
	

def getSoftclipStats(samfileGss,workdirGss):
	'''
	generate .sft.statistics files
	'''
	os.chdir(workdirGss)
	
	###################################
	#### re compilation begin here ####
	softre=re.compile(r'.*?(\d+)S.*?')											# soft number
	insre=re.compile(r'.*?(\d+)I.*?')											# insertion number
	matre=re.compile(r'.*?(\d+)M.*?')											# match number
	padre=re.compile(r'.*?(\d+)P.*?')											# padding number
	mismre=re.compile(r'.*?(\d+)X.*?')											# mismatch number
	mat2re=re.compile(r'.*?(\d+)=.*?')											# type 2 match number
	
	softreb=re.compile(r'^(\d+)S')												# soft beginning
	softree=re.compile(r'(\d+)S$')												# soft ending
	###################################
	
	samfiled=workdirGss+"/"+samfileGss
	softclipping=workdirGss+"/"+samfileGss+".sft.statistics"
	
	stsout=open(softclipping,'w')
	
	samin=open(samfiled,'r')
	
	for line in samin.xreadlines():
		infoar=line.strip().split('\t')
		# if infoar[4] > qualthres:
		# mapping quality could be present for it would filter out mappings
		# basic info listed below
		
		# stsout.write('info\tflag\tchr\tpos\tqual\tCIGAR\traw_seq\tseq_qual\ttotal_soft_ratio\thead_soft_ratio\thead_soft_seq\ttail_soft_ratio\ttail_soft_seq\n')
		# no heading because there would be awk in reTmap
		
		insar=insre.findall(infoar[5])
		matar=matre.findall(infoar[5])
		padar=padre.findall(infoar[5])
		misar=mismre.findall(infoar[5])
		mat2ar=mat2re.findall(infoar[5])
		
		softar=softre.findall(infoar[5])
		softnbr=sum([int(i) for i in softar])
		
		########### set soft-clipping length filtration here ###########
		if softnbr<20:
			continue
		########### set soft-clipping length filtration here ###########
		
		stsout.write('\t'.join(infoar[0:6])+'\t')
		stsout.write('\t'.join(infoar[9:11])+'\t')
		
		insnbr=sum([int(i) for i in insar])
		matnbr=sum([int(i) for i in matar])
		padnbr=sum([int(i) for i in padar])
		misnbr=sum([int(i) for i in misar])
		mat2nbr=sum([int(i) for i in mat2ar])
		
		totalnbr=sum([softnbr,insnbr,matnbr,padnbr,misnbr,mat2nbr])
		
		# total soft clipping ratio
		stsout.write(str(round(100*softnbr/totalnbr,3))+"\t")
		# soft clipping in beginning
		# note that pat should be searched in info[5]
		softbear=softreb.findall(infoar[5])
		softenar=softree.findall(infoar[5])
		
		if softbear:
			stsout.write(str(round(100*int(softbear[0])/totalnbr,3))+"\t")
			stsout.write(infoar[9][:int(softbear[0])]+"\t")
		else:
			stsout.write("NA\t-\t")
			
		# soft clipping in ending
		if softenar:
			stsout.write(str(round(100*int(softenar[0])/totalnbr,3))+"\t")
			stsout.write(infoar[9][-int(softenar[-1]):]+"\n")
		else:
			stsout.write("NA\t-\n")
	
	samin.close()
	stsout.close()
	return

def reTmap(stsfile,chrxfaRt,workdirRt):
	'''
	recieve .sft.statistics and generate tmp fa files for retmapping
	'''
	os.chdir(workdirRt)
	
	infh=open(stsfile,'r')
	stsremappedfile=workdirRt+'/reTmapdir/'+stsfile+'.retmapped'
	outfh=open(stsremappedfile,'w')
	
	# info    flag    chr     pos     qual    CIGAR   raw_seq seq_qual        total_soft_ratio        head_soft_ratio head_soft_seq   tail_soft_ratio tail_soft_seq
	
	# raw sequence
	rawtmpfa=stsfile+'.tmp.fa'
	cmd1="awk '{print \">\"$1,$7}' OFS='\n' %s > %s" % (stsfile, rawtmpfa)
	sp.call(cmd1,shell=True)
	
	# head soft sequence
	softhtmpfa=stsfile+'.sfth.tmp.fa'
	cmd2="awk '{print \">\"$1,$11}' OFS='\n' %s > %s" % (stsfile, softhtmpfa)
	sp.call(cmd2,shell=True)
	
	# tail soft sequence
	softttmpfa=stsfile+'.sftt.tmp.fa'
	cmd3="awk '{print \">\"$1,$13}' OFS='\n' %s > %s" % (stsfile, softttmpfa)
	sp.call(cmd3,shell=True)
	
	# raw sequence tmap
	tmpsam=stsfile+'.rawtmp.sam'
	tmapcmd1="tmap mapall -f %s -r %s -s %s stage1 map4" % (chrxfaRt, rawtmpfa, tmpsam)
	sp.call(tmapcmd1,shell=True)
	os.remove(rawtmpfa)
	
	# head soft tmap
	softhsam=stsfile+'.sfth.sam'
	tmapcmd2="tmap mapall -f %s -r %s -s %s stage1 map4" % (chrxfaRt, softhtmpfa, softhsam)
	sp.call(tmapcmd2,shell=True)
	os.remove(softhtmpfa)
	
	# tail soft tmap
	softtsam=stsfile+'.sftt.sam'
	tmapcmd3="tmap mapall -f %s -r %s -s %s stage1 map4" % (chrxfaRt,softttmpfa,softtsam)
	sp.call(tmapcmd3,shell=True)
	os.remove(softttmpfa)
	
	# rawdct
	i=0
	rawsamdict={}
	for line in open(tmpsam,'r').xreadlines():
		if line.startswith('@'):
			pass
		else:
			i+=1
			samar=line.strip().split('\t')
			# careful for there would be 0 for non-matches
			if samar[3]=='0':
				tmp=[samar[j] for j in (1,2,3,4,5)]
			else:
				tmp=[samar[1],'chrX',str(31132343+int(samar[3])),samar[4],samar[5]]
			rawsamdict[i]=tmp
	# shdct
	i=0
	softhdict={}
	for line in open(softhsam,'r').xreadlines():
		
		if line.startswith('@'):
			pass
		else:
			i+=1
			samar=line.strip().split('\t')
			if samar[3]=='0':
				tmp=[samar[j] for j in (1,2,3,4,5)]
			else:
				tmp=[samar[1],'chrX',str(31132343+int(samar[3])),samar[4],samar[5]]
			# softhdict[i]=[samar[j] for j in (1,2,3,4,5)]
			softhdict[i]=tmp
	# stdct
	i=0
	softtdict={}
	for line in open(softtsam,'r').xreadlines():
		if line.startswith('@'):
			pass
		else:
			i+=1
			samar=line.strip().split('\t')
			if samar[3]=='0':
				tmp=[samar[j] for j in (1,2,3,4,5)]
			else:
				tmp=[samar[1],'chrX',str(31132343+int(samar[3])),samar[4],samar[5]]
			# softtdict[i]=[samar[j] for j in (1,2,3,4,5)]
			softtdict[i]=tmp
	
	# final file generation
	outfh.write('info\tflag\tchr\tpos\tqual\tCIGAR\traw_seq\tseq_qual\ttotal_soft_ratio\thead_soft_ratio\thead_soft_seq\ttail_soft_ratio\ttail_soft_seq\trp1\trp2\trp3\trp4\trp5\tsh1\tsh2\tsh3\tsh4\tsh5\tst1\tst2\tst3\tst4\tst5\n')
	i=0
	for line in infh.xreadlines():
		i+=1
		# print(line+"\t%s\t%s\t%s\n" % ('\t'.join(rawsamdict[i]), '\t'.join(softhdict[i]), '\t'.join(softtdict[i])))
		try:
			outfh.write(line.strip()+"\t%s\t%s\t%s\n" % ('\t'.join(rawsamdict[i]), '\t'.join(softhdict[i]), '\t'.join(softtdict[i])))
		except Exception:
			pass
	infh.close()
	outfh.close()
	
	os.remove(tmpsam)
	os.remove(softhsam)
	os.remove(softtsam)
		
def hierachyMain(filename):
	infh=open(filename,'r')
	outfh=open(filename+'.posReshaped','w')
	
	headsft=re.compile(r'^(\d+)S')
	tailsft=re.compile(r'.*?(\d+)S$')
	
	outfh.write('Info\toriMap\tsfHmap\tsfTmap\tori5\tori3\tsfH5\tsfH3\tsfT5\tsfT3\toriCIGAR\treMapCIGAR\tFLAG\n')
	
	# consider H as 5' T as 3', which gives H always < T
	for line in infh.xreadlines():
		linear=line.strip().split('\t')
		# soft head clip
		if not linear[10]=='-' and linear[1]==linear[18]=='0':
			# if headsft.match(linear[5]):
			headSftLen=int(headsft.findall(linear[5])[0])
			oriH=int(linear[3])
			oriT=int(linear[3])+len(linear[6])-headSftLen
			remapH=int(linear[20])
			remapT=int(linear[20])+headSftLen
			outfh.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'	%	(linear[0],linear[3],linear[20],linear[25],oriH,oriT,remapH,remapT,'NA','NA',linear[5],linear[22],0))
		elif not linear[10]=='-' and linear[1]=='16' and linear[18]=='0':
			# if headsft.match(linear[5]):
			headSftLen=int(headsft.findall(linear[5])[0])
			oriH=int(linear[3])
			oriT=int(linear[3])+len(linear[6])-headSftLen
			remapH=int(linear[20])
			remapT=int(linear[20])+headSftLen
			outfh.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'	%	(linear[0],linear[3],linear[20],linear[25],oriH,oriT,remapH,remapT,'NA','NA',linear[5],linear[22],16))
			
		# soft tail clip
		elif not linear[12]=='-' and linear[1]==linear[23]=='0':
			tailSftLen=int(tailsft.findall(linear[5])[0])
			oriH=int(linear[3])
			oriT=int(linear[3])+len(linear[6])-tailSftLen
			remapH=int(linear[25])
			remapT=int(linear[25])+tailSftLen
			outfh.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'	%	(linear[0],linear[3],linear[20],linear[25],oriH,oriT,'NA','NA',remapH,remapT,linear[5],linear[27],0))
		elif not linear[12]=='-' and linear[1]=='16' and linear[23]=='0':
			tailSftLen=int(tailsft.findall(linear[5])[0])
			oriH=int(linear[3])
			oriT=int(linear[3])+len(linear[6])-tailSftLen
			remapH=int(linear[25])
			remapT=int(linear[25])+tailSftLen
			outfh.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'	%	(linear[0],linear[3],linear[20],linear[25],oriH,oriT,'NA','NA',remapH,remapT,linear[5],linear[27],16))
	
	infh.close()
	outfh.close()
	
	# potentFiltra(reshapedFile=filename+'.posReshaped')
	# classifier(filenameClas=filename+'.posReshaped.potentialDel')
	# classifier(filenameClas=filename+'.posReshaped.potentialDup')
	# getMedian(filenameGm=filename+'.posReshaped.potentialDel.mergedInfo')
	# getMedian(filenameGm=filename+'.posReshaped.potentialDup.mergedInfo')
	# mapRegion(filenameMapR=filename+'.posReshaped.potentialDel.mergedInfo.InfoSum')
	# mapRegion(filenameMapR=filename+'.posReshaped.potentialDup.mergedInfo.InfoSum')
	
def potentFiltra(reshapedFile):
	infh=open(reshapedFile,'r')
	potdel=open(reshapedFile+'.potentialDel','w')
	potdup=open(reshapedFile+'.potentialDup','w')
	
	potdelList=[]
	potdupList=[]
	
	ind=1
	for line in infh.xreadlines():
		if line.startswith('Info'):
			continue
		
		linear=line.strip().split()
		# no soft head
		## axrange=[31137344-500, 33229673+500]
		if linear[7]=='NA':
			
			if int(linear[5])<int(linear[8]) and 31137344<int(linear[8])<33229673 and 31137344<int(linear[5])<33229673:
				potdelList.append((int(linear[5]),int(linear[8]),ind))
				
			elif int(linear[5])>int(linear[8]) and 31137344<int(linear[8])<33229673 and 31137344<int(linear[5])<33229673:
				potdupList.append((int(linear[8]),int(linear[5]),ind))
				
			
		# no soft tail
		elif linear[9]=='NA':
			
			if int(linear[7])<int(linear[4]) and 31137344<int(linear[7])<33229673 and 31137344<int(linear[4])<33229673:
				potdelList.append((int(linear[7]),int(linear[4]),ind))
				
			elif int(linear[7])>int(linear[4]) and 31137344<int(linear[7])<33229673 and 31137344<int(linear[4])<33229673:
				potdupList.append((int(linear[4]),int(linear[7]),ind))
				
		ind+=1
			
	potdelList.sort()
	potdupList.sort()
	
	potdel.write('\n'.join('\t'.join(str(i) for i in potdelList[j]) for j in range(len(potdelList))))
	potdel.write('\n')
	
	potdup.write('\n'.join('\t'.join(str(i) for i in potdupList[j]) for j in range(len(potdupList))))
	potdup.write('\n')
	
	potdel.close()
	potdup.close()
	
	

class node:
	def __init__(self,x,y,index):
		self.x=int(x)
		self.y=int(y)
		self.index=index
		self.xp=int(x)//100
		self.yp=int(y)//100

def classifier(filenameClas,gap=50):
	infh=open(filenameClas,'r')
	
	outfh=open(filenameClas+'.mergedInfo','w')
	
	rawpat=re.compile(r'(.*?posReshaped).*?')
	rawname=rawpat.findall(filenameClas)[0]
	
	nodear=[]
	xpar=[]
	
	for line in infh.xreadlines():
		linear=line.strip().split('\t')
		newnode=node(linear[0],linear[1],linear[2])
		nodear.append(newnode)
		xpar.append(newnode.xp)
		
	infh.close()
	
	commonxpar=Counter(xpar).most_common()
	rawPosreshFh=open(rawname,'r')
	
	ind=1
	univdct={}
	for line in rawPosreshFh.xreadlines():
		if line.startswith('Info'):
			continue
		univdct[str(ind)]=line.strip()
		ind+=1
	rawPosreshFh.close()
	
	clsind=1
	for candxpar in commonxpar:
		tmpnodear=[i for i in nodear if i.xp==candxpar[0]]
		candxar=[i.x for i in tmpnodear]
		candx=Counter(candxar).most_common(1)[0][0]
		
		candyar=[i.y for i in tmpnodear if i.x==candx]
		candy=Counter(candyar).most_common(1)[0][0]
		
		tmpcount=candxpar[1]
		
		for subnode in tmpnodear:
			outfh.write(str(clsind)+'\t'+str(candx)+'\t'+str(candy)+'\t'+str(tmpcount)+'\t')
			outfh.write(univdct[subnode.index]+'\n')
			
		clsind+=1
		
	outfh.close()
	
	
def getMedian(filenameGm):
	'''
	get classifier center using median axis
	'''
	infh=open(filenameGm,'r')
	outfh=open(filenameGm+'.InfoSum','w')
	
	tmpar=[]
	# xtmpar=[]
	ind='1'
	counts=''
	
	outfh.write('supportedReads\tBreakPoint1\tBreakPoint2\tgapLength\tInfo\n')
	
	for line in infh.xreadlines():
		linear=line.strip().split('\t')
		if linear[0]==ind:
			tmpar.append((linear[1],linear[2],linear[16]))
			counts=linear[3]
			# xtmpar.append(linear[1])
		else:
			xtmpar=[i[0] for i in tmpar]
			commonx=Counter(xtmpar).most_common(1)[0][0]
			ytmpar=[i[1] for i in tmpar if i[0]==commonx]
			commony=Counter(ytmpar).most_common(1)[0][0]
			
			# exoninfo=getRegionInfo(start=int(commonx),stop=int(commony),dmdDct=exondct)
			
			xyintv=abs(int(commonx)-int(commony))
			
			outfh.write('%s\t%s\t%s\t%s\t'	%	(counts,str(commonx),str(commony),str(xyintv)))
			
			for item in tmpar:
				outfh.write('%s\t'	%	(str(item[2])+':'+item[0]+'-'+item[1]))
				
			outfh.write('\n')
			
			ind=linear[0]
			counts=linear[3]
			tmpar=[]
			tmpar.append((linear[1],linear[2],linear[16]))
	
	infh.close()
	outfh.close()
	
def mapRegion(filenameMapR):
	infh=open(filenameMapR,'r')
	outfh=open(filenameMapR+'.exonAnno','w')
	
	exonfh=open('/home/ljzhang/data/dmdExons/dmd.exons.txt','r')
	exondct={}
	# ind=0
	introname='INTRON-80-'
	introlist=[31135344]
	for line in exonfh.xreadlines():
		linear=line.strip().split('\t')
		exondct['EXON-'+linear[0].lstrip('EXON')]=[linear[1],linear[2]]
		
		introname=introname+linear[0].lstrip('EXON')
		introlist.append(linear[1])
		exondct[introname]=introlist
		introlist=[]
		
		introname='INTRON-'+linear[0].lstrip('EXON')+'-'
		introlist.append(linear[2])
		
	exonfh.close()
	
	outfh.write('exonAnno\tsupportedReads\tBreakPoint1\tBreakPoint2\tgapLength\tInfo\n')
	for line in infh.xreadlines():
		if line.startswith('supportedReads'):
			continue
		
		linear=line.strip().split('\t')
		exoninfo=getRegionInfo(start=int(linear[1]),stop=int(linear[2]),dmdDct=exondct)
		outfh.write('%s\t%s'	%	(exoninfo,line))
		
	infh.close()
	outfh.close()
	
def getRegionInfo(start,stop,dmdDct):
	
	exon1='NA'
	exon2='NA'
	
	for key in dmdDct.keys():
		rg1=int(dmdDct[key][0])-1
		rg2=int(dmdDct[key][1])+1
		
		if rg1<start<rg2:
			exon1=key
			
		if rg1<stop<rg2:
			exon2=key
		
	return exon1+'__'+exon2

# mainExec(filename='IonXpress_006_rawlib.basecaller.chrX.sorted.sam.sft.statistics.retmapped')

# def sumDept(stt,end,dct):
	
	# stt=int(stt)
	# end=int(end)
	
	# if stt>end:
		# tp=stt
		# stt=end
		# end=tp
	# else:
		# pass
	
	# depth=0
	# while stt<=end:
		# sttp=str(stt)
		# if dct.has_key(sttp):
			# depth+=int(dct[sttp])
		# stt+=1
	# return depth

# def bisearch(tgt,poslistb):
	##put bisearch here
	# headb=0
	# tailb=len(poslistb)
	
	# while tailb>headb:
		# mid=(headb+tailb)//2
		# try:
			# if poslistb[mid]<tgt:
				# if poslistb[mid+1]>tgt:
					# return mid
				# elif poslistb[mid+1]==tgt:
					# return mid+1
				# else:
					# headb=mid+1
			# elif poslistb[mid]>tgt:
				# if poslistb[mid-1]<tgt:
					# return mid-1
				# elif poslistb[mid-1]==tgt:
					# return mid-1
				# else:
					# tailb=mid-1
			# else:
				# return mid
		# except Exception as error:
			# print error
			# break
	# return -1
	
	
if __name__ == '__main__':
	main()
	
	
	
	

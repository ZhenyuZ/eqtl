#!/usr/bin/python -O
# Anuar Konkashbaev
# Modified one line to call plate with female samples only
'''

  Recalling multi plate GWAS dataset with different CQC thresholds

  Usage:
        anCallCrlmm.py \\
                        --cel-files CELFILES \\
                        --anno AFILE \\
                        --out OFILE \\
                        cel1 [cel2 [...]]
'''
################################################################################
# importing libraries
import sys,os,re
import struct,csv
import gzip,time
from optparse import OptionParser
################################################################################
# globals
cerr=sys.stderr.write
cout=sys.stdout.write
dels=re.compile(r'[\t ]')
bases=re.compile(r'[ACGT]')
basematch={'A':re.compile(r'A'),'C':re.compile(r'C'),'G':re.compile(r'G'),
           'T':re.compile(r'T')}
compl={"A":"T","T":"A","G":"C","C":"G"}
chr2int={'1':1, '2':2, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9, '10':10,
         '11':11, '12':12, '13':13, '14':14, '15':15, '16':16, '17':17, '18':18,
         '19':19, '20':20, '21':21, '22':22, '23':23, '24':24, '25':25, '26':26,
         'X':23, 'Y':24, 'XY':25, 'YX':25, 'MT':26, '0':0, 'NA':0, None:0,
         'chr1':1,'chr2':2, 'chr3':3, 'chr4':4, 'chr5':5, 'chr6':6, 'chr7':7,
         'chr8':8, 'chr9':9, 'chr10':10,'chr11':11, 'chr12':12, 'chr13':13,
         'chr14':14, 'chr15':15, 'chr16':16, 'chr17':17, 'chr18':18,'chr19':19,
         'chr20':20, 'chr21':21, 'chr22':22, 'chr23':23, 'chr24':24, 'chr25':25,
         'chr26':26,'chrX':23, 'chrY':24, 'chrXY':25, 'chrYX':25, 'chrMT':26,
         'chr0':0}
snpflip={'A':'T','T':'A','C':'G','G':'C'}
ab={'1':'A A','2':'A B','3':'B B'}
################################################################################
# classes
################################################################################
class anBedFileWriter(object):
    '''
      anBedFileWriter - write SNP-major bed file snp by snp
    '''
    def __init__(self,filename):
        self.open(filename)
    def open(self,filename):
        '''
          Initializing bed file writer
        '''
        self.fName=filename
        self.fout=open(filename,'w')
        self.fout.write(struct.pack('B',108))
        self.fout.write(struct.pack('B',27))
        self.fout.write(struct.pack('B',1))
        self.isit=True
    def close(self):
        '''
          Closing file
        '''
        self.fout.close()
        self.isit=False
    def writeSnp(self,data):
        '''
          writes SNP data to filename
          data - an array of genotypes encoded as 0 -> aa,
                                                  1 -> ab,
                                                  2 -> bb,
                                                  3 -> nn
        '''
        num=0
        nsubj=len(data)
        dgt=[0,2,3,1]
        for i in xrange(len(data)):
            if(data[i]<0 or data[i]>3):
                self.crap('wrong genotype at %d %d'%(i,data[i]))
            data[i]=dgt[data[i]]
        nbytes=int(nsubj/4)
        nextra=nsubj%4
        # pack into bytes
        i=0
        while(i<nbytes):
            x4=data[4*i]
            x3=data[4*i+1]
            x2=data[4*i+2]
            x1=data[4*i+3]
            # write a byte to a file
            c=(x1<<6)+(x2<<4)+(x3<<2)+x4
            self.fout.write(struct.pack("B",c))
            i=i+1
        if(nextra==1):
            x4=data[4*i]
            x3=0
            x2=0
            x1=0
            # write a byte to a file
            c=(x1<<6)+(x2<<4)+(x3<<2)+x4
            self.fout.write(struct.pack("B",c))
        elif(nextra==2):
            x4=data[4*i]
            x3=data[4*i+1]
            x2=0
            x1=0
            # write a byte to a file
            c=(x1<<6)+(x2<<4)+(x3<<2)+x4
            self.fout.write(struct.pack("B",c))
        elif(nextra==3):
            x4=data[4*i]
            x3=data[4*i+1]
            x2=data[4*i+2]
            x1=0
            # write a byte to a file
            c=(x1<<6)+(x2<<4)+(x3<<2)+x4
            self.fout.write(struct.pack("B",c))
        return num
    def crap(self,m):
        self.close()
        cerr(m+'\n')
        sys.exit(1)
################################################################################
# functions
################################################################################
def removeFile(fname):
    print 'Removing '+fname+'\n'
    if(os.path.isfile(fname)):
        os.unlink(fname)
    return True
################################################################################
def runCrlmm(cels,prefix,cdfName):
    print 'Calling %d cel-files\n'%(len(cels))
    rfile=prefix+'runcrlmm.r'
    fout=open(rfile,'w')
    #echo 'library(crlmm)' > runCrlmm.$celfiles.r
    #fout.write('.libPaths("/userhome/genegateRPackages/lib.2.15.1")\n')
    fout.write('library(crlmm)\n')
    #echo 'celFiles=as.matrix(read.table("'$celfiles'",header=FALSE))' >>
    #runCrlmm.$celfiles.r
    fout.write('celFiles=c("')
    fout.write('","'.join(cels))
    fout.write('")\n')
    #echo 'system.time(crlmmResult <- crlmm(celFiles, verbose = TRUE, cdfName=
    #"GenomeWideSNP_5"))' >> runCrlmm.$celfiles.r
    
    fout.write(r'crlmmResult <- crlmm(celFiles,gender=rep(2L, length(celFiles)), verbose=TRUE,cdfName="'+cdfName+'")'+'\n')
    #echo 'system.time(write.table(calls(crlmmResult),file="'$celfiles.calls'",
    #sep="\t"))' >> runCrlmm.$celfiles.r
    fout.write(r'write.table(calls(crlmmResult),file="')
    fout.write(prefix+r'.calls",sep="\t")'+'\n')
    #echo 'system.time(write.table(round(confs(crlmmResult)*10000)/10000.0,file=
    #"'$celfiles.confs'",sep="\t"))' >> runCrlmm.$celfiles.r
    fout.write(r'write.table(round(confs(crlmmResult)*10000)/10000.0,file="')
    fout.write(prefix+r'.confs",sep="\t")'+'\n')
    #echo 'system.time(write.table(crlmmResult[["batchQC"]],"'$celfiles.batchqc'"
    #,sep="\t"))' >> runCrlmm.$celfiles.r
    fout.write(r'write.table(crlmmResult[["batchQC"]],"')
    fout.write(prefix+r'.batchqc",sep="\t")'+'\n')
    #echo 'system.time(write.table(crlmmResult[["SNR"]],"'$celfiles.snr'",
    #sep="\t"))' >> runCrlmm.$celfiles.r
    fout.write(r'write.table(crlmmResult[["SNR"]],"')
    fout.write(prefix+r'.snr",sep="\t")'+'\n')
    fout.close()
    print 'Output prefix is '+prefix+'\n'
    #run r
    cmd='R --vanilla CMD BATCH '+rfile
    if(os.system(cmd)):
        sys.stderr.write('running R failed, check '+rfile+'.Rout\n')
        return True
    else:
        removeFile('.RData')
        removeFile(rfile+'.Rout')
    return False
################################################################################
def readAnnotation(fname):
    snpinfo={}
    # open file
    fin=open(fname,'r')
    # skip header
    line=fin.readline()
    while(line.startswith('#')):
        line=fin.readline()
    # read header
    heads=line.strip().replace('"','').split(',')
    # read SNPs
    annoReader=csv.reader(fin)
    rsnumlist={}
    for arrs in annoReader:
        # get basic info
        snpid=arrs[0]
        if(snpid.startswith('SNP_A') and chr2int.has_key(arrs[2])):
            rsnum=arrs[1]
            if(rsnumlist.has_key(rsnum)):
                continue
            else:
                rsnumlist[rsnum]=True
            # check on chromosome
            chrom=str(chr2int[arrs[2]])
            if(arrs[21]=='1' or arrs[5]=='1'):
                chrom='25'
            pos=arrs[3]
            # get genetic position
            try:
                gpos=str(int(float(arrs[11].split()[0])*1000)/1000.0)
            except:
                gpos='0.0'
            # make dict
            snpinfo[snpid]={'rsmap':'\t'.join([chrom,rsnum,gpos,pos])}
            if(arrs[4]=='+'):
                snpinfo[snpid]['a']=arrs[8]
                snpinfo[snpid]['b']=arrs[9]
                snpinfo[snpid]['1']=arrs[8]+' '+arrs[8]
                snpinfo[snpid]['2']=arrs[8]+' '+arrs[9]
                snpinfo[snpid]['3']=arrs[9]+' '+arrs[9]
            else:
                snpinfo[snpid]['a']=snpflip[arrs[8]]
                snpinfo[snpid]['b']=snpflip[arrs[9]]
                snpinfo[snpid]['1']=snpflip[arrs[8]]+' '+snpflip[arrs[8]]
                snpinfo[snpid]['2']=snpflip[arrs[8]]+' '+snpflip[arrs[9]]
                snpinfo[snpid]['3']=snpflip[arrs[9]]+' '+snpflip[arrs[9]]
            snpinfo[snpid]['0']='0 0'
            # print snpinfo[snpid]
            # break;
    # close file
    fin.close()
    return snpinfo
################################################################################
def crlmm2plink(fname,anno,gtThresh=0.999):
    # check arguments
    pfiles=[]
    # process calls and confs file
    fnamecalls=fname+'.calls'
    fnameconfs=fname+'.confs'
    ## open files
    fcalls=open(fnamecalls,'r')
    fconfs=open(fnameconfs,'r')
    # open output bed file
    fbim=open(fname+'.bim','w')
    fbed=anBedFileWriter(fname+'.bed')
    # reading calls and confidences
    dudes=map(lambda x:x.replace('"',''),fcalls.readline().strip().split())
    confdudes=map(lambda x:x.replace('"',''),fconfs.readline().strip().split())
    for line in fcalls:
        calls=line.strip().split()
        confs=fconfs.readline().strip().split()
        snpid=calls[0].replace('"','')
        if(anno.has_key(snpid)):
            if(snpid!=confs[0].replace('"','')):
                return []
            fbim.write(anno[snpid]['rsmap']+'\t'+anno[snpid]['a']+'\t'+anno[snpid]['b']+'\n')
            # write snp to bed
            snpd=[]
            for i in xrange(1,len(calls)):
                if(float(confs[i])>gtThresh):
                    snpd.append(int(calls[i])-1)
                else:
                    snpd.append(3)
            fbed.writeSnp(snpd)
        else:
            if(snpid!=confs[0].replace('"','')):
                return []
            fbim.write('0\t'+snpid+'\t0\t0\ta\tb\n')
            # write snp to bed
            snpd=[]
            for i in xrange(1,len(calls)):
                if(float(confs[i])>gtThresh):
                    snpd.append(int(calls[i])-1)
                else:
                    snpd.append(3)
            fbed.writeSnp(snpd)
    # close output bim and bed files
    fbed.close()
    fbim.close()
    # close input files
    fconfs.close()
    fcalls.close()
    # writing fam
    fout=open(fname+'.fam','w')
    fout.write('\n'.join([x+'\t'+x+'\t0'*4 for x in dudes])+'\n')
    fout.close()
    # removing files
    removeFile(fnamecalls)
    removeFile(fnameconfs)
    return [fname+'.bed',fname+'.bim',fname+'.fam']
################################################################################
# MAIN
if(__name__=='__main__'):
    # processing options
    ## populating options
    parser=OptionParser(__doc__)
    parser.add_option('-c','--cel-files',dest='celFiles',default=None,
                      help='a file with cel-file list',type='string')
    parser.add_option('-a','--anno',dest='afile',default=None,type='string',
                      help='An Affimetrix annotation file')
    parser.add_option('-o','--out',dest='ofile',default=None,type='string',
                      help='An output prefix')
    opt,celfiles=parser.parse_args()
    ## checing arguments
    if(opt.celFiles!=None and len(celfiles)<1):
        fin=open(opt.celFiles,'r')
        celfiles=[x.strip().split()[0] for x in fin if(not x.startswith('cel_files'))]
        fin.close()
    elif(len(celfiles)<1 and opt.celFiles==None):
        parser.error('CEL file list have to be specified either with --cel-files or on command line')
    elif(len(celfiles)>1 and opt.celFiles!=None):
        parser.error('Only one CEL file list should be specified')
    elif(opt.afile==None):
        parser.error('At this point annotation file is required ...')
    elif(not os.path.isfile(opt.afile)):
        parser.error('Annotation file '+opt.afile+'does not exist')
    cout('Calling using %d CEL files\n'%(len(celfiles)))
    ## run CRLMM
    if(opt.ofile!=None):
        prefix=opt.ofile+'.crlmm'
    else:
        if(opt.celFiles!=None):
            prefix=opt.celFiles+'.crlmm'
        else:
            prefix='crlmm'
    cerr('\n'.join(celfiles)+'\n')
    if(runCrlmm(celfiles,prefix,'GenomeWideSNP_6')):
        cerr('CRLMM failed\n')
        sys.exit(1)
    ## convert to binary and rename duplicates
    cout('Converting to plink\n')
    # read annotation
    cout('Reading annotation from '+opt.afile+'\n')
    anno=readAnnotation(opt.afile)
    bfiles=crlmm2plink(prefix,anno)
    if(len(bfiles)!=3):
        cerr('CRLMM to plink failed\n')
        sys.exit(1)
    sys.exit()

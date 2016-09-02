import smtplib
import os
import fnmatch
from email.mime.text import MIMEText
from xml.dom.minidom import parse

DATA_DIRS=[
    '/data1/illumina_data/',
    #'/data1/illumina_data/NextSeqData/'
]


#################FUNCTIONS###################

def ConvertBcl(run_dir):

    LOG=run_dir+'/conversionLog.txt';
    ERROR=run_dir+'/conversionError.txt';

    FCID = os.path.split(run_dir)[-1].split("_")[-1]

    os.system("touch "+run_dir+"/ConversionDone.txt")
    os.system("date >>"+run_dir+'/conversionLog.txt')
    os.system("/data/fastqconversion/bcl2fastq_v2.17.1.14/bin/bcl2fastq --runfolder-dir %s -r 18 -d 18 -p 18 -w 4 -l WARNING 1>>%s 2>>%s" % (run_dir, LOG, ERROR))


    #ADD FLOWCELL ID TO FASTQ NAME
    for root, dirs, files in os.walk(run_dir+"/Data/Intensities/BaseCalls"):
        for name in files:
            if fnmatch.fnmatch(name, '*fastq.gz'):
                name_parts = name.split('_')
                if name_parts[1] != FCID:
                    newname = name.split('_')[0]+'_'+FCID+'_'+'_'.join(name.split('_')[1:])
                    os.rename( os.path.join(root, name), os.path.join(root, newname) )




    os.system("date >>"+run_dir+'/conversionLog.txt')
    run_parameters = parse(run_dir+"/runParameters.xml")
    experiment_name = run_parameters.getElementsByTagName('ExperimentName')[0].firstChild.nodeValue

    #msg = MIMEText('You can find the data for'+experiment_name+' in '+run_dir)
    #msg['Subject'] = 'Conversion of '+experiment_name+' finished'
    #msg['From'] = 'S.vanLieshout@hartwigmedicalfoundation.nl'
    #msg['To'] = 'S.vanLieshout@hartwigmedicalfoundation.nl'
    #s = smtplib.SMTP('smtp.gmail.com')
    #s.set_debuglevel(1)
    #s.sendmail('S.vanLieshout@hartwigmedicalfoundation.nl',['S.vanLieshout@hartwigmedicalfoundation.nl'], msg.as_string())
    #s.quit()
#




def SubDirPath (d):
    return filter(os.path.isdir, [os.path.join(d,f) for f in os.listdir(d)])

#ConvertBcl('/data_hodor/illumina_data/NEXTSEQ500_HU02/150529_NS500414_0094_HCWVTBGXX')

#############ZE MAGIC#####################
for d in DATA_DIRS:
    for run_dir in SubDirPath(d):
        if(os.path.isfile(run_dir+'/SampleSheet.csv') and os.path.isfile(run_dir+'/RTAComplete.txt') and not os.path.isfile(run_dir+'/ConversionDone.txt')):
            ConvertBcl(run_dir)

# coding: utf-8

# Contact L.Samoylova <liubov.samoylova@xfel.eu>, A.Buzmakov <buzmakov@gmail.com>
# SPB S2E simulation project, European XFEL Hamburg <www.xfel.eu>
# May 2014
# Wave optics software is based on 
# SRW core library <https://github.com/ochubar/SRW>, and 
# WPG framework <https://github.com/samoylv/WPG>

# This is a fake pulses generator for Demo version only. Realy return only 3fs_20 pulse (AB)
# In[ ]:


import shutil
import os
import errno

def mkdir_p(path):
    """
    Create directory tree, if not exists (mkdir -p)
    
    :param path: Path to be created
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-o", "--output-dir", dest="out_dir", help="Output directory")
    

    (options, args) = parser.parse_args()
       
    if not options.out_dir:   # if output directory not given
        parser.error('Output directory is not specified')
    else:
        out_dir=options.out_dir

    mkdir_p(out_dir)    
    shutil.copyfile('FELsource_out_0000001.h5', 
        os.path.join(out_dir,'FELsource_out_0000001.h5'))
# In[ ]:

if __name__ == '__main__':
    main()





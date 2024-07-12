import numpy as np
import copy
import os
from .matter_powerspectrum import *
import cosmopower

base_dir = os.path.join(os.path.dirname(__file__))

print(base_dir)
def download_data(download_dir):
    """
    Download the data needed for the emulators to the specified directory.

    Parameters
    ----------
    download_dir : str
        the data will be downloaded to this directory
    """
    from six.moves import urllib
    import shutil
    import gzip
    import glob
    import gdown

    #url = 'https://drive.google.com/file/d/1xup9pe2Ro2vGKWih2Gx9zuXZKn1WRaLU/view?usp=sharing' #v1
    #url = 'https://drive.google.com/file/d/14tF2GKYOSo5RVRGCSmW3WwxS4bz5XaZv/view?usp=sharing' #v2
    #url = 'https://drive.google.com/file/d/1Bmr4UpA69RTt15Vbu69I7D7ghUUEVSMP/view?usp=sharing' #v3
    filenames = 'models_v3.zip'
    version = '_v3'
    out_filenames = 'models'+version
    file_path = os.path.join(download_dir, filenames)
    final_path = os.path.join(download_dir, out_filenames)

    #url = 'https://drive.google.com/uc?id=1xup9pe2Ro2vGKWih2Gx9zuXZKn1WRaLU' #v1
    #url = 'https://drive.google.com/uc?id=14tF2GKYOSo5RVRGCSmW3WwxS4bz5XaZv' #v2
    url = 'https://drive.google.com/uc?id=1Bmr4UpA69RTt15Vbu69I7D7ghUUEVSMP' #v3

    # do not re-download
    if not os.path.exists(final_path):

        print("\n As it is the first instance of the emulator, "
                    "we need to download some data, it can take a few "
                    "seconds...\n")

        print("Downloading %s...\n" % out_filenames)

        #file_path, _ = urllib.request.urlretrieve(url=url,
        #                                            filename=file_path,
        #                                            reporthook=None)
        gdown.download(url, file_path, quiet=False)

        print("Download finished. Extracting files.")

        # unzip the file
        shutil.unpack_archive(
            filename=file_path, extract_dir=download_dir)
        os.remove(file_path)
        print("Done.\n")


download_data(base_dir)

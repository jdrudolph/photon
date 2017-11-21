import urllib.request
import gzip
import tarfile
from shutil import copyfile, copytree
import pytest
import tempfile
import os
import os.path

@pytest.fixture()
def clean_dir_with_data(request):
    newpath = tempfile.mkdtemp()
    copyfile(request.param, os.path.join(newpath, 'data.csv'))
    copytree('static', os.path.join(newpath, 'static'))
    copytree('templates', os.path.join(newpath, 'templates'))
    if not os.path.isfile('db.tar.gz'):
        url = 'http://cs.tau.ac.il/~jandanielr/db.tar.gz'
        print('downloading data')
        urllib.request.urlretrieve(url, 'db.tar.gz')
    db_path = os.path.abspath('db.tar.gz')
    os.chdir(newpath)
    with tarfile.open(db_path) as uncompressed:
            print('extracting data')
            uncompressed.extractall()

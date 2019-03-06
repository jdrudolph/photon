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
    copyfile('db.tar.gz', os.path.join(newpath, 'db.tar.gz'))
    os.chdir(newpath)
    with tarfile.open('db.tar.gz') as uncompressed:
            print('extracting data')
            uncompressed.extractall()

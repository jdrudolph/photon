#!/usr/bin/env python
import os
import io
import zipfile
from glob import glob
from flask import Flask, url_for, redirect, render_template, request, jsonify, send_from_directory, send_file
from phos.defaults import ROOT_DIR, WORK_DIR, parameters, descriptions

app = Flask(__name__.split('.')[0],
        template_folder=os.path.join(ROOT_DIR, 'templates'),
        static_folder=os.path.join(ROOT_DIR, 'static'))
app.config['CELERY_BROKER_URL'] = 'redis://localhost:6379/0'
app.config['CELERY_RESULT_BACKEND'] = 'redis://localhost:6379/0'
app.config['CELERY_ACCEPT_CONTENT'] = ['json']
app.config['CELERY_TASK_SERIALIZER'] = 'json'
app.config['CELERY_RESULT_SERIALIZER'] = 'json'

from celery import Celery
from celery.utils import uuid
celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)

import phos.pipeline

@celery.task()
def phos_task(task_id, data, _parameters):
    phos.util.deep_update(parameters, _parameters)
    return phos.pipeline.run(task_id, data, parameters)

@app.route('/', methods=['GET', 'POST'])
def index():
    parameters['anat']['anchor'] = None
    return render_template('index.html', parameters=parameters,
            descriptions=descriptions)

@app.route('/submit', methods=['POST'])
def submit():
    task_id = uuid()
    TASK_DIR = os.path.join(WORK_DIR, task_id)
    os.makedirs(TASK_DIR)
    file = request.files['data']
    filename = os.path.join(TASK_DIR, 'data.tab')
    file.save(filename)
    form = request.form
    _parameters = {k:{} for k in parameters.keys()}
    for _k,v in form.items():
        k, kk = _k.split(':')
        _parameters[k][kk] = phos.util.num(v)
    phos_task.apply_async(args=(task_id, filename, _parameters), task_id=task_id)
    return render_template('submit.html', task_id=task_id)

@app.route('/get/<task_id>')
def show_result(task_id):
    return send_from_directory(os.path.join(WORK_DIR, task_id),
            filename='result.html')

@app.route('/download/<task_id>')
def download_result(task_id):
    memory_zip = io.BytesIO()
    with zipfile.ZipFile(memory_zip, 'w') as zf:
        for path in glob(os.path.join(WORK_DIR, task_id, '*')):
            rpath = os.path.relpath(path)
            arcname = os.path.join(*os.path.split(rpath)[1:])
            zf.write(rpath, arcname=arcname)
    memory_zip.seek(0)
    return send_file(memory_zip, attachment_filename='{}.zip'.format(task_id),
            as_attachment=True)

@app.route('/status/<task_id>')
def status(task_id):
    task = phos_task.AsyncResult(task_id)
    return jsonify({'state' : task.state})

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='PHOTON server')
    parser.add_argument('--debug', action='store_const',
            const=True, default=False, help='start server in debug mode. \
                    Do not set this on a production system. \
                    Allows for the execution of arbitrary code on your system')
    parser.add_argument('--host', type=str, default='0.0.0.0',
            help="server host. Change to '127.0.0.1' to be visible only on your computer")
    args = parser.parse_args()
    app.run(debug=args.debug, host=args.host)
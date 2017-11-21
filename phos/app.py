#!/usr/bin/env python
import os
from os.path import join, abspath, dirname
import io
import zipfile
from glob import glob
from flask import Flask, url_for, redirect, render_template, request, jsonify, send_from_directory, send_file
from phos.defaults import make_defaults, descriptions
import phos.util
defaults = make_defaults(abspath(dirname(dirname(__file__))))
template_dir = abspath(join(defaults['root'], 'templates'))
app = Flask(__name__,
        template_folder=template_dir,
        static_folder=abspath(join(defaults['root'], 'static')))
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

import json
import redis
def task_key(task_id):
    return 'task:{}:progress'.format(task_id)

redis_db = redis.StrictRedis()
def make_set_progress(task_id):
    def set_progress(message):
        redis_db.set(task_key(task_id), json.dumps(message))
    return set_progress

def get_progress(task_id):
    raw_message = redis_db.get(task_key(task_id))
    if not raw_message:
        return "no message found"
    message = json.loads(raw_message.decode('utf-8'))
    return message

@celery.task()
def phos_task(task_id, data, _parameters):
    parameters = phos.util.read_json(abspath(join(defaults['root'], 'parameters.json')))
    phos.util.deep_update(parameters, _parameters)
    return phos.pipeline.run(task_id, data, parameters, work_dir=defaults['work'], db=defaults['db'], template_dir=template_dir, set_progress = make_set_progress(task_id))

@app.route('/', methods=['GET', 'POST'])
def index():
    parameters = phos.util.read_json(abspath(join(defaults['root'], 'parameters.json')))
    parameters['anat']['anchor'] = None
    return render_template('index.html', parameters=parameters,
            descriptions=descriptions)

@app.route('/submit', methods=['POST'])
def submit():
    task_id = uuid()
    TASK_DIR = os.path.join(defaults['work'], task_id)
    os.makedirs(TASK_DIR)
    file = request.files['data']
    filename = os.path.join(TASK_DIR, 'data.tab')
    file.save(filename)
    form = request.form
    parameters = phos.util.read_json(abspath(join(defaults['root'], 'parameters.json')))
    _parameters = {k:{} for k in parameters.keys()}
    for _k,v in form.items():
        k, kk = _k.split(':')
        _parameters[k][kk] = phos.util.num(v)
    phos_task.apply_async(args=(task_id, filename, _parameters), task_id=task_id)
    return render_template('submit.html', task_id=task_id)

@app.route('/get/<task_id>')
def show_result(task_id):
    return send_from_directory(os.path.join(defaults['work'], task_id),
            filename='result.html')

@app.route('/download/<task_id>')
def download_result(task_id):
    memory_zip = io.BytesIO()
    with zipfile.ZipFile(memory_zip, 'w') as zf:
        for path in glob(os.path.join(defaults['work'], task_id, '*')):
            rpath = os.path.relpath(path)
            arcname = os.path.join(*os.path.split(rpath)[1:])
            zf.write(rpath, arcname=arcname)
    memory_zip.seek(0)
    return send_file(memory_zip, attachment_filename='{}.zip'.format(task_id),
            as_attachment=True)

@app.route('/status/<task_id>')
def status(task_id):
    task = phos_task.AsyncResult(task_id)
    info = str(task.info) if task.state == 'FAILURE' else ''
    message = get_progress(task_id)
    print('status info', info)
    return jsonify({'state' : task.state, 'info': info, 'message' : message})

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

import os
import re
import uuid
import time
import copy
import numpy
from threading import Thread
from flask import Flask, jsonify, render_template, request, session, Markup
from simplejson import JSONEncoder

import pyd2d

SESSION_LIFETIME = 3000
DEBUG = True
ROUND = 0

__author__ = 'Clemens Blank'
app = Flask(__name__)

class JSONEncoder_ignore(JSONEncoder):
    def __init__(self, **kwargs):
        kwargs['ignore_nan'] = True
        super(JSONEncoder_ignore, self).__init__(**kwargs)

try:
    from flask.ext.compress import Compress
    app.config['COMPRESS_LEVEL'] = 4
    Compress(app)
except:
    pass

@app.route('/', methods=['GET'])
def d2d_presenter():
    data = {}
    data['xvalue'] = [0.342300000, float('nan'), float('inf'), float('-inf'), 7]
    return jsonify(data=data)

if __name__ == "__main__":
    app.secret_key = os.urandom(24)
    app.permanent_session_lifetime = SESSION_LIFETIME
    app.debug = DEBUG
    app.json_encoder = JSONEncoder_ignore
    app.run(threaded=True,
            host="127.0.0.1",
            port=int("5000")
            )

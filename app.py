import os

from flask import Flask, render_template, request, jsonify
from werkzeug.utils import secure_filename

app = Flask(__name__)

basedir = os.path.abspath(os.path.dirname(__file__))
app.config['ALLOWED_EXTENSIONS'] = {'txt'}


@app.route('/')
def hello_world():
    return render_template('index.html')


@app.route('/upload_file', methods=['POST'])
def upldfile():
    # https://github.com/moremorefor/flask-fileupload-ajax-example
    if request.method == 'POST':
        files = request.files['file']
        if files and allowed_file(files.filename):
            filename = secure_filename(files.filename)
            app.logger.info('FileName: ' + filename)
            updir = os.path.join(basedir, 'upload/')
            files.save(os.path.join(updir, filename))
            return jsonify(filename=filename)
        else:
            return jsonify(filename='wrong extension')


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']


if __name__ == '__main__':
    app.run()

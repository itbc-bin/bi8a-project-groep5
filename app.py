import os

from flask import Flask, render_template, request, jsonify, send_file
from werkzeug.utils import secure_filename

app = Flask(__name__)

basedir = os.path.abspath(os.path.dirname(__file__))
app.config['ALLOWED_EXTENSIONS'] = {'txt'}

test_data = [
    {'zoekwoord': 'google', 'aantal': 5, 'link': 'https://google.com'},
    {'zoekwoord': 'facebook', 'aantal': 3,
     'link': 'https://facebook.com'},
    {'zoekwoord': 'twitter', 'aantal': 8,
     'link': 'https://twitter.com'}
]


@app.route('/')
def hello_world():
    global test_data

    return render_template('index.html', test_data=test_data)


@app.route('/download', methods=['GET'])
def download():
    global test_data
    with open('data.tsv', 'w') as file:
        file.write('zoekwoord\taantal\tlink\n')
        for data in test_data:
            file.write(
                f'{data["zoekwoord"]}\t{data["aantal"]}\t{data["link"]}\n')
    return send_file('data.tsv',
                     mimetype='text/csv',
                     attachment_filename='data.tsv',
                     as_attachment=True)


@app.route('/upload_file', methods=['POST'])
def upldfile():
    # https://github.com/moremorefor/flask-fileupload-ajax-example
    if request.method == 'POST':
        files = request.files['file']
        if files and allowed_file(files.filename):
            filename = secure_filename(files.filename)
            app.logger.info('FileName: ' + filename)
            updir = os.path.join(basedir, 'upload/')
            old_file = get_old_gen_panel_file()
            file = os.path.join('upload', old_file)
            old_dir = os.path.join(os.getcwd(), f'upload{os.path.sep}old_files')
            os.rename(os.path.join(os.getcwd(), file), os.path.join(old_dir, old_file))
            files.save(os.path.join(updir, filename))
            return jsonify(filename=filename)
        else:
            return jsonify(filename='wrong extension')


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']


def get_old_gen_panel_file():
    for file in os.listdir(os.path.join(os.getcwd(), 'upload')):
        if file.endswith('.txt'):
            return file


if __name__ == '__main__':
    app.run(debug=True)

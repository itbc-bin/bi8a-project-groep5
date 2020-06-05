import os
import random
import re
import string
import threading
import time

import mysql.connector
from flask import Flask, render_template, request, jsonify, send_file, json
from werkzeug.utils import secure_filename

from co_occurrence_algorithm.co_occurrence import CoOccurrence
from pubmed_search.pubmed_searcher import PubmedSearch

app = Flask(__name__)
app.config['SECRET_KEY'] = 'CnzOd54-fbuNu_X3_-PDzQ'
app.config['ALLOWED_EXTENSIONS'] = {'txt'}

basedir = os.path.abspath(os.path.dirname(__file__))

result_ids = []
articles = []


@app.errorhandler(404)
def error_page_not_found(error):
    return render_template('error_pages/404.html')


@app.route('/')
def home_page():
    term = ''
    symbols = None
    first_visit = True
    results = []
    ids_data = []
    if request.args.get("pheno_input"):
        first_visit = False
        term = request.args.get("pheno_input")
        words = request.args.getlist("symbols_input")
        symbols = [symbol.strip() for symbol in words[0].split()]
        date = request.args.get('calendar_input')
        if date:
            search = PubmedSearch(search_word=term,
                                  gene_symbols=words[0],
                                  date=date)
        else:
            search = PubmedSearch(search_word=term,
                                  gene_symbols=words[0])
        search.search_pubmed()
        results = search.results
        ids_data = search.ids_data
        thread = threading.Thread(target=parse_results, args=(search,))
        thread.daemon = True
        thread.start()
    return render_template('home.html', results=results,
                           first_visit=first_visit, ids_data=ids_data,
                           term=term, symbols=symbols)


@app.route('/about_project')
def about_project_page():
    return render_template('about_project.html')


@app.route('/about_makers')
def about_makers_page():
    return render_template('about_makers.html')


@app.route('/download', methods=['GET'])
def download():
    articles_data = request.args.get('results')
    articles_data = articles_data.replace("\'", "\"")
    articles_data = json.loads(articles_data)
    download_dir = os.path.join(basedir, 'data_files')
    with open(os.path.join(download_dir, 'data.tsv'), 'w') as file:
        file.write('Search Word\tAmount of hits\tlink\n')
        for article in articles_data:
            file.write(
                f'{article["search_word"]}\t{article["amount_hits"]}'
                f'\t{article["link"]}\n')
    return send_file(os.path.join(download_dir, 'data.tsv'),
                     mimetype='text/csv',
                     attachment_filename=os.path.join(download_dir,
                                                      'data.tsv'),
                     as_attachment=True)


@app.route('/upload_file', methods=['POST'])
def upload_file():
    # https://github.com/moremorefor/flask-fileupload-ajax-example
    if request.method == 'POST':
        files = request.files['file']
        if files and allowed_file(files.filename):
            filename = secure_filename(files.filename)
            app.logger.info('FileName: ' + filename)
            updir = os.path.join(basedir, 'data_files')
            old_file = get_old_gen_panel_file()
            if old_file:
                file = os.path.join('data_files', old_file)
                old_dir = os.path.join(basedir,
                                       os.path.join('data_files', 'old_files'))
                os.rename(os.path.join(basedir, file),
                          os.path.join(old_dir, old_file))

            files.save(os.path.join(updir, filename))
            return jsonify(filename=filename)
        else:
            return jsonify(filename='Wrong extension')


@app.route('/results')
def results_page():
    url_results = get_all_previous_results()
    return render_template('results.html', url_results=url_results)


@app.route('/results/<result_id>', methods=['GET'])
def individual_result(result_id):
    result_list, title = get_algorithm_results(result_id)
    table_list = get_results(result_list)
    if result_list:
        return render_template('individual_result.html', title=title,
                               result_list=result_list, table_list=table_list)
    else:
        return render_template('error_pages/404.html'), 404


@app.route('/results/<result_id>', methods=['POST'])
def do_algorithm(result_id):
    global articles

    # wait until thread is ready inserting articles into database
    while not articles:
        time.sleep(2)
    # if result_id alreay taken, generate a new one
    while result_id in result_ids:
        result_id = f'''_{"".join(random.choice(
            string.ascii_lowercase + string.digits) for _ in range(9))}'''
    result_ids.append(result_id)
    url = f'/results/{result_id}'

    options = {'Title': False, 'Sentence': False, 'Abstract': False,
               'Multiple Abstracts': False}

    data = json.loads(request.form['data'])
    selected_options = data['options']
    job_title = data['jobTitle']
    term = data['term']

    for selected_option in selected_options:
        options[selected_option] = True
    co_occurrence = CoOccurrence(data=articles, url_id=result_id, term=term,
                                 title=job_title, in_title=options['Title'],
                                 in_sentence=options['Sentence'],
                                 in_abstract=options['Abstract'],
                                 in_multiple_abstracts=options[
                                     'Multiple Abstracts'])
    co_occurrence.pre_process_data()
    co_occurrence.calculate_co_occurrence()
    co_occurrence.get_co_occurence()
    co_occurrence.save_to_db()
    return json.dumps({'status': 'OK', 'url': url})


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']


def get_old_gen_panel_file():
    for file in os.listdir(os.path.join(basedir, 'data_files')):
        if file.endswith('.txt'):
            return file


def get_algorithm_results(result_id):
    results_list = []
    connction = connection_database()
    cursor = connction.cursor()
    cursor.execute("select combination, amount, pmid, title from results inner"
                   " join algorithm_results ar on results.id = ar.results_id "
                   "inner join pmid_result pr on ar.id = "
                   "pr.algorithm_results_id where url_id = '{}'"
                   .format(result_id))
    all_results = cursor.fetchall()
    if not all_results:
        return {}, ''
    title = all_results[0][-1]
    for result in all_results:
        results_list.append(
            {'combination': result[0],
             'amount': result[1],
             'PMID': result[2]}
        )
    return results_list, title


def parse_results(search):
    global articles
    search.parse_results()
    search.insert_to_database()
    articles = search.get_articles()


def get_results(result_list):
    connection = connection_database()
    cursor = connection.cursor()
    table_list = []

    for info in result_list:
        cursor.execute(
            "select title, abstract, keywords, authors, publication_year, "
            "article_link, pubmed_id from articles where "
            "pubmed_id = '{}'".format(info['PMID']))
        data = cursor.fetchall()[0]
        result = {
            'Title': get_clean_html(data[0]),
            'Abstract': get_clean_html(data[1]),
            'Keywords': data[2],
            'Authors': data[3],
            'Publication_year': data[4],
            'Link': data[5],
            'pmid': data[6]
        }
        table_list.append(result)
    connection.close()
    return table_list


def get_all_previous_results():
    connection = connection_database()
    cursor = connection.cursor()
    cursor.execute("select url_id, title, creation_date from results")
    url_results = cursor.fetchall()
    all_urls = [{
        'url': f'results/{result[0]}',
        'title': result[1],
        'date': result[2]
    } for result in url_results]
    connection.close()
    return all_urls


def connection_database():
    """Make connection to the database
    :return: The connection
    """
    connection = mysql.connector.connect(
        host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
        user='owe7_pg5@hannl-hlo-bioinformatica-mysqlsrv',
        database='owe7_pg5',
        password='blaat1234')

    return connection


def get_clean_html(raw_html):
    html_pattern = re.compile(r'<.*?>')
    clean_text = re.sub(html_pattern, '', raw_html)
    return clean_text


if __name__ == '__main__':
    app.run(debug=True)

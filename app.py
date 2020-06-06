import os
import random
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
    """
    Custom error handler page.
    :param error: Type of HTTP error
    :return: 404.html template.
    """
    return render_template('error_pages/404.html')


@app.route('/')
def home_page():
    """
    Homepage, gets user input and performs pubmed search. Starts thread to
    add results to database.
    :return: home.html template.
    """
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
            pubmed_search = PubmedSearch(search_word=term,
                                         gene_symbols=words[0],
                                         date=date)
        else:
            pubmed_search = PubmedSearch(search_word=term,
                                         gene_symbols=words[0])
        pubmed_search.search_pubmed()
        results = pubmed_search.results
        ids_data = pubmed_search.ids_data
        thread = threading.Thread(target=parse_results, args=(pubmed_search,))
        thread.daemon = True
        thread.start()
    return render_template('home.html', results=results,
                           first_visit=first_visit, ids_data=ids_data,
                           term=term, symbols=symbols)


@app.route('/about_project')
def about_project_page():
    """
    About project page.
    :return: about_project.html template.
    """
    return render_template('about_project.html')


@app.route('/about_makers')
def about_makers_page():
    """
    About makers page.
    :return: about_makers.html template.
    """
    return render_template('about_makers.html')


@app.route('/results')
def results_page():
    """
    Page with all of the algorithm results which are retrieved from database.
    :return: results.html template.
    """
    url_results = get_all_previous_results()
    return render_template('results.html', url_results=url_results)


@app.route('/results/<result_id>', methods=['GET'])
def individual_result(result_id):
    """
    Get an individual algorithm result based on id in url.
    :param result_id: id of the algorithm searched, which is generated at each
    search and added to database.
    :return: individual_results.html template.
    """
    gene_symbols = get_gene_symbols()
    result_list, title = get_algorithm_results(result_id)
    table_list = get_results(result_list)
    processed_results_list = process_results(result_list, gene_symbols)
    if result_list:
        return render_template('individual_result.html', title=title,
                               result_list=result_list, table_list=table_list,
                               processed_results_list=processed_results_list)
    else:
        return render_template('error_pages/404.html'), 404


@app.route('/results/<result_id>', methods=['POST'])
def do_algorithm(result_id):
    """
    Perform the co-occurrence algorithm and add results to database based on
    random generated id from url.
    :param result_id: Random generated id (via Javascript) for the search
    :return: json object with the status and url.
    """
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


@app.route('/download', methods=['GET'])
def download():
    """
    Download the data from the table. Request is made via Ajax in Javascript
    when user presses download button.
    :return: tsv file containing the data.
    """
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
    """
    Upload new gene panel file, and move old gene_panel file to old_files
    directory.
    :return: JSON object to let the user (via Javascript) know if the upload
    was succesful.
    """
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


def parse_results(pubmed_search):
    """
    Parse the articles of the pubmed search, and add to database. This is
    performed in another thread, so the user does not have to wait.
    :param pubmed_search: The PubmedSearch object containing the data.
    """
    global articles
    pubmed_search.parse_results()
    pubmed_search.insert_to_database()
    articles = pubmed_search.get_articles()


def get_all_previous_results():
    """
    Get all performed co-occurrence algorithm results from the databse.
    :return: A list with dictionaries containing the url, title and date of
    the algorithm search.
    """
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


def get_gene_symbols():
    """
    Get all of the gene symbols from the gene panel file. First check where
    the gene panel file is saved. If there is no file in data)_files, check the
    old_files directory.
    :return: List containing the gene symbols, if there is no file present,
    return empty list.
    """
    data_dir = os.path.join(basedir, 'data_files')
    gen_panel_file = None
    for _file in os.listdir(data_dir):
        if _file.endswith('.txt') and 'GenPanel' in _file:
            gen_panel_file = os.path.join(data_dir, _file)
            break
    if not gen_panel_file:
        for _file in os.listdir(os.path.join(data_dir, 'old_files')):
            if _file.endswith('.txt') and 'GenPanel' in _file:
                gen_panel_file = os.path.join(
                    os.path.join(data_dir, 'old_files'), _file)

    gene_symbols = []
    if gen_panel_file:
        with open(gen_panel_file, 'r') as symbols_file:
            symbols_file.readline()
            for line in symbols_file:
                gene_symbols.append(line.split('\t')[0].strip())

    return gene_symbols


def get_algorithm_results(result_id):
    """
    Get the combination, amount, pmid and title of the co-occurrence algorithm
    search from the database, and return the results.
    :param result_id: The id of the algorithm search.
    :return: List containing dictionarir with the combination, amount en pubmed
    id. And a string which is the title of the search.
    """
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


def get_results(result_list):
    """
    Get the data of the articles from the database, from a list containing
    multiple pubmed ids.
    :param result_list: A list with dictionaries, containing the the pubmed id.
    :return: List that can be used in individual_results.html with the data
    of the articles.
    """
    connection = connection_database()
    cursor = connection.cursor()
    table_list = []
    all_combinations = [info['combination'] for info in result_list]
    all_combinations = list(set(all_combinations))

    for info in result_list:
        cursor.execute(
            "select title, abstract, keywords, authors, publication_year, "
            "article_link, pubmed_id from articles where "
            "pubmed_id = '{}'".format(info['PMID']))
        data = cursor.fetchall()[0]
        title = make_bold(data[0], all_combinations)
        abstract = make_bold(data[1], all_combinations)
        result = {
            'Title': title,
            'Abstract': abstract,
            'Keywords': data[2],
            'Authors': data[3],
            'Publication_year': data[4],
            'Link': data[5],
            'pmid': data[6]
        }
        table_list.append(result)
    connection.close()
    return table_list


def make_bold(part, all_combinations):
    """
    Make the gene symbol and phenotype bold in the text. Which can be eiter
    the title or the abstract.
    :param part: Text of title or abstract.
    :param all_combinations: List of combinations that should be made bold in
    the text.
    :return: The text with the gene symbols and phenotypes replaced to be bold.
    """
    for combi in all_combinations:
        gene_symbol = combi.split(', ')[0]
        phenotype = combi.split(', ')[1]
        part = part.replace(gene_symbol, f'<strong>{gene_symbol}</strong>')
        part = part.replace(phenotype, f'<strong>{phenotype}</strong>')
    return part


def process_results(result_list, gene_symbols):
    """
    Process the results of the algorithm. Add pubmed ids which belong to the
    same combination to a new list of dictionaries where the pubmed ids are
    merged. Also add if the gene symbol from the combination is already in the
    genepanel file.
    :param result_list: List containing the combinations, amounts, and pubmed
    ids from the co-occurrence algorithm.
    :param gene_symbols: A list containing all the gene symbols from the
    gene panel file.
    :return: A list with dictionaries containing the combination, its amount,
    its pubmed ids and True/False whether or not in occurs in the genepanel
    file.
    """
    new_list = []
    for item in result_list:
        if combi_in_list(item['combination'], new_list):
            append_pmid(new_list, item['combination'], item['PMID'])
        else:
            data = {item['combination']: {
                'amount': item['amount'],
                'PMIDS': [str(item['PMID'])],
            }
            }
            if gene_symbols:
                data[item['combination']]['in_genepanel'] = \
                    item['combination'].split(',')[0] in gene_symbols

            new_list.append(data)
    return new_list


def combi_in_list(combination, _list):
    """
    Return true if a comnbination is already in the new list, false when the
    list does not contain the combination.
    :param combination: The combination.
    :param _list: The list containing all the combinations already added.
    :return: Boolean whether or not the current combination already occures in
    the list.
    """
    for element in _list:
        if combination in element:
            return True
    return False


def append_pmid(_list, combination, pmid):
    """
    Append a pubmed id to its corresponding combination in the list.
    :param _list: The list containing the combinations and its data.
    :param combination: The combination to which the pubmed id should be added.
    :param pmid: The pubmed that should be added.
    """
    for element in _list:
        try:
            element[combination]['PMIDS'].append(str(pmid))
        except KeyError:
            pass


def allowed_file(filename):
    """
    Helper function for upload file, only allow txt files to be uploaded.
    :param filename: The name of the file.
    :return: boolean whether or not the file ends with .txt.
    """
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']


def get_old_gen_panel_file():
    """
    Get the old gen_panel file after user uploaded new one, so it can be moved
    to the old_files directory.
    :return: The name of the old file.
    """
    for file in os.listdir(os.path.join(basedir, 'data_files')):
        if file.endswith('.txt'):
            return file


def connection_database():
    """Make connection to the database.
    :return: The connection.
    """
    connection = mysql.connector.connect(
        host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
        user='owe7_pg5@hannl-hlo-bioinformatica-mysqlsrv',
        database='owe7_pg5',
        password='blaat1234')

    return connection


if __name__ == '__main__':
    app.run(debug=True)

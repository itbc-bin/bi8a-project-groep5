import os
import random
import string
import threading

import mysql.connector
from flask import Flask, render_template, request, jsonify, send_file, json, \
    session
from werkzeug.utils import secure_filename

from co_occurrence_algorithm.co_occurrence import CoOccurrence
from pubmed_searcher import PubmedSearch

# TODO clean up code
# TODO write documentation
# TODO handle bugs
# TODO visualise co occurrence data
# TODO cleanup index.html
# TODO make better 404 error page


app = Flask(__name__)
app.config['SECRET_KEY'] = 'CnzOd54-fbuNu_X3_-PDzQ'
app.config['ALLOWED_EXTENSIONS'] = {'txt'}

basedir = os.path.abspath(os.path.dirname(__file__))

result_ids = []
results = []
test_data = [
    {'zoekwoord': 'google', 'aantal': 5, 'link': 'https://google.com'},
    {'zoekwoord': 'facebook', 'aantal': 3,
     'link': 'https://facebook.com'},
    {'zoekwoord': 'twitter', 'aantal': 8,
     'link': 'https://twitter.com'}
]
e_mail = ''


@app.route('/')
def home_page():
    global results
    global e_mail
    first_time = True
    articles = []
    if request.args.get("pheno_input"):
        first_time = False
        term = request.args.get("pheno_input")
        words = request.args.getlist("symbols_input")
        e_mail = request.args.get("input_mail")
        date = request.args.get('calendar_input')
        if date:
            search = PubmedSearch(e_mail=e_mail, search_word=term,
                                  gene_symbols=words[0],
                                  date=date)
        else:
            search = PubmedSearch(e_mail=e_mail, search_word=term,
                                  gene_symbols=words[0])
        search.search_pubmed()
        articles_data = search.articles_data
        articles = articles_data
        x = threading.Thread(target=parse_results, args=(search,))
        x.daemon = True
        x.start()
    session['articles'] = articles
    return render_template('index.html', articles=articles,
                           results=results, first_time=first_time)


@app.route('/download', methods=['GET'])
def download():
    articles_data = session.get('articles')
    with open('data_files/data.tsv', 'w') as file:
        file.write('zoekwoord\taantal\tlink\n')
        for article in articles_data:
            file.write(
                f'{article["zoekwoord"]}\t{article["aantal_hits"]}'
                f'\t{article["link"]}\n')
    return send_file('data_files/data.tsv',
                     mimetype='text/csv',
                     attachment_filename='data_files/data.tsv',
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
            old_dir = os.path.join(os.getcwd(),
                                   f'upload{os.path.sep}old_files')
            os.rename(os.path.join(os.getcwd(), file),
                      os.path.join(old_dir, old_file))
            files.save(os.path.join(updir, filename))
            return jsonify(filename=filename)
        else:
            return jsonify(filename='wrong extension')


@app.route('/results/<result_id>', methods=['GET'])
def render_results(result_id):
    result_dict = get_algorithm_results(result_id)
    if result_dict:
        print(result_dict)
        return render_template('results.html', title='test',
                               results=result_dict)
    else:
        return render_template('error_pages/404.html'), 404


@app.route('/results/<result_id>', methods=['POST'])
def do_algorithm(result_id):
    global results
    global e_mail
    print(result_id)
    options = {'Title': False, 'Sentence': False, 'Abstract': False,
               'Multiple Abstracts': False}
    test_data1 = [
        {
            'Title': "PRRT2 gene variant in a child with dysmorphic features, congenital microcephaly, and severe epileptic seizures: genotype-phenotype correlation?",
            'Authors': "Pavone Piero P, Corsello Giovanni G, Cho Sung Yoon SY, Pappalardo Xena Giada XG, Ruggieri Martino M, Marino Simona Domenica SD, Jin Dong Kyu DK, Marino Silvia S, Falsaperla Raffaele R",
            'Publication_year': "2019",
            'Keywords': "Dysmorphic features, Epileptic encephalopathy, Microcephaly, PRRT2 mutation",
            'Abstract': "Mutations in Proline-rich Transmembrane Protein 2 (PRRT2) have been primarily associated with individuals presenting with infantile epilepsy, including benign familial infantile epilepsy, benign infantile epilepsy, and benign myoclonus of early infancy, and/or with dyskinetic paroxysms such as paroxysmal kinesigenic dyskinesia, paroxysmal non-kinesigenic dyskinesia, and exercise-induced dyskinesia. However, the clinical manifestations of this disorder vary widely. PRRT2 encodes a protein expressed in the central nervous system that is mainly localized in the pre-synaptic neurons and is involved in the modulation of synaptic neurotransmitter release. The anomalous function of this gene has been proposed to cause dysregulation of neuronal excitability and cerebral disorders. We hereby report on a young child followed-up for three years who presents with a spectrum of clinical manifestations such as congenital microcephaly, dysmorphic features, severe intellectual disability, and drug-resistant epileptic encephalopathy in association with a synonymous variant in PRRT2 gene (c.501C > T; p.Thr167Ile) of unknown clinical significance variant (VUS) revealed by diagnostic exome sequencing. Several hypotheses have been advanced on the specific role that PRRT2 gene mutations play to cause the clinical features of affected patients. To our knowledge, the severe phenotype seen in this case has never been reported in association with any clinically actionable variant, as the missense substitution detected in PRRT2 gene. Intriguingly, the same mutation was reported in the healthy father: the action of modifying factors in the affected child may be hypothesized. The report of similar observations could extend the spectrum of clinical manifestations linked to this mutation.",
            'PMID': "31801583",
            'Link': "https://pubmed.ncbi.nlm.nih.gov/31801583/"
        },
        {
            'Abstract': "We studied the presence of benign infantile epilepsy (BIE), paroxysmal kinesigenic dyskinesia (PKD), and PKD with infantile convulsions (PKD/IC) in patients with a 16p11.2 deletion including PRRT2 or with a PRRT2 loss-of-function sequence variant. Index patients were recruited from seven Dutch university hospitals. The presence of BIE, PKD and PKD/IC was retrospectively evaluated using questionnaires and medical records. We included 33 patients with a 16p11.2 deletion: three (9%) had BIE, none had PKD or PKD/IC. Twelve patients had a PRRT2 sequence variant: BIE was present in four (p = 0.069), PKD in six (p < 0.001) and PKD/IC in two (p = 0.067). Most patients with a deletion had undergone genetic testing because of developmental problems (87%), whereas all patients with a sequence variant were tested because of a movement disorder (55%) or epilepsy (45%). BIE, PKD and PKD/IC clearly showed incomplete penetrance in patients with 16p11.2 deletions, but were found in all and 95% of patients with a PRRT2 sequence variant in our study and a large literature cohort, respectively. Deletions and sequence variants have the same underlying loss-of-function disease mechanism. Thus, differences in ascertainment have led to overestimating the frequency of BIE, PKD and PKD/IC in patients with a PRRT2 sequence variant. This has important implications for counseling if genome-wide sequencing shows such variants in patients not presenting the PRRT2-related phenotypes.",
            'Authors': "Vlaskamp Danique R M DRM, Callenbach Petra M C PMC, Rump Patrick P, Giannini Lucia A A LAA, Brilstra Eva H EH, Dijkhuizen Trijnie T, Vos Yvonne J YJ, van der Kevie-Kersemaekers Anne-Marie F AF, Knijnenburg Jeroen J, de Leeuw Nicole N, van Minkelen Rick R, Ruivenkamp Claudia A L CAL, Stegmann Alexander P A APA, Brouwer Oebele F OF, van Ravenswaaij-Arts Conny M A CMA",
            'Keywords': "Benign infantile epilepsy, Microarray, Movement disorder, Seizure, Sequencing",
            'Link': "https://pubmed.ncbi.nlm.nih.gov/30125676/",
            'PMID': "30125676",
            'Publication_year': "2019",
            'Title': "PRRT2-related phenotypes in patients with a 16p11.2 deletion."
        }
    ]
    while result_id in result_ids:
        result_id = ''.join(
            random.choice(string.ascii_lowercase + string.digits) for _ in
            range(22))
    result_ids.append(result_id)
    url = f'/results/{result_id}'
    data = json.loads(request.form['data'])
    selected_options = data['options']
    job_title = data['jobTitle']
    send_mail = data['sendMail']
    for selected_option in selected_options:
        options[selected_option] = True
    # print(results)
    # print(options)
    co_occurrence = CoOccurrence(data=test_data1, url_id=result_id,
                                 title=job_title, notify=send_mail,
                                 email=e_mail,
                                 in_sentence=options['Sentence'],
                                 in_abstract=options['Abstract'],
                                 multiple_times_in_abstract=options[
                                     'Multiple Abstracts'])
    # co_occurrence.pre_process_data()
    # print('pre processed')
    # print(co_occurrence.tokenized)
    # co_occurrence.calculate_co_occurrence()
    # print('calculated')
    # co_occurr = co_occurrence.get_co_occurence()
    # co_occurrence.save_to_db()

    # print('getting co-occurrentce')
    # print(co_occurr)
    # co_occur = {'PRRT2, Paroxysmal non-kinesigenic dyskinesia': 7,
    #             'PRRT2, Paroxysmal kinesigenic dyskinesia': 7}
    # print(co_occur)
    co_occur = ''
    return json.dumps({'status': 'OK', 'url': url, 'data': co_occur})


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']


def get_old_gen_panel_file():
    for file in os.listdir(os.path.join(os.getcwd(), 'upload')):
        if file.endswith('.txt'):
            return file


def get_algorithm_results(result_id):
    results_dict = {}
    connction = connection_database()
    cursor = connction.cursor()
    cursor.execute("select combination, amount from results join "
                   "algorithm_results ar on results.id = ar.results_id "
                   "where url_id = '{}'".format(result_id))
    for result in cursor.fetchall():
        results_dict[result[0]] = result[1]
    return results_dict


def parse_results(search):
    global results
    search.parse_results()
    search.insert_to_database()
    results = search.get_results()
    print('ready')


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


if __name__ == '__main__':
    app.run(debug=True)

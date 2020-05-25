import os
import random
import string

from flask import Flask, render_template, request, jsonify, send_file, json
from werkzeug.utils import secure_filename

from co_occurrence_algorithm.co_occurrence import CoOccurrence
from pubmed_searcher import PubmedSearch

app = Flask(__name__)

basedir = os.path.abspath(os.path.dirname(__file__))
app.config['ALLOWED_EXTENSIONS'] = {'txt'}

result_ids = []
test_data = [
    {'zoekwoord': 'google', 'aantal': 5, 'link': 'https://google.com'},
    {'zoekwoord': 'facebook', 'aantal': 3,
     'link': 'https://facebook.com'},
    {'zoekwoord': 'twitter', 'aantal': 8,
     'link': 'https://twitter.com'}
]


@app.route('/')
def home_page():
    # print(request.args.get("pheno_input"))
    # print(request.args.get("symbols_input"))
    # print("zo dus: ", request.args.get("calendar_input"))

    results = []
    articles_data = []
    if request.args.get("pheno_input"):
        term = request.args.get("pheno_input")
        words = request.args.getlist("symbols_input")
        email = request.args.get("input_mail")
        search = PubmedSearch(e_mail=email, search_word=term,
                              gene_symbols=words[0],
                              year=2015)
        search.search_pubmed()
        search.parse_results()
        results = search.results
        # search.insert_to_database()
        articles_data = search.articles_data
    print(articles_data)

    # print(request.args.getlist("focus_input"))
    return render_template('index.html', articles=articles_data,
                           results=results)


@app.route('/indextest')
def hello_world():
    global test_data
    return render_template('index.html', test_data=test_data)


@app.route('/download', methods=['GET'])
def download():
    print('mandje')
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
            old_dir = os.path.join(os.getcwd(),
                                   f'upload{os.path.sep}old_files')
            os.rename(os.path.join(os.getcwd(), file),
                      os.path.join(old_dir, old_file))
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


@app.route('/test')
def test():
    TEST_DATA = [{
        'title': 'PRRT2 gene variant in a child with dysmorphic features, congenital microcephaly, and severe epileptic seizures: genotype-phenotype correlation?',
        'abstract': 'BACKGROUND: Mutations in Proline-rich Transmembrane Protein 2 (PRRT2) have been primarily associated with individuals presenting with infantile epilepsy, including benign familial infantile epilepsy, benign infantile epilepsy, and benign myoclonus of early infancy, and/or with dyskinetic paroxysms such as paroxysmal kinesigenic dyskinesia, paroxysmal non-kinesigenic dyskinesia, and exercise-induced dyskinesia. However, the clinical manifestations of this disorder vary widely. PRRT2 encodes a protein expressed in the central nervous system that is mainly localized in the pre-synaptic neurons and is involved in the modulation of synaptic neurotransmitter release. The anomalous function of this gene has been proposed to cause dysregulation of neuronal excitability and cerebral disorders. CASE PRESENTATION: We hereby report on a young child followed-up for three years who presents with a spectrum of clinical manifestations such as congenital microcephaly, dysmorphic features, severe intellectual disability, and drug-resistant epileptic encephalopathy in association with a synonymous variant in PRRT2 gene (c.501C > T; p.Thr167Ile) of unknown clinical significance variant (VUS) revealed by diagnostic exome sequencing. CONCLUSION: Several hypotheses have been advanced on the specific role that PRRT2 gene mutations play to cause the clinical features of affected patients. To our knowledge, the severe phenotype seen in this case has never been reported in association with any clinically actionable variant, as the missense substitution detected in PRRT2 gene. Intriguingly, the same mutation was reported in the healthy father: the action of modifying factors in the affected child may be hypothesized. The report of similar observations could extend the spectrum of clinical manifestations linked to this mutation.'},
        {
            'title': 'PRRT2-related phenotypes in patients with a 16p11.2 deletion.',
            'abstract': 'We studied the presence of benign infantile epilepsy (BIE), paroxysmal kinesigenic dyskinesia (PKD), and PKD with infantile convulsions (PKD/IC) in patients with a 16p11.2 deletion including PRRT2 or with a PRRT2 loss-of-function sequence variant. Index patients were recruited from seven Dutch university hospitals. The presence of BIE, PKD and PKD/IC was retrospectively evaluated using questionnaires and medical records. We included 33 patients with a 16p11.2 deletion: three (9%) had BIE, none had PKD or PKD/IC. Twelve patients had a PRRT2 sequence variant: BIE was present in four (p = 0.069), PKD in six (p < 0.001) and PKD/IC in two (p = 0.067). Most patients with a deletion had undergone genetic testing because of developmental problems (87%), whereas all patients with a sequence variant were tested because of a movement disorder (55%) or epilepsy (45%). BIE, PKD and PKD/IC clearly showed incomplete penetrance in patients with 16p11.2 deletions, but were found in all and 95% of patients with a PRRT2 sequence variant in our study and a large literature cohort, respectively. Deletions and sequence variants have the same underlying loss-of-function disease mechanism. Thus, differences in ascertainment have led to overestimating the frequency of BIE, PKD and PKD/IC in patients with a PRRT2 sequence variant. This has important implications for counseling if genome-wide sequencing shows such variants in patients not presenting the PRRT2-related phenotypes.'}]

    test_co_occr = CoOccurrence(TEST_DATA)
    test_co_occr.pre_process_data()
    test_co_occr.calculate_co_occurrence()
    co_occurr = test_co_occr.get_co_occurence()
    print(co_occurr)
    return 'test'


@app.route('/results/<result_id>', methods=['GET'])
def render_results(result_id):
    global result_ids
    if result_id in result_ids:
        return f'hoi {result_id}'
    else:
        return 'neen'


@app.route('/results/<result_id>', methods=['POST'])
def do_algorithm(result_id):
    print(result_id)
    options = {'Title': False, 'Sentence': False, 'Abstract': False,
               'Multiple Abstracts': False}
    global result_ids
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

    data = json.loads(request.get_data())
    selected_options = data['options']
    results = data['results']
    # for selected_option in selected_options:
    #     options[selected_option] = True
    # print(results)
    # print(options)
    # co_occurrence = CoOccurrence(data=test_data1, title=options['Title'],
    #                              in_sentence=options['Sentence'],
    #                              in_abstract=options['Abstract'],
    #                              multiple_times_in_abstract=options[
    #                                  'Multiple Abstracts'])
    # co_occurrence.pre_process_data()
    # print('pre processed')
    # print(co_occurrence.tokenized)
    # co_occurrence.calculate_co_occurrence()
    # print('calculated')
    # co_occurr = co_occurrence.get_co_occurence()
    # print('getting co-occurrentce')
    # print(co_occurr)
    co_occur = {'PRRT2, Paroxysmal non-kinesigenic dyskinesia': 7, 'PRRT2, Paroxysmal kinesigenic dyskinesia': 7}
    # print(co_occur)
    return json.dumps({'status': 'OK', 'url': url, 'data': co_occur})


if __name__ == '__main__':
    app.run(debug=True)

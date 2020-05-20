from flask import Flask, render_template

from co_occurrence_algorithm import co_occurrence

app = Flask(__name__)


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/test')
def test():
    TEST_DATA = [{
        'title': 'PRRT2 gene variant in a child with dysmorphic features, congenital microcephaly, and severe epileptic seizures: genotype-phenotype correlation?',
        'abstract': 'BACKGROUND: Mutations in Proline-rich Transmembrane Protein 2 (PRRT2) have been primarily associated with individuals presenting with infantile epilepsy, including benign familial infantile epilepsy, benign infantile epilepsy, and benign myoclonus of early infancy, and/or with dyskinetic paroxysms such as paroxysmal kinesigenic dyskinesia, paroxysmal non-kinesigenic dyskinesia, and exercise-induced dyskinesia. However, the clinical manifestations of this disorder vary widely. PRRT2 encodes a protein expressed in the central nervous system that is mainly localized in the pre-synaptic neurons and is involved in the modulation of synaptic neurotransmitter release. The anomalous function of this gene has been proposed to cause dysregulation of neuronal excitability and cerebral disorders. CASE PRESENTATION: We hereby report on a young child followed-up for three years who presents with a spectrum of clinical manifestations such as congenital microcephaly, dysmorphic features, severe intellectual disability, and drug-resistant epileptic encephalopathy in association with a synonymous variant in PRRT2 gene (c.501C > T; p.Thr167Ile) of unknown clinical significance variant (VUS) revealed by diagnostic exome sequencing. CONCLUSION: Several hypotheses have been advanced on the specific role that PRRT2 gene mutations play to cause the clinical features of affected patients. To our knowledge, the severe phenotype seen in this case has never been reported in association with any clinically actionable variant, as the missense substitution detected in PRRT2 gene. Intriguingly, the same mutation was reported in the healthy father: the action of modifying factors in the affected child may be hypothesized. The report of similar observations could extend the spectrum of clinical manifestations linked to this mutation.'},
        {
            'title': 'PRRT2-related phenotypes in patients with a 16p11.2 deletion.',
            'abstract': 'We studied the presence of benign infantile epilepsy (BIE), paroxysmal kinesigenic dyskinesia (PKD), and PKD with infantile convulsions (PKD/IC) in patients with a 16p11.2 deletion including PRRT2 or with a PRRT2 loss-of-function sequence variant. Index patients were recruited from seven Dutch university hospitals. The presence of BIE, PKD and PKD/IC was retrospectively evaluated using questionnaires and medical records. We included 33 patients with a 16p11.2 deletion: three (9%) had BIE, none had PKD or PKD/IC. Twelve patients had a PRRT2 sequence variant: BIE was present in four (p = 0.069), PKD in six (p < 0.001) and PKD/IC in two (p = 0.067). Most patients with a deletion had undergone genetic testing because of developmental problems (87%), whereas all patients with a sequence variant were tested because of a movement disorder (55%) or epilepsy (45%). BIE, PKD and PKD/IC clearly showed incomplete penetrance in patients with 16p11.2 deletions, but were found in all and 95% of patients with a PRRT2 sequence variant in our study and a large literature cohort, respectively. Deletions and sequence variants have the same underlying loss-of-function disease mechanism. Thus, differences in ascertainment have led to overestimating the frequency of BIE, PKD and PKD/IC in patients with a PRRT2 sequence variant. This has important implications for counseling if genome-wide sequencing shows such variants in patients not presenting the PRRT2-related phenotypes.'}]

    test_co_occr = co_occurrence.CoOccurrence(TEST_DATA)
    test_co_occr.pre_process_data()
    test_co_occr.calculate_co_occurrence()
    co_occurr = test_co_occr.get_co_occurence()
    return co_occurr


if __name__ == '__main__':
    app.run()
